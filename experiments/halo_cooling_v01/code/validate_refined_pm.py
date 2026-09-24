#!/usr/bin/env python3
"""Validate refined periodic PM forces without running a collapse simulation.

Execution requires the SHA256 of the already frozen protocol.  Numerical
operator checks have fixed thresholds; particle-lattice and spherical-force
errors are diagnostics, not assertions of continuum accuracy.  Every execution
attempt is retained, including failures.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import time
import traceback
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import scipy

import pm_refined as refined
import spherical_collapse as sc


ROOT = Path(__file__).resolve().parents[1]
TOLERANCES = {
    "mass_absolute": 1e-12,
    "uniform_force_absolute": 1e-12,
    "self_force_absolute": 1e-12,
    "net_force_normalized": 1e-11,
    "rfft_fullfft_relative": 1e-11,
    "direct_dft_relative": 1e-11,
}


def file_sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def json_safe(value):
    """Retain a failed nonfinite value as text instead of emitting invalid JSON."""
    if isinstance(value, dict):
        return {str(key): json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(item) for item in value]
    if isinstance(value, np.ndarray):
        return json_safe(value.tolist())
    if isinstance(value, (np.bool_,)):
        return bool(value)
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (float, np.floating)):
        return float(value) if np.isfinite(value) else str(value)
    return value


def particle_centers(nparticle):
    axis = (np.arange(nparticle, dtype=float) + 0.5) / nparticle
    return np.stack(np.meshgrid(axis, axis, axis, indexing="ij"), axis=-1).reshape(-1, 3)


def relative_l2(actual, expected):
    denominator = np.linalg.norm(expected)
    if denominator == 0:
        return float(np.linalg.norm(actual - expected))
    return float(np.linalg.norm(actual - expected) / denominator)


def classify_checks(checks):
    """Classify unchanged outcomes; a diagnostic failure remains a failure.

    The compensated-sphere sign/finite checks describe physical discretization.
    They are not frozen implementation-accuracy gates.  Nonfinite force errors
    raised by the solver still appear as implementation execution exceptions.
    """
    scopes = ("implementation_protocol", "physical_discretization_diagnostic")
    for check in checks:
        inferred = scopes[1] if check["name"].startswith("compensated_core_") else scopes[0]
        check.setdefault("scope", inferred)
        if check["scope"] not in scopes:
            raise ValueError(f"Unrecognized validation scope: {check['scope']}")
    implementation = [check for check in checks if check["scope"] == scopes[0]]
    physical = [check for check in checks if check["scope"] == scopes[1]]
    implementation_failures = [check for check in implementation if not check["passed"]]
    physical_failures = [check for check in physical if not check["passed"]]
    return {
        "implementation_status": "FAIL" if implementation_failures else "PASS",
        "implementation_failures": implementation_failures,
        "implementation_checks_passed": len(implementation) - len(implementation_failures),
        "implementation_checks_total": len(implementation),
        "physical_diagnostic_status": "FAIL" if physical_failures else "PASS",
        "physical_diagnostic_failures": physical_failures,
        "physical_diagnostic_checks_passed": len(physical) - len(physical_failures),
        "physical_diagnostic_checks_total": len(physical),
    }


def explicit_cic_mass(positions, nmesh):
    """Independent small-grid scalar loop, used only by the direct DFT check."""
    mass = np.zeros((nmesh,) * 3, dtype=float)
    for particle in positions:
        scaled = np.mod(particle, 1.0) * nmesh
        lower = np.floor(scaled).astype(int)
        fraction = scaled - lower
        for dx in range(2):
            for dy in range(2):
                for dz in range(2):
                    offsets = (dx, dy, dz)
                    weight = 1.0
                    for axis, offset in enumerate(offsets):
                        weight *= fraction[axis] if offset else 1.0 - fraction[axis]
                    index = tuple((lower[axis] + offsets[axis]) % nmesh for axis in range(3))
                    mass[index] += weight / len(positions)
    return mass


def explicit_cic_interpolate(grid, positions):
    nmesh = grid.shape[0]
    values = np.zeros((len(positions), 3), dtype=float)
    for number, particle in enumerate(positions):
        scaled = np.mod(particle, 1.0) * nmesh
        lower = np.floor(scaled).astype(int)
        fraction = scaled - lower
        for dx in range(2):
            for dy in range(2):
                for dz in range(2):
                    offsets = (dx, dy, dz)
                    weight = 1.0
                    for axis, offset in enumerate(offsets):
                        weight *= fraction[axis] if offset else 1.0 - fraction[axis]
                    index = tuple((lower[axis] + offsets[axis]) % nmesh for axis in range(3))
                    values[number] += weight * grid[index]
    return values


def direct_dft_force(delta):
    """Explicit finite Fourier sums: no FFT or production Poisson helper.

    The reference is the declared discrete spectral Poisson operator, including
    its zero DC and zero component derivative on each Nyquist plane.  Agreement
    with this sum does not establish agreement with continuum gravity.
    """
    nmesh = delta.shape[0]
    if delta.shape != (nmesh,) * 3 or nmesh > 8 or nmesh % 2:
        raise ValueError("Direct DFT reference is restricted to small even cubic grids")
    indices = np.arange(nmesh, dtype=int)
    integers = np.where(indices < nmesh // 2, indices, indices - nmesh)
    coordinates = np.stack(np.meshgrid(indices / nmesh, indices / nmesh,
                                      indices / nmesh, indexing="ij"), axis=-1).reshape(-1, 3)
    modes = np.stack(np.meshgrid(integers, integers, integers,
                                indexing="ij"), axis=-1).reshape(-1, 3)
    forward = np.exp(-2j * np.pi * (coordinates @ modes.T))
    coefficients = forward.T @ np.asarray(delta).ravel()
    k = 2 * np.pi * modes
    k2 = np.sum(k * k, axis=1)
    inverse_k2 = np.zeros_like(k2)
    inverse_k2[k2 > 0] = 1.0 / k2[k2 > 0]
    force = np.empty((nmesh**3, 3), dtype=float)
    imaginary = 0.0
    for axis in range(3):
        gradient = k[:, axis].copy()
        gradient[modes[:, axis] == -nmesh // 2] = 0.0
        transformed = forward.conj() @ (1j * gradient * inverse_k2 * coefficients) / nmesh**3
        force[:, axis] = transformed.real
        imaginary = max(imaginary, float(np.max(np.abs(transformed.imag))))
    return force.reshape((nmesh,) * 3 + (3,)), imaginary


def run_checks(checks, detail, workers=4, large=False):
    def upper(name, value, tolerance_name):
        threshold = TOLERANCES[tolerance_name]
        checks.append({"name": name, "passed": bool(np.isfinite(value) and value <= threshold),
                       "value": float(value), "inclusive_upper_threshold": threshold,
                       "tolerance_name": tolerance_name})

    def predicate(name, condition, **extra):
        checks.append({"name": name, "passed": bool(condition), **extra})

    rng = np.random.default_rng(20260924)
    detail["interpretation"] = {
        "sinc_reference_used_for_unequal_grids": False,
        "physical_force_errors": "Diagnostic only; no machine-precision continuum claim",
        "lattice_aliases": "Uniform particles may produce nonzero mesh delta when Nforce>Nparticle",
        "twofold_force_grid": "For a regular load, all Bragg gradient components are DC/Nyquist, hence mesh force is zero with this operator",
        "fourfold_force_grid": "Non-Nyquist Bragg modes may produce force between particles; particle forces still cancel by symmetry",
        "phase_response": "A force-only test, not a fitted growth reference or a collapse trajectory",
    }

    # Check the revised transform against the unchanged full-complex transform.
    fft_rows = []
    for nmesh in (8, 16, 32):
        mesh = refined.PMGrid(nmesh, workers=workers)
        old_mesh = refined.legacy.PMGrid(nmesh)
        delta = rng.standard_normal((nmesh,) * 3)
        delta -= delta.mean()
        actual, imaginary = mesh.mesh_force(delta)
        expected, old_imaginary = old_mesh.mesh_force(delta)
        error = relative_l2(actual, expected)
        upper(f"rfft_fullfft_mesh:nf={nmesh}", error, "rfft_fullfft_relative")
        random_positions = rng.uniform(-1.0, 2.0, size=(257, 3))
        random_positions[:4] = [[0, 0, 0], [1, 1, 1], [-1e-15, 0.5, 0.25], [1 + 1e-15, 0, 0.75]]
        actual_particle, diagnostics = mesh.force(random_positions)
        expected_particle, _ = old_mesh.force(random_positions)
        particle_error = relative_l2(actual_particle, expected_particle)
        upper(f"rfft_fullfft_particles:nf={nmesh}", particle_error, "rfft_fullfft_relative")
        mass = refined.cic_deposit(random_positions, nmesh)
        upper(f"random_raw_mass:nf={nmesh}", abs(float(mass.sum()) - 1), "mass_absolute")
        force_scale = float(np.mean(np.linalg.norm(actual_particle, axis=1)))
        net_ratio = float(np.linalg.norm(np.mean(actual_particle, axis=0)) / max(force_scale, 1e-300))
        upper(f"random_net_force:nf={nmesh}", net_ratio, "net_force_normalized")
        for number, position in enumerate(((0.5, 0.5, 0.5), (0.17351, 0.29783, 0.73111),
                                           (0.5 / nmesh,) * 3)):
            self_force, _ = mesh.force(np.asarray([position], dtype=float))
            upper(f"single_self_force:nf={nmesh}:case={number}",
                  np.max(np.abs(self_force)), "self_force_absolute")
        fft_rows.append({"nforce": nmesh, "mesh_relative_l2": error,
                         "particle_relative_l2": particle_error,
                         "rfft_imaginary_diagnostic": imaginary,
                         "fullfft_imaginary_diagnostic": old_imaginary,
                         "random_diagnostics": diagnostics})
    detail["rfft_fullfft"] = fft_rows

    # This reference uses explicit Fourier sums and separate scalar CIC loops.
    direct_n = 8
    mesh = refined.PMGrid(direct_n, workers=workers)
    delta = rng.standard_normal((direct_n,) * 3)
    delta -= delta.mean()
    direct, direct_imaginary = direct_dft_force(delta)
    actual, _ = mesh.mesh_force(delta)
    direct_error = relative_l2(actual, direct)
    upper("explicit_DFT_mesh:nf=8", direct_error, "direct_dft_relative")
    positions = rng.uniform(-0.3, 1.3, size=(17, 3))
    positions[:3] = [[0, 0, 0], [1, 0.5, 0.25], [-1e-15, 0.25, 0.75]]
    explicit_mass = explicit_cic_mass(positions, direct_n)
    explicit_grid_force, particle_direct_imaginary = direct_dft_force(explicit_mass * direct_n**3 - 1)
    direct_particle = explicit_cic_interpolate(explicit_grid_force, positions)
    actual_particle, _ = mesh.force(positions)
    direct_particle_error = relative_l2(actual_particle, direct_particle)
    upper("explicit_DFT_and_CIC_particles:nf=8", direct_particle_error, "direct_dft_relative")
    detail["explicit_dft"] = {"nforce": direct_n, "mesh_relative_l2": direct_error,
                              "particles_relative_l2": direct_particle_error,
                              "mesh_imaginary_max": direct_imaginary,
                              "particles_mesh_imaginary_max": particle_direct_imaginary,
                              "reference": "Explicit forward/inverse DFT sums; independent scalar-loop CIC for particles"}

    shifts = ((0.0, 0.0, 0.0), (0.25, 0.375, 0.5), (0.5, 0.5, 0.5))
    configurations = [(8, 8), (8, 16), (8, 32)]
    if large:
        configurations += [(64, 128), (64, 256)]
    uniform_rows = []
    for nparticle, nmesh in configurations:
        mesh = refined.PMGrid(nmesh, workers=workers)
        q = particle_centers(nparticle)
        for shift in shifts:
            x = np.mod(q + np.asarray(shift) / nmesh, 1.0)
            mass = refined.cic_deposit(x, nmesh)
            delta = mass * nmesh**3 - 1
            grid_force, imaginary = mesh.mesh_force(delta)
            force = refined.cic_gather(grid_force, x)
            mass_error = abs(float(mass.sum()) - 1.0)
            force_error = float(np.max(np.abs(force)))
            label = f"np={nparticle}:nf={nmesh}:shift={shift}"
            upper(f"uniform_raw_mass:{label}", mass_error, "mass_absolute")
            upper(f"uniform_particle_force:{label}", force_error, "uniform_force_absolute")
            uniform_rows.append({"nparticle": nparticle, "nforce": nmesh,
                                 "shift_force_cells": shift, "mass_absolute_error": mass_error,
                                 "mesh_delta_rms": float(np.sqrt(np.mean(delta**2))),
                                 "mesh_force_max_abs": float(np.max(np.abs(grid_force))),
                                 "particle_force_max_abs": force_error,
                                 "imaginary_diagnostic": imaginary})
            del mass, delta, grid_force, force
        del mesh
    detail["uniform_lattices"] = uniform_rows

    # Perturbation reversal and translation expose the CIC node/alias issue.
    # No sinc value, fitted growth exponent, or amplitude-dependent pass bound.
    response_rows = []
    q = particle_centers(8)
    for nmesh in (8, 16, 32):
        mesh = refined.PMGrid(nmesh, workers=workers)
        for shift in shifts:
            x0 = np.mod(q + np.asarray(shift) / nmesh, 1.0)
            baseline, _ = mesh.force(x0)
            for mode in ((1, 0, 0), (1, 1, 0), (1, 1, 1)):
                direction = np.asarray(mode, dtype=float)
                direction /= np.linalg.norm(direction)
                shape = np.sin(2 * np.pi * (q @ np.asarray(mode)))[:, None] * direction
                for fraction in (1.0, 0.5):
                    amplitude = fraction * 2e-5 / nmesh
                    displacement = amplitude * shape
                    positive, _ = mesh.force(np.mod(x0 + displacement, 1.0))
                    negative, _ = mesh.force(np.mod(x0 - displacement, 1.0))
                    positive -= baseline
                    negative -= baseline
                    gain = float(np.sum(positive * displacement) / np.sum(displacement**2))
                    residual = positive - gain * displacement
                    transverse = positive - (positive @ direction)[:, None] * direction
                    nonodd = float(np.linalg.norm(positive + negative)
                                   / max(np.linalg.norm(positive) + np.linalg.norm(negative), 1e-300))
                    row = {"nparticle": 8, "nforce": nmesh, "shift_force_cells": shift,
                           "mode": mode, "amplitude_box": amplitude,
                           "force_projection_on_displacement": gain,
                           "orthogonal_residual_relative_l2": float(np.linalg.norm(residual)
                                                                      / max(np.linalg.norm(positive), 1e-300)),
                           "transverse_force_relative_l2": float(np.linalg.norm(transverse)
                                                                   / max(np.linalg.norm(positive), 1e-300)),
                           "reversal_nonodd_relative_l2": nonodd,
                           "interpretation": "Diagnostic of finite CIC response; not a universal analytic growth factor"}
                    predicate(f"phase_response_finite:nf={nmesh}:shift={shift}:mode={mode}:amplitude={amplitude}",
                              np.all(np.isfinite([gain, row["orthogonal_residual_relative_l2"],
                                                 row["transverse_force_relative_l2"], nonodd])))
                    response_rows.append(row)
    detail["phase_response"] = response_rows

    # Uniform compensated core: shell theorem gives g=-(delta_nl_core/3)*r.
    # This tests physical discretization separately from transform equivalence.
    reference_background = sc.load_backgrounds()["reference"]
    detail["sphere_background_input"] = {"path": str(reference_background.path),
                                         "sha256": file_sha256(reference_background.path)}
    sphere_rows = []
    sphere_baseline = None
    sphere_configs = [(32, 32, (0.0, 0.0, 0.0)), (32, 64, (0.0, 0.0, 0.0)),
                      (32, 128, (0.0, 0.0, 0.0)), (32, 128, (0.25, 0.375, 0.5))]
    for nparticle, nmesh, shift in sphere_configs:
        initial = refined.compensated_initial(nparticle, delta_i=0.05,
                                               reference_background=reference_background,
                                               shift_cells=shift, force_mesh=nmesh)
        x = np.asarray(initial["positions"], dtype=float)
        q = np.asarray(initial["q"], dtype=float)
        q_radius = np.asarray(initial["q_radius"], dtype=float)
        center = np.asarray(initial["center"], dtype=float)
        metadata = initial["metadata"]
        core_radius = float(metadata["core_lagrangian_radius_box"])
        core_delta = float(metadata["delta_nl_core"])
        scale_core = float(metadata["scale_core"])
        radial_vector = np.mod(x - center + 0.5, 1.0) - 0.5
        selection = (q_radius >= 0.25 * core_radius) & (q_radius <= 0.75 * core_radius)
        if not np.any(selection):
            raise RuntimeError("No particles in the predeclared analytic core-force interval")
        mesh = refined.PMGrid(nmesh, workers=workers)
        actual, diagnostics = mesh.force(x)
        expected = -(core_delta / 3) * radial_vector
        selected_force = actual[selection]
        selected_radial = radial_vector[selection]
        selected_radius = np.linalg.norm(selected_radial, axis=1)
        unit_radial = selected_radial / selected_radius[:, None]
        radial_component = np.sum(selected_force * unit_radial, axis=1)
        tangent = selected_force - radial_component[:, None] * unit_radial
        expected_norm = float(np.linalg.norm(expected[selection]))
        inward_fraction = float(np.mean(radial_component < 0))
        mean_radial_component = float(np.mean(radial_component))
        mapping_error = float(np.max(np.abs(radial_vector[selection] - scale_core * q[selection])))
        label = f"np={nparticle}:nf={nmesh}:shift={shift}"
        predicate(f"compensated_core_force_finite:{label}", np.all(np.isfinite(selected_force)))
        predicate(f"compensated_core_mean_force_inward:{label}",
                  np.isfinite(mean_radial_component) and mean_radial_component < 0,
                  mean_radial_component=mean_radial_component,
                  particles_with_inward_force_fraction=inward_fraction,
                  qualification="Only the mean sign is a gate; individual-force accuracy is diagnostic")
        row = {"nparticle": nparticle, "nforce": nmesh, "shift_force_cells": shift,
               "delta_nl_core": core_delta, "scale_core": scale_core,
               "core_lagrangian_radius_box": core_radius,
               "core_q_interval_fraction": [0.25, 0.75], "selected_particles": int(selection.sum()),
               "core_mapping_max_abs_box": mapping_error,
               "force_relative_l2_to_continuum_sphere": relative_l2(selected_force, expected[selection]),
               "tangential_force_relative_l2_to_analytic": float(np.linalg.norm(tangent) / expected_norm),
               "radial_force_relative_l2_to_analytic": relative_l2(radial_component,
                                                                   -(core_delta / 3) * selected_radius),
               "mean_radial_component": mean_radial_component,
               "particles_with_inward_force_fraction": inward_fraction,
               "diagnostics": diagnostics,
               "interpretation": "Finite particle/mesh error, diagnostic only; not evidence of a resolved halo"}
        if nmesh == 128 and shift == (0.0, 0.0, 0.0):
            sphere_baseline = selected_force.copy()
        elif nmesh == 128:
            row["force_relative_l2_to_unshifted_configuration"] = relative_l2(selected_force, sphere_baseline)
        sphere_rows.append(row)
    detail["compensated_sphere"] = sphere_rows


def write_report(result):
    lines = ["# 细化 PM 网格的独立验证", "",
             f"状态：**{result['status']}，{result['checks_passed']}/{result['checks_total']} 项通过**。",
             f"第 {result['attempt']} 次记录；历次尝试及失败保留。", "",
             f"实现检查：**{result['implementation_status']}，{result['implementation_checks_passed']}/{result['implementation_checks_total']} 项通过**；物理离散诊断：**{result['physical_diagnostic_status']}，{result['physical_diagnostic_checks_passed']}/{result['physical_diagnostic_checks_total']} 项通过**。", "",
             f"冻结协议 SHA256：`{result['protocol_sha256']}`。", "",
             "本脚本只检查力算子与初始条件，不执行晕坍缩生产计算。RFFT 与旧 fullFFT、显式 DFT 的相符检查，验证的是声明的离散算子；它们不代表连续物理误差达到机器精度。", "",
             "质量、自力、净力和规则粒子处力使用预先规定阈值。补偿球的连续力误差、粒子格相位误差单独记录，不据结果重设精度门槛。", "",
             "每项检查以 scope 区分 implementation_protocol 与 physical_discretization_diagnostic。总状态保留所有失败；生产脚本仅以顶层 implementation_failures 是否为空判断力算法实现检查是否阻断。物理离散诊断的失败仍是失败，必须与随后六个固定配置的收敛结果一并解释。算法实现通过不能视为物理收敛。", "",
             "| 检查 | 类别 | 状态 | 数值或说明 |", "| --- | --- | --- | --- |"]
    for check in result["checks"]:
        value = check.get("value", check.get("mean_radial_component", "布尔条件，详见 JSON"))
        rendered = f"{value:.8g}" if isinstance(value, (int, float)) else str(value)
        lines.append(f"| {check['name']} | {check['scope']} | {'PASS' if check['passed'] else 'FAIL'} | {rendered} |")
    lines += ["", "## 规则粒子格与力网格", "",
              "规则粒子格在更细网格上通常产生非零沉积密度，并非加入了真实物理扰动。粒子处的力仍应由对称性抵消。二倍网格的 Bragg 模只有零或 Nyquist 分量，本实现的 Nyquist 安全梯度使基态网格力为零；四倍网格则可以在粒子之间存在网格力。", "",
              "| 每方向粒子数 | 力网格 | 平移（力网格单元） | 网格密度 RMS | 网格力最大值 | 粒子力最大值 |",
              "| --- | --- | --- | --- | --- | --- |"]
    for row in result["detail"].get("uniform_lattices", []):
        lines.append(f"| {row['nparticle']} | {row['nforce']} | {row['shift_force_cells']} | {row['mesh_delta_rms']:.6g} | {row['mesh_force_max_abs']:.6g} | {row['particle_force_max_abs']:.6g} |")
    lines += ["", "## CIC 相位与别名响应", "",
              "此处未把粒子居中且与网格一一对应时的 sinc 公式用于不等分辨率。JSON 保留正负小扰动、两种幅度及三种整体平移下的力投影、非奇对称部分与横向分量。CIC 节点的分段导数和粒子格别名都可能导致离散响应差异；这些结果不能反拟合成解析增长参照后再用于自证。", "",
              "## 补偿球的连续物理参照", "",
              "均匀核心使用实际初始非线性过密度 δ=0.05；其严格映射为 r=q/(1+δ)^(1/3)，解析初始力为 g=−δr/3。固定选择初始拉格朗日半径 q/qL∈[0.25,0.75]，避免球心及核心边缘。有限性及平均径向力向内保留为物理诊断的通过／失败条件，不作为算法实现通过线；下表的离散力误差是诊断值，未预设它应达到高精度或随单独细化力网格单调改善。", "",
              "| 粒子／力网格 | 平移 | 核心粒子数 | 对连续球力相对 L2 | 横向力相对 L2 | 个体力向内比例 | 映射最大误差 |",
              "| --- | --- | --- | --- | --- | --- | --- |"]
    for row in result["detail"].get("compensated_sphere", []):
        lines.append(f"| {row['nparticle']}³ / {row['nforce']}³ | {row['shift_force_cells']} | {row['selected_particles']} | {row['force_relative_l2_to_continuum_sphere']:.6g} | {row['tangential_force_relative_l2_to_analytic']:.6g} | {row['particles_with_inward_force_fraction']:.6g} | {row['core_mapping_max_abs_box']:.6g} |")
    lines += ["", "## 失败记录", ""]
    if result["failures"]:
        lines += [f"- {row['name']}：{json.dumps(json_safe(row), ensure_ascii=False)}" for row in result["failures"]]
    else:
        lines.append("本次无未通过检查。既往尝试文件及哈希保留在 JSON 中。")
    lines += ["", "细化力网格不会增加粒子采样信息。本验证不等同于晕质量、密度剖面、气体冷却或星系形成收敛；后续仍需独立改变粒子数、力网格、时间步与网格相位。", ""]
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--protocol-sha256", required=True,
                        help="SHA256 supplied only after protocol.md has been frozen")
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--large", action="store_true",
                        help="Also check uniform 64^3 particles on 128^3 and 256^3 force grids")
    args = parser.parse_args()
    if args.workers < 1:
        parser.error("--workers must be positive")
    expected_hash = args.protocol_sha256.lower()
    if len(expected_hash) != 64 or any(character not in "0123456789abcdef" for character in expected_hash):
        parser.error("--protocol-sha256 must contain exactly 64 hexadecimal characters")
    protocol_path = ROOT / "protocol.md"
    actual_hash = file_sha256(protocol_path)
    if actual_hash != expected_hash:
        raise RuntimeError("Protocol hash mismatch; refusing to run an unfrozen or changed protocol")

    started = time.perf_counter()
    result_directory = ROOT / "results"
    report_directory = ROOT / "reports"
    result_directory.mkdir(parents=True, exist_ok=True)
    report_directory.mkdir(parents=True, exist_ok=True)
    attempts = sorted(result_directory.glob("refined_pm_validation_attempt_*.json"))
    attempt = max([int(path.stem.rsplit("_", 1)[1]) for path in attempts], default=0) + 1
    checks, detail = [], {}
    try:
        run_checks(checks, detail, workers=args.workers, large=args.large)
    except Exception:
        checks.append({"name": "execution_exception", "passed": False,
                       "detail": traceback.format_exc()})
    failures = [row for row in checks if not row["passed"]]
    classification = classify_checks(checks)
    code_paths = {"pm_refined.py": ROOT / "code" / "pm_refined.py",
                  "validate_refined_pm.py": Path(__file__),
                  "spherical_collapse.py": ROOT / "code" / "spherical_collapse.py",
                  "legacy_pm.py": Path(refined.legacy.__file__)}
    result = {"status": "FAIL" if failures else "PASS", "attempt": attempt,
              "created_utc": datetime.now(timezone.utc).isoformat(),
              "checks_passed": len(checks) - len(failures), "checks_total": len(checks),
              "checks": checks, "failures": failures, "detail": detail,
              **classification,
              "scope_classification": {
                  "schema_version": 1,
                  "rule": "compensated_core_* checks are physical_discretization_diagnostic; all remaining checks are implementation_protocol",
                  "production_gate": "implementation_failures must be empty; status retains every failure",
                  "outcomes_changed_by_classification": False,
              },
              "protocol_sha256": actual_hash, "thresholds": TOLERANCES,
              "command_configuration": {"workers": args.workers, "large": args.large},
              "code_sha256": {name: file_sha256(path) for name, path in code_paths.items()},
              "environment": {"python": platform.python_version(), "numpy": np.__version__,
                              "scipy": scipy.__version__},
              "previous_attempts": [{"path": path.name, "sha256": file_sha256(path)} for path in attempts],
              "wall_seconds": time.perf_counter() - started,
              "scope": "Force-operator checks and initial compensated-sphere discretization; no production collapse or halo claim"}
    encoded = json.dumps(json_safe(result), indent=2, ensure_ascii=False, allow_nan=False) + "\n"
    attempt_path = result_directory / f"refined_pm_validation_attempt_{attempt:03d}.json"
    with attempt_path.open("x") as handle:
        handle.write(encoded)
    (result_directory / "refined_pm_validation.json").write_text(encoded)
    (report_directory / "refined_pm_validation_cn.md").write_text(write_report(result))
    print(json.dumps(json_safe({"status": result["status"], "attempt": attempt,
                               "checks_passed": result["checks_passed"],
                               "checks_total": result["checks_total"],
                               "implementation_status": result["implementation_status"],
                               "implementation_failures": result["implementation_failures"],
                               "physical_diagnostic_status": result["physical_diagnostic_status"],
                               "physical_diagnostic_failures": result["physical_diagnostic_failures"],
                               "wall_seconds": result["wall_seconds"], "failures": failures}),
                     ensure_ascii=False, indent=2), flush=True)
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
