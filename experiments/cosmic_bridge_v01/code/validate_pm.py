#!/usr/bin/env python3
"""Execute the frozen small PM checks; preserve every attempt and failure."""

from __future__ import annotations

import json
import platform
import time
import traceback
import tempfile
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import scipy

from pm import (Background, PMGrid, advance_kdk, cic_deposit, file_sha256,
                generate_initial_modes, initial_conditions, integrate)


ROOT = Path(__file__).resolve().parents[1]
PROTOCOL_SHA256 = "74a2cab351d53f98169b885542d23d18b4d425c6fd6b75b91b49f5c53422e2e1"


def plane_initial(n, epsilon=2e-6, a_initial=0.02):
    centers = (np.arange(n) + 0.5) / n
    q = np.stack(np.meshgrid(centers, centers, centers, indexing="ij"), axis=-1).reshape((-1, 3))
    displacement = np.zeros_like(q)
    displacement[:, 0] = -epsilon / (2 * np.pi) * np.sin(2 * np.pi * q[:, 0])
    return {"positions": np.remainder(q + displacement, 1.0),
            "momenta": np.sqrt(a_initial) * displacement,
            "lagrangian": q, "displacement": displacement,
            "metadata": {"a_initial": a_initial, "epsilon_initial": epsilon, "f_initial": 1.0,
                         "E_initial": a_initial**(-1.5), "initial_condition": "single x-axis fundamental mode"}}


def pm_growth_reference(n, ratio=50.0):
    kh = 2 * np.pi / n
    response = np.sin(kh) / kh
    rplus = (-1 + np.sqrt(1 + 24 * response)) / 4
    rminus = (-1 - np.sqrt(1 + 24 * response)) / 4
    growing = (1 - rminus) / (rplus - rminus)
    decaying = (rplus - 1) / (rplus - rminus)
    value = growing * ratio**rplus + decaying * ratio**rminus
    return {"F": float(response), "r_plus": float(rplus), "r_minus": float(rminus),
            "growing_weight": float(growing), "decaying_weight": float(decaying),
            "growth": float(value), "continuum_growth": float(ratio),
            "derivation": "CIC center load, node mesh, spectral gradient: F=sin(kh)/(kh); no fitted force correction"}


def run_checks(checks, detail):
    def upper(name, value, threshold, **extra):
        checks.append({"name": name, "passed": bool(np.isfinite(value) and value < threshold),
                       "value": float(value), "strict_upper_threshold": threshold, **extra})

    rng = np.random.default_rng(20260924)
    # Additive pre-production input checks requested by internal review; no
    # scientific parameter or frozen numerical threshold is changed.
    with tempfile.TemporaryDirectory(prefix="rsm_pm_input_check_") as temporary:
        path = Path(temporary) / "nonfinite_power.csv"
        path.write_text("k_h_mpc,pk_z49_mpc_over_h3\n0.001,1\n1,nan\n10,1\n")
        try:
            generate_initial_modes(path)
            rejected = False
        except ValueError:
            rejected = True
        checks.append({"name": "reject_nonfinite_power_input", "passed": rejected})
    for condition in ("nonfinite", "nonhermitian"):
        coefficient = np.zeros((4, 4, 4), dtype=complex)
        coefficient[1, 0, 0] = np.nan if condition == "nonfinite" else 1.0
        try:
            initial_conditions({"coefficients": coefficient, "metadata": {}}, 4, 0.02, 1.0, 1.0)
            rejected = False
        except ValueError:
            rejected = True
        checks.append({"name": f"reject_{condition}_mode_input", "passed": rejected})
    for n in (16, 32):
        mesh = PMGrid(n)
        q = plane_initial(n)["lagrangian"]
        for shift in (0.0, 0.137):
            force, diag = mesh.force(np.remainder(q + shift, 1.0))
            upper(f"uniform_force:n={n}:shift={shift}", np.max(np.abs(force)), 1e-12)
            upper(f"uniform_raw_mass:n={n}:shift={shift}", diag["mass_relative_error"], 1e-12)
        random_positions = rng.uniform(-1.5, 2.5, size=(257, 3))
        random_positions[:4] = [[0, 0, 0], [1, 1, 1], [-1e-15, 0.5, 0.25], [1 + 1e-15, 0, 0.75]]
        force, diag = mesh.force(random_positions)
        upper(f"random_raw_mass:n={n}", diag["mass_relative_error"], 1e-12)
        upper(f"random_total_force:n={n}", diag["total_force_normalized"], 1e-12)
        single = np.array([[0.17351, 0.29783, 0.73111]])
        mass = cic_deposit(single, n)
        force_single, _ = mesh.force(single)
        grid_force, _ = mesh.mesh_force(mass * n**3 - 1)
        scale = float(np.max(np.linalg.norm(grid_force, axis=-1)))
        upper(f"single_particle_self_force:n={n}", np.linalg.norm(force_single) / scale, 1e-12)
        grid = np.stack(np.meshgrid(*([np.arange(n) / n] * 3), indexing="ij"), axis=-1)
        for index in ((1, 0, 0), (1, 2, 0), (2, 1, 3)):
            k = 2 * np.pi * np.asarray(index, dtype=float)
            phase = np.einsum("...i,i->...", grid, k)
            delta = np.cos(phase)
            actual, _ = mesh.mesh_force(delta)
            expected = -np.sin(phase)[..., None] * k / np.dot(k, k)
            relative = np.linalg.norm(actual - expected) / np.linalg.norm(expected)
            upper(f"analytic_mesh_fourier_force:n={n}:mode={index}", relative, 1e-12)

    eds = Background.eds()
    for a1, a2 in ((0.02, 0.04), (0.02, 1.0), (0.73, 1.0), (0.1, 0.02)):
        kick, drift = eds.factors(a1, a2)
        expected_kick = 2 * (np.sqrt(a2) - np.sqrt(a1))
        expected_drift = 2 * (a1**(-0.5) - a2**(-0.5))
        upper(f"eds_kick_factor:{a1}:{a2}", abs(kick / expected_kick - 1), 1e-11)
        upper(f"eds_drift_factor:{a1}:{a2}", abs(drift / expected_drift - 1), 1e-11)

    initial = plane_initial(16)
    mesh = PMGrid(16)
    x, p = initial["positions"].copy(), initial["momenta"].copy()
    edges = np.geomspace(0.02, 0.024, 9)
    for a1, a2 in zip(edges[:-1], edges[1:]):
        x, p, _, _ = advance_kdk(x, p, a1, a2, eds, mesh)
    for a2, a1 in zip(edges[:0:-1], edges[-2::-1]):
        x, p, _, _ = advance_kdk(x, p, a2, a1, eds, mesh)
    displacement = np.remainder(x - initial["positions"] + 0.5, 1.0) - 0.5
    upper("short_time_reversal_position_absolute", np.max(np.abs(displacement)), 1e-10)
    upper("short_time_reversal_momentum_absolute", np.max(np.abs(p - initial["momenta"])), 1e-10)

    plane_results = {}
    for n in (16, 32):
        initial = plane_initial(n)
        reference = pm_growth_reference(n)
        values = []
        for steps in (128, 256, 512):
            started = time.perf_counter()
            result = integrate(initial, eds, n, steps=steps, a_final=1, a_outputs=())
            s = np.remainder(result["final_positions"] - initial["lagrangian"] + 0.5, 1.0) - 0.5
            s0 = initial["displacement"]
            growth = float(np.sum(s * s0) / np.sum(s0**2))
            residual = s - growth * s0
            transverse = float(np.max(np.abs(s[:, 1:])))
            value = {"steps": steps, "growth": growth,
                     "relative_error_discrete_linear_reference": abs(growth / reference["growth"] - 1),
                     "relative_error_continuum_EdS": abs(growth / 50.0 - 1),
                     "nonfundamental_displacement_relative_l2": float(np.linalg.norm(residual) / np.linalg.norm(s)),
                     "transverse_displacement_max_abs": transverse,
                     "diagnostics": result["diagnostics"], "wall_seconds": time.perf_counter() - started}
            values.append(value)
            upper(f"plane_mass:n={n}:steps={steps}", result["diagnostics"]["max_mass_relative_error"], 1e-12)
            print(json.dumps({"stage": "EdS plane", "nmesh": n, **value}, ensure_ascii=False), flush=True)
        upper(f"fine_step_discrete_growth:n={n}", values[-1]["relative_error_discrete_linear_reference"], 5e-4)
        difference1 = abs(values[0]["growth"] - values[1]["growth"])
        difference2 = abs(values[1]["growth"] - values[2]["growth"])
        order = float(np.log2(difference1 / difference2))
        checks.append({"name": f"KDK_time_difference_order:n={n}", "passed": bool(1.7 <= order <= 2.3),
                       "value": order, "inclusive_interval": [1.7, 2.3],
                       "coarse_medium_difference": difference1, "medium_fine_difference": difference2})
        plane_results[str(n)] = {"reference": reference, "runs": values, "time_order": order}
    detail["plane_mode"] = plane_results
    error16 = plane_results["16"]["runs"][-1]["relative_error_continuum_EdS"]
    error32 = plane_results["32"]["runs"][-1]["relative_error_continuum_EdS"]
    detail["spatial_comparison"] = {
        "continuum_relative_error_16": error16, "continuum_relative_error_32": error32,
        "error_reduction": error16 / error32,
        "qualification": "Two resolutions diagnose the expected CIC spatial bias; they do not establish a general convergence theorem."}


def main():
    started = time.perf_counter()
    result_dir = ROOT / "results"
    report_dir = ROOT / "reports"
    result_dir.mkdir(parents=True, exist_ok=True)
    report_dir.mkdir(parents=True, exist_ok=True)
    actual_protocol_hash = file_sha256(ROOT / "protocol.md")
    if actual_protocol_hash != PROTOCOL_SHA256:
        raise RuntimeError("Frozen protocol hash changed; refusing an unrecorded protocol mutation")
    attempts = sorted(result_dir.glob("pm_validation_attempt_*.json"))
    attempt = len(attempts) + 1
    checks, detail = [], {}
    failure_exception = None
    try:
        run_checks(checks, detail)
    except Exception:
        failure_exception = traceback.format_exc()
        checks.append({"name": "execution_exception", "passed": False, "detail": failure_exception})
    failures = [row for row in checks if not row["passed"]]
    result = {"status": "PASS" if not failures else "FAIL", "attempt": attempt,
              "created_utc": datetime.now(timezone.utc).isoformat(),
              "checks_passed": len(checks) - len(failures), "checks_total": len(checks),
              "checks": checks, "failures": failures, "detail": detail,
              "protocol_sha256": actual_protocol_hash,
              "code_sha256": {"pm.py": file_sha256(ROOT / "code/pm.py"),
                              "validate_pm.py": file_sha256(__file__)},
              "environment": {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__},
              "previous_attempts": [{"path": path.name, "sha256": file_sha256(path)} for path in attempts],
              "implementation_review_note": "After attempt 1 passed, additive finite/Hermitian input guards were introduced before production; no numerical algorithm, physical parameter, or frozen threshold was changed.",
              "wall_seconds": time.perf_counter() - started,
              "scope": "Small deterministic numerical validation; no observed-data fit, halo classification, or inflation calculation"}
    encoded = json.dumps(result, indent=2, ensure_ascii=False, allow_nan=False) + "\n"
    attempt_path = result_dir / f"pm_validation_attempt_{attempt:03d}.json"
    attempt_path.write_text(encoded)
    (result_dir / "pm_validation.json").write_text(encoded)
    lines = ["# 周期宇宙学PM的独立验证", "",
             f"状态：**{result['status']}，{result['checks_passed']}/{result['checks_total']}项通过**；第{attempt}次记录，历次失败不覆盖。", "",
             "冻结协议：`protocol.md`，哈希`" + actual_protocol_hash + "`。阈值未依据计算结果调整。", "",
             "此实现采用网格节点j/N、粒子初值(j+1/2)/N、CIC沉积/回插、谱Poisson及Nyquist安全梯度。", "",
             "| 检查 | 状态 | 结果 |", "| --- | --- | --- |"]
    for row in checks:
        description = f"{row['value']:.8g}" if "value" in row else "见JSON异常记录"
        lines.append(f"| {row['name']} | {'PASS' if row['passed'] else 'FAIL'} | {description} |")
    lines += ["", "## EdS增长的两种参照", "",
              "小幅平面模式初始密度幅度为2×10⁻⁶，a=.02→1，初速度采用连续增长支f=1。", "",
              "有限CIC网格的独立解析力响应为F=sin(kh)/(kh)，r±=(−1±sqrt(1+24F))/4。因为初速度不是网格自身的纯增长支，解析参照保留增长与衰减两支。", "",
              "时间收敛比较该网格解析增长；同时报告对连续EdS增长a/ai的空间偏差，不能把这两个误差混为一项。", "",
              "| 网格 | 步数 | 增长 | 对网格解析相对差 | 对连续EdS相对差 |", "| --- | --- | --- | --- | --- |"]
    for n, block in detail.get("plane_mode", {}).items():
        for row in block["runs"]:
            lines.append(f"| {n}³ | {row['steps']} | {row['growth']:.10g} | {row['relative_error_discrete_linear_reference']:.6g} | {row['relative_error_continuum_EdS']:.6g} |")
    lines += ["", "## 失败与修正记录", ""]
    if failures:
        lines += [f"- {row['name']}：{json.dumps(row, ensure_ascii=False)}" for row in failures]
    else:
        lines.append("本次无未通过检查；既往尝试及其哈希保存在JSON中。")
    lines += ["", "PM计算验证的是引入标准牛顿引力后的无碰撞物质接口，并未从黎曼结构导出引力。背景膨胀与有限网格平滑都保留；没有要求宇宙学物理能量在膨胀中恒定，也不因KDK名称宣称整个CIC粒子力已严格辛。", ""]
    (report_dir / "pm_validation_cn.md").write_text("\n".join(lines))
    print(json.dumps({"status": result["status"], "attempt": attempt,
                      "checks_passed": result["checks_passed"], "checks_total": result["checks_total"],
                      "wall_seconds": result["wall_seconds"], "failures": failures}, ensure_ascii=False, indent=2), flush=True)
    return 0 if not failures else 1


if __name__ == "__main__":
    raise SystemExit(main())
