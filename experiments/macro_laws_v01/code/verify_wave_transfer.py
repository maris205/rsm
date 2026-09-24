#!/usr/bin/env python3
"""Small, deterministic verification of a separately defined 3D wave lattice.

This implementation imports no CMB project code and loads no observations.
Read ../protocol_wave.md for the pre-recorded parameters and thresholds.
The model is externally scheduled, unlike RSM's autonomous clock coupling.
"""

from __future__ import annotations

import hashlib
import json
import platform
from datetime import datetime, timezone
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
STEPS = 45
OFFSET = 10.0
A2 = 0.45
SEED = 2026092404
SOURCE_COMMIT = "f42db51748003fc83be064dade9bb52de5eb8319"
SOURCE_HASHES = {
    "code/lattice_transfer.py": "6ed9effd4f9f734715907d4ea3193f31e8169542a378c4fc7a7c52adb570c5b4",
    "paper/main.tex": "e513fb52ca871f01b932f2504d4a305123f1d2e507b192c79dc031ab75fb2dcd",
    "reports/theory_derivation.md": "cd2ea373103c79c8d1dfd6a066c80c85274c2d68e4922517b2e85132d28b2cd3",
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def lattice_operator(q: np.ndarray) -> np.ndarray:
    """Normalized six-neighbour operator on the first three periodic axes."""
    result = -6.0 * q.copy()
    for axis in range(3):
        result += np.roll(q, 1, axis=axis)
        result += np.roll(q, -1, axis=axis)
    return result / 6.0


def evolve(q0: np.ndarray, coupling: np.ndarray, means: bool = False):
    previous, current = q0.copy(), q0.copy()
    mean_values = [float(np.mean(current))] if means else []
    for c2 in coupling:
        following = 2.0 * current - previous + c2 * lattice_operator(current)
        previous, current = current, following
        if means:
            mean_values.append(float(np.mean(current)))
    return current, np.asarray(mean_values)


def mode_transfer(lam: np.ndarray, coupling: np.ndarray) -> np.ndarray:
    previous = np.ones_like(lam, dtype=float)
    current = previous.copy()
    for c2 in coupling:
        previous, current = current, (2.0 - c2 * lam) * current - previous
    return current


def make_schedules():
    indices = np.arange(1, STEPS + 1, dtype=float)
    logs = np.log(indices + OFFSET)
    target = np.sqrt(A2) * np.sum(1.0 / logs)
    a0 = float((target / STEPS) ** 2)
    return indices, np.array([a0, A2]), np.stack([
        np.full(STEPS, a0), A2 / logs**2,
    ])


def lambda_stable(k: np.ndarray) -> np.ndarray:
    return (2.0 / 3.0) * np.sum(np.sin(k / 2.0)**2, axis=-1)


def relative_norm(x: np.ndarray, y: np.ndarray) -> float:
    return float(np.linalg.norm(x - y) / np.linalg.norm(y))


def main() -> int:
    results_dir = ROOT / "results"
    reports_dir = ROOT / "reports"
    results_dir.mkdir(parents=True, exist_ok=True)
    reports_dir.mkdir(parents=True, exist_ok=True)
    checks = []

    def check(name, passed, detail):
        checks.append({"name": name, "passed": bool(passed), "detail": detail})

    def upper(name, observed, threshold, **extra):
        check(name, np.isfinite(observed) and observed < threshold,
              {"observed": float(observed), "strict_upper_threshold": threshold, **extra})

    indices, amplitudes, schedules = make_schedules()
    phases = np.sum(np.sqrt(schedules), axis=1) / np.sqrt(6.0)
    upper("matched_long_wave_distance", abs(phases[0] - phases[1]), 1e-13,
          distances_lattice_sites=phases.tolist())

    # Coordinate-space and FFT routes share only their fixed model definition.
    size = 12
    rng = np.random.default_rng(SEED)
    initial = rng.standard_normal((size, size, size))
    k = 2.0 * np.pi * np.fft.fftfreq(size)
    wave_vectors = np.stack(np.meshgrid(k, k, k, indexing="ij"), axis=-1)
    eigenvalue_grid = 1.0 - np.sum(np.cos(wave_vectors), axis=-1) / 3.0
    initial_fft = np.fft.fftn(initial)
    coordinate_errors = []
    mean_errors = []
    for p, schedule in zip((0, 2), schedules):
        direct, means = evolve(initial, schedule, means=True)
        fft_output = np.fft.ifftn(initial_fft * mode_transfer(eigenvalue_grid, schedule))
        error = relative_norm(fft_output.real, direct)
        coordinate_errors.append(error)
        mean_error = float(np.max(np.abs(means - means[0])))
        mean_errors.append(mean_error)
        upper(f"coordinate_fft_relative_l2:p={p}", error, 1e-11)
        upper(f"linear_mean_conservation:p={p}", mean_error, 1e-13)

    coordinates = np.stack(np.meshgrid(*([np.arange(size)] * 3), indexing="ij"), axis=-1)
    mode_indices = [(0, 0, 0), (1, 0, 0), (1, 2, 3), (6, 6, 6)]
    mode_operator_errors = []
    for index in mode_indices:
        mode_k = 2.0 * np.pi * np.asarray(index) / size
        mode = np.exp(1j * np.einsum("...i,i->...", coordinates, mode_k))
        lam = float(1.0 - np.mean(np.cos(mode_k)))
        mode_operator_errors.append(float(np.max(np.abs(lattice_operator(mode) + lam * mode))))
    upper("coordinate_fourier_eigenvalue", max(mode_operator_errors), 1e-13,
          mode_indices=mode_indices, errors=mode_operator_errors)

    # Exact closed form for the constant schedule, with the stated backward start.
    lambda_axis = np.linspace(0.0, 2.0, 4097)
    transfers = np.stack([mode_transfer(lambda_axis, row) for row in schedules])
    omega = np.arccos(1.0 - amplitudes[0] * lambda_axis / 2.0)
    closed = np.cos((STEPS + 0.5) * omega) / np.cos(omega / 2.0)
    closed[0] = 1.0
    upper("constant_transfer_closed_form", np.max(np.abs(transfers[0] - closed)), 2e-11)

    J = np.array([[0.0, 1.0], [-1.0, 0.0]])
    single_symplectic = []
    product_symplectic = []
    product_determinant = []
    for schedule in schedules:
        for lam in (0.0, 0.2, 1.0, 2.0):
            product = np.eye(2)
            for c2 in schedule:
                matrix = np.array([[2.0 - c2 * lam, -1.0], [1.0, 0.0]])
                single_symplectic.append(float(np.max(np.abs(matrix.T @ J @ matrix - J))))
                product = matrix @ product
            product_symplectic.append(float(np.max(np.abs(product.T @ J @ product - J))))
            product_determinant.append(float(abs(np.linalg.det(product) - 1.0)))
    upper("single_mode_symplectic", max(single_symplectic), 1e-14)
    upper("45_step_mode_symplectic", max(product_symplectic), 1e-10)
    upper("45_step_mode_determinant", max(product_determinant), 1e-10)
    check("selected_frozen_oscillatory_upper_bound", np.max(schedules) * 2.0 < 4.0,
          {"max_c2_lambda": float(np.max(schedules) * 2.0), "strict_upper_bound": 4.0,
           "qualification": "The zero mode is neutral; this is not a proof for arbitrary time-dependent products."})

    # The cubic fourth-order term is direction dependent, even at equal |k|.
    direction_names = np.array(["axis", "face_diagonal", "body_diagonal"])
    directions = np.array([[1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [1.0, 1.0, 1.0]])
    directions /= np.linalg.norm(directions, axis=1, keepdims=True)
    quartic = -np.sum(directions**4, axis=1) / 72.0
    small_r = 0.01
    small_lambda = lambda_stable(small_r * directions)
    leading_error = np.abs(small_lambda / small_r**2 - 1.0 / 6.0) / (1.0 / 6.0)
    quartic_measured = (small_lambda - small_r**2 / 6.0) / small_r**4
    upper("long_wave_coefficient_one_sixth", np.max(leading_error), 1e-5,
          r=small_r, errors_by_direction=leading_error.tolist())
    upper("quartic_direction_coefficient", np.max(np.abs(quartic_measured - quartic)), 1e-7,
          expected=quartic.tolist(), measured=quartic_measured.tolist(), r=small_r)
    test_r = np.array([0.1, 0.05])
    test_lambda = lambda_stable(directions[:, None, :] * test_r[None, :, None])
    leading = test_r[None, :]**2 / 6.0
    corrected = leading + quartic[:, None] * test_r[None, :]**4
    errors2 = np.abs(test_lambda - leading) / test_lambda
    errors4 = np.abs(test_lambda - corrected) / test_lambda
    orders2 = np.log2(errors2[:, 0] / errors2[:, 1])
    orders4 = np.log2(errors4[:, 0] / errors4[:, 1])
    check("relative_long_wave_error_order_2", np.all((orders2 >= 1.95) & (orders2 <= 2.05)),
          {"orders": orders2.tolist(), "inclusive_interval": [1.95, 2.05], "radii": test_r.tolist()})
    check("relative_corrected_error_order_4", np.all((orders4 >= 3.95) & (orders4 <= 4.05)),
          {"orders": orders4.tolist(), "inclusive_interval": [3.95, 4.05], "radii": test_r.tolist()})
    radius_axis = np.geomspace(0.01, np.pi, 241)
    curve_k = directions[:, None, :] * radius_axis[None, :, None]
    lambda_curves = lambda_stable(curve_k)
    cosine_curves = 1.0 - np.sum(np.cos(curve_k), axis=-1) / 3.0
    upper("sine_cosine_eigenvalue_identity", np.max(np.abs(lambda_curves - cosine_curves)), 1e-14)
    leading_curves = np.broadcast_to(radius_axis[None, :]**2 / 6.0, lambda_curves.shape).copy()
    corrected_curves = leading_curves + quartic[:, None] * radius_axis[None, :]**4

    # Full covariance: evolve every column of a deterministic covariance square root.
    # The direct route never draws a Monte Carlo ensemble or averages radial bins.
    cov_size = 4
    dimension = cov_size**3
    identity = np.eye(dimension)
    basis_volumes = identity.reshape((cov_size, cov_size, cov_size, dimension))
    Lmatrix = lattice_operator(basis_volumes).reshape((dimension, dimension))
    lam, U = np.linalg.eigh(-Lmatrix)
    upper("dense_laplacian_symmetry", np.max(np.abs(Lmatrix - Lmatrix.T)), 1e-13)
    upper("dense_eigendecomposition", np.max(np.abs((-Lmatrix) @ U - U * lam)), 1e-12)
    cov_rng = np.random.default_rng(SEED + 1)
    B = cov_rng.standard_normal((dimension, dimension))
    root_cov = np.concatenate([B / np.sqrt(dimension), np.sqrt(0.2) * identity], axis=1)
    K0 = root_cov @ root_cov.T
    dense_outputs = []
    dense_transfers = []
    covariance_errors = []
    for p, schedule in zip((0, 2), schedules):
        direct_basis, _ = evolve(basis_volumes, schedule)
        Fdirect = direct_basis.reshape((dimension, dimension))
        transfer = mode_transfer(lam, schedule)
        Fmodal = (U * transfer) @ U.T
        root_volume = root_cov.reshape((cov_size, cov_size, cov_size, 2 * dimension))
        evolved_root, _ = evolve(root_volume, schedule)
        evolved_root = evolved_root.reshape((dimension, 2 * dimension))
        Kdirect = evolved_root @ evolved_root.T
        Kformula = Fmodal @ K0 @ Fmodal.T
        cov_error = relative_norm(Kdirect, Kformula)
        covariance_errors.append(cov_error)
        upper(f"coordinate_modal_operator:p={p}", relative_norm(Fdirect, Fmodal), 1e-10)
        upper(f"complete_covariance_propagation:p={p}", cov_error, 1e-10)
        dense_outputs.append(Fdirect)
        dense_transfers.append(transfer)
    dense_transfers = np.asarray(dense_transfers)

    retained = np.all(np.abs(dense_transfers) > 0.1, axis=0)
    check("retained_mode_subspace_nonempty", np.any(retained),
          {"retained_modes": int(np.sum(retained)), "total_modes": dimension,
           "rule": "Both schedules have |T_N| > 0.1 before division."})
    Ur = U[:, retained]
    target_power = 0.1 + 1.0 / (1.0 + lam[retained])
    target_cov = (Ur * target_power) @ Ur.T
    degenerate_outputs = []
    prior_statistics = []
    for p, Fdirect, transfer in zip((0, 2), dense_outputs, dense_transfers):
        input_power = target_power / transfer[retained]**2
        input_cov = (Ur * input_power) @ Ur.T
        final_cov = Fdirect @ input_cov @ Fdirect.T
        degenerate_outputs.append(final_cov)
        upper(f"prior_absorption_target_covariance:p={p}", relative_norm(final_cov, target_cov), 1e-10)
        prior_statistics.append({"p": p, "minimum_retained_abs_T": float(np.min(np.abs(transfer[retained]))),
                                 "maximum_initial_mode_power": float(np.max(input_power)),
                                 "minimum_initial_mode_power": float(np.min(input_power))})
    upper("prior_absorption_equal_final_covariances", relative_norm(*degenerate_outputs), 1e-10)

    # Finite sums are a diagnostic of the asymptotic formula, not its proof.
    sum_sizes = np.array([45, 100, 1000, 10000, 100000, 1000000], dtype=int)
    sum_exponents = np.array([0.0, 1.0, 2.0, 3.0])
    all_n = np.arange(1, sum_sizes[-1] + 1, dtype=float)
    log_all = np.log(all_n + OFFSET)
    distances = []
    asymptotic_ratios = []
    for p in sum_exponents:
        partial = np.cumsum(log_all**(-p / 2.0))[sum_sizes - 1]
        actual = np.sqrt(A2 / 6.0) * partial
        leading_asymptotic = np.sqrt(A2 / 6.0) * sum_sizes / np.log(sum_sizes)**(p / 2.0)
        distances.append(actual)
        asymptotic_ratios.append(actual / leading_asymptotic)
    distances = np.asarray(distances)
    asymptotic_ratios = np.asarray(asymptotic_ratios)
    check("finite_propagation_sum_diagnostics", np.all(np.isfinite(asymptotic_ratios)) and np.all(distances > 0),
          {"qualification": "No convergence rate, finite horizon, or asymptotic theorem is proved by these finite sums."})
    upper("constant_schedule_sum_ratio", np.max(np.abs(asymptotic_ratios[0] - 1.0)), 2e-14)

    metadata = {
        "model": "Externally scheduled normalized three-dimensional linear wave lattice",
        "units": {"n": "dimensionless update index", "lambda": "dimensionless eigenvalue of -L",
                  "k_radius": "radians per lattice spacing", "q": "arbitrary scalar field units",
                  "mode_power": "dimensionless |T_N|^2", "propagation_distance": "lattice sites"},
        "steps": STEPS, "coefficient_indices": "1,...,45", "initial_pair": "q[-1] = q[0]",
        "offset": OFFSET, "exponents": [0, 2], "amplitudes": amplitudes.tolist(),
        "amplitude_matching": "Equal sum sqrt(c_n^2), which only matches a long-wave phase proxy",
        "coordinate_grid_size": size, "covariance_grid_size": cov_size,
        "seed_coordinate": SEED, "seed_covariance": SEED + 1,
        "source_commit": SOURCE_COMMIT, "source_sha256_verified_before_implementation": SOURCE_HASHES,
        "protocol_sha256": sha256(ROOT / "protocol_wave.md"),
        "script_sha256": sha256(Path(__file__)),
        "not_included": ["observations", "sky fitting", "autonomous clock", "drag", "clipping", "nonlinear force"],
    }
    curve_path = results_dir / "wave_transfer_curves.npz"
    np.savez_compressed(
        curve_path, metadata_json=np.array(json.dumps(metadata, ensure_ascii=False, sort_keys=True)),
        coefficient_indices=indices, exponents=np.array([0, 2]), amplitudes=amplitudes,
        coupling_schedule=schedules, propagation_distance_sites=phases,
        lambda_axis=lambda_axis, mode_transfer=transfers, mode_power=transfers**2,
        lambda_grid=lambda_axis, transfer_squared_p0=transfers[0]**2, transfer_squared_p2=transfers[1]**2,
        constant_closed_transfer=closed, constant_closed_power=closed**2,
        k_radius=radius_axis, direction_names=direction_names, directions=directions,
        lambda_exact=lambda_curves, lambda_leading=leading_curves, lambda_corrected=corrected_curves,
        k_grid=radius_axis, lambda_axis_exact=lambda_curves[0], lambda_diagonal_exact=lambda_curves[2],
        lambda_second=leading_curves[0], lambda_axis_fourth=corrected_curves[0],
        lambda_diagonal_fourth=corrected_curves[2],
        relative_error_leading=np.abs(lambda_curves - leading_curves) / lambda_curves,
        relative_error_corrected=np.abs(lambda_curves - corrected_curves) / lambda_curves,
        sum_sizes=sum_sizes, sum_exponents=sum_exponents, cumulative_distance_sites=distances,
        cumulative_asymptotic_ratio=asymptotic_ratios,
    )
    passed = sum(row["passed"] for row in checks)
    result = {
        "status": "PASS" if passed == len(checks) else "FAIL",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "checks_passed": passed, "checks_total": len(checks), "checks": checks,
        "metadata": metadata,
        "runtime": {"python": platform.python_version(), "numpy": np.__version__},
        "summary": {
            "coordinate_fft_relative_l2": coordinate_errors, "mean_max_abs_errors": mean_errors,
            "covariance_relative_errors": covariance_errors,
            "matched_distance_lattice_sites": phases.tolist(),
            "long_wave_relative_error_orders": orders2.tolist(),
            "corrected_relative_error_orders": orders4.tolist(),
            "prior_absorption": {"retained_modes": int(np.sum(retained)), "excluded_modes": int(np.sum(~retained)),
                                 "total_modes": dimension, "schedule_statistics": prior_statistics},
        },
        "finite_sum_diagnostic": {"p": sum_exponents.tolist(), "N": sum_sizes.tolist(),
                                  "A_all_p": A2, "offset": OFFSET, "ratios": asymptotic_ratios.tolist(),
                                  "not_a_proof": True},
        "artifacts": {"results/wave_transfer_curves.npz": sha256(curve_path)},
    }
    result_path = results_dir / "wave_transfer_checks.json"
    result_path.write_text(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False) + "\n")

    report = [
        "# 三维线性格点的微观—宏观连接：小型确定性验证", "",
        f"运行状态：**{result['status']}，{passed}/{len(checks)} 项通过**。协议先记录后运行；本地检查不构成外部注册或同行评审。", "",
        "## 模型与检验范围", "",
        "这是与主 RSM 自主时钟模型分开的三维线性候选。周期六邻居算子除以6，主递推的初态为 q[-1] = q[0]，执行45步；不加载天空数据、不作观测拟合、没有阻尼、非线性或裁剪。单位均为格点/迭代单位。", "",
        f"平方对数分支 A2={A2}、c0={OFFSET:g}；常系数分支 A0={amplitudes[0]:.12g} 匹配累计长波距离代理量，双方均为 {phases[0]:.12g} 个格点。该条件不匹配每个 Fourier 模态的精确相位。", "",
        f"坐标/FFT检查网格为{size}³；完整协方差检查网格为{cov_size}³，直接演化确定性协方差平方根的全部列，无大样本ensemble及径向平均。", "",
        "## 检查结果", "", "| 检查 | 状态 | 数值或说明 |", "| --- | --- | --- |",
    ]
    for row in checks:
        detail = row["detail"]
        if "observed" in detail:
            description = f"{detail['observed']:.6g} < {detail['strict_upper_threshold']:.3g}"
        elif "orders" in detail:
            description = ", ".join(f"{x:.6f}" for x in detail["orders"])
        elif "retained_modes" in detail:
            description = f"保留 {detail['retained_modes']}/{detail['total_modes']} 个模式"
        elif "max_c2_lambda" in detail:
            description = f"max(c²λ) = {detail['max_c2_lambda']:.6g} < 4；仅冻结系数条件"
        else:
            description = "有限且为正；只作有限求和诊断"
        report.append(f"| {row['name']} | {'PASS' if row['passed'] else 'FAIL'} | {description} |")
    report += ["", "## 初谱退化与有限求和", "",
               f"按事先固定的双方 |T|>0.1 规则保留 {int(np.sum(retained))}/{dimension} 个模式，排除 {int(np.sum(~retained))} 个；没有除以传递零点。各历史输入功率由同一目标功率除以该历史 T² 给出。保留子空间上可得到同一完整终态协方差；这显示初态自由度，不能解释成观测对任一历史的支持。", "",
               "有限和与渐近主项之比：", "", "| p | N=45 | N=100 | N=1,000 | N=10,000 | N=100,000 | N=1,000,000 |",
               "| --- | --- | --- | --- | --- | --- | --- |"]
    for p, row in zip(sum_exponents, asymptotic_ratios):
        report.append(f"| {p:g} | " + " | ".join(f"{x:.8f}" for x in row) + " |")
    report += ["", "此处所有p单独固定A=0.45，A在比值中抵消；不是45步匹配设计的新增模型比较。有限和不能证明渐近等价或给出有限声学视界。p=2的数学主项为N/ln N，仍发散。", "",
               "## 解释边界", "",
               "- 长波展开的1/6与四阶方向项来自三维归一化邻居几何；现有参数尚未定标为物理光速、格距或宇宙年龄。",
               "- 初始相等位置是零后向差分，不是连续零初速度的二阶起步；常系数闭式保留(N+1/2)相位及分母。",
               "- 辛性与固定模态振荡条件没有证明所有时变系统的稳定性或自主能量守恒。",
               "- 平均场保持是本线性周期模型的性质，尚未识别为物质质量守恒。",
               "- 功率是场的二阶统计。不同历史都可产生振荡特征，不能由峰形单独推出平方对数规律、原子能级或CMB光子能谱。",
               "- 协方差吸收允许每个保留模态自由调整初始方差；固定低维初谱、传递零点或额外共同观测可能限制此自由度。", "",
               "## 失败记录与来源", ""]
    failures = [row for row in checks if not row["passed"]]
    if failures:
        report += [f"- {row['name']}：{json.dumps(row['detail'], ensure_ascii=False)}" for row in failures]
    else:
        report.append("本次协议内检查无失败。全部数值、阈值与条件保存在JSON中；没有在运行后调整阈值。")
    report += ["", f"来源提交：`{SOURCE_COMMIT}`。实现前以 `git show` 核对三份来源SHA256；本脚本不导入相邻项目代码。", "",
               f"协议 SHA256：`{metadata['protocol_sha256']}`。", "",
               f"脚本 SHA256：`{metadata['script_sha256']}`。", "",
               f"曲线 NPZ SHA256：`{sha256(curve_path)}`。", "",
               "曲线文件同时保存轴、方向、单位、参数JSON；读取时无需允许pickle。完整协方差检查是有限矩阵恒等式核对，不是新观测或物理机制发现。", ""]
    (reports_dir / "wave_transfer_checks_cn.md").write_text("\n".join(report))
    print(json.dumps({"status": result["status"], "checks_passed": passed, "checks_total": len(checks),
                      "summary": result["summary"], "failures": failures}, ensure_ascii=False, indent=2))
    return 0 if not failures else 1


if __name__ == "__main__":
    raise SystemExit(main())
