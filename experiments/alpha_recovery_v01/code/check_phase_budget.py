#!/usr/bin/env python3
"""Analytic adiabatic and frozen-integrator budget; no observational fitting.

Run from any directory. Only the two sibling result/report files are written.
The redshift grid is a deterministic evaluation grid, not an observation set.
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
YEAR_SECONDS = 365.25 * 86400.0
MPC_KM = 3.0856775814913673e19
H0 = 67.4
OMEGA_M = 0.315
TSTAR_SECONDS = 5.391247e-44
NU0_HZ = 5.0e14
OMEGA0_PER_SECOND = 2.0 * math.pi * NU0_HZ
Z_MAX = 4.2
D0_VALUES = (0.0, -2.5e-19, 2.5e-19, -1.25e-18, 1.25e-18)


def age_years(z):
    """Flat matter+Lambda age: inherited fixed background, no radiation term."""
    z = np.asarray(z, dtype=float)
    return 2.0 / (3.0 * H0 / MPC_KM * YEAR_SECONDS * np.sqrt(1.0 - OMEGA_M)) * np.arcsinh(
        np.sqrt((1.0 - OMEGA_M) / OMEGA_M) / (1.0 + z) ** 1.5
    )


T0_YEARS = float(age_years(0.0))
L0 = math.log(T0_YEARS * YEAR_SECONDS / TSTAR_SECONDS)


def transfer_years(t_years, s):
    """Stable T_s, including its exact zero at today."""
    ell = np.log(np.asarray(t_years) / T0_YEARS)
    return -T0_YEARS * L0 / s * np.expm1(-s * np.log1p(ell / L0))


def transfer_derivative(t_years, s):
    t_years = np.asarray(t_years)
    log_ratio = np.log(t_years / T0_YEARS)
    return (T0_YEARS / t_years) * (L0 / (L0 + log_ratio)) ** (s + 1)


def phase_bias(x):
    """2 asin(x/2)/x - 1, avoiding cancellation for small x."""
    if not 0.0 < x < 2.0:
        raise ValueError("Strict frozen stability requires 0 < x < 2.")
    if x < 0.01:
        return x*x/24.0 + 3.0*x**4/640.0 + 5.0*x**6/7168.0
    return 2.0 * math.asin(x / 2.0) / x - 1.0


def drift_factor_excess(x):
    """d ln(omega_num)/d ln(omega) - 1 at FIXED step size."""
    if x < 0.01:
        return x*x/12.0 + 11.0*x**4/720.0 + 191.0*x**6/60480.0
    return x / (2.0 * math.asin(x / 2.0) * math.sqrt(1.0 - x*x/4.0)) - 1.0


def budget():
    redshifts = np.linspace(0.0, Z_MAX, 4201)
    times = age_years(redshifts)
    cases = []
    checks = {}
    for s in (1, 2, 3):
        transfer = transfer_years(times, s)
        derivative = transfer_derivative(times, s)
        for d0 in D0_VALUES:
            delta = d0 * transfer
            a = 1.0 + delta
            for kappa in (-7.0, 7.0):
                # D0 is per year, while the optical angular frequency is per second.
                omega = OMEGA0_PER_SECOND * np.exp(kappa * np.log1p(delta))
                drift_per_year = d0 * derivative / a
                epsilon = np.abs(kappa * drift_per_year) / (YEAR_SECONDS * omega)
                maximum = int(np.argmax(epsilon))
                cases.append({
                    "s": s, "D0_per_year": d0, "kappa": kappa,
                    "delta_alpha_zmax": float(delta[-1]),
                    "A_min": float(a.min()), "A_max": float(a.max()),
                    "A_positive_over_interval_by_monotonicity": bool(a[0] > 0 and a[-1] > 0),
                    "max_adiabatic_epsilon_on_grid": float(epsilon[maximum]),
                    "z_of_grid_maximum": float(redshifts[maximum]),
                    "epsilon_today": float(epsilon[0]),
                    "max_abs_D_per_year": float(np.max(np.abs(drift_per_year))),
                })

    x_values = (0.8, 0.4, 0.2, 0.1, 0.05, 0.025, 0.0125, 0.00625)
    rows = []
    for i, x in enumerate(x_values):
        rows.append({
            "x_equals_dt_omega": x,
            "dt_seconds_at_nu0": x / OMEGA0_PER_SECOND,
            "fractional_frozen_frequency_bias": phase_bias(x),
            "fractional_drift_amplification_minus_one": drift_factor_excess(x),
            "bias_ratio_to_next_halving": phase_bias(x)/phase_bias(x/2.0),
        })

    target_bias = 1.0e-18
    low, high = 1.0e-12, 1.0e-5
    for _ in range(100):
        middle = (low + high) / 2.0
        if phase_bias(middle) > target_bias:
            high = middle
        else:
            low = middle
    x_target = (low + high) / 2.0
    age_seconds = T0_YEARS * YEAR_SECONDS

    # Independent finite differences/eigenvalues check the analytic implementation.
    derivative_errors = []
    for s in (1, 2, 3):
        for t in age_years(np.array([0.0, 0.5, 1.0, 2.0, Z_MAX])):
            dt = t * 1.0e-5
            finite_difference = (transfer_years(t + dt, s) - transfer_years(t - dt, s)) / (2.0 * dt)
            derivative_errors.append(abs(float(finite_difference/transfer_derivative(t, s)-1.0)))
    eigenphase_errors = []
    drift_errors = []
    for x in x_values:
        matrix = np.array([[1.0 - x*x, x], [-x, 1.0]])  # Omega=1, h=x
        theta_eigenvalue = max(np.angle(np.linalg.eigvals(matrix)))
        eigenphase_errors.append(abs(float(theta_eigenvalue - 2.0*math.asin(x/2.0))))
        step = 1.0e-5
        log_plus = math.log(2.0*math.asin(x*math.exp(step)/2.0))
        log_minus = math.log(2.0*math.asin(x*math.exp(-step)/2.0))
        finite_difference = (log_plus-log_minus)/(2.0*step)
        drift_errors.append(abs(finite_difference - 1.0 - drift_factor_excess(x)))
    checks = {
        "transfer_today_exact_zero": all(float(transfer_years(T0_YEARS, s)) == 0.0 for s in (1, 2, 3)),
        "transfer_derivative_max_relative_error": max(derivative_errors),
        "transfer_derivative_check_pass": max(derivative_errors) < 1.0e-8,
        "all_zero_injections_exactly_zero": all(c["delta_alpha_zmax"] == 0.0 and c["max_adiabatic_epsilon_on_grid"] == 0.0 for c in cases if c["D0_per_year"] == 0.0),
        "all_injections_A_positive": all(c["A_positive_over_interval_by_monotonicity"] for c in cases),
        "frozen_matrix_eigenphase_max_absolute_error": max(eigenphase_errors),
        "frozen_matrix_eigenphase_check_pass": max(eigenphase_errors) < 1.0e-12,
        "frozen_drift_derivative_max_absolute_error": max(drift_errors),
        "frozen_drift_derivative_check_pass": max(drift_errors) < 1.0e-8,
        "asymptotic_phase_bias_halving_ratio": rows[-1]["bias_ratio_to_next_halving"],
        "phase_bias_second_order_halving_check_pass": abs(rows[-1]["bias_ratio_to_next_halving"]-4.0) < 1.0e-4,
    }
    assert all(v for k, v in checks.items() if isinstance(v, bool)), checks
    return {
        "scope": "Deterministic analytic effective-probe budget; no observation fitting and no cosmic optical-cycle integration.",
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "parameters": {
            "H0_km_s_Mpc": H0, "Omega_m": OMEGA_M,
            "background": "fixed flat matter+Lambda, radiation omitted",
            "tstar_seconds_fixed_reference": TSTAR_SECONDS,
            "seconds_per_year": YEAR_SECONDS, "nu0_hz_illustrative": NU0_HZ,
            "kappa_illustrative_values": [-7.0, 7.0],
            "D0_per_year_injections": list(D0_VALUES), "s_values": [1, 2, 3],
            "z_min": 0.0, "z_max": Z_MAX, "grid_count": len(redshifts),
            "t0_years": T0_YEARS, "t_at_zmax_years": float(times[-1]), "L0": L0,
        },
        "cosmological_probe_cases": cases,
        "frozen_phase_convergence": rows,
        "illustrative_absolute_phase_bias_requirement": {
            "target_fractional_frequency_bias": target_bias,
            "x_required_max": x_target, "dt_required_seconds_at_nu0": x_target/OMEGA0_PER_SECOND,
            "steps_over_t0_at_that_dt": age_seconds*OMEGA0_PER_SECOND/x_target,
            "optical_cycles_over_t0_at_constant_nu0": age_seconds*NU0_HZ,
            "interpretation": "Absolute uncalibrated numerical frequency-bias illustration, not a requirement for differential clock data or an observed drift.",
        },
        "checks": checks,
    }


def write_report(result):
    p = result["parameters"]
    cases = result["cosmological_probe_cases"]
    extreme = max(cases, key=lambda case: case["max_adiabatic_epsilon_on_grid"])
    required = result["illustrative_absolute_phase_bias_requirement"]
    lines = [
        "# 有效探针的绝热与数值相位预算",
        "",
        "## 材料与范围",
        "",
        "本记录由本地确定性脚本生成。输入是已定义的最小有效模型、固定背景和人为注入幅度；没有增加观测点，没有重新拟合原子钟，也没有把经典探针当成真实原子微观推导。",
        "",
        "来源：[最小公设](../../../draft_v01/minimal_framework_v01.md)、[数学审查](../../../reports/minimal_dynamics_v01_review.md)、[已有钟约束与符号](../../../reports/alpha_interface_v01_review.md)。",
        "",
        "## 参数与公式",
        "",
        r"采用 $H_0=67.4\,\mathrm{km\,s^{-1}\,Mpc^{-1}}$、$\Omega_m=0.315$ 的平直物质加 $\Lambda$ 背景，略去辐射，固定 $t_*=5.391247\times10^{-44}\,\mathrm s$。这里只在 $0\le z\le4.2$ 使用它。参考尺度不解释为数值步长。",
        "",
        f"得到 $t_0={p['t0_years']/1e9:.6f}$ Gyr、$t(4.2)={p['t_at_zmax_years']/1e9:.6f}$ Gyr、$L_0={p['L0']:.6f}$。对数指数分别取 $s=1,2,3$。",
        "",
        r"注入 $D_0=0,\pm2.5\times10^{-19},\pm1.25\times10^{-18}\,\mathrm{yr}^{-1}$；后两种绝对幅度是既有压缩钟估计误差的 1 倍和 5 倍，用于测试。5 倍误差注入不是声称观测已允许或已发现该信号。取示意探针 $\nu_0=5\times10^{14}\,\mathrm{Hz}$ 和 $\kappa=\pm7$；它们不是新增测量，也不把单探针 $\kappa$ 冒充实际钟比敏感度。",
        "",
        r"$$A=1+D_0\mathcal T_s(t),\quad Q_s(t)=\mathcal T_s'(t)=\frac{t_0}{t}\left(\frac{L_0}{L}\right)^{s+1},\quad D(t)=\frac{D_0Q_s}{A}.$$",
        "",
        r"$$\omega=2\pi\nu_0 A^\kappa,\qquad \epsilon_{\rm ad}=\frac{|\dot\omega|}{\omega^2}=\frac{|\kappa D_0|Q_s}{(365.25\times86400)\,2\pi\nu_0 A^{\kappa+1}}.$$",
        "",
        "此处 $D_0$ 按年计，角频率按秒计，公式中的年秒换算不可遗漏。$Q_s>0$，所以 $A$ 随时间单调；两个端点为正即保证本区间内为正。",
        "",
        "## 慢变尺度",
        "",
        r"下表每行合并正负注入和两种 $\kappa$ 符号。最大绝热参数是在明确列出的 4201 点网格上求得；正性由端点与单调性保证。",
        "",
        r"| $s$ | $\lvert D_0\rvert$ / yr$^{-1}$ | $\lvert\delta_\alpha(4.2)\rvert$ / ppm | 区间最小 $A$ | 最大 $\epsilon_{\rm ad}$ |",
        "| --- | ---: | ---: | ---: | ---: |",
    ]
    for s in (1, 2, 3):
        for amplitude in (2.5e-19, 1.25e-18):
            selected = [c for c in cases if c["s"] == s and abs(c["D0_per_year"]) == amplitude]
            lines.append(f"| {s} | {amplitude:.3e} | {abs(selected[0]['delta_alpha_zmax'])*1e6:.6f} | {min(c['A_min'] for c in selected):.12f} | {max(c['max_adiabatic_epsilon_on_grid'] for c in selected):.3e} |")
    lines += [
        "",
        rf"零注入严格给出 $A=1$、$\epsilon_{{\rm ad}}=0$。全组最大值为 **{extreme['max_adiabatic_epsilon_on_grid']:.3e}**，出现在 $z={extreme['z_of_grid_maximum']:.1f}$。这些指定光学探针的驱动尺度远慢于振荡尺度，因此使用瞬时有效频率有很大的尺度分离。局部绝热指标不是对任意多模态系统、长期总相位或所有频率估计器的严格误差上界。",
        "",
        "## 冻结积分器的相位偏差与步长",
        "",
        r"对严格冻结的单模式定义 $x=h\Omega=\Delta t\,\omega$，稳定振荡域为 $0<x<2$。辛 Euler 的本征相位及其频率偏差是",
        "",
        r"$$\theta=2\arcsin(x/2),\qquad b_\nu(x)=\frac{\theta}{x}-1=\frac{x^2}{24}+O(x^4).$$",
        "",
        "在固定物理步长下，若冻结系数被缓慢改变，将瞬时数值频率当作读数还会改变漂移响应：",
        "",
        r"$$\frac{d\ln\omega_{\rm num}}{dt}=B(x)\frac{d\ln\omega}{dt},\qquad B(x)=\frac{x}{2\arcsin(x/2)\sqrt{1-x^2/4}}=1+\frac{x^2}{12}+O(x^4).$$",
        "",
        r"这是冻结色散关系对参数的导数；完整时变轨迹还含绝热与估计误差。常数 $\omega$ 和固定步长只产生固定频率偏置，不凭空产生时间漂移。实际频率比还可能消除共同偏置，其残余要按两个探针分别计算；自适应变步长则需另行分析。",
        "",
        r"| $x$ | $\Delta t$ / s（示意光频） | $b_\nu$ | $B-1$ | 再减半步长时偏差比 |",
        "| ---: | ---: | ---: | ---: | ---: |",
    ]
    for row in result["frozen_phase_convergence"]:
        lines.append(f"| {row['x_equals_dt_omega']:.5f} | {row['dt_seconds_at_nu0']:.3e} | {row['fractional_frozen_frequency_bias']:.3e} | {row['fractional_drift_amplification_minus_one']:.3e} | {row['bias_ratio_to_next_halving']:.6f} |")
    lines += [
        "",
        "冻结相位偏差随步长平方收敛；这不把辛 Euler 完整状态更新的一阶精度改成二阶，也不是时变乘积长期稳定性证明。接近 $x=2$ 时偏差及其导数增大，不能只满足稳定阈值就称为高精度。",
        "",
        rf"若示意性地要求未校准绝对数值频率偏置小于 $10^{{-18}}$，此积分器需 $x\lesssim{required['x_required_max']:.3e}$，即 $\Delta t\lesssim{required['dt_required_seconds_at_nu0']:.3e}$ s。按恒定示意频率跨越整个 $t_0$ 约有 {required['optical_cycles_over_t0_at_constant_nu0']:.3e} 个周期、该步长约需 {required['steps_over_t0_at_that_dt']:.3e} 步。这只是说明逐周期宇宙积分不合适；不是实际钟比实验必须满足的未校准绝对偏置要求。",
        "",
        "## 对本轮恢复实验的处理",
        "",
        r"观测模拟应直接计算 $\delta_\alpha=D_0\mathcal T_s$ 与 $d\ln R/dt=\Delta\kappa\,D(t)$，无需积分原子光学载波。积分器步长是计算选择，不是物理漂移参数；未来若从轨迹提取频率，必须单独做步长、绝热和读数误差收敛。当前预算既没有发现新的物理变化，也没有从绝热性证明逆对数平方优于其他指数。",
        "",
        "## 复现与核验",
        "",
        "```bash",
        "python riemann_model/experiments/alpha_recovery_v01/code/check_phase_budget.py",
        "```",
        "",
        "依赖 Python 与 NumPy。输出：[数值 JSON](../results/phase_budget.json)、本报告。脚本用稳定 `log1p/expm1` 计算今天附近的传递函数，用小量级数避免相位偏差的灾难性相消；源码散列记录在 JSON 中。",
        "",
        f"全部布尔核验通过：传递函数今天为零、解析导数与有限差分一致、零注入严格为零、各注入保持正性、冻结矩阵本征相位一致、冻结色散导数一致、相位步长比趋近 4。传递导数最大相对差为 {result['checks']['transfer_derivative_max_relative_error']:.3e}，矩阵本征相位最大绝对差为 {result['checks']['frozen_matrix_eigenphase_max_absolute_error']:.3e}。这些是实现与公式检查。没有执行有限区间非自治玩具轨道、能量积分或观测显著性检验。",
        "",
    ]
    return "\n".join(lines)


def main():
    result = budget()
    (ROOT / "results").mkdir(parents=True, exist_ok=True)
    (ROOT / "reports").mkdir(parents=True, exist_ok=True)
    (ROOT / "results" / "phase_budget.json").write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n")
    (ROOT / "reports" / "phase_budget_cn.md").write_text(write_report(result))
    print(json.dumps({"cases": len(result["cosmological_probe_cases"]), "checks": result["checks"]}, indent=2))


if __name__ == "__main__":
    main()
