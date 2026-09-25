#!/usr/bin/env python3
"""Independent direct high-precision determinants for the frozen SUGRA trial.

No cosmological ODE is integrated here.  The charged spectrum is reconstructed
from the defining functions, rather than the production small-splitting series.
Decimal background samples are fixed inputs, not high-precision observations.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
from datetime import datetime, timezone
from pathlib import Path
import platform

import mpmath as mp

ROOT = Path(__file__).resolve().parents[1]
REPO = ROOT.parents[1]
MODELS = ("global", "canonical", "shift")
INDICES = (0, 1, 256, 512, 768, 1024)


def number(x):
    return mp.mpf(str(x))


def decimal(x):
    return mp.nstr(x, 110, min_fixed=-5, max_fixed=8)


def relative(a, b):
    return abs(a-b)/max(abs(a), abs(b), mp.mpf("1e-1000"))


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def elementary(model, chi, inputs):
    """Independent physical mass and tree functions, all in eV units."""
    chi = number(chi)
    ci = number(inputs["chi_i"])
    eps = number(inputs["epsilon"])
    fsq = number(inputs["F_squared_eV2"])
    mpl2 = number(inputs["Mpl_eV"])**2
    mi = number(inputs.get("holomorphic_mass_i_eV", inputs.get("protected_mass_eV", inputs.get("mass_eV", "1e11"))))
    unit = number(inputs["potential_unit_eV4"])
    u = chi-ci
    g = (1+(ci/chi)**2)/2
    gp = -ci**2/chi**3
    mass_hol = mi*mp.sqrt(g)
    r = gp/(2*g)
    U = unit*number(inputs["W_i"])*mp.exp(-2*u)
    if model == "global":
        k = mp.mpf(0)
        s = mass_hol**2
        d = mp.mpf(0)
        B = mp.sqrt(2*U/fsq)*mass_hol*r
        V = U
    elif model == "canonical":
        k = eps*chi**2/2
        a = 1-eps*chi/2
        s = mp.exp(k)*mass_hol**2
        d = mp.exp(k)*U*(a*a-eps)/mpl2
        B = mp.exp(k)*mp.sqrt(2*U/fsq)*mass_hol*(a*(r+eps*chi/2)+eps/2)
        V = mp.exp(k)*U*(a*a-3*eps/2)
    elif model == "shift":
        k = mp.mpf(0)
        s = mass_hol**2
        d = U*(1-eps)/mpl2
        B = mp.sqrt(2*U/fsq)*mass_hol*(r+eps/2)
        V = (1-3*eps/2)*U
    else:
        raise ValueError(model)
    return dict(chi=chi, u=u, g=g, r=r, k=k, Mhol=mass_hol,
                physical_mass_eV=mp.sqrt(s), s=s, d=d, B=B,
                plus=s+d+abs(B), minus=s+d-abs(B), h2=B*B, U=U, V=V)


def raw_potential(model, chi, inputs):
    p = elementary(model, chi, inputs)
    mu2 = elementary(model, number(inputs["chi_i"]), inputs)["s"]
    def determinant(x):
        return x*x*(mp.log(x/mu2)-mp.mpf("1.5"))
    return (determinant(p["plus"])+determinant(p["minus"])
            -2*determinant(p["s"]))/ (32*mp.pi**2)


def evaluate(model, chi, today, inputs):
    chi = number(chi)
    p = elementary(model, chi, inputs)
    # Differentiate the full definitions. d(chi) is never held constant.
    for key in ("s", "d", "h2", "V"):
        p[key+"p"] = mp.diff(lambda c: elementary(model, c, inputs)[key], chi)
    p["Bp"] = mp.diff(lambda c: elementary(model, c, inputs)["B"], chi)
    p["U1"] = raw_potential(model, chi, inputs)
    p["U1p"] = mp.diff(lambda c: raw_potential(model, c, inputs), chi)
    p["mu2"] = elementary(model, number(inputs["chi_i"]), inputs)["s"]
    # An independent chain-rule derivative of the three unsimplified terms.
    def determinant_prime(x):
        return 2*x*(mp.log(x/p["mu2"])-1)
    hprime = mp.sign(p["B"])*p["Bp"]
    p["U1p_chain"] = (
        determinant_prime(p["plus"])*(p["sp"]+p["dp"]+hprime)
        +determinant_prime(p["minus"])*(p["sp"]+p["dp"]-hprime)
        -2*determinant_prime(p["s"])*p["sp"])/(32*mp.pi**2)
    p["loop_force_ratio"] = abs(p["U1p"]/p["Vp"])
    p["split_ratio"] = abs(p["B"])/p["s"]
    p["gravity_ratio"] = abs(p["d"])/p["s"]
    now = elementary(model, number(today), inputs)
    weighted = (mp.log(p["plus"]/now["plus"])
                +mp.log(p["minus"]/now["minus"])
                +4*mp.log(p["s"]/now["s"]))/3
    c = number(inputs["alpha_reference"])*weighted/(4*mp.pi)
    p["delta_alpha"] = c/(1-c)
    # The physical mass trend is also an independent threshold-sign check.
    p["dlog_s"] = p["sp"]/p["s"]
    return p


def load_primary():
    path = ROOT/"code/run_sugra.py"
    spec = importlib.util.spec_from_file_location("sugra_primary_precision", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def trajectory_rows(path):
    with path.open() as handle:
        rows = list(csv.DictReader(handle))
    grouped = {}
    for model in MODELS:
        group = [r for r in rows if r.get("model", r.get("case")) == model
                 and float(r.get("y_initial", r.get("y_i", r.get("initial_y", "0")))) == 0]
        if len(group) != 1025:
            raise ValueError(f"Expected 1025 axis rows for {model}, got {len(group)}")
        grouped[model] = group
    return grouped


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dps", type=int, default=320)
    parser.add_argument("--confirm-dps", type=int, default=400)
    args = parser.parse_args()
    if args.dps < 320 or args.confirm_dps <= args.dps:
        raise ValueError("Require at least 320 digits and a higher-precision repeat")
    (ROOT/"results").mkdir(parents=True, exist_ok=True)
    (ROOT/"reports").mkdir(parents=True, exist_ok=True)
    inputs_path = ROOT/"results/inputs.json"
    trajectory_path = ROOT/"results/axis_trajectories.csv"
    if not trajectory_path.exists():
        trajectory_path = ROOT/"results/trajectories.csv"
    inputs = json.loads(inputs_path.read_text(), parse_float=str)
    float_inputs = json.loads(inputs_path.read_text())
    primary = load_primary()
    rows = trajectory_rows(trajectory_path)
    mp.mp.dps = args.confirm_dps
    records = []
    comparison_max = {}
    repeat_max = mp.mpf(0)
    chain_max = mp.mpf(0)
    threshold_max = mp.mpf(0)
    checks = []
    compared = ("s", "sp", "d", "dp", "h2", "h2p", "U1", "U1p", "loop_force_ratio")
    for model in MODELS:
        today = rows[model][-1]["chi"]
        for index in INDICES:
            row = rows[model][index]
            chi = row["chi"]
            values = []
            for digits in (args.dps, args.confirm_dps):
                with mp.workdps(digits):
                    values.append(evaluate(model, chi, today, inputs))
            reference = values[1]
            repeat_errors = {key: relative(values[0][key], reference[key]) for key in compared}
            repeat_max = max(repeat_max, *repeat_errors.values())
            chain_error = relative(reference["U1p"], reference["U1p_chain"])
            chain_max = max(chain_max, chain_error)
            got = primary.particle_budget(model, float(chi), float_inputs)
            errors = {}
            for key in compared:
                if key not in got:
                    raise KeyError(f"Production particle_budget missing {key}; actual keys: {list(got)}")
                error = relative(number(got[key]), reference[key])
                errors[key] = error
                comparison_max[key] = max(comparison_max.get(key, mp.mpf(0)), error)
            threshold_error = relative(number(row["delta_alpha"]), reference["delta_alpha"])
            threshold_max = max(threshold_max, threshold_error)
            records.append(dict(model=model, index=index, N=row["N"], chi=chi,
                values={k: decimal(v) for k, v in reference.items()},
                precision_errors={k: decimal(v) for k, v in repeat_errors.items()},
                primary_errors={k: decimal(v) for k, v in errors.items()},
                chain_derivative_relative_error=decimal(chain_error),
                threshold_relative_error=decimal(threshold_error),
                production_delta_alpha=row["delta_alpha"]))
            checks.append(dict(name=f"{model}_{index}_spectrum_positive", passed=reference["minus"]>0))
            checks.append(dict(name=f"{model}_{index}_small_splitting", passed=
                max(reference["split_ratio"], reference["gravity_ratio"])<mp.mpf(".01")))
    checks.extend([
        dict(name="320_400_digit_direct_supertrace_agreement", passed=repeat_max<mp.mpf("1e-90")),
        dict(name="mp_diff_vs_full_chain_derivative", passed=chain_max<mp.mpf("1e-90")),
        dict(name="primary_series_vs_direct_determinants", passed=all(v<mp.mpf("1e-8") for v in comparison_max.values())),
        dict(name="physical_mass_weighted_threshold", passed=threshold_max<mp.mpf("1e-8")),
    ])
    for model in MODELS:
        selected = [r for r in records if r["model"] == model]
        trends = [number(r["values"]["dlog_s"]) for r in selected]
        expected_sign = 1 if model == "canonical" else -1
        checks.append(dict(name=f"{model}_physical_mass_trend_sign", passed=all(mp.sign(x)==expected_sign for x in trends)))
    files = [Path(__file__), ROOT/"protocol.md", inputs_path, trajectory_path, ROOT/"code/run_sugra.py"]
    report = dict(created_utc=datetime.now(timezone.utc).isoformat(),
        scope="Independent charged-particle spectrum, flat-space one-loop budget, and threshold on saved tree backgrounds. No quantum-feedback integration, full SUGRA quantum calculation, or observational test.",
        arithmetic_digits=[args.dps,args.confirm_dps],
        input_precision_note="Samples and parameter decimals are fixed inputs; extra arithmetic digits do not improve physical precision.",
        renormalization="One fixed mu^2=s(chi_i) per model. No field-dependent renormalization scale or potential subtraction.",
        python=platform.python_version(), mpmath=mp.__version__,
        source_sha256={str(p.relative_to(REPO)):sha(p) for p in files},
        sample_count=len(records), threshold_count=len(records),
        max_precision_relative_error=decimal(repeat_max),
        max_chain_derivative_relative_error=decimal(chain_max),
        max_primary_relative_errors={k:decimal(v) for k,v in comparison_max.items()},
        max_threshold_relative_error=decimal(threshold_max),
        checks=checks, all_checks_pass=all(c["passed"] for c in checks), records=records)
    out = ROOT/"results/particle_precision.json"
    out.write_text(json.dumps(report,ensure_ascii=False,indent=2,allow_nan=False)+"\n")
    brief = {k:report[k] for k in ("sample_count","max_precision_relative_error", "max_chain_derivative_relative_error", "max_primary_relative_errors", "max_threshold_relative_error","all_checks_pass")}
    print(json.dumps(brief,indent=2))
    text = ["# 带电粒子谱与一圈预算独立高精度复核", "", "2026-09-25。材料仅为冻结协议、程序和已存树级背景；没有新增观测。", "", "使用独立质量函数重建两个复标量和一个 Dirac 场的三个行列式，在 320 位及 400 位直接相减，并直接对完整 χ 依赖求导。引力共同分裂 d(χ) 的导数包含在内；μ² 固定为每个模型初始物理质量平方。", "", f"共 {len(records)} 个预定点，全部检查通过：{report['all_checks_pass']}。详细数据和 SHA-256 在 results/particle_precision.json。", "", "| 核验 | 最大相对误差 |", "|---|---:|", f"| 320 / 400 位重复 | {decimal(repeat_max)} |", f"| 直接微分 / 完整链式微分 | {decimal(chain_max)} |", f"| α 的物理质量阈值 | {decimal(threshold_max)} |"]
    text += [f"| 主程序 {key} | {decimal(value)} |" for key,value in comparison_max.items()]
    text += ["", "三个模型的抽样标量谱均为正，分裂在预定小分裂域内。canonical 的物理质量随 χ 增长；global 和 shift 的物理质量随 χ 下降。因而不能把仅使用全纯质量得到的阈值方向直接移植到 canonical 模型。", "", "执行命令：`timeout 180s python experiments/sugra_shift_v01/code/check_particle_precision.py`。首次执行的全部数值检查通过，但最后写 Markdown 时 reports/ 目录不存在，进程返回非零；原始日志保留为 results/particle_precision_attempt_001.log。增加输出目录初始化后重新执行，未改变物理公式、参数、样本或验收门槛。", "", "本复核只检查已定义理论的数值实现。主背景使用树级势，一圈势力是在轨迹上重新估算；这些检查不表示已经积分全量子反馈，不包括曲率有效作用量、引力子 / gravitino 循环、可见部门或观测约束。", ""]
    (ROOT/"reports/particle_precision_cn.md").write_text("\n".join(text))
    if not report["all_checks_pass"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
