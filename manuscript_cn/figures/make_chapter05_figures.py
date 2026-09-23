#!/usr/bin/env python3
"""Replot Chapter 5 from frozen summary-data results; never refit.

Only reports/chapter05_sources/ is read.  Raw King rows and existing fit
parameters are copied unchanged from the earlier alpha project; the plotted
clock band is a conditional transfer of today's drift uncertainty.  Hubble
statistics are existing fixed-clock model comparisons, not new likelihoods.
"""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

from make_overview_figures import (
    BLUE, FONT, GOLD, GRAY, INK, PALE_BLUE, PALE_GOLD, PALE_TEAL, TEAL,
)


OUT = Path(__file__).resolve().parent
REPO = OUT.parents[1]
SOURCES = REPO / "reports/chapter05_sources"
plt.rcParams.update({"axes.labelcolor": INK, "xtick.color": INK, "ytick.color": INK,
                     "axes.spines.top": False, "axes.spines.right": False})


def read_inputs():
    # Check the same bytes that were archived, regardless of sibling-workspace
    # changes.  The manifest can include extra audit/source files not plotted.
    manifest = json.loads((SOURCES / "manifest.json").read_text())
    files = manifest["files"]
    for row in files:
        rel = row.get("snapshot_path", row.get("path"))
        if rel is None:
            raise ValueError("Manifest must record a snapshot_path for each input")
        path = SOURCES / rel
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        if digest != row["sha256"]:
            raise ValueError(f"Frozen source hash mismatch: {path}")
    alpha = json.loads((SOURCES / "alpha/results/results.json").read_text())
    table = [r.split() for r in (SOURCES / "alpha/raw/King2012_tablea1.dat").read_text().splitlines() if r.strip()]
    if len(table) != 295 or any(len(row) != 9 for row in table):
        raise ValueError("Unexpected King catalogue shape")
    selected = [r for r in table if int(r[8]) == 0]
    if len(selected) != 293:
        raise ValueError("Expected original authors' retained 293 rows")
    data = {
        "z": np.array([float(r[3]) for r in selected]),
        "y": np.array([float(r[4]) * 10 for r in selected]),
        "sigma": np.hypot([float(r[5]) * 10 for r in selected],
                          [{1: 0, 2: 17.43, 3: 9.05}[int(r[7])] for r in selected]),
        "telescope": np.array([r[6] for r in selected]),
    }
    with (SOURCES / "hubble/comparison.csv").open() as stream:
        hubble = list(csv.DictReader(stream))
    null = json.loads((SOURCES / "hubble/null_profile.json").read_text())["best"]
    return alpha, data, hubble, null


def save(fig, name):
    for ext in ("pdf", "svg", "png"):
        target = OUT / f"{name}.{ext}"
        fig.savefig(target, dpi=240,
                    metadata={"Creator": Path(__file__).name} if ext == "pdf" else None)
        if ext == "svg":
            target.write_text("\n".join(line.rstrip() for line in target.read_text().splitlines()) + "\n")
    plt.close(fig)


def transfer(z, background):
    """Delta-alpha/alpha divided by present drift; years, same saved background."""
    om = background["Omega_m"]
    t0 = background["t0_years"]
    l0 = background["L0"]
    ratio = np.arcsinh(np.sqrt((1 - om) / om) / (1 + np.asarray(z)) ** 1.5) / np.arcsinh(np.sqrt((1 - om) / om))
    logratio = np.log(ratio)
    if np.any(l0 + logratio <= 0):
        raise ValueError("Response outside logarithmic domain")
    return -t0 * l0 / 2 * np.expm1(-2 * np.log1p(logratio / l0))


def alpha_constraints(alpha, data):
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(9.5, 4.75), gridspec_kw={"width_ratios": [1.15, 1]})
    fig.subplots_adjust(left=.095, right=.98, bottom=.29, top=.76, wspace=.32)
    fig.suptitle("精细结构常数：历史测量与今日漂移的尺度对照", y=.965, fontsize=14, color=INK)
    fig.text(.5, .895, "旧类星体测量与原子钟条件外推分开显示；两栏使用不同纵轴尺度", ha="center", fontsize=10.3, color=GRAY)
    z = np.linspace(0, 4.2, 401)
    response = transfer(z, alpha["background"])
    for name, color in [("Keck", BLUE), ("VLT", GOLD)]:
        mask = data["telescope"] == name
        ax.errorbar(data["z"][mask], data["y"][mask], yerr=data["sigma"][mask],
                    fmt="o", ms=2.5, lw=.6, alpha=.53, color=color,
                    label=f"{name} ({sum(mask)})", capsize=0)
    fit = next(r for r in alpha["fits"]["king"] if r["model"] == "log2")
    ax.plot(z, 1e-13 * response * fit["d"], color="#9D3B40", lw=1.6, label="类星体单独拟合")
    ax.axhline(0, color=GRAY, lw=.8, ls="--")
    lo, hi = np.min(data["y"] - data["sigma"]), np.max(data["y"] + data["sigma"])
    margin = .06 * (hi - lo)
    ax.set_ylim(lo - margin, hi + margin)
    ax.set_xlim(-.07, 4.3)
    ax.set_title("(a) 293 条记录及完整误差棒", fontsize=10.5, pad=9)
    ax.set_xlabel(r"吸收体红移 $z$")
    ax.set_ylabel(r"$\Delta\alpha/\alpha_0\;(\mathrm{ppm})$")
    ax.legend(fontsize=8.8, loc="upper left", frameon=False)
    ax.grid(axis="y", color="#E5E9ED", lw=.6)

    d, sigma = alpha["clock"]
    mean = 1e-13 * response * d
    halfwidth = 1.96e-13 * np.abs(response) * sigma
    bx.fill_between(z, mean - halfwidth, mean + halfwidth, color=TEAL, alpha=.20,
                    label="95% 逐点区间")
    bx.plot(z, mean, color=TEAL, lw=1.8, label="原子钟条件中心值")
    bx.axhline(0, color=GRAY, lw=.8, ls="--")
    bx.set_xlim(-.07, 4.3)
    bx.set_title("(b) 共同时间响应下的外推", fontsize=10.5, pad=9)
    bx.set_xlabel(r"红移 $z$")
    bx.set_ylabel(r"$\Delta\alpha/\alpha_0\;(\mathrm{ppm})$")
    bx.ticklabel_format(axis="y", style="plain", useOffset=False)
    bx.grid(axis="y", color="#E5E9ED", lw=.6)
    bx.legend(fontsize=8.8, loc="lower left", frameon=False)
    fig.text(.095, .123,
             r"固定背景与 $t_*$；左：已存档的逆对数平方拟合。右：同质、无屏蔽假设下的原子钟约束。",
             ha="left", fontsize=9.5, color=INK)
    fig.text(.095, .068, "左侧沿用作者的 293 行选择及既定额外散布预算；右侧曲线不是高红移测量或检出。",
             ha="left", fontsize=9.5, color=GRAY)
    save(fig, "chapter05_alpha_constraints_cn")


def response_extrapolation(alpha):
    bg = alpha["background"]
    l0 = bg["L0"]
    ratio = np.logspace(-2, 12, 601)
    fig, ax = plt.subplots(figsize=(9, 4.25))
    fig.subplots_adjust(left=.12, right=.98, bottom=.30, top=.76)
    fig.suptitle("候选响应的归一化延拓", fontsize=14, color=INK, y=.963)
    fig.text(.5, .875, r"$r_s(t)=[L_0/(L_0+\ln(t/t_0))]^s,\qquad L_0=\ln(t_0/t_*)$",
             ha="center", fontsize=12.5, color=INK)
    ax.axvspan(1, ratio[-1], facecolor=PALE_GOLD, zorder=0)
    for s, color, style in [(1, GRAY, "--"), (2, BLUE, "-"), (3, TEAL, "-.")]:
        ax.plot(ratio, (l0 / (l0 + np.log(ratio))) ** s,
                color=color, lw=2.1 if s == 2 else 1.5, ls=style, label=rf"$s={s}$")
    ax.axvline(1, color=GRAY, lw=.8, ls=":")
    ax.scatter([1], [1], color=BLUE, s=23, zorder=4)
    ax.text(1.6, 1.026, r"今天：$r_s(t_0)=1$", fontsize=9.7, color=INK)
    ax.text(3e5, 1.086, "未来：函数延拓，非数据拟合", ha="center", fontsize=10, color=GOLD)
    ax.set_xscale("log")
    ax.set_xlim(ratio[0], ratio[-1])
    ax.set_ylim(.56, 1.14)
    ax.set_xticks([1e-2, 1, 1e4, 1e8, 1e12], [rf"$10^{{{p}}}$" for p in (-2, 0, 4, 8, 12)])
    ax.set_xlabel(r"时间比 $t/t_0$（对数坐标）")
    ax.set_ylabel(r"设定响应的剩余比例 $r_s(t)$")
    ax.grid(axis="y", color="#E1E5E8", lw=.6)
    ax.legend(loc="lower left", ncol=3, frameon=False, fontsize=10)
    mantissa, exponent = f"{bg['tstar_s']:.6e}".split("e")
    fig.text(.12, .122, rf"同源参考：$t_0={bg['t0_years']/1e9:.3f}$ 十亿年，$t_*={mantissa}\times10^{{{int(exponent)}}}\,\mathrm{{s}}$，$L_0={l0:.3f}$。",
             fontsize=9.7, color=INK)
    fig.text(.12, .061, "曲线仅显示预先设定的响应函数，不表示宇宙完成比例、剩余寿命或冻结年代。",
             fontsize=9.7, color=GRAY)
    save(fig, "chapter05_response_extrapolation_cn")


def card(ax, x, y, w, h, title, lines=(), *, color=BLUE, fill=PALE_BLUE, dashed=False):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=.035,rounding_size=.07",
                               linewidth=1.05, edgecolor=color, facecolor=fill,
                               linestyle=(0, (4, 3)) if dashed else "solid", zorder=3))
    ax.text(x + w / 2, y + h - .21, title, ha="center", va="center", fontsize=11.5, color=color, zorder=4)
    for i, line in enumerate(lines):
        ty = y + (h - .42) * (1 - (i + .5) / len(lines))
        ax.text(x + w / 2, ty, line, ha="center", va="center", fontsize=10.3, color=INK, zorder=4)


def arrow(ax, start, end, color=GOLD):
    ax.add_patch(FancyArrowPatch(start, end, arrowstyle="-|>", mutation_scale=12,
                               linewidth=1.2, color=color, linestyle=(0, (4, 3)), zorder=5))


def future_conditions():
    fig = plt.figure(figsize=(9, 4.15))
    ax = fig.add_axes([.025, .04, .95, .93])
    ax.set_xlim(0, 9); ax.set_ylim(0, 4.15); ax.axis("off")
    ax.text(4.5, 3.98, "从有限数据到未来极限：需要哪些条件", ha="center", va="center", fontsize=14, color=INK)
    ax.text(4.5, 3.61, "观测区间内的拟合与无限未来的物理结论分开建立", ha="center", va="center", fontsize=10.3, color=GRAY)
    card(ax, .16, 1.74, 2.41, 1.41, "有限红移拟合", ["距离与校准似然", r"约束所选归一化 $g(t)$", "记录模型与参数退化"])
    card(ax, 3.29, 1.74, 2.41, 1.41, "未来延拓条件", ["正性与函数定义域", "完整背景与物质演化", "扰动、稳定性与适用范围"], color=GOLD, fill=PALE_GOLD, dashed=True)
    card(ax, 6.42, 1.74, 2.41, 1.41, "条件极限", [r"膨胀归一化极限 $g_\infty$", r"真实膨胀率极限 $H_\infty$", "分别求得，不能直接等同"], color=TEAL, fill=PALE_TEAL)
    arrow(ax, (2.63, 2.43), (3.23, 2.43))
    arrow(ax, (5.76, 2.43), (6.36, 2.43))
    ax.text(4.5, 1.28, r"若 $H(t)=H_0g(t)E_{\rm ref}[z(t)]$，则 $H_\infty=H_0\lim_{t\to\infty}\{g(t)E_{\rm ref}[z(t)]\}$。",
            ha="center", va="center", fontsize=12, color=INK)
    ax.text(4.5, .79, "还需证明极限存在，并验证参考背景在未来延拓中的适用性。", ha="center", va="center", fontsize=10.6, color=INK)
    ax.text(4.5, .33, "本图是推断流程；未给出终值、宇宙寿命或计算冻结时刻。", ha="center", va="center", fontsize=10.2, color=GRAY)
    save(fig, "chapter05_future_conditions_cn")


def hubble_comparison(rows, null):
    chosen = [next(r for r in rows if r["model"] == model and r["case"] == "1.0_future_positive")
              for model in ("logarithmic", "exponential", "power_law")]
    q = float(chosen[0]["q"])
    for r in chosen:
        if abs(float(r["q"]) - q) > 1e-13:
            raise ValueError("Comparison clock differs")
        if abs(null["chi2"] - float(r["chi2"]) - float(r["delta_chi2"])) > 1e-8:
            raise ValueError("Common null inconsistency")
        if abs(float(r["delta_AIC"]) - (2 - float(r["delta_chi2"]))) > 1e-10:
            raise ValueError("AIC comparison not one added parameter")
    vals = [[0] + [float(r[field]) for r in chosen] for field in ("delta_chi2", "delta_AIC")]
    colors = [GRAY, BLUE, GOLD, TEAL]
    labels = ["常数归一化", "逆对数平方", "指数", "幂律"]
    fig, axes = plt.subplots(1, 2, figsize=(9, 4.6), sharey=True)
    fig.subplots_adjust(left=.15, right=.97, bottom=.32, top=.75, wspace=.25)
    fig.suptitle("固定时钟下的校准背景模型比较", fontsize=14, color=INK, y=.965)
    fig.text(.5, .88, r"同一校准似然、固定 $q=0.951245$；选择允许未来归一化为正的参数分支",
             ha="center", fontsize=10.4, color=GRAY)
    for i, (ax, values) in enumerate(zip(axes, vals)):
        yy = np.arange(4)
        ax.barh(yy, values, height=.56, color=colors, alpha=.87)
        ax.axvline(0, color=GRAY, lw=.8)
        ax.scatter([0], [0], color=GRAY, s=16, zorder=4)
        for y, value in zip(yy, values):
            ax.text(value + (.13 if value >= 0 else -.13), y, f"{value:+.3f}" if i else f"{value:.3f}",
                    ha="left" if value >= 0 else "right", va="center", fontsize=10, color=INK)
        ax.set_yticks(yy, labels)
        ax.grid(axis="x", color="#E6EAED", lw=.6)
        ax.set_axisbelow(True)
        ax.set_ylim(3.55, -.55)
        if i == 0:
            ax.set_xlim(-.12, 5.95)
            ax.set_title("(a) 拟合改善：较大更好", fontsize=10.7)
            ax.set_xlabel(r"$\Delta\chi^2=\chi^2_{\rm null}-\chi^2_{\rm model}$")
        else:
            ax.set_xlim(-4.1, 2.3)
            ax.set_title("(b) 参数代价后：较小更好", fontsize=10.7)
            ax.set_xlabel(r"$\Delta\mathrm{AIC}=\mathrm{AIC}_{\rm model}-\mathrm{AIC}_{\rm null}$")
    fig.text(.15, .15, r"$N=1671$ 个观测向量分量（$1657+14$）；完整协方差；常数模型 $k=4$，其余 $k=5$。",
             fontsize=9.4, color=INK)
    fig.text(.15, .097, "指数最优解为负斜率；此处不强制归一化随时间不减。", fontsize=9.4, color=INK)
    fig.text(.15, .044, "记录共享校准，N 不等于独立天体数；图示为作者条件模型比较，未转换为发现显著性。",
             fontsize=9.4, color=GRAY)
    save(fig, "chapter05_hubble_comparison_cn")


if __name__ == "__main__":
    alpha, data, hubble, null = read_inputs()
    # Re-evaluate published transfer values as a plotting consistency check,
    # without adjusting a model, solving an optimizer or rerunning a fit.
    for row in alpha["clock_transfer"]:
        check = float(transfer(row["z"], alpha["background"]))
        if not np.isclose(check, row["T_years"], rtol=1e-12, atol=1e-5):
            raise ValueError("Transfer formula disagrees with archived values")
    alpha_constraints(alpha, data)
    response_extrapolation(alpha)
    future_conditions()
    hubble_comparison(hubble, null)
    print(f"Font: {FONT}; King rows: {len(data['z'])}; no fits run")
    print(f"Saved four Chapter 5 diagrams as PDF/SVG/PNG in {OUT}")
