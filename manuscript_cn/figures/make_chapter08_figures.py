#!/usr/bin/env python3
"""Replot Chapter 8 from frozen results; no fitting or simulation.

Inputs and their SHA-256 digests are listed below and in captions_chapter08.md.
Only summary statistics are transformed for display.  The alpha panel selects
the fixed-s=2 rejection test, never the scanned-exponent test.  Arithmetic
ranges remain the frozen empirical central ranges, not confidence intervals.
Run from any directory: python /path/to/make_chapter08_figures.py
"""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import PercentFormatter
import numpy as np

from make_overview_figures import BLUE, GRAY, INK, PALE_GOLD, TEAL


OUT = Path(__file__).resolve().parent
REPO = OUT.parents[1]
ALPHA = "experiments/alpha_recovery_v01"
ARITH = "experiments/arithmetic_residual_v01"
MAIN = "king_clock_offsets"
INPUT_SHA256 = {
    f"{ALPHA}/results/recovery.csv": "70c5749a40e1c3e03563ff6d5e8c9a137bdef45d53154f355d0cf5eda1f057a7",
    f"{ALPHA}/results/summary.json": "c3806224d57c034218093da24e80c8c76a0f39791078ec61eaba3b08a2e25a6a",
    f"{ALPHA}/results/shape_geometry.csv": "e15b14f65cb721a202bef546ce03949d8b001f5989ead8cb975064363209e3bf",
    f"{ALPHA}/protocol.md": "3e0b6ecff2c73a9e830a4b7b481d07753e5953b86572c8c02ac56c73a1d1e100",
    f"{ARITH}/results/comparison_blocks.csv": "b1234f062c46209efb0acbf452c3ec65036b66b5d0b24994abbed26c8dc8e7e8",
    f"{ARITH}/results/comparison_summary.json": "278a9ae61df1d4fba46b4db4f240c5adc34e139213a5bd976b15904942d2f610",
    f"{ARITH}/protocol.md": "952804fa8aaf7ffa2092449010dc8a8c5d045ed319df1a7812a993ee565e9d56",
}

plt.rcParams.update({
    "axes.labelcolor": INK, "xtick.color": INK, "ytick.color": INK,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.edgecolor": GRAY, "axes.labelsize": 10.1,
})


def read_csv(relative):
    with (REPO / relative).open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def check_sources():
    for relative, expected in INPUT_SHA256.items():
        actual = hashlib.sha256((REPO / relative).read_bytes()).hexdigest()
        if actual != expected:
            raise ValueError(f"Frozen input hash mismatch: {relative}")


def select(recovery, true_s, d):
    records = [row for row in recovery if row["design"] == MAIN
               and int(row["true_s"]) == true_s and float(row["d_injected"]) == d]
    if len(records) != 1:
        raise ValueError((true_s, d, len(records)))
    return records[0]


def save(fig, stem):
    for extension in ("pdf", "svg", "png"):
        path = OUT / f"{stem}.{extension}"
        fig.savefig(path, dpi=240,
                    metadata={"Creator": Path(__file__).name} if extension == "pdf" else None)
        if extension == "svg":
            path.write_text("\n".join(line.rstrip() for line in path.read_text().splitlines()) + "\n")
        print(path.relative_to(REPO))
    plt.close(fig)


def recovery_figure():
    rows = read_csv(f"{ALPHA}/results/recovery.csv")
    summary = json.loads((REPO / ALPHA / "results/summary.json").read_text())
    geometry = read_csv(f"{ALPHA}/results/shape_geometry.csv")
    assert summary["all_checks_passed"] and summary["n_evaluation"] == 20000
    records = [select(rows, 0, 0)] + [select(rows, 2, d) for d in (2.5, 5, 12.5)]
    x = np.array([float(r["d_injected"]) for r in records]) / summary["current_clock_sigma_d"]
    y = np.array([float(r["fixed_s2_detection"]) for r in records])
    err = 1.96 * np.array([float(r["detection_mc_se"]) for r in records])
    theoretical = np.array([float(r["fixed_s2_analytic_power"]) for r in records])
    assert np.allclose(x, [0, 1, 2, 5])
    assert np.allclose(y, [.04795, .1696, .5176, .99885])
    # This audit uses the archived no-noise residuals but does not display them.
    nearest = min(float(r["delta_chi2_at_d12p5"]) for r in geometry
                  if r["design"] == MAIN and r["true_s"] == "2")
    assert np.isclose(nearest, 3.4086648410236295e-9, rtol=1e-12, atol=0)

    fig, (ax, bx) = plt.subplots(1, 2, figsize=(9.6, 4.8),
                                gridspec_kw={"width_ratios": [1, 1.08]})
    fig.subplots_adjust(left=.087, right=.975, bottom=.30, top=.76, wspace=.42)
    fig.suptitle("合成恢复：检出漂移幅度与识别时间形状", y=.964, fontsize=14, color=INK)
    fig.text(.5, .894, "King 293 条记录 + 原子钟；Keck、VLT 各一个自由偏移", ha="center", fontsize=10.1, color=GRAY)

    ax.plot(x, theoretical, ls="none", marker="x", ms=8, mew=1.3,
            color=GRAY, label="同一注入点的解析功效", zorder=2)
    ax.errorbar(x, y, yerr=err, fmt="o", color=BLUE, markersize=5,
                capsize=3, elinewidth=1.2, label="20,000 次合成恢复", zorder=3)
    ax.axhline(.05, color=GRAY, ls=":", lw=1)
    ax.text(4.95, .095, "零假设名义拒绝率 5%", color=GRAY, fontsize=8.5, ha="right")
    for xx, yy, label, offset, align in zip(
        x, y, ["4.795%", "16.960%", "51.760%", "99.885%"],
        [(11, 9), (0, 10), (0, 10), (-3, -17)], ["left", "center", "center", "right"]
    ):
        ax.annotate(label, (xx, yy), xytext=offset, textcoords="offset points",
                    ha=align, fontsize=9.2, color=INK)
    ax.set_title(r"$\mathrm{(a)}$ 固定 $s=2$，检验零漂移", fontsize=10.6, pad=12)
    ax.set(xlim=(-.18, 5.25), ylim=(0, 1.07), xticks=[0, 1, 2, 3, 4, 5],
           xlabel=r"注入漂移 $D_0/\sigma_{\rm clock}$", ylabel="拒绝零漂移的比例")
    ax.yaxis.set_major_formatter(PercentFormatter(1, decimals=0))
    ax.grid(axis="y", color="#E4E8EC", lw=.6)
    ax.legend(loc="upper left", frameon=False, fontsize=8.1)

    candidates = [select(rows, s, 12.5) for s in (1, 2, 3)]
    keys = ["aic_constant", "aic_s1", "aic_s2", "aic_s3"]
    percentages = 100 * np.array([[float(r[k]) for k in keys] for r in candidates])
    assert np.allclose(percentages.sum(axis=1), 100)
    assert np.all(percentages == percentages[0])
    bx.imshow(percentages, vmin=0, vmax=55, cmap="Blues", aspect="auto")
    for (iy, ix), value in np.ndenumerate(percentages):
        bx.text(ix, iy, f"{value:.3f}%", ha="center", va="center", fontsize=9.5,
                color="white" if value > 35 else INK)
    bx.set_title(r"$\mathrm{(b)\ AIC}$ 选择频率，$D_0=5\sigma_{\rm clock}$", fontsize=10.6, pad=12)
    bx.set(xticks=range(4), xticklabels=["常数", "$s=1$", "$s=2$", "$s=3$"],
           yticks=range(3), yticklabels=["$s=1$", "$s=2$", "$s=3$"],
           xlabel="选中的候选模型", ylabel="注入的真实指数")
    bx.set_xticks(np.arange(-.5, 4), minor=True)
    bx.set_yticks(np.arange(-.5, 3), minor=True)
    bx.grid(which="minor", color="white", lw=2)
    bx.tick_params(which="minor", bottom=False, left=False)
    for spine in bx.spines.values():
        spine.set_visible(False)

    fig.text(.087, .157, r"$\sigma_{\rm clock}=2.5\times10^{-19}\,\mathrm{yr}^{-1}$；固定背景和参考尺度。左侧误差棒：约 $95\%$ 蒙特卡洛误差。", fontsize=9.5, color=INK)
    fig.text(.087, .105, "右侧各行共用噪声实现；选择频率不是理论成立概率。常数候选保留，未重新归一化。", fontsize=9.5, color=INK)
    fig.text(.087, .054, "这些是既定观测设计上的合成试验，不是新观测检出；五倍误差注入不表示与今日测量相容。", fontsize=9.3, color=GRAY)
    save(fig, "chapter08_recovery_cn")


def arithmetic_figure():
    rows = read_csv(f"{ARITH}/results/comparison_blocks.csv")
    summary = json.loads((REPO / ARITH / "results/comparison_summary.json").read_text())
    baseline = summary["baseline_T"]
    assert np.isclose(baseline, .04513942199602207, rtol=0, atol=1e-15)
    assert len(rows) == 16 and all(r["inside_central_95"] == "True" for r in rows)
    assert all(r["surrogates"] == "99" for r in rows)
    fig, ax = plt.subplots(figsize=(9.2, 6.05))
    fig.subplots_adjust(left=.16, right=.963, bottom=.245, top=.79)
    fig.suptitle("算术残差：固定终点的能量交换", y=.968, fontsize=14, color=INK)
    fig.text(.5, .909, r"$T=[E_R(20)-E_R(0)]/\mathcal{H}(0)$，$D=T-T_0$；模型时间与能量均无量纲", ha="center", fontsize=11, color=INK)
    fig.text(.5, .86, "八个真实零点残差块，与每块两族各 99 条替代输入比较", ha="center", fontsize=10, color=GRAY)

    ax.axvline(0, color=GRAY, lw=1, ls="--", zorder=1)
    ax.axhspan(-1.50, -.50, color=PALE_GOLD, zorder=0)
    ax.axhline(-.5, color="#D5D9DD", lw=.8)
    positions = list(range(7, -1, -1)) + [-1]
    for family, color, offset in [("iaaft", BLUE, .13), ("exact_spectrum", TEAL, -.13)]:
        ordered = sorted([r for r in rows if r["family"] == family], key=lambda r: int(r["block"]))
        records = ordered + [summary["aggregate"][family]]
        for row, position in zip(records, positions):
            values = [(float(row[k]) - baseline) * 1e4 for k in
                      ("null_q025", "null_q25", "null_median", "null_q75", "null_q975")]
            yy = position + offset
            ax.hlines(yy, values[0], values[4], color=color, lw=1.5, zorder=2)
            ax.vlines([values[0], values[4]], yy-.045, yy+.045, color=color, lw=1)
            ax.hlines(yy, values[1], values[3], color=color, lw=3.2, alpha=.75, zorder=3)
            ax.plot(values[2], yy, "o", ms=3.6, color=color, zorder=4)
    real = sorted([r for r in rows if r["family"] == "iaaft"], key=lambda r: int(r["block"]))
    real_values = [float(r["real_D"]) * 1e4 for r in real]
    real_values.append(summary["aggregate"]["iaaft"]["real_D"] * 1e4)
    assert np.isclose(np.mean(real_values[:-1]), real_values[-1], rtol=0, atol=1e-11)
    ax.plot(real_values, positions, "D", ms=4.4, color="#111111", zorder=5)
    ax.set(yticks=positions, yticklabels=[f"评价块 {j}" for j in range(1, 9)] + ["八块等权均值"],
           xlabel=r"相对纯包络基线的差 $10^4 D$", xlim=(-9.6, 10.1), ylim=(-1.5, 7.6))
    ax.set_xticks([-8, -4, 0, 4, 8])
    ax.grid(axis="x", color="#E4E8EC", lw=.6, zorder=0)
    handles = [Line2D([], [], marker="D", color="#111111", ls="none", markersize=4.5, label="真实输入"),
               Line2D([], [], marker="o", color=BLUE, lw=1.5, markersize=3.6, label="IAAFT 替代族"),
               Line2D([], [], marker="o", color=TEAL, lw=1.5, markersize=3.6, label="配对的精确节点谱族")]
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(.5, -.14), frameon=False, ncol=3, fontsize=9.1)
    fig.text(.07, .111, "细条：经验 2.5%—97.5% 范围；粗条：四分位范围；圆点：中位数。不是参数置信区间。", fontsize=9.4, color=INK)
    fig.text(.07, .067, "汇总范围按同一替代序号跨块平均后计算。两族相互配对；四个节点质量失败仍保留。", fontsize=9.4, color=INK)
    fig.text(.07, .023, "窗化连续驱动及其导数的谱未严格匹配；图示结果未建立黎曼特异效应，也不证明等价。", fontsize=9.4, color=GRAY)
    save(fig, "chapter08_arithmetic_cn")


def main():
    check_sources()
    recovery_figure()
    arithmetic_figure()


if __name__ == "__main__":
    main()
