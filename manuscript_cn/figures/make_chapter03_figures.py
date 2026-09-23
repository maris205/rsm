#!/usr/bin/env python3
"""Draw the Chapter 3 update/readout schematic; no new experiment is run.

The K-D-K stages represent the separable Hamiltonian benchmark already defined
in experiments/log_clock_coupling_v01/reports/derivation_and_limits.md, Sec. 9.
Every kick evaluates all forces at the same configuration for that stage.
The shared R is a global reduced degree of freedom, not a local causal field.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

from make_overview_figures import (
    BLUE, FONT, GOLD, GRAY, INK, PALE_BLUE, PALE_GOLD, PALE_TEAL, TEAL,
)


OUT = Path(__file__).resolve().parent


def card(ax, x, y, w, h, title, lines=(), *, color=BLUE, fill=PALE_BLUE, dashed=False):
    ax.add_patch(FancyBboxPatch(
        (x, y), w, h, boxstyle="round,pad=.035,rounding_size=.07",
        linewidth=1.05, edgecolor=color, facecolor=fill,
        linestyle=(0, (4, 3)) if dashed else "solid", zorder=3,
    ))
    ax.text(x + w / 2, y + h - .20, title, ha="center", va="center", fontsize=10.7,
            color=color, zorder=4)
    for i, line in enumerate(lines):
        ty = y + (h - .39) * (1 - (i + .5) / len(lines))
        ax.text(x + w / 2, ty, line, ha="center", va="center", fontsize=9.8,
                color=INK, zorder=4)


def arrow(ax, start, end, *, color=BLUE, dashed=False):
    ax.add_patch(FancyArrowPatch(
        start, end, arrowstyle="-|>", mutation_scale=10.5, linewidth=1.15,
        color=color, linestyle=(0, (4, 3)) if dashed else "solid", zorder=5,
    ))


def draw():
    fig = plt.figure(figsize=(9.0, 5.55))
    ax = fig.add_axes([.025, .035, .95, .94])
    ax.set_xlim(0, 9)
    ax.set_ylim(0, 5.55)
    ax.axis("off")

    ax.text(4.5, 5.37, "共同更新规则与独立物理读数", ha="center", va="center", fontsize=14, color=INK)
    ax.text(4.5, 5.04, "计算流程定义变量怎样更新；物理接口还需定义变量怎样被观测", ha="center", va="center", fontsize=10, color=GRAY)
    ax.text(4.5, 4.63,
            r"$H=T+V$；纯包络 $M^2(R)=M_\infty^2+b/\ln^2(R/R_*)$ 同时进入场力与时钟反馈",
            ha="center", va="center", fontsize=10.8, color=INK)

    card(ax, .15, 3.12, 1.99, 1.14, "完整状态", [
        r"$Z^k=(q_\ell,p_\ell,R,P_R)^k$", "场与时钟共同进入状态",
    ])
    card(ax, 2.39, 3.12, 1.99, 1.14, "同一配置评估", [
        r"$-\partial_{q_\ell}V,\;-\partial_RV$", "每一子步使用同一配置",
    ])
    card(ax, 4.63, 3.12, 1.99, 1.14, "完整更新", [
        r"$Z^{k+1}=\Phi_h(Z^k)$", "场与时钟均被推进",
    ], color=TEAL, fill=PALE_TEAL)
    card(ax, 6.87, 3.12, 1.99, 1.14, "独立读数", [
        r"$\mathcal{O}=\mathcal{R}[Z]$", "指定测量与校准规则",
    ], color=GOLD, fill=PALE_GOLD, dashed=True)
    arrow(ax, (2.19, 3.68), (2.34, 3.68))
    arrow(ax, (4.43, 3.68), (4.58, 3.68))
    arrow(ax, (6.67, 3.68), (6.82, 3.68), color=GOLD, dashed=True)

    ax.text(.17, 2.77, "完整步的分解", fontsize=10.3, color=TEAL, ha="left", va="center")
    for x, title, detail in [
        (1.64, r"动量半步 $K(h/2)$", r"同一 $(q^k,R^k)$ 的力"),
        (4.12, r"坐标整步 $D(h)$", r"共同推进 $(q,R)$"),
        (6.60, r"动量半步 $K(h/2)$", r"同一 $(q^{k+1},R^{k+1})$ 的力"),
    ]:
        card(ax, x, 2.27, 2.22, .72, title, [detail], color=TEAL, fill=PALE_TEAL)
    arrow(ax, (3.91, 2.62), (4.07, 2.62), color=TEAL)
    arrow(ax, (6.39, 2.62), (6.55, 2.62), color=TEAL)
    ax.text(4.5, 1.98, r"$h$ 是数值步长，不是已证实的基本时间节拍；辛更新仍需误差与能量核验。",
            ha="center", va="center", fontsize=9.8, color=GRAY)

    card(ax, .15, .72, 4.24, .91, "计算依赖图", [
        "哪些变量进入同一个演化方程右端", "箭头表示求值依赖与数据流",
    ])
    card(ax, 4.63, .72, 4.23, .91, "物理接口", [
        "空间、传播、钟读数与测量概率需定义", "物理因果关系需独立建立与检验",
    ], color=GOLD, fill=PALE_GOLD, dashed=True)
    ax.text(4.5, .32, r"共享 $R$ 是当前全局约化；物理局域驱动场及其传播规律尚待建立。",
            ha="center", va="center", fontsize=10, color=INK)

    for ext in ("pdf", "svg", "png"):
        path = OUT / f"chapter03_update_and_readout_cn.{ext}"
        fig.savefig(path, dpi=240,
                    metadata={"Creator": Path(__file__).name} if ext == "pdf" else None)
        if ext == "svg":
            # Matplotlib puts trailing spaces in path data; keep versioned text clean.
            path.write_text("\n".join(line.rstrip() for line in path.read_text().splitlines()) + "\n")
    plt.close(fig)


if __name__ == "__main__":
    draw()
    print(f"Font: {FONT}")
    print(f"Saved Chapter 3 diagram as PDF/SVG/PNG in {OUT}")
