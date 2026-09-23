#!/usr/bin/env python3
"""Draw the Chapter 4 lattice/long-wave schematic; no numerical experiment.

Equations follow manuscript_v01/chapters/03_autonomous_dynamics.md and
04_continuum_limit.md.  Canonical normalization and the common clock source
are retained.  The response, including optional arithmetic residual coupling,
is a chosen input; the diagram does not derive it uniquely from Riemann data.
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


def card(ax, x, y, w, h, title, lines=(), *, color=BLUE, fill=PALE_BLUE,
         dashed=False, textsize=10.8, titlesize=11.2):
    ax.add_patch(FancyBboxPatch(
        (x, y), w, h, boxstyle="round,pad=.035,rounding_size=.07",
        linewidth=1.05, edgecolor=color, facecolor=fill,
        linestyle=(0, (4, 3)) if dashed else "solid", zorder=3,
    ))
    ax.text(x + w / 2, y + h - .20, title, ha="center", va="center",
            fontsize=titlesize, color=color, zorder=4)
    for i, line in enumerate(lines):
        ty = y + (h - .40) * (1 - (i + .5) / len(lines))
        ax.text(x + w / 2, ty, line, ha="center", va="center", fontsize=textsize,
                color=INK, zorder=4)


def arrow(ax, start, end, *, color=BLUE, dashed=False):
    ax.add_patch(FancyArrowPatch(
        start, end, arrowstyle="-|>", mutation_scale=11, linewidth=1.2,
        color=color, linestyle=(0, (4, 3)) if dashed else "solid", zorder=5,
    ))


def draw():
    fig = plt.figure(figsize=(9.0, 5.8))
    ax = fig.add_axes([.025, .035, .95, .94])
    ax.set_xlim(0, 9)
    ax.set_ylim(0, 5.8)
    ax.axis("off")
    ax.text(4.5, 5.62, "同一相互作用：从格点到长波场", ha="center", va="center", fontsize=14, color=INK)
    ax.text(4.5, 5.28, "指定正则耦合后推导宏观方程；逆对数平方与算术耦合仍是输入假设",
            ha="center", va="center", fontsize=10.2, color=GRAY)

    card(ax, .16, 4.26, 8.68, .76, "算术研究启发的候选响应", [
        r"$M^2(R)=M_{\rm env}^2(R)+\varepsilon_{\rm ar}w(R)S(R),\qquad"
        r" M_{\rm env}^2(R)=M_\infty^2+b/\ln^2(R/R_*)$",
    ], color=GOLD, fill=PALE_GOLD, dashed=True, textsize=11.4)
    arrow(ax, (1.54, 4.20), (1.54, 3.96), color=GOLD, dashed=True)

    card(ax, .16, 2.08, 2.66, 1.82, "正则格点与共享时钟", [
        r"$p_\ell=a\dot q_\ell,\qquad P_R=I\dot R$",
        r"$V_{\rm int}=\frac{1}{2} M^2(R)Q_2$",
        r"$Q_2=a\sum_\ell q_\ell^2$",
    ], textsize=12)
    card(ax, 3.17, 2.08, 2.66, 1.82, "同一 Hamilton 量的方程", [
        r"$\dot q_\ell=\partial_{p_\ell}H$",
        r"$\dot p_\ell=-\partial_{q_\ell}H$",
        r"$I\ddot R=-\frac{1}{2}[M^2]'(R)Q_2$",
    ], textsize=12)
    card(ax, 6.18, 2.08, 2.66, 1.82, "长波场与积分反馈", [
        r"$\phi_{\tau\tau}=v^2\phi_{xx}-M^2(R)\phi-g\phi^3$",
        r"$I\ddot R=-\frac{1}{2}[M^2]'(R)\int_0^L\phi^2\,\mathrm{d}x$",
        r"$q_\ell\simeq\phi(x_\ell,\tau),\quad p_\ell/a\simeq\phi_\tau$",
    ], color=TEAL, fill=PALE_TEAL, textsize=10.6)
    arrow(ax, (2.88, 2.99), (3.11, 2.99))
    arrow(ax, (5.89, 2.99), (6.12, 2.99), color=TEAL)

    ax.text(4.5, 1.79,
            r"固定 $L,I$ 及模型系数；在光滑、受控的有限时窗内，格点求和对应同一场积分。",
            ha="center", va="center", fontsize=10.2, color=INK)

    card(ax, .16, .47, 2.66, 1.03, "空间尺度：格距", [
        r"$a=L/N,\quad |\kappa|a\ll1$",
        r"首修正 $+(v^2a^2/12)\phi_{xxxx}$",
    ], textsize=10.8)
    card(ax, 3.17, .47, 2.66, 1.03, "计算尺度：时间步长", [
        r"$h$ 控制时间离散误差",
        "独立加密；与空间修正分开",
    ], textsize=10.6)
    card(ax, 6.18, .47, 2.66, 1.03, "统计平均：另需闭合", [
        r"一般 $\langle\phi^3\rangle\neq\langle\phi\rangle^3$",
        "需定义统计态与约化规则",
    ], color=GOLD, fill=PALE_GOLD, dashed=True, textsize=10.6)
    ax.text(4.5, .13, r"共享 $R$ 仍是均匀模式：连续描述保留全局积分反馈，局域驱动场尚待建立。",
            ha="center", va="center", fontsize=10.2, color=GRAY)

    for ext in ("pdf", "svg", "png"):
        target = OUT / f"chapter04_micro_macro_bridge_cn.{ext}"
        fig.savefig(target, dpi=240,
                    metadata={"Creator": Path(__file__).name} if ext == "pdf" else None)
        if ext == "svg":
            # Matplotlib emits whitespace-only indentation lines; keep tracked
            # SVGs clean under git diff --check without changing their content.
            target.write_text("\n".join(line.rstrip() for line in target.read_text().splitlines()) + "\n")
    plt.close(fig)


if __name__ == "__main__":
    draw()
    print(f"Font: {FONT}")
    print(f"Saved Chapter 4 diagram as PDF/SVG/PNG in {OUT}")
