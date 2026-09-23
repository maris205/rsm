#!/usr/bin/env python3
"""Reproduce the three Chinese Chapter 2 schematics (not new fit results).

Run from any directory.  Uses the same fonts and palette as the front matter.
The prime indicator is exact for the displayed integers.  Parameter landmarks
are reference values for the autonomous quadratic family, not a new scan or a
fit of arithmetic data.  Dashed gold links denote added physical hypotheses.
"""

from __future__ import annotations

import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle

from make_overview_figures import (
    BLUE, FONT, GOLD, GRAY, INK, PALE_BLUE, PALE_GOLD, PALE_TEAL, TEAL,
)


OUT = Path(__file__).resolve().parent


def canvas(title, subtitle, height=5.6):
    fig = plt.figure(figsize=(9, height))
    ax = fig.add_axes([0.025, 0.035, 0.95, 0.94])
    ax.set_xlim(0, 9)
    ax.set_ylim(0, height)
    ax.axis("off")
    ax.text(4.5, height - .16, title, ha="center", va="center", fontsize=14, color=INK)
    ax.text(4.5, height - .49, subtitle, ha="center", va="center", fontsize=9.7, color=GRAY)
    return fig, ax


def card(ax, x, y, w, h, title, lines=(), *, color=BLUE, fill=PALE_BLUE, dashed=False):
    ax.add_patch(FancyBboxPatch(
        (x, y), w, h, boxstyle="round,pad=.035,rounding_size=.075",
        linewidth=1.05, edgecolor=color, facecolor=fill,
        linestyle=(0, (4, 3)) if dashed else "solid", zorder=3,
    ))
    ax.text(x + w / 2, y + h - .20, title, ha="center", va="center", fontsize=11,
            color=color, zorder=4)
    for i, line in enumerate(lines):
        ty = y + (h - .41) * (1 - (i + .5) / len(lines))
        ax.text(x + w / 2, ty, line, ha="center", va="center", fontsize=10,
                color=INK, zorder=4)


def arrow(ax, start, end, *, color=BLUE, dashed=False):
    ax.add_patch(FancyArrowPatch(
        start, end, arrowstyle="-|>", mutation_scale=11, linewidth=1.2,
        color=color, linestyle=(0, (4, 3)) if dashed else "solid", zorder=5,
    ))


def save(fig, name):
    for ext in ("pdf", "svg", "png"):
        fig.savefig(OUT / f"{name}.{ext}", dpi=240,
                    metadata={"Creator": Path(__file__).name} if ext == "pdf" else None)
    plt.close(fig)


def is_prime(n):
    return n >= 2 and all(n % d for d in range(2, math.isqrt(n) + 1))


def arithmetic_encoding():
    fig, ax = canvas("从算术筛选到符号动力学", "两个各自明确的编码对象；建立对应需要指定规则并进行检验")
    card(ax, .18, 2.78, 8.64, 1.96, "算术输入：素数指示序列")
    ax.text(.44, 4.14, r"整数 $m$", fontsize=10, color=INK, va="center", zorder=4)
    ax.text(.44, 3.64, r"$b_m$", fontsize=11, color=INK, va="center", zorder=4)
    for m in range(1, 19):
        x = 1.55 + (m - 1) * .397
        ax.text(x, 4.14, str(m), ha="center", va="center", color=INK, fontsize=10, zorder=4)
        prime = is_prime(m)
        ax.add_patch(Rectangle((x - .157, 3.48), .314, .32,
                               facecolor=TEAL if prime else "#DFE5EA", edgecolor="none", zorder=4))
        ax.text(x, 3.64, str(int(prime)), ha="center", va="center",
                color="white" if prime else INK, fontsize=10, zorder=5)
    ax.text(4.5, 3.11, r"$b_m=1$：素数；$b_m=0$：非素数。$m=1$ 既不是素数，也不是合数。",
            ha="center", va="center", fontsize=10, color=INK, zorder=4)
    arrow(ax, (4.5, 2.72), (4.5, 2.21))
    ax.text(4.69, 2.48, "对应／统计检验", va="center", fontsize=10.5, color=BLUE)
    card(ax, .18, .43, 8.64, 1.72, "动力学输入：二次映射轨道与符号分区")
    ax.text(2.10, 1.15, r"$x_{n+1}=1-u_nx_n^2$", ha="center", va="center", fontsize=17, color=INK, zorder=4)
    arrow(ax, (3.90, 1.16), (4.47, 1.16))
    for x, label, sign, color in [(5.02, "L", r"$x_n<0$", TEAL),
                                  (6.48, "R", r"$x_n>0$", BLUE),
                                  (7.93, "C", r"$x_n=0$", GRAY)]:
        ax.text(x, 1.23, label, fontsize=16, color=color, ha="center", va="center", zorder=4)
        ax.text(x, .83, sign, fontsize=11, color=INK, ha="center", va="center", zorder=4)
    ax.text(4.5, .09, "素数与左右分区之间的标签配对属于建模选择；有限数值对应需另行核验。",
            ha="center", va="center", fontsize=9.6, color=GRAY)
    save(fig, "chapter02_arithmetic_encoding_cn")


def parameter_landmarks():
    fig, ax = canvas("二次映射中的两个参数位置", r"自主参数族 $f_u(x)=1-ux^2$；参数标记示意，不表示新扫描或算术拟合")
    x0, x1 = .80, 8.20
    umin, umax = 1.30, 1.65
    to_x = lambda u: x0 + (u - umin) / (umax - umin) * (x1 - x0)
    line_y = 2.29
    card(ax, .28, 3.00, 3.88, 1.55, "倍周期累积点", [r"$u_F\simeq 1.401155$", "倍周期级联的极限位置"], color=TEAL, fill=PALE_TEAL)
    card(ax, 4.83, 3.00, 3.88, 1.55, "双带合并点", [r"$u_M\simeq 1.543689$", "两条混沌带合并为一条"], color=GOLD, fill=PALE_GOLD)
    ax.plot([x0, x1], [line_y, line_y], color=INK, lw=1.2)
    for u in (1.30, 1.35, 1.40, 1.45, 1.50, 1.55, 1.60, 1.65):
        x = to_x(u)
        ax.plot([x, x], [line_y - .05, line_y + .05], color=GRAY, lw=.8)
        ax.text(x, line_y - .19, f"{u:.2f}", ha="center", va="top", color=GRAY, fontsize=9)
    for u, color in [(1.401155, TEAL), (1.543689, GOLD)]:
        x = to_x(u)
        ax.plot([x, x], [line_y, 2.91], color=color, lw=1.4)
        ax.scatter([x], [line_y], s=42, c=color, zorder=5)
    ax.text(8.42, line_y, r"$u$", va="center", fontsize=12, color=INK)
    ax.text(4.5, 1.66, "周期窗口与带结构需要分别识别，不能用单一阈值划分整个参数区间。",
            ha="center", va="center", fontsize=10, color=INK)
    card(ax, .70, .39, 7.60, .99, "与常见 Logistic 参数的换算", [
        r"$y_{n+1}=r y_n(1-y_n),\quad u=r(r-2)/4$",
        r"$r_F\simeq3.569946,\qquad r_M\simeq3.678574$",
    ])
    ax.text(4.5, .08, "两处标记指向不同的动力学结构；本框架必须说明具体采用哪一个参数基准。",
            ha="center", va="center", fontsize=9.5, color=GRAY)
    save(fig, "chapter02_parameter_landmarks_cn")


def physical_bridge():
    fig, ax = canvas("从算术数值对应到物理检验", "保留研究动机，同时逐层给出新增假设与可计算对象", height=7.2)
    card(ax, .19, 5.19, 4.05, 1.07, "算术来源与谱统计", ["素数结构；黎曼零点", "GUE 统计联系提供启发"])
    card(ax, 4.76, 5.19, 4.05, 1.07, "已发表的数值研究", ["经标定的非自治动力学模型", "有限区间的数值谱对应"], color=TEAL, fill=PALE_TEAL)
    arrow(ax, (2.21, 5.13), (3.58, 4.66), color=GOLD, dashed=True)
    arrow(ax, (6.79, 5.13), (5.42, 4.66), color=GOLD, dashed=True)
    card(ax, 1.00, 3.63, 7.00, .97, "新增假设：内部时钟与响应规律", [
        r"$M^2(R)$ 含 $1/\ln^2(R/R_*)$；需定义索引、内部时钟与物理时间的关系",
    ], color=GOLD, fill=PALE_GOLD, dashed=True)
    arrow(ax, (4.5, 3.57), (4.5, 3.15))
    card(ax, 1.00, 2.25, 7.00, .84, "给定模型：时钟—物质的 Hamilton 耦合", [
        "明确自由度、辛结构、能量交换与连续近似",
    ])
    arrow(ax, (4.5, 2.19), (4.5, 1.80), color=GOLD, dashed=True)
    card(ax, .49, .88, 8.02, .86, "指定物理读出与独立检验", [
        r"原子钟／光谱的 $\alpha$、膨胀史 $H(z)$、$\mathrm{CMB}$：标定、常规模型与留出数据对照",
    ], color=GOLD, fill=PALE_GOLD, dashed=True)
    ax.text(4.5, .54, r"$\mathrm{GUE}$ 统计联系不等同于物理同构；算术索引 $n$ 不自动等同于宇宙时间 $t$。",
            ha="center", va="center", fontsize=10, color=INK)
    ax.text(4.5, .16, "实线：给定模型内的构造；虚线：新增假设或待建立的物理对应。",
            ha="center", va="center", fontsize=9.8, color=GRAY)
    save(fig, "chapter02_physical_bridge_cn")


if __name__ == "__main__":
    arithmetic_encoding()
    parameter_landmarks()
    physical_bridge()
    print(f"Font: {FONT}")
    print(f"Saved 3 diagrams as PDF/SVG/PNG in {OUT}")
