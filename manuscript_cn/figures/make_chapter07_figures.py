#!/usr/bin/env python3
"""Draw Chapter 7's evidence pipeline; no simulation or fit is performed.

The solid paths organize the existing log_clock_coupling_v01 and
arithmetic_residual_v01 scripts.  Their independent DOP853 implementations
consume the same defined cases, not the main integrator's implementation.
The dashed ECS/Gym and physical-readout interfaces are proposed work.
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
         dashed=False):
    ax.add_patch(FancyBboxPatch(
        (x, y), w, h, boxstyle="round,pad=.035,rounding_size=.07",
        linewidth=1.05, edgecolor=color, facecolor=fill,
        linestyle=(0, (4, 3)) if dashed else "solid", zorder=3,
    ))
    ax.text(x + w / 2, y + h - .20, title, ha="center", va="center",
            fontsize=10.8, color=color, zorder=4)
    for i, line in enumerate(lines):
        ty = y + (h - .40) * (1 - (i + .5) / len(lines))
        ax.text(x + w / 2, ty, line, ha="center", va="center", fontsize=9.7,
                color=INK, zorder=4)


def arrow(ax, start, end, *, color=BLUE, dashed=False):
    ax.add_patch(FancyArrowPatch(
        start, end, arrowstyle="-|>", mutation_scale=10.5, linewidth=1.10,
        color=color, linestyle=(0, (4, 3)) if dashed else "solid", zorder=2,
    ))


def draw():
    fig = plt.figure(figsize=(9.0, 5.8))
    ax = fig.add_axes([.025, .035, .95, .94])
    ax.set_xlim(0, 9)
    ax.set_ylim(0, 5.8)
    ax.axis("off")

    ax.text(4.5, 5.59, "黎曼引擎：从固定输入到可核查产物",
            ha="center", va="center", fontsize=14, color=INK)
    ax.text(4.5, 5.26, "现有脚本的逻辑整理；尚非已发布的统一软件包",
            ha="center", va="center", fontsize=10, color=GRAY)

    # Shared case definitions feed two separate implementations.  No arrow
    # runs from the KDK solver into the independent DOP853 implementation.
    top = 3.64
    height = 1.13
    card(ax, .18, top, 1.78, height, "固定实验协议", [
        "参数、初值、种子", "比较指标、输入哈希",
    ])
    card(ax, 2.32, top, 2.18, height, "输入与耦合状态", [
        r"$M^2(R)$ 与 $\partial_R M^2$",
        r"$(q_\ell,\dot q_\ell,R,U)$；$U=\dot R$",
    ])
    card(ax, 4.89, top, 1.63, height, "主积分 KDK", [
        "共同推进场与时钟", r"输出状态轨迹 $Z(\tau)$",
    ])
    card(ax, 6.95, top, 1.87, height, "主结果诊断", [
        "模式、能量与能流", "有效域与控制组",
    ], color=TEAL, fill=PALE_TEAL)
    arrow(ax, (2.02, 4.21), (2.26, 4.21))
    arrow(ax, (4.56, 4.21), (4.83, 4.21))
    arrow(ax, (6.58, 4.21), (6.89, 4.21))

    lower = 1.92
    card(ax, 2.32, lower, 2.18, 1.04, "独立参考 DOP853", [
        "相同输入，分别实现右端", "独立步长与精度控制",
    ])
    card(ax, 4.89, lower, 1.63, 1.04, "主解—参考解", [
        "轨迹误差、收敛", "失败与敏感性记录",
    ], color=TEAL, fill=PALE_TEAL)
    card(ax, 6.95, lower, 1.87, 1.04, "归档与复现", [
        "结果、图表与检查", "版本、配置及哈希",
    ], color=TEAL, fill=PALE_TEAL)
    arrow(ax, (3.41, 3.58), (3.41, 3.02))
    ax.text(2.58, 3.27, "同一案例", fontsize=9.5, color=GRAY,
            ha="center", va="center")
    arrow(ax, (4.56, 2.44), (4.83, 2.44))
    arrow(ax, (5.70, 3.58), (5.70, 3.02), color=TEAL)
    arrow(ax, (6.58, 2.44), (6.89, 2.44), color=TEAL)
    arrow(ax, (7.88, 3.58), (7.88, 3.02), color=TEAL)

    ax.text(.27, 2.72, "验证分支", fontsize=10.7, color=BLUE, ha="left")
    ax.text(.27, 2.36, "独立实现检验数值误差", fontsize=9.5, color=GRAY,
            ha="left")
    ax.text(.27, 2.02, "物理对应另需观测检验", fontsize=9.5, color=GRAY,
            ha="left")

    # Future software adapters do not imply a physical quantum engine.
    card(ax, .18, .56, 4.15, .83, "待实现：ECS／Gym 适配层", [
        "状态封装、统一调度与交互接口",
    ], color=GOLD, fill=PALE_GOLD, dashed=True)
    card(ax, 4.67, .56, 4.15, .83, "待建立：物理读出", [
        "量纲、宇宙时间标定与真实相互作用",
    ], color=GOLD, fill=PALE_GOLD, dashed=True)
    ax.text(4.5, 1.61,
            r"$h$：数值步长；$\tau$：模型演化时间；$R$：动态内部坐标，尚非宇宙年龄。",
            ha="center", va="center", fontsize=9.7, color=INK)
    ax.text(4.5, .20, "实线：已有脚本的数据与检验流程　　虚线框：后续工程或物理接口",
            ha="center", va="center", fontsize=9.5, color=GRAY)

    for ext in ("pdf", "svg", "png"):
        path = OUT / f"chapter07_engine_pipeline_cn.{ext}"
        fig.savefig(path, dpi=240,
                    metadata={"Creator": Path(__file__).name} if ext == "pdf" else None)
        if ext == "svg":
            path.write_text("\n".join(line.rstrip() for line in path.read_text().splitlines()) + "\n")
    plt.close(fig)


if __name__ == "__main__":
    draw()
    print(f"Font: {FONT}")
    print(f"Saved Chapter 7 diagram as PDF/SVG/PNG in {OUT}")
