#!/usr/bin/env python3
"""Reproduce the two conceptual front-matter diagrams for the Chinese RSM draft.

These diagrams contain no newly calculated scientific data.  The limited
numerical finding shown in figure 2 is taken from the archived v0.5 study.
Run from any directory: python /path/to/make_overview_figures.py
Requires Python 3, Matplotlib, and a Chinese font (override with RSM_CN_FONT).
"""

from __future__ import annotations

import os
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch


OUT = Path(__file__).resolve().parent
FONT_CANDIDATES = [
    os.environ.get("RSM_CN_FONT", ""),
    "/usr/share/fonts/truetype/droid/DroidSansFallbackFull.ttf",
    "/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc",
    "/usr/share/fonts/truetype/arphic-gbsn00lp/gbsn00lp.ttf",
]
FONT_PATH = next((Path(p) for p in FONT_CANDIDATES if p and Path(p).is_file()), None)
if FONT_PATH is None:
    raise SystemExit("Chinese font missing: set RSM_CN_FONT to a CJK font file.")
font_manager.fontManager.addfont(str(FONT_PATH))
FONT = font_manager.FontProperties(fname=str(FONT_PATH)).get_name()
matplotlib.rcParams.update(
    {
        "font.family": [FONT, "DejaVu Sans"],
        "font.size": 10.5,
        "mathtext.fontset": "stix",
        "axes.unicode_minus": False,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "path",
        "figure.facecolor": "white",
        "savefig.facecolor": "white",
    }
)

INK = "#223344"
BLUE = "#416B89"
TEAL = "#26776B"
GOLD = "#A87935"
GRAY = "#6B7580"
PALE_BLUE = "#F0F5F9"
PALE_TEAL = "#EFF7F4"
PALE_GOLD = "#FBF6ED"


def canvas(title: str, subtitle: str):
    fig = plt.figure(figsize=(9.0, 6.2))
    ax = fig.add_axes([0.025, 0.035, 0.95, 0.94])
    ax.set_xlim(0, 9)
    ax.set_ylim(0, 6.2)
    ax.axis("off")
    ax.text(4.5, 6.06, title, ha="center", va="center", fontsize=14, color=INK)
    ax.text(4.5, 5.73, subtitle, ha="center", va="center", fontsize=9.6, color=GRAY)
    return fig, ax


def box(ax, x, y, w, h, title, lines=(), *, color=BLUE, fill=PALE_BLUE, dashed=False):
    patch = FancyBboxPatch(
        (x, y), w, h,
        boxstyle="round,pad=0.045,rounding_size=0.075",
        linewidth=1.05, edgecolor=color, facecolor=fill,
        linestyle=(0, (4, 3)) if dashed else "solid", zorder=3,
    )
    ax.add_patch(patch)
    ax.text(x + w / 2, y + h - 0.19, title, ha="center", va="center", fontsize=11, color=color, zorder=4)
    if lines:
        available = h - 0.38
        for i, line in enumerate(lines):
            ty = y + available * (1 - (i + 0.5) / len(lines))
            ax.text(x + w / 2, ty, line, ha="center", va="center", fontsize=10.1, color=INK, zorder=4)


def arrow(ax, start, end, *, color=BLUE, dashed=False, rad=0):
    ax.add_patch(FancyArrowPatch(
        start, end, arrowstyle="-|>", mutation_scale=11,
        linewidth=1.15, color=color,
        linestyle=(0, (4, 3)) if dashed else "solid",
        connectionstyle=f"arc3,rad={rad}", zorder=2,
    ))


def save(fig, name):
    for ext in ("pdf", "svg", "png"):
        fig.savefig(OUT / f"{name}.{ext}", dpi=240, metadata={"Creator": "make_overview_figures.py"} if ext == "pdf" else None)
    plt.close(fig)


def architecture():
    fig, ax = canvas(
        "RSM 总体构架：从动力学假设到可检验接口",
        "沿用原文构建顺序；底层模型、描述层与物理对应分别陈述",
    )
    box(ax, 0.18, 4.64, 3.64, 0.78, "算术与非自治动力学启发 · 第2章", ["零点谱统计；已发表的映射数值对应"])
    box(ax, 4.30, 4.64, 4.52, 0.78, "候选响应与内部演化 · 第2—4章", [r"规定 $1/\ln^2(R/R_*)$；内部坐标 $R$ 参与反馈"], color=GOLD, fill=PALE_GOLD, dashed=True)
    arrow(ax, (3.88, 5.03), (4.22, 5.03), color=GOLD, dashed=True)

    boxes = [(0.18, "离散状态", [r"格点自由度 $(q_\ell,p_\ell)$", "状态变量与物理实体需区分"]),
             (3.18, "同步更新 · 第3章", ["统一数值更新规则", "步长不直接等同普朗克时间"]),
             (6.18, "辛／Hamilton 演化 · 第4章", ["共同 Hamilton 量与正则方程", "能量交换、辛积分与误差核验"])]
    for x, title, lines in boxes:
        box(ax, x, 3.26, 2.64, 1.00, title, lines)
    arrow(ax, (6.56, 4.59), (7.50, 4.32))
    arrow(ax, (4.98, 4.59), (4.50, 4.32))
    arrow(ax, (4.72, 4.59), (1.50, 4.32))

    box(ax, 0.18, 1.93, 4.10, 0.93, "逻辑描述：状态关联 · 第3、6章", ["状态标识、关联关系与访问规则", "对应量子关联的机制尚待建立"], color=GOLD, fill=PALE_GOLD, dashed=True)
    box(ax, 4.72, 1.93, 4.10, 0.93, "空间描述：局部传播 · 第3—6章", ["近邻相互作用与长波连续描述", "对应真实时空与物质的机制尚待建立"], color=GOLD, fill=PALE_GOLD, dashed=True)
    arrow(ax, (1.5, 3.20), (2.23, 2.92), color=GOLD, dashed=True)
    arrow(ax, (7.5, 3.20), (6.77, 2.92), color=GOLD, dashed=True)

    box(ax, 0.18, 0.55, 5.20, 0.99, "已完成的模型基准 · 第4、8章", ["内部时钟—非线性格点双向反馈；长波连续近似", "真实算术残差与替代输入对照"], color=TEAL, fill=PALE_TEAL)
    box(ax, 5.83, 0.55, 2.99, 0.99, "ECS 实现方案 · 第7—8章", ["实体、组件、系统的工程分工", "拟支持模拟、复现与观测接口"])
    # Route around the descriptive hypotheses: the numerical baseline does not
    # depend on their still-unestablished interpretation as real physics.
    ax.plot([8.92, 8.92, 5.56, 5.56], [3.72, 1.72, 1.72, 1.26], color=TEAL, lw=1.0, zorder=1)
    ax.plot([8.85, 8.92], [3.72, 3.72], color=TEAL, lw=1.0, zorder=1)
    arrow(ax, (5.56, 1.26), (5.45, 1.26), color=TEAL)
    arrow(ax, (5.77, 0.90), (5.44, 0.90), color=BLUE)
    ax.text(4.5, 0.19, "蓝：模型定义／工程组织    绿：已计算基准    虚线金色：启发与待建立的物理对应", ha="center", va="center", fontsize=9.2, color=GRAY)
    save(fig, "rsm_architecture_cn")


def reasoning():
    fig, ax = canvas(
        "RSM 研究路径：从数值启发到物理检验",
        "逆对数平方是可检验的候选假设；谱统计相似本身不推出物理常数漂移",
    )
    box(ax, 0.18, 4.66, 4.10, 0.78, "黎曼零点谱统计背景", ["与随机矩阵统计的联系；Hilbert–Pólya 思路"])
    box(ax, 4.72, 4.66, 4.10, 0.78, "作者已发表的非自治映射研究", ["算术序列与映射的数值对应，尚非物理同构"])
    box(ax, 1.18, 3.46, 6.64, 0.82, "候选响应假设 · 第2—4章", [r"$M^2(R)$ 含 $1/\ln^2(R/R_*)$；索引 $j$ 与时间 $t$ 的对应需另建"], color=GOLD, fill=PALE_GOLD, dashed=True)
    arrow(ax, (2.23, 4.60), (3.48, 4.34), color=GOLD, dashed=True)
    arrow(ax, (6.77, 4.60), (5.52, 4.34), color=GOLD, dashed=True)

    box(ax, 0.18, 2.00, 4.10, 1.04, "给定模型后的推导 · 第4章", ["自主 Hamilton 系统；时钟—格点双向反馈", "条件晚时尾律与长波连续近似"], color=TEAL, fill=PALE_TEAL)
    box(ax, 4.72, 2.00, 4.10, 1.04, "算术输入对照 · 第8章", ["真实零点残差与固定替代集合比较", "八个真实块未检出特殊终点能量响应"], color=TEAL, fill=PALE_TEAL)
    arrow(ax, (3.35, 3.40), (2.23, 3.10), color=TEAL)
    arrow(ax, (5.65, 3.40), (6.77, 3.10), color=TEAL)

    box(ax, 0.18, 0.66, 4.10, 0.92, "待建立的物理桥接 · 第5—6章", ["真实物质耦合；内部演化与宇宙时间标定", "选择可观测量，给出误差与适用范围"], color=GOLD, fill=PALE_GOLD, dashed=True)
    box(ax, 4.72, 0.66, 4.10, 0.92, "观测检验与模型筛选 · 第5、8章", [r"原子钟／光谱的 $\alpha$；膨胀史的 $H(z)$", "先检验额外变化，再比较时间依赖形式"], color=GOLD, fill=PALE_GOLD, dashed=True)
    arrow(ax, (2.23, 1.94), (2.23, 1.64), color=GOLD, dashed=True)
    arrow(ax, (6.77, 1.94), (3.72, 1.64), color=GOLD, dashed=True)
    arrow(ax, (4.34, 1.12), (4.66, 1.12), color=GOLD, dashed=True)
    ax.text(4.5, 0.29, "实线：给定模型内已定义、推导或计算的连接；虚线：启发性假设与待建立的物理连接", ha="center", va="center", fontsize=9.2, color=GRAY)
    ax.text(4.5, 0.04, "“未检出”限于当前检验量与对照设计，不等于排除全部算术效应。", ha="center", va="center", fontsize=9.2, color=GRAY)
    save(fig, "rsm_reasoning_cn")


if __name__ == "__main__":
    architecture()
    reasoning()
    print(f"Font: {FONT} ({FONT_PATH})")
    print(f"Saved PDF/SVG/PNG (9.0 × 6.2 in, PNG 240 dpi) in {OUT}")
