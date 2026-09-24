#!/usr/bin/env python3
"""Redraw four archived CMB-conditioned lattice states and saved diagnostics.

Input is the small, locally frozen plotting snapshot. No neighboring project,
HEALPix, lattice integration, sky fitting, or cosmic-age assignment is used.
All heatmaps preserve the stored field values and share one linear colour
scale; the curve below uses full-volume variance, not the plotted slice.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.ticker import NullFormatter
import numpy as np

from make_overview_figures import BLUE, FONT, GOLD, GRAY, INK, PALE_GOLD, TEAL


OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[1]
SOURCE_DIR = ROOT / "reports/chapter04_timeline_sources"
SOURCE = SOURCE_DIR / "timeline_snapshot.npz"
COMMIT = "f42db51748003fc83be064dade9bb52de5eb8319"

plt.rcParams.update({
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.edgecolor": GRAY, "axes.labelcolor": INK,
    "text.color": INK, "xtick.color": INK, "ytick.color": INK,
    "axes.axisbelow": True, "svg.hashsalt": "rsm-chapter04-timeline-20260924",
})


def load_snapshot():
    manifest = json.loads((SOURCE_DIR / "manifest.json").read_text())
    expected = manifest["local_files"][SOURCE.name]
    if (SOURCE.stat().st_size != expected["bytes"]
            or hashlib.sha256(SOURCE.read_bytes()).hexdigest() != expected["sha256"]):
        raise ValueError("The plotting snapshot differs from its manifest.")
    with np.load(SOURCE, allow_pickle=False) as archive:
        data = {key: archive[key].copy() for key in archive.files}
    metadata = json.loads(str(data["metadata_json"]))
    if metadata != manifest["extraction"] or metadata["source_commit"] != COMMIT:
        raise ValueError("Snapshot provenance differs from the fixed source.")
    if (data["slices_z48"].shape != (4, 96, 96)
            or not np.array_equal(data["snapshot_steps"], [0, 45, 70, 95])
            or data["histories"].shape != (4, 96, 5)):
        raise ValueError("Unexpected saved state or history shape.")
    return data, metadata


def draw(data, metadata):
    fig = plt.figure(figsize=(9.4, 6.8))
    fig.suptitle("从CMB约束初态到非线性场域", y=.972, fontsize=14)
    fig.text(.5, .918, r"同一 $96^3$ 格点上的实际存档状态；前 $45$ 步相同，从 $46$ 步起延续不同后期规则。",
             ha="center", fontsize=10.0, color=GRAY)

    titles = ("拟合初态", "CMB约束锚点", "计算中间态", "后期场域")
    image_lefts = (.066, .282, .498, .714)
    norm = Normalize(*metadata["shared_colour_limits"])
    for i, (left, title, step) in enumerate(zip(image_lefts, titles, data["snapshot_steps"])):
        ax = fig.add_axes([left, .563, .19, .263])
        ax.imshow(data["slices_z48"][i], origin="lower", interpolation="nearest",
                  cmap="RdBu_r", norm=norm)
        ax.set_title(f"({chr(97+i)}) {title}\n" + rf"$n={step}$", fontsize=10.2, pad=7)
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)
        variance = data["histories"][0, step, 1]
        ax.set_xlabel(rf"体积方差 $={variance:.4g}$", fontsize=9.7, labelpad=6)
        if i < 3:
            fig.text(left + .203, .692, "→", fontsize=15, color=GRAY,
                     ha="center", va="center")
    cax = fig.add_axes([.926, .563, .012, .263])
    colour_bar = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap="RdBu_r"), cax=cax,
                             ticks=[-3, 0, 3])
    colour_bar.set_label(r"标量场 $q$", fontsize=10, labelpad=3)
    colour_bar.ax.tick_params(labelsize=9)
    fig.text(.49, .496, r"上排均为中央切片 $z=48$，共用原始 $q$ 单位与色标；方差取完整三维体积。",
             ha="center", fontsize=9.8, color=GRAY)

    ax = fig.add_axes([.095, .207, .85, .218])
    ax.set_title("(e) 从同一锚点延续：每一步的三维体积方差", fontsize=10.6, pad=9)
    steps, histories = data["steps"], data["histories"]
    ax.axvspan(45, 95, color=PALE_GOLD, alpha=.75, zorder=0)
    ax.plot(steps[:46], histories[0, :46, 1], color=INK, lw=1.7, zorder=5)
    lines = {}
    styles = {
        0: (GOLD, "-", 2., r"$p=2$：双阱与阻尼（$4.282$）"),
        1: (BLUE, "-", 1.6, r"$p=2$：去掉局部力（$0.0388$）"),
        2: (GRAY, (0, (4, 2)), 1.65, "仅晚期常系数：同双阱与阻尼（4.273）"),
        3: (TEAL, (0, (1, 1.8)), 1.8, r"$p=2$：双阱，无晚期阻尼（$3.181$）"),
    }
    for i, (color, linestyle, width, label) in styles.items():
        lines[i], = ax.plot(steps[45:], histories[i, 45:, 1], color=color,
                           ls=linestyle, lw=width, label=label, zorder=3+i)
    ax.axvline(45, color=GRAY, lw=.9, ls=(0, (2, 2)))
    ax.text(43, 5.25, r"$n=45$ 锚点", ha="right", color=GRAY, fontsize=9.3)
    ax.text(3, .15, "共同线性前段", color=INK, fontsize=9.3)
    ax.set_yscale("log")
    ax.set_xlim(0, 95)
    ax.set_ylim(.026, 8)
    ax.set_xticks([0, 15, 30, 45, 60, 75, 95])
    ax.set_yticks([.03, .1, .3, 1, 3], ["0.03", "0.1", "0.3", "1", "3"])
    ax.yaxis.set_minor_formatter(NullFormatter())
    ax.grid(axis="y", which="major", color="#DCE2E6", lw=.55)
    ax.set_xlabel(r"迭代步 $n$（无量纲）", fontsize=10.1, labelpad=3)
    ax.set_ylabel(r"$\operatorname{Var}(q)$", fontsize=11)
    ax.legend(handles=[lines[i] for i in (0, 2, 1, 3)],
              loc="upper center", bbox_to_anchor=(.5, -.305), ncol=2,
              frameon=False, fontsize=9.0, handlelength=2.6, columnspacing=2.1,
              borderaxespad=0., labelspacing=.45)
    fig.text(.5, .02,
             "图例括号为末态方差。局部势与阻尼为规定输入；模型步数尚未标定宇宙年龄。",
             ha="center", fontsize=9.6, color=GRAY)

    for extension in ("pdf", "svg", "png"):
        path = OUT / f"chapter04_timeline_cn.{extension}"
        file_metadata = None
        if extension == "pdf":
            file_metadata = {"Creator": Path(__file__).name, "CreationDate": None, "ModDate": None}
        elif extension == "svg":
            file_metadata = {"Creator": Path(__file__).name, "Date": None}
        fig.savefig(path, dpi=240, metadata=file_metadata)
        if extension == "svg":
            path.write_text("\n".join(line.rstrip() for line in path.read_text().splitlines()) + "\n")
    plt.close(fig)


if __name__ == "__main__":
    data, metadata = load_snapshot()
    draw(data, metadata)
    print(f"Font: {FONT}")
    print(f"Input: {SOURCE.relative_to(ROOT)}")
    print(f"Input SHA256: {hashlib.sha256(SOURCE.read_bytes()).hexdigest()}")
    print(f"Saved chapter04_timeline_cn.pdf/png/svg in {OUT}")
