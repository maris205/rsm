#!/usr/bin/env python3
"""Plot Chapter 4's deterministic wave-transfer and long-wave diagnostics.

Read only the saved curves produced by
experiments/macro_laws_v01/code/verify_wave_transfer.py. This script performs
no evolution, fit, data selection by fit quality, or observational comparison.
The constant and inverse-log-square histories have the same cumulative
long-wave propagation proxy; their exact modal phases need not coincide.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
import numpy as np

from make_overview_figures import BLUE, FONT, GOLD, GRAY, INK, TEAL


OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[1]
SOURCE = ROOT / "experiments/macro_laws_v01/results/wave_transfer_curves.npz"

plt.rcParams.update({
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.edgecolor": GRAY,
    "axes.labelcolor": INK,
    "text.color": INK,
    "xtick.color": INK,
    "ytick.color": INK,
    "axes.axisbelow": True,
    "svg.hashsalt": "rsm-chapter04-macro-wave-20260924",
})


def load_curves():
    """Read checked curves without recalculating their generating model."""
    with np.load(SOURCE, allow_pickle=False) as archive:
        data = {key: archive[key].copy() for key in archive.files}
    metadata = json.loads(str(data["metadata_json"]))
    if (metadata["steps"] != 45 or metadata["exponents"] != [0, 2]
            or metadata["initial_pair"] != "q[-1] = q[0]"):
        raise ValueError("Saved model differs from this figure's stated design.")
    if not np.allclose(data["propagation_distance_sites"],
                       data["propagation_distance_sites"][0], rtol=1e-13, atol=0):
        raise ValueError("The two cumulative propagation proxies do not match.")
    return data


def draw(data):
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.9))
    fig.subplots_adjust(left=.083, right=.982, bottom=.255, top=.78, wspace=.36)
    fig.suptitle("从格点迭代到波动统计：可计算的尺度规律", y=.967, fontsize=14)
    fig.text(.5, .879, "相同的长波传播尺度，允许不同的精确传递；长波近似另有方向与尺度边界。",
             ha="center", fontsize=10, color=GRAY)

    ax, bx = axes
    lam = data["lambda_grid"]
    p0 = data["transfer_squared_p0"]
    p2 = data["transfer_squared_p2"]
    ax.plot(lam, p0, color=BLUE, lw=1.65, label=r"常系数 $p=0$")
    ax.plot(lam, p2, color=GOLD, lw=1.65, ls=(0, (4, 2)), label=r"逆对数平方 $p=2$")
    ax.set_title(r"$\mathrm{(a)}$ 两种历史均产生峰谷", fontsize=10.6, pad=11)
    ax.set_xlabel(r"格点本征值 $\lambda$", fontsize=10.5)
    ax.set_ylabel(r"模态功率传递 $|T_{45}(\lambda)|^2$", fontsize=10.5)
    ax.set_xlim(0, 2)
    ax.set_ylim(0, 1.35 * max(float(p0.max()), float(p2.max())))
    ax.set_xticks([0, .5, 1, 1.5, 2])
    ax.grid(axis="y", color="#DCE2E6", lw=.6)
    ax.legend(loc="upper left", frameon=False, fontsize=9.4, ncol=2,
              handlelength=2.1, columnspacing=1.3)

    k = data["k_grid"]
    use = (k >= .05) & (k <= np.pi)
    for direction, color, label in (("axis", BLUE, "轴向"), ("diagonal", TEAL, "体对角")):
        exact = data[f"lambda_{direction}_exact"]
        for order, key, ls, text in ((2, "lambda_second", (0, (4, 2)), "二阶"),
                                     (4, f"lambda_{direction}_fourth", "solid", "四阶")):
            error = np.abs(data[key][use] / exact[use] - 1)
            if not np.all(np.isfinite(error)) or not np.all(error > 0):
                raise ValueError(f"Invalid saved dispersion curve: {direction}, order {order}")
            bx.loglog(k[use], error, color=color, ls=ls, lw=1.6, label=f"{text} · {label}")
    bx.set_title(r"$\mathrm{(b)}$ 长波色散的截断误差", fontsize=10.6, pad=11)
    bx.set_xlabel(r"径向波数 $k$（格点单位）", fontsize=10.5)
    bx.set_ylabel(r"$|\lambda_{\rm approx}-\lambda|/\lambda$", fontsize=11.5)
    bx.set_xlim(.05, np.pi)
    bx.set_ylim(1e-9, 3)
    bx.set_xticks([.05, .1, .3, 1, 3], ["0.05", "0.1", "0.3", "1", "3"])
    bx.set_yticks([10.**power for power in (-9, -7, -5, -3, -1, 0)],
                 [rf"$10^{{{power}}}$" for power in (-9, -7, -5, -3, -1, 0)])
    bx.xaxis.set_minor_formatter(NullFormatter())
    bx.yaxis.set_minor_formatter(NullFormatter())
    bx.grid(which="major", color="#DCE2E6", lw=.6)
    bx.legend(loc="lower right", frameon=False, fontsize=9.0, ncol=2,
              handlelength=2.0, columnspacing=1.1)

    fig.text(.083, .123,
             r"$N=45$；两种系数匹配 $\sum_n\sqrt{c_n^2}$。右图：$\lambda_2=k^2/6$，"
             r"$\lambda_4=k^2/6-\sum_i k_i^4/72$。",
             fontsize=10, color=INK)
    fig.text(.083, .052,
             "仅验证规定模型的传递与长波展开，不含观测拟合；峰谷本身不构成对平方对数律的支持。",
             fontsize=9.7, color=GRAY)

    for extension in ("pdf", "svg", "png"):
        path = OUT / f"chapter04_wave_spectrum_cn.{extension}"
        metadata = None
        if extension == "pdf":
            metadata = {"Creator": Path(__file__).name, "CreationDate": None, "ModDate": None}
        elif extension == "svg":
            metadata = {"Creator": Path(__file__).name, "Date": None}
        fig.savefig(path, dpi=240, metadata=metadata)
        if extension == "svg":
            path.write_text("\n".join(line.rstrip() for line in path.read_text().splitlines()) + "\n")
    plt.close(fig)


if __name__ == "__main__":
    draw(load_curves())
    print(f"Font: {FONT}")
    print(f"Input: {SOURCE.relative_to(ROOT)}")
    print(f"Input SHA256: {hashlib.sha256(SOURCE.read_bytes()).hexdigest()}")
    print(f"Saved chapter04_wave_spectrum_cn.pdf/png/svg in {OUT}")
