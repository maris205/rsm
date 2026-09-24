#!/usr/bin/env python3
"""Draw saved physical-background and matched-PM bridge results.

No solver, CAMB calculation, initial-condition generation or fit is executed.
All six main panels use existing tables/snapshots. Density maps share one
log-density colour scale; zero projections are masked rather than repaired.
The companion spectrum figure retains the actual small paired responses.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.ticker import FuncFormatter, NullFormatter, ScalarFormatter
import numpy as np

from make_overview_figures import BLUE, FONT, GOLD, GRAY, INK, PALE_BLUE, TEAL


OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[1]
EXPERIMENT = ROOT / "experiments/cosmic_bridge_v01"
RESULTS = EXPERIMENT / "results"
DIAGNOSTICS = EXPERIMENT / "figures"
PROTOCOL_SHA = "74a2cab351d53f98169b885542d23d18b4d425c6fd6b75b91b49f5c53422e2e1"
PAIRS = ((64, 256), (64, 128), (64, 512), (32, 256))
plt.rcParams.update({
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.edgecolor": GRAY, "axes.labelcolor": INK, "text.color": INK,
    "xtick.color": INK, "ytick.color": INK, "axes.axisbelow": True,
    "svg.hashsalt": "rsm-chapter04-cosmic-bridge-20260924",
})


def sha(path):
    hasher = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024*1024), b""):
            hasher.update(block)
    return hasher.hexdigest()


def load_npz(path):
    with np.load(path, allow_pickle=False) as archive:
        return {key: archive[key].copy() for key in archive.files}


def load_inputs():
    paths = [EXPERIMENT / "protocol.md", RESULTS / "background_eps0.csv",
             RESULTS / "background_eps1e-4.csv", RESULTS / "background_summary.json",
             RESULTS / "clustering_summary.json",
             EXPERIMENT / "inputs/reference_linear_spectrum.npz",
             EXPERIMENT / "inputs/reference_spectrum_metadata.json"]
    paths += [RESULTS / f"{label}_n{n}_s{steps}.npz"
              for n, steps in PAIRS for label in ("reference", "clock")]
    for path in paths:
        if not path.is_file():
            raise FileNotFoundError(f"Required saved input is not ready: {path}")
    hashes = {str(path.relative_to(ROOT)): {"sha256": sha(path), "bytes": path.stat().st_size}
              for path in paths}
    if sha(EXPERIMENT / "protocol.md") != PROTOCOL_SHA:
        raise ValueError("Protocol changed from the frozen figure design.")
    summaries = {name: json.loads((RESULTS / f"{name}_summary.json").read_text())
                 for name in ("background", "clustering")}
    reference_meta = json.loads((EXPERIMENT / "inputs/reference_spectrum_metadata.json").read_text())
    if (summaries["background"]["metadata"]["protocol_sha256"] != PROTOCOL_SHA
            or summaries["clustering"]["protocol_sha256"] != PROTOCOL_SHA
            or reference_meta["protocol_sha256"] != PROTOCOL_SHA):
        raise ValueError("Saved inputs do not use the same frozen protocol.")
    backgrounds = {name: np.genfromtxt(RESULTS / f"background_{suffix}.csv", delimiter=",", names=True)
                   for name, suffix in (("reference", "eps0"), ("clock", "eps1e-4"))}
    if not np.array_equal(backgrounds["reference"]["a"], backgrounds["clock"]["a"]):
        raise ValueError("Matched growth tables have different scale-factor axes.")
    runs = {}
    for n, steps in PAIRS:
        for label in ("reference", "clock"):
            name = f"{label}_n{n}_s{steps}"
            path = RESULTS / f"{name}.npz"
            if sha(path) != summaries["clustering"]["runs"][name]["npz_sha256"]:
                raise ValueError(f"Snapshot is not the version bound by the summary: {name}")
            runs[name] = load_npz(path)
        if not np.array_equal(runs[f"reference_n{n}_s{steps}"]["k"],
                              runs[f"clock_n{n}_s{steps}"]["k"]):
            raise ValueError("A paired spectrum uses different k bins.")
    reference_path = EXPERIMENT / "inputs/reference_linear_spectrum.npz"
    if sha(reference_path) != reference_meta["output_files"][reference_path.name]["sha256"]:
        raise ValueError("External reference spectrum hash differs from its metadata.")
    return backgrounds, runs, summaries, load_npz(reference_path), hashes


def save_figure(fig, stem, extensions):
    stem.parent.mkdir(parents=True, exist_ok=True)
    for ext in extensions:
        target = stem.with_suffix("." + ext)
        metadata = None
        if ext == "pdf":
            metadata = {"Creator": Path(__file__).name, "CreationDate": None, "ModDate": None}
        elif ext == "svg":
            metadata = {"Creator": Path(__file__).name, "Date": None}
        fig.savefig(target, dpi=240, metadata=metadata)
        if ext == "svg":
            target.write_text("\n".join(line.rstrip() for line in target.read_text().splitlines()) + "\n")
    plt.close(fig)


def log_ticks(axis, which, ticks, labels):
    if which == "x":
        axis.set_xticks(ticks, labels)
        axis.xaxis.set_minor_formatter(NullFormatter())
    else:
        axis.set_yticks(ticks, labels)
        axis.yaxis.set_minor_formatter(NullFormatter())


def main_figure(backgrounds, runs, summaries):
    reference, clock = backgrounds["reference"], backgrounds["clock"]
    source = runs["clock_n64_s256"]
    box = float(summaries["clustering"]["box_mpc_over_h"])
    fig = plt.figure(figsize=(10., 7.2))
    fig.suptitle("带物理时间的最小连接：背景、线性增长与无碰撞聚类", y=.975, fontsize=14)
    fig.text(.5, .931, r"外部标准线性谱与共同粒子初态；主例 $\varepsilon=10^{-4}$，$64^3$ 网格，$256$ 步。",
             ha="center", fontsize=10, color=GRAY)
    lefts = (.075, .385, .695)
    axes = [fig.add_axes([left, .64, .23, .225]) for left in lefts]
    ax, bx, cx = axes
    ax.plot(reference["t_Gyr"], reference["a"], color=BLUE, lw=1.7,
            label=r"$\varepsilon=0$")
    ax.plot(clock["t_Gyr"], clock["a"], color=GOLD, lw=1.7, ls=(0, (4, 2)),
            label=r"$\varepsilon=10^{-4}$")
    ax.set_title("(a) 物理时间与尺度因子", fontsize=10.5, pad=9)
    ax.set_xlabel(r"宇宙年龄 $t\;(\mathrm{Gyr})$", fontsize=10)
    ax.set_ylabel(r"$a_{\rm cos}$", fontsize=11)
    ax.set_xlim(0, max(reference["t_Gyr"].max(), clock["t_Gyr"].max()))
    ax.set_ylim(0, 1.04)
    ax.legend(frameon=False, fontsize=9.0, loc="lower right")

    response = (clock["chi"][-1] / clock["chi"])**2
    bx.semilogx(clock["t_Gyr"], response, color=GOLD, lw=1.8)
    bx.set_title("(b) 指定的逆对数响应读数", fontsize=10.5, pad=9)
    bx.set_xlabel(r"宇宙年龄 $t\;(\mathrm{Gyr})$", fontsize=10)
    bx.set_ylabel(r"$\chi_0^2/\chi^2(t)$", fontsize=11)
    log_ticks(bx, "x", [.05, .5, 5], ["0.05", "0.5", "5"])
    bx.set_xlim(clock["t_Gyr"][0], clock["t_Gyr"][-1])
    bx.axhline(1, color=GRAY, lw=.65, ls=":")

    matched = clock["D_matched"] / reference["D_matched"] - 1
    cx.plot(clock["a"], matched, color=TEAL, lw=1.8)
    cx.axhline(0, color=GRAY, lw=.65)
    cx.set_title("(c) 同初动量的线性增长差", fontsize=10.5, pad=9)
    cx.set_xlabel(r"尺度因子 $a_{\rm cos}$", fontsize=10)
    cx.set_ylabel(r"$D_{\rm clock}/D_{\rm ref}-1$", fontsize=11)
    cx.set_xlim(.02, 1)
    small_response_ticks = ScalarFormatter(useMathText=False)
    small_response_ticks.set_scientific(True)
    small_response_ticks.set_powerlimits((-4, -4))
    cx.yaxis.set_major_formatter(small_response_ticks)
    cx.yaxis.get_offset_text().set_fontfamily("DejaVu Sans")
    cx.yaxis.get_offset_text().set_fontsize(9)
    for axis in axes:
        axis.grid(axis="y", color="#DCE2E6", lw=.55)

    maps, snapshot_info = [], []
    for value in (.1, .5, 1.):
        where = np.flatnonzero(np.isclose(source["a"], value, rtol=0, atol=1e-12))
        if len(where) != 1:
            raise ValueError(f"Missing stored density snapshot at a={value}")
        projection = np.asarray(source["projection"][where[0]], dtype=float)
        if not np.isfinite(projection).all() or np.any(projection < 0):
            raise ValueError("Projected density must be finite and nonnegative.")
        mean = float(np.mean(projection))
        if mean <= 0:
            raise ValueError("Projected density has zero mean.")
        shown = np.ma.log10(np.ma.masked_less_equal(projection / mean, 0))
        age = float(np.interp(np.log(value), clock["N"], clock["t_Gyr"]))
        maps.append(shown)
        snapshot_info.append({"a": value, "z": 1/value-1, "age_Gyr": age,
                              "projection_mean_before_normalization": mean,
                              "zero_pixels_masked": int(np.count_nonzero(projection == 0)),
                              "log_projection_min": float(shown.min()), "log_projection_max": float(shown.max())})
    maximum = max(float(np.ma.max(np.ma.abs(shown))) for shown in maps)
    if maximum <= 0:
        raise ValueError("All projected fields are uniform; no nonzero colour range is available.")
    norm = Normalize(-maximum, maximum)
    cmap = plt.colormaps["RdBu_r"].copy()
    cmap.set_bad("#D6D6D6")
    for i, (left, shown, info) in enumerate(zip(lefts, maps, snapshot_info)):
        axis = fig.add_axes([left, .18, .23, .32])
        axis.imshow(shown.T, origin="lower", interpolation="nearest", cmap=cmap, norm=norm,
                    extent=(0, box, 0, box))
        axis.set_title(rf"$\mathrm{{({chr(100+i)})}}\ a={info['a']:g},\ z={info['z']:g}$"
                       + "\n" + rf"$t={info['age_Gyr']:.3f}\ \mathrm{{Gyr}}$", fontsize=10.5, pad=8)
        axis.set_xticks([0, box/2, box])
        axis.set_yticks([0, box/2, box])
        axis.set_xlabel(r"$x\;(\mathrm{Mpc}/h)$", fontsize=10, labelpad=3)
        if i == 0:
            axis.set_ylabel(r"$y\;(\mathrm{Mpc}/h)$", fontsize=10)
    cax = fig.add_axes([.255, .075, .49, .015])
    colorbar = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax, orientation="horizontal")
    colorbar.set_label(r"共同色标：$\log_{10}[\Sigma/\langle\Sigma\rangle]$", fontsize=10, labelpad=3)
    colorbar.ax.xaxis.set_label_position("top")
    colorbar.ax.tick_params(labelsize=9)
    fig.text(.5, .012, "下排为无碰撞物质的全盒密度投影；未模拟气体、恒星形成或可识别星系。",
             ha="center", fontsize=9.7, color=GRAY)
    save_figure(fig, OUT / "chapter04_cosmic_bridge_cn", ("pdf", "png", "svg"))
    return {"main_density_snapshots": snapshot_info, "common_log_colour_limits": [-maximum, maximum],
            "response_initial_and_final": [float(response[0]), float(response[-1])],
            "matched_linear_growth_fractional_final": float(matched[-1]),
            "age_assignment": "linear interpolation in saved ln(a) table; no background re-integration",
            "projection_orientation": "array.T puts the original x-index horizontally and y-index vertically"}


def power_figure(runs, linear_reference):
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(10., 4.8))
    fig.subplots_adjust(left=.085, right=.975, bottom=.30, top=.80, wspace=.32)
    fig.suptitle("共同初态下的物质功率与微小配对响应", y=.967, fontsize=14)
    fig.text(.5, .881, "同一随机实现、初始位移与实际动量；原始网格密度估计，没有观测拟合。",
             ha="center", fontsize=10, color=GRAY)
    reference, clock = runs["reference_n64_s256"], runs["clock_n64_s256"]
    k = reference["k"]
    ax.loglog(k, reference["pk"][-1], color=BLUE, marker="o", ms=3, lw=1.7,
              label=r"参考背景，$64^3/256$")
    ax.loglog(k, clock["pk"][-1], color=GOLD, ls=(0, (4, 2)), lw=1.7,
              label=r"时钟背景，$64^3/256$")
    use = (linear_reference["k_h_mpc"] >= k.min()) & (linear_reference["k_h_mpc"] <= k.max())
    ax.loglog(linear_reference["k_h_mpc"][use], linear_reference["pk_z0_mpc_over_h3"][use],
              color=GRAY, ls=":", lw=1.5, label=r"外部 $\mathrm{CAMB}$ 线性参考，$z=0$")
    ax.set_title(r"$\mathrm{(a)}$ 末态 $a=1$ 的绝对功率", fontsize=10.6, pad=9)
    ax.set_xlabel(r"$k\;(h/\mathrm{Mpc})$", fontsize=10.5)
    ax.set_ylabel(r"$P(k)\;[(\mathrm{Mpc}/h)^3]$", fontsize=11)
    ax.legend(frameon=False, fontsize=9.0, loc="lower left")
    ax.yaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value:g}"))
    ax.yaxis.set_minor_formatter(NullFormatter())

    styles = ((GOLD, "-", "o"), (BLUE, (0, (4, 2)), "v"),
              (TEAL, (0, (1, 1.8)), "s"), (GRAY, (0, (6, 2)), "D"))
    response_info = {}
    for (n, steps), (color, linestyle, marker) in zip(PAIRS, styles):
        rr, cc = runs[f"reference_n{n}_s{steps}"], runs[f"clock_n{n}_s{steps}"]
        response = cc["pk"][-1]/rr["pk"][-1]-1
        if not np.isfinite(response).all():
            raise ValueError("A paired power response is nonfinite; report it before plotting.")
        bx.semilogx(rr["k"], response, color=color, ls=linestyle, marker=marker,
                    lw=1.4, ms=3.3, label=rf"${n}^3$，${steps}$ 步")
        response_info[f"n{n}_s{steps}"] = {"k_h_mpc": rr["k"].tolist(), "fractional_response": response.tolist()}
    bx.axhline(0, color=GRAY, lw=.7)
    bx.set_title(r"$\mathrm{(b)}$ 配对功率比：保留实际幅度", fontsize=10.6, pad=9)
    bx.set_xlabel(r"$k\;(h/\mathrm{Mpc})$", fontsize=10.5)
    bx.set_ylabel(r"$P_{\rm clock}/P_{\rm ref}-1$", fontsize=11)
    bx.yaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value:.1e}"))
    bx.legend(frameon=False, fontsize=9.0, loc="lower left")
    for axis in (ax, bx):
        axis.axvspan(k.min(), .3, color=PALE_BLUE, alpha=.8, zorder=0)
        axis.axvline(.3, color=GRAY, lw=.75, ls=":")
        axis.set_xlim(k.min()*.93, k.max()*1.05)
        log_ticks(axis, "x", [.15, .3, .6, 1, 2], ["0.15", "0.3", "0.6", "1", "2"])
        axis.grid(axis="y", color="#DCE2E6", lw=.55)
    fig.text(.085, .169, r"浅蓝：按频箱平均波数选择 $\bar{k}\leq0.3\ h/\mathrm{Mpc}$。网格窗口未去卷积，未减去粒子散粒噪声。",
             fontsize=9.7, color=INK)
    fig.text(.085, .104, "外部CAMB曲线未经过同一有限盒、带限和网格窗口；偏离不是观测误差或纯非线性效应。",
             fontsize=9.7, color=GRAY)
    fig.text(.085, .039, "分辨率差为数值诊断，不是观测误差条或检测显著性；一盒实现没有宇宙方差估计。",
             fontsize=9.7, color=GRAY)
    save_figure(fig, DIAGNOSTICS / "matter_power_diagnostics", ("pdf", "png"))
    return response_info


if __name__ == "__main__":
    backgrounds, runs, summaries, linear_reference, input_hashes = load_inputs()
    reading = main_figure(backgrounds, runs, summaries)
    responses = power_figure(runs, linear_reference)
    for relative, record in input_hashes.items():
        if sha(ROOT / relative) != record["sha256"]:
            raise RuntimeError(f"Input changed during rendering: {relative}")
    outputs = [OUT / f"chapter04_cosmic_bridge_cn.{suffix}" for suffix in ("pdf", "png", "svg")]
    outputs += [DIAGNOSTICS / f"matter_power_diagnostics.{suffix}" for suffix in ("pdf", "png")]
    manifest = {"scope": "Saved conditional background and PM results; no simulation or fit during plotting.",
                "protocol_sha256": PROTOCOL_SHA, "generator_sha256": sha(__file__), "font": FONT,
                "inputs": input_hashes, "readout": reading, "paired_power_responses": responses,
                "outputs": {str(path.relative_to(ROOT)): {"sha256": sha(path), "bytes": path.stat().st_size} for path in outputs}}
    (DIAGNOSTICS / "plot_inputs.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False, allow_nan=False) + "\n")
    print(f"Saved Chapter 4 cosmic bridge and matter-power diagnostics; {len(input_hashes)} input hashes recorded.")
