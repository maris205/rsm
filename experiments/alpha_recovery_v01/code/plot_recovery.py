#!/usr/bin/env python3
"""Plot the existing synthetic-recovery results; do not rerun or alter fits.

Run from any directory. Inputs are ../results/{summary.json,recovery.csv,
shape_geometry.csv}; outputs are ../figures/*.png and *.pdf.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator, NullFormatter, PercentFormatter
import numpy as np
from scipy.stats import norm


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
FIGURES = ROOT / "figures"
BLUE, ORANGE, GREEN, GREY = "#0072B2", "#D55E00", "#009E73", "#666666"
MAIN = "king_clock_offsets"


def read_csv(name):
    with (RESULTS / name).open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def select(rows, *, design, true_s, d):
    selected = [
        row for row in rows
        if row["design"] == design
        and int(row["true_s"]) == true_s
        and float(row["d_injected"]) == d
    ]
    assert len(selected) == 1, (design, true_s, d, len(selected))
    return selected[0]


def finish(fig, stem):
    FIGURES.mkdir(parents=True, exist_ok=True)
    for extension in ("png", "pdf"):
        path = FIGURES / f"{stem}.{extension}"
        fig.savefig(path, dpi=220, facecolor="white", bbox_inches="tight")
        print(path.relative_to(ROOT))
    plt.close(fig)


def recovery_capability(summary, recovery, designs):
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.9),
                             gridspec_kw={"width_ratios": [1.08, 1]})
    fig.subplots_adjust(left=0.075, right=0.975, top=0.90,
                        bottom=0.27, wspace=0.57)
    ax, right = axes
    sigma_clock = summary["current_clock_sigma_d"]
    sigma_joint = designs[MAIN]["sigma_d"]["2"]
    threshold = np.sqrt(designs[MAIN]["fixed_s2_5pct_threshold"])
    x = np.linspace(0, 5.2, 300)
    mu = x * sigma_clock / sigma_joint
    analytic = norm.sf(threshold - mu) + norm.cdf(-threshold - mu)
    ax.plot(x, analytic, color=GREY, lw=1.7, label="Gaussian analytic power")
    for sign, color, marker, label in (
        (1, BLUE, "o", r"Simulation, $d>0$"),
        (-1, ORANGE, "s", r"Simulation, $d<0$"),
    ):
        records = [select(recovery, design=MAIN, true_s=2, d=sign*d)
                   for d in (2.5, 5, 12.5)]
        xx = np.array([abs(float(r["d_injected"])) for r in records]) / sigma_clock
        yy = np.array([float(r["fixed_s2_detection"]) for r in records])
        err = np.array([float(r["detection_mc_se"]) for r in records])
        # A small horizontal displacement makes the two injection signs visible.
        ax.errorbar(xx + sign*0.025, yy, yerr=1.96*err, color=color,
                    marker=marker, ms=5, linestyle="none", capsize=3,
                    label=label, zorder=4)
    null = select(recovery, design=MAIN, true_s=0, d=0)
    ax.errorbar(0, float(null["fixed_s2_detection"]),
                yerr=1.96*float(null["detection_mc_se"]),
                color="black", marker="D", ms=4.5, capsize=3,
                label="Simulation, null", zorder=5)
    ax.axhline(0.05, color=GREY, linestyle=":", lw=1)
    ax.text(3.0, 0.079, "Nominal 5% false-positive rate", color=GREY, fontsize=8.5)
    ax.set(title=r"(a) Detecting amplitude, fixed $s=2$",
           xlabel=r"Injected $|d|/2.5$ (clock standard errors)",
           ylabel="Fraction rejecting zero drift",
           xlim=(-0.12, 5.25), ylim=(0, 1.055), xticks=range(6))
    ax.yaxis.set_major_formatter(PercentFormatter(1, decimals=0))
    ax.legend(loc="upper left", frameon=False, fontsize=8.7)
    ax.grid(axis="y", color="#dddddd", lw=0.7)

    values = [designs["king_only_offsets"]["sigma_D0_per_year"],
              designs[MAIN]["sigma_D0_per_year"], sigma_clock*1e-19]
    positions = [2, 1, 0]
    colors = [ORANGE, BLUE, GREY]
    for yy, value, color in zip(positions, values, colors):
        right.plot(value, yy, marker="o", color=color, ms=8, zorder=3)
        exponent = int(np.floor(np.log10(value)))
        mantissa = value / 10**exponent
        right.annotate(rf"${mantissa:.2f}\times10^{{{exponent}}}$", (value, yy),
                       xytext=(0, 12), textcoords="offset points", ha="center",
                       fontsize=10, color=color)
    right.set_xscale("log")
    right.set(title=r"(b) Drift precision, fixed $s=2$",
              xlabel=r"Expected $1\sigma$ error on $D_0$ (yr$^{-1}$)",
              xlim=(1.05e-19, 5e-16), ylim=(-0.8, 2.7),
              yticks=positions,
              yticklabels=["King + offsets", "King + offsets\n+ clock", "Clock alone"])
    right.xaxis.set_major_locator(LogLocator(base=10, numticks=5))
    right.xaxis.set_minor_formatter(NullFormatter())
    right.grid(axis="x", color="#dddddd", lw=0.7)
    information = designs[MAIN]["clock_information_fraction"]*100
    right.text(0.97, 0.03, f"Clock fraction of joint information: {information:.4f}%",
               transform=right.transAxes, ha="right", va="bottom", fontsize=8.8)
    fig.text(0.075, 0.14,
             "Synthetic data at the current observing design; 20,000 repetitions per injection. "
             "King: 293 measurements, with two instrument offsets.", fontsize=9)
    fig.text(0.075, 0.095,
             "Bars in (a): approximate 95% Monte Carlo uncertainty; signs are slightly displaced horizontally. "
             "Measured means are not used.", fontsize=9)
    fig.text(0.075, 0.05,
             r"$d=D_0/(10^{-19}\,{\rm yr}^{-1})$. Detecting drift does not establish its time dependence.",
             fontsize=9, color=GREY)
    finish(fig, "recovery_capability")


def shape_identifiability(summary, recovery, geometry):
    fig, axes = plt.subplots(1, 2, figsize=(12.0, 5.6),
                             gridspec_kw={"width_ratios": [0.9, 1.1]})
    fig.subplots_adjust(left=0.075, right=0.975, top=0.87,
                        bottom=0.32, wspace=0.37)
    ax, right = axes
    records = [select(recovery, design=MAIN, true_s=s, d=12.5) for s in (1, 2, 3)]
    keys = ("aic_constant", "aic_s1", "aic_s2", "aic_s3")
    values = np.array([[float(r[k]) for k in keys] for r in records])*100
    assert np.allclose(values.sum(axis=1), 100)
    ax.imshow(values, cmap="viridis", vmin=0, vmax=55, aspect="auto")
    for (iy, ix), value in np.ndenumerate(values):
        ax.text(ix, iy, f"{value:.3f}%", ha="center", va="center",
                fontsize=11, color="white" if value < 30 else "#111111")
    ax.set(title="(a) AIC selection frequencies",
           xticks=range(4), xticklabels=["Constant", r"$s=1$", r"$s=2$", r"$s=3$"],
           yticks=range(3), yticklabels=[r"$s=1$", r"$s=2$", r"$s=3$"],
           xlabel="Selected candidate", ylabel="Injected exponent")
    ax.set_xticks(np.arange(-0.5, 4, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, 3, 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=1.5)
    ax.tick_params(which="minor", bottom=False, left=False)
    ax.text(0.5, 1.04, r"$d=+12.5$: five clock standard errors",
            transform=ax.transAxes, ha="center", fontsize=9, color=GREY)

    curve = summary["ideal_precision_curve"]
    xx = np.array([r["error_reduction"] for r in curve])
    for field, d, color, style in (
        ("min_delta_chi2_d2p5", 2.5, BLUE, "-"),
        ("min_delta_chi2_d12p5", 12.5, ORANGE, "--"),
    ):
        yy = np.array([r[field] for r in curve])
        right.loglog(xx, yy, color=color, ls=style, lw=2,
                     label=rf"Injected $s=2$, $d={d:g}$")
        # Check that the current-precision curve starts at the CSV geometry.
        expected = min(float(r[f"delta_chi2_at_d{str(d).replace('.', 'p')}"])
                       for r in geometry if r["design"] == MAIN and r["true_s"] == "2")
        assert np.isclose(yy[0], expected, rtol=1e-8, atol=0)
    for target in (1, 9):
        right.axhline(target, color=GREY, ls=":" if target == 1 else "-.", lw=1)
        right.text(1.6, target*1.3, rf"$\Delta\chi^2={target}$", fontsize=9, color=GREY)
    for row in summary["ideal_precision_budget"]:
        color = BLUE if row["d_injected"] == 2.5 else ORANGE
        right.plot(row["astronomical_error_reduction"], row["achieved_min_delta_chi2"],
                   "o", color=color, ms=4)
    right.set(title="(b) Ideal precision needed for shape",
              xlabel="Reduction factor in astronomical standard errors",
              ylabel=r"Nearest alternative: noise-free $\Delta\chi^2$",
              xlim=(1, 1e7), ylim=(1e-11, 1e2))
    right.xaxis.set_major_locator(LogLocator(base=10, numticks=8))
    right.yaxis.set_major_locator(LogLocator(base=10, numticks=8))
    right.grid(which="major", color="#dddddd", lw=0.7)
    right.legend(loc="lower right", frameon=False, fontsize=9)
    fig.text(0.075, 0.21,
             "King + clock, with Keck/VLT offsets and a fixed Planck-time reference scale. "
             "Panel (a): 20,000 synthetic repetitions per row.", fontsize=9)
    fig.text(0.075, 0.165,
             "Near-collinear histories; winner frequencies are not evidence probabilities. "
             "A selected candidate is not an amplitude detection.", fontsize=9)
    fig.text(0.075, 0.12,
             "Panel (b): refit amplitude and offsets for each alternative (s = 1 or 3); "
             "clock precision is unchanged. All astronomical errors scale ideally.", fontsize=9)
    fig.text(0.075, 0.075,
             r"The reference levels 1 and 9 are residual-size benchmarks, not discovery significances. "
             "This is not an instrument forecast.", fontsize=9, color=GREY)
    finish(fig, "shape_identifiability")


def main():
    with (RESULTS / "summary.json").open(encoding="utf-8") as stream:
        summary = json.load(stream)
    assert summary["all_checks_passed"]
    assert summary["n_evaluation"] == 20000
    recovery = read_csv("recovery.csv")
    geometry = read_csv("shape_geometry.csv")
    designs = {row["name"]: row for row in summary["designs"]}
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 10,
        "axes.titlesize": 11, "axes.titlepad": 18,
        "axes.spines.top": False, "axes.spines.right": False,
        "pdf.fonttype": 42, "ps.fonttype": 42,
        "savefig.pad_inches": 0.12,
    })
    recovery_capability(summary, recovery, designs)
    shape_identifiability(summary, recovery, geometry)


if __name__ == "__main__":
    main()
