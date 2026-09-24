#!/usr/bin/env python3
"""Plot the saved controlled-collapse / conditional-cooling experiment.

This module does not import a simulation, ODE or cooling solver.  All six PM
runs are mandatory, including their explicit failure/missing-event records.
No missing PM curve is replaced by an illustration or extrapolation.  The
cooling grid is a prescribed-state diagnostic, not hydrodynamics or a galaxy
catalogue.  Run --check-inputs to check readiness without creating a figure.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from matplotlib.patches import Rectangle
from matplotlib.ticker import FuncFormatter, NullFormatter
import numpy as np

from make_overview_figures import BLUE, FONT, GOLD, GRAY, INK, TEAL


OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[1]
EXPERIMENT = ROOT / "experiments/halo_cooling_v01"
RESULTS = EXPERIMENT / "results"
PROTOCOL_SHA = "f3a0cc5abf21a4374439db756cbf8c4a3103bcabfa2bee6f9bd1921b28e608f2"
CASES = ("ref_force128", "ref_main", "ref_halfstep", "ref_particles128", "ref_shift", "clock_main")
CASE_LABELS = ("粗力网格", "参考主例", "半步长", "加粒子", "相位平移", "时钟主例")
BRANCHES = ("reference", "clock")
BRANCH_COLOURS = {"reference": BLUE, "clock": GOLD}
plt.rcParams.update({
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.edgecolor": GRAY, "axes.labelcolor": INK, "text.color": INK,
    "xtick.color": INK, "ytick.color": INK, "axes.axisbelow": True,
    "hatch.linewidth": .55,
})


def sha(path):
    hasher = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            hasher.update(block)
    return hasher.hexdigest()


def read_json(path):
    return json.loads(Path(path).read_text())


def read_csv(path):
    with Path(path).open(newline="") as stream:
        return list(csv.DictReader(stream))


def column(rows, name):
    return np.asarray([float(row[name]) for row in rows], dtype=float)


def truth(value):
    if isinstance(value, bool):
        return value
    if value in ("True", "true", "1"):
        return True
    if value in ("False", "false", "0"):
        return False
    raise ValueError(f"Unrecognized saved boolean: {value!r}")


def required_paths():
    names = ["sphere_summary.json", "sphere_reference.csv", "sphere_clock.csv",
             "cooling_grid.csv", "cooling_events.json", "cooling_summary.json",
             "spherical_pm_summary.json"]
    summary_path = RESULTS / "spherical_pm_summary.json"
    runs = read_json(summary_path).get("runs", {}) if summary_path.is_file() else {}
    for case in CASES:
        suffix = "_execution_failure" if runs.get(case, {}).get("status") == "execution_failed" else ""
        names.append(f"pm_{case}{suffix}.json")
    return [EXPERIMENT / "protocol.md"] + [RESULTS / name for name in names]


def load_inputs():
    paths = required_paths()
    missing = [str(path.relative_to(ROOT)) for path in paths if not path.is_file()]
    if missing:
        raise FileNotFoundError("Saved results are not all ready:\n" + "\n".join(missing))
    if sha(EXPERIMENT / "protocol.md") != PROTOCOL_SHA:
        raise ValueError("The frozen protocol changed.")
    sphere = read_json(RESULTS / "sphere_summary.json")
    cooling = read_json(RESULTS / "cooling_summary.json")
    events = read_json(RESULTS / "cooling_events.json")
    pm_summary = read_json(RESULTS / "spherical_pm_summary.json")
    if sphere["protocol_sha256"] != PROTOCOL_SHA or cooling["manifest"]["protocol_sha256"] != PROTOCOL_SHA:
        raise ValueError("ODE or cooling results do not use the frozen protocol.")
    if pm_summary["protocol_sha256"] != PROTOCOL_SHA:
        raise ValueError("PM results do not use the frozen protocol.")
    if pm_summary["input_sha256"]["sphere_summary.json"] != sha(RESULTS / "sphere_summary.json"):
        raise ValueError("PM calibration differs from the loaded sphere summary.")
    if events["status"] != "event_readouts_computed":
        raise ValueError("The two saved ODE cooling-event readouts are not ready.")
    if events["input_sha256"] != sha(RESULTS / "sphere_summary.json"):
        raise ValueError("Cooling-event input differs from the loaded sphere summary.")
    for name in ("cooling_grid.csv", "cooling_events.json"):
        if cooling["output_sha256"][name] != sha(RESULTS / name):
            raise ValueError(f"Cooling result changed after its summary: {name}")
    ode = {branch: read_csv(RESULTS / f"sphere_{branch}.csv") for branch in BRANCHES}
    pm = {}
    for case in CASES:
        saved = pm_summary["runs"].get(case)
        if saved is None or saved["status"] not in ("completed", "execution_failed"):
            raise ValueError(f"PM case is not completed or explicitly failed: {case}")
        suffix = "_execution_failure" if saved["status"] == "execution_failed" else ""
        pm[case] = read_json(RESULTS / f"pm_{case}{suffix}.json")
    for case, run in pm.items():
        if run["status"] == "execution_failed":
            if run != pm_summary["runs"][case] or run["protocol_sha256"] != PROTOCOL_SHA:
                raise ValueError(f"Explicit execution-failure record differs from the matrix summary: {case}")
            continue
        if run["metadata"]["protocol_sha256"] != PROTOCOL_SHA:
            raise ValueError(f"PM run uses a different protocol: {case}")
        if sha(RESULTS / f"pm_{case}.json") != pm_summary["runs"][case]["json_sha256"]:
            raise ValueError(f"Saved PM run differs from its completed matrix summary: {case}")
        if not isinstance(run.get("history"), list) or not run["history"]:
            raise ValueError(f"Saved PM history missing for {case}; do not substitute a curve.")
        if "event200" not in run:
            raise ValueError(f"Saved PM event status is absent for {case}.")
    grid = read_csv(RESULTS / "cooling_grid.csv")
    hashes = {str(path.relative_to(ROOT)): {"sha256": sha(path), "bytes": path.stat().st_size}
              for path in paths}
    return sphere, ode, pm, pm_summary, grid, events, cooling, hashes


def radius_panel(ax, sphere, ode, pm):
    h = float(sphere["inputs"]["h"])
    ql = float(sphere["inputs"]["R_L_Mpc_h"])
    readout = {}
    for branch, case in (("reference", "ref_main"), ("clock", "clock_main")):
        color = BRANCH_COLOURS[branch]
        event = sphere["models"][branch]["events"]["200"]
        rows = ode[branch]
        mask = column(rows, "a") < event["a"]
        age = np.r_[column(rows, "t_Gyr")[mask], event["t_Gyr"]]
        radius = np.r_[column(rows, "physical_radius_Mpc_h")[mask], event["physical_radius_Mpc_h"]] / h
        branch_label = "参考" if branch == "reference" else "时钟"
        ax.plot(age, radius, color=color, lw=1.75, ls="-" if branch == "reference" else (0, (5, 2)),
                label=branch_label + r" $\mathrm{ODE}$")
        if pm[case]["status"] == "completed":
            history = pm[case]["history"]
            pm_age = column(history, "t_Gyr")
            pm_radius = column(history, "a") * ql * column(history, "y_median") / h
            ax.plot(pm_age, pm_radius, color=color, lw=1.0, ls=":", marker="o" if branch == "reference" else "s",
                    ms=2.4, markevery=max(1, len(history)//17), fillstyle="none", label=branch_label + r" $\mathrm{PM}$")
            pm_record = {"pm_last_age_Gyr": float(pm_age[-1]), "pm_last_radius_proxy_Mpc": float(pm_radius[-1]),
                         "pm_history_points": len(history)}
        else:
            ax.text(.03, .85 if branch == "reference" else .77, branch_label + " PM 执行失败：无轨迹",
                    transform=ax.transAxes, fontsize=8.5, color="#A44838")
            pm_record = {"execution_failed": True, "failure_record": pm[case]}
        turn = sphere["models"][branch]["turnaround"]
        ax.scatter([turn["t_Gyr"]], [turn["physical_radius_Mpc_h"]/h], color=color, marker="v", s=22, zorder=5)
        ax.scatter([event["t_Gyr"]], [event["physical_radius_Mpc_h"]/h], color=color, marker="D", s=19, zorder=5)
        readout[branch] = {"ode_turnaround_age_Gyr": turn["t_Gyr"], "ode_Delta200_age_Gyr": event["t_Gyr"], **pm_record}
    ax.set_title("(a) 物理半径轨迹", fontsize=10.5, pad=9)
    ax.set_xlabel(r"宇宙年龄 $t\;(\mathrm{Gyr})$", fontsize=9.5)
    ax.set_ylabel(r"$R=a\,q_L y\;(\mathrm{Mpc})$", fontsize=10.5)
    ax.set_xticks([0, 2, 4, 6])
    ax.set_yticks([0, .2, .4, .6])
    ax.set_ylim(bottom=0)
    ax.scatter([], [], color=GRAY, marker="v", s=20, label="转向")
    ax.scatter([], [], color=GRAY, marker="D", s=18, label=r"$\Delta=200$")
    handles, labels = ax.get_legend_handles_labels()
    order = [0, 1, 4, 2, 3, 5]
    ax.legend([handles[i] for i in order], [labels[i] for i in order], frameon=False, fontsize=9.0,
              loc="lower center", ncol=2, columnspacing=.55, handlelength=1.3,
              borderpad=.15, labelspacing=.3, handletextpad=.4)
    ax.grid(axis="y", color="#DCE2E6", lw=.55)
    return readout


def error_panel(ax, sphere, pm, pm_summary):
    readout = {}
    event_errors, radius_errors = [], []
    for case in CASES:
        if pm[case]["status"] == "execution_failed":
            event_errors.append(None)
            radius_errors.append(None)
            readout[case] = {"execution_failed": True, "failure_record": pm[case]}
            continue
        branch = "clock" if case == "clock_main" else "reference"
        ode_event = sphere["models"][branch]["events"]["200"]["a"]
        event = pm[case]["event200"]
        history = pm[case]["history"]
        if event is None:
            stop = min(float(history[-1]["a"]), ode_event)
            event_error = None
        else:
            stop = min(float(event["a"]), ode_event)
            event_error = float(event["a"])/ode_event-1
        kept = [row for row in history if float(row["a"]) <= stop*(1+1e-13)
                and row.get("relative_y_error_vs_ode") is not None]
        if not kept:
            raise ValueError(f"No stored common pre-event trajectory errors for {case}")
        radius_error = max(abs(float(row["relative_y_error_vs_ode"])) for row in kept)
        diagnostics = pm[case].get("diagnostics", {})
        saved = diagnostics.get("max_radius_relative_error_vs_ode_common_pre_event")
        if saved is not None and not np.isclose(radius_error, saved, rtol=1e-10, atol=1e-13):
            raise ValueError(f"Stored/reconstructed common-interval radius errors disagree: {case}")
        event_errors.append(None if event_error is None else 100*abs(event_error))
        radius_errors.append(100*radius_error)
        readout[case] = {"event_relative_error_vs_own_ode": event_error,
                         "max_common_pre_event_radius_relative_error": radius_error,
                         "sample_count": len(kept), "last_comparison_a": float(kept[-1]["a"]),
                         "pm_event_missing": event is None,
                         "saved_protocol_failures": pm[case].get("failures", [])}
    for j, (event_error, radius_error) in enumerate(zip(event_errors, radius_errors)):
        if event_error is not None:
            ax.scatter(event_error, j-.12, color=BLUE, marker="o", s=28)
        if radius_error is not None:
            ax.scatter(radius_error, j+.12, color=TEAL, marker="s", s=25)
    maximum = max([3.] + [v for v in event_errors + radius_errors if v is not None])
    positive = [v for v in event_errors + radius_errors if v is not None and v > 0]
    minimum = min([1.] + positive)/1.3
    if any(v == 0 for v in event_errors + radius_errors if v is not None):
        ax.set_xscale("symlog", linthresh=minimum)
    else:
        ax.set_xscale("log")
    ax.set_xlim(minimum, maximum*1.25)
    ax.set_xticks([1, 3, 10, 30, 100], ["1", "3", "10", "30", "100"])
    ax.xaxis.set_minor_formatter(NullFormatter())
    for j, event_error in enumerate(event_errors):
        if event_error is None:
            label = "执行失败" if radius_errors[j] is None else "事件缺失"
            ax.text(maximum*1.13, j-.12, label, ha="right", va="center", color="#A44838", fontsize=8)
    ax.axvline(2, color=BLUE, lw=.8, ls=":")
    ax.axvline(3, color=TEAL, lw=.8, ls=(0, (4, 2)))
    matrix_failures = pm_summary.get("failures", [])
    implicated = {case for case in CASES if any(case in str(item.get("name", "")) for item in matrix_failures)}
    if any("_event_vs_main" in str(item.get("name", "")) for item in matrix_failures):
        implicated.add("ref_main")
    labels = [label + (r" $\dagger$" if pm[case].get("failures") or case in implicated
                       or pm[case]["status"] == "execution_failed" else "") for case, label in zip(CASES, CASE_LABELS)]
    ax.set_yticks(range(len(CASES)), labels, fontsize=9.1)
    ax.invert_yaxis()
    ax.set_title(r"$\mathrm{(b)}$ 对各自 $\mathrm{ODE}$ 的误差", fontsize=10.5, pad=9)
    ax.set_xlabel("绝对相对误差（%，对数轴）", fontsize=9.2)
    ax.scatter([], [], color=BLUE, marker="o", s=24, label=r"事件 $a$")
    ax.scatter([], [], color=TEAL, marker="s", s=24, label="半径最大值")
    ax.legend(frameon=False, fontsize=9.0, loc="upper center", bbox_to_anchor=(.5, -.18), ncol=2,
              borderpad=.1, labelspacing=.25, handletextpad=.25, columnspacing=.6, handlelength=.9)
    ax.grid(axis="x", color="#DCE2E6", lw=.55)
    return {"runs": readout, "matrix_protocol_failures": matrix_failures,
            "matrix_comparisons": pm_summary.get("comparisons", {})}


def cooling_panel(fig, ax, rows):
    tables = {}
    fields = ("Tchar_K", "tcool_over_tdyn", "tauchem_over_tcool", "tauchem_over_tdyn")
    for branch in BRANCHES:
        tables[branch] = {(float(row["M200m_Msun_phys"]), float(row["z"])): row
                          for row in rows if row["background"] == branch}
    if tables["reference"].keys() != tables["clock"].keys():
        raise ValueError("Cooling backgrounds do not share the same M,z grid.")
    for key, ref in tables["reference"].items():
        other = tables["clock"][key]
        for field in fields:
            if ref[field] != other[field]:
                raise ValueError(f"Same-M,z cooling changed between backgrounds: {key}, {field}")
    masses = sorted({key[0] for key in tables["reference"]})
    redshifts = sorted({key[1] for key in tables["reference"]})
    if len(tables["reference"]) != len(masses)*len(redshifts):
        raise ValueError("Cooling grid is incomplete; missing cells may not be interpolated.")
    xx = np.log10(masses)
    if not np.allclose(np.diff(xx), .25, rtol=0, atol=1e-12):
        raise ValueError("Unexpected mass-grid spacing.")
    shape = (len(redshifts), len(masses))
    values = np.full(shape, np.nan)
    cold, hot, sensitive, cie = (np.zeros(shape, dtype=bool) for _ in range(4))
    raw_cells = []
    for j, redshift in enumerate(redshifts):
        for i, mass in enumerate(masses):
            row = tables["reference"][(mass, redshift)]
            ratio = float(row["tcool_over_tdyn"])
            temperature = float(row["Tchar_K"])
            cold[j, i] = temperature < 1e4
            hot[j, i] = temperature > 1e8
            sensitive[j, i] = 1e4 <= temperature < 10**4.5
            cie[j, i] = truth(row["CIE_relaxation_slower_than_cooling"]) or truth(row["CIE_relaxation_slower_than_dynamics"])
            if ratio > 0 and np.isfinite(ratio):
                values[j, i] = np.log10(ratio)
            elif not (cold[j, i] or hot[j, i]):
                raise ValueError("Nonfinite cooling ratio in declared atomic domain; report before plotting.")
            raw_cells.append({"log10_M_Msun": float(xx[i]), "z": redshift, "Tchar_K": temperature,
                              "log10_tcool_over_tdyn": None if not np.isfinite(values[j, i]) else float(values[j, i]),
                              "below_atomic_range": bool(cold[j, i]), "above_validated_range": bool(hot[j, i]),
                              "sensitive_temperature": bool(sensitive[j, i]),
                              "local_CIE_warning": bool(cie[j, i])})
    shown = np.ma.masked_where(cold | hot | ~np.isfinite(values), values)
    limit = max(1., float(np.ceil(np.max(np.abs(shown)))))
    cmap = plt.colormaps["RdBu_r"].copy()
    cmap.set_bad("#D9DEE2")
    image = ax.imshow(shown, origin="lower", aspect="auto", interpolation="nearest", cmap=cmap,
                      norm=TwoSlopeNorm(vmin=-limit, vcenter=0, vmax=limit),
                      extent=(xx[0]-.125, xx[-1]+.125, -.5, len(redshifts)-.5))
    for j in range(shape[0]):
        for i in range(shape[1]):
            if cold[j, i] or hot[j, i] or sensitive[j, i]:
                ax.add_patch(Rectangle((xx[i]-.125, j-.5), .25, 1, fill=False, edgecolor="#566571",
                                       lw=0, hatch="///" if cold[j, i] or hot[j, i] else ".."))
            if cie[j, i]:
                ax.plot(xx[i], j, marker="x", ms=3, mew=.55, color="#202A33")
    ax.set_yticks(range(len(redshifts)), [f"{z:g}" for z in redshifts])
    ax.set_xticks([8, 10, 12, 14])
    ax.set_xlabel(r"$\log_{10}(M/M_\odot)$", fontsize=10)
    ax.set_ylabel(r"给定红移 $z$", fontsize=9.5)
    ax.set_title("(c) 条件冷却网格", fontsize=10.5, pad=9)
    cax = fig.add_axes([.467, .232, .014, .223])
    bar = fig.colorbar(image, cax=cax)
    bar.set_label(r"$\log_{10}(t_{\rm cool}/t_{\rm dyn})$", fontsize=9.8, labelpad=3)
    bar.set_ticks([-4, -2, 0, 2, 4])
    bar.ax.tick_params(labelsize=9)
    return {"colour_limits": [-limit, limit], "cold_cells_gray_hatched": int(cold.sum()),
            "hot_cells_gray_hatched": int(hot.sum()),
            "sensitive_cells_dotted": int(sensitive.sum()), "local_CIE_warning_cells": int(cie.sum()),
            "grid_interpolation": "none; redshifts shown as discrete equal-height sampled rows",
            "backgrounds_identical_at_fixed_M_z": True, "cells": raw_cells}


def timescale_panel(ax, events, cooling):
    readouts = {item["label"]: item for item in events["readouts"]}
    if set(readouts) != set(BRANCHES):
        raise ValueError("Require exactly one saved conditional event per background.")
    gyr_s = float(cooling["manifest"]["constants"]["GYR_S"])
    fields = ("tcool_Gyr", "tdyn_Gyr", "tff_Gyr", "tau_H_s", "tau_He_s")
    labels = (r"$t_{\rm cool}$", r"$t_{\rm dyn}$", r"$t_{\rm ff}$", r"$\tau_{\rm H}$", r"$\tau_{\rm He}$")
    manifest = {}
    for branch, offset in (("reference", -.08), ("clock", .08)):
        item = readouts[branch]
        micro = item["microphysics"]
        values = np.asarray([float(micro[field])/(gyr_s if field.endswith("_s") else 1) for field in fields])
        if not np.isfinite(values).all() or np.any(values <= 0):
            raise ValueError("Event contains a nonfinite/nonpositive timescale; report before plotting.")
        ax.semilogy(np.arange(5)+offset, values, color=BRANCH_COLOURS[branch], marker="o" if branch == "reference" else "s",
                    ms=4.5, linestyle="none", fillstyle="none" if branch == "clock" else "full",
                    label="参考事件" if branch == "reference" else "时钟事件")
        manifest[branch] = {"time_Gyr": values.tolist(), "fields": list(fields),
                            "Tchar_K": micro["Tchar_K"], "M_Msun_phys": micro["M200m_Msun_phys"],
                            "tcool_over_tdyn": micro["tcool_over_tdyn"], "source_event": item["source_event"]}
    reference = manifest["reference"]
    ax.set_title(r"$\mathrm{(d)}$ 阈值事件的条件时间尺度", fontsize=10.0, pad=9)
    ax.set_xticks(range(5), labels, fontsize=11)
    ax.set_xlim(-.5, 4.5)
    ax.set_ylabel(r"时间尺度（$\mathrm{Gyr}$）", fontsize=9.5)
    ax.set_yticks([1e-5, 1e-3, .1, 1, 10], [r"$10^{-5}$", r"$10^{-3}$", r"$0.1$", r"$1$", r"$10$"])
    ax.yaxis.set_minor_formatter(NullFormatter())
    ax.grid(axis="y", color="#DCE2E6", lw=.55)
    ax.legend(frameon=False, fontsize=9., loc="upper right", borderpad=.1, labelspacing=.25, handletextpad=.4)
    ax.text(.04, .05, rf"$M={reference['M_Msun_phys']/1e12:.2f}\times10^{{12}}M_\odot$" + "\n"
            + rf"$T_{{\rm ref}}={reference['Tchar_K']/1e6:.2f}\times10^6\,\mathrm{{K}}$", transform=ax.transAxes,
            fontsize=8.8, va="bottom", color=GRAY)
    return manifest


def draw(sphere, ode, pm, pm_summary, grid, events, cooling):
    fig = plt.figure(figsize=(7., 7.5))
    fig.suptitle("受控高密度事件与原初气体的条件冷却", y=.976, fontsize=13)
    fig.text(.5, .940, "标准引力与固定原子速率；六例保存结果，球基准尚未通过。",
             ha="center", fontsize=9.2, color=GRAY)
    ax = fig.add_axes([.102, .603, .35, .271])
    bx = fig.add_axes([.648, .603, .325, .271])
    cx = fig.add_axes([.102, .239, .303, .249])
    dx = fig.add_axes([.648, .239, .325, .249])
    readout = {"radius": radius_panel(ax, sphere, ode, pm), "pm_errors": error_panel(bx, sphere, pm, pm_summary),
               "cooling_grid": cooling_panel(fig, cx, grid), "event_timescales": timescale_panel(dx, events, cooling)}
    fig.axes[-1].set_position([.423, .249, .016, .228])
    for axis in (ax, bx, cx, dx):
        axis.tick_params(labelsize=9)
    particle_shift = 100*pm_summary["comparisons"]["ref_particles128_event_relative_to_main"]
    phase_shift = 100*pm_summary["comparisons"]["ref_shift_event_relative_to_main"]
    fig.text(.102, .156, r"斜线：$T<10^4$ 或 $T>10^8\,\mathrm{K}$；点纹：$10^4\leq T<10^{4.5}\,\mathrm{K}$。", fontsize=8.9)
    fig.text(.102, .124, "叉号：局部平衡松弛较慢；无UV/金属/分子，光学薄需另检验。", fontsize=8.9, color=GRAY)
    fig.text(.102, .092, f"加粒子/相位事件相对主例：{particle_shift:+.2f}% / {phase_shift:+.2f}%，均未过1%门限。", fontsize=8.9)
    fig.text(.102, .060, r"$\mathrm{PM}$ 为标签半径；$\dagger$ 门限未过，虚线 $2\%/3\%$；微小时钟差未分辨。", fontsize=8.9, color=GRAY)
    fig.text(.102, .028, "冷却为规定气体态；未模拟气体热化、恒星形成或反馈。", fontsize=8.9, color=GRAY)
    for suffix in ("pdf", "png"):
        metadata = {"Creator": Path(__file__).name, "CreationDate": None, "ModDate": None} if suffix == "pdf" else None
        fig.savefig(OUT / f"chapter04_halo_cooling_cn.{suffix}", dpi=240, metadata=metadata)
    plt.close(fig)
    return readout


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check-inputs", action="store_true", help="Report readiness; do not draw or write a manifest.")
    args = parser.parse_args()
    if args.check_inputs:
        missing = [str(p.relative_to(ROOT)) for p in required_paths() if not p.is_file()]
        print(json.dumps({"ready": not missing, "missing": missing}, ensure_ascii=False, indent=2))
        return
    sphere, ode, pm, pm_summary, grid, events, cooling, hashes = load_inputs()
    readout = draw(sphere, ode, pm, pm_summary, grid, events, cooling)
    for relative, record in hashes.items():
        if sha(ROOT / relative) != record["sha256"]:
            raise RuntimeError(f"Input changed during drawing: {relative}")
    outputs = [OUT / f"chapter04_halo_cooling_cn.{suffix}" for suffix in ("pdf", "png")]
    manifest = {"scope": "Figure from saved controlled-collapse and prescribed-state cooling outputs; no numerical evolution during plotting.",
                "protocol_sha256": PROTOCOL_SHA, "generator_sha256": sha(__file__), "font": FONT,
                "inputs": hashes, "readout": readout,
                "interpretation": "No observed objects, no galaxy formation, no hydrodynamics. PM failure and missing-event records retained.",
                "outputs": {str(p.relative_to(ROOT)): {"sha256": sha(p), "bytes": p.stat().st_size} for p in outputs}}
    folder = EXPERIMENT / "figures"
    folder.mkdir(parents=True, exist_ok=True)
    (folder / "plot_inputs.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2, allow_nan=False)+"\n")
    print(f"Saved Chapter 4 controlled-collapse/cooling figure; {len(hashes)} input hashes recorded.")


if __name__ == "__main__":
    main()
