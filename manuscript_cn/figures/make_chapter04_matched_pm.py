#!/usr/bin/env python3
"""Plot saved matched particle/force-scale collapse diagnostics (Figure 4.6).

No PM kernel, integrator, initial-condition generator or fitting code is
imported. Both frozen new cases are required, including explicit execution
failures or resource stops. The preselected old 64/256 halfstep case is only
context; the complete old matrix remains in Figure 4.5. Missing final results
produce a readiness report and no figure. --check-inputs never draws.
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
from matplotlib.ticker import NullFormatter
import numpy as np

from make_overview_figures import BLUE, FONT, GOLD, GRAY, INK, PALE_BLUE, TEAL


OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[1]
EXPERIMENT = ROOT / "experiments/pm_matched_scale_v01"
RESULTS = EXPERIMENT / "results"
ARCHIVE = ROOT / "experiments/halo_cooling_v01"
PROTOCOL_SHA = "4409835c2fd93047d5df7bcb9675b0771e72fc0bfdfb0dd69e8a89c54cae9bf0"
OLD_PROTOCOL_SHA = "f3a0cc5abf21a4374439db756cbf8c4a3103bcabfa2bee6f9bd1921b28e608f2"
NEW_CASES = ("ref_matched64", "ref_matched128")
OLD_CASE = "ref_halfstep"
DISPLAY_CASES = (*NEW_CASES, OLD_CASE)
LABELS = {"ref_matched64": r"$64/64$", "ref_matched128": r"$128/128$", OLD_CASE: r"旧 $64/256$"}
COLOURS = {"ref_matched64": BLUE, "ref_matched128": TEAL, OLD_CASE: GOLD}
STYLES = {"ref_matched64": "-", "ref_matched128": (0, (5, 2)), OLD_CASE: ":"}
TERMINAL_STATES = ("completed", "stopped_resource_limit", "execution_failed")
plt.rcParams.update({
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.edgecolor": GRAY, "axes.labelcolor": INK, "text.color": INK,
    "xtick.color": INK, "ytick.color": INK, "axes.axisbelow": True,
})


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024*1024), b""):
            digest.update(block)
    return digest.hexdigest()


def json_read(path):
    return json.loads(Path(path).read_text())


def csv_read(path):
    with Path(path).open(newline="") as stream:
        return list(csv.DictReader(stream))


def column(rows, key):
    return np.asarray([float(row[key]) for row in rows])


def inventory():
    """Status-sensitive inventory: a saved failed case is terminal, not missing."""
    summary_path = RESULTS / "matched_pm_summary.json"
    summary = json_read(summary_path) if summary_path.is_file() else None
    paths = [EXPERIMENT / "protocol.md", summary_path, ARCHIVE / "protocol.md",
             ARCHIVE / "results/sphere_summary.json", ARCHIVE / "results/sphere_reference.csv",
             ARCHIVE / "results/spherical_pm_summary.json", ARCHIVE / f"results/pm_{OLD_CASE}.json"]
    statuses = {}
    for case in NEW_CASES:
        record = {} if summary is None else summary.get("runs", {}).get(case, {})
        status = record.get("status", "pending")
        statuses[case] = status
        suffix = "_execution_failure" if status == "execution_failed" else ""
        paths.append(RESULTS / f"pm_{case}{suffix}.json")
        if status != "execution_failed":
            paths.append(RESULTS / f"pm_{case}.npz")
    missing = [str(path.relative_to(ROOT)) for path in paths if not path.is_file()]
    ready = not missing and all(status in TERMINAL_STATES for status in statuses.values())
    return {"ready": ready, "case_status": statuses, "missing": missing}, paths, summary


def load_inputs(paths, summary):
    if sha(EXPERIMENT / "protocol.md") != PROTOCOL_SHA or summary["protocol_sha256"] != PROTOCOL_SHA:
        raise ValueError("New experiment protocol identity changed.")
    if sha(ARCHIVE / "protocol.md") != OLD_PROTOCOL_SHA:
        raise ValueError("Archived protocol identity changed.")
    sphere_path = ARCHIVE / "results/sphere_summary.json"
    sphere = json_read(sphere_path)
    old_summary = json_read(ARCHIVE / "results/spherical_pm_summary.json")
    old = json_read(ARCHIVE / f"results/pm_{OLD_CASE}.json")
    if old["status"] != "completed" or not old["metadata"]["halfstep"]:
        raise ValueError("The preselected archived halfstep context is unavailable.")
    if (old["metadata"]["nparticle"], old["metadata"]["nmesh"]) != (64, 256):
        raise ValueError("Archived context has the wrong particle/force dimensions.")
    if old["metadata"]["protocol_sha256"] != OLD_PROTOCOL_SHA or sphere["protocol_sha256"] != OLD_PROTOCOL_SHA:
        raise ValueError("Archived inputs disagree on protocol identity.")
    if sha(ARCHIVE / f"results/pm_{OLD_CASE}.json") != old_summary["runs"][OLD_CASE]["json_sha256"]:
        raise ValueError("Archived context differs from its matrix summary.")
    if summary["archived_source_sha256"]["halo_cooling_v01/results/sphere_summary.json"] != sha(sphere_path):
        raise ValueError("New cases use a different initial spherical calibration.")
    runs, profiles = {OLD_CASE: old}, {}
    for case, n in zip(NEW_CASES, (64, 128)):
        saved = summary["runs"][case]
        suffix = "_execution_failure" if saved["status"] == "execution_failed" else ""
        path = RESULTS / f"pm_{case}{suffix}.json"
        run = json_read(path)
        runs[case] = run
        if saved["status"] == "execution_failed":
            if run != saved or run["protocol_sha256"] != PROTOCOL_SHA:
                raise ValueError(f"Explicit failure record differs from summary: {case}")
            continue
        if run["status"] != saved["status"] or sha(path) != saved["json_sha256"]:
            raise ValueError(f"Terminal history differs from summary: {case}")
        meta = run["metadata"]
        if (meta["nparticle"], meta["nmesh"], meta["halfstep"]) != (n, n, True):
            raise ValueError(f"Saved case differs from the frozen matched-scale matrix: {case}")
        if meta["background"] != "reference" or meta["protocol_sha256"] != PROTOCOL_SHA:
            raise ValueError(f"Saved background/protocol differs: {case}")
        if meta["spherical_calibration_sha256"] != sha(sphere_path):
            raise ValueError(f"Saved initial calibration differs: {case}")
        array_path = RESULTS / f"pm_{case}.npz"
        if sha(array_path) != run["array_sha256"] or run["array_sha256"] != saved["array_sha256"]:
            raise ValueError(f"Saved terminal profile differs from run/summary: {case}")
        with np.load(array_path, allow_pickle=False) as archive:
            profiles[case] = {key: archive[key].copy() for key in
                              ("q_edges_box", "q_shell_count", "q_shell_y_percentiles", "final_a")}
        if not run["history"] or "event200" not in run:
            raise ValueError(f"Missing saved history or explicit event status: {case}")
    hashes = {str(path.relative_to(ROOT)): {"sha256": sha(path), "bytes": path.stat().st_size} for path in paths}
    return sphere, csv_read(ARCHIVE / "results/sphere_reference.csv"), runs, profiles, hashes


def grid(ax):
    ax.tick_params(labelsize=9)
    ax.grid(axis="y", color="#DCE2E6", lw=.55)


def radii_panel(ax, sphere, sphere_rows, runs):
    event = sphere["models"]["reference"]["events"]["200"]
    h, ql = (float(sphere["inputs"][key]) for key in ("h", "R_L_Mpc_h"))
    before = column(sphere_rows, "a") < event["a"]
    times = np.r_[column(sphere_rows, "t_Gyr")[before], event["t_Gyr"]]
    radius = np.r_[column(sphere_rows, "physical_radius_Mpc_h")[before], event["physical_radius_Mpc_h"]]/h
    ax.plot(times, radius, color=GRAY, lw=1.7, label="连续球")
    readout = {}
    for index, case in enumerate(DISPLAY_CASES):
        run = runs[case]
        if run["status"] == "execution_failed":
            ax.text(.04, .93-index*.1, LABELS[case]+" 执行失败", transform=ax.transAxes, fontsize=9, color=COLOURS[case])
            readout[case] = {"status": run["status"], "failure_record": run}
            continue
        history = run["history"]
        age = column(history, "t_Gyr")
        shown = column(history, "a")*ql*column(history, "y_median")/h
        ax.plot(age, shown, color=COLOURS[case], ls=STYLES[case], lw=1.5, label=LABELS[case])
        ax.plot(age[-1], shown[-1], marker="o", color=COLOURS[case], ms=3)
        readout[case] = {"status": run["status"], "points": len(history),
                         "last_a": float(history[-1]["a"]), "last_t_Gyr": float(age[-1]),
                         "last_radius_proxy_Mpc": float(shown[-1])}
    ax.set_title("(a) 同初态的半径与年龄", fontsize=10.0, pad=9)
    ax.set_xlabel(r"宇宙年龄 $t\;(\mathrm{Gyr})$", fontsize=9.5)
    ax.set_ylabel(r"$a\,q_L y\;(\mathrm{Mpc})$", fontsize=10.5)
    ax.set_ylim(bottom=0)
    ax.legend(frameon=False, fontsize=9, loc="lower center", ncol=2, handlelength=1.3,
              columnspacing=.6, handletextpad=.4, borderpad=.1, labelspacing=.3)
    ax.set_xticks([0, 2, 4, 6])
    grid(ax)
    return readout


def error_readout(run):
    if run["status"] == "execution_failed":
        return {"status": run["status"], "failure_record": run}
    event, ode_event = run["event200"], run["ode_event200"]
    history = run["history"]
    endpoint = min(float(ode_event["a"]), float(event["a"] if event is not None else history[-1]["a"]))
    common = [row for row in history if float(row["a"]) <= endpoint*(1+1e-13)
              and row.get("relative_y_error_vs_ode") is not None]
    if not common:
        raise ValueError("Saved run has no common pre-event radius samples.")
    maximum = max(abs(float(row["relative_y_error_vs_ode"])) for row in common)
    if not np.isclose(maximum, run["diagnostics"]["max_radius_relative_error_vs_ode_common_pre_event"], rtol=1e-10, atol=1e-13):
        raise ValueError("Reconstructed radius maximum differs from the saved diagnostic.")
    event_error = None if event is None else float(event["a"])/float(ode_event["a"])-1
    saved = run["diagnostics"]["event_relative_error_vs_ode"]
    if event_error is None:
        if saved is not None:
            raise ValueError("Missing event has a non-null saved error.")
    elif not np.isclose(event_error, saved, rtol=1e-10, atol=1e-13):
        raise ValueError("Reconstructed event error differs from saved diagnostic.")
    return {"status": run["status"], "event_relative_error_vs_ode": event_error,
            "max_radius_relative_error_common_pre_event": maximum, "common_pre_event_samples": len(common),
            "common_interval_endpoint_a": endpoint, "last_sample_a": float(common[-1]["a"]),
            "saved_scientific_failures": run["failures"]}


def errors_panel(ax, runs, summary):
    result = {case: error_readout(runs[case]) for case in DISPLAY_CASES}
    values, missing = [], []
    for j, case in enumerate(DISPLAY_CASES):
        record = result[case]
        if record["status"] == "execution_failed":
            missing.append((j, "执行失败"))
            continue
        e = record["event_relative_error_vs_ode"]
        r = 100*record["max_radius_relative_error_common_pre_event"]
        if e is None:
            missing.append((j-.12, "事件缺失"))
        else:
            ax.scatter(100*abs(e), j-.12, marker="o", color=BLUE, s=27)
            values.append(100*abs(e))
        ax.scatter(r, j+.12, marker="s", color=TEAL, s=25)
        values.append(r)
    limit = max(100., max([4.] + values)*1.3)
    ax.set_xscale("symlog", linthresh=.01, linscale=1)
    ax.set_xlim(0, limit)
    ticks = [value for value in (0, .01, .1, 1, 10, 100, 1000) if value <= limit]
    ax.set_xticks(ticks, [f"{x:g}" for x in ticks])
    ax.xaxis.set_minor_formatter(NullFormatter())
    for j, label in missing:
        ax.text(.99, 1-(j+.5)/3, label, transform=ax.transAxes, fontsize=9,
                ha="right", va="center", color="#A44838")
    ax.axvline(2, color=BLUE, lw=.85, ls=":")
    ax.axvline(3, color=TEAL, lw=.85, ls=(0, (4, 2)))
    matrix_failed = bool(summary.get("failures"))
    labels = [LABELS[case] + (r" $\dagger$" if runs[case].get("failures")
                            or runs[case]["status"] != "completed"
                            or case in NEW_CASES and matrix_failed else "") for case in DISPLAY_CASES]
    ax.set_yticks(range(3), labels)
    ax.set_ylim(2.5, -.5)
    ax.set_title(r"$\mathrm{(b)}$ 相对连续球的误差", fontsize=10.2, pad=9)
    ax.set_xlabel("绝对相对误差（%，近零线性）", fontsize=9.3)
    ax.scatter([], [], marker="o", color=BLUE, s=23, label=r"事件 $a$")
    ax.scatter([], [], marker="s", color=TEAL, s=23, label="半径最大值")
    ax.legend(frameon=False, fontsize=9, ncol=2, loc="upper center", bbox_to_anchor=(.5, -.18),
              columnspacing=.6, handlelength=.9, handletextpad=.2, borderpad=.1)
    ax.tick_params(labelsize=9)
    ax.grid(axis="x", color="#DCE2E6", lw=.55)
    return {"runs": result, "matrix_comparisons": summary["comparisons"], "matrix_failures": summary["failures"],
            "horizontal_axis": {"scale": "symlog", "linear_threshold_percent": .01, "limits": [0, limit]},
            "context_selection": "Preselected old ref_halfstep only; complete old matrix is retained in Figure 4.5."}


def homology_panel(ax, runs):
    result = {}
    ax.axhline(10, color=GRAY, lw=.9, ls=":", label=r"$10\%$ 门限")
    for index, case in enumerate(NEW_CASES):
        run = runs[case]
        if run["status"] == "execution_failed":
            ax.text(.04, .93-index*.1, LABELS[case]+" 执行失败", transform=ax.transAxes, fontsize=9, color=COLOURS[case])
            result[case] = {"status": run["status"], "failure_record": run}
            continue
        history = run["history"]
        scatter = column(history, "nonhomology_fraction")
        p16, p50, p84 = (column(history, field) for field in ("y_p16", "y_median", "y_p84"))
        reconstructed = (p84-p16)/(2*p50)
        if not np.allclose(scatter, reconstructed, rtol=1e-10, atol=1e-13):
            raise ValueError(f"Saved nonhomology definition disagrees with percentiles: {case}")
        maximum = float(scatter.max())
        if not np.isclose(maximum, run["diagnostics"]["max_nonhomology_fraction"], rtol=1e-10, atol=1e-13):
            raise ValueError(f"Saved maximum nonhomology disagrees with full history: {case}")
        ax.plot(column(history, "a"), 100*scatter, color=COLOURS[case], ls=STYLES[case], lw=1.6, label=LABELS[case])
        result[case] = {"status": run["status"], "history_points": len(history), "maximum_fraction": maximum,
                        "maximum_at_a": float(history[int(np.argmax(scatter))]["a"]),
                        "last_a": float(history[-1]["a"]), "threshold_fraction": .1}
    ax.set_title("(c) 全轨迹的同调散布", fontsize=10.2, pad=9)
    ax.set_xlabel(r"尺度因子 $a$", fontsize=9.5)
    ax.set_ylabel(r"$(P_{84}-P_{16})/(2P_{50})\;(\%)$", fontsize=9.5)
    ax.set_xlim(.02, .55)
    ax.set_ylim(bottom=0)
    ax.set_xticks([.1, .3, .5])
    ax.legend(frameon=False, fontsize=9, loc="upper left", bbox_to_anchor=(0, .95),
              borderpad=.1, handlelength=1.4, labelspacing=.3)
    grid(ax)
    return result


def profiles_panel(ax, sphere, runs, profiles):
    result = {}
    ql = float(sphere["inputs"]["R_L_Mpc_h"])
    ax.axvspan(.25, .75, color=PALE_BLUE, alpha=1., zorder=0)
    for index, case in enumerate(NEW_CASES):
        run = runs[case]
        if run["status"] == "execution_failed":
            ax.text(.04, .93-index*.1, LABELS[case]+" 无终态剖面", transform=ax.transAxes, fontsize=9, color=COLOURS[case])
            result[case] = {"status": run["status"], "failure_record": run}
            continue
        data = profiles[case]
        edges, counts, quartiles = (np.asarray(data[key]) for key in ("q_edges_box", "q_shell_count", "q_shell_y_percentiles"))
        if edges.ndim != 1 or counts.shape != (len(edges)-1,) or quartiles.shape != (len(counts), 3):
            raise ValueError(f"Unexpected saved terminal profile dimensions: {case}")
        if np.any(np.diff(edges) <= 0) or np.any(counts < 0):
            raise ValueError(f"Invalid saved profile bins: {case}")
        populated = counts > 0
        if np.any(~np.isfinite(quartiles[populated])) or np.any(np.diff(quartiles[populated], axis=1) < 0):
            raise ValueError(f"Saved percentiles are nonfinite or unordered: {case}")
        box = float(run["metadata"]["box_mpc_h"])
        center = (edges[1:]+edges[:-1])/2 * box / ql
        shown = np.ma.masked_where(np.broadcast_to(~populated[:, None], quartiles.shape), quartiles)
        final_a = float(data["final_a"])
        if not np.isclose(final_a, run["history"][-1]["a"], rtol=0, atol=1e-13):
            raise ValueError(f"Terminal profile and history end at different a: {case}")
        ax.fill_between(center, shown[:, 0], shown[:, 2], where=populated, color=COLOURS[case], alpha=.18, linewidth=0)
        ax.plot(center, shown[:, 1], color=COLOURS[case], ls=STYLES[case], lw=1.5,
                label=LABELS[case] + rf"，$a_f={final_a:.4f}$")
        result[case] = {"status": run["status"], "final_a": final_a,
                        "q_edges_over_qL": (edges*box/ql).tolist(), "q_shell_count": counts.tolist(),
                        "q_shell_y_percentiles": quartiles.tolist(), "empty_bins_masked": int((~populated).sum())}
    ax.set_title("(d) 各自终态的标签壳剖面", fontsize=10.2, pad=9)
    ax.set_xlabel(r"原始标签半径 $q/q_L$", fontsize=9.5)
    ax.set_ylabel(r"径向缩放 $y=r/(a\,q)$", fontsize=9.5)
    ax.set_xlim(0, 1.6)
    ax.set_ylim(bottom=0)
    ax.set_xticks([0, .5, 1, 1.5])
    ax.legend(frameon=False, fontsize=8.8, loc="upper left", handlelength=1.3, handletextpad=.3, borderpad=.1)
    grid(ax)
    return result


def draw(sphere, sphere_rows, runs, profiles, summary):
    fig = plt.figure(figsize=(7., 7.5))
    fig.suptitle("匹配粒子与力网格尺度的受控塌缩检验", y=.976, fontsize=13)
    fig.text(.5, .940, r"两档新配置与预先指定旧参照；标签为每边 $N_p/N_f$，均为参考背景。",
             ha="center", fontsize=9.2, color=GRAY)
    ax = fig.add_axes([.102, .603, .35, .271])
    bx = fig.add_axes([.648, .603, .325, .271])
    cx = fig.add_axes([.102, .245, .35, .243])
    dx = fig.add_axes([.648, .245, .325, .243])
    readout = {"radius": radii_panel(ax, sphere, sphere_rows, runs),
               "errors": errors_panel(bx, runs, summary),
               "nonhomology": homology_panel(cx, runs),
               "terminal_profiles": profiles_panel(dx, sphere, runs, profiles)}
    fig.text(.102, .157, r"误差线为新例门限 $2\%/3\%$；$\dagger$ 表示逐例或双档门限未过。", fontsize=8.9)
    fig.text(.102, .124, r"剖面带为 $P_{16}$ 至 $P_{84}$；浅蓝为统计核心 $0.25\leq q/q_L\leq0.75$。", fontsize=8.9)
    fig.text(.102, .091, "三例均将三项步长上限减半，实际步长可不同；旧矩阵见图4.5。", fontsize=8.9, color=GRAY)
    fig.text(.102, .058, "两档比较不能给出一般收敛阶，也不检验时钟微小响应或真实星系。", fontsize=8.9, color=GRAY)
    fine = readout["errors"]["runs"]["ref_matched128"].get("max_radius_relative_error_common_pre_event")
    fine_text = "128档无半径误差记录" if fine is None else (
        f"128档半径：{100*fine:.6f}%（{'未过' if fine > .03 else '通过'}3%）")
    pair = summary["comparisons"]["matched128_event_relative_to_matched64"]
    pair_text = "双档事件比较缺失" if pair is None else (
        f"双档事件差：{100*pair:+.3f}%（{'未过' if abs(pair) > .01 else '通过'}1%）")
    fig.text(.102, .025, fine_text + "；" + pair_text + "。", fontsize=8.9, color=INK)
    for extension in ("pdf", "png"):
        metadata = {"Creator": Path(__file__).name, "CreationDate": None, "ModDate": None} if extension == "pdf" else None
        fig.savefig(OUT / f"chapter04_matched_pm_cn.{extension}", dpi=240, metadata=metadata)
    plt.close(fig)
    return readout


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check-inputs", action="store_true", help="Check availability without drawing or writing.")
    args = parser.parse_args()
    readiness, paths, summary = inventory()
    if args.check_inputs or not readiness["ready"]:
        print(json.dumps(readiness, ensure_ascii=False, indent=2))
        return
    sphere, sphere_rows, runs, profiles, hashes = load_inputs(paths, summary)
    readout = draw(sphere, sphere_rows, runs, profiles, summary)
    for relative, item in hashes.items():
        if sha(ROOT / relative) != item["sha256"]:
            raise RuntimeError(f"Input changed while drawing: {relative}")
    outputs = [OUT / f"chapter04_matched_pm_cn.{extension}" for extension in ("pdf", "png")]
    manifest = {"scope": "Saved matched-scale sphere diagnostics only; no solver or fit during plotting.",
                "protocol_sha256": PROTOCOL_SHA, "generator_sha256": sha(__file__), "font": FONT,
                "figure_size_points": [504, 540], "inputs": hashes, "readout": readout,
                "outputs": {str(path.relative_to(ROOT)): {"sha256": sha(path), "bytes": path.stat().st_size} for path in outputs}}
    destination = EXPERIMENT / "figures"
    destination.mkdir(parents=True, exist_ok=True)
    (destination / "plot_inputs.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2, allow_nan=False)+"\n")
    print(f"Saved Figure 4.6 and its manifest; {len(hashes)} input hashes recorded.")


if __name__ == "__main__":
    main()
