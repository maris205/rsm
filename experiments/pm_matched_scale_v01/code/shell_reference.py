#!/usr/bin/env python3
"""Frozen continuous shell reference and exactly two initial PM force calls.

No production PM trajectory, old y-equation, initial-amplitude recalibration,
clock branch, smoothing, or force renormalization is performed here.
"""
from __future__ import annotations

import csv
import gc
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import platform
import resource
import sys
import time
import traceback
from datetime import datetime, timezone

# Importing archived modules must not create/change their bytecode caches.
sys.dont_write_bytecode = True
import numpy as np
import scipy
from scipy.integrate import quad, solve_ivp
from scipy.interpolate import PchipInterpolator
from scipy.sparse import bmat, diags


ROOT = Path(__file__).resolve().parents[1]
REPO = ROOT.parents[1]
OUT = ROOT/"results"
PROTOCOL_SHA256 = "4409835c2fd93047d5df7bcb9675b0771e72fc0bfdfb0dd69e8a89c54cae9bf0"
DEPENDENCIES = {
    "experiments/halo_cooling_v01/code/pm_refined.py": "f8e5bcbc32d726a0345282e228458894ac3b5e75baf0051f65cdbfb05a268602",
    "experiments/halo_cooling_v01/code/spherical_collapse.py": "4e0afba9b783575ec02ad60ec2f6acfaf87da8a0824f6a978b115e89ba522e47",
    "experiments/halo_cooling_v01/results/sphere_summary.json": "eff75c2d2a7144effad867560d54fafc0a3bec3e344f493031e236639b5bc275",
    "experiments/halo_cooling_v01/results/sphere_reference.csv": "f96b0572abb16132aa7425143fdc498445bf2ddf45f9d971767d8c26075ee8a1",
    "experiments/cosmic_bridge_v01/code/pm.py": "1b2607ac0b7ba35c5df9ceb7375910ba75b7bed65a1edec79c0583f067a5d982",
    "experiments/cosmic_bridge_v01/results/background_eps0.csv": "1efdc1e2ef28dccd06cf9d0d0aa83396ab4db2a6a7dc9bc51b48b7ef6a1091a0",
    "experiments/cosmic_bridge_v01/results/background_summary.json": "2bd9d4d15cb60432543bdb485816767e335bda4b47783cfbaa06b4c34a808307",
}
BOX_MPC_H = 10.0
QL, QC = 2.5/BOX_MPC_H, 4.0/BOX_MPC_H
AI, A_MAX = .02, .55
DELTA = .07295182597218097
OM, OR = .315, 9.2e-5
Q_MPC_H = np.concatenate(([.0625], np.arange(1, 33)*.125))
Q_BOX = Q_MPC_H/BOX_MPC_H
RTOL, ATOL, MAX_STEP_N = 2e-10, 2e-12, .01
MASS_TOLERANCE, SHELL_TOLERANCE = 1e-9, 1e-6
FFT_WORKERS = 4


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def array_sha(value):
    return hashlib.sha256(np.ascontiguousarray(value).view(np.uint8)).hexdigest()


def plain(value):
    if isinstance(value, dict):
        return {str(k): plain(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [plain(v) for v in value]
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    return value


def save_json(path, value):
    path.write_text(json.dumps(plain(value), indent=2, ensure_ascii=False, allow_nan=False)+"\n")


class ReferenceBackground:
    """Read-only PCHIP background; smooth acceleration assembled independently."""
    def __init__(self):
        source = json.loads((REPO/"experiments/cosmic_bridge_v01/results/background_summary.json").read_text())
        reference = source["models"]["eps0"]
        self.omega_lambda = float(reference["omega_lambda"])
        self.f_reference = float(reference["initial_pm"]["f_regular"])
        table = np.genfromtxt(REPO/"experiments/cosmic_bridge_v01/results/background_eps0.csv",
                              delimiter=",", names=True)
        self.n = np.log(table["a"])
        self.log_e = PchipInterpolator(self.n, np.log(table["E"]), extrapolate=False)
        self.age = PchipInterpolator(self.n, table["t_Gyr"], extrapolate=False)

    def E(self, n):
        values = np.asarray(n)
        if np.any(values < self.n[0]-2e-12) or np.any(values > self.n[-1]+2e-12):
            raise ValueError("Background extrapolation forbidden")
        return np.exp(self.log_e(np.clip(values, self.n[0], self.n[-1])))

    def evaluate_a(self, a):
        return {"a": a, "E": float(self.E(np.log(a)))}

    def smooth_acceleration(self, n):
        return self.omega_lambda-OR*np.exp(-4*n)


def window(q):
    """Independent factored C2 enclosed-contrast profile; all radii in box units."""
    s = np.clip((np.asarray(q)-QL)/(QC-QL), 0, 1)
    return (1-s)**3*(1+3*s+6*s*s)


def window_q(q):
    s = np.clip((np.asarray(q)-QL)/(QC-QL), 0, 1)
    return -30*s*s*(1-s)**2/(QC-QL)


def initial_map(q, f_reference):
    q = np.asarray(q)
    d = DELTA*window(q)
    y = (1+d)**(-1/3)
    yq = -DELTA*window_q(q)/(3*(1+d)**(4/3))
    yn = -f_reference*d*y/(3*(1+d))
    jacobian = y*y*(y+q*yq)
    return y, yn, jacobian


def mass_integral(q):
    """Independent quadrature of local density excess in Lagrange coordinates."""
    points = [0.]+[x for x in (QL, QC) if x < q]+[float(q)]
    value, error = 0., 0.
    for left, right in zip(points[:-1], points[1:]):
        if right <= left:
            continue
        segment, uncertainty = quad(lambda s: (1-float(initial_map(s, 0)[2]))*s*s,
                                    left, right, epsabs=1e-15, epsrel=1e-11)
        value += segment
        error += uncertainty
    return value, error


def continuous_initial_force(q_vectors):
    """Fast exact mass expression, separately audited against local-density quadrature."""
    q = np.linalg.norm(q_vectors, axis=1)
    y = initial_map(q, 0)[0]
    x = q*y
    integral = (q**3-x**3)/3
    magnitude = np.zeros_like(q)
    np.divide(-integral, x*x, out=magnitude, where=x > 0)
    unit = np.zeros_like(q_vectors)
    np.divide(q_vectors, q[:, None], out=unit, where=q[:, None] > 0)
    return magnitude[:, None]*unit, magnitude, unit


def profile_audit(bg, check, condition):
    radii = np.unique(np.concatenate((np.linspace(0, QC, 513), Q_BOX, [QL, .45, .5])))
    mass_scale = QL**3*DELTA/(3*(1+DELTA))
    rows = []
    for q in radii:
        y, yn, jac = initial_map(q, bg.f_reference)
        x = q*y
        integral, uncertainty = mass_integral(q)
        exact = (q**3-x**3)/3
        rows.append({"q_box": q, "q_Mpc_h": q*BOX_MPC_H, "x_initial_box": x,
                     "y_initial": y, "y_N_initial": yn, "Jacobian": jac,
                     "local_density_ratio": 1/jac,
                     "I_quadrature_box3": integral, "I_closed_form_box3": exact,
                     "I_quadrature_estimated_error_box3": uncertainty,
                     "normalized_integral_difference": (integral-exact)/mass_scale,
                     "g_radial_quadrature_box": -integral/(x*x) if x > 0 else 0.,
                     "g_radial_closed_form_box": -exact/(x*x) if x > 0 else 0.})
    path = OUT/"shell_reference_initial_profile.csv"
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    maximum = max(abs(row["normalized_integral_difference"]) for row in rows)
    check("independent_initial_mass_integrals", maximum, MASS_TOLERANCE)
    condition("initial_Jacobian_positive", all(row["Jacobian"] > 0 for row in rows))
    outer_integral = mass_integral(QC)[0]
    check("outer_boundary_mass_compensation", outer_integral/mass_scale, MASS_TOLERANCE)
    check("outer_boundary_displacement", float(initial_map(QC, 0)[0])-1, MASS_TOLERANCE)
    return {"sample_count": len(rows), "maximum_normalized_mass_error": maximum,
            "normalization_maximum_core_excess_box3": mass_scale,
            "outer_integral_box3": outer_integral,
            "minimum_initial_Jacobian": min(row["Jacobian"] for row in rows),
            "analytic_Jacobian_lower_bound": 1/(1+DELTA),
            "coordinates": "All q,x and g use box units; I has box-volume units",
            "force_evaluation": "Closed expression is used for particles only after independent mass quadrature check; no force fitting"}


def force_region(force, reference, unit, mask):
    actual, ideal, directions = force[mask], reference[mask], unit[mask]
    denominator = float(np.linalg.norm(ideal))
    radial = np.sum(actual*directions, axis=1)
    radial_reference = np.sum(ideal*directions, axis=1)
    transverse = actual-radial[:, None]*directions
    return {"count": int(mask.sum()), "reference_L2": denominator,
            "reference_RMS": denominator/np.sqrt(mask.sum()),
            "vector_relative_L2": float(np.linalg.norm(actual-ideal)/denominator),
            "radial_relative_L2": float(np.linalg.norm(radial-radial_reference)/denominator),
            "transverse_relative_L2": float(np.linalg.norm(transverse)/denominator),
            "mean_radial_actual": float(np.mean(radial)),
            "mean_radial_reference": float(np.mean(radial_reference)),
            "particles_with_inward_force_fraction": float(np.mean(radial < 0)),
            "classification": "descriptive discretization diagnostics; no added force-accuracy pass threshold"}


def initial_force_case(nparticle, bg, pm, check, condition):
    label = f"ref_matched{nparticle}"
    started = time.perf_counter()
    initial = pm.compensated_initial(nparticle, DELTA, bg, force_mesh=nparticle)
    q_vectors, q = initial["q"], initial["q_radius"]
    y, yn, jac = initial_map(q, bg.f_reference)
    expected_x = np.remainder(.5+q_vectors*y[:, None], 1)
    expected_p = AI**2*float(bg.E(np.log(AI)))*q_vectors*yn[:, None]
    displacement_error = np.remainder(initial["positions"]-expected_x+.5, 1)-.5
    check(label+"_independent_initial_position", float(np.max(np.abs(displacement_error))), 1e-12)
    check(label+"_independent_initial_momentum", float(np.max(np.abs(initial["momenta"]-expected_p))), 1e-12)
    condition(label+"_initial_Jacobian_positive", bool(np.all(jac > 0)))
    core = (q >= .25*QL)&(q <= .75*QL)
    condition(label+"_fixed_core_labels", bool(np.array_equal(core, initial["core_mask"])))
    actual_radius = np.linalg.norm(np.remainder(initial["positions"]-.5+.5, 1)-.5, axis=1)
    check(label+"_initial_radius_over_q", float(np.max(np.abs(actual_radius/q-y))), 1e-12)
    mesh = pm.PMGrid(nparticle, FFT_WORKERS)
    force_started = time.perf_counter()
    force, operator_diagnostics = mesh.force(initial["positions"])
    force_seconds = time.perf_counter()-force_started
    # The sole initial force call for this case ends above. Save all values,
    # including exterior forces and any physically inaccurate components.
    force_path = OUT/f"shell_reference_initial_force{nparticle}.npy"
    np.save(force_path, force, allow_pickle=False)
    reference, radial_reference, unit = continuous_initial_force(q_vectors)
    transition, outer = (q > QL)&(q < QC), q >= QC
    regions = {"fixed_core": force_region(force, reference, unit, core),
               "transition": force_region(force, reference, unit, transition)}
    core_norm = regions["fixed_core"]["reference_L2"]
    core_rms = regions["fixed_core"]["reference_RMS"]
    outer_norm = float(np.linalg.norm(force[outer]))
    outer_rms = outer_norm/np.sqrt(outer.sum())
    exterior = {"count": int(outer.sum()), "reference_force": 0.,
                "actual_absolute_L2_box": outer_norm, "actual_RMS_box": outer_rms,
                "actual_max_norm_box": float(np.linalg.norm(force[outer], axis=1).max()),
                "L2_over_core_reference_L2": outer_norm/core_norm,
                "RMS_over_core_reference_RMS": outer_rms/core_rms,
                "classification": "descriptive; no division by individual zero-force references"}
    initial_q4 = float(np.mean(np.sum(unit[core]**4, axis=1))-.6)
    check(label+"_deposited_mass", operator_diagnostics["mass_relative_error"], 1e-12)
    check(label+"_net_force", operator_diagnostics["total_force_normalized"], 1e-11)
    condition(label+"_finite_force", bool(np.all(np.isfinite(force))))
    # Record every radial bin, including empty or zero-reference bins, with
    # absolute RMS values; relative scientific diagnostics use whole regions.
    edges = np.linspace(0, np.sqrt(3)/2, 65)
    bins = np.minimum(np.searchsorted(edges, q, side="right")-1, len(edges)-2)
    radial_actual = np.sum(force*unit, axis=1)
    transverse = force-radial_actual[:, None]*unit
    sums = {
        "reference_squared": np.sum(reference*reference, axis=1),
        "vector_error_squared": np.sum((force-reference)**2, axis=1),
        "radial_error_squared": (radial_actual-radial_reference)**2,
        "transverse_squared": np.sum(transverse*transverse, axis=1),
        "radial_actual": radial_actual, "radial_reference": radial_reference,
    }
    counts = np.bincount(bins, minlength=64)
    binned = {key: np.bincount(bins, weights=values, minlength=64) for key, values in sums.items()}
    rows = []
    for i, count in enumerate(counts):
        row = {"q_left_Mpc_h": edges[i]*BOX_MPC_H, "q_right_Mpc_h": edges[i+1]*BOX_MPC_H,
               "count": int(count)}
        for key, values in binned.items():
            value = float(values[i]/count) if count else None
            row[key.replace("squared", "RMS")] = (np.sqrt(max(value, 0)) if "squared" in key and value is not None else value)
        rows.append(row)
    with (OUT/f"shell_reference_initial_force{nparticle}_bins.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    output = {"label": label, "nparticle_per_side": nparticle, "nmesh_per_side": nparticle,
              "force_calls": 1, "force_seconds": force_seconds, "FFT_workers": FFT_WORKERS,
              "preparation_and_diagnostics_seconds": time.perf_counter()-started,
              "timing_scope": "One initial force call, not a production benchmark or full-run estimate",
              "initial_Q4": initial_q4, "initial_Jacobian_min": float(jac.min()),
              "regions": regions, "exterior": exterior, "operator_diagnostics": operator_diagnostics,
              "initial_metadata": initial["metadata"],
              "initial_positions_array_sha256": array_sha(initial["positions"]),
              "initial_momenta_array_sha256": array_sha(initial["momenta"]),
              "force_array_sha256": array_sha(force), "force_file_sha256": sha(force_path),
              "process_peak_RSS_KiB": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    print(json.dumps({"case": label, "force_seconds": force_seconds,
                      "core_vector_relative_L2": regions["fixed_core"]["vector_relative_L2"],
                      "core_transverse_relative_L2": regions["fixed_core"]["transverse_relative_L2"],
                      "transition_vector_relative_L2": regions["transition"]["vector_relative_L2"]}), flush=True)
    return output


def shells(bg, previous_summary, check, condition):
    size = len(Q_BOX)
    y, yn, jac = initial_map(Q_BOX, bg.f_reference)
    e_i = float(bg.E(np.log(AI)))
    r0 = AI*Q_BOX*y
    v0 = AI*e_i*Q_BOX*(y+yn)
    state0 = np.concatenate((r0, v0))
    core = Q_BOX <= QL
    def rhs(n, state):
        r, v = state[:size], state[size:]
        e = bg.E(n)
        return np.concatenate((v/e, (-OM*Q_BOX**3/(2*r*r)+bg.smooth_acceleration(n)*r)/e))

    def derivative(n, state):
        r = state[:size]
        e = bg.E(n)
        upper = diags(np.full(size, 1/e))
        lower = diags((OM*Q_BOX**3/r**3+bg.smooth_acceleration(n))/e)
        return bmat([[None, upper], [lower, None]], format="csc")

    def core200(n, state):
        # In the continuum all core shells are homologous. Define the first
        # sampled core shell threshold explicitly rather than choosing a shell
        # after viewing numerical event times.
        core_y = state[:size][core]/(np.exp(n)*Q_BOX[core])
        return float(np.min(core_y)-200**(-1/3))
    core200.direction, core200.terminal = -1, True

    def crossing(n, state):
        return float(np.min(np.diff(state[:size])))
    crossing.direction, crossing.terminal = -1, True
    started = time.perf_counter()
    solved = solve_ivp(rhs, (np.log(AI), np.log(A_MAX)), state0, method="Radau",
                       rtol=RTOL, atol=ATOL, max_step=MAX_STEP_N, jac=derivative,
                       events=(core200, crossing), dense_output=True)
    condition("physical_shell_solver_success", solved.success, message=solved.message)
    event_n = float(solved.t_events[0][0]) if len(solved.t_events[0]) else None
    crossing_n = float(solved.t_events[1][0]) if len(solved.t_events[1]) else None
    condition("core200_event_exists", event_n is not None)
    condition("no_sampled_shell_crossing_before_endpoint", crossing_n is None)
    end_n = float(solved.t[-1])
    # Exactly 512 log-spaced pre-endpoint samples, followed by the exact endpoint.
    grid = np.concatenate((np.linspace(np.log(AI), end_n, 512, endpoint=False), [end_n]))
    states = solved.sol(grid)
    radius, velocity = states[:size], states[size:]
    scale, expansion = np.exp(grid), bg.E(grid)
    shell_y = radius/(Q_BOX[:, None]*scale[None, :])
    shell_yn = velocity/(Q_BOX[:, None]*scale[None, :]*expansion[None, :])-shell_y
    minimum_gap = float(np.min(np.diff(radius, axis=0)))
    condition("sampled_shell_order_positive", minimum_gap > 0, minimum_gap_physical_box=minimum_gap)
    check("core_homology_relative", float(np.max(np.abs(shell_y[core]/shell_y[np.flatnonzero(core)[-1]]-1))), SHELL_TOLERANCE)
    old = np.genfromtxt(REPO/"experiments/halo_cooling_v01/results/sphere_reference.csv", delimiter=",", names=True, dtype=None, encoding="utf-8")
    common_end = min(end_n, float(old["N"][-1]))
    common_n = np.linspace(np.log(AI), common_end, 513)
    common = solved.sol(common_n)
    common_y = common[:size]/(Q_BOX[:, None]*np.exp(common_n)[None, :])
    common_yn = common[size:]/(Q_BOX[:, None]*np.exp(common_n)[None, :]*bg.E(common_n)[None, :])-common_y
    comparisons = {}
    for field, values in (("y", common_y), ("y_N", common_yn)):
        reference = PchipInterpolator(old["N"], old[field], extrapolate=False)(common_n)
        error = float(np.max(np.abs(values[core]/reference[None, :]-1)))
        comparisons[field] = error
        check("core_vs_saved_sphere_"+field, error, SHELL_TOLERANCE)
    event_a = float(np.exp(event_n)) if event_n is not None else None
    reference_a = previous_summary["models"]["reference"]["events"]["200"]["a"]
    check("core_event_vs_saved_sphere", None if event_a is None else event_a/reference_a-1, SHELL_TOLERANCE)
    check("outer_shell_remains_background", float(np.max(np.abs(shell_y[-1]-1))), SHELL_TOLERANCE)
    np.savez(OUT/"shell_reference_trajectories.npz", q_Mpc_h=Q_MPC_H, q_box=Q_BOX,
             N=grid, a=scale, t_Gyr=bg.age(grid), R_physical_box=radius,
             V_dR_dH0t=velocity, y=shell_y, y_N=shell_yn, Delta_enclosed=shell_y**-3,
             common_comparison_N=common_n, common_comparison_y=common_y,
             common_comparison_y_N=common_yn)
    return {"status": "complete" if solved.success else "failed", "solver_message": solved.message,
            "method": "Radau", "rtol": RTOL, "atol": ATOL, "max_step_N": MAX_STEP_N,
            "q_Mpc_h": Q_MPC_H, "number_of_shells": size,
            "saved_time_count": len(grid), "pre_endpoint_sample_count": 512,
            "initial_R_physical_box": r0, "initial_V_dR_dH0t": v0,
            "initial_Jacobian": jac, "core200_event_a": event_a,
            "core200_event_t_Gyr": float(bg.age(event_n)) if event_n is not None else None,
            "crossing_event_a": float(np.exp(crossing_n)) if crossing_n is not None else None,
            "endpoint_a": float(scale[-1]), "minimum_sampled_shell_gap_physical_box": minimum_gap,
            "comparison_to_saved_core": comparisons, "comparison_end_a": float(np.exp(common_end)),
            "nfev": solved.nfev, "njev": solved.njev, "wall_seconds": time.perf_counter()-started,
            "scope": "Sampled non-crossing spherical shells; does not establish three-dimensional stability or resolved halos"}


def run():
    started = time.perf_counter()
    OUT.mkdir(exist_ok=True)
    existing = list(OUT.glob("shell_reference*"))
    if existing:
        raise RuntimeError("Refusing automatic rerun/overwrite; existing shell_reference outputs must be preserved")
    output = {"created_utc": datetime.now(timezone.utc).isoformat(),
              "command": "OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 python experiments/pm_matched_scale_v01/code/shell_reference.py",
              "argv": sys.argv, "working_directory": str(Path.cwd()),
              "protocol_sha256": sha(ROOT/"protocol.md"), "code_sha256": sha(Path(__file__)),
              "input_sha256": {path: sha(REPO/path) for path in DEPENDENCIES},
              "software": {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__},
              "resources": {"FFT_workers": FFT_WORKERS, "OMP_NUM_THREADS": os.getenv("OMP_NUM_THREADS"),
                            "OPENBLAS_NUM_THREADS": os.getenv("OPENBLAS_NUM_THREADS")},
              "inputs": {"delta_i": DELTA, "a_initial": AI, "L_Mpc_h": BOX_MPC_H,
                         "qL_Mpc_h": QL*BOX_MPC_H, "qc_Mpc_h": QC*BOX_MPC_H,
                         "epsilon": 0., "mass_integral_tolerance": MASS_TOLERANCE,
                         "core_comparison_tolerance": SHELL_TOLERANCE},
              "checks": [], "initial_force_cases": {}, "exceptions": []}
    checks = output["checks"]
    def check(name, value, threshold, **detail):
        passed = value is not None and np.isfinite(value) and abs(value) <= threshold
        checks.append({"name": name, "value": value, "threshold": threshold,
                       "passed": bool(passed), **detail})
    def condition(name, passed, **detail):
        checks.append({"name": name, "passed": bool(passed), **detail})
    try:
        condition("frozen_protocol_identity", output["protocol_sha256"] == PROTOCOL_SHA256)
        for path, expected in DEPENDENCIES.items():
            condition("dependency_identity:"+path, output["input_sha256"][path] == expected)
        if not all(row["passed"] for row in checks):
            raise RuntimeError("Frozen protocol or archived input identity mismatch")
        previous = json.loads((REPO/"experiments/halo_cooling_v01/results/sphere_summary.json").read_text())
        condition("unchanged_initial_amplitude", previous["delta_i"] == DELTA)
        if previous["delta_i"] != DELTA:
            raise RuntimeError("Stored initial amplitude differs from frozen input")
        bg = ReferenceBackground()
        output["profile"] = profile_audit(bg, check, condition)
        spec = importlib.util.spec_from_file_location("archived_refined_pm_for_shell_diagnostic", REPO/"experiments/halo_cooling_v01/code/pm_refined.py")
        pm = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(pm)
        for nparticle in (64, 128):
            try:
                result = initial_force_case(nparticle, bg, pm, check, condition)
                output["initial_force_cases"][result["label"]] = result
                save_json(OUT/f"shell_reference_initial_force{nparticle}.json", result)
            except Exception as error:
                output["exceptions"].append({"stage": f"initial_force_{nparticle}", "exception": repr(error), "traceback": traceback.format_exc()})
                condition(f"initial_force_{nparticle}_completed", False)
            gc.collect()
        try:
            output["shells"] = shells(bg, previous, check, condition)
        except Exception as error:
            output["exceptions"].append({"stage": "shells", "exception": repr(error), "traceback": traceback.format_exc()})
            condition("shells_completed", False)
        condition("exactly_two_completed_initial_force_calls", len(output["initial_force_cases"]) == 2)
        current_hashes = {path: sha(REPO/path) for path in DEPENDENCIES}
        condition("archived_inputs_unchanged", current_hashes == output["input_sha256"])
    except Exception as error:
        output["exceptions"].append({"stage": "setup", "exception": repr(error), "traceback": traceback.format_exc()})
        condition("setup_complete", False)
    output["failures"] = [row for row in checks if not row["passed"]]
    output["summary"] = {"checks_total": len(checks), "checks_passed": sum(row["passed"] for row in checks),
                         "failures": len(output["failures"]), "exceptions": len(output["exceptions"])}
    output["wall_seconds"] = time.perf_counter()-started
    output["limitations"] = [
        "Force comparisons use whole-region L2 norms; exterior zero-reference forces have absolute and core-normalized diagnostics.",
        "Initial PM force errors are descriptive, not newly fitted acceptance thresholds.",
        "Isolated and periodically copied continuous compensated spheres agree only while non-overlapping and exactly compensated; finite particles and non-spherical evolution need not.",
        "No production PM trajectory, cooling, gas, clock response or observational constraint is calculated here."]
    save_json(OUT/"shell_reference_summary.json", output)
    artifacts = {str(path.relative_to(REPO)): sha(path) for path in sorted(OUT.glob("shell_reference*")) if path.is_file()}
    save_json(OUT/"shell_reference_manifest.json", {"created_utc": datetime.now(timezone.utc).isoformat(),
               "protocol_sha256": PROTOCOL_SHA256, "code_sha256": sha(Path(__file__)),
               "input_sha256": output["input_sha256"], "artifact_sha256": artifacts})
    print(json.dumps(output["summary"]), flush=True)
    return int(bool(output["failures"] or output["exceptions"]))


if __name__ == "__main__":
    raise SystemExit(run())
