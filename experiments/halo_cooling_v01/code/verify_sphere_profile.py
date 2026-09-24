#!/usr/bin/env python3
"""Independent continuous initial-profile and saved-artifact audit.

This supplementary audit does not import or rerun spherical_collapse.py,
does not run PM, and does not alter the frozen production protocol/results.
Its implementation tolerances are stated here before its first execution;
they are not additional preregistered scientific acceptance thresholds.
Run from the repository root:
    python experiments/halo_cooling_v01/code/verify_sphere_profile.py
"""
from __future__ import annotations

import hashlib
import json
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import scipy
from scipy.integrate import quad


EXPERIMENT = Path(__file__).resolve().parents[1]
REPOSITORY = EXPERIMENT.parents[1]
EXPECTED_PROTOCOL = "f3a0cc5abf21a4374439db756cbf8c4a3103bcabfa2bee6f9bd1921b28e608f2"
LIMITS = {
    "cartesian_jacobian_relative": 1e-7,
    "compensation_mass_relative_to_background_sphere": 1e-10,
    "boundary_absolute": 1e-12,
    "saved_table_relative": 1e-12,
    "jacobian_lower_bound_slack": 1e-14,
}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def relative_path(path):
    return str(Path(path).relative_to(REPOSITORY))


def run():
    checks = []

    def check(name, value, threshold, **details):
        checks.append({"name": name, "value": float(value), "threshold": float(threshold),
                       "passed": bool(np.isfinite(value) and abs(value) <= threshold),
                       **details})

    def condition(name, passed, **details):
        checks.append({"name": name, "passed": bool(passed), **details})

    summary_path = EXPERIMENT/"results/sphere_summary.json"
    summary = json.loads(summary_path.read_text())
    delta = float(summary["delta_i"])
    q_l = float(summary["inputs"]["R_L_Mpc_h"])
    q_c = float(summary["inputs"]["R_comp_Mpc_h"])
    condition("positive_delta_and_ordered_radii", delta > 0 and 0 < q_l < q_c)

    # Independent implementation, using a factored polynomial different from
    # production's expanded expression. W'=-30 s^2 (1-s)^2/(qc-qL).
    def window(q):
        q = np.asarray(q, dtype=float)
        s = np.clip((q-q_l)/(q_c-q_l), 0.0, 1.0)
        return (1-s)**3*(1+3*s+6*s*s)

    def window_q(q):
        q = np.asarray(q, dtype=float)
        s = np.clip((q-q_l)/(q_c-q_l), 0.0, 1.0)
        return -30*s*s*(1-s)**2/(q_c-q_l)

    def window_qq(q):
        q = np.asarray(q, dtype=float)
        s = np.clip((q-q_l)/(q_c-q_l), 0.0, 1.0)
        return -60*s*(1-s)*(1-2*s)/(q_c-q_l)**2

    def dilation(q):
        return (1+delta*window(q))**(-1/3)

    def dilation_q(q):
        return -delta*window_q(q)*(1+delta*window(q))**(-4/3)/3

    def jacobian(q):
        y = dilation(q)
        return y*y*(y+np.asarray(q)*dilation_q(q))

    def cartesian_map(vector):
        return dilation(np.linalg.norm(vector))*vector

    def numerical_jacobian(vector):
        # Fourth-order central Cartesian differences: independent of the
        # spherical eigenvalue formula used in jacobian().
        step = 2e-5*max(1.0, float(np.linalg.norm(vector)))
        matrix = np.empty((3, 3))
        for axis in range(3):
            displacement = np.eye(3)[axis]*step
            matrix[:, axis] = (-cartesian_map(vector+2*displacement)
                               +8*cartesian_map(vector+displacement)
                               -8*cartesian_map(vector-displacement)
                               +cartesian_map(vector-2*displacement))/(12*step)
        return float(np.linalg.det(matrix))

    for name, value in (("W_core", window(q_l)-1), ("W_outer", window(q_c)),
                        ("Wq_core", window_q(q_l)), ("Wq_outer", window_q(q_c)),
                        ("Wqq_core", window_qq(q_l)), ("Wqq_outer", window_qq(q_c))):
        check(name, value, LIMITS["boundary_absolute"])

    dense_q = np.linspace(0, 1.25*q_c, 100001)
    dense_j = jacobian(dense_q)
    dense_density = 1/dense_j
    lower_bound = 1/(1+delta)
    condition("sampled_positive_jacobian", np.all(dense_j > 0), minimum=float(dense_j.min()))
    check("analytic_jacobian_lower_bound", max(0.0, lower_bound-float(dense_j.min())),
          LIMITS["jacobian_lower_bound_slack"], analytic_lower_bound=lower_bound)
    condition("sampled_nonincreasing_window", np.all(window_q(dense_q) <= 0))
    condition("sampled_window_between_zero_and_one", np.all((window(dense_q) >= 0)&(window(dense_q) <= 1)))

    directions = np.array([[1., 0, 0], [0, 1., 0], [0, 0, 1.],
                           [1., 1., 1.], [1., -2., 3.], [-3., 2., 1.]])
    directions /= np.linalg.norm(directions, axis=1)[:, None]
    radii = [0, .25, .625, 1.25, 2., 2.4, q_l-.01, q_l, q_l+.01,
             2.75, 3., 3.3, 3.7, q_c-.01, q_c, q_c+.01, 4.5]
    jacobian_samples = []
    for radius in radii:
        for direction_index, direction in enumerate(directions):
            exact = float(jacobian(radius))
            numerical = numerical_jacobian(radius*direction)
            error = numerical/exact-1
            jacobian_samples.append({"q_Mpc_h": radius, "direction_index": direction_index,
                                     "analytic": exact, "cartesian_fd": numerical,
                                     "relative_difference": error})
    max_j_error = max(abs(row["relative_difference"]) for row in jacobian_samples)
    check("independent_cartesian_jacobian", max_j_error, LIMITS["cartesian_jacobian_relative"],
          sample_count=len(jacobian_samples))

    # In Eulerian comoving volume dV=4 pi J q^2 dq and rho/rhobar=1/J.
    # Excess mass / (4 pi rho_m0) = integral (1-J) q^2 dq: rho_m0 is
    # the constant mean comoving mass density, not the physical rho(a_i).
    def excess_integrand(q):
        return (1-float(jacobian(q)))*q*q

    core_integral, core_error = quad(excess_integrand, 0, q_l, epsabs=1e-12, epsrel=1e-12)
    transition_integral, transition_error = quad(excess_integrand, q_l, q_c,
                                                 epsabs=1e-12, epsrel=1e-12)
    background_volume_factor = q_c**3/3
    core_exact = q_l**3*delta/(3*(1+delta))
    check("core_excess_mass_integral", (core_integral-core_exact)/background_volume_factor,
          LIMITS["compensation_mass_relative_to_background_sphere"])
    check("transition_net_deficit_integral", (transition_integral+core_exact)/background_volume_factor,
          LIMITS["compensation_mass_relative_to_background_sphere"])
    check("total_compensated_excess_mass_integral", (core_integral+transition_integral)/background_volume_factor,
          LIMITS["compensation_mass_relative_to_background_sphere"])
    for radius in (q_c, 1.125*q_c, 1.25*q_c):
        check("outer_displacement_q_"+str(radius), radius*(float(dilation(radius))-1),
              LIMITS["boundary_absolute"])
        check("outer_enclosed_excess_fraction_q_"+str(radius), 1-float(dilation(radius))**3,
              LIMITS["boundary_absolute"])
        check("outer_local_density_q_"+str(radius), 1/float(jacobian(radius))-1,
              LIMITS["boundary_absolute"])

    # Hash records distinguish verification against the production manifest
    # from newly archived current hashes (CSVs had no prior published hashes).
    tracked_paths = [EXPERIMENT/"protocol.md", EXPERIMENT/"code/spherical_collapse.py",
                     summary_path, EXPERIMENT/"results/sphere_reference.csv",
                     EXPERIMENT/"results/sphere_clock.csv"]
    tracked_paths += [REPOSITORY/path for path in summary["input_sha256"]]
    hashes_before = {relative_path(path): sha(path) for path in tracked_paths}
    condition("protocol_matches_frozen_hash", sha(EXPERIMENT/"protocol.md") == EXPECTED_PROTOCOL)
    condition("production_recorded_protocol_matches", summary["protocol_sha256"] == EXPECTED_PROTOCOL)
    condition("production_code_matches_recorded_hash",
              sha(EXPERIMENT/"code/spherical_collapse.py") == summary["code_sha256"])
    for path, expected in summary["input_sha256"].items():
        condition("production_input_hash_"+path, sha(REPOSITORY/path) == expected,
                  expected_sha256=expected, current_sha256=sha(REPOSITORY/path))

    tables = {}
    for label in ("reference", "clock"):
        path = EXPERIMENT/f"results/sphere_{label}.csv"
        table = np.genfromtxt(path, delimiter=",", names=True, dtype=None, encoding="utf-8")
        condition(label+"_trajectory_row_count", len(table) == summary["inputs"]["trajectory_points"],
                  rows=len(table))
        condition(label+"_trajectory_monotonic_a", np.all(np.diff(table["a"]) > 0))
        condition(label+"_trajectory_label", np.all(table["label"] == label))
        condition(label+"_trajectory_finite_numeric_values", all(np.all(np.isfinite(table[field]))
                  for field in table.dtype.names if field != "label"))
        initial, event = summary["initial"][label], summary["models"][label]["events"]["200"]
        for field, expected in (("a", initial["a_i"]), ("y", initial["y_i"]),
                                ("y_N", initial["y_N_i"]), ("E", initial["E_i"])):
            check(label+"_saved_initial_"+field, float(table[field][0])/expected-1,
                  LIMITS["saved_table_relative"])
        for field in ("a", "y", "y_N", "E", "Delta", "t_Gyr"):
            check(label+"_saved_final_"+field, float(table[field][-1])/event[field]-1,
                  LIMITS["saved_table_relative"])
        tables[label] = {"rows": len(table), "first_a": float(table["a"][0]),
                         "last_a": float(table["a"][-1])}
    hashes_after = {relative_path(path): sha(path) for path in tracked_paths}
    condition("production_inputs_and_artifacts_unchanged_during_audit", hashes_before == hashes_after)

    output = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "scope": "Supplementary independent continuous initial-profile and saved-artifact audit; no ODE/PM production rerun",
        "threshold_status": "Implementation-audit tolerances fixed in this script before its first execution; no change to the frozen production protocol",
        "command": "python experiments/halo_cooling_v01/code/verify_sphere_profile.py",
        "argv": sys.argv, "working_directory": str(Path.cwd()),
        "software": {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__},
        "audit_code_sha256": sha(Path(__file__)), "protocol_sha256": EXPECTED_PROTOCOL,
        "limits": LIMITS,
        "inputs": {"delta_i": delta, "qL_Mpc_h": q_l, "q_comp_Mpc_h": q_c},
        "analytic_result": {
            "W_factorization": "(1-s)^3 (1+3s+6s^2)",
            "W_prime": "-30 s^2 (1-s)^2/(q_comp-qL)",
            "Jacobian": "J=y^2(y+q dy/dq)=(1+delta W)^(-1)[1-delta q W_prime/(3(1+delta W))]",
            "positive_lower_bound": lower_bound,
            "local_density_ratio": "rho/rho_background=1/J; delta_i W is enclosed mean, not local contrast",
            "mass_compensation": "Delta M/(4 pi rho_m0) = integral_0^qc (1-J) q^2 dq = [q^3-(q y)^3]_0^qc/3 = 0; rho_m0 is constant mean comoving mass density",
            "continuum_scope": "Proof concerns the prescribed continuous initial map; it does not remove discrete PM anisotropy, force error or later shell crossing",
        },
        "sampled_extrema": {"Jacobian_min": float(dense_j.min()), "Jacobian_max": float(dense_j.max()),
                            "local_density_ratio_min": float(dense_density.min()),
                            "local_density_ratio_max": float(dense_density.max()),
                            "local_density_min_q_Mpc_h": float(dense_q[np.argmin(dense_density)])},
        "mass_integrals_in_units_4pi_rho_m0_comoving": {
            "core_excess": core_integral, "transition_net_deficit": transition_integral,
            "total_excess": core_integral+transition_integral,
            "core_analytic": core_exact,
            "quad_absolute_error_estimates": [core_error, transition_error]},
        "cartesian_jacobian_samples": jacobian_samples,
        "saved_tables": tables,
        "artifact_current_sha256": hashes_after,
        "hash_interpretation": "Protocol, production code and external inputs are checked against pre-existing recorded hashes. The current summary/CSV hashes are newly archived here; before/after equality checks that this audit did not modify them.",
        "checks": checks,
        "summary": {"checks_total": len(checks), "checks_passed": sum(row["passed"] for row in checks),
                    "failures": sum(not row["passed"] for row in checks)},
        "failures": [row for row in checks if not row["passed"]],
    }
    target = EXPERIMENT/"results/sphere_profile_check.json"
    if target.exists():
        previous = target.read_bytes()
        archived = EXPERIMENT/f"results/sphere_profile_check_previous_{hashlib.sha256(previous).hexdigest()[:12]}.json"
        archived.write_bytes(previous)
    target.write_text(json.dumps(output, ensure_ascii=False, indent=2, allow_nan=False)+"\n")
    print(json.dumps(output["summary"]))
    return int(bool(output["failures"]))


if __name__ == "__main__":
    raise SystemExit(run())
