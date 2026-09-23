#!/usr/bin/env python3
"""Independent adaptive reference for the arithmetic-residual pulse experiment.

The source residual samples are shared inputs. The natural cubic interpolation,
mass derivative, lattice dynamics, energy and work diagnostics are independently
implemented here, without importing the production integrator or its splines.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
import time
from pathlib import Path

import numpy as np
import scipy
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline

ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "data" / "inputs.npz"
OUT = ROOT / "results"
LENGTH = 2.0 * np.pi
INERTIA = 100.0
R_STAR = np.exp(-2.0)
EPSILON = 0.02
TIMES = np.linspace(0.0, 20.0, 201)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def constitutive(R, spline):
    """Return M^2 and dM^2/dR, including a compactly windowed residual."""
    radius = np.asarray(R)
    chi = np.log(radius / R_STAR)
    mass = 1.0 + 2.0 * (chi ** -2 - 0.25)
    derivative = -4.0 / (radius * chi ** 3)
    inside = (radius > 1.0) & (radius < 3.0)
    # Clip interpolation arguments so that extrapolated values never enter the
    # evaluation outside the pulse, even though its window there is zero.
    coordinate = np.clip(radius, 1.0, 3.0)
    phase = np.pi * (coordinate - 1.0) / 2.0
    sine, cosine = np.sin(phase), np.cos(phase)
    window = np.where(inside, sine ** 4, 0.0)
    window_prime = np.where(inside, 2.0 * np.pi * sine ** 3 * cosine, 0.0)
    sample = spline(coordinate)
    sample_prime = spline(coordinate, 1)
    return (
        mass + EPSILON * window * sample,
        derivative + EPSILON * (window_prime * sample + window * sample_prime),
    )


def solve_case(N, spline, *, rtol=1e-12, atol=1e-14, max_step=0.005):
    spacing = LENGTH / N
    x = np.arange(N) * spacing
    initial = np.r_[0.5 * np.cos(x), np.zeros(N), 1.0, 0.1, 0.0]

    def rhs(_, state):
        field = state[:N]
        speed = state[N : 2 * N]
        radius, radius_speed = state[2 * N : 2 * N + 2]
        mass, mass_prime = constitutive(radius, spline)
        laplacian = (np.roll(field, 1) + np.roll(field, -1) - 2.0 * field) / spacing ** 2
        acceleration = laplacian - mass * field - field ** 3
        square_norm = spacing * np.dot(field, field)
        force = -0.5 * square_norm * mass_prime
        # Work is the change in field energy due to the internal coordinate.
        power = -force * radius_speed
        return np.r_[speed, acceleration, radius_speed, force / INERTIA, power]

    started = time.monotonic()
    solution = solve_ivp(
        rhs,
        (TIMES[0], TIMES[-1]),
        initial,
        t_eval=TIMES,
        method="DOP853",
        rtol=rtol,
        atol=atol,
        max_step=max_step,
    )
    if not solution.success:
        raise RuntimeError(solution.message)
    field = solution.y[:N].T
    speed = solution.y[N : 2 * N].T
    radius, radius_speed, work = solution.y[2 * N :]
    mass, _ = constitutive(radius, spline)
    difference = np.roll(field, -1, axis=1) - field
    field_energy = spacing * np.sum(
        0.5 * speed ** 2
        + 0.5 * (difference / spacing) ** 2
        + 0.5 * mass[:, None] * field ** 2
        + 0.25 * field ** 4,
        axis=1,
    )
    clock_energy = 0.5 * INERTIA * radius_speed ** 2
    total_energy = field_energy + clock_energy
    modes = 2.0 * field @ np.cos(x[:, None] * np.array([1, 3])) / N
    trajectory = dict(
        q=field,
        velocity=speed,
        R=radius,
        Rdot=radius_speed,
        mass_squared=mass,
        field_energy=field_energy,
        clock_energy=clock_energy,
        total_energy=total_energy,
        work=work,
        C1=modes[:, 0],
        C3=modes[:, 1],
    )
    energy_scale = total_energy[0]
    diagnostics = {
        "N": N,
        "rtol": rtol,
        "atol": atol,
        "max_step": max_step,
        "nfev": solution.nfev,
        "elapsed_seconds": time.monotonic() - started,
        "primary_T": float((clock_energy[-1] - clock_energy[0]) / energy_scale),
        "initial_total_energy": float(energy_scale),
        "final_R": float(radius[-1]),
        "min_sampled_R": float(radius.min()),
        "min_sampled_Rdot": float(radius_speed.min()),
        "min_sampled_mass_squared": float(mass.min()),
        "max_sampled_total_energy_relative_error": float(np.max(np.abs(total_energy - energy_scale)) / energy_scale),
        "max_sampled_field_work_relative_error": float(np.max(np.abs(field_energy - field_energy[0] - work)) / energy_scale),
        "max_sampled_clock_work_relative_error": float(np.max(np.abs(clock_energy - clock_energy[0] + work)) / energy_scale),
    }
    return trajectory, diagnostics


def derivative_checks(nodes, splines, caseids):
    # Midpoints avoid the cubic knots; the finite-difference stencil is orders
    # of magnitude smaller than their separation. Check outside/endpoints too.
    probe = np.r_[0.8, 1.0, (nodes[:-1] + nodes[1:]) / 2.0, 3.0, 3.2]
    step = 1e-6
    rows = []
    for caseid in caseids:
        spline = splines[int(caseid)]
        _, analytic = constitutive(probe, spline)
        numerical = (
            constitutive(probe + step, spline)[0]
            - constitutive(probe - step, spline)[0]
        ) / (2.0 * step)
        endpoint = np.array([0.8, 1.0, 3.0, 3.2])
        endpoint_mass, endpoint_derivative = constitutive(endpoint, spline)
        chi = np.log(endpoint / R_STAR)
        endpoint_pulse = endpoint_mass - (1.0 + 2.0 * (chi ** -2 - 0.25))
        endpoint_pulse_derivative = endpoint_derivative + 4.0 / (endpoint * chi ** 3)
        error = float(np.max(np.abs(analytic - numerical) / np.maximum(1.0, np.abs(analytic))))
        rows.append({
            "caseid": int(caseid),
            "max_scaled_finite_difference_error": error,
            "max_endpoint_or_outside_pulse": float(np.max(np.abs(endpoint_pulse))),
            "max_endpoint_or_outside_pulse_derivative": float(np.max(np.abs(endpoint_pulse_derivative))),
            "passed": bool(error < 2e-7 and np.all(endpoint_pulse == 0.0) and np.all(endpoint_pulse_derivative == 0.0)),
        })
    return rows


def select_caseids(family, block, surrogate):
    selected = np.flatnonzero((family == "baseline") | (family == "real")
                             | ((block == 0) & (surrogate == 1)))
    if len(selected) != 11:
        raise ValueError(f"Expected baseline, eight real and two first surrogates; got {len(selected)}")
    if sorted(family[selected].tolist()) != sorted(["baseline"] + ["real"] * 8 + ["iaaft", "exact_spectrum"]):
        raise ValueError("Unexpected selected case families")
    return selected


def save_group(N, selected, splines, family, block, surrogate):
    arrays, rows = [], []
    for caseid in selected:
        print(f"reference N={N} case={caseid} family={family[caseid]} block={block[caseid]}", flush=True)
        trajectory, diagnostics = solve_case(N, splines[int(caseid)])
        diagnostics.update(caseid=int(caseid), family=str(family[caseid]), block=int(block[caseid]), surrogate=int(surrogate[caseid]))
        rows.append(diagnostics)
        arrays.append(trajectory)
    stacked = {key: np.stack([row[key] for row in arrays]) for key in arrays[0]}
    np.savez_compressed(
        OUT / f"reference_N{N}_selected.npz",
        caseids=selected, t=TIMES, x=np.arange(N) * LENGTH / N,
        case_family=family[selected], case_block=block[selected], case_surrogate=surrogate[selected],
        **stacked,
    )
    return stacked, rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    with np.load(INPUT, allow_pickle=False) as inputs:
        eta = np.array(inputs["eta"])
        nodes = np.array(inputs["nodes"])
        family = np.asarray(inputs["case_family"]).astype(str)
        block = np.asarray(inputs["case_block"])
        surrogate = np.asarray(inputs["case_surrogate"])
    if eta.shape != (1593, 128) or nodes.shape != (128,):
        raise ValueError(f"Unexpected source shape {eta.shape}, {nodes.shape}")
    if not np.allclose(nodes, np.linspace(1.0, 3.0, 128), rtol=0.0, atol=1e-14):
        raise ValueError("Unexpected clock nodes")
    selected = select_caseids(family, block, surrogate)
    splines = {int(i): CubicSpline(nodes, eta[i], bc_type="natural") for i in selected}
    summary = {
        "method": "Independent scipy solve_ivp DOP853, separate natural cubic reconstruction and Hamiltonian diagnostics",
        "input_sha256": sha256(INPUT),
        "code_sha256": sha256(Path(__file__)),
        "versions": {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__},
        "parameters": {"L": LENGTH, "I": INERTIA, "Rstar": R_STAR, "b": 2.0, "epsilon": EPSILON, "duration": 20.0, "output_step": 0.1},
        "selected_caseids": selected.tolist(),
        "precision_refinement": {
            "reason": "Initial rtol=1e-10, atol=1e-12, max_step=.01 narrowly failed energy/work and tight-primary guards; all initial artifacts retained as reference_initial_*.",
            "initial_summary": "reference_initial_summary.json",
            "initial_summary_sha256": sha256(OUT / "reference_initial_summary.json") if (OUT / "reference_initial_summary.json").exists() else None,
        },
        "derivative_checks": derivative_checks(nodes, splines, selected),
        "groups": {},
    }
    groups = {}
    for N in [64, 128]:
        groups[N], rows = save_group(N, selected, splines, family, block, surrogate)
        summary["groups"][str(N)] = rows
        (OUT / "reference_summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    tight_caseid = int(np.flatnonzero((family == "real") & (block == 0))[0])
    position = int(np.flatnonzero(selected == tight_caseid)[0])
    print(f"reference tight N=64 case={tight_caseid}", flush=True)
    tight, tight_diagnostics = solve_case(64, splines[tight_caseid], rtol=1e-13, atol=1e-15, max_step=0.0025)
    np.savez_compressed(OUT / "reference_N64_real0_tight.npz", caseids=np.array([tight_caseid]), t=TIMES, x=np.arange(64) * LENGTH / 64, **tight)
    normal_diagnostics = summary["groups"]["64"][position]
    summary["tight_check"] = {
        "caseid": tight_caseid,
        "diagnostics": tight_diagnostics,
        "absolute_primary_T_change": abs(tight_diagnostics["primary_T"] - normal_diagnostics["primary_T"]),
        "max_field_RMS_change": float(np.max(np.sqrt(np.mean((tight["q"] - groups[64]["q"][position]) ** 2, axis=1)))),
        "max_clock_change": float(np.max(np.abs(tight["R"] - groups[64]["R"][position]))),
    }

    prior_path = ROOT.parent / "log_clock_coupling_v01" / "results" / "reference_full_N64_I100.npz"
    if prior_path.exists():
        baseline_position = int(np.flatnonzero(family[selected] == "baseline")[0])
        with np.load(prior_path) as prior:
            positions = np.argmin(np.abs(prior["t"][:, None] - TIMES[None, :]), axis=0)
            if not np.allclose(prior["t"][positions], TIMES, rtol=0.0, atol=1e-12):
                raise ValueError("Old baseline times cannot be aligned exactly")
            baseline = groups[64]
            summary["v04_baseline_regression"] = {
                "source": str(prior_path.relative_to(ROOT.parent.parent)),
                "source_sha256": sha256(prior_path),
                "max_field_RMS_difference": float(np.max(np.sqrt(np.mean((baseline["q"][baseline_position] - prior["q"][positions]) ** 2, axis=1)))),
                "max_clock_difference": float(np.max(np.abs(baseline["R"][baseline_position] - prior["R"][positions]))),
                "max_clock_speed_difference": float(np.max(np.abs(baseline["Rdot"][baseline_position] - prior["Rdot"][positions]))),
            }
    else:
        summary["v04_baseline_regression"] = {"available": False}

    summary["checks"] = {
        "all_derivative_checks_passed": all(row["passed"] for row in summary["derivative_checks"]),
        "positive_sampled_branch_all_selected": all(row["min_sampled_R"] > R_STAR and row["min_sampled_Rdot"] > 0.0 and row["min_sampled_mass_squared"] > 0.0 for rows in summary["groups"].values() for row in rows),
        "initial_mass_is_one_all_selected": all(float(constitutive(1.0, spline)[0]) == 1.0 for spline in splines.values()),
        "energy_and_work_relative_errors_below_1e_8": all(max(row["max_sampled_total_energy_relative_error"], row["max_sampled_field_work_relative_error"], row["max_sampled_clock_work_relative_error"]) < 1e-8 for rows in summary["groups"].values() for row in rows),
        "tight_primary_change_below_1e_9": summary["tight_check"]["absolute_primary_T_change"] < 1e-9,
    }
    summary["output_sha256"] = {p.name: sha256(p) for p in sorted(OUT.glob("reference_*.npz"))}
    (OUT / "reference_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({"checks": summary["checks"], "tight_check": summary["tight_check"], "regression": summary["v04_baseline_regression"]}, indent=2), flush=True)


if __name__ == "__main__":
    main()
