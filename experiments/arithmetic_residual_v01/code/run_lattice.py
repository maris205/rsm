#!/usr/bin/env python3
"""Fixed-parameter arithmetic-residual Hamiltonian experiment.

Simultaneous velocity Verlet is canonical kick-drift-kick since p=a*qdot
and P_R=I*Rdot have constant kinetic masses. Input cases are batched but
never coupled to one another. Only results/lattice_* are written.
"""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import platform
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results"
L = 2 * np.pi
RSTAR = float(np.exp(-2))
I = 100.0
EPSILON = 0.02
FINAL_TIME = 20.0
OUTPUT_STEP = 0.1
ENERGY_GATE = 1e-4
WORK_GATE = 1e-4
PRIMARY_REFINEMENT_GATE = 1e-6


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def dump_json(path, obj):
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2,
                              default=lambda x: x.tolist() if isinstance(x, np.ndarray)
                              else x.item() if isinstance(x, np.generic) else str(x)) + "\n")


def mass_and_derivative(R, coefficients, nodes):
    """Natural cubic spline in R, with compact C3 sin^4 window."""
    chi = np.log(R / RSTAR)
    interval = np.searchsorted(nodes, R, side="right") - 1
    interval = np.clip(interval, 0, len(nodes) - 2)
    local_R = R - nodes[interval]
    c = coefficients[np.arange(len(R)), :, interval]
    eta = ((c[:, 0] * local_R + c[:, 1]) * local_R + c[:, 2]) * local_R + c[:, 3]
    eta_prime = (3 * c[:, 0] * local_R + 2 * c[:, 1]) * local_R + c[:, 2]
    angle = np.pi * (R - 1) / 2
    sine = np.sin(angle)
    cosine = np.cos(angle)
    inside = (R > 1) & (R < 3)
    window = np.where(inside, sine**4, 0.0)
    window_prime = np.where(inside, 2 * np.pi * sine**3 * cosine, 0.0)
    mass2 = 1 + 2 * (chi**-2 - 0.25) + EPSILON * window * eta
    derivative = -4 / (R * chi**3) + EPSILON * (window_prime * eta + window * eta_prime)
    return mass2, derivative


def choose_selected(family, block, surrogate):
    selected = [int(np.flatnonzero(family == "baseline")[0])]
    for block_id in sorted(set(block[family == "real"].tolist())):
        for name in ("real", "iaaft", "exact_spectrum"):
            matches = np.flatnonzero((family == name) & (block == block_id))
            if len(matches) == 0:
                raise ValueError(f"Missing block {block_id}, family {name}")
            if name == "real":
                if len(matches) != 1:
                    raise ValueError("Expected one real residual per block")
                selected.append(int(matches[0]))
            else:
                selected.append(int(matches[np.argmin(surrogate[matches])]))
    return np.asarray(selected, dtype=np.int64)


def evolve(coefficients, nodes, case_indices, selected_case_indices, dt, N=64, label="main"):
    began = time.perf_counter()
    count = len(case_indices)
    a = L / N
    stride = round(OUTPUT_STEP / dt)
    if not np.isclose(stride * dt, OUTPUT_STEP, rtol=0, atol=1e-14):
        raise ValueError("dt must exactly divide output interval")
    steps = round(FINAL_TIME / dt)
    output_count = steps // stride + 1
    t = np.arange(output_count) * OUTPUT_STEP
    x = np.arange(N) * a
    q = np.repeat((0.5 * np.cos(x))[None, :], count, axis=0)
    velocity = np.zeros_like(q)
    R = np.ones(count)
    Rdot = np.full(count, 0.1)
    acceleration = np.empty_like(q)
    coefficients = coefficients[case_indices]
    index_lookup = {int(case_id): index for index, case_id in enumerate(case_indices)}
    selected_local = np.asarray([index_lookup[int(i)] for i in selected_case_indices], dtype=np.int64)
    scalar_keys = ("R", "Rdot", "M2", "Q2", "Hfield", "Hclock", "Htotal", "C1", "C3",
                   "power_field", "integrated_power_field", "field_work_residual", "clock_work_residual")
    history = {key: np.empty((count, output_count)) for key in scalar_keys}
    q_selected = np.empty((len(selected_local), output_count, N))
    velocity_selected = np.empty_like(q_selected)
    work = np.zeros(count)
    min_R = R.copy()
    min_Rdot = Rdot.copy()
    min_mass = np.full(count, np.inf)
    max_mass = np.full(count, -np.inf)
    max_abs_derivative = np.zeros(count)
    max_abs_q = np.max(np.abs(q), axis=1)
    cosine1 = 2 * np.cos(x) / N
    cosine3 = 2 * np.cos(3 * x) / N

    def forces():
        if np.any(~np.isfinite(R)) or np.any(R <= RSTAR):
            raise FloatingPointError("Clock left R > Rstar domain")
        mass2, mass_derivative = mass_and_derivative(R, coefficients, nodes)
        if np.any(~np.isfinite(mass2)) or np.any(mass2 <= 0):
            raise FloatingPointError("Nonpositive or nonfinite mass squared")
        q2 = a * np.einsum("ij,ij->i", q, q)
        acceleration[:, 1:-1] = q[:, 2:] - 2 * q[:, 1:-1] + q[:, :-2]
        acceleration[:, 0] = q[:, 1] - 2 * q[:, 0] + q[:, -1]
        acceleration[:, -1] = q[:, 0] - 2 * q[:, -1] + q[:, -2]
        acceleration[:] /= a * a
        acceleration[:] -= mass2[:, None] * q + q * q * q
        Racc = -0.5 * q2 * mass_derivative / I
        np.minimum(min_R, R, out=min_R)
        np.minimum(min_mass, mass2, out=min_mass)
        np.maximum(max_mass, mass2, out=max_mass)
        np.maximum(max_abs_derivative, np.abs(mass_derivative), out=max_abs_derivative)
        np.maximum(max_abs_q, np.max(np.abs(q), axis=1), out=max_abs_q)
        return mass2, mass_derivative, q2, Racc

    mass2, derivative, q2, Racc = forces()

    def save(column):
        dq = np.roll(q, -1, axis=1) - q
        qsq = q * q
        Hfield = a * np.sum(0.5 * velocity * velocity + 0.5 * (dq / a)**2
                            + 0.5 * mass2[:, None] * qsq + 0.25 * qsq * qsq, axis=1)
        Hclock = 0.5 * I * Rdot * Rdot
        if column == 0:
            history["Hfield"][:, 0] = Hfield
            history["Hclock"][:, 0] = Hclock
        values = dict(R=R, Rdot=Rdot, M2=mass2, Q2=q2, Hfield=Hfield, Hclock=Hclock,
                      Htotal=Hfield + Hclock, C1=np.einsum("ij,j->i", q, cosine1),
                      C3=np.einsum("ij,j->i", q, cosine3),
                      power_field=0.5 * q2 * derivative * Rdot,
                      integrated_power_field=work,
                      field_work_residual=Hfield - history["Hfield"][:, 0] - work,
                      clock_work_residual=Hclock - history["Hclock"][:, 0] + work)
        for key, value in values.items():
            history[key][:, column] = value
        q_selected[:, column] = q[selected_local]
        velocity_selected[:, column] = velocity[selected_local]

    save(0)
    for step in range(1, steps + 1):
        old_power = 0.5 * q2 * derivative * Rdot
        velocity += 0.5 * dt * acceleration
        Rdot += 0.5 * dt * Racc
        q += dt * velocity
        R += dt * Rdot
        mass2, derivative, q2, Racc = forces()
        velocity += 0.5 * dt * acceleration
        Rdot += 0.5 * dt * Racc
        np.minimum(min_Rdot, Rdot, out=min_Rdot)
        new_power = 0.5 * q2 * derivative * Rdot
        work += 0.5 * dt * (old_power + new_power)
        if step % stride == 0:
            save(step // stride)
        if step % 2000 == 0:
            print(json.dumps(dict(run=label, step=step, steps=steps,
                                  elapsed_seconds=time.perf_counter()-began)), flush=True)
    total0 = history["Htotal"][:, 0]
    energy_deviation = np.max(np.abs(history["Htotal"] - total0[:, None]), axis=1) / total0
    field_work_relative = np.max(np.abs(history["field_work_residual"]), axis=1) / total0
    clock_work_relative = np.max(np.abs(history["clock_work_residual"]), axis=1) / total0
    primary = (history["Hclock"][:, -1] - history["Hclock"][:, 0]) / total0
    diagnostics = dict(min_R=min_R, min_Rdot=min_Rdot, min_mass_squared=min_mass,
                       max_mass_squared=max_mass, max_abs_mass_derivative=max_abs_derivative,
                       max_abs_q=max_abs_q, sampled_max_relative_energy_deviation=energy_deviation,
                       sampled_max_relative_field_work_residual=field_work_relative,
                       sampled_max_relative_clock_work_residual=clock_work_relative,
                       primary_T=primary)
    checks = [dict(name="finite_trajectories", passed=all(np.all(np.isfinite(v)) for v in history.values())),
              dict(name="R_greater_than_Rstar", passed=np.all(min_R > RSTAR)),
              dict(name="positive_mass_squared", passed=np.all(min_mass > 0)),
              dict(name="identical_initial_total_energy", passed=np.ptp(total0) < 1e-13),
              dict(name="initial_window_zero_mass_one", passed=np.all(history["M2"][:, 0] == 1)),
              dict(name="energy_accuracy", passed=np.max(energy_deviation) < ENERGY_GATE,
                   value=np.max(energy_deviation), threshold=ENERGY_GATE),
              dict(name="field_work_accuracy", passed=np.max(field_work_relative) < WORK_GATE,
                   value=np.max(field_work_relative), threshold=WORK_GATE),
              dict(name="clock_work_accuracy", passed=np.max(clock_work_relative) < WORK_GATE,
                   value=np.max(clock_work_relative), threshold=WORK_GATE)]
    parameters = dict(N=N, L=L, a=a, I=I, b=2.0, epsilon=EPSILON, Rstar=RSTAR,
                      R_initial=1.0, Rdot_initial=0.1, field_amplitude=0.5,
                      dt=dt, final_time=FINAL_TIME, output_step=OUTPUT_STEP)
    payload = dict(t=t, x=x, case_indices=case_indices, selected_case_indices=selected_case_indices,
                   selected_q=q_selected, selected_velocity=velocity_selected, **history, **diagnostics,
                   parameters_json=json.dumps(parameters))
    summary = dict(parameters=parameters, cases=count, selected_cases=len(selected_local),
                   runtime_seconds=time.perf_counter()-began, checks=checks,
                   all_checks_passed=all(bool(c["passed"]) for c in checks),
                   ranges={key: dict(min=float(np.min(value)), max=float(np.max(value)))
                           for key, value in diagnostics.items()})
    return payload, summary


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=ROOT / "data" / "inputs.npz")
    parser.add_argument("--smoke", action="store_true", help="Zero-input full-duration development check")
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    began = time.perf_counter()
    if args.smoke:
        nodes = np.linspace(1, 3, 128)
        payload, summary = evolve(np.zeros((1, 4, 127)), nodes, np.array([0]), np.array([0]), .00125,
                                  label="smoke_zero_residual")
        np.savez_compressed(OUT / "lattice_smoke_zero.npz", **payload)
        dump_json(OUT / "lattice_smoke_zero.json", summary)
        print(json.dumps(summary, default=lambda x: x.item() if isinstance(x, np.generic) else str(x)), flush=True)
        return
    inputs = np.load(args.input, allow_pickle=False)
    coefficients = inputs["spline_coefficients"]
    nodes = inputs["nodes"]
    family = inputs["case_family"]
    block = inputs["case_block"]
    surrogate = inputs["case_surrogate"]
    count = len(family)
    if coefficients.shape != (count, 4, len(nodes)-1):
        raise ValueError("Unexpected spline coefficient layout")
    selected = choose_selected(family, block, surrogate)
    main, summary_main = evolve(coefficients, nodes, np.arange(count), selected, .00125, label="main")
    metadata = dict(case_family=family, case_block=block, case_surrogate=surrogate)
    np.savez_compressed(OUT / "lattice_main.npz", **main, **metadata)
    dump_json(OUT / "lattice_main.json", summary_main)
    refined, summary_refined = evolve(coefficients, nodes, selected, selected, .000625, label="refined_selected")
    np.savez_compressed(OUT / "lattice_refined_selected.npz", **refined,
                        **{key: value[selected] for key, value in metadata.items()})
    dump_json(OUT / "lattice_refined_selected.json", summary_refined)
    baseline = int(np.flatnonzero(family == "baseline")[0])
    baseline_mode = main["C1"][baseline]
    c1_relative_rms = np.sqrt(np.mean((main["C1"]-baseline_mode[None, :])**2, axis=1)
                              / np.mean(baseline_mode**2))
    primary_difference = refined["primary_T"] - main["primary_T"][selected]
    baseline_selected_local = int(np.flatnonzero(selected == baseline)[0])
    main_contrast = main["primary_T"][selected] - main["primary_T"][baseline]
    fine_contrast = refined["primary_T"] - refined["primary_T"][baseline_selected_local]
    q_difference = refined["selected_q"] - main["selected_q"]
    q_relative_rms = np.sqrt(np.mean(q_difference**2, axis=(1, 2))
                             / np.mean(refined["selected_q"]**2, axis=(1, 2)))
    refinement = dict(case_indices=selected, primary_T_fine_minus_main=primary_difference,
                      primary_contrast_main=main_contrast, primary_contrast_fine=fine_contrast,
                      primary_contrast_fine_minus_main=fine_contrast-main_contrast,
                      q_global_relative_rms_difference=q_relative_rms,
                      max_absolute_primary_difference=float(np.max(np.abs(primary_difference))),
                      max_q_global_relative_rms_difference=float(np.max(q_relative_rms)))
    checks = [dict(run="main", **c) for c in summary_main["checks"]]
    checks += [dict(run="refined", **c) for c in summary_refined["checks"]]
    checks.append(dict(name="primary_refinement_accuracy", run="comparison",
                       value=refinement["max_absolute_primary_difference"], threshold=PRIMARY_REFINEMENT_GATE,
                       passed=refinement["max_absolute_primary_difference"] < PRIMARY_REFINEMENT_GATE))
    rows = []
    diagnostic_keys = ("min_R", "min_Rdot", "min_mass_squared", "max_mass_squared", "max_abs_mass_derivative",
                       "sampled_max_relative_energy_deviation", "sampled_max_relative_field_work_residual",
                       "sampled_max_relative_clock_work_residual")
    for index in range(count):
        rows.append(dict(case_index=index, family=str(family[index]), block=int(block[index]),
                         surrogate=int(surrogate[index]), primary_T=float(main["primary_T"][index]),
                         primary_T_minus_baseline=float(main["primary_T"][index]-main["primary_T"][baseline]),
                         C1_relative_rms_difference_from_baseline=float(c1_relative_rms[index]),
                         R_final=float(main["R"][index, -1]), Rdot_final=float(main["Rdot"][index, -1]),
                         M2_final=float(main["M2"][index, -1]),
                         **{key: float(main[key][index]) for key in diagnostic_keys}))
    with (OUT / "lattice_cases.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    refinement_rows = [dict(case_index=int(case_id), family=str(family[case_id]), block=int(block[case_id]),
                            surrogate=int(surrogate[case_id]), primary_T_main=float(main["primary_T"][case_id]),
                            primary_T_fine=float(refined["primary_T"][i]),
                            primary_T_fine_minus_main=float(primary_difference[i]),
                            contrast_main=float(main_contrast[i]), contrast_fine=float(fine_contrast[i]),
                            contrast_fine_minus_main=float(fine_contrast[i]-main_contrast[i]),
                            q_global_relative_rms_difference=float(q_relative_rms[i]))
                       for i, case_id in enumerate(selected)]
    with (OUT / "lattice_refinement.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(refinement_rows[0]))
        writer.writeheader()
        writer.writerows(refinement_rows)
    previous_path = ROOT.parent / "log_clock_coupling_v01" / "results" / "reference_full_N64_I100.npz"
    previous_comparison = None
    if previous_path.exists():
        previous = np.load(previous_path, allow_pickle=False)
        previous_time_indices = np.searchsorted(previous["t"], main["t"]-1e-12)
        if not np.allclose(previous["t"][previous_time_indices], main["t"], rtol=0, atol=1e-10):
            raise ValueError("Previous-stage reference sampling mismatch")
        new_baseline_q = main["selected_q"][int(np.flatnonzero(selected == baseline)[0])]
        previous_q = previous["q"][previous_time_indices]
        previous_T = ((previous["clock_energy"][-1]-previous["clock_energy"][0])
                      / previous["total_energy"][0])
        previous_comparison = dict(path=str(previous_path.relative_to(ROOT.parent)), sha256=sha256(previous_path),
                                   primary_T_reference=float(previous_T),
                                   primary_T_main_minus_reference=float(main["primary_T"][baseline]-previous_T),
                                   q_global_relative_rms=float(np.sqrt(np.mean((new_baseline_q-previous_q)**2)
                                                                      / np.mean(previous_q**2))),
                                   R_max_absolute_difference=float(np.max(np.abs(main["R"][baseline]
                                                                     -previous["R"][previous_time_indices]))),
                                   note="Prior v0.4 independent DOP853 N64, same zero-residual model; no fit or input retuning.")
    combined = dict(scope="Dimensionless arithmetic-residual response sensitivity experiment; no observations or physical inference.",
                    input_path=str(args.input.relative_to(ROOT)), input_sha256=sha256(args.input),
                    code_sha256=sha256(__file__),
                    software=dict(python=platform.python_version(), numpy=np.__version__),
                    runtime_seconds=time.perf_counter()-began,
                    main=summary_main, refined_selected=summary_refined, refinement=refinement,
                    baseline_previous_v04_reference=previous_comparison,
                    primary="T=(Hclock(final)-Hclock(initial))/Htotal(initial)",
                    secondary="RMS(C1_case-C1_baseline)/RMS(C1_baseline), over 201 output times",
                    checks=checks, checks_count=len(checks), all_checks_passed=all(bool(c["passed"]) for c in checks),
                    failed_checks=[c for c in checks if not c["passed"]],
                    data_layout="Scalar histories (case,time), selected_q and selected_velocity (selected_case,time,grid)",
                    selected_case_indices=selected)
    dump_json(OUT / "lattice_summary.json", combined)
    print(json.dumps(dict(cases=count, checks=len(checks), all_checks_passed=combined["all_checks_passed"],
                          max_relative_energy_deviation=summary_main["ranges"]["sampled_max_relative_energy_deviation"]["max"],
                          max_absolute_primary_refinement_difference=refinement["max_absolute_primary_difference"],
                          runtime_seconds=combined["runtime_seconds"])), flush=True)
    if not combined["all_checks_passed"]:
        raise SystemExit("A numerical gate failed; preserve and inspect the saved results")


if __name__ == "__main__":
    main()
