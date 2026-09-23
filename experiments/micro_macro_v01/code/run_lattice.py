#!/usr/bin/env python3
"""Periodic phi^4 Hamiltonian lattice, dimensionless benchmark, no observations.

Only writes results/lattice_* under this experiment. The canonical momentum is
p_j = a * velocity_j. Saved velocities are at integer KDK steps, not half steps.
"""
from pathlib import Path
import csv
import hashlib
import json
import platform
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results"
L = 2 * np.pi
SPEED = 1.
OMEGA0 = 1.
AMPLITUDE = .5
FINAL_TIME = 20.
OUTPUT_STEP = .02


def write_json(path, content):
    path.write_text(json.dumps(content, ensure_ascii=False, indent=2,
        default=lambda x: x.item() if isinstance(x, np.generic) else x.tolist()) + "\n")


def energy(q, velocity, spacing, g):
    dq = np.roll(q, -1) - q
    return float(spacing * np.sum(
        .5*velocity**2 + .5*SPEED**2*(dq/spacing)**2
        + .5*OMEGA0**2*q**2 + .25*g*q**4))


def run_case(name, n, g, mode, requested_dt):
    spacing = L/n
    stride = int(round(OUTPUT_STEP/requested_dt))
    if stride < 1 or not np.isclose(stride*requested_dt, OUTPUT_STEP, rtol=0, atol=1e-14):
        raise ValueError("The requested step does not align with the output grid")
    dt = OUTPUT_STEP/stride
    count = int(round(FINAL_TIME/OUTPUT_STEP))+1
    times = np.arange(count)*OUTPUT_STEP
    x = np.arange(n)*spacing
    q = AMPLITUDE*np.cos(mode*x)
    velocity = np.zeros(n)
    acceleration = np.empty(n)
    all_q = np.empty((count,n))
    all_velocity = np.empty((count,n))
    all_energy = np.empty(count)
    max_abs_q = float(np.max(np.abs(q)))

    def force():
        acceleration[1:-1] = q[2:]-2*q[1:-1]+q[:-2]
        acceleration[0] = q[1]-2*q[0]+q[-1]
        acceleration[-1] = q[0]-2*q[-1]+q[-2]
        acceleration[:] *= SPEED**2/spacing**2
        acceleration[:] -= OMEGA0**2*q
        if g:
            acceleration[:] -= g*q*q*q

    start = time.perf_counter()
    force()
    for output_index in range(count):
        all_q[output_index] = q
        all_velocity[output_index] = velocity
        all_energy[output_index] = energy(q,velocity,spacing,g)
        if output_index == count-1:
            break
        for _ in range(stride):
            velocity += .5*dt*acceleration
            q += dt*velocity
            force()
            velocity += .5*dt*acceleration
            max_abs_q = max(max_abs_q,float(np.max(np.abs(q))))
    elapsed = time.perf_counter()-start
    selected_modes = sorted(set((1,3,5,mode)))
    cos_basis = np.cos(np.outer(x,selected_modes))
    cosine_q = (2./n)*all_q@cos_basis
    cosine_velocity = (2./n)*all_velocity@cos_basis
    e0 = all_energy[0]
    frequency = np.sqrt(OMEGA0**2 + 4*SPEED**2/spacing**2*np.sin(mode*spacing/2)**2)
    initial_energy_expected = L*AMPLITUDE**2*frequency**2/4
    if g:
        # Used only for initial mode 1, resolved by every requested grid.
        initial_energy_expected += 3*L*g*AMPLITUDE**4/32
    max_curvature_frequency_bound = np.sqrt(
        OMEGA0**2 + 4*SPEED**2/spacing**2 + 3*g*max_abs_q**2)
    checks = []
    def check(label, passed, **evidence):
        checks.append(dict(name=label,passed=bool(passed),**evidence))
    check("initial_field", np.array_equal(all_q[0],AMPLITUDE*np.cos(mode*x)))
    check("initial_velocity",np.all(all_velocity[0]==0.))
    check("initial_energy_analytic",abs(e0/initial_energy_expected-1)<5e-13,
          relative_error=abs(e0/initial_energy_expected-1))
    check("positive_saved_energy",np.all(all_energy>0),minimum=float(np.min(all_energy)))
    check("finite_saved_fields",np.all(np.isfinite(all_q)) and np.all(np.isfinite(all_velocity)))
    check("strict_local_curvature_step_bound",dt*max_curvature_frequency_bound<2,
          dt_times_frequency_bound=dt*max_curvature_frequency_bound,
          interpretation="Conservative instantaneous curvature scale; not a nonlinear stability theorem.")
    check("same_final_time",times[-1]==FINAL_TIME and
          np.isclose((count-1)*stride*dt,FINAL_TIME,rtol=0,atol=1e-12))
    parameters = dict(name=name,n=n,spacing=spacing,L=L,speed=SPEED,omega0=OMEGA0,
        g=g,initial_amplitude=AMPLITUDE,initial_mode=mode,canonical_momentum="p_j=a*velocity_j",
        requested_dt=requested_dt,actual_dt=dt,steps_per_output=stride,
        internal_steps=(count-1)*stride,output_step=OUTPUT_STEP,final_time=FINAL_TIME,
        output_count=count,selected_cosine_modes=selected_modes,
        method="Störmer–Verlet kick-drift-kick; reuse force after last kick",
        saved_velocity_convention="integer-step velocity after both kicks")
    summary = dict(parameters=parameters,wall_seconds=elapsed,
        initial_energy=e0,final_energy=float(all_energy[-1]),
        maximum_saved_relative_energy_deviation=float(np.max(np.abs(all_energy/e0-1))),
        signed_saved_relative_energy_min=float(np.min(all_energy/e0-1)),
        signed_saved_relative_energy_max=float(np.max(all_energy/e0-1)),
        energy_sampling="Output times only; does not assert exact energy conservation.",
        max_abs_q_over_internal_steps=max_abs_q,
        final_cosine_q={str(k):float(cosine_q[-1,j]) for j,k in enumerate(selected_modes)},
        maximum_absolute_cosine_q={str(k):float(np.max(np.abs(cosine_q[:,j])))
                                  for j,k in enumerate(selected_modes)},
        checks=checks)

    if g==0:
        theta=2*np.arcsin(.5*dt*frequency)
        step_numbers=np.arange(count,dtype=np.int64)*stride
        phase=step_numbers*theta
        spatial=AMPLITUDE*np.cos(mode*x)
        analytic_q=np.cos(phase)[:,None]*spatial
        analytic_velocity=(-frequency*np.sqrt(1-(.5*dt*frequency)**2)
                           *np.sin(phase)[:,None]*spatial)
        analytic_energy_relative=-(.5*dt*frequency)**2*np.sin(phase)**2
        q_error=float(np.max(np.abs(all_q-analytic_q)))
        velocity_error=float(np.max(np.abs(all_velocity-analytic_velocity)))
        energy_error=float(np.max(np.abs(all_energy/e0-1-analytic_energy_relative)))
        check("linear_discrete_time_field",q_error<2e-11,max_abs_error=q_error)
        check("linear_discrete_time_velocity",velocity_error<3e-10,max_abs_error=velocity_error)
        check("linear_energy_oscillation_formula",energy_error<2e-12,
              max_abs_relative_energy_error=energy_error)
        summary["linear_analytic"]=dict(discrete_space_frequency=frequency,
            continuum_frequency=float(np.sqrt(OMEGA0**2+SPEED**2*mode**2)),
            verlet_phase_per_step=theta,verlet_frequency=theta/dt,
            max_abs_q_error=q_error,max_abs_velocity_error=velocity_error,
            max_abs_relative_energy_formula_error=energy_error,
            linear_expected_max_possible_relative_energy_deviation=(.5*dt*frequency)**2,
            comparison="Exact discrete-space and discrete-time KDK trajectory; not continuum PDE.")
    summary["all_checks_passed"]=all(c["passed"] for c in checks)
    np.savez_compressed(OUT/f"lattice_{name}.npz",q=all_q,velocity=all_velocity,
        energy=all_energy,t=times,x=x,cosine_modes=np.array(selected_modes),
        cosine_q=cosine_q,cosine_velocity=cosine_velocity,
        parameters_json=json.dumps(parameters,ensure_ascii=False))
    write_json(OUT/f"lattice_{name}.json",dict(summary=summary,t=times,
        cosine_q={str(k):cosine_q[:,j] for j,k in enumerate(selected_modes)},
        cosine_velocity={str(k):cosine_velocity[:,j] for j,k in enumerate(selected_modes)}))
    mode_rows=[]
    for j,t in enumerate(times):
        row=dict(case=name,n=n,g=g,initial_mode=mode,dt=dt,t=t,energy=all_energy[j],
                 relative_energy_deviation=all_energy[j]/e0-1)
        for k in (1,3,5):
            index=selected_modes.index(k)
            row[f"cosine_q_{k}"]=cosine_q[j,index]
            row[f"cosine_velocity_{k}"]=cosine_velocity[j,index]
        index=selected_modes.index(mode)
        row["cosine_q_initial_mode"]=cosine_q[j,index]
        row["cosine_velocity_initial_mode"]=cosine_velocity[j,index]
        mode_rows.append(row)
    print(json.dumps(dict(case=name,wall_seconds=elapsed,
        max_relative_energy_deviation=summary["maximum_saved_relative_energy_deviation"],
        checks_passed=summary["all_checks_passed"]),ensure_ascii=False),flush=True)
    return summary,mode_rows


def main():
    OUT.mkdir(exist_ok=True,parents=True)
    start=time.perf_counter()
    summaries=[]
    rows=[]
    for n in (64,128,256,512):
        a=L/n
        dt=OUTPUT_STEP/int(np.ceil(OUTPUT_STEP/(.5*a*a)))
        summary,mode_rows=run_case(f"spatial_N{n}",n,1.,1,dt)
        summaries.append(summary)
        rows.extend(mode_rows)
    for dt in (.01,.005,.0025,.00125):
        label=str(dt).replace(".","p")
        summary,mode_rows=run_case(f"time_N128_dt{label}",128,1.,1,dt)
        summaries.append(summary)
        rows.extend(mode_rows)
    for mode in (1,8,24):
        summary,mode_rows=run_case(f"linear_N64_k{mode}",64,0.,mode,.0025)
        summaries.append(summary)
        rows.extend(mode_rows)
    with (OUT/"lattice_modes.csv").open("w",newline="") as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    checks=[dict(case=s["parameters"]["name"],**c) for s in summaries for c in s["checks"]]
    combined=dict(scope="Dimensionless periodic phi4 lattice benchmark; no cosmological calibration or observations.",
        code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        software=dict(python=platform.python_version(),numpy=np.__version__),
        runtime_wall_seconds=time.perf_counter()-start,
        Hamiltonian="H=sum_j a[velocity_j^2/2+speed^2((q_j+1-q_j)/a)^2/2+omega0^2 q_j^2/2+g q_j^4/4]",
        canonical_momentum="p_j=a*velocity_j",
        spatial_policy="dt=.02/ceil(.02/(.5*a^2)); time error is asymptotically subleading for smooth fixed modes.",
        cases=summaries,checks_count=len(checks),all_checks_passed=all(c["passed"] for c in checks),
        failed_checks=[c for c in checks if not c["passed"]],
        files=dict(modes_csv="lattice_modes.csv",
             trajectory_npz=[f"lattice_{s['parameters']['name']}.npz" for s in summaries],
             per_case_json=[f"lattice_{s['parameters']['name']}.json" for s in summaries]))
    write_json(OUT/"lattice_summary.json",combined)
    print(json.dumps(dict(cases=len(summaries),checks=len(checks),
        all_checks_passed=combined["all_checks_passed"],
        runtime_seconds=combined["runtime_wall_seconds"]),ensure_ascii=False),flush=True)
    if not combined["all_checks_passed"]:
        raise SystemExit("Lattice checks failed; inspect lattice_summary.json.")


if __name__=="__main__":
    main()
