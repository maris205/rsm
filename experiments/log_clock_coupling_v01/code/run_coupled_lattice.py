#!/usr/bin/env python3
"""Autonomous field + positive internal clock Hamiltonian, dimensionless.

The same simultaneous kick-drift-kick map advances canonical p=a*velocity
and P=I*Rdot, stored as velocities because both kinetic masses are constant.
Only writes this experiment's results/lattice_* artifacts. No observations.
"""
from pathlib import Path
import csv
import hashlib
import json
import platform
import time

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"results"
L=2*np.pi
N=128
SPEED=1.
M0=1.
G=1.
INITIAL_AMPLITUDE=.5
RI=1.
RSTAR=float(np.exp(-2.))
CHI_I=2.
TD=10.
INITIAL_RDOT=1/TD
FINAL_TIME=20.
OUTPUT_STEP=.02


def write_json(path,obj):
    path.write_text(json.dumps(obj,ensure_ascii=False,indent=2,
        default=lambda x:x.item() if isinstance(x,np.generic) else x.tolist())+"\n")


def run_case(name,inertia,dt,b=2.,amplitude=INITIAL_AMPLITUDE,roles=()):
    spacing=L/N
    stride=int(round(OUTPUT_STEP/dt))
    if stride<1 or not np.isclose(stride*dt,OUTPUT_STEP,rtol=0,atol=1e-14):
        raise ValueError("Step must divide the output interval exactly.")
    dt=OUTPUT_STEP/stride
    count=int(round(FINAL_TIME/OUTPUT_STEP))+1
    t=np.arange(count)*OUTPUT_STEP
    x=np.arange(N)*spacing
    q=amplitude*np.cos(x)
    velocity=np.zeros(N)
    R=RI
    Rdot=INITIAL_RDOT
    acceleration=np.empty(N)
    history={key:np.empty(count) for key in
        ("R","Rdot","chi","M2","Q2","Hfield","Hclock","Htotal",
         "power_field","integrated_power_field","field_work_residual")}
    saved_q=np.empty((count,N))
    saved_velocity=np.empty((count,N))
    min_R=R
    min_chi=CHI_I
    min_mass2=M0*M0
    max_mass2=M0*M0
    max_abs_q=abs(amplitude)
    work=0.

    def forces(current_R):
        nonlocal min_R,min_chi,min_mass2,max_mass2
        if not np.isfinite(current_R) or current_R<=RSTAR:
            raise FloatingPointError("Clock left the R > Rstar, chi > 0 domain.")
        chi=float(np.log(current_R/RSTAR))
        mass2=M0*M0+b*(chi**-2-CHI_I**-2)
        if not np.isfinite(mass2) or mass2<=0:
            raise FloatingPointError("Nonpositive or nonfinite field mass squared.")
        q2=float(spacing*np.dot(q,q))
        acceleration[1:-1]=q[2:]-2*q[1:-1]+q[:-2]
        acceleration[0]=q[1]-2*q[0]+q[-1]
        acceleration[-1]=q[0]-2*q[-1]+q[-2]
        acceleration[:]*=SPEED**2/spacing**2
        acceleration[:]-=mass2*q
        acceleration[:]-=G*q*q*q
        clock_acceleration=b*q2/(inertia*current_R*chi**3)
        min_R=min(min_R,current_R)
        min_chi=min(min_chi,chi)
        min_mass2=min(min_mass2,mass2)
        max_mass2=max(max_mass2,mass2)
        return clock_acceleration,chi,mass2,q2

    began=time.perf_counter()
    Racc,chi,mass2,q2=forces(R)
    initial_field_energy=None
    for j in range(count):
        dq=np.roll(q,-1)-q
        hfield=float(spacing*np.sum(
            .5*velocity**2+.5*SPEED**2*(dq/spacing)**2+
            .5*mass2*q**2+.25*G*q**4))
        hclock=.5*inertia*Rdot*Rdot
        power=-b*Rdot*q2/(R*chi**3)
        if initial_field_energy is None:
            initial_field_energy=hfield
        saved_q[j]=q
        saved_velocity[j]=velocity
        for key,value in dict(R=R,Rdot=Rdot,chi=chi,M2=mass2,Q2=q2,
            Hfield=hfield,Hclock=hclock,Htotal=hfield+hclock,
            power_field=power,integrated_power_field=work,
            field_work_residual=hfield-initial_field_energy-work).items():
            history[key][j]=value
        if j==count-1:
            break
        for _ in range(stride):
            previous_power=-b*Rdot*q2/(R*chi**3)
            # Simultaneous half kicks at the same old configuration.
            velocity+=.5*dt*acceleration
            Rdot+=.5*dt*Racc
            # Simultaneous drifts using the half-step momenta.
            q+=dt*velocity
            R+=dt*Rdot
            # Both new forces are evaluated at the same new configuration.
            Racc,chi,mass2,q2=forces(R)
            velocity+=.5*dt*acceleration
            Rdot+=.5*dt*Racc
            new_power=-b*Rdot*q2/(R*chi**3)
            work+=.5*dt*(previous_power+new_power)
            max_abs_q=max(max_abs_q,float(np.max(np.abs(q))))
    elapsed=time.perf_counter()-began
    cosine_numbers=np.array((1,3,5))
    basis=np.cos(np.outer(x,cosine_numbers))
    cosine_q=2/N*saved_q@basis
    cosine_velocity=2/N*saved_velocity@basis
    initial_total=history["Htotal"][0]
    fractional_energy=history["Htotal"]/initial_total-1
    initial_omega_squared=M0*M0+4*SPEED**2/spacing**2*np.sin(spacing/2)**2
    expected_initial_field=L*amplitude**2*initial_omega_squared/4+3*L*G*amplitude**4/32
    expected_initial_total=expected_initial_field+.5*inertia*INITIAL_RDOT**2
    checks=[]
    def check(label,passed,**evidence):
        checks.append(dict(name=label,passed=bool(passed),**evidence))
    check("initial_field_and_velocity",
          np.array_equal(saved_q[0],amplitude*np.cos(x)) and np.all(saved_velocity[0]==0.))
    check("initial_clock",history["R"][0]==RI and history["Rdot"][0]==INITIAL_RDOT)
    check("initial_total_energy",abs(initial_total/expected_initial_total-1)<5e-13,
          relative_error=abs(initial_total/expected_initial_total-1))
    check("positive_clock_and_mass_domain",min_R>RSTAR and min_chi>0 and min_mass2>0,
          min_R=min_R,min_chi=min_chi,min_M2=min_mass2)
    check("finite_fields_and_clock",np.all(np.isfinite(saved_q)) and
          np.all(np.isfinite(saved_velocity)) and
          all(np.all(np.isfinite(h)) for h in history.values()))
    check("positive_total_energy",np.all(history["Htotal"]>0) and
          np.all(history["Hfield"]>=0) and np.all(history["Hclock"]>0))
    check("positive_clock_direction",np.all(history["Rdot"]>0))
    check("nonnegative_clock_acceleration_branch",np.min(np.diff(history["Rdot"]))>=-1e-14,
          minimum_output_velocity_increment=float(np.min(np.diff(history["Rdot"]))))
    check("nonpositive_field_exchange_power",np.max(history["power_field"])<=1e-14,
          maximum_power=float(np.max(history["power_field"])))
    check("time_alignment",t[-1]==FINAL_TIME and
          np.isclose((count-1)*stride*dt,FINAL_TIME,rtol=0,atol=1e-12))
    if b==0 or amplitude==0:
        r_error=float(np.max(np.abs(history["R"]-(RI+INITIAL_RDOT*t))))
        check("uncoupled_clock_linear",r_error<1e-11,max_abs_R_error=r_error)
        check("uncoupled_clock_velocity_constant",np.all(history["Rdot"]==INITIAL_RDOT))
        check("uncoupled_exchange_zero",np.all(history["power_field"]==0.) and
              np.all(history["integrated_power_field"]==0.))
    if amplitude==0:
        check("zero_field_exact",np.all(saved_q==0.) and np.all(saved_velocity==0.)
              and np.all(history["Hfield"]==0.))
    if b==0:
        old=ROOT.parent/"micro_macro_v01/results/lattice_time_N128_dt0p00125.npz"
        if old.exists():
            with np.load(old,allow_pickle=False) as reference:
                error=float(np.max(np.abs(saved_q-reference["q"])))
                velocity_error=float(np.max(np.abs(saved_velocity-reference["velocity"])))
            check("b0_matches_previous_autonomous_phi4",max(error,velocity_error)<5e-13,
                  max_abs_q_error=error,max_abs_velocity_error=velocity_error,
                  source=str(old.relative_to(ROOT.parent)),source_sha256=hashlib.sha256(old.read_bytes()).hexdigest())
    parameters=dict(name=name,roles=list(roles),N=N,L=L,spacing=spacing,
        speed=SPEED,m0=M0,g=G,initial_field_amplitude=amplitude,
        initial_field_velocity=0.,b=b,I=inertia,R_initial=RI,
        Rstar=RSTAR,chi_initial=CHI_I,Td=TD,Rdot_initial=INITIAL_RDOT,
        dt=dt,output_dt=OUTPUT_STEP,final_elapsed_time=FINAL_TIME,
        output_count=count,steps_per_output=stride,internal_steps=(count-1)*stride,
        canonical_variables="p_j=a*velocity_j; P=I*Rdot",
        method="Simultaneous separable-Hamiltonian kick-drift-kick Störmer-Verlet",
        field_power_definition="-b*Rdot*Q2/(R*chi**3); Q2=a*sum(q**2)",
        units="Dimensionless benchmark; elapsed time is not calibrated cosmic time.")
    summary=dict(parameters=parameters,wall_seconds=elapsed,
        initial_total_energy=initial_total,final_total_energy=history["Htotal"][-1],
        max_sampled_relative_total_energy_deviation=float(np.max(np.abs(fractional_energy))),
        energy_sampling_note="Energy deviation measured at saved output times; exact conservation is not claimed.",
        initial_field_energy=history["Hfield"][0],final_field_energy=history["Hfield"][-1],
        initial_clock_energy=history["Hclock"][0],final_clock_energy=history["Hclock"][-1],
        clock_energy_gain=history["Hclock"][-1]-history["Hclock"][0],
        field_energy_change=history["Hfield"][-1]-history["Hfield"][0],
        final_R=history["R"][-1],final_Rdot=history["Rdot"][-1],
        final_chi=history["chi"][-1],final_M2=history["M2"][-1],
        external_linear_clock_R_final=RI+INITIAL_RDOT*FINAL_TIME,
        fractional_R_departure_from_external_final=history["R"][-1]/(RI+INITIAL_RDOT*FINAL_TIME)-1,
        min_R_over_steps=min_R,min_chi_over_steps=min_chi,
        min_M2_over_steps=min_mass2,max_M2_over_steps=max_mass2,
        max_abs_q_over_steps=max_abs_q,
        field_linear_curvature_step_scale=dt*np.sqrt(
            max_mass2+4*SPEED**2/spacing**2+3*G*max_abs_q**2),
        curvature_note="Field block diagnostic only; not a proof of nonlinear coupled stability.",
        integrated_field_power_final=history["integrated_power_field"][-1],
        max_abs_field_work_residual=float(np.max(np.abs(history["field_work_residual"]))),
        work_note="Trapezoidal integration at internal steps; independent continuous-time reference is separate.",
        final_cosine_q={str(k):float(cosine_q[-1,j]) for j,k in enumerate(cosine_numbers)},
        max_absolute_cosine_q={str(k):float(np.max(np.abs(cosine_q[:,j])))
                               for j,k in enumerate(cosine_numbers)},
        checks=checks,all_checks_passed=all(c["passed"] for c in checks))
    np.savez_compressed(OUT/f"lattice_{name}.npz",q=saved_q,velocity=saved_velocity,
        t=t,x=x,**history,cosine_modes=cosine_numbers,cosine_q=cosine_q,
        cosine_velocity=cosine_velocity,parameters_json=json.dumps(parameters,ensure_ascii=False))
    write_json(OUT/f"lattice_{name}.json",summary)
    rows=[]
    for j,time_value in enumerate(t):
        row=dict(case=name,I=inertia,b=b,dt=dt,t=time_value,
            R=history["R"][j],Rdot=history["Rdot"][j],M2=history["M2"][j],
            Q2=history["Q2"][j],Hfield=history["Hfield"][j],Hclock=history["Hclock"][j],
            Htotal=history["Htotal"][j],power_field=history["power_field"][j])
        for index,k in enumerate(cosine_numbers):
            row[f"cosine_q_{k}"]=cosine_q[j,index]
            row[f"cosine_velocity_{k}"]=cosine_velocity[j,index]
        rows.append(row)
    print(json.dumps(dict(case=name,seconds=elapsed,R_final=history["R"][-1],
        mass2_final=history["M2"][-1],
        max_relative_energy_deviation=summary["max_sampled_relative_total_energy_deviation"],
        checks_passed=summary["all_checks_passed"]),ensure_ascii=False),flush=True)
    return summary,rows


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    began=time.perf_counter()
    summaries=[]
    all_rows=[]
    for inertia in (10.,100.,1000.):
        roles=("main","time_refinement") if inertia==100 else ("main",)
        summary,rows=run_case(f"main_I{int(inertia)}",inertia,.00125,roles=roles)
        summaries.append(summary)
        all_rows.extend(rows)
    time_map={"0.00125":"main_I100"}
    for dt in (.01,.005,.0025):
        name=f"time_I100_dt{str(dt).replace('.','p')}"
        summary,rows=run_case(name,100.,dt,roles=("time_refinement",))
        summaries.append(summary)
        all_rows.extend(rows)
        time_map[str(dt)]=name
    for name,inertia,dt,b,amplitude in (
        ("null_b0_I100",100.,.00125,0.,INITIAL_AMPLITUDE),
        ("driver_q0_I100",100.,.01,2.,0.)):
        summary,rows=run_case(name,inertia,dt,b=b,amplitude=amplitude,roles=("null_control",))
        summaries.append(summary)
        all_rows.extend(rows)
    with (OUT/"lattice_modes.csv").open("w",newline="") as stream:
        writer=csv.DictWriter(stream,fieldnames=list(all_rows[0]))
        writer.writeheader()
        writer.writerows(all_rows)
    checks=[dict(case=s["parameters"]["name"],**c) for s in summaries for c in s["checks"]]
    combined=dict(scope="Dimensionless autonomous nonlinear lattice coupled to a canonical logarithmic internal clock. No observations or cosmic-time calibration.",
        code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        software=dict(python=platform.python_version(),numpy=np.__version__),
        runtime_seconds=time.perf_counter()-began,
        Hamiltonian="H=I*Rdot^2/2+a*sum[velocity^2/2+speed^2*(delta q/a)^2/2+M2(R)*q^2/2+g*q^4/4]",
        time_refinement_cases=time_map,
        main_cases=["main_I10","main_I100","main_I1000"],
        cases=summaries,checks_count=len(checks),
        all_checks_passed=all(c["passed"] for c in checks),
        failed_checks=[c for c in checks if not c["passed"]],
        trajectory_files=[f"lattice_{s['parameters']['name']}.npz" for s in summaries],
        modes_csv="lattice_modes.csv")
    write_json(OUT/"lattice_summary.json",combined)
    print(json.dumps(dict(cases=len(summaries),checks=len(checks),
        all_checks_passed=combined["all_checks_passed"],
        runtime_seconds=combined["runtime_seconds"]),ensure_ascii=False),flush=True)
    if not combined["all_checks_passed"]:
        raise SystemExit("A coupled-lattice check failed; inspect the saved summary.")


if __name__=="__main__":
    main()
