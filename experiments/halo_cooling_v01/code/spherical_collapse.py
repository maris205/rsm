#!/usr/bin/env python3
"""Controlled fixed-mass spherical overdensity in saved smooth backgrounds.

This is a collapse-event reference, not a virialization or halo finder.
Run only after this experiment's protocol.md is frozen:
    python experiments/halo_cooling_v01/code/spherical_collapse.py

API:
    backgrounds = load_backgrounds()
    delta_i, trace = calibrate_delta(backgrounds['reference'])
    initial = initial_conditions(delta_i, backgrounds)
    ref = integrate_sphere(backgrounds['reference'], initial['reference'])
    ref.evaluate_a(a), ref.events, ref.turnaround

Main trajectories solve y=r/(a R_L), d y/d ln(a). The independent Radau
implementation solves physical r/R_L and dr/(H0 R_L dt) with a separately
assembled smooth-fluid acceleration. No production PM or cooling is run here.
"""
from __future__ import annotations

import csv
import hashlib
import json
import platform
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import scipy
from scipy.integrate import solve_ivp
from scipy.interpolate import PchipInterpolator
from scipy.optimize import brentq


ROOT = Path(__file__).resolve().parents[1]
BRIDGE = ROOT.parent/"cosmic_bridge_v01"
AI = 0.02
TARGET_A200 = 0.5
R_L_MPC_H = 2.5
R_COMP_MPC_H = 4.0
OMEGA_M0 = 0.315
OMEGA_R0 = 9.2e-5
H0_KM_S_MPC = 67.4
H = H0_KM_S_MPC/100
H0_GYR = H0_KM_S_MPC/3.0856775814913673e19*(365.25*86400)*1e9
RHO_CRIT0_UNITS = 2.77536627e11  # (Msun/h)/(Mpc/h)^3
THRESHOLDS = (100.0,200.0,500.0)
RTOL = 2e-10
ATOL = 2e-12
MAX_STEP_N = 0.01
CHECK_THRESHOLD = 1e-6
OUTPUT_COUNT = 4097
EXPECTED_PROTOCOL_SHA256 = "f3a0cc5abf21a4374439db756cbf8c4a3103bcabfa2bee6f9bd1921b28e608f2"


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def scalar(value):
    array = np.asarray(value)
    return float(array) if array.ndim == 0 else array


def serializable(value):
    if isinstance(value,dict):
        return {key:serializable(item) for key,item in value.items()}
    if isinstance(value,(list,tuple)):
        return [serializable(item) for item in value]
    if isinstance(value,np.ndarray):
        return value.tolist()
    if isinstance(value,np.generic):
        return value.item()
    return value


class BackgroundTable:
    def __init__(self, label, epsilon, table_path, omega_lambda, f_reference):
        self.label,self.epsilon = str(label),float(epsilon)
        self.omega_m0,self.omega_r0 = OMEGA_M0,OMEGA_R0
        self.omega_lambda = float(omega_lambda)
        self.f_reference = float(f_reference)
        self.path = Path(table_path)
        table = np.genfromtxt(self.path,delimiter=",",names=True)
        self.n = np.log(table["a"])
        for column in ("a","E","dlnE","t_Gyr","Omega_chi"):
            if np.any(~np.isfinite(table[column])):
                raise ValueError(f"Nonfinite required background column: {column}")
        if np.any(np.diff(self.n)<=0) or np.any(table["E"]<=0):
            raise ValueError("Ordered positive expansion table required")
        self._log_e = PchipInterpolator(self.n,np.log(table["E"]),extrapolate=False)
        self._f = PchipInterpolator(self.n,table["dlnE"],extrapolate=False)
        self._age = PchipInterpolator(self.n,table["t_Gyr"],extrapolate=False)
        if self.epsilon > 0:
            # Direct acceleration uses rho_chi+3p_chi instead of F=E_N/E.
            if np.any(~np.isfinite(table["w_chi"])):
                raise ValueError("Finite scalar pressure required for positive epsilon")
            pressure_combination = table["Omega_chi"]*(1+3*table["w_chi"])
            self._pressure = PchipInterpolator(self.n,pressure_combination,extrapolate=False)
        else:
            self._pressure = None

    def evaluate_n(self,n):
        n = np.asarray(n,dtype=float)
        if np.any(n<self.n[0]-2e-12) or np.any(n>self.n[-1]+2e-12):
            raise ValueError("Sphere requested background extrapolation")
        n = np.clip(n,self.n[0],self.n[-1])
        e = np.exp(self._log_e(n))
        output = {"a":np.exp(n),"N":n,"E":e,"dlnE":self._f(n),
                  "Omega_m":self.omega_m0*np.exp(-3*n)/(e*e),
                  "t_Gyr":self._age(n)}
        return {key:scalar(value) for key,value in output.items()}

    def evaluate_a(self,a):
        return self.evaluate_n(np.log(a))

    def smooth_acceleration(self,n):
        """Return (a_ddot/a)/H0^2 after subtracting homogeneous dust."""
        values = self.evaluate_n(n)
        pressure = 0.0 if self._pressure is None else self._pressure(np.asarray(n))
        return scalar(self.omega_lambda-self.omega_r0*np.exp(-4*np.asarray(n))
                      -0.5*np.asarray(values["E"])**2*pressure)


class EdSBackground:
    label = "EdS_exact_benchmark"
    epsilon = 0.0
    omega_m0 = 1.0
    omega_r0 = 0.0
    omega_lambda = 0.0
    f_reference = 1.0

    def evaluate_n(self,n):
        n = np.asarray(n,dtype=float)
        return {"a":scalar(np.exp(n)),"N":scalar(n),"E":scalar(np.exp(-1.5*n)),
                "dlnE":scalar(np.full_like(n,-1.5)),
                "Omega_m":scalar(np.ones_like(n)),
                "t_Gyr":scalar((2/3)*np.exp(1.5*n)/H0_GYR)}

    def evaluate_a(self,a):
        return self.evaluate_n(np.log(a))

    def smooth_acceleration(self,n):
        return scalar(np.zeros_like(np.asarray(n),dtype=float))


def load_backgrounds():
    summary = json.loads((BRIDGE/"results/background_summary.json").read_text())
    reference_f = summary["models"]["eps0"]["initial_pm"]["f_regular"]
    backgrounds = {}
    for label,key,epsilon in (("reference","eps0",0.0),("clock","eps1e-4",1e-4)):
        backgrounds[label] = BackgroundTable(label,epsilon,BRIDGE/f"results/background_{key}.csv",
                                             summary["models"][key]["omega_lambda"],reference_f)
    return backgrounds


def initial_conditions(delta_i,backgrounds):
    if delta_i <= 0:
        raise ValueError("Positive initial enclosed overdensity required")
    ref = backgrounds["reference"]
    e_ref = ref.evaluate_a(AI)["E"]
    f_ref = ref.f_reference
    y_i = (1+delta_i)**(-1/3)
    yn_ref = -f_ref*delta_i*y_i/(3*(1+delta_i))
    output = {}
    for label,bg in backgrounds.items():
        e = bg.evaluate_a(AI)["E"]
        yn_i = yn_ref*e_ref/e
        output[label] = {"a_i":AI,"delta_i":float(delta_i),"y_i":float(y_i),
                         "y_N_i":float(yn_i),"E_i":float(e),"E_reference_i":float(e_ref),
                         "f_reference_i":float(f_ref),
                         "common_peculiar_momentum_per_RL":float(AI**2*e*yn_i),
                         "velocity_convention":"finite-ai delta_N=f_reference*delta; matched physical p, not exact nonlinear pure growing mode"}
    return output


def compensation_window(q_mpc_h):
    """C2 enclosed-overdensity window in Lagrangian radius, not local density."""
    q = np.asarray(q_mpc_h,dtype=float)
    s = np.clip((q-R_L_MPC_H)/(R_COMP_MPC_H-R_L_MPC_H),0.0,1.0)
    return scalar(1-10*s**3+15*s**4-6*s**5)


def initial_shell_map(q_mpc_h,delta_i,f_reference,E_reference,E_model=None):
    """Return y, y_N and common p/q for the specified compensated mapping."""
    if E_model is None:
        E_model = E_reference
    enclosed_delta = delta_i*np.asarray(compensation_window(q_mpc_h))
    y = (1+enclosed_delta)**(-1/3)
    yn = -f_reference*enclosed_delta*y/(3*(1+enclosed_delta))*E_reference/E_model
    return {"W":compensation_window(q_mpc_h),"y":scalar(y),"y_N":scalar(yn),
            "momentum_per_q":scalar(AI**2*E_model*yn)}


class SphereSolution:
    def __init__(self,background,initial,solution,coordinate):
        self.background,self.initial,self.solution,self.coordinate = background,initial,solution,coordinate
        self.nmin,self.nmax = float(solution.t[0]),float(solution.t[-1])
        self.events = {}
        for i,threshold in enumerate(THRESHOLDS):
            if len(solution.t_events[i]):
                event_n = float(solution.t_events[i][0])
                self.events[str(int(threshold))] = self.evaluate_n(event_n)
            else:
                self.events[str(int(threshold))] = None
        self.turnaround = (self.evaluate_n(float(solution.t_events[-1][0]))
                           if len(solution.t_events[-1]) else None)

    def evaluate_n(self,n):
        n = np.asarray(n,dtype=float)
        if np.any(n<self.nmin-2e-12) or np.any(n>self.nmax+2e-12):
            raise ValueError("Sphere evaluation outside integrated interval")
        n = np.clip(n,self.nmin,self.nmax)
        bg = self.background.evaluate_n(n)
        first,second = self.solution.sol(n)
        a,e = np.exp(n),np.asarray(bg["E"])
        if self.coordinate == "comoving_ratio":
            y,yn = first,second
        else:
            y = first/a
            yn = second/(e*a)-y
        result = {**bg,"y":scalar(y),"y_N":scalar(yn),"Delta":scalar(y**-3),
                  "delta":scalar(y**-3-1),"comoving_radius_Mpc_h":scalar(R_L_MPC_H*y),
                  "physical_radius_Mpc_h":scalar(R_L_MPC_H*a*y),
                  "physical_radial_velocity_km_s":scalar(100*R_L_MPC_H*a*e*(y+yn)),
                  "peculiar_radial_velocity_km_s":scalar(100*R_L_MPC_H*a*e*yn)}
        return result

    def evaluate_a(self,a):
        return self.evaluate_n(np.log(a))


def integrate_sphere(background,initial,method="DOP853",rtol=RTOL,atol=ATOL,max_step=MAX_STEP_N,
                     coordinate="comoving_ratio",a_final=1.0):
    n_i = float(np.log(initial["a_i"]))
    if coordinate == "comoving_ratio":
        state0 = (initial["y_i"],initial["y_N_i"])

        def rhs(n,state):
            y,yn = state
            if y <= 0:
                raise FloatingPointError("Reached nonpositive sphere radius before terminal event")
            bg = background.evaluate_n(n)
            return yn,-(2+bg["dlnE"])*yn-0.5*bg["Omega_m"]*(y**-3-1)*y

        def y_of_state(n,state):
            return state[0]

        def turnaround(n,state):
            return state[0]+state[1]
    elif coordinate == "physical_radius":
        a,e = initial["a_i"],initial["E_i"]
        state0 = (a*initial["y_i"],a*e*(initial["y_i"]+initial["y_N_i"]))

        def rhs(n,state):
            radius,velocity = state
            if radius <= 0:
                raise FloatingPointError("Reached nonpositive physical radius before terminal event")
            bg = background.evaluate_n(n)
            acceleration = (-background.omega_m0/(2*radius**2)
                            +background.smooth_acceleration(n)*radius)
            return velocity/bg["E"],acceleration/bg["E"]

        def y_of_state(n,state):
            return state[0]*np.exp(-n)

        def turnaround(n,state):
            return state[1]
    else:
        raise ValueError("Unknown sphere coordinates")
    events = []
    for threshold in THRESHOLDS:
        def event(n,state,threshold=threshold):
            return y_of_state(n,state)-threshold**(-1/3)
        event.direction = -1
        event.terminal = threshold == THRESHOLDS[-1]
        events.append(event)
    turnaround.direction = -1
    turnaround.terminal = False
    events.append(turnaround)
    solution = solve_ivp(rhs,(n_i,np.log(a_final)),state0,method=method,
                         rtol=rtol,atol=atol,max_step=max_step,events=events,dense_output=True)
    if not solution.success:
        raise RuntimeError(solution.message)
    return SphereSolution(background,initial,solution,coordinate)


def calibrate_delta(reference,target_a=TARGET_A200):
    """Calibrate a numerical initial condition, never an observation fit."""
    trace = []

    def residual(delta):
        initial = initial_conditions(delta,{"reference":reference})["reference"]
        solution = integrate_sphere(reference,initial)
        event = solution.events["200"]
        # No event by a=1 is a censored positive bracket, not a false event.
        value = (1.0 if event is None else event["a"])-target_a
        trace.append({"delta_i":float(delta),"a200":None if event is None else event["a"],
                      "no_event_by_a1":event is None,"root_residual":float(value)})
        return value
    delta = brentq(residual,1e-3,0.5,xtol=2e-12,rtol=2e-12)
    return float(delta),trace


def cycloid_delta(theta):
    return 4.5*(theta-np.sin(theta))**2/(1-np.cos(theta))**3


def cycloid_setup():
    theta200 = brentq(lambda theta:cycloid_delta(theta)-200,np.pi,2*np.pi-1e-5,
                     xtol=1e-14,rtol=1e-14)
    tau200 = (2/3)*TARGET_A200**1.5
    b = tau200/(theta200-np.sin(theta200))
    theta_i = brentq(lambda theta:b*(theta-np.sin(theta))-(2/3)*AI**1.5,
                    1e-4,theta200,xtol=1e-14,rtol=1e-14)
    yi = cycloid_delta(theta_i)**(-1/3)
    yn_i = yi*(1.5*np.sin(theta_i)*(theta_i-np.sin(theta_i))/(1-np.cos(theta_i))**2-1)
    background = EdSBackground()
    initial = {"a_i":AI,"delta_i":cycloid_delta(theta_i)-1,"y_i":yi,"y_N_i":yn_i,
               "E_i":background.evaluate_a(AI)["E"],"velocity_convention":"exact EdS cycloid growing initial state"}
    events = {}
    for threshold in THRESHOLDS:
        theta = brentq(lambda theta:cycloid_delta(theta)-threshold,np.pi,2*np.pi-1e-5,
                       xtol=1e-14,rtol=1e-14)
        tau = b*(theta-np.sin(theta))
        a = (1.5*tau)**(2/3)
        events[str(int(threshold))] = {"theta":theta,"a":a,"t_Gyr":tau/H0_GYR,
                                      "Delta":threshold,"delta_linear":0.6*(0.75*(theta-np.sin(theta)))**(2/3)}
    return background,initial,{"B_H0":b,"theta_i":theta_i,"events":events,
        "turnaround":{"a":(1.5*b*np.pi)**(2/3),"t_Gyr":b*np.pi/H0_GYR,"Delta":9*np.pi**2/16},
        "formal_collapse_delta_linear":0.6*(1.5*np.pi)**(2/3),
        "formal_collapse_a":(3*b*np.pi)**(2/3),
        "interpretation":"Exact EdS pre-shell-crossing benchmark. No virialization is inferred from Delta=200."}


def write_trajectory(path,solution):
    event = solution.events["200"]
    if event is None:
        raise RuntimeError("Cannot export requested trajectory without Delta=200 event")
    grid = np.linspace(np.log(AI),event["N"],OUTPUT_COUNT)
    values = solution.evaluate_n(grid)
    fields = ["a","N","E","dlnE","t_Gyr","y","y_N","Delta","delta",
              "comoving_radius_Mpc_h","physical_radius_Mpc_h",
              "physical_radial_velocity_km_s","peculiar_radial_velocity_km_s"]
    with Path(path).open("w",newline="") as stream:
        writer = csv.DictWriter(stream,fieldnames=["label","epsilon"]+fields)
        writer.writeheader()
        for i in range(len(grid)):
            writer.writerow({"label":solution.background.label,"epsilon":solution.background.epsilon,
                             **{key:float(values[key][i]) for key in fields}})


def run():
    started = time.perf_counter()
    protocol = ROOT/"protocol.md"
    if not protocol.exists():
        raise RuntimeError("Freeze protocol.md before executing sphere calculations")
    if sha(protocol) != EXPECTED_PROTOCOL_SHA256:
        raise RuntimeError("Sphere protocol bytes differ from the explicitly frozen SHA256")
    (ROOT/"results").mkdir(exist_ok=True)
    checks = []

    def check(name,value,threshold=CHECK_THRESHOLD,**detail):
        checks.append({"name":name,"value":float(value),"threshold":float(threshold),
                       "passed":bool(np.isfinite(value) and abs(value)<=threshold),**serializable(detail)})

    def failure(name,error):
        checks.append({"name":name,"passed":False,"exception":repr(error)})

    output = {"created_utc":datetime.now(timezone.utc).isoformat(),
        "scope":"Fixed-mass high-density collapse events; not virialization, halo identification or observation fitting",
        "protocol_sha256":sha(protocol),"code_sha256":sha(Path(__file__)),
        "execution":{"command":"python experiments/halo_cooling_v01/code/spherical_collapse.py",
                     "argv":sys.argv,"python_executable":sys.executable,"working_directory":str(Path.cwd())},
        "software":{"python":platform.python_version(),"numpy":np.__version__,"scipy":scipy.__version__},
        "inputs":{"initial_a":AI,"target_reference_a200":TARGET_A200,"R_L_Mpc_h":R_L_MPC_H,
                  "R_comp_Mpc_h":R_COMP_MPC_H,"Omega_m0":OMEGA_M0,"Omega_r0":OMEGA_R0,
                  "H0_km_s_Mpc":H0_KM_S_MPC,"h":H,"thresholds_mean_density":THRESHOLDS,
                  "rtol":RTOL,"atol":ATOL,"max_step_N":MAX_STEP_N,"check_threshold":CHECK_THRESHOLD,
                  "delta_root_bracket":[1e-3,0.5],"trajectory_points":OUTPUT_COUNT,
                  "core_mass_Msun_h":4*np.pi*R_L_MPC_H**3*RHO_CRIT0_UNITS*OMEGA_M0/3,
                  "profile":"W=1 for q<=2.5; 1-10s^3+15s^4-6s^5 for 2.5<q<4; zero outside",
                  "profile_meaning":"delta_i W is enclosed mean overdensity in the initial Lagrangian mapping, not local shell density"},
        "models":{},"benchmarks":{},"event_rows":[]}
    try:
        backgrounds = load_backgrounds()
        output["input_sha256"] = {str(bg.path.relative_to(ROOT.parent.parent)):sha(bg.path)
                                  for bg in backgrounds.values()}
        output["input_sha256"]["experiments/cosmic_bridge_v01/results/background_summary.json"] = sha(BRIDGE/"results/background_summary.json")
        delta_i,trace = calibrate_delta(backgrounds["reference"])
        initials = initial_conditions(delta_i,backgrounds)
        output.update(delta_i=delta_i,initial=initials,calibration_trace=trace)
        check("matched_initial_physical_momentum",
              initials["clock"]["common_peculiar_momentum_per_RL"]
              /initials["reference"]["common_peculiar_momentum_per_RL"]-1)
        solutions = {}
        for label,bg in backgrounds.items():
            try:
                solution = integrate_sphere(bg,initials[label])
                solutions[label] = solution
                output["models"][label] = {"epsilon":bg.epsilon,"initial":initials[label],
                    "events":solution.events,"turnaround":solution.turnaround,
                    "solver_evaluations":solution.solution.nfev}
                all_events = [("turnaround",solution.turnaround)]+[("Delta"+threshold,event)
                                                                   for threshold,event in solution.events.items()]
                mass_msun_h = output["inputs"]["core_mass_Msun_h"]
                for event_name,event in all_events:
                    if event is not None:
                        output["event_rows"].append({
                            "label":label,"eps":bg.epsilon,"epsilon":bg.epsilon,"event":event_name,
                            "a_event":event["a"],"z_event":1/event["a"]-1,"t_Gyr":event["t_Gyr"],
                            "Delta":event["Delta"],"M_Msun_h":mass_msun_h,
                            "M_solar_physical":mass_msun_h/H,"qL_Mpc_h":R_L_MPC_H,
                            "qL_Mpc_physical_comoving":R_L_MPC_H/H,
                            "R_physical_Mpc_h":event["physical_radius_Mpc_h"],
                            "R_physical_Mpc":event["physical_radius_Mpc_h"]/H,
                            "physical_radial_velocity_km_s":event["physical_radial_velocity_km_s"],
                            "peculiar_radial_velocity_km_s":event["peculiar_radial_velocity_km_s"],
                            "interpretation":"fixed-mass spherical event; not an equilibrium/virialized halo"})
                for threshold,event in solution.events.items():
                    if event is None:
                        failure(label+"_event_"+threshold,RuntimeError("Event not reached before a=1"))
                    else:
                        check(label+"_threshold_"+threshold,event["Delta"]/float(threshold)-1)
                if solution.turnaround is None:
                    failure(label+"_turnaround",RuntimeError("No turnaround event"))
                else:
                    ta=solution.turnaround
                    check(label+"_turnaround_radial_derivative",ta["y_N"]+ta["y"])
                if label == "reference":
                    check("reference_target_a200",solution.events["200"]["a"]/TARGET_A200-1)
                grid = np.linspace(np.log(AI),solution.events["200"]["N"],1025)
                base = solution.evaluate_n(grid)
                for name,kwargs in (("tight",dict(rtol=2e-12,atol=2e-14,max_step=0.005)),
                                    ("physical_R_Radau",dict(method="Radau",coordinate="physical_radius",
                                                             rtol=2e-11,atol=2e-13,max_step=0.005))):
                    alternate = integrate_sphere(bg,initials[label],**kwargs)
                    # Slightly different event locations require a common interval.
                    common_end = min(solution.events["200"]["N"],alternate.events["200"]["N"])
                    common = np.linspace(np.log(AI),common_end,1025)
                    left,right = solution.evaluate_n(common),alternate.evaluate_n(common)
                    for field in ("y","y_N"):
                        check(label+"_"+name+"_"+field,float(np.max(np.abs(np.asarray(right[field])/left[field]-1))))
                    for threshold in ("100","200","500"):
                        check(label+"_"+name+"_event_a"+threshold,
                              alternate.events[threshold]["a"]/solution.events[threshold]["a"]-1)
                    check(label+"_"+name+"_turnaround_a",
                          alternate.turnaround["a"]/solution.turnaround["a"]-1)
                    output["models"][label][name] = {"events":alternate.events,"turnaround":alternate.turnaround,
                                                       "solver_evaluations":alternate.solution.nfev}
                write_trajectory(ROOT/f"results/sphere_{label}.csv",solution)
            except Exception as error:
                failure(label+"_execution",error)
        if all(label in solutions for label in ("reference","clock")):
            output["paired_event_shifts"] = {}
            for threshold in ("100","200","500"):
                ref,clock = solutions["reference"].events[threshold],solutions["clock"].events[threshold]
                if ref is not None and clock is not None:
                    output["paired_event_shifts"][threshold] = {
                        "delta_a":clock["a"]-ref["a"],"fractional_delta_a":clock["a"]/ref["a"]-1,
                        "delta_t_Myr":1000*(clock["t_Gyr"]-ref["t_Gyr"]),
                        "reference_a":ref["a"],"clock_a":clock["a"],
                        "reference_t_Gyr":ref["t_Gyr"],"clock_t_Gyr":clock["t_Gyr"]}
    except Exception as error:
        failure("main_execution",error)

    try:
        eds,initial,analytic = cycloid_setup()
        solution = integrate_sphere(eds,initial)
        independent = integrate_sphere(eds,initial,coordinate="physical_radius",method="Radau",
                                       rtol=2e-11,atol=2e-13,max_step=0.005)
        for threshold in ("100","200","500"):
            check("EdS_cycloid_a"+threshold,solution.events[threshold]["a"]/analytic["events"][threshold]["a"]-1)
            check("EdS_physical_R_Radau_a"+threshold,independent.events[threshold]["a"]/analytic["events"][threshold]["a"]-1)
        check("EdS_cycloid_turnaround_a",solution.turnaround["a"]/analytic["turnaround"]["a"]-1)
        check("EdS_cycloid_turnaround_Delta",solution.turnaround["Delta"]/analytic["turnaround"]["Delta"]-1)
        output["benchmarks"]["EdS"] = {"initial":initial,"analytic":analytic,
                                          "events":solution.events,"turnaround":solution.turnaround,
                                          "independent_events":independent.events}
    except Exception as error:
        failure("EdS_benchmark",error)

    output["checks"] = checks
    output["failures"] = [item for item in checks if not item["passed"]]
    output["summary"] = {"checks_total":len(checks),"checks_passed":sum(item["passed"] for item in checks),
                         "failures":len(output["failures"]),"models_completed":len(output["models"])}
    output["execution"]["wall_seconds"] = time.perf_counter()-started
    output["limitations"] = [
        "Spherical smooth-background equations neglect tides, shell crossing and clustered scalar/radiation perturbations.",
        "Delta=100/200/500 are fixed-mean-density collapse events; Delta=200 does not establish virialization.",
        "Only EdS has an analytic cycloid benchmark here; no EdS collapse/virial overdensity is imposed on general backgrounds.",
        "Reference delta_i is calibrated to a prescribed numerical event at a=.5, not an observed halo.",
        "Physical ages inherit the declared prior background and its idealized early tail.",
        "PM comparison must separately measure finite-particle core mass, radial estimator, anisotropy and periodic compensation errors."]
    path=ROOT/"results/sphere_summary.json"
    if path.exists():
        previous=path.read_bytes()
        previous_hash=hashlib.sha256(previous).hexdigest()
        (ROOT/f"results/sphere_summary_previous_{previous_hash[:12]}.json").write_bytes(previous)
    path.write_text(json.dumps(serializable(output),ensure_ascii=False,indent=2,allow_nan=False)+"\n")
    print(json.dumps(output["summary"],ensure_ascii=False))
    return 1 if output["failures"] or len(output["models"])!=2 else 0


if __name__ == "__main__":
    raise SystemExit(run())
