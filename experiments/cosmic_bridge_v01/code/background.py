#!/usr/bin/env python3
"""Radiation+dust+Lambda+canonical exponential-clock reference backgrounds.

Run only after protocol.md has been frozen:
    python experiments/cosmic_bridge_v01/code/background.py

Public API:
    ref = Background(epsilon=0.0)
    bg = Background(epsilon=1e-4)
    bg.evaluate_a(a) / bg.evaluate_n(N) -> vectorized dictionary
    bg.E(a), bg.dlnE_dN(a)
    bg.growth_regular() -> GrowthSolution, normalized D(0.02)=1
    bg.growth_matched(reference=ref) -> GrowthSolution with matched physical p_i
    bg.pm_interpolators() -> dictionary of PCHIP callables, argument a

The regular growth branch neglects radiation/scalar perturbations. It is an
initial-condition convention, not a primordial transfer function. The scalar
is the older canonical exponential-potential candidate, not the fixed-box R.
"""
from __future__ import annotations

import csv
import hashlib
import json
import platform
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import scipy
from scipy.integrate import quad, solve_ivp
from scipy.interpolate import PchipInterpolator
from scipy.optimize import brentq


ROOT = Path(__file__).resolve().parents[1]
H0_KM_S_MPC = 67.4
MPC_KM = 3.0856775814913673e19
YEAR_S = 365.25 * 86400.0
H0_S = H0_KM_S_MPC / MPC_KM
H0_GYR = H0_S * YEAR_S * 1e9
OMEGA_M0 = 0.315
OMEGA_R0 = 9.2e-5
TSTAR_RAD_S = 5.391247e-44
AI = 1e-7
A_PM = 0.02
EPSILONS = (0.0, 1e-4, 1e-2)
EPSILON_LABELS = {0.0:"eps0",1e-4:"eps1e-4",1e-2:"eps1e-2"}
DEFAULT_RTOL = 1e-10
DEFAULT_ATOL = 2e-12
DEFAULT_MAX_STEP = 0.04
RELATIVE_THRESHOLD = 1e-6
CLOSURE_THRESHOLD = 1e-10
PM_GRID_SIZE = 4097


def _scalar_or_array(value):
    value = np.asarray(value)
    return float(value) if value.ndim == 0 else value


def _jsonable(value):
    if isinstance(value, dict):
        return {key: _jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(item) for item in value]
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    return value


def _relative_difference(a, b):
    a, b = np.asarray(a), np.asarray(b)
    if np.any(b == 0):
        raise ValueError("Relative comparison requires nonzero reference values")
    return float(np.max(np.abs(a / b - 1.0)))


class GrowthSolution:
    def __init__(self, background, solution, normalization_a, metadata):
        self.background = background
        self.solution = solution
        self.normalization_a = float(normalization_a)
        self.metadata = metadata
        self.nmin = float(solution.t[0])
        self.nmax = float(solution.t[-1])
        self.log_normalization = float(solution.sol(np.log(normalization_a))[0])

    def evaluate_n(self, n):
        n = np.asarray(n, dtype=float)
        if np.any(n < self.nmin - 2e-12) or np.any(n > self.nmax + 2e-12):
            raise ValueError("Growth evaluation outside its integrated interval")
        logd, f = self.solution.sol(np.clip(n, self.nmin, self.nmax))
        return {"N": _scalar_or_array(n), "a": _scalar_or_array(np.exp(n)),
                "D": _scalar_or_array(np.exp(logd - self.log_normalization)),
                "f": _scalar_or_array(f)}

    def evaluate_a(self, a):
        a = np.asarray(a, dtype=float)
        if np.any(a <= 0):
            raise ValueError("Scale factor must be positive")
        return self.evaluate_n(np.log(a))

    def D(self, a):
        return self.evaluate_a(a)["D"]

    def f(self, a):
        return self.evaluate_a(a)["f"]


class Background:
    """Stable N=ln(a) integration with state (u, q=dchi/dN, y=Ht).

    epsilon=0 omits the scalar entirely; chi, q and scalar w are then None.
    Setting omega_lambda skips the present-day closure shoot for benchmarks.
    The declared tstar_rad fixes Vs/Fchi^2=1/(4*tstar_rad^2) in every run.
    """
    def __init__(self, epsilon=1e-4, ai=AI, omega_m0=OMEGA_M0,
                 omega_r0=OMEGA_R0, omega_lambda=None,
                 tstar_rad_s=TSTAR_RAD_S, rtol=DEFAULT_RTOL,
                 atol=DEFAULT_ATOL, max_step=DEFAULT_MAX_STEP):
        self.epsilon = float(epsilon)
        self.ai, self.ni = float(ai), float(np.log(ai))
        self.omega_m0, self.omega_r0 = float(omega_m0), float(omega_r0)
        self.tstar_rad_s = float(tstar_rad_s)
        self.rtol, self.atol, self.max_step = float(rtol), float(atol), float(max_step)
        if not (0 <= self.epsilon < 1 and 0 < self.ai < 1):
            raise ValueError("Require 0 <= epsilon < 1 and 0 < ai < 1")
        if min(self.omega_m0, self.omega_r0) < 0:
            raise ValueError("Negative reference fluid density")
        self.log_potential_amplitude = -np.log(4.0) - 2*np.log(H0_S*self.tstar_rad_s)
        self._set_initial_asymptotics()
        self.closure_history = []
        if omega_lambda is None:
            self.omega_lambda = float(brentq(self._closure, 0.0, 0.9,
                                            xtol=2e-13, rtol=2e-13))
            self.closure_was_shot = True
        else:
            self.omega_lambda = float(omega_lambda)
            self.closure_was_shot = False
        self.solution = self._integrate(self.omega_lambda)
        self._regular_growth = None

    def _set_initial_asymptotics(self):
        eps = self.epsilon
        if self.omega_r0 > 0:
            self.initial_branch = "radiation tracking with first-order dust correction"
            self.initial_x = self.omega_m0*self.ai/self.omega_r0
            if self.initial_x >= 0.01:
                raise ValueError("Radiation asymptotic start requires matter/radiation < 0.01")
            amplitude = self.omega_r0/(1-eps)
            coefficient = (1-eps)/(1-2*eps/3)
            correction = coefficient*self.initial_x
            self.q_i = 2-correction/2
            self.W_i = amplitude*self.ai**-4*(1+correction)
            self.tau_i = self.ai**2/(2*np.sqrt(amplitude))*(1-correction/3)
            self.tail_order = "O((Omega_m0*ai/Omega_r0)^2); Lambda omitted in early tail"
        elif self.omega_m0 > 0:
            self.initial_branch = "exact matter tracking when Lambda=0"
            self.initial_x = None
            initial_e2 = self.omega_m0*self.ai**-3/(1-3*eps/4)
            self.q_i = 1.5
            self.W_i = 9*initial_e2/8
            self.tau_i = 2/(3*np.sqrt(initial_e2))
            self.tail_order = "exact pure-matter tracking tail; Lambda omitted in early tail"
        else:
            raise ValueError("At least one reference radiation or dust component is required")
        if eps > 0:
            self.chi_i = 0.5*(self.log_potential_amplitude-np.log(self.W_i))
            if self.chi_i <= 0:
                raise ValueError("Chosen response branch requires chi_i > 0")
        else:
            self.chi_i = None

    def _terms(self, n, u, q, omega_lambda):
        n = np.asarray(n)
        radiation = self.omega_r0*np.exp(-4*n)
        matter = self.omega_m0*np.exp(-3*n)
        if self.epsilon > 0:
            potential = self.W_i*np.exp(-2*np.asarray(u))
            denominator = 1-self.epsilon*np.asarray(q)**2/6
            if np.any(denominator <= 0):
                raise FloatingPointError("Friedmann kinetic denominator crossed zero")
            e2 = (radiation+matter+omega_lambda+self.epsilon*potential/3)/denominator
            scalar_kinetic_fraction = self.epsilon*np.asarray(q)**2/2
        else:
            potential = np.zeros_like(n, dtype=float)
            denominator = np.ones_like(n, dtype=float)
            e2 = radiation+matter+omega_lambda
            scalar_kinetic_fraction = 0.0
        if np.any(e2 <= 0) or np.any(~np.isfinite(e2)):
            raise FloatingPointError("Nonpositive or nonfinite Friedmann value")
        dln_e = -2*radiation/e2-1.5*matter/e2-scalar_kinetic_fraction
        return e2, dln_e, potential, denominator, radiation/e2, matter/e2

    def _integrate(self, omega_lambda):
        if self.epsilon > 0:
            initial_e2 = self._terms(self.ni, 0.0, self.q_i, omega_lambda)[0]
            initial = (0.0, self.q_i, np.sqrt(initial_e2)*self.tau_i)

            def rhs(n, state):
                u, q, y = state
                e2, f, potential, *_ = self._terms(n,u,q,omega_lambda)
                return q, -(3+f)*q+2*potential/e2, 1+f*y
        else:
            initial_e2 = self._terms(self.ni, 0.0, 0.0, omega_lambda)[0]
            initial = (np.sqrt(initial_e2)*self.tau_i,)

            def rhs(n, state):
                f = self._terms(n,0.0,0.0,omega_lambda)[1]
                return (1+f*state[0],)
        result = solve_ivp(rhs, (self.ni,0.0), initial, method="DOP853",
                           rtol=self.rtol, atol=self.atol, max_step=self.max_step,
                           dense_output=True)
        if not result.success:
            raise RuntimeError(result.message)
        return result

    def _closure(self, omega_lambda):
        if self.epsilon == 0:
            value = self.omega_r0+self.omega_m0+omega_lambda-1.0
        else:
            solution = self._integrate(omega_lambda)
            u, q, _ = solution.y[:,-1]
            value = float(self._terms(0.0,u,q,omega_lambda)[0]-1.0)
        self.closure_history.append({"omega_lambda":float(omega_lambda),"E0_squared_minus_one":value})
        return value

    def evaluate_n(self, n):
        n = np.asarray(n,dtype=float)
        if np.any(n < self.ni-2e-12) or np.any(n > 2e-12):
            raise ValueError("Background evaluation outside its integrated interval")
        n = np.clip(n,self.ni,0.0)
        state = self.solution.sol(n)
        if self.epsilon > 0:
            u,q,y = state
        else:
            u,q,y = np.zeros_like(n),np.zeros_like(n),state[0]
        e2,f,w,denominator,orad,omat = self._terms(n,u,q,self.omega_lambda)
        e = np.sqrt(e2)
        potential_over_e2 = w/e2
        scalar_energy_over_e2 = q*q/2+potential_over_e2
        output = {
            "N":n,"a":np.exp(n),"E":e,"E2":e2,"dlnE":f,"Ht":y,
            "tau":y/e,"t_Gyr":y/e/H0_GYR,
            "Omega_r":orad,"Omega_m":omat,"Omega_lambda":self.omega_lambda/e2,
            "Omega_chi":self.epsilon*scalar_energy_over_e2/3,
            "kinetic_denominator":denominator,
        }
        output = {key:_scalar_or_array(value) for key,value in output.items()}
        if self.epsilon > 0:
            output.update(chi=_scalar_or_array(self.chi_i+u),u=_scalar_or_array(u),
                          q=_scalar_or_array(q),W=_scalar_or_array(w),
                          w_chi=_scalar_or_array((q*q/2-potential_over_e2)/scalar_energy_over_e2),
                          scalar_energy=_scalar_or_array(e2*scalar_energy_over_e2))
        else:
            output.update(chi=None,u=None,q=None,W=None,w_chi=None,scalar_energy=None)
        return output

    def evaluate_a(self, a):
        a = np.asarray(a,dtype=float)
        if np.any(a <= 0):
            raise ValueError("Scale factor must be positive")
        return self.evaluate_n(np.log(a))

    def E(self, a):
        return self.evaluate_a(a)["E"]

    def dlnE_dN(self, a):
        return self.evaluate_a(a)["dlnE"]

    def _growth(self, n_initial, logd_initial, f_initial, normalization_a, metadata):
        def rhs(n,state):
            values = self.evaluate_n(n)
            f_growth = state[1]
            return f_growth, 1.5*values["Omega_m"]-f_growth*f_growth-(2+values["dlnE"])*f_growth
        solution = solve_ivp(rhs,(n_initial,0.0),(logd_initial,f_initial),
                             method="DOP853",rtol=self.rtol,atol=self.atol,
                             max_step=self.max_step,dense_output=True)
        if not solution.success:
            raise RuntimeError(solution.message)
        return GrowthSolution(self,solution,normalization_a,metadata)

    def growth_regular(self, normalization_a=A_PM):
        if self._regular_growth is not None and self._regular_growth.normalization_a == normalization_a:
            return self._regular_growth
        if self.omega_r0 > 0:
            correction = 1.5*(1-self.epsilon)*self.initial_x
            logd_i = np.log1p(correction)
            f_i = correction/(1+correction)
            branch = "regular smooth-radiation dust-growth branch; not primordial transfer"
        else:
            f_i = (np.sqrt(25-18*self.epsilon)-1)/4
            logd_i = 0.0
            branch = "matter-tracking growing power law"
        metadata = {"branch":branch,"initial_a":self.ai,"initial_f":float(f_i),
                    "normalization_a":float(normalization_a),
                    "radiation_and_scalar_perturbations":"omitted; physical use requires subhorizon smooth-field regime"}
        growth = self._growth(self.ni,logd_i,f_i,normalization_a,metadata)
        if normalization_a == A_PM:
            self._regular_growth = growth
        return growth

    def growth_matched(self, reference=None, f_initial=None, a_initial=A_PM):
        if (reference is None) == (f_initial is None):
            raise ValueError("Specify exactly one of reference or f_initial")
        metadata = {"branch":"late dust growth with matched initial D and physical canonical momentum",
                    "initial_a":float(a_initial),"normalization_a":float(a_initial)}
        if reference is not None:
            reference_f = float(reference.growth_regular().f(a_initial))
            reference_e = float(reference.E(a_initial))
            model_e = float(self.E(a_initial))
            f_initial = reference_f*reference_e/model_e
            metadata.update(reference_epsilon=reference.epsilon,reference_f_initial=reference_f,
                            reference_E_initial=reference_e,model_E_initial=model_e,
                            matching_relation="f_model_i * E_model_i = f_reference_i * E_reference_i")
        metadata["initial_f"] = float(f_initial)
        return self._growth(float(np.log(a_initial)),0.0,float(f_initial),a_initial,metadata)

    def pm_interpolators(self, a_initial=A_PM, count=PM_GRID_SIZE):
        """Return named callables E(a), dlnE_dN(a), Ht(a), t_Gyr(a)."""
        n = np.linspace(np.log(a_initial),0.0,int(count))
        values = self.evaluate_n(n)
        result = {"N_grid":n,"a_grid":np.exp(n)}
        for output_name,input_name in (("E","E"),("dlnE_dN","dlnE"),("Ht","Ht"),("t_Gyr","t_Gyr")):
            interpolation = PchipInterpolator(n,values[input_name],extrapolate=False)
            result[output_name] = lambda a, interpolation=interpolation: interpolation(np.log(a))
        return result

    def transfer_chi_Gyr(self, a):
        if self.epsilon == 0:
            raise ValueError("No physical clock response is assigned to epsilon=0 control")
        values, today = self.evaluate_a(a), self.evaluate_a(1.0)
        delta = np.asarray(values["u"])-today["u"]
        return _scalar_or_array(-today["chi"]/(2*H0_GYR*today["q"])
                                *np.expm1(-2*np.log1p(delta/today["chi"])))


def run():
    protocol_path = ROOT/"protocol.md"
    if not protocol_path.exists():
        raise RuntimeError("Freeze protocol.md before executing background checks")
    output_dir = ROOT/"results"
    output_dir.mkdir(exist_ok=True)
    checks = []

    def record(name, value, threshold=RELATIVE_THRESHOLD, **details):
        finite = bool(np.isfinite(value))
        checks.append({"name":name,"value":float(value),"threshold":float(threshold),
                       "passed":finite and float(value)<=threshold,**_jsonable(details)})

    def exception(name, error):
        checks.append({"name":name,"passed":False,"exception":repr(error)})

    models, matched, rows, model_summaries = {},{},[],[]
    data = {"metadata":{
        "run_utc":datetime.now(timezone.utc).isoformat(),
        "python":platform.python_version(),"numpy":np.__version__,"scipy":scipy.__version__,
        "script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "protocol_sha256":hashlib.sha256(protocol_path.read_bytes()).hexdigest(),
        "scope":"declared reference backgrounds and conditional dust growth; no observation fitting",
    },"inputs":{
        "H0_km_s_Mpc":H0_KM_S_MPC,"Omega_m0":OMEGA_M0,"Omega_r0":OMEGA_R0,
        "epsilons":EPSILONS,"epsilon_1e_minus_2":"diagnostic only; not asserted observationally allowed",
        "initial_a":AI,"pm_initial_a":A_PM,"pm_grid_count":PM_GRID_SIZE,
        "tstar_rad_seconds":TSTAR_RAD_S,
        "tstar_role":"declared potential/response anchor; not a derived fundamental time",
        "potential_ratio":"Vs/Fchi^2 = 1/(4*tstar_rad^2), held fixed across initial-a checks",
        "rtol":DEFAULT_RTOL,"atol":DEFAULT_ATOL,"max_step_N":DEFAULT_MAX_STEP,
        "comparison_threshold":RELATIVE_THRESHOLD,"closure_threshold":CLOSURE_THRESHOLD,
        "growth_matching":"same initial D and physical p; f_i is rescaled by E_ref/E_model",
    }}
    for eps in EPSILONS:
        label = f"eps_{eps:g}"
        try:
            models[eps] = Background(epsilon=eps)
        except Exception as error:
            exception(label+"_construction",error)

    reference = models.get(0.0)
    n_full = np.linspace(np.log(AI),0.0,4097)
    n_pm = np.linspace(np.log(A_PM),0.0,PM_GRID_SIZE)
    n_table = np.unique(np.r_[n_full,n_pm])
    for eps,bg in models.items():
        label = f"eps_{eps:g}"
        try:
            values = bg.evaluate_n(n_table)
            today = bg.evaluate_a(1.0)
            growth = bg.growth_regular()
            gv = growth.evaluate_n(n_table)
            record(label+"_present_Friedmann_closure",abs(today["E"]-1),CLOSURE_THRESHOLD)
            fraction_sum = values["Omega_r"]+values["Omega_m"]+values["Omega_lambda"]+values["Omega_chi"]
            record(label+"_density_sum",float(np.max(np.abs(fraction_sum-1))),1e-10)
            record(label+"_positive_background",0 if np.all(values["E"]>0) and np.all(values["Ht"]>0) else 1,0)
            record(label+"_positive_kinetic_denominator",0 if np.all(values["kinetic_denominator"]>0) else 1,0)
            if eps > 0:
                record(label+"_positive_chi",0 if np.all(values["chi"]>0) else 1,0)
                record(label+"_canonical_scalar_w",max(0.0,float(np.max(np.abs(values["w_chi"])-1))),1e-10)

            age_comparisons = []
            for n in np.linspace(bg.ni,0.0,19):
                integral = 0.0 if n == bg.ni else quad(lambda x:1/bg.evaluate_n(x)["E"],bg.ni,n,
                                                     epsabs=1e-300,epsrel=3e-11,limit=200)[0]
                independent_tau = bg.tau_i+integral
                actual_tau = bg.evaluate_n(n)["tau"]
                age_comparisons.append({"N":float(n),"tau_state":actual_tau,"tau_quadrature":independent_tau,
                                        "relative_error":abs(actual_tau/independent_tau-1)})
            record(label+"_independent_age_integral",max(item["relative_error"] for item in age_comparisons),
                   comparisons=age_comparisons)

            nd = np.linspace(bg.ni+0.002,-0.002,257)
            dh = 2e-4
            sampled = [np.asarray(bg.evaluate_n(nd+j*dh)["E2"]) for j in (-2,-1,1,2)]
            derivative = (sampled[0]-8*sampled[1]+8*sampled[2]-sampled[3])/(12*dh)
            center = bg.evaluate_n(nd)
            expected = -4*bg.omega_r0*np.exp(-4*nd)-3*bg.omega_m0*np.exp(-3*nd)
            if eps > 0:
                expected -= eps*center["E2"]*center["q"]**2
            record(label+"_independent_Friedmann_derivative",_relative_difference(derivative,expected))

            if eps > 0:
                edges = np.linspace(bg.ni,0.0,13)
                energy_errors = []
                for left,right in zip(edges[:-1],edges[1:]):
                    change = bg.evaluate_n(right)["scalar_energy"]-bg.evaluate_n(left)["scalar_energy"]
                    loss = quad(lambda x:-3*bg.evaluate_n(x)["E2"]*bg.evaluate_n(x)["q"]**2,
                                left,right,epsabs=1e-300,epsrel=3e-10,limit=200)[0]
                    energy_errors.append(abs(change/loss-1))
                record(label+"_integrated_scalar_energy",max(energy_errors),interval_relative_errors=energy_errors)

            initial_momentum_comparison = None
            matched_values = None
            if reference is not None:
                matched[eps] = bg.growth_matched(reference=reference)
                matched_values = matched[eps].evaluate_n(n_pm)
                initial_momentum_comparison = (matched[eps].f(A_PM)*bg.E(A_PM)
                                                /(reference.growth_regular().f(A_PM)*reference.E(A_PM))-1)
                record(label+"_matched_initial_physical_momentum",abs(initial_momentum_comparison),1e-12)
                record(label+"_matched_initial_D",abs(matched[eps].D(A_PM)-1),1e-12)
            else:
                exception(label+"_matched_growth",RuntimeError("epsilon=0 reference unavailable"))

            interpolated = bg.pm_interpolators()
            n_mid = (n_pm[:-1]+n_pm[1:])/2
            mids = bg.evaluate_n(n_mid)
            for field,function in (("E","E"),("dlnE","dlnE_dN"),("t_Gyr","t_Gyr")):
                record(label+"_PM_interpolation_"+field,
                       _relative_difference(interpolated[function](np.exp(n_mid)),mids[field]))

            summary = {"epsilon":eps,"omega_lambda":bg.omega_lambda,"today":today,
                       "initial_chi":bg.chi_i,"initial_W":bg.W_i if eps>0 else None,
                       "early_tail_tau":bg.tau_i,"initial_branch":bg.initial_branch,
                       "early_tail_order":bg.tail_order,"closure_history":bg.closure_history,
                       "regular_growth":growth.metadata,"regular_D_today_normalized_at_PM":growth.D(1.0),
                       "regular_f_PM":growth.f(A_PM),"regular_f_today":growth.f(1.0),
                       "matched_growth":matched[eps].metadata if eps in matched else None,
                       "matched_D_today":matched[eps].D(1.0) if eps in matched else None,
                       "initial_pm":{"a":A_PM,"E":bg.E(A_PM),"f_regular":growth.f(A_PM),
                                     "f_matched":matched[eps].f(A_PM) if eps in matched else None},
                       "age_samples":[{"a":a,**bg.evaluate_a(a)} for a in (AI,1/1101,1/101,A_PM,0.5,1.0)]}
            model_summaries.append(summary)
            for i,n in enumerate(n_table):
                row = {"epsilon":eps}
                for field in ("a","N","E","dlnE","Ht","t_Gyr","chi","q","Omega_chi","w_chi","Omega_m","Omega_r"):
                    row[field] = None if values[field] is None else float(values[field][i])
                row.update(D_regular_norm_PM=float(gv["D"][i]),f_regular=float(gv["f"][i]),
                           D_matched=None,f_matched=None)
                if n >= np.log(A_PM)-1e-12 and eps in matched:
                    late = matched[eps].evaluate_n(n)
                    row.update(D_matched=late["D"],f_matched=late["f"])
                rows.append(row)
        except Exception as error:
            exception(label+"_checks_or_output",error)

        # Keep the declared physical potential ratio fixed in both alternatives.
        for alternative_name,kwargs in (
                ("tighter_tolerance",{"rtol":2e-12,"atol":2e-14,"max_step":0.02}),
                ("earlier_start",{"ai":1e-8})):
            try:
                alternate = Background(epsilon=eps,**kwargs)
                na = np.linspace(bg.ni,0.0,1025)
                base,other = bg.evaluate_n(na),alternate.evaluate_n(na)
                for field in ("E","t_Gyr"):
                    record(label+"_"+alternative_name+"_"+field,_relative_difference(other[field],base[field]))
                if eps > 0:
                    record(label+"_"+alternative_name+"_chi",_relative_difference(other["chi"],base["chi"]))
                ga,gb = alternate.growth_regular().evaluate_n(n_pm),bg.growth_regular().evaluate_n(n_pm)
                for field in ("D","f"):
                    record(label+"_"+alternative_name+"_growth_"+field,_relative_difference(ga[field],gb[field]))
                if reference is not None:
                    # The comparison holds the same reference physical p_i fixed.
                    ma = alternate.growth_matched(reference=reference).evaluate_n(n_pm)
                    mb = bg.growth_matched(reference=reference).evaluate_n(n_pm)
                    for field in ("D","f"):
                        record(label+"_"+alternative_name+"_matched_"+field,_relative_difference(ma[field],mb[field]))
            except Exception as error:
                exception(label+"_"+alternative_name,error)

    benchmark_summaries = []
    for eps in (1e-4,1e-2):
        for branch,om,orad in (("radiation",0.0,OMEGA_R0),("matter",OMEGA_M0,0.0)):
            name = f"pure_{branch}_eps_{eps:g}"
            try:
                bg = Background(epsilon=eps,omega_m0=om,omega_r0=orad,omega_lambda=0.0)
                n = np.linspace(bg.ni,0.0,1025)
                values = bg.evaluate_n(n)
                if branch == "radiation":
                    q_expected,y_expected,omega_expected,w_expected = 2.0,0.5,eps,1/3
                    e_expected = np.sqrt(orad/(1-eps))*np.exp(-2*n)
                    tstar = TSTAR_RAD_S
                else:
                    q_expected,y_expected,omega_expected,w_expected = 1.5,2/3,0.75*eps,0.0
                    e_expected = np.sqrt(om/(1-0.75*eps))*np.exp(-1.5*n)
                    tstar = np.sqrt(2)*TSTAR_RAD_S
                for field,expected in (("q",q_expected),("Ht",y_expected),("Omega_chi",omega_expected)):
                    record(name+"_"+field,_relative_difference(values[field],expected))
                record(name+"_w_chi_absolute",float(np.max(np.abs(values["w_chi"]-w_expected))))
                record(name+"_E",_relative_difference(values["E"],e_expected))
                log_time = np.log(values["tau"]/(H0_S*tstar))
                record(name+"_log_time_absolute",float(np.max(np.abs(values["chi"]-log_time))))
                if branch == "matter":
                    exponent = (np.sqrt(25-18*eps)-1)/4
                    growth = bg.growth_regular().evaluate_n(n_pm)
                    record(name+"_growth_power",_relative_difference(growth["D"],(np.exp(n_pm)/A_PM)**exponent))
                benchmark_summaries.append({"name":name,"expected_q":q_expected,"expected_Ht":y_expected,
                                            "expected_Omega_chi":omega_expected,"tstar_seconds":tstar})
            except Exception as error:
                exception(name,error)

    for name,om,orad in (("Meszaros",OMEGA_M0,OMEGA_R0),("EdS",1.0,0.0)):
        try:
            bg = Background(epsilon=0.0,omega_m0=om,omega_r0=orad,omega_lambda=0.0)
            growth = bg.growth_regular().evaluate_n(n_pm)
            if orad > 0:
                expected = (1+1.5*om*np.exp(n_pm)/orad)/(1+1.5*om*A_PM/orad)
                f_expected = (1.5*om*np.exp(n_pm)/orad)/(1+1.5*om*np.exp(n_pm)/orad)
            else:
                expected = np.exp(n_pm)/A_PM
                f_expected = np.ones_like(n_pm)
            record(name+"_growth_D",_relative_difference(growth["D"],expected))
            record(name+"_growth_f",_relative_difference(growth["f"],f_expected))
        except Exception as error:
            exception(name,error)

    data["models"] = {EPSILON_LABELS[item["epsilon"]]:item for item in model_summaries}
    data["analytic_benchmarks"] = benchmark_summaries
    data["checks"] = checks
    data["failures"] = [item for item in checks if not item["passed"]]
    data["summary"] = {"total_checks":len(checks),"passed":sum(item["passed"] for item in checks),
                       "failed":len(data["failures"]),"background_models_completed":len(model_summaries),
                       "grid_rows":len(rows),"all_requested_models_completed":len(model_summaries)==len(EPSILONS)}
    with (output_dir/"background_grid.csv").open("w",newline="") as stream:
        if rows:
            writer = csv.DictWriter(stream,fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
    for eps,bg in models.items():
        if eps not in matched:
            continue
        values = bg.evaluate_n(n_pm)
        regular_values = bg.growth_regular().evaluate_n(n_pm)
        matched_values = matched[eps].evaluate_n(n_pm)
        pm_rows = []
        for i in range(len(n_pm)):
            row = {"epsilon":eps}
            for field in ("a","N","E","dlnE","Ht","t_Gyr","chi","q","Omega_chi","w_chi","Omega_m","Omega_r"):
                row[field] = None if values[field] is None else float(values[field][i])
            row.update(D_regular_norm_PM=float(regular_values["D"][i]),
                       f_regular=float(regular_values["f"][i]),
                       D_matched=float(matched_values["D"][i]),
                       f_matched=float(matched_values["f"][i]))
            pm_rows.append(row)
        with (output_dir/f"background_{EPSILON_LABELS[eps]}.csv").open("w",newline="") as stream:
            writer = csv.DictWriter(stream,fieldnames=list(pm_rows[0]))
            writer.writeheader()
            writer.writerows(pm_rows)
    with (output_dir/"background_summary.json").open("w") as stream:
        json.dump(_jsonable(data),stream,ensure_ascii=False,indent=2,allow_nan=False)
        stream.write("\n")
    print(json.dumps(data["summary"],ensure_ascii=False))
    return 1 if data["failures"] or not data["summary"]["all_requested_models_completed"] else 0


if __name__ == "__main__":
    raise SystemExit(run())
