#!/usr/bin/env python3
"""Independent proper-time and direct-K/W audit of the SUGRA clock experiment.

No production implementation enters the independent differential equations.
Backgrounds are solved against tau=H_ref*t with physical velocities, rather than
the production N equations.  The production potential is loaded only after all
independent integrations, to compare it with direct high-precision K/W derivatives.
The one-dimensional monotone N(tau) map is inverted only for comparison.  The
real potential used in the ODE is symbolically differentiated here; an additional
mpmath audit reconstructs it directly from complex K and W before differentiating.
"""
from __future__ import annotations

import csv
import hashlib
import importlib.util
import json
import math
from pathlib import Path

import mpmath as mp
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
from scipy.optimize import brentq


ROOT = Path(__file__).resolve().parents[1]
MODELS = ("global", "canonical", "shift")
FIELDS = ("V", "Vchi", "Vy", "Vchichi", "Vchiy", "Vyy")


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def symbolic_potentials():
    """Own derivation: covariant W derivatives, then ordinary real derivatives."""
    x, y, ci, eps, wi = sp.symbols("x y ci eps wi", real=True)
    tree = wi * sp.exp(-2 * (x - ci))
    # These are |D_Z W|^2 / U and K/Mpl^2, with the same holomorphic W.
    dz2 = {
        "global": sp.Integer(1),
        "canonical": (1 - eps * x / 2) ** 2 + eps**2 * y**2 / 4,
        "shift": 1 + eps**2 * y**2,
    }
    kahler = {"global": 0, "canonical": eps * (x*x + y*y) / 2,
              "shift": eps * y*y}
    out = {}
    for name in MODELS:
        v = tree if name == "global" else tree * sp.exp(kahler[name]) * (
            dz2[name] - 3 * eps / 2)
        exprs = [v, sp.diff(v, x), sp.diff(v, y), sp.diff(v, x, 2),
                 sp.diff(v, x, y), sp.diff(v, y, 2)]
        out[name] = sp.lambdify((x, y, ci, eps, wi), exprs, "numpy", cse=True)
    return out


POTENTIALS = symbolic_potentials()


class ProperTimeModel:
    def __init__(self, pars, name):
        self.p = pars
        self.name = name
        self.ci = float(pars["chi_i"])
        self.eps = float(pars["epsilon"])
        self.wi = float(pars["W_i"])

    def potential(self, u, y):
        vals = POTENTIALS[self.name](self.ci + u, y, self.ci, self.eps, self.wi)
        return dict(zip(FIELDS, map(float, vals)))

    def hubble(self, n, u, y, wx, wy):
        v = self.potential(u, y)["V"]
        fluid = (self.p["omega_r"] * math.exp(-4*n)
                 + self.p["omega_m"] * math.exp(-3*n)
                 + self.p["omega_lambda"])
        e2 = fluid + self.eps * (0.5 * (wx*wx + wy*wy) + v) / 3
        if not np.isfinite(e2) or e2 <= 0:
            raise FloatingPointError("Nonpositive/nonfinite proper-time Friedmann H^2")
        return math.sqrt(e2)

    def rhs(self, tau, state):
        n, u, y, wx, wy, loss = state
        v = self.potential(u, y)
        e = self.hubble(n, u, y, wx, wy)
        return [e, wx, wy, -3*e*wx-v["Vchi"], -3*e*wy-v["Vy"],
                3*e*(wx*wx+wy*wy)]

    def solve(self, yi):
        p = self.p
        ini = [p["N_i"], 0.0, yi, p["wchi_i"], 0.0, 0.0]

        def today(t, state):
            return state[0]

        def displacement_domain(t, state):
            return 20-abs(state[1])

        def transverse_domain(t, state):
            return 1-abs(state[2])

        for ev in (today, displacement_domain, transverse_domain):
            ev.terminal = True
        today.direction = 1
        displacement_domain.direction = -1
        transverse_domain.direction = -1
        sol = solve_ivp(self.rhs, (p["tau_i"], p["tau_i"]+3), ini,
                        method="DOP853", rtol=2e-12, atol=2e-14,
                        max_step=5e-4, dense_output=True,
                        events=(today, displacement_domain, transverse_domain))
        if not sol.success or len(sol.t_events[0]) != 1:
            raise RuntimeError(f"{self.name}, yi={yi}: {sol.message}, events={sol.t_events}")
        return sol

    def invert_n(self, sol, n_grid):
        lo, hi = sol.t[0], sol.t[-1]
        values = []
        for n in n_grid:
            if abs(n-sol.y[0, 0]) < 2e-13:
                tau = lo
            elif abs(n) < 2e-13:
                tau = hi
            else:
                tau = brentq(lambda t: sol.sol(t)[0]-n, lo, hi,
                             xtol=3e-15, rtol=9e-16)
            values.append(tau)
        return np.array(values)

    def modes(self, background, k):
        """Both bases are integrated in cosmic time, no N-friction equations."""
        e_i = self.hubble(*background.y[:5, 0])
        ini = [1.0, 0.0, 0.0, e_i]

        def rhs(tau, state):
            n, u, y, wx, wy = background.sol(tau)[:5]
            e = self.hubble(n, u, y, wx, wy)
            omega2 = k*k*math.exp(-2*n)+self.potential(u, y)["Vyy"]
            return [state[1], -3*e*state[1]-omega2*state[0],
                    state[3], -3*e*state[3]-omega2*state[2]]

        result = solve_ivp(rhs, (background.t[0], background.t[-1]), ini,
                           method="DOP853", rtol=2e-12, atol=2e-14,
                           max_step=5e-4, dense_output=True)
        if not result.success:
            raise RuntimeError(result.message)
        return result


def direct_kw_potential(name, chi, y, pars):
    """Dimensionless V/F^2 H_ref^2 from K^{Z Zbar} D W D Wbar-3|W|^2.

    z=Z/F, k=K/Mpl^2, w=W/(F sqrt(U_i)).  z and zbar are independent
    during the Wirtinger derivatives.  Conjugacy is imposed after differentiation.
    This uses no simplified real potential or derivative from the ODE.
    """
    eps, ci, wi = [mp.mpf(str(pars[key])) for key in ("epsilon", "chi_i", "W_i")]
    z = (chi+1j*y)/mp.sqrt(2)
    zb = mp.conj(z)

    def w(zz):
        return -mp.exp(-mp.sqrt(2)*(zz-ci/mp.sqrt(2)))/mp.sqrt(2)

    dw = mp.diff(w, z)
    if name == "global":
        return wi * abs(dw)**2

    def kahler(zz, zbc):
        if name == "canonical":
            return eps*zz*zbc
        if name == "shift":
            return -eps*(zz-zbc)**2/2
        raise ValueError(name)

    k_z = mp.diff(lambda zz: kahler(zz, zb), z)
    k_zzbar = mp.diff(lambda zzb: mp.diff(lambda zz: kahler(zz, zzb), z), zb)
    cov_dw = dw+k_z*w(z)
    value = wi*mp.exp(kahler(z, zb)) * (
        eps/k_zzbar*abs(cov_dw)**2 - 3*eps*abs(w(z))**2)
    if abs(mp.im(value)) > mp.mpf("1e-50"):
        raise ArithmeticError("Real potential has nontrivial imaginary part")
    return mp.re(value)


def direct_kw_values(name, chi, y, pars):
    f = lambda x, yy: direct_kw_potential(name, x, yy, pars)
    return [f(chi, y), mp.diff(lambda x: f(x, y), chi),
            mp.diff(lambda yy: f(chi, yy), y),
            mp.diff(lambda x: f(x, y), chi, 2),
            mp.diff(lambda x: mp.diff(lambda yy: f(x, yy), y), chi),
            mp.diff(lambda yy: f(chi, yy), y, 2)]


def csv_groups(path, keys):
    groups = {}
    with path.open(newline="") as stream:
        for row in csv.DictReader(stream):
            groups.setdefault(tuple(row[key] for key in keys), []).append(row)
    return groups


def normalized_error(actual, expected, floor=1.0):
    """Maximum absolute error divided by max(floor,max(abs(expected)))."""
    actual, expected = np.asarray(actual), np.asarray(expected)
    scale = max(float(floor), float(np.max(np.abs(expected))))
    return float(np.max(np.abs(actual-expected))/scale)


def main():
    # Archive adapters below are the only coupling to production artifacts.
    pars = json.loads((ROOT/"results/inputs.json").read_text())
    if "physical_wchi_i" in pars:
        pars["wchi_i"] = pars["physical_wchi_i"]
    elif "wchi_i" not in pars:
        eps, q, ni = pars["epsilon"], pars["q_i"], pars["N_i"]
        fluid = (pars["omega_r"]*math.exp(-4*ni)+pars["omega_m"]*math.exp(-3*ni)
                 +pars["omega_lambda"])
        eref2 = (fluid+eps*pars["W_i"]/3)/(1-eps*q*q/6)
        pars["wchi_i"] = q*math.sqrt(eref2)
    groups = csv_groups(ROOT/"results/trajectories.csv", ["model", "y_initial"])
    checks, output, backgrounds = [], [], {}
    for (name, yi_text), rows in groups.items():
        model = ProperTimeModel(pars, name)
        yi = float(yi_text)
        sol = model.solve(yi)
        n = np.array([float(row["N"]) for row in rows])
        tau = model.invert_n(sol, n)
        states = sol.sol(tau)
        e = np.array([model.hubble(*s[:5]) for s in states.T])
        candidate = {"u": states[1], "y": states[2], "wchi": states[3],
                     "wy": states[4], "E": e, "tau": tau}
        for field, values in candidate.items():
            if field == "wchi":
                expected = np.array([float(row["qx"])*float(row["E"]) for row in rows])
            elif field == "wy":
                expected = np.array([float(row["qy"])*float(row["E"]) for row in rows])
            else:
                expected = np.array([float(row[field]) for row in rows])
            err = normalized_error(values, expected)
            checks.append({"kind": "proper_time_background", "model": name,
                           "y_i": yi, "field": field,
                           "error_definition": "max_abs_difference/max(1,max_abs_primary)",
                           "normalized_max_error": err, "limit": 1e-6, "passed": err<1e-6})
        rho = .5*(states[3]**2+states[4]**2)+np.array([
            model.potential(uu, yy)["V"] for uu, yy in states[[1, 2]].T])
        balance = rho+states[5]-rho[0]
        rel_balance = float(max(abs(balance))/max(1, abs(rho[0])))
        checks.append({"kind": "proper_time_energy_identity", "model": name,
                       "y_i": yi, "normalized_max_error": rel_balance,
                       "limit": 1e-8, "passed": rel_balance<1e-8})
        for j in range(len(n)):
            output.append({"model": name, "y_i": yi, "N": n[j],
                           **{key: val[j] for key, val in candidate.items()}})
        if yi == 0:
            backgrounds[name] = (model, sol)

    mode_output = []
    mode_groups = csv_groups(ROOT/"results/linear_modes.csv", ["model", "kappa"])
    for (name, k_text), rows in mode_groups.items():
        model, background = backgrounds[name]
        k = float(k_text)
        mode = model.modes(background, k)
        n = np.array([float(row["N"]) for row in rows])
        tau = model.invert_n(background, n)
        vals = mode.sol(tau)
        bg = background.sol(tau)
        e = np.array([model.hubble(*state[:5]) for state in bg.T])
        candidate = {"y1": vals[0], "q1": vals[1]/e,
                     "y2": vals[2], "q2": vals[3]/e}
        field_map = {"y1": "D0", "q1": "D0_N", "y2": "D1", "q2": "D1_N"}
        for key, value in candidate.items():
            reference = np.array([float(row[field_map[key]]) for row in rows])
            err = normalized_error(value, reference)
            checks.append({"kind": "proper_time_linear_modes", "model": name,
                           "k_Href": k, "field": key,
                           "error_definition": "max_abs_difference/max(1,max_abs_primary)",
                           "normalized_max_error": err, "limit": 1e-6, "passed": err<1e-6})
        for j in range(len(n)):
            mode_output.append({"model": name, "k": k, "N": n[j],
                                **{key: val[j] for key, val in candidate.items()}})

    mp.mp.dps = 70
    potential_output = []
    spec = importlib.util.spec_from_file_location("production_potential_for_comparison_only",
                                                 ROOT/"code/run_sugra.py")
    production = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(production)
    primary_fields = dict(zip(FIELDS, ("V", "Vx", "Vy", "Vxx", "Vxy", "Vyy")))
    # All points are preset and include off-axis values, not selected for success.
    for name in MODELS:
        model = ProperTimeModel(pars, name)
        for u in (0.0, .75, 2.0):
            for yy in (0.0, .001, .1):
                chi, y = mp.mpf(str(pars["chi_i"]))+mp.mpf(str(u)), mp.mpf(str(yy))
                direct = direct_kw_values(name, chi, y, pars)
                symbolic = model.potential(u, yy)
                primary = production.potential(name, float(chi), float(y), pars)
                for field, exact in zip(FIELDS, direct):
                    actual = symbolic[field]
                    err = float(abs(mp.mpf(actual)-exact)/max(mp.mpf(1), abs(exact)))
                    record = {"kind": "direct_complex_KW_vs_symbolic_potential", "model": name,
                              "u": u, "y": yy, "field": field,
                              "direct_high_precision": mp.nstr(exact, 55),
                              "symbolic_float": actual,
                              "normalized_max_error": err, "limit": 2e-11,
                              "passed": err<2e-11}
                    checks.append(record)
                    potential_output.append(record)
                    primary_error = float(abs(mp.mpf(float(primary[primary_fields[field]]))-exact)/max(
                        mp.mpf(1), abs(exact)))
                    checks.append({"kind": "direct_complex_KW_vs_primary_potential",
                                   "model": name, "u": u, "y": yy, "field": field,
                                   "normalized_max_error": primary_error,
                                   "limit": 2e-11, "passed": primary_error<2e-11})

    for filename, rows in (("independent_trajectories.csv", output),
                           ("independent_linear_modes.csv", mode_output)):
        with (ROOT/"results"/filename).open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
    summary = {"protocol_sha256": sha256(ROOT/"protocol.md"),
               "independent_code_sha256": sha256(Path(__file__)),
               "inputs_sha256": sha256(ROOT/"results/inputs.json"),
               "primary_trajectories_sha256": sha256(ROOT/"results/trajectories.csv"),
               "primary_linear_modes_sha256": sha256(ROOT/"results/linear_modes.csv"),
               "primary_code_sha256": sha256(ROOT/"code/run_sugra.py"),
               "method": "Independent cosmic time, physical velocities, direct complex K/W",
               "error_normalization": "max_abs_difference/max(1,max_abs_primary); dimensionless units",
               "checks_count": len(checks), "checks_passed": sum(c["passed"] for c in checks),
               "all_passed": all(c["passed"] for c in checks), "checks": checks}
    (ROOT/"results/independent_checks.json").write_text(json.dumps(summary, indent=2)+"\n")
    print(json.dumps({k: v for k, v in summary.items() if k!="checks"}, indent=2))
    if not summary["all_passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
