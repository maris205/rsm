#!/usr/bin/env python3
"""Independent proper-time/canonical-momentum audit of the threshold toy model.

The production solver uses N as its independent variable.  This implementation
instead integrates N, chi-chi_i, canonical p=A*dchi/dtau, and dissipated scalar
energy against tau=H_ref*t, then inverts the monotone N coordinate.  It does not
import production code.  Input values and comparison trajectories are archived
production artifacts, not fitted parameters.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq


ROOT = Path(__file__).resolve().parents[1]


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


class Model:
    def __init__(self, pars, kinetic, ell):
        self.pars = pars
        self.kinetic = kinetic
        self.ell = float(ell)
        self.ci = float(pars["chi_i"])
        self.f2 = float(pars["F2_eV2"])
        self.h2 = float(pars["H_ref_eV"]) ** 2
        self.mi2 = math.sqrt(32 * math.pi**2 * self.ell * self.f2 * self.h2)

    def terms(self, u):
        chi = self.ci + u
        g = 0.5 * (1 + (self.ci / chi) ** 2)
        gp = -self.ci**2 / chi**3
        gpp = 3 * self.ci**2 / chi**4
        logg = math.log(g)
        tree = self.pars["W_i"] * math.exp(-2 * u)
        # Independent differentiation of m^4[log(m^2/m_i^2)-3/2].
        cw = self.ell * (g * g * (logg - 1.5) + 1.5)
        cwp = 2 * self.ell * g * gp * (logg - 1)
        if self.kinetic == "chi":
            k, kp = 1.0, 0.0
        elif self.kinetic == "R":
            k = math.exp(2 * u)
            kp = 2 * k
        else:
            raise ValueError(f"Unknown kinetic completion {self.kinetic}")
        factor = self.mi2 / (96 * math.pi**2 * self.f2)
        da = factor * gp**2 / g
        dap = factor * (2 * gp * gpp / g - gp**3 / g**2)
        return {
            "chi": chi, "g": g, "gp": gp, "A": k + da,
            "Ap": kp + dap, "V": tree + cw, "Vp": -2 * tree + cwp,
            "cw": cw, "cwp": cwp, "da": da, "dap": dap,
        }

    def background(self, n, u, p):
        out = self.terms(u)
        rad = self.pars["omega_r"] * math.exp(-4 * n)
        mat = self.pars["omega_m"] * math.exp(-3 * n)
        numerator = rad + mat + self.pars["omega_lambda"] + self.pars["epsilon"] * out["V"] / 3
        rho = p*p / (2 * out["A"]) + out["V"]
        e2 = rad + mat + self.pars["omega_lambda"] + self.pars["epsilon"] * rho / 3
        if e2 <= 0:
            raise FloatingPointError(f"Nonpositive H^2={e2} in independent solve")
        out.update(E=math.sqrt(e2), rho=rho, numerator=numerator)
        return out

    def rhs(self, tau, state):
        n, u, p, lost = state
        terms = self.background(n, u, p)
        velocity = p / terms["A"]
        return [
            terms["E"],
            velocity,
            -3 * terms["E"] * p + 0.5 * terms["Ap"] * velocity**2 - terms["Vp"],
            3 * terms["E"] * p*p / terms["A"],
        ]

    def solve(self):
        n0, q0 = self.pars["N_i"], self.pars["q_i"]
        initial = self.terms(0)
        fluid = (self.pars["omega_r"] * math.exp(-4*n0)
                 + self.pars["omega_m"] * math.exp(-3*n0)
                 + self.pars["omega_lambda"])
        e02 = (fluid + self.pars["epsilon"] * initial["V"] / 3) / (
            1 - self.pars["epsilon"] * initial["A"] * q0*q0 / 6)
        p0 = initial["A"] * q0 * math.sqrt(e02)
        y0 = [n0, 0.0, p0, 0.0]

        def today(tau, y):
            return y[0]

        def chi_domain(tau, y):
            return 20 - abs(y[1])

        def mass_domain(tau, y):
            return 0.25 - abs(math.log(self.terms(y[1])["g"]))

        def kinetic_domain(tau, y):
            val = self.background(*y[:3])
            q = y[2] / (val["A"] * val["E"])
            return 1 - self.pars["epsilon"] * val["A"] * q*q / 6 - 1e-4

        def energy_domain(tau, y):
            return self.background(*y[:3])["numerator"] - 1e-10

        events = [today, chi_domain, mass_domain, kinetic_domain, energy_domain]
        for event in events:
            event.terminal = True
        today.direction = 1
        for event in events[1:]:
            event.direction = -1
        tau0 = self.pars["tau_i"]
        result = solve_ivp(
            self.rhs, (tau0, tau0 + 10.0), y0, method="DOP853",
            rtol=2e-11, atol=2e-13, max_step=0.002,
            dense_output=True, events=events,
        )
        if not result.success:
            raise RuntimeError(result.message)
        result.initial_rho = p0*p0 / (2 * initial["A"]) + initial["V"]
        return result


def load_input(path):
    """Only this adapter depends on the production archive's JSON schema."""
    raw = json.loads(Path(path).read_text())
    # The production author and independent auditor agree on physical input
    # names; neither imports the other's numerical functions.
    pars = dict(raw)
    pars["F2_eV2"] = raw["F_squared_eV2"]
    pars["alpha_ref"] = raw["alpha_reference"]
    return pars


def read_trajectories(path):
    groups = {}
    with Path(path).open(newline="") as handle:
        for row in csv.DictReader(handle):
            groups.setdefault(row["case"], []).append(row)
    return groups


def finite_difference_checks(model):
    errors = []
    kinetic_errors = []
    for u in (-1.0, 0.0, 1.0):
        h = 0.01
        vals = [model.terms(u + i*h) for i in (-2, -1, 1, 2)]
        fd = (vals[0]["cw"] - 8*vals[1]["cw"]
              + 8*vals[2]["cw"] - vals[3]["cw"]) / (12*h)
        exact = model.terms(u)["cwp"]
        if model.ell:
            errors.append(abs(fd-exact) / max(abs(exact), 1e-300))
            fd_da = (vals[0]["da"] - 8*vals[1]["da"]
                     + 8*vals[2]["da"] - vals[3]["da"]) / (12*h)
            exact_da = model.terms(u)["dap"]
            kinetic_errors.append(abs(fd_da-exact_da)/max(abs(exact_da),1e-300))
    return {"CW_gradient_fd_relative": max(errors, default=0.0),
            "delta_A_gradient_fd_relative": max(kinetic_errors, default=0.0)}


def number(row, name):
    return float(row[name])


def audit_case(name, rows, pars):
    case = next(item for item in pars["case_definitions"] if item["case"] == name)
    kinetic = case["kinetic"]
    ell = case["ell"]
    model = Model(pars, kinetic, ell)
    sol = model.solve()
    target_n = np.asarray([number(row, "N") for row in rows])
    initial_tau, last_tau = sol.t[0], sol.t[-1]
    achieved_n = float(sol.y[0, -1])
    if target_n[-1] > achieved_n + 1e-9:
        raise RuntimeError(f"Independent domain ended at N={achieved_n}; production target N={target_n[-1]}")
    times = []
    left = initial_tau
    for n in target_n:
        if abs(n-pars["N_i"]) < 1e-12:
            t = initial_tau
        elif abs(n-achieved_n) < 1e-10:
            t = last_tau
        else:
            t = brentq(lambda t: float(sol.sol(t)[0])-n, left, last_tau,
                       xtol=3e-15, rtol=1e-14)
        times.append(t)
        left = t
    times = np.asarray(times)
    states = sol.sol(times)
    terms = [model.background(*states[:3, i]) for i in range(len(times))]
    chi = np.asarray([item["chi"] for item in terms])
    energy = np.asarray([item["rho"] for item in terms])
    e = np.asarray([item["E"] for item in terms])
    g = np.asarray([item["g"] for item in terms])
    # Compute alpha through inverse alpha itself rather than reuse the
    # production delta/(1-delta) expression.
    log_mass_ratio = np.log(g / g[-1])
    inverse_alpha = 1/pars["alpha_ref"] - log_mass_ratio/(12*math.pi)
    delta_alpha = 1/(pars["alpha_ref"] * inverse_alpha) - 1
    prod = {key: np.asarray([number(row, key) for row in rows])
            for key in ("chi", "E", "tau", "delta_alpha", "rho", "q")}
    q = states[2] / (np.asarray([item["A"] for item in terms]) * e)
    displacement_scale = max(np.max(abs(prod["chi"]-pars["chi_i"])), 1.0)
    alpha_scale = max(np.max(abs(prod["delta_alpha"])), 1e-12)
    errors = {
        "chi_displacement_normalized": float(np.max(abs(chi-prod["chi"]))/displacement_scale),
        "H_relative": float(np.max(abs(e-prod["E"])/abs(prod["E"]))),
        "age_relative": float(np.max(abs(times-prod["tau"])/abs(prod["tau"]))),
        "alpha_normalized": float(np.max(abs(delta_alpha-prod["delta_alpha"]))/alpha_scale),
        "scalar_energy_initial_normalized": float(np.max(abs(energy-prod["rho"]))/abs(sol.initial_rho)),
        "clock_rate_normalized": float(np.max(abs(q-prod["q"]))/max(np.max(abs(prod["q"])),1.0)),
        "scalar_continuity_initial_energy_normalized": float(
            np.max(abs(energy + states[3] - sol.initial_rho)) / abs(sol.initial_rho)),
    }
    errors.update(finite_difference_checks(model))
    passed = all(value < 1e-6 for value in errors.values())
    event_names = ("today", "chi_domain", "mass_domain", "kinetic_domain", "energy_domain")
    terminated = [event_names[i] for i, ts in enumerate(sol.t_events) if len(ts)]
    record = {
        "case": name, "passed": passed, "errors": errors,
        "end_event": terminated, "independent_end_N": achieved_n,
        "solver_function_evaluations": sol.nfev,
        "independent_end_age_tau": float(last_tau),
        "initial_scalar_energy": float(sol.initial_rho),
        "alpha_normalization": alpha_scale,
        "chi_displacement_normalization": displacement_scale,
        "minimum_q": float(np.min(q)),
        "clock_reversal": bool(np.min(q) < 0),
        "delta_alpha_at_initial": float(delta_alpha[0]),
    }
    audit_rows = [{"case": name, "N": float(n), "tau": float(t),
                   "chi": float(ch), "E": float(ev), "delta_alpha": float(da),
                   "rho_scalar": float(rho), "cumulative_dissipation": float(q)}
                  for n,t,ch,ev,da,rho,q in zip(target_n,times,chi,e,delta_alpha,energy,states[3])]
    return record, audit_rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--inputs", type=Path, default=ROOT/"results/inputs.json")
    parser.add_argument("--trajectories", type=Path, default=ROOT/"results/trajectories.csv")
    args = parser.parse_args()
    pars = load_input(args.inputs)
    trajectories = read_trajectories(args.trajectories)
    results, output_rows = [], []
    for name, rows in trajectories.items():
        try:
            result, case_rows = audit_case(name, rows, pars)
            results.append(result)
            output_rows.extend(case_rows)
            print(f"{name}: passed={result['passed']} errors={result['errors']}", flush=True)
        except Exception as exc:
            result = {"case": name, "passed": False,
                      "exception_type": type(exc).__name__, "exception": str(exc)}
            results.append(result)
            print(json.dumps(result), flush=True)
    report = {
        "method": "proper time, canonical momentum; independent RHS, no production imports",
        "parameters": pars, "cases": results,
        "passed": len(results) == 8 and all(case["passed"] for case in results),
        "threshold": 1e-6,
        "input_hashes": {str(args.inputs):sha256(args.inputs),
                         str(args.trajectories):sha256(args.trajectories),
                         "protocol.md":sha256(ROOT/"protocol.md"),
                         "independent_check.py":sha256(__file__)},
    }
    (ROOT/"results/independent_checks.json").write_text(json.dumps(report, indent=2)+"\n")
    if output_rows:
        with (ROOT/"results/independent_trajectories.csv").open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(output_rows[0]))
            writer.writeheader()
            writer.writerows(output_rows)
    summary = ["# 独立数值复核", "", "本检查不导入主求解器；以宇宙时间和正则动量重新求解同一冻结模型，然后反解 N 将两个求解器的读数对齐。", "",
               "- 状态：" + ("8 个案例全部通过。" if report["passed"] else "存在失败；详见 JSON，不据此宣布通过。"),
               "- 阈值：轨迹与能量误差使用明确非零尺度归一，均要求小于 10⁻⁶。",
               "- 检查：χ 增量、H、年龄、低能 α 阈值、标量能量、时钟速度，以及 ∆ρ + ∫3H Aχ̇²dt=0；CW 斜率和动能修正斜率另作五点差分复核。",
               "- 这些检查验证两种数值实现一致，不验证新带电粒子的现实可行性，也不代表谱线观测支持模型。", "",
               "| 案例 | 通过 | 最大误差 |", "|---|---|---:|"]
    for result in results:
        err = max(result.get("errors", {}).values(), default=float("nan"))
        summary.append(f"| {result['case']} | {result['passed']} | {err:.3e} |")
    (ROOT/"reports/independent_audit.md").write_text("\n".join(summary)+"\n")
    raise SystemExit(0 if report["passed"] else 1)


if __name__ == "__main__":
    main()
