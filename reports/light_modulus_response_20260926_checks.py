#!/usr/bin/env python3
"""Independent retained-chiral EFT algebra; no imports from the main experiment."""
from __future__ import annotations

import hashlib
import json
import platform
from pathlib import Path

import mpmath as mp
import sympy as sp

mp.mp.dps = 100
checks = []


def check(name, passed, detail=None):
    row = {"name": name, "passed": bool(passed)}
    if detail is not None:
        row["detail"] = detail
    checks.append(row)


def close(a, b, tol=mp.mpf("1e-75")):
    return abs(a - b) <= tol * max(mp.mpf(1), abs(a), abs(b))


for n in (1, 2, 3, 5):
    for q in (2, 3, 5):
        Q = sp.zeros(n + 1, n)
        Q[0, 0] = q
        for a in range(1, n):
            Q[a, a - 1], Q[a, a] = 1, -q
        Q[n, n - 1] = 1
        w = sp.Matrix([1] + [-q**a for a in range(1, n + 1)])
        W = (w.T * w)[0]
        P = Q * (Q.T * Q).inv() * Q.T
        tag = f"q{q}_n{n}"
        check(f"{tag}: null direction", Q.T * w == sp.zeros(n, 1))
        check(f"{tag}: exact projector", P == sp.eye(n + 1) - w * w.T / W)
        check(f"{tag}: projector rank", P.rank() == n)
        tau = sp.cancel(q * P[0, n] / P[0, 0])
        expected = sp.Rational(q * q - 1, q) / (sp.Integer(q)**n - sp.Integer(q)**(-n))
        check(f"{tag}: cross matching", tau == expected)
        check(f"{tag}: self matching", sp.cancel(q*q * P[n, n] / P[0, 0]) == 1)
        beta_product_hat = sp.cancel(4*q*w[0]*w[n] / (W * P[0, 0]))
        beta_x_squared_hat = sp.cancel(4*q*q*w[n]**2 / (W * P[0, 0]))
        check(f"{tag}: current product", beta_product_hat == -4*tau)
        check(f"{tag}: hidden current unsuppressed", beta_x_squared_hat == 4*(q*q-1)/(1-sp.Integer(q)**(-2*n)))
        j = sp.Matrix(sp.symbols(f"j0:{n+1}"))
        ell_light = sp.Symbol("s") * w
        ell = ell_light - P*j
        # v=1: v² L²+2v² J·L after massive-vector elimination.
        actual = (ell.T*ell)[0] + 2*(j.T*ell)[0]
        expected_k = (ell_light.T*ell_light)[0] + 2*(j.T*ell_light)[0] - (j.T*P*j)[0]
        check(f"{tag}: retained-coordinate completion", sp.expand(actual-expected_k) == 0)

# Local two-chiral inverse metric with hidden scalar kept, not frozen.
sig, bx, lam, S, R = sp.symbols("sig bx Lambda S R", real=True, positive=True)
metric_hidden_schur = 1 + bx*sig - 4*R/lam**2 - bx**2*R/2
potential = S/metric_hidden_schur
hidden_local_mass = sp.diff(potential, R).subs({sig: 0, R: 0})
check("hidden Schur mass includes light mixing", sp.simplify(hidden_local_mass-S*(4/lam**2+bx**2/2)) == 0)

# The nonstationary clock Hessian must not be presented as vacuum stability.
bi, tau, y = sp.symbols("bi tau y", real=True)
V = S/(1+bx*sig-2*tau*y*y/lam**2)
gradient_sig = sp.diff(V, sig).subs({sig: 0, y: 0})
ordinary_yy = sp.diff(V, y, 2).subs({sig: 0, y: 0})
Gamma_sig_yy = bi/2
covariant_yy = ordinary_yy-Gamma_sig_yy*gradient_sig
check("nonzero light gradient", sp.simplify(gradient_sig+S*bx) == 0)
check("ordinary clock Hessian", sp.simplify(ordinary_yy-4*tau*S/lam**2) == 0)
check("covariant clock Hessian", sp.simplify(covariant_yy.subs(bi, -4*tau/(bx*lam**2))-2*tau*S/lam**2) == 0)
sig_new = sp.symbols("sig_new", real=True)
transformed = V.subs(sig, sig_new+2*tau*y*y/(bx*lam**2))
check("coordinate ordinary Hessian can vanish", sp.diff(sp.simplify(transformed), y) == 0)

# Independent exact normalization for the frozen benchmark, evaluated at high precision.
q, n = 3, 127
L = mp.mpf("2e12")
g = mp.mpf(1)
Pgrav = mp.mpf("5.929225e54")
mg = mp.mpf("1e-15")
source = 3*Pgrav*mg**2
W = (mp.mpf(q)**(2*n+2)-1)/(q*q-1)
u = [1/mp.sqrt(W)] + [-mp.mpf(q)**a/mp.sqrt(W) for a in range(1, n+1)]
p00 = 1-1/W
tau_num = (mp.mpf(q)-1/mp.mpf(q))/(mp.mpf(q)**n-mp.mpf(q)**(-n))
v = L*mp.sqrt(p00)/(2*g*q)
ci = 1/(v*L*mp.sqrt(p00))
cx = q*ci
beta_i = 2*v*ci*u[0]
beta_x = 2*v*cx*u[-1]
b = -2*v*v*cx*u[-1]
check("benchmark quartic normalization", close(v*v*ci*ci*p00, 1/L**2))
check("benchmark product identity", close(beta_i*beta_x*L**2/tau_num, -4))
check("benchmark betaX squared", close(beta_x**2*L**2, 4*(q*q-1)/(1-mp.mpf(q)**(-2*n))))
check("benchmark light gradient not tiny-overlap suppressed", abs(beta_x*L) > 1)
check("benchmark origin nonstationary", abs(source*beta_x) > 0)

# Independent homogeneous Minkowski time from energy conservation, not an ODE solve.
def B(r):
    return sum(ua**2/mp.sqrt(1+ua**2*r*r) for ua in u)

target_metric_change = mp.mpf("0.1")
r_target = target_metric_change/b
# r=t² removes the integrable square-root endpoint singularity.
def transformed_time_integrand(t):
    if t == 0:
        return mp.sqrt(2/b)
    r = t*t
    return mp.sqrt(2*B(r)*(1+b*r)/b)

s_target = mp.quad(transformed_time_integrand, [0, mp.sqrt(r_target)])
hbar_eVs = mp.mpf("6.582119569e-16")
time_unit_s = hbar_eVs*v/mp.sqrt(source)
rprime = mp.sqrt(2*(1-1/(1+b*r_target))/B(r_target))
check("D-flat metric canonical at origin", close(B(0), 1))
check("target energy conservation", close(B(r_target)*rprime**2/2+1/(1+b*r_target), 1))
check("target hidden metric positive", 1+b*r_target > 0)
check("energy integral positive finite", mp.isfinite(s_target) and s_target > 0)

def num(x):
    return mp.nstr(x, 45)

benchmark = {
    "q": q, "n": n, "Lambda_eV": num(L), "g": num(g),
    "m_gravitino_eV": num(mg), "S_eV4": num(source), "v_eV": num(v),
    "u_left": num(u[0]), "u_right": num(u[-1]), "tau": num(tau_num),
    "c_I_eV_minus2": num(ci), "c_X_eV_minus2": num(cx),
    "beta_I_eV_minus1": num(beta_i), "beta_X_eV_minus1": num(beta_x),
    "beta_X_times_Lambda": num(beta_x*L),
    "dV_dsigma_eV3": num(-source*beta_x),
    "ordinary_clock_Hessian_eV2": num(4*tau_num*source/L**2),
    "covariant_clock_Hessian_eV2": num(2*tau_num*source/L**2),
    "local_hidden_Hessian_eV2": num(source*(4/L**2+beta_x**2/2)),
    "sqrt_local_hidden_Hessian_eV_NOT_vacuum_mass": num(mp.sqrt(source*(4/L**2+beta_x**2/2))),
    "local_saxion_Hessian_eV2": num(2*source*beta_x**2),
    "source_over_v4": num(source/v**4),
}
transient = {
    "scope": "Fixed-source homogeneous Minkowski D-flat truncation, initially at rest; not a cosmological solution or SUGRA lifetime.",
    "b": num(b), "relative_hidden_metric_change": num(target_metric_change),
    "r_target": num(r_target), "dimensionless_time_energy_quadrature": num(s_target),
    "time_unit_seconds": num(time_unit_s), "target_time_seconds": num(time_unit_s*s_target),
    "B_target": num(B(r_target)), "dimensionless_rprime_target": num(rprime),
}
path = Path(__file__)
out = {
    "scope": "Own analytic derivation and independent implementation, same model family; no external empirical data or external review.",
    "versions": {"python": platform.python_version(), "sympy": sp.__version__, "mpmath": mp.__version__},
    "precision_decimal_digits": mp.mp.dps, "numeric_comparison_tolerance": "1e-75 * max(1, abs(a), abs(b))",
    "script_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
    "passed": sum(c["passed"] for c in checks), "total": len(checks),
    "checks": checks, "benchmark": benchmark, "conditional_transient": transient,
}
path.with_suffix(".json").write_text(json.dumps(out, indent=2, ensure_ascii=False)+"\n")
print(json.dumps({"passed": out["passed"], "total": out["total"], "benchmark": benchmark, "conditional_transient": transient}, indent=2))
if out["passed"] != out["total"]:
    raise SystemExit(1)
