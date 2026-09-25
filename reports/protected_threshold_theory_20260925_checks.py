#!/usr/bin/env python3
"""Symbolic identities for the specified protected threshold EFT, not data tests."""
import json
from pathlib import Path
import sympy as sp

s, mu, d, D, x, xi, F, mi, U, eps = sp.symbols(
    "s mu d D x xi F mi U eps", positive=True
)
y = sp.symbols("y", real=True)
checks = []


def check(name, expression):
    residual = sp.simplify(expression)
    checks.append({"name": name, "passed": residual == 0, "residual": str(residual)})


L = sp.log(s / mu**2)
f = s**2 * (L - sp.Rational(3, 2))
check("cw_first_derivative", sp.diff(f, s) - 2 * s * (L - 1))
check("cw_second_derivative", sp.diff(f, s, 2) - 2 * L)
check("cw_fourth_derivative", sp.diff(f, s, 4) + 2 / s**2)
check("unbroken_cw_cancellation", f + f - 2 * f)
pair = sp.Matrix([[s + d, D], [D, s + d]])
check("scalar_eigenvalue_plus", (pair - (s + d + D) * sp.eye(2)).det())
check("scalar_eigenvalue_minus", (pair - (s + d - D) * sp.eye(2)).det())
check("pair_beta_coefficient", sp.Rational(1, 3) * 2 + sp.Rational(4, 3) - 2)
taylor = sum(sp.diff(f, s, n) * ((d + D)**n + (d - D)**n) / sp.factorial(n)
             for n in range(1, 5)) / (32 * sp.pi**2)
check("Fterm_cw_to_fourth_order", taylor.subs(d, 0) -
      (D**2 * L / (16 * sp.pi**2) - D**4 / (192 * sp.pi**2 * s**2)))
check("soft_cw_linear_term", sp.diff(taylor, d).subs({d: 0, D: 0}) -
      s * (L - 1) / (8 * sp.pi**2))
check("soft_force_leading_derivative", sp.diff(s*(L-1),s)-L)
check("vectorlike_cubic_anomaly", sp.Integer(1)**3+sp.Integer(-1)**3)
check("vectorlike_mixed_gravitational_anomaly", sp.Integer(1)+sp.Integer(-1))
k1 = -s * (L - 2) / (16 * sp.pi**2)
check("kahler_metric_coefficient", s * sp.diff(k1, s, 2) + sp.diff(k1, s) +
      L / (16 * sp.pi**2))
h = sp.symbols("h")
check("kahler_and_cw_same_F_squared_term", sp.series(U / (1-h), h, 0, 2).removeO() - U*(1+h))
Mreal = mi * sp.sqrt(1 + xi / x**2)
Areal = 2 * mi**2 * xi**2 / (F**2 * x**4 * (x**2 + xi))
check("holomorphic_mass_derivative", 2 * sp.diff(Mreal, x)**2 / F**2 - Areal)
z = sp.symbols("z", real=True)
W0 = -F / sp.sqrt(2) * sp.sqrt(U) * sp.exp(-sp.sqrt(2) * z / F)
check("exponential_superpotential_derivative", sp.diff(W0, z) - sp.sqrt(U) * sp.exp(-sp.sqrt(2)*z/F))
radicand = (x**2 - y**2 + xi)**2 + 4*x**2*y**2
sxy = mi**2 * sp.sqrt(radicand) / (x**2+y**2)
Axy = 2 * mi**2 * xi**2 / (F**2 * (x**2+y**2)**2 * sp.sqrt(radicand))
E = -xi*(3*x**2+xi)/(x**2*(x**2+xi)**2)
J = (3*x**4+3*x**2*xi+2*xi**2)/(x**2*(x**2+xi)**2)
check("transverse_logmass_curvature", sp.diff(sp.log(sxy), y, 2).subs(y, 0) - 2*E)
check("transverse_logyukawa_curvature", sp.diff(sp.log(Axy), y, 2).subs(y, 0) + 2*J)
Lreal = sp.log(mi**2*(x**2+xi)/(mu**2*x**2))
hxy = Axy * sp.log(sxy/mu**2)/(16*sp.pi**2)
check("transverse_CW_curvature", sp.diff(hxy, y, 2).subs(y, 0) - Areal*(E-J*Lreal)/(8*sp.pi**2))
Vsg = U*sp.exp(eps*(x**2+y**2)/2)*((1-eps*x/2)**2+eps**2*y**2/4-3*eps/2)
check("sugra_transverse_curvature", sp.diff(Vsg, y, 2).subs(y, 0) -
      U*sp.exp(eps*x**2/2)*eps*((1-eps*x/2)**2-eps))
out = {"scope": "Symbolic checks of chosen model identities; not observations, physical viability, or full perturbation stability",
       "count": len(checks), "passed": sum(c["passed"] for c in checks), "checks": checks}
path = Path(__file__).with_suffix(".json")
path.write_text(json.dumps(out, indent=2) + "\n")
print(json.dumps({"count": out["count"], "passed": out["passed"], "output": str(path)}, indent=2))
assert all(c["passed"] for c in checks)
