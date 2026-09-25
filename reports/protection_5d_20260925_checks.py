"""Independent coefficient and convention checks for the cited flat 5D sector.

This is not a rerun of the paper's full supergravity diagram calculation.
Run from repository root: python reports/protection_5d_20260925_checks.py
"""

from pathlib import Path
import hashlib
import json

import mpmath as mp
import sympy as sp


mp.mp.dps = 80
rows = []


def record(name, passed, **detail):
    rows.append({"name": name, "passed": bool(passed), **detail})


def numerical(name, actual, expected, tol=mp.mpf("1e-65")):
    error = abs(actual - expected) / max(abs(expected), mp.mpf(1))
    record(name, error < tol, relative_to_max_one=str(error))


def identity(name, expr):
    record(name, sp.simplify(expr) == 0)


# y = 6 M5 L x in RSS Eq. (2.7). Differentiate analytically at zero
# boundary backgrounds before numerical quadrature, avoiding finite differences.
cuts = [0, mp.mpf("0.1"), 1, 5, 20, mp.inf]
zeta3 = mp.zeta(3)
i0 = mp.quad(lambda y: y * mp.log1p(-mp.exp(-y)), cuts)
i1 = mp.quad(lambda y: -y**2 / (3 * mp.expm1(y)), cuts)
i2 = mp.quad(
    lambda y: -y**3 * mp.exp(-y) / (9 * (-mp.expm1(-y))**2), cuts
)
numerical("RSS_integral_constant", i0, -zeta3)
numerical("RSS_integral_one_boundary_derivative", i1, -2 * zeta3 / 3)
numerical("RSS_integral_cross_boundary_derivative", i2, -2 * zeta3 / 3)
numerical("RSS_A", -i0 / (4 * mp.pi**2), zeta3 / (4 * mp.pi**2))
numerical("RSS_B", -i1 / (4 * mp.pi**2), zeta3 / (6 * mp.pi**2))
numerical("RSS_C", -i2 / (4 * mp.pi**2), zeta3 / (6 * mp.pi**2))

P, L0, t, z, Q2, X2, mg, R, eps = sp.symbols(
    "P L0 t z Q2 X2 mg R eps", positive=True
)
L = L0 * t
m5cube = 2 * P / L0
lam = z / (12 * sp.pi**2 * P * L0**2)
frame0 = -sp.Rational(3, 2) * m5cube * L + Q2 + X2
loop = z / (4 * sp.pi**2 * L**2)
loop += z * (Q2 + X2) / (6 * sp.pi**2 * m5cube * L**3)
loop += z * Q2 * X2 / (6 * sp.pi**2 * m5cube**2 * L**4)
omega = -(frame0 + loop) / (3 * P)
expected = t - lam / t**2
expected -= (1 + lam / t**3) * (Q2 + X2) / (3 * P)
expected -= lam * Q2 * X2 / (6 * P**2 * t**4)
identity("normalized_frame_complete_displayed_terms", omega - expected)
identity("tree_Einstein_mass", m5cube * L0 / 2 - P)
identity("LO_RSS_M5_cube_ratio", m5cube / (P / L0) - 2)
identity("loop_log_sign", (-(loop.subs({Q2: 0, X2: 0})) / (3 * P)) + lam / t**2)
identity("canonical_boundary_metric_factor", -3 * P * sp.diff(omega, Q2).subs(X2, 0) - (1 + lam / t**3))
identity("cross_boundary_contact", sp.diff(omega, Q2, X2) + lam / (6 * P**2 * t**4))
identity("radion_second_derivative", sp.diff(t - lam / t**2, t, 2) + 6 * lam / t**4)
identity("Casimir_LO_RSS", -18 * lam * mg**2 * P + 3 * z * mg**2 / (2 * sp.pi**2 * L0**2))
identity("RSS_radius_convention_Eq6_3", -(2 * sp.pi * eps)**2 * 6 / (2 * sp.pi * R)**4 + 3 * eps**2 / (2 * sp.pi**2 * R**4))
identity("KK_conversion", lam.subs(L0, 2 * sp.pi / R) - z * R**2 / (48 * sp.pi**4 * P))

# Scaling of the frame-function terms under masses -> scale * masses and
# length -> length / scale: every term has mass dimension two.
scale = sp.symbols("scale", positive=True)
scaling = {P: scale**2 * P, L0: L0 / scale, Q2: scale**2 * Q2, X2: scale**2 * X2}
identity("frame_function_mass_dimension_two", (frame0 + loop).subs(scaling, simultaneous=True) - scale**2 * (frame0 + loop))
identity("lambda_dimensionless", lam.subs(scaling, simultaneous=True) - lam)
identity("vacuum_mass_dimension_four", (-18 * lam * mg**2 * P).subs({**scaling, mg: scale * mg}, simultaneous=True) + 18 * scale**4 * lam * mg**2 * P)

mp_ev = mp.mpf("2.435e27")
benchmarks = []
for mkk_gev in [100, 1000, 10000, 1000000]:
    mkk = mp.mpf(mkk_gev) * 10**9
    length = 2 * mp.pi / mkk
    m5 = (2 * mp_ev**2 / length) ** (mp.mpf(1) / 3)
    coef = zeta3 / (12 * mp.pi**2 * mp_ev**2 * length**2)
    mq_ratio = mp.mpf("1e11") / mkk
    numerical(f"Planck_conversion_{mkk_gev}GeV", m5**3 * length / (2 * mp_ev**2), 1)
    numerical(f"lambda_conversion_{mkk_gev}GeV", coef, zeta3 / (48 * mp.pi**4) * (mkk / mp_ev)**2)
    benchmarks.append({
        "mKK_GeV": mkk_gev, "M5_RSS_GeV": str(m5 / 10**9),
        "lambda5": str(coef), "M_Q_over_mKK": str(mq_ratio),
        "mKK_over_M5": str(mkk / m5),
        "qualifier": "hierarchy illustration only; M_Q = mKK point is not controlled four-dimensional threshold matching",
    })

record("positive_gravity_frame_coefficient", zeta3 / (4 * mp.pi**2) > 0)
record("negative_pure_gravity_potential", -3 * zeta3 / (2 * mp.pi**2) < 0)
record("all_benchmarks_KK_below_M5", all(mp.mpf(b["mKK_over_M5"]) < 1 for b in benchmarks))

out = {
    "scope": "Coefficient-integral, convention, and dimensional checks; not original full 5D diagram reproduction or observational evidence",
    "mpmath_dps": mp.mp.dps,
    "checks": rows,
    "passed": sum(r["passed"] for r in rows),
    "total": len(rows),
    "benchmarks": benchmarks,
    "kappa5_zeta3_over_48pi4": str(zeta3 / (48 * mp.pi**4)),
    "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
}
dest = Path(__file__).with_suffix(".json")
dest.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n")
print(json.dumps({"passed": out["passed"], "total": out["total"], "output": str(dest)}))
if out["passed"] != out["total"]:
    raise SystemExit(1)
