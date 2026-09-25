#!/usr/bin/env python3
"""Symbolic identities for a conditional constant–clock bridge; no simulation."""
from pathlib import Path
import hashlib
import json
import sympy as s

R, Rs, b, I, V, P, Q = s.symbols('R Rs b I V P Q', positive=True)
mi = s.symbols('mi', real=True)
chi = s.log(R / Rs)
f = mi + b / chi**2
fp, fpp = s.diff(f, R), s.diff(f, R, 2)
Pi = P / fp
checks = []


def check(name, expression):
    residual = s.simplify(expression)
    checks.append({'name': name, 'passed': residual == 0, 'residual': str(residual)})


check('canonical_point_transform', s.diff(f, R) * s.diff(Pi, P) - 1)
check('kinetic_energy_preserved', fp**2 * Pi**2 / (2 * I) - P**2 / (2 * I))
check('inverse_log_first_derivative', fp + 2 * b / (R * chi**3))
check('inverse_log_second_derivative', fpp - 2 * b * (chi + 3) / (R**2 * chi**4))
Rdot, Pdot = P / I, -Q * fp / 2
Pi_chain = s.diff(Pi, R) * Rdot + s.diff(Pi, P) * Pdot
check('transformed_momentum_equation', Pi_chain - (-fpp * Pi**2 / I - Q / 2))
check('clock_sign_and_volume', -s.diff(P, P) * s.diff(-R / V, R) - 1 / V)
check('lattice_clock_energy_exchange', Q * fp * Rdot / 2 + P * Pdot / I)

g1, g2, x, sa, sb, xdot = s.symbols('g1 g2 x sa sb xdot', real=True)
check('rank_one_cross_epoch', (g1 * sa) * (g2 * sb) - (g1 * sb) * (g2 * sa))
D1, D2 = g1 * xdot / (1 + g1 * x), g2 * xdot / (1 + g2 * x)
check('logarithmic_drift_cross_product', g2 * (1 + g1 * x) * D1 - g1 * (1 + g2 * x) * D2)
c0, v = s.symbols('c0 v', positive=True)
response = 1 + g1 * (chi**-2 - c0**-2)
check('alpha_logarithmic_drift', s.diff(response, R) * v / response + 2 * g1 * v / (R * chi**3 * response))

B, BR, A2, Bmag2 = s.symbols('B BR A2 Bmag2', positive=True)
hem = A2 / (2 * B) + B * Bmag2 / 2
F2 = 2 * (Bmag2 - A2 / B**2)
check('EM_derivative_at_fixed_canonical_variables', s.diff(hem, B) * BR - BR * F2 / 4)
t, C, Rinf = s.symbols('t C Rinf', positive=True)
Rmatter = Rinf - C / t
check('free_R_matter_background_counterexample', s.diff(Rmatter, t, 2) + 3 * (s.Rational(2, 3) / t) * s.diff(Rmatter, t))

ct, eps = s.symbols('ct eps', real=True)
delta_threshold = ct * s.log(1 + eps) / (1 - ct * s.log(1 + eps))
check('threshold_linear_coefficient', s.diff(delta_threshold, eps).subs(eps, 0) - ct)
check('threshold_quadratic_coefficient', s.diff(delta_threshold, eps, 2).subs(eps, 0) / 2 - (ct**2 - ct / 2))
alpha0, bem, eta, massinf2, mass02 = s.symbols('alpha0 bem eta massinf2 mass02', positive=True)
mass2 = massinf2 + eta / chi**2
alpha_ir = 1 / (1 / alpha0 - bem / (4 * s.pi) * s.log(mass2 / mass02))
check('moving_threshold_drift_sign_and_factor', s.diff(alpha_ir, R) * v / alpha_ir + bem * alpha_ir * eta * v / (2 * s.pi * mass2 * R * chi**3))

root = Path(__file__).resolve().parent
result = {
    'scope': 'Exact symbolic identities on the stated nonsingular branch; no observational or dynamical convergence test.',
    'sympy_version': s.__version__,
    'domain': 'R>Rs>0, b,I,V>0; denominators nonzero; B>0. Some checks also hold on larger domains.',
    'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    'checks': checks, 'passed': all(c['passed'] for c in checks),
    'n_checks': len(checks), 'n_passed': sum(c['passed'] for c in checks),
}
(root / 'derivation_checks.json').write_text(json.dumps(result, ensure_ascii=False, indent=2) + '\n')
print(json.dumps({k: result[k] for k in ['passed', 'n_checks', 'n_passed']}))
raise SystemExit(0 if result['passed'] else 1)
