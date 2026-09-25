#!/usr/bin/env python3
"""Exact, local symmetry/soft-term audit; no import of the main experiment.

All barred holomorphic variables are independent symbols. Unit phases are
represented by an algebraic nonzero u and the barred variable transforms by
u**(-1). This tests polynomial identities, not numerical phase sampling.
"""

import hashlib
import json
from pathlib import Path

import sympy as s

checks = []


def equal(name, lhs, rhs=0):
    residual = s.factor(s.simplify(lhs - rhs))
    checks.append({"name": name, "passed": residual == 0,
                   "residual": str(residual)})


def nonzero(name, expression):
    residual = s.factor(s.simplify(expression))
    checks.append({"name": name, "passed": residual != 0,
                   "nonzero_witness": str(residual)})


P = s.symbols("P", positive=True)
T, Tb, Z, Zb = s.symbols("T Tb Z Zb")
qp, qpb, qm, qmb = s.symbols("Qp Qpb Qm Qmb")
e, eb = s.symbols("e eb")
alpha, shift = s.symbols("alpha shift", real=True)
k = -(Z - Zb) ** 2 / 2
matter = qp * qpb + qm * qmb
Y = T + Tb - (k + matter) / (3 * P)
S = T + Z ** 2 / (6 * P)
Sb = Tb + Zb ** 2 / (6 * P)
equal("clock_holomorphic_coordinate_to_standard_coset", Y,
      S + Sb - (Z * Zb + matter) / (3 * P))
heis = {qp: qp + e, qpb: qpb + eb,
        T: T + eb * qp / (3 * P) + e * eb / (6 * P) + s.I * alpha,
        Tb: Tb + e * qpb / (3 * P) + e * eb / (6 * P) - s.I * alpha}
equal("finite_heisenberg_Y_invariance", Y.subs(heis, simultaneous=True), Y)
equal("imaginary_T_shift_Y_invariance",
      Y.subs({T: T + s.I * alpha, Tb: Tb - s.I * alpha}, simultaneous=True), Y)
equal("real_clock_shift_k_invariance",
      k.subs({Z: Z + shift, Zb: Zb + shift}, simultaneous=True), k)

a, b, eta = s.symbols("a b eta", positive=True)
delta, v = s.symbols("delta v", real=True)
quartic = a * delta ** 4 + b * v ** 4
nonzero("b_breaks_axionic_T_translation", quartic.subs(v, v + 2 * alpha) - quartic)
equal("quartic_preserves_CP", quartic.subs(v, -v), quartic)
rotation_generator = -v * s.diff(quartic, delta) + delta * s.diff(quartic, v)
equal("centered_continuous_rotation_breaking_polynomial", rotation_generator,
      4 * delta * v * (b * v ** 2 - a * delta ** 2))
nonzero("positive_quartics_not_centered_U1", rotation_generator)
e_real = s.symbols("e_real", real=True)
equal("Heisenberg_breaking_witness_at_center",
      (a * (T + Tb - 1) ** 4).subs(heis, simultaneous=True)
      .subs({T: s.Rational(1, 2), Tb: s.Rational(1, 2), qp: 0, qpb: 0,
             e: e_real, eb: e_real, alpha: 0}), a * e_real ** 8 / (81 * P ** 4))
nonzero("quartic_Heisenberg_breaking_witness_nonzero", a * e_real ** 8 / (81 * P ** 4))
lam = s.symbols("lam", positive=True)
scale = {T: lam * T, Tb: lam * Tb, Z: s.sqrt(lam) * Z,
         Zb: s.sqrt(lam) * Zb, qp: s.sqrt(lam) * qp,
         qpb: s.sqrt(lam) * qpb, qm: s.sqrt(lam) * qm, qmb: s.sqrt(lam) * qmb}
equal("ideal_homogeneous_scaling", Y.subs(scale, simultaneous=True), lam * Y)
t = s.symbols("t", real=True)
nonzero("centered_quartic_breaks_homogeneous_scaling",
        a * (lam * t - 1) ** 4 + b * lam ** 4 * v ** 4
        - lam * (a * (t - 1) ** 4 + b * v ** 4))

M, Wc = s.symbols("M Wc")
equal("mass_superpotential_breaks_Q_translation",
      (Wc + M * qp * qm).subs(qp, qp + e) - (Wc + M * qp * qm), M * e * qm)
nonzero("mass_breaking_nonzero", M * e * qm)
g, A = s.symbols("g A")
equal("charged_translation_covariant_derivative_mismatch", -s.I * g * A * (qp + e)
      - (-s.I * g * A * qp), -s.I * g * A * e)
nonzero("gauge_breaking_nonzero", -s.I * g * A * e)
r = s.symbols("r", positive=True)
nonzero("fixed_Wc_and_mass_not_required_scale_weight",
        Wc + r ** 2 * M * qp * qm - r ** 3 * (Wc + M * qp * qm))
kap = s.symbols("kap", positive=True)
w = s.exp(-kap * Z)
equal("clock_exponential_shift_factor", w.subs(Z, Z + shift), s.exp(-kap * shift) * w)
nonzero("constant_and_clock_not_same_shift_character",
        Wc + s.exp(-kap * shift) * w - s.exp(-kap * shift) * (Wc + w))
xi, minf = s.symbols("xi minf", positive=True)
equal("clock_mass_squared_shift_derivative",
      s.diff(minf ** 2 * (1 + xi / Z ** 2), Z), -2 * minf ** 2 * xi / Z ** 3)

tau, taub = s.symbols("tau taub")
u, up, um = s.symbols("u up um", nonzero=True)
contact = tau * taub * matter
phase = {tau: u * tau, taub: taub / u, qp: up * qp, qpb: qpb / up,
         qm: um * qm, qmb: qmb / um}
equal("dangerous_contact_invariant_under_all_linear_phases_including_R",
      contact.subs(phase, simultaneous=True), contact)
gauge = {qp: u * qp, qpb: qpb / u, qm: qm / u, qmb: u * qmb}
equal("contact_U1_gauge_invariant", contact.subs(gauge, simultaneous=True), contact)
equal("W_mass_U1_gauge_invariant", (qp * qm).subs(gauge, simultaneous=True), qp * qm)
equal("contact_centered_parity_even",
      contact.subs({tau: -tau, taub: -taub}, simultaneous=True), contact)
equal("contact_CP_even", contact.subs({tau: taub, taub: tau, qp: qpb, qpb: qp,
      qm: qmb, qmb: qm}, simultaneous=True), contact)
equal("quadratic_hidden_deformation_phase_invariant",
      (tau * taub).subs(phase, simultaneous=True), tau * taub)
cubic = tau ** 2 * taub + tau * taub ** 2
equal("cubic_hidden_deformation_CP_even",
      cubic.subs({tau: taub, taub: tau}, simultaneous=True), cubic)
equal("cubic_odd_under_hypothetical_centered_parity",
      cubic.subs({tau: -tau, taub: -taub}, simultaneous=True), -cubic)
nonzero("centered_parity_not_actual_base_K_invariance",
        (1 + tau + taub).subs({tau: -tau, taub: -taub}, simultaneous=True)
        - (1 + tau + taub))
c = s.symbols("c", real=True)
ell = c * tau * taub
equal("contact_mixed_log_metric_curvature", s.diff(ell, tau, taub), c)
h = s.Function("h")
hb = s.Function("hb")
equal("holomorphic_rescaling_cannot_change_mixed_curvature",
      s.diff(ell + h(tau) + hb(taub), tau, taub), c)
X = s.symbols("X", positive=True)
equal("soft_mass_contact_shift", -X * s.diff(ell, tau, taub), -c * X)

# Direct one-complex-modulus SUGRA potential with constant W, not a mass ansatz.
c2, c3 = s.symbols("c2 c3", real=True)
f = 1 + tau + taub + a * (tau + taub) ** 4 + b * (tau - taub) ** 4
f += c2 * tau * taub + c3 * (tau ** 2 * taub + tau * taub ** 2)
K = -3 * P * s.log(f)
Ktt = s.diff(K, tau, taub)
Wsq = s.symbols("Wsq", positive=True)
potential = f ** (-3) * Wsq / P * (s.diff(K, tau) * s.diff(K, taub) / (P * Ktt) - 3)
q = s.diff(f, tau, taub)
D = s.diff(f, tau) * s.diff(f, taub) - f * q
equal("direct_SUGRA_potential_identity", potential, 3 * q * Wsq / (P * f ** 2 * D))
center = {tau: 0, taub: 0}
equal("hidden_c2_center_potential", potential.subs(center), 3 * c2 * Wsq / (P * (1 - c2)))
equal("hidden_c2_center_metric", Ktt.subs(center), 3 * P * (1 - c2))
equal("hidden_c2_compensator_over_Wbar_P",
      (1 - s.diff(K, tau) * s.diff(K, taub) / (3 * P * Ktt)).subs(center),
      -c2 / (1 - c2))
equal("hidden_c2_common_soft_mass", 2 * potential.subs(center) / (3 * P),
      2 * c2 * Wsq / (P ** 2 * (1 - c2)))
equal("hidden_c3_center_potential_zero", potential.subs(c2, 0).subs(center), 0)
equal("hidden_c3_real_gradient", (s.diff(potential, tau) + s.diff(potential, taub))
      .subs(c2, 0).subs(center) / 2, 6 * c3 * Wsq / P)
Vbaseline = potential.subs({c2: 0, c3: 0})
equal("quartic_real_hessian", (s.diff(Vbaseline, tau, 2) + 2 * s.diff(Vbaseline, tau, taub)
      + s.diff(Vbaseline, taub, 2)).subs(center) / 4, 72 * a * Wsq / P)
equal("quartic_imaginary_hessian", (-s.diff(Vbaseline, tau, 2) + 2 * s.diff(Vbaseline, tau, taub)
      - s.diff(Vbaseline, taub, 2)).subs(center) / 4, 72 * b * Wsq / P)
equal("cubic_leading_displacement", -(6 * c3 * Wsq / P) / (72 * a * Wsq / P), -c3 / (12 * a))

# Nonzero K alone would not rule out a holomorphic+antiholomorphic Kähler
# transformation. A mixed derivative is a direct metric counterexample.
Kbase = K.subs({c2: 0, c3: 0})
axion_metric_change = (s.diff(Kbase, tau, taub).subs(
    {tau: s.I / 8, taub: -s.I / 8, a: 1, b: 1, P: 1})
    - s.diff(Kbase, tau, taub).subs({tau: 0, taub: 0, a: 1, b: 1, P: 1}))
nonzero("quartic_axion_breaking_not_a_Kahler_transformation", axion_metric_change)
qY = 1 - qp * qpb / (3 * P)
e0 = s.Rational(1, 2)
dshift = (e0 * (qp + qpb) + e0 ** 2) / (3 * P)
vshift = -s.I * e0 * (qp - qpb) / (3 * P)
dK_heis = -3 * P * s.log((qY + a * dshift ** 4 + b * vshift ** 4) / qY)
nonzero("quartic_Q_translation_breaking_not_a_Kahler_transformation",
        s.diff(dK_heis, qp, qpb).subs({qp: 0, qpb: 0, a: 1, b: 1, P: 1}))

script = Path(__file__).resolve()
result = {"scope": "Exact finite transformations, permitted local operators, constant-W SUGRA derivatives; not a classification of every possible UV symmetry",
          "sympy_version": s.__version__, "script_sha256": hashlib.sha256(script.read_bytes()).hexdigest(),
          "checks": checks, "total": len(checks), "passed": sum(x["passed"] for x in checks),
          "physical_inference": "The tested existing symmetries do not force c=0. Permitted does not imply a computed nonzero loop coefficient."}
target = script.with_suffix(".json")
target.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n")
print(json.dumps({"passed": result["passed"], "total": result["total"], "output": str(target)}))
raise SystemExit(0 if all(x["passed"] for x in checks) else 1)
