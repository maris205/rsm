#!/usr/bin/env python3
"""Independent exact-algebra tests of prescribed endpoint bypasses.

No main experiment or prior check implementation is imported. These are
sensitivity tests, not a computation of radiative generation or a loop floor.
Run from any directory; outputs the adjacent *_checks.json file.
"""
from pathlib import Path
import hashlib
import json
import sys

import mpmath as mp
import sympy as sp

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
checks = []


def check(name, passed, detail=None):
    row = {"name": name, "passed": bool(passed)}
    if detail is not None:
        row["detail"] = detail
    checks.append(row)


def same(a, b):
    return sp.cancel(a - b) == 0


# Symbolic endpoint inverse and source transformations, independent of n.
f, t, d, eps = sp.symbols("f t d eps", real=True)
C = f * sp.Matrix([[1, t], [t, 1]])
E = sp.Matrix([[0, 1], [1, 0]])
D = 1 + 2*d*f*t - d*d*f*f*(1-t*t)
td = t - d*f*(1-t*t)
Cd = (C.inv() + d*E).inv()
for i in range(2):
    for j in range(2):
        check(f"symbolic_mass_endpoint_{i}{j}",
              same(Cd[i, j], f/D * (1 if i == j else td)))
L = sp.Matrix([[1, eps], [0, 1]])
Ce = L.T * C * L
check("symbolic_current_cross", same(Ce[0, 1]/f, t+eps))
check("symbolic_current_self", same(Ce[1, 1]/f, 1+2*eps*t+eps**2))
check("symbolic_current_determinant", same(Ce.det(), f*f*(1-t*t)))
Cboth = L.T * Cd * L
check("symbolic_combined_cross", same(Cboth[0, 1]/Cboth[0, 0], td+eps))
check("symbolic_combined_self", same(Cboth[1, 1]/Cboth[0, 0], 1+2*eps*td+eps**2))
check("symbolic_updated_determinant", same(Cd.det(), f*f*(1-t*t)/D))
Zp, X, product, vev2, lam, kap = sp.symbols("Zp X product vev2 lam kap")
check("neutral_link_W_field_redefinition",
      sp.expand(lam*(Zp-kap*X/lam)*(product-vev2)+kap*X*product
                -lam*Zp*(product-vev2)-kap*vev2*X) == 0)
for label, delta, target in [("lower", -t/(f*(1-t*t)), 2*t),
                              ("upper", t/(2*f*(1-t*t)), t/2),
                              ("zero", t/(f*(1-t*t)), 0)]:
    check("mass_factor_two_"+label, same(td.subs(d, delta), target))
for label, e, target in [("lower", -t/2, t/2), ("upper", t, 2*t), ("zero", -t, 0)]:
    check("current_factor_two_"+label, same((t+eps).subs(eps, e), target))

# Direct rational full matrix inversion; no Woodbury formula in comparison.
for q in (2, 3, 5):
    for n in (2, 3, 5, 8):
        h = sp.diag(*([q*q+1]*n))
        for i in range(n-1):
            h[i, i+1] = h[i+1, i] = -q
        inv = h.inv()
        ff, tt = inv[0, 0], inv[0, n-1]/inv[0, 0]
        check(f"direct_q{q}_n{n}_endpoint",
              same(tt, sp.Rational((q*q-1)*q**(n-1), q**(2*n)-1)))
        # Verify a constructive vectorlike Higgs charge factorization.
        charge = sp.zeros(n+1, n)
        charge[0, 0] = q
        for i in range(n-1):
            charge[i+1, i], charge[i+1, i+1] = 1, -q
        charge[n, n-1] = 1
        check(f"higgs_charge_gram_q{q}_n{n}", charge.T*charge == h)
        check(f"higgs_charge_rank_q{q}_n{n}", charge.rank() == n)
        null_phase = sp.Matrix([1]+[-q**j for j in range(1,n+1)])
        check(f"higgs_one_uneaten_chiral_q{q}_n{n}",
              charge.T*null_phase == sp.zeros(n,1) and n+1-charge.rank() == 1)
        # Exact anomaly cancellation for each pair of opposite charge vectors.
        check(f"higgs_vectorlike_linear_q{q}_n{n}",
              all(x+(-x) == 0 for x in charge))
        check(f"higgs_vectorlike_cubic_q{q}_n{n}",
              all(charge[a,i]*charge[a,j]*charge[a,k]
                  +(-charge[a,i])*(-charge[a,j])*(-charge[a,k]) == 0
                  for a in range(n+1) for i in range(n)
                  for j in range(n) for k in range(n)))
        for ratio in (sp.Rational(-1), sp.Rational(1,2), sp.Rational(2)):
            dd = ratio*tt/(ff*(1-tt*tt))
            hp = h.copy()
            hp[0, n-1] += dd
            hp[n-1, 0] += dd
            invp = hp.inv()
            denom = 1+2*dd*ff*tt-dd*dd*ff*ff*(1-tt*tt)
            tag = f"direct_q{q}_n{n}_r{ratio}"
            check(tag+"_self", same(invp[0,0], ff/denom))
            check(tag+"_cross", same(invp[0,n-1]/invp[0,0], tt-dd*ff*(1-tt*tt)))
            check(tag+"_determinant_lemma", same(hp.det()/h.det(), denom))
            # Sylvester's criterion compared to exact two-channel bound.
            sylvester = all(hp[:k,:k].det() > 0 for k in range(1, n+1))
            bound = -1/(ff*(1+tt)) < dd < 1/(ff*(1-tt))
            check(tag+"_SPD", sylvester == bound)

# A tiny cross element must be checked with sufficient precision; a standard
# float endpoint inverse or visual spectrum plot is not a protection test.
mp.mp.dps = 180
q, n = mp.mpf(3), 127
ff = (q**(2*n)-1)/(q**(2*n+2)-1)
tt = (q-q**-1)/(q**n-q**-n)
h = mp.matrix(n)
for i in range(n):
    h[i,i] = q*q+1
    if i+1 < n:
        h[i,i+1] = h[i+1,i] = -q


def close(a, b, tol=mp.mpf("1e-110")):
    return abs(a-b) <= tol*max(abs(a), abs(b), mp.mpf("1e-170"))


# Direct dense LU of the full 127x127 deformed matrix, not endpoint updating.
for r in (mp.mpf("0"), mp.mpf("0.5"), mp.mpf("2")):
    dd = r*tt/(ff*(1-tt*tt))
    hp = h.copy()
    hp[0,n-1] += dd
    hp[n-1,0] += dd
    b = mp.matrix(n,1)
    b[0] = 1
    sol = mp.lu_solve(hp, b)
    denom = 1+2*dd*ff*tt-dd*dd*ff*ff*(1-tt*tt)
    check(f"full127_r{r}_self", close(sol[0], ff/denom))
    check(f"full127_r{r}_cross", close(sol[n-1]/sol[0], tt-dd*ff*(1-tt*tt)))
    check(f"full127_r{r}_residual", mp.norm(hp*sol-b, p=mp.inf) < mp.mpf("1e-170"))

P = mp.mpf("2.435e27")**2
A = mp.zeta(3)/(48*mp.pi**4)*mp.mpf("1e14")**2/P
Lambda = mp.mpf("2e12")
tau_critical = 2*A*Lambda**2/P
bound = tt/(ff*(1-tt*tt))
delta_restore = (tt-tau_critical)/(ff*(1-tt*tt))
# Holding M fixed changes c0 by 1/D; solve that precise threshold separately.
def restored_difference(x):
    dd = x*bound
    denom = 1+2*dd*ff*tt-dd*dd*ff*ff*(1-tt*tt)
    return (tt-dd*ff*(1-tt*tt))/tt - tau_critical/tt*denom
fixed_m_root_fraction = mp.findroot(restored_difference, (mp.mpf("0.1"), mp.mpf("0.2")))
fixed_m_root = fixed_m_root_fraction*bound
check("representative_rH_above_two", tt*P/(A*Lambda**2) > 2)
check("factor_two_lower_not_enough_for_restoration", tt*P/(2*A*Lambda**2) < 2)
check("fixed_Lambda_restoration_boundary", close(tt-delta_restore*ff*(1-tt*tt), tau_critical))
check("fixed_M_restoration_boundary", abs(restored_difference(fixed_m_root_fraction)) < mp.mpf("1e-170"))
check("fixed_M_and_fixed_Lambda_equal_through_110_relative_digits",
      close(fixed_m_root, delta_restore))
check("epsilon_zero_cross_keeps_positive_self", 1-tt*tt > 0)
check("mass_sign_flip_does_not_destroy_SPD", 2*bound < 1/(ff*(1-tt)))

def val(x):
    return mp.nstr(x, 65)

examples = []
for e in (mp.mpf("-1e-50"), -tt, -tt/2, mp.mpf(0), tt, mp.mpf("1e-50")):
    examples.append({"type":"current_leakage", "epsilon":val(e),
                     "tau_eff":val(tt+e), "tau_ratio":val((tt+e)/tt),
                     "zeta_minus_one":val(2*e*tt+e*e)})
for dd in (-bound, bound/2, bound, 2*bound, mp.mpf("1e-50")):
    denom = 1+2*dd*ff*tt-dd*dd*ff*ff*(1-tt*tt)
    tnew = tt-dd*ff*(1-tt*tt)
    examples.append({"type":"mass_bypass", "delta":val(dd),
                     "tau_normalized":val(tnew), "tau_ratio":val(tnew/tt),
                     "D_minus_one":val(denom-1),
                     "is_SPD":bool(-1/(ff*(1+tt)) < dd < 1/(ff*(1-tt)))})

inputs = [ROOT/"experiments/mediator_chain_v01/protocol.md",
          ROOT/"reports/mediator_chain_sources_20260926.md",
          ROOT/"reports/mediator_chain_algebra_20260926.md"]
out = {"scope":"Prescribed endpoint mass/current deformations; no computed radiative generation or universal loop floor",
       "precision_decimal_digits":180, "check_count":len(checks),
       "passed":sum(c["passed"] for c in checks),
       "failed":sum(not c["passed"] for c in checks), "checks":checks,
       "representative":{"q":3,"n":127,"Lambda_eV":"2e12",
                         "f":val(ff),"tau":val(tt),"r_H":val(tt*P/(A*Lambda**2)),
                         "tau_critical_leading_restoration":val(tau_critical),
                         "epsilon_factor_two_lower":val(-tt/2),
                         "epsilon_factor_two_upper":val(tt),
                         "delta_factor_two_lower":val(-bound),
                         "delta_factor_two_upper":val(bound/2),
                         "delta_SPD_lower":val(-1/(ff*(1+tt))),
                         "delta_SPD_upper":val(1/(ff*(1-tt))),
                         "epsilon_restoration_lower_strict":val(tau_critical-tt),
                         "delta_restoration_upper_fixed_Lambda_strict":val(delta_restore),
                         "delta_restoration_upper_fixed_M_strict":val(fixed_m_root),
                         "relative_negative_tau_margin":val(1-tau_critical/tt)},
       "examples":examples,
       "input_hashes":{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},
       "script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
       "interpretation":"Exact algebra and high precision matrix checks only; no statistical or observational claim."}
target = HERE/"chain_locality_breaking_20260926_checks.json"
target.write_text(json.dumps(out, indent=2, ensure_ascii=False)+"\n")
print(json.dumps({"passed":out["passed"],"total":out["check_count"],"failed":out["failed"],"output":str(target)}, indent=2))
sys.exit(1 if out["failed"] else 0)
