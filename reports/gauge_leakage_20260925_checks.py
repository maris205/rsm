#!/usr/bin/env python3
"""Independent leading-log gauge leakage budget, no background reintegration.

Reads prior archived inputs and the soft_0 trajectory only. Does not import
prior production modules. Fixed RG endpoint=100 GeV; fixed CW scale=100 GeV.
This isolates the common scalar soft contribution, not complete gauge mediation.
"""
from pathlib import Path
from datetime import datetime, timezone
import csv
import hashlib
import json
import math
import platform
import time

import mpmath as mp
from scipy.integrate import solve_ivp

ROOT = Path(__file__).resolve().parents[1]
REPORTS = ROOT / 'reports'
PARENT = ROOT / 'experiments/protected_threshold_v01/results'


def digest(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    started = time.monotonic()
    inputs = json.loads((PARENT / 'inputs.json').read_text())
    budget = json.loads((PARENT / 'budgets.json').read_text())
    with (PARENT / 'trajectories.csv').open(newline='') as f:
        rows = [r for r in csv.DictReader(f) if r['case'] == 'soft_0']
    alpha = inputs['alpha_reference']
    g2 = 4 * math.pi * alpha
    uv, ir, mi = 1e13, 1e11, 1e11  # all eV
    L = math.log(uv / ir)
    kappa = 8 * g2 * L / (16 * math.pi**2)
    c = math.sqrt(kappa)
    old = next(r for r in budget['records']
               if r['mass_eV'] == mi and r['mu_over_mi'] == 1)
    scalar_limit = old['msoft_for_force_point1_eV']
    gaugino_limit = scalar_limit / c
    Mlam = 1e9
    d = kappa * Mlam**2
    # Dimensionless RGE with |Mlambda|=1: two charge-conjugate soft masses.
    def rge(t, x):
        S = x[0] - x[1]
        common = -8 * g2 / (16 * math.pi**2)
        return [common + 2 * g2 * S / (16 * math.pi**2),
                common - 2 * g2 * S / (16 * math.pi**2)]
    sol = solve_ivp(rge, [0, -L], [0, 0], method='DOP853',
                    rtol=1e-12, atol=1e-14)
    # Independently reconstruct masses and tree forces from original formula.
    x0 = inputs['chi_i']
    unit = inputs['potential_unit_eV4']
    Wi = inputs['W_i']
    force_per_d = []
    exact_large_force = []
    d_over_s = []
    mp.mp.dps = 90
    md = mp.mpf(str(d))
    for row in rows:
        chi = float(row['chi'])
        s = 0.5 * mi**2 * (1 + (x0 / chi)**2)
        sp = -mi**2 * x0**2 / chi**3
        W = Wi * math.exp(-2 * (chi-x0))
        ls = math.log(s / mi**2)
        tree = 2 * unit * W
        force_per_d.append(abs(sp * ls / (8 * math.pi**2)) / tree)
        ms, msp = mp.mpf(str(s)), mp.mpf(str(sp))
        mmu2 = mp.mpf(str(mi))**2
        # Exact common-soft determinant derivative, no small-d subtraction.
        exact = msp * ((ms+md)*(mp.log((ms+md)/mmu2)-1)
                      - ms*(mp.log(ms/mmu2)-1)) / (8*mp.pi**2)
        exact_large_force.append(float(abs(exact) / mp.mpf(str(tree))))
        d_over_s.append(d/s)
    max_per_d = max(force_per_d)
    reconstructed_limit = math.sqrt(0.1 / max_per_d)
    leading_force = max_per_d * d
    exact_force = max(exact_large_force)
    checks = []
    def check(name, value, tolerance):
        checks.append(dict(name=name, value=float(value), tolerance=tolerance,
                           passed=bool(value <= tolerance)))
    check('coefficient_from_alpha_normalization',
          abs(kappa - 2*alpha*L/math.pi), 1e-16)
    check('rge_scalar_plus_relative', abs(sol.y[0,-1]/kappa-1), 1e-11)
    check('rge_scalar_minus_relative', abs(sol.y[1,-1]/kappa-1), 1e-11)
    check('charge_trace_zero', abs(sol.y[0,-1]-sol.y[1,-1]), 1e-15)
    check('old_force_budget_reconstruction_relative',
          abs(reconstructed_limit/scalar_limit-1), 1e-9)
    check('gaugino_limit_maps_to_force_point1',
          abs(max_per_d*kappa*gaugino_limit**2-0.1), 1e-10)
    # With d/s~2e-6 the exact common CW force is well approximated at its maximum.
    check('large_case_exact_vs_leading_max_relative',
          abs(exact_force/leading_force-1), 1e-3)
    check('large_case_d_over_s_small', max(d_over_s), 1e-4)
    check('gaugino_below_pair_threshold', Mlam/ir, 0.1)
    out = dict(
        run_utc=datetime.now(timezone.utc).isoformat(),
        elapsed_seconds=time.monotonic()-started,
        python=platform.python_version(), mpmath_precision=mp.mp.dps,
        scope='Frozen fixed-scale common-soft gauge channel budget; not full RG-improved potential, no new cosmic ODE or observed limit',
        input_sha256={str(p.relative_to(ROOT)):digest(p) for p in
                      [PARENT/'inputs.json', PARENT/'budgets.json',
                       PARENT/'trajectories.csv', Path(__file__)]},
        settings=dict(q_abs=1, alpha=alpha,g_squared=g2,
                      rg_uv_eV=uv,rg_ir_eV=ir, CW_mu_eV=mi,
                      scalar_soft_UV_eV2=0, fixed_g_and_gaugino=True,
                      pair_initial_supersymmetric_mass_eV=mi),
        results=dict(kappa_d_over_gaugino_mass2=kappa,
                     msoft_over_gaugino_mass=c,
                     prior_msoft_force_point1_limit_eV=scalar_limit,
                     independently_reconstructed_msoft_limit_eV=reconstructed_limit,
                     conditional_gaugino_force_point1_limit_eV=gaugino_limit,
                     illustrative_gaugino_eV=Mlam,
                     illustrative_msoft_eV=c*Mlam,
                     illustrative_common_d_eV2=d,
                     illustrative_leading_force_over_tree=leading_force,
                     illustrative_exact_common_force_over_tree=exact_force,
                     max_d_over_s=max(d_over_s),
                     scalar_force_max_index=force_per_d.index(max_per_d),
                     archived_points=len(rows),
                     generated_b_over_Mlambda_M_leadinglog=-kappa),
        uncomputed=['running g and Mlambda','physical moving decoupling scale',
                    'complete two-loop potential and local counterterms',
                    'induced holomorphic b contribution to total CW force',
                    'full hidden/visible matter embedding and experimental constraints'],
        checks=checks, checks_passed=sum(c['passed'] for c in checks),
        check_count=len(checks))
    target=REPORTS/'gauge_leakage_20260925_checks.json'
    target.write_text(json.dumps(out,indent=2,ensure_ascii=False)+'\n')
    print(json.dumps(out['results'],indent=2))
    print(f"Checks {out['checks_passed']}/{out['check_count']}")
    if out['checks_passed'] != out['check_count']:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
