#!/usr/bin/env python3
"""Finite-interval clock stabilization and the contacts of a shared boundary current.

Selected Wilsonian sectors only; no complete heavy spectrum or cosmic solution.
"""
from pathlib import Path
from datetime import datetime, timezone
import csv
import hashlib
import json
import platform
import numpy as np
import scipy
from scipy.special import zeta

ROOT = Path(__file__).resolve().parents[1]
REPO = ROOT.parents[1]
OUT = ROOT / 'results'
PHASES = {'zero': (1., 0.), 'quadrature': (0., 1.), 'pi': (-1., 0.)}
MG = np.logspace(-12, 0, 49)
MV = np.logspace(11, 14, 49)
MI = 1e11
C5 = float(zeta(3, 1) / (48 * np.pi ** 4))


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def dump(p, obj):
    p.write_text(json.dumps(obj, indent=2, ensure_ascii=False, allow_nan=False) + '\n')


def csvout(p, rows):
    with p.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def inputs():
    folder = REPO / 'experiments/sugra_shift_v01/results'
    p = json.loads((folder / 'inputs.json').read_text())
    with (folder / 'trajectories.csv').open() as stream:
        chi = np.array([float(r['chi']) for r in csv.DictReader(stream) if r['case'] == 'shift_axis'])
    P = p['Mpl_eV'] ** 2
    p['P'] = P
    p['rho'] = 3 * p['omega_lambda'] * P * p['H_ref_eV'] ** 2
    U = p['potential_unit_eV4'] * p['W_i'] * np.exp(-2 * (chi - p['chi_i']))
    return p, chi, U


def kernel(mg, phase, p, chi, U):
    co, si = PHASES[phase]
    P, F2 = p['P'], p['F_squared_eV2']
    F = np.sqrt(F2)
    wm = F * np.sqrt(U / 2)
    rw = wm / (mg * P)
    X = mg ** 2 * ((co - rw) ** 2 + si ** 2)
    dX = mg ** 2 * np.stack((2 * rw * (co - rw), 2 * rw * si / 3))
    rel = p['xi'] * (chi ** -2 - p['chi_i'] ** -2) / (1 + p['xi'] / p['chi_i'] ** 2)
    s = MI ** 2 * (1 + rel)
    L = np.log1p(rel)
    r = -p['xi'] / (chi * (chi ** 2 + p['xi']))
    ds = np.stack((2 * r * s, np.zeros_like(s)))
    kc = -(X * L * ds + s * (L - 1) * dX) / (8 * np.pi ** 2)
    ft = mg * (co - rw - 1j * si)
    eta = p['rho'] / (6 * mg ** 2 * P)
    Btree = np.sqrt(s) * (np.sqrt(2 * U) / F * r + 2 * eta * ft)
    rr = p['rho'] / (3 * P * mg ** 2)
    mass0 = (4 * (1 + rr) * U / (3 * P) - 2 * rr * mg * wm * co / F2
             + 2 * rr * mg ** 2 * (2 + rr) * ((co - rw) ** 2 + si ** 2))
    mass0_x = (-8 * (1 + rr) * U / (3 * P) + 2 * rr * mg * wm * co / F2
               + 4 * rr * mg ** 2 * (2 + rr) * rw * (co - rw))
    return dict(kc=kc, s=s, L=L, ds=ds, X=X, ft=ft, Btree=Btree,
                mass0=mass0, mass0_x=mass0_x, sin=si)


def normratio(v, U):
    return np.hypot(v[0], v[1]) / (2 * U)


def evaluate(mg, mv, kk, sigma, phase, p, chi, U, f):
    P, F2 = p['P'], p['F_squared_eV2']
    F = np.sqrt(F2)
    A = C5 * kk ** 2 / P
    c = (12 * A + 3 * A ** 2) / (1 + A) ** 2
    beta = -6 * A / (1 + A)
    lam = mv  # Declared gZ=sqrt(2), with physical Stueckelberg mass normalization.
    mass = f['mass0'] - 24 * A * mg ** 2 + 12 * U / lam ** 2
    stable = bool(np.min(mass) > 0)
    J = -6 * np.sqrt(2) * A * F * np.sqrt(U) * mg * f['sin']
    d = 2 * sigma * U / lam ** 2
    dx = -2 * d
    dy = 2 * sigma * J / lam ** 2
    dgrad = np.stack((dx, dy))
    contact_force = (f['L'] * f['ds'] * d + f['s'] * (f['L'] - 1) * dgrad) / (8 * np.pi ** 2)
    new_By = 4j * np.sqrt(2) * sigma * np.sqrt(f['s']) * F * np.sqrt(U) / lam ** 2
    new_Btree_force_y = 2 * np.real(np.conj(f['Btree']) * new_By) * f['L'] / (16 * np.pi ** 2)
    # B contains -2 F^T d_t ln(hQ); beta=+2 d_t ln(hQ).
    B5d = -beta * np.sqrt(f['s']) * f['ft']
    new_B5d_force_y = 2 * np.real(np.conj(B5d) * new_By) * f['L'] / (16 * np.pi ** 2)
    cr = normratio(c * f['kc'], U)
    qr = normratio(contact_force, U)
    br = np.abs(new_Btree_force_y) / (2 * U)
    b5r = np.abs(new_B5d_force_y) / (2 * U)
    if stable:
        ystar = -J / (F2 * mass)
        val = -J ** 2 / (2 * F2 * mass)
        # The large 12U/Lambda² terms cancel analytically in 2m²+(m²)'.
        residual = 2 * f['mass0'] + f['mass0_x'] - 48 * A * mg ** 2
        val_x = J ** 2 * residual / (2 * F2 * mass ** 2)
        vr = np.abs(val_x) / (2 * U)
        field_fraction = float(np.max(np.sqrt(6) * F * np.abs(ystar) / lam))
        valley_force = float(np.max(vr))
        valley_energy = float(np.max(np.abs(val) / U))
        feedback = float(np.max(cr + qr + br + b5r + vr))
    else:
        field_fraction = valley_force = valley_energy = feedback = None
    M5 = (P * kk / np.pi) ** (1 / 3)
    auxiliary = max(np.sqrt(3) * p['Mpl_eV'] * mg, float(np.sqrt(np.max(U)))) / mv ** 2
    hierarchy = max(MI / mv, mv / kk, kk / M5, mg / mv,
                    np.sqrt(48) * mg / mv, np.sqrt(np.max(np.abs(mass))) / mv, p['H_ref_eV'] / mv)
    loop_parameter = 2 * max(1, sigma ** 2) / (16 * np.pi ** 2)
    Ucrit = (24 * A * mg ** 2 - 4 * p['rho'] / (3 * P)) / (12 / lam ** 2 + 4 / (3 * P))
    chicrit = float(p['chi_i'] + .5 * np.log(U[0] / Ucrit)) if Ucrit > 0 else None
    necessary = (stable and feedback <= .1 and field_fraction <= .1
                 and auxiliary <= .1 * (1 + 1e-14) and hierarchy <= .1 * (1 + 1e-14)
                 and loop_parameter <= .1)
    return dict(mG_eV=float(mg), M_vector_eV=float(mv), mKK_eV=float(kk),
                sigma=float(sigma), phase=phase, Lambda_eV=float(lam), A_5d=float(A),
                c_Q=float(c), beta_Q=float(beta),
                min_my2_eV2=float(np.min(mass)), max_my2_eV2=float(np.max(mass)),
                min_my2_over_Href2=float(np.min(mass) / p['H_ref_eV'] ** 2),
                clockoff_my2_leading_eV2=float(4 * p['rho'] / (3 * P) - 24 * A * mg ** 2),
                RSS_c_2loop_force_ratio=float(np.max(cr)),
                contact_d_1loop_force_ratio=float(np.max(qr)),
                contact_Btree_1loop_force_ratio=float(np.max(br)),
                contact_B5d_2loop_force_ratio=float(np.max(b5r)),
                quadratic_valley_force_ratio=valley_force,
                quadratic_valley_energy_over_U=valley_energy,
                valley_metric_boundary_fraction=field_fraction,
                selected_feedback_sum_ratio=feedback,
                known_axis_feedback_without_valley=float(np.max(cr + qr + br + b5r)),
                charged_over_vector=float(MI / mv), vector_over_KK=float(mv / kk),
                KK_over_M5=float(kk / M5), max_hierarchy_ratio=float(hierarchy),
                FTcanonical_over_vector2=float(np.sqrt(3) * p['Mpl_eV'] * mg / mv ** 2),
                max_auxiliary_over_vector2=float(auxiliary),
                max_coupling_loop_parameter=float(loop_parameter),
                local_positive_curvature=stable, hierarchy_pass=bool(hierarchy <= .1 * (1 + 1e-14)),
                auxiliary_pass=bool(auxiliary <= .1 * (1 + 1e-14)),
                perturbative_coupling_pass=bool(loop_parameter <= .1),
                necessary_conditions_pass=bool(necessary),
                Ucrit_leading_eV4=float(Ucrit) if Ucrit > 0 else None,
                chi_critical_leading=chicrit,
                chi_gap_after_last=chicrit - float(chi[-1]) if chicrit is not None else None)


def main():
    p, chi, U = inputs()
    OUT.mkdir(exist_ok=True)
    grid, reps, sensitivity, checks = [], [], [], []

    def check(name, passed, detail=''):
        checks.append(dict(name=name, passed=bool(passed), detail=detail))

    for mg in MG:
        for phase in ('zero', 'quadrature'):
            f = kernel(mg, phase, p, chi, U)
            for mv in MV:
                grid.append(evaluate(mg, mv, 1e14, 1, phase, p, chi, U, f))
    for mg, mv, kk in [(1., 1e9, 1e12), (1., 1e12, 1e12),
                       (1e-6, 1e12, 1e14), (1e-5, 1e12, 1e14),
                       (1e-4, 1e12, 1e14), (1e-6, 1e13, 1e14)]:
        for phase in PHASES:
            f = kernel(mg, phase, p, chi, U)
            for sigma in (-1, 0, 1):
                reps.append(evaluate(mg, mv, kk, sigma, phase, p, chi, U, f))
    for sigma in np.linspace(-3, 3, 25):
        for phase in ('zero', 'quadrature'):
            f = kernel(1e-6, phase, p, chi, U)
            sensitivity.append(evaluate(1e-6, 1e12, 1e14, sigma, phase, p, chi, U, f))
    frontier = []
    for kk in (1e12, 1e13, 1e14):
        A = C5 * kk ** 2 / p['P']
        for mv in MV:
            frontier.append(dict(mKK_eV=kk, M_vector_eV=float(mv),
                                  local_curvature_mG_max_eV=float(np.sqrt((4 * (np.min(U) + p['rho']) / (3 * p['P']) + 12 * np.min(U) / mv ** 2) / (24 * A))),
                                  auxiliary_mG_max_eV=float(.1 * mv ** 2 / (np.sqrt(3) * p['Mpl_eV'])),
                                  two_mass_gaps_pass=bool(MI / mv <= .1 * (1 + 1e-14) and mv / kk <= .1 * (1 + 1e-14))))

    check('grid_count', len(grid) == 4802)
    check('representative_count', len(reps) == 54)
    check('sigma_count', len(sensitivity) == 50)
    check('both_outcomes_retained', {r['necessary_conditions_pass'] for r in grid} == {False, True})
    for r in reps:
        tag = f"{r['mG_eV']}_{r['M_vector_eV']}_{r['mKK_eV']}_{r['sigma']}_{r['phase']}"
        if r['sigma'] == 0:
            check('zero_contact_' + tag, r['contact_d_1loop_force_ratio'] == r['contact_Btree_1loop_force_ratio'] == r['contact_B5d_2loop_force_ratio'] == 0)
        if r['phase'] in ('zero', 'pi'):
            check('aligned_B_interference_' + tag, r['contact_Btree_1loop_force_ratio'] == r['contact_B5d_2loop_force_ratio'] == 0)
        check('declared_mass_normalization_' + tag, r['Lambda_eV'] == r['M_vector_eV'])
    sample = next(r for r in reps if r['mG_eV'] == 1e-6 and r['M_vector_eV'] == 1e12 and r['sigma'] == 1 and r['phase'] == 'zero')
    quadrature = next(r for r in reps if r['mG_eV'] == 1e-6 and r['M_vector_eV'] == 1e12 and r['sigma'] == 1 and r['phase'] == 'quadrature')
    check('finite_interval_example_passes', sample['necessary_conditions_pass'])
    check('generic_phase_example_passes', quadrature['necessary_conditions_pass'])
    check('future_failure_retained', sample['clockoff_my2_leading_eV2'] < 0 and sample['chi_gap_after_last'] > 0)
    old = [r for r in reps if r['mG_eV'] == 1]
    check('original_one_eV_examples_fail_joint_conditions', not any(r['necessary_conditions_pass'] for r in old))
    check('higher_sigma_perturbative_control', not sensitivity[-1]['perturbative_coupling_pass'])
    csvout(OUT / 'scan.csv', grid)
    csvout(OUT / 'representatives.csv', reps)
    csvout(OUT / 'sigma_sensitivity.csv', sensitivity)
    csvout(OUT / 'conditional_frontiers.csv', frontier)
    dump(OUT / 'inputs.json', dict(parent_parameters=p, chi_points=len(chi), chi_first=float(chi[0]),
                                  chi_last=float(chi[-1]), Umin_eV4=float(np.min(U)), Umax_eV4=float(np.max(U)),
                                  masses_eV=MG.tolist(), vector_masses_eV=MV.tolist(), C5=C5,
                                  gZ=float(np.sqrt(2)), sigma_sensitivity=np.linspace(-3, 3, 25).tolist(),
                                  primary_KK_eV=1e14, hierarchy_budget=.1, auxiliary_budget=.1, force_budget=.1,
                                  status='Conditional finite-interval Wilsonian model, not full UV or observed universe.'))
    summary = dict(created_utc=datetime.now(timezone.utc).isoformat(), python=platform.python_version(),
                   numpy=np.__version__, scipy=scipy.__version__, checks=checks,
                   checks_passed=sum(c['passed'] for c in checks), check_count=len(checks),
                   grid_rows=len(grid), representative_rows=len(reps), sample_aligned=sample,
                   sample_quadrature=quadrature,
                   pass_counts={phase: sum(r['necessary_conditions_pass'] for r in grid if r['phase'] == phase)
                                for phase in ('zero', 'quadrature')},
                   source_hashes={str(path.relative_to(REPO)): sha(path) for path in [Path(__file__), ROOT / 'protocol.md',
                                  REPO / 'experiments/sugra_shift_v01/results/inputs.json',
                                  REPO / 'experiments/sugra_shift_v01/results/trajectories.csv']},
                   scope=['specified local slices and quadratic valley only', 'selected CW sectors with fixed finite boundary',
                          'heavy auxiliary hierarchy is a conservative diagnostic',
                          'quartic protection disappears as U approaches zero',
                          'no derivation of physical Riemann interactions or inverse-log exponent'])
    dump(OUT / 'summary.json', summary)
    print('Main checks', summary['checks_passed'], '/', summary['check_count'])
    print('Rows', len(grid), len(reps), 'pass counts', summary['pass_counts'])
    print('Aligned example', json.dumps(sample))
    print('Quadrature example', json.dumps(quadrature))
    for c in checks:
        if not c['passed']:
            print('FAILED', c)
    raise SystemExit(0 if summary['checks_passed'] == summary['check_count'] else 1)


if __name__ == '__main__':
    main()
