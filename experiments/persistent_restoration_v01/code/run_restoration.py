#!/usr/bin/env python3
"""Persistent local restoration, correlated soft terms, and conditional late tails."""
from pathlib import Path
from datetime import datetime, timezone
import csv
import hashlib
import json
import platform
import numpy as np
import scipy
from scipy.special import zeta
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parents[1]
REPO = ROOT.parents[1]
OUT = ROOT / 'results'
PHASES = {'zero': (1., 0.), 'quadrature': (0., 1.), 'pi': (-1., 0.)}
MG = np.logspace(-12, -3, 37)
RH = np.logspace(-1, 12, 53)
MI, MV, KK = 1e11, 1e12, 1e14
C5 = float(zeta(3, 1) / (48 * np.pi ** 4))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def dump(path, data):
    path.write_text(json.dumps(data, indent=2, ensure_ascii=False, allow_nan=False) + '\n')


def csvout(path, rows):
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def inputs():
    folder = REPO / 'experiments/sugra_shift_v01/results'
    p = json.loads((folder / 'inputs.json').read_text())
    with (folder / 'trajectories.csv').open() as stream:
        chi = np.array([float(r['chi']) for r in csv.DictReader(stream) if r['case'] == 'shift_axis'])
    p['P'] = p['Mpl_eV'] ** 2
    p['rho'] = 3 * p['omega_lambda'] * p['P'] * p['H_ref_eV'] ** 2
    p['Ui'] = p['potential_unit_eV4'] * p['W_i']
    return p, chi


def kernel(mg, phase, p, chi):
    co, si = PHASES[phase]
    P, F2 = p['P'], p['F_squared_eV2']
    F = np.sqrt(F2)
    U = p['Ui'] * np.exp(-2 * (chi - p['chi_i']))
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
    m0 = (4 * (1 + rr) * U / (3 * P) - 2 * rr * mg * wm * co / F2
          + 2 * rr * mg ** 2 * (2 + rr) * ((co - rw) ** 2 + si ** 2))
    m0x = (-8 * (1 + rr) * U / (3 * P) + 2 * rr * mg * wm * co / F2
           + 4 * rr * mg ** 2 * (2 + rr) * rw * (co - rw))
    return dict(U=U, rw=rw, X=X, kc=kc, s=s, L=L, ds=ds, ft=ft, Btree=Btree,
                mass0=m0, mass0_x=m0x, cos=co, sin=si, rr=rr)


def ratios(v, U):
    return np.hypot(v[0], v[1]) / (2 * U)


def arrays(mg, rh, sigma, p, f):
    P, F2 = p['P'], p['F_squared_eV2']
    F, U = np.sqrt(F2), f['U']
    A = C5 * KK ** 2 / P
    tau = rh * A * MV ** 2 / P
    c = (12 * A + 3 * A ** 2) / (1 + A) ** 2
    D = 6 * sigma * rh * A
    # Stable even at sigma*rH=2; this residual is not full higher-order matching.
    net_c = 6 * A * (2 - sigma * rh) - 3 * A ** 2 * (7 + 4 * A) / (1 + A) ** 2
    beta = -6 * A / (1 + A)
    shift = (rh - 2) + rh * (-2 * f['cos'] * f['rw'] + f['rw'] ** 2)
    mass = f['mass0'] + 12 * U / MV ** 2 + 12 * A * mg ** 2 * shift
    fX2 = 3 * P * f['X']
    radial = 48 * mg ** 2 + (4 * tau ** 2 * fX2 + 2 * tau * U) / MV ** 2
    mixed2 = 16 * tau ** 2 * fX2 * U / MV ** 4
    schur = mass - mixed2 / radial
    J = -6 * np.sqrt(2) * A * F * np.sqrt(U) * mg * f['sin']
    dU = 2 * sigma * U / MV ** 2
    dy_extra = -12 * np.sqrt(2) * sigma * A * (2 + tau) * F * np.sqrt(U) * mg * f['sin'] / MV ** 2
    contact = (f['L'] * f['ds'] * dU + f['s'] * (f['L'] - 1) * np.stack((-2 * dU, dy_extra))) / (8 * np.pi ** 2)
    new_By = 4j * np.sqrt(2) * sigma * np.sqrt(f['s']) * F * np.sqrt(U) / MV ** 2
    B5d = -beta * np.sqrt(f['s']) * f['ft']
    btree_y = 2 * np.real(np.conj(f['Btree']) * new_By) * f['L'] / (16 * np.pi ** 2)
    b5_y = 2 * np.real(np.conj(B5d) * new_By) * f['L'] / (16 * np.pi ** 2)
    source_r = ratios(c * f['kc'], U)
    hidden_r = ratios(-D * f['kc'], U)
    contact_r = ratios(contact, U)
    axis_budget = source_r + hidden_r + contact_r + (np.abs(btree_y) + np.abs(b5_y)) / (2 * U)
    signed = net_c * f['kc'] + contact
    signed[1] += btree_y + b5_y
    good = mass > 0
    valley_x = np.zeros_like(U)
    field_fraction = np.full_like(U, np.nan)
    val_energy = np.full_like(U, np.nan)
    if np.any(good):
        # Analytically remove the U/Λ² terms and near-rH=2 cancellation.
        residual = (2 * f['mass0'] + f['mass0_x']
                    + 24 * A * mg ** 2 * ((rh - 2) - rh * f['cos'] * f['rw']))
        valley_x[good] = J[good] ** 2 * residual[good] / (2 * F2 * mass[good] ** 2)
        field_fraction[good] = np.sqrt(6) * np.abs(J[good]) / (F * mass[good] * MV)
        val_energy[good] = J[good] ** 2 / (2 * F2 * mass[good] * U[good])
    return dict(A=A, tau=tau, c=c, D=D, net_c=net_c, mass=mass, radial=radial, schur=schur,
                source_r=source_r, hidden_r=hidden_r, contact_r=contact_r,
                btree_r=np.abs(btree_y) / (2 * U), b5_r=np.abs(b5_y) / (2 * U),
                axis_budget=axis_budget, budget=axis_budget + np.abs(valley_x) / (2 * U),
                signed_ratio=ratios(signed, U), signed_chi_ratio=np.abs(signed[0]) / (2 * U),
                valley_r=np.abs(valley_x) / (2 * U), field_fraction=field_fraction,
                val_energy=val_energy, soft_fraction=np.abs(dU + D * f['X'] - c * f['X']) / f['s'])


def record(mg, rh, sigma, phase, p, f):
    a = arrays(mg, rh, sigma, p, f)
    stable = bool(np.all(a['mass'] > 0) and np.all(a['schur'] > 0) and np.all(a['radial'] > 0))
    off = 2 * f['rr'] * mg ** 2 * (2 + f['rr']) + 12 * a['A'] * mg ** 2 * (rh - 2)
    aux = np.sqrt(3) * p['Mpl_eV'] * mg / MV ** 2
    M5 = (p['P'] * KK / np.pi) ** (1 / 3)
    hierarchy = max(MI / MV, MV / KK, KK / M5, mg / MV,
                    np.sqrt(np.max(np.abs(a['radial']))) / MV,
                    np.sqrt(np.max(np.abs(a['mass']))) / MV, p['H_ref_eV'] / MV)
    coupling = 2 * max(1, sigma ** 2, a['tau'] ** 2) / (16 * np.pi ** 2)
    budget = float(np.max(a['budget'])) if stable else None
    field = float(np.max(a['field_fraction'])) if stable else None
    passed = (stable and off > 0 and budget <= .1 and field <= .1 and aux <= .1
              and hierarchy <= .1 * (1 + 1e-14) and coupling <= .1)
    return dict(mG_eV=float(mg), rH=float(rh), sigma=float(sigma), phase=phase,
                M_vector_eV=MV, mKK_eV=KK, tau=float(a['tau']), A=float(a['A']),
                source_c=float(a['c']), persistent_soft_coefficient=float(a['D']),
                formal_net_c=float(a['net_c']), min_my2_eV2=float(np.min(a['mass'])),
                min_transverse_schur_eV2=float(np.min(a['schur'])),
                clockoff_my2_leading_eV2=float(off),
                max_radial_mass_eV=float(np.sqrt(np.max(np.abs(a['radial'])))),
                source_CW_force_ratio=float(np.max(a['source_r'])),
                persistent_CW_force_ratio=float(np.max(a['hidden_r'])),
                clock_contact_CW_force_ratio=float(np.max(a['contact_r'])),
                contact_Btree_force_ratio=float(np.max(a['btree_r'])),
                contact_B5d_force_ratio=float(np.max(a['b5_r'])),
                selected_signed_axis_ratio=float(np.max(a['signed_ratio'])),
                selected_signed_chi_ratio=float(np.max(a['signed_chi_ratio'])),
                selected_axis_moduli_sum=float(np.max(a['axis_budget'])),
                selected_feedback_sum_ratio=budget,
                quadratic_valley_force_ratio=float(np.max(a['valley_r'])) if stable else None,
                max_valley_energy_over_U=float(np.max(a['val_energy'])) if stable else None,
                valley_boundary_fraction=field, FT_over_MV2=float(aux),
                max_hierarchy_ratio=float(hierarchy), max_coupling_loop_parameter=float(coupling),
                max_soft_over_charged_mass2=float(np.max(a['soft_fraction'])),
                finite_interval_local_pass=stable, clockoff_local_pass=bool(off > 0),
                hierarchy_pass=bool(hierarchy <= .1 * (1 + 1e-14)),
                auxiliary_pass=bool(aux <= .1), perturbative_pass=bool(coupling <= .1),
                stated_conditions_pass=bool(passed),
                leading_cancellation_branch=bool(sigma * rh == 2))


def first_crossing(fn, start, stop):
    xs = np.linspace(start, stop, 321)
    last = fn(xs[0])
    if last >= 0:
        return float(start)
    for left, right in zip(xs[:-1], xs[1:]):
        current = fn(right)
        if current >= 0:
            return float(brentq(fn, left, right, xtol=1e-11))
        last = current
    return None


def future(mg, rh, sigma, phase, p, start):
    stop = start + 80

    def evaluate_at(x):
        f = kernel(mg, phase, p, np.array([x]))
        return arrays(mg, rh, sigma, p, f)

    massroot = first_crossing(lambda x: -evaluate_at(x)['mass'][0], start, stop)
    valid_stop = stop if massroot is None else max(start, massroot - 1e-9)
    budgetroot = first_crossing(lambda x: np.log10(max(evaluate_at(x)['budget'][0], 1e-300)) + 1, start, valid_stop)
    signedroot = first_crossing(lambda x: np.log10(max(evaluate_at(x)['signed_chi_ratio'][0], 1e-300)) + 1, start, valid_stop)
    roots = [r for r in (massroot, budgetroot) if r is not None]
    return dict(mG_eV=mg, rH=rh, sigma=sigma, phase=phase,
                chi_start=start, chi_search_end=stop, chi_local_mass_failure=massroot,
                chi_moduli_budget_10percent=budgetroot,
                chi_signed_chi_budget_10percent=signedroot,
                chi_first_stated_failure=min(roots) if roots else None,
                ordering_note=('curvature_or_feedback_within_chi_1e-9' if massroot is not None and phase == 'quadrature' and budgetroot is None
                               else 'Only positive-curvature interval searched for feedback roots.'),
                interpretation='Fixed-parameter coordinate extrapolation only; selected matching, no cosmic time prediction.')


def main():
    p, chi = inputs()
    grid, reps, limits, curves, checks = [], [], [], [], []

    def check(name, ok, detail=''):
        checks.append(dict(name=name, passed=bool(ok), detail=detail))

    for mg in MG:
        for phase in ('zero', 'quadrature'):
            f = kernel(mg, phase, p, chi)
            for rh in RH:
                grid.append(record(mg, rh, 1, phase, p, f))
    A = C5 * KK ** 2 / p['P']
    for rh in (-4, 0, 2, 4, 1 / A, p['P'] / (A * MV ** 2)):
        for sigma in (-1, 0, 1):
            for phase in PHASES:
                reps.append(record(1e-6, rh, sigma, phase, p, kernel(1e-6, phase, p, chi)))
    xs = np.linspace(float(chi[-1]), float(chi[-1]) + 30, 601)
    for rh in (0, 2, 4, 1e6):
        for phase in ('zero', 'quadrature'):
            limits.append(future(1e-6, rh, 1, phase, p, float(chi[-1])))
        a = arrays(1e-6, rh, 1, p, kernel(1e-6, 'zero', p, xs))
        for i, x in enumerate(xs):
            curves.append(dict(rH=rh, chi=float(x), local_my2_eV2=float(a['mass'][i]),
                               moduli_budget=float(a['budget'][i]),
                               signed_chi_budget=float(a['signed_chi_ratio'][i])))
    sample = next(r for r in reps if r['rH'] == 4 and r['sigma'] == 1 and r['phase'] == 'zero')
    sampleq = next(r for r in reps if r['rH'] == 4 and r['sigma'] == 1 and r['phase'] == 'quadrature')
    check('grid_count', len(grid) == 3922)
    check('representative_count', len(reps) == 54)
    check('future_limits_count', len(limits) == 8)
    check('both_grid_outcomes', {r['stated_conditions_pass'] for r in grid} == {True, False})
    check('sample_persistent_local_pass', sample['stated_conditions_pass'])
    check('sample_generic_phase_pass', sampleq['stated_conditions_pass'])
    for r in reps:
        tag = f"{r['rH']}_{r['sigma']}_{r['phase']}"
        check('contact_to_source_ratio_' + tag,
              np.isclose(r['persistent_CW_force_ratio'], abs(r['sigma'] * r['rH'] / 2) * r['source_CW_force_ratio'], rtol=1e-12, atol=0))
        if r['rH'] <= 0:
            check('no_persistence_' + tag, not r['clockoff_local_pass'])
        if r['rH'] >= 1 / A and r['sigma'] != 0:
            check('large_coupling_feedback_failure_' + tag, r['selected_axis_moduli_sum'] > .1)
    check('sample_future_radiative_failure', next(r for r in limits if r['rH'] == 4 and r['phase'] == 'zero')['chi_moduli_budget_10percent'] is not None)
    for name, rows in [('scan', grid), ('representatives', reps), ('future_limits', limits), ('future_curves', curves)]:
        csvout(OUT / (name + '.csv'), rows)
    dump(OUT / 'inputs.json', dict(parent_parameters=p, chi_first=float(chi[0]), chi_last=float(chi[-1]),
                                  chi_points=len(chi), mG_grid_eV=MG.tolist(), rH_grid=RH.tolist(),
                                  M_vector_eV=MV, mKK_eV=KK, Mi_eV=MI, C5=C5,
                                  all_matching_conditions='Fixed mu=Mi, inherited nu0=0, frozen canonical current coefficients.'))
    summary = dict(created_utc=datetime.now(timezone.utc).isoformat(), python=platform.python_version(),
                   numpy=np.__version__, scipy=scipy.__version__, checks=checks,
                   checks_passed=sum(c['passed'] for c in checks), check_count=len(checks),
                   sample_aligned=sample, sample_quadrature=sampleq, future_limits=limits,
                   grid_pass_counts={phase: sum(r['stated_conditions_pass'] for r in grid if r['phase'] == phase)
                                     for phase in ('zero', 'quadrature')},
                   source_hashes={str(q.relative_to(REPO)): sha(q) for q in [Path(__file__), ROOT / 'protocol.md',
                                  REPO / 'experiments/sugra_shift_v01/results/inputs.json',
                                  REPO / 'experiments/sugra_shift_v01/results/trajectories.csv']},
                   scope=['The tiny positive tau is an unfixed input, not a derived naturalness result.',
                          'Persistent transverse restoration does not remove the selected quantum clock slope.',
                          'Only declared local conditions and finite matching sectors; no physical time extrapolation.'])
    dump(OUT / 'summary.json', summary)
    print('Main', summary['checks_passed'], '/', summary['check_count'], summary['grid_pass_counts'])
    print('Aligned', json.dumps(sample))
    print('Future', json.dumps(limits))
    for c in checks:
        if not c['passed']:
            print('FAILED', c)
    raise SystemExit(0 if summary['checks_passed'] == summary['check_count'] else 1)


if __name__ == '__main__':
    main()
