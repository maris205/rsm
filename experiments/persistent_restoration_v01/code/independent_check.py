#!/usr/bin/env python3
"""Independent implementation of the declared persistent-restoration approximation.

Reads the protocol and historical inputs, never the current main source.  The
comparison mode additionally reads main CSV outputs.  This is an internal AI
audit of a conditional model, not an observation or external peer review.
"""
from pathlib import Path
from datetime import datetime, timezone
import argparse
import csv
import hashlib
import json
import platform

import mpmath as mp
import numpy as np
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parents[1]
REPO = ROOT.parents[1]
OUT = ROOT / 'results'
PARENT = REPO / 'experiments/sequestering_test_v01/results/inputs.json'
COORD = REPO / 'experiments/sugra_shift_v01/results/trajectories.csv'
PHASE = {'zero': (1, 0), 'quadrature': (0, 1), 'pi': (-1, 0)}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def dump(name, value):
    (OUT / name).write_text(json.dumps(value, indent=2, ensure_ascii=False) + '\n')


def csvwrite(name, rows):
    with (OUT / name).open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator='\n')
        w.writeheader()
        w.writerows(rows)


def parameters(high=False):
    raw = json.loads(PARENT.read_text())['parent_parameters']
    cast = (lambda x: mp.mpf(str(x))) if high else float
    p = {k: cast(v) for k, v in raw.items() if isinstance(v, (int, float))}
    p['P'] = p['Mpl_eV'] ** 2
    p['F'] = (mp.sqrt if high else np.sqrt)(p['F_squared_eV2'])
    p['Ui'] = p['potential_unit_eV4'] * p['W_i']
    p['rho'] = p['rho_lambda_eV4']
    p['C5'] = mp.zeta(3) / (48 * mp.pi ** 4) if high else float(mp.zeta(3) / (48 * mp.pi ** 4))
    with COORD.open() as f:
        chi = [cast(r['chi']) for r in csv.DictReader(f) if r['case'] == 'shift_axis']
    return p, chi


def point(chi, mg, rh, sigma, phase, p, high=False, mv=1e12, kk=1e14):
    """Jets from the declared finite CW polynomial and local quadratic valley."""
    b = mp if high else np
    pi, exp, sqrt, log = b.pi, b.exp, b.sqrt, b.log
    co, si = PHASE[phase]
    P, F, F2 = p['P'], p['F'], p['F_squared_eV2']
    A = p['C5'] * kk ** 2 / P
    tau = rh * A * mv ** 2 / P
    U = p['Ui'] * exp(-2 * (chi - p['chi_i']))
    s = p['holomorphic_mass_i_eV'] ** 2 * (1 + p['xi'] / chi ** 2) / (1 + p['xi'] / p['chi_i'] ** 2)
    M = sqrt(s)
    r = -p['xi'] / (chi * (chi ** 2 + p['xi']))
    sx = 2 * r * s
    L = log(s / p['holomorphic_mass_i_eV'] ** 2)
    c = (12 * A + 3 * A ** 2) / (1 + A) ** 2
    beta = -6 * A / (1 + A)
    D = 6 * sigma * rh * A
    dc = 6 * A * (sigma * rh - 2) + 3 * A ** 2 * (7 + 4 * A) / (1 + A) ** 2
    wm = F * sqrt(U / 2)
    rw = wm / (mg * P)
    fre, fim = mg * (co - rw), -mg * si
    X = fre ** 2 + fim ** 2
    Xx = 2 * fre * mg * rw
    Xy = 2 * mg * wm * si / (3 * P)
    pp, px, py = s * X, sx * X + s * Xx, s * Xy
    kernelx = -s * (2 * r * X * L + Xx * (L - 1)) / (8 * pi ** 2)
    kernely = -py * (L - 1) / (8 * pi ** 2)
    sourcex, sourcey = c * kernelx, c * kernely
    persistx, persisty = -D * kernelx, -D * kernely
    d = 2 * sigma * U / mv ** 2
    dx = -2 * d
    dy = -12 * sqrt(2) * sigma * A * (2 + tau) * F * sqrt(U) * mg * si / mv ** 2
    contactx = (sx * d * L + s * dx * (L - 1)) / (8 * pi ** 2)
    contacty = s * dy * (L - 1) / (8 * pi ** 2)
    eta = p['rho'] / (6 * P * mg ** 2)
    B0re = M * (sqrt(2 * U) * r / F + 2 * eta * fre)
    B0im = M * 2 * eta * fim
    dByim = 4 * sqrt(2) * sigma * M * F * sqrt(U) / mv ** 2
    Btree = B0im * dByim * L / (8 * pi ** 2)
    B5d = (-beta * M * fim) * dByim * L / (8 * pi ** 2)
    rr = p['rho'] / (3 * P * mg ** 2)
    b1 = 4 * (1 + rr) * U / (3 * P)
    b2 = -2 * rr * mg * wm * co / F2
    b3 = 2 * rr * mg ** 2 * (2 + rr) * ((co - rw) ** 2 + si ** 2)
    base = b1 + b2 + b3
    restoration = 12 * U / mv ** 2
    net_mass = 12 * A * mg ** 2 * ((rh - 2) + rh * (-2 * co * rw + rw ** 2))
    mass2 = base + restoration + net_mass
    massx = (-2 * b1 - b2 + 2 * rr * mg ** 2 * (2 + rr) * (2 * co * rw - 2 * rw ** 2)
             - 2 * restoration + 12 * A * mg ** 2 * rh * (2 * co * rw - 2 * rw ** 2))
    Jy = -6 * sqrt(2) * A * F * sqrt(U) * mg * si
    residual = b2 + 4 * rr * mg ** 2 * (2 + rr) * (1 - co * rw) + 24 * A * mg ** 2 * ((rh - 2) - rh * co * rw)
    with np.errstate(divide='ignore', invalid='ignore'):
        valley_y = -Jy / (F2 * mass2)
        valley_v = -Jy ** 2 / (2 * F2 * mass2)
        valley_x = 36 * A ** 2 * mg ** 2 * si ** 2 * U * residual / mass2 ** 2
    delta_mX2 = (12 * tau ** 2 * P * X + 2 * tau * U) / mv ** 2
    return locals()


def summarize(mg, rh, sigma, phase, p, chis, mv=1e12, kk=1e14):
    x = point(np.asarray(chis), mg, rh, sigma, phase, p, mv=mv, kk=kk)
    U = x['U']
    ratios = {k: np.hypot(x[k + 'x'], x[k + 'y']) / (2 * U) for k in ('source', 'persist', 'contact')}
    ratios['Btree'] = np.abs(x['Btree']) / (2 * U)
    ratios['B5d'] = np.abs(x['B5d']) / (2 * U)
    axis_budget = sum(ratios.values())
    radial = 48 * mg ** 2 + x['delta_mX2']
    schur = x['mass2'] - 48 * x['tau'] ** 2 * p['P'] * x['X'] * U / (mv ** 4 * radial)
    stable = bool(np.all(x['mass2'] > 0) and np.all(radial > 0) and np.all(schur > 0))
    ratios['valley'] = np.abs(x['valley_x']) / (2 * U)
    budget = axis_budget + ratios['valley']
    netx = -x['dc'] * x['kernelx'] + x['contactx']
    nety = -x['dc'] * x['kernely'] + x['contacty'] + x['Btree'] + x['B5d']
    signed = np.hypot(netx, nety) / (2 * U)
    M5 = (p['P'] * kk / np.pi) ** (1 / 3)
    mass_hierarchy = max(p['holomorphic_mass_i_eV'] / mv, mv / kk, kk / M5,
                         np.sqrt(np.max(radial)) / mv,
                         mg / mv, np.sqrt(np.max(np.abs(x['mass2']))) / mv, p['H_ref_eV'] / mv)
    FT = np.sqrt(3) * p['Mpl_eV'] * mg / mv ** 2
    FZ = np.sqrt(np.max(U)) / mv ** 2
    perturb = 2 * max(1., sigma ** 2, x['tau'] ** 2) / (16 * np.pi ** 2)
    metric = float(np.max(np.sqrt(6) * p['F'] * np.abs(x['valley_y']) / mv)) if stable else None
    rr = p['rho'] / (3 * p['P'] * mg ** 2)
    asymptotic = 2 * rr * mg ** 2 * (2 + rr) + 12 * x['A'] * mg ** 2 * (rh - 2)
    fullbudget = float(np.max(budget)) if stable else None
    out = dict(mG_eV=mg, rH=rh, sigma=sigma, phase=phase, M_vector_eV=mv, mKK_eV=kk,
               tau=float(x['tau']), A=float(x['A']), c=float(x['c']), D=float(x['D']),
               min_mass2=float(np.min(x['mass2'])), max_mass2=float(np.max(x['mass2'])),
               axis_budget=float(np.max(axis_budget)), budget=fullbudget,
               signed=float(np.max(signed)), signed_chi=float(np.max(np.abs(netx) / (2 * U))),
               signed_with_valley=float(np.max(np.hypot(netx + x['valley_x'], nety) / (2 * U))) if stable else None,
               valley_energy=float(np.max(np.abs(x['valley_v']) / U)) if stable else None,
               radial_mass=float(np.sqrt(np.max(radial))), min_schur=float(np.min(schur)),
               formal_net_c=float(-x['dc']), soft_fraction=float(np.max(np.abs(x['dc'] * x['X'] + x['d']) / x['s'])),
               metric_fraction=metric, FT_ratio=float(FT), FZ_ratio=float(FZ),
               hierarchy=float(mass_hierarchy), perturbativity=float(perturb),
               stable=stable, asymptotic_mass2=float(asymptotic), asymptotic_positive=bool(asymptotic > 0),
               mass_scale=float(np.max(np.abs(x['base']) + x['restoration'] + abs(24 * x['A'] * mg ** 2) + abs(12 * rh * x['A'] * x['X']))))
    for k, v in ratios.items():
        out[k + '_ratio'] = float(np.max(v)) if (k != 'valley' or stable) else None
    out['necessary'] = bool(stable and asymptotic > 0 and fullbudget <= .1 and max(FT, FZ) <= .1
                            and mass_hierarchy <= .1 * (1 + 1e-14) and metric <= .1 and perturb <= .1)
    return out


def future(mg, rh, sigma, phase, p, lastchi, stop=None):
    """First local-curvature or conservative-budget failure in the stated interval."""
    stop = lastchi + 80 if stop is None else stop
    def vals(chi):
        v = point(float(chi), mg, rh, sigma, phase, p)
        axis = (np.hypot(v['sourcex'], v['sourcey']) + np.hypot(v['persistx'], v['persisty'])
                + np.hypot(v['contactx'], v['contacty']) + abs(v['Btree']) + abs(v['B5d'])) / (2 * v['U'])
        total = axis + abs(v['valley_x']) / (2 * v['U']) if v['mass2'] > 0 else np.nan
        return v['mass2'], total
    points = np.linspace(lastchi, stop, 1601)
    previous = vals(points[0])
    if previous[0] <= 0:
        return dict(first_chi=lastchi, cause='already_nonpositive_curvature', root_residual=0.)
    if previous[1] >= .1:
        return dict(first_chi=lastchi, cause='already_feedback_exceeds_limit', root_residual=0.)
    for a, b in zip(points[:-1], points[1:]):
        now = vals(b)
        roots = []
        if now[0] <= 0:
            root = brentq(lambda xx: vals(xx)[0], a, b, xtol=1e-12)
            roots.append((root, 'curvature_or_feedback_within_root_resolution', float(vals(root)[0])))
            # The local quadratic valley can diverge arbitrarily close to the
            # mass zero. Never turn its excluded side into an infinite budget.
            b = root - 1e-10
            now = vals(b) if b > a else (0., np.nan)
        if np.isfinite(now[1]) and now[1] >= .1:
            root = brentq(lambda xx: vals(xx)[1] - .1, a, b, xtol=1e-12)
            roots.append((root, 'conservative_feedback', vals(root)[1] - .1))
        if roots:
            root, cause, residual = min(roots)
            return dict(first_chi=float(root), cause=cause, root_residual=float(residual))
        previous = now
    return dict(first_chi=None, cause='none_in_search_interval', root_residual=None)


def derive():
    checks, rows, earlier = [], [], {}
    def check(label, error, tolerance):
        checks.append(dict(name=label, error=mp.nstr(error, 35), tolerance=str(tolerance), passed=bool(error <= tolerance)))
    for digits in (160, 220):
        mp.mp.dps = digits
        p, chis = parameters(True)
        for mt in ('1e-12', '1e-6', '1e-3'):
            mg = mp.mpf(mt)
            for rh in (0, 2, 4, 100):
                for phase in PHASE:
                    for index in (0, 1024):
                        chi = chis[index]
                        v = point(chi, mg, rh, 1, phase, p, True, mp.mpf('1e12'), mp.mpf('1e14'))
                        label = f'{digits}:{mt}:{rh}:{phase}:{index}'
                        # Direct differentiation of each independent finite potential.
                        for name, factor in (('source', -v['c']), ('persist', v['D'])):
                            def potential(xx):
                                w = point(xx, mg, rh, 1, phase, p, True, mp.mpf('1e12'), mp.mpf('1e14'))
                                return factor * w['s'] * w['X'] * (w['L'] - 1) / (8 * mp.pi ** 2)
                            numerical = mp.diff(potential, chi)
                            analytical = v[name + 'x']
                            check(label + ':' + name + '_potential_derivative', abs(numerical - analytical) / (abs(analytical) + mp.mpf('1e-200')), mp.mpf('1e-70'))
                        def contactpotential(xx):
                            w = point(xx, mg, rh, 1, phase, p, True, mp.mpf('1e12'), mp.mpf('1e14'))
                            return w['s'] * w['d'] * (w['L'] - 1) / (8 * mp.pi ** 2)
                        direct = mp.diff(contactpotential, chi)
                        check(label + ':contact_potential_derivative', abs(direct - v['contactx']) / abs(v['contactx']), mp.mpf('1e-70'))
                        directmass = mp.diff(lambda xx: point(xx, mg, rh, 1, phase, p, True, mp.mpf('1e12'), mp.mpf('1e14'))['mass2'], chi)
                        check(label + ':mass_derivative', abs(directmass - v['massx']) / (abs(directmass) + mp.mpf('1e-200')), mp.mpf('1e-70'))
                        if v['mass2'] > 0 and phase == 'quadrature':
                            direct = mp.diff(lambda xx: point(xx, mg, rh, 1, phase, p, True, mp.mpf('1e12'), mp.mpf('1e14'))['valley_v'], chi)
                            check(label + ':quadratic_valley_derivative', abs(direct - v['valley_x']) / (abs(direct) + mp.mpf('1e-200')), mp.mpf('1e-70'))
                        # Check stable cancellation identities with a component denominator.
                        residual = 2 * v['mass2'] + directmass
                        check(label + ':valley_residual_identity', abs(residual - v['residual']) / (abs(v['base']) + abs(v['restoration']) + abs(v['net_mass'])), mp.mpf('1e-70'))
                        check(label + ':stable_net_soft_identity', abs(v['D'] - v['c'] - v['dc']) / (abs(v['D']) + v['c']), mp.mpf('1e-70'))
                        key = (mt, rh, phase, index)
                        values = [v[k] for k in ('mass2', 'sourcex', 'persistx', 'contactx', 'valley_x')]
                        if digits == 160:
                            earlier[key] = values
                        else:
                            error = max(abs(a - b) / (abs(a) + abs(b) + mp.mpf('1e-200')) for a, b in zip(values, earlier[key]))
                            check(label + ':precision_agreement', error, mp.mpf('1e-65'))
                        rows.append(dict(digits=digits, mG_eV=mt, rH=rh, phase=phase, chi_index=index,
                                         **{k: mp.nstr(v[k], 45) for k in ('mass2', 'sourcex', 'persistx', 'contactx', 'valley_x', 'dc')}))
        # Asymptotic chi^-3 coefficient: no extrapolated cosmological clock.
        for rh in (0, 2, 4):
            chi, mg = mp.mpf('1e8'), mp.mpf('1e-6')
            v = point(chi, mg, rh, 1, 'zero', p, True, mp.mpf('1e12'), mp.mpf('1e14'))
            sinf = p['holomorphic_mass_i_eV'] ** 2 / (1 + p['xi'] / p['chi_i'] ** 2)
            Linf = mp.log(sinf / p['holomorphic_mass_i_eV'] ** 2)
            expected = -p['xi'] * sinf * v['dc'] * mg ** 2 * Linf / (4 * mp.pi ** 2)
            actual = -v['dc'] * v['kernelx'] * chi ** 3
            check(f'{digits}:asymptotic_tail_rH{rh}', abs(actual / expected - 1), mp.mpf('1e-10'))
    # Explicit dimensional exponents in eV, dimensionless chi and couplings.
    dimensions = {'A': 2 - 2, 'tau': 2 - 2, 'D': 0, 'persistent_mass2': 0 + 2,
                  'persistent_soft2': 0 + 2, 'CW_force4': 2 + 2,
                  'transverse_force4': 1 + 2 + 1, 'valley_force4': 2 + 4 + 2 - 4}
    expected_dimensions = {'A': 0, 'tau': 0, 'D': 0, 'persistent_mass2': 2,
                           'persistent_soft2': 2, 'CW_force4': 4, 'transverse_force4': 4, 'valley_force4': 4}
    for k, d in dimensions.items():
        check('dimensions:' + k, abs(d - expected_dimensions[k]), 0)
    p, chis = parameters()
    reps = []
    for mg in (1e-12, 1e-6, 1e-3):
        for rh in (0, 2, 4, 100):
            for sigma in (-1, 0, 1):
                for phase in PHASE:
                    reps.append(summarize(mg, rh, sigma, phase, p, chis))
    futures = []
    for mg in (1e-12, 1e-6, 1e-3):
        for rh in (0, 2, 4, 100):
            for phase in ('zero', 'quadrature'):
                v = future(mg, rh, 1, phase, p, chis[-1])
                futures.append(dict(mG_eV=mg, rH=rh, sigma=1, phase=phase, start_chi=chis[-1], stop_chi=chis[-1] + 80, **v))
                if v['cause'] == 'conservative_feedback':
                    check(f'future_root:{mg}:{rh}:{phase}', abs(v['root_residual']), 1e-9)
    csvwrite('independent_high_precision.csv', rows)
    csvwrite('independent_representatives.csv', reps)
    csvwrite('independent_future.csv', futures)
    hashes = {str(x.relative_to(REPO)): sha(x) for x in (Path(__file__), ROOT / 'protocol.md', PARENT, COORD)}
    result = dict(status='pass' if all(c['passed'] for c in checks) else 'fail',
                  checks_passed=sum(c['passed'] for c in checks), checks_total=len(checks), checks=checks,
                  high_precision_points=len(rows), representative_points=len(reps), future_points=len(futures),
                  scope='Declared leading jets only; local curvature and a quadratic valley, not full covariant stability or a UV completion.',
                  dimensions=dimensions, hashes=hashes, precision_digits=[160, 220],
                  python=platform.python_version(), mpmath=mp.__version__, numpy=np.__version__,
                  utc=datetime.now(timezone.utc).isoformat())
    dump('independent_summary.json', result)
    print('derive', result['checks_passed'], '/', result['checks_total'], flush=True)
    for c in checks:
        if not c['passed']:
            print('FAIL', c, flush=True)
    return result['status'] == 'pass'


def compare():
    p, chis = parameters()
    fields = {'tau': 'tau', 'A': 'A', 'source_c': 'c', 'persistent_soft_coefficient': 'D',
              'formal_net_c': 'formal_net_c', 'min_my2_eV2': 'min_mass2', 'min_transverse_schur_eV2': 'min_schur',
              'clockoff_my2_leading_eV2': 'asymptotic_mass2', 'max_radial_mass_eV': 'radial_mass',
              'source_CW_force_ratio': 'source_ratio', 'persistent_CW_force_ratio': 'persist_ratio',
              'clock_contact_CW_force_ratio': 'contact_ratio', 'contact_Btree_force_ratio': 'Btree_ratio',
              'contact_B5d_force_ratio': 'B5d_ratio', 'selected_signed_axis_ratio': 'signed',
              'selected_signed_chi_ratio': 'signed_chi', 'selected_axis_moduli_sum': 'axis_budget',
              'selected_feedback_sum_ratio': 'budget', 'quadratic_valley_force_ratio': 'valley_ratio',
              'max_valley_energy_over_U': 'valley_energy', 'valley_boundary_fraction': 'metric_fraction',
              'FT_over_MV2': 'FT_ratio', 'max_hierarchy_ratio': 'hierarchy',
              'max_coupling_loop_parameter': 'perturbativity', 'max_soft_over_charged_mass2': 'soft_fraction'}
    bools = {'finite_interval_local_pass': 'stable', 'clockoff_local_pass': 'asymptotic_positive',
             'stated_conditions_pass': 'necessary'}
    maxima, worst, differences, counts, hashes = {k: 0. for k in fields}, {}, [], {}, {}
    for name in ('scan.csv', 'representatives.csv'):
        path = OUT / name
        hashes[str(path.relative_to(REPO))] = sha(path)
        with path.open() as f:
            rows = list(csv.DictReader(f))
        counts[name] = len(rows)
        for index, row in enumerate(rows):
            mg, rh, sigma, mv, kk = [float(row[k]) for k in ('mG_eV', 'rH', 'sigma', 'M_vector_eV', 'mKK_eV')]
            v = summarize(mg, rh, sigma, row['phase'], p, chis, mv, kk)
            for mainkey, key in fields.items():
                expected = v[key]
                if expected is None:
                    if row[mainkey] != '':
                        differences.append(dict(file=name, row=index, field=mainkey, expected=None, got=row[mainkey]))
                    continue
                if row[mainkey] == '':
                    differences.append(dict(file=name, row=index, field=mainkey, expected=expected, got=''))
                    continue
                got = float(row[mainkey])
                scale = v['mass_scale'] if mainkey in ('min_my2_eV2', 'min_transverse_schur_eV2') else abs(expected) or 1.
                error = abs(got - expected) / scale
                if error > maxima[mainkey]:
                    maxima[mainkey] = error
                    worst[mainkey] = dict(file=name, row=index, expected=expected, got=got)
            for mainkey, key in bools.items():
                if (row[mainkey] == 'True') != v[key]:
                    differences.append(dict(file=name, row=index, field=mainkey, expected=v[key], got=row[mainkey]))
        print('compared', name, len(rows), flush=True)
    futurepath = OUT / 'future_limits.csv'
    hashes[str(futurepath.relative_to(REPO))] = sha(futurepath)
    future_max = 0.
    with futurepath.open() as f:
        rows = list(csv.DictReader(f))
    counts[futurepath.name] = len(rows)
    for index, row in enumerate(rows):
        v = future(float(row['mG_eV']), float(row['rH']), float(row['sigma']), row['phase'], p,
                   float(row['chi_start']), float(row['chi_search_end']))
        got = float(row['chi_first_stated_failure']) if row['chi_first_stated_failure'] else None
        if (v['first_chi'] is None) != (got is None):
            differences.append(dict(file=futurepath.name, row=index, field='chi_first_stated_failure', expected=v, got=got))
        elif got is not None:
            future_max = max(future_max, abs(got - v['first_chi']))
    tolerance = 1e-10
    checks = [dict(name=k, max_normalized_error=e, tolerance=tolerance, passed=bool(e <= tolerance)) for k, e in maxima.items()]
    checks.append(dict(name='first_future_failure_chi', maximum_absolute_error=future_max, tolerance=2e-9, passed=future_max <= 2e-9))
    checks.append(dict(name='classification_and_missing_value_agreement', difference_count=len(differences), passed=not differences))
    hashes.update({str(x.relative_to(REPO)): sha(x) for x in (Path(__file__), ROOT / 'protocol.md', PARENT, COORD)})
    result = dict(status='pass' if all(c['passed'] for c in checks) else 'fail', checks=checks,
                  checks_passed=sum(c['passed'] for c in checks), checks_total=len(checks), row_counts=counts,
                  maximum_normalized_error=max(maxima.values()), maximum_future_chi_error=future_max,
                  differences=differences[:100], worst=worst, hashes=hashes,
                  scope='Independent formulas; main outputs read solely for comparison, no main source read or imported.',
                  utc=datetime.now(timezone.utc).isoformat())
    dump('independent_comparison.json', result)
    print('compare', result['checks_passed'], '/', result['checks_total'], 'max', result['maximum_normalized_error'], flush=True)
    for c in checks:
        if not c['passed']:
            print('FAIL', c, flush=True)
    return result['status'] == 'pass'


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--derive', action='store_true')
    ap.add_argument('--compare', action='store_true')
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    okay = not args.derive or derive()
    if args.compare:
        okay = compare() and okay
    raise SystemExit(0 if okay else 1)
