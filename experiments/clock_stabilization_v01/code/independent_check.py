#!/usr/bin/env python3
"""Independent analytic/numerical audit of a declared local EFT.

No main-code import or execution.  --derive reads only protocols and historical
inputs; --compare additionally reads final main CSV data, never main code.
High precision checks use 180 and 240 decimal digits. This is an internal AI
reimplementation, not independent observational evidence or external review.
"""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import platform
from datetime import datetime, timezone

import mpmath as mp
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
REPO = ROOT.parents[1]
OUT = ROOT / 'results'
OLD = REPO / 'experiments/sugra_shift_v01/results'
PARENT = REPO / 'experiments/sequestering_test_v01/results/inputs.json'
PHASE = {'zero': (1, 0), 'quadrature': (0, 1), 'pi': (-1, 0)}


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def dump(name, obj):
    (OUT / name).write_text(json.dumps(obj, indent=2, ensure_ascii=False) + '\n')


def write_csv(name, rows):
    with (OUT / name).open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator='\n')
        w.writeheader()
        w.writerows(rows)


def load():
    raw = json.loads(PARENT.read_text())['parent_parameters']
    with (OLD / 'trajectories.csv').open() as f:
        chis = [r['chi'] for r in csv.DictReader(f) if r['case'] == 'shift_axis']
    return raw, chis


def parameters(high_precision=False):
    raw, chis = load()
    cast = (lambda x: mp.mpf(str(x))) if high_precision else float
    p = {k: cast(v) for k, v in raw.items() if isinstance(v, (int, float))}
    p['P'] = p['Mpl_eV'] ** 2
    p['F'] = (mp.sqrt if high_precision else np.sqrt)(p['F_squared_eV2'])
    p['Ui'] = p['potential_unit_eV4'] * p['W_i']
    p['rho'] = p['rho_lambda_eV4']
    p['C5'] = mp.zeta(3) / (48 * mp.pi ** 4) if high_precision else float(mp.zeta(3) / (48 * mp.pi ** 4))
    return p, [cast(c) for c in chis]


def point(chi, mg, mv, kk, sigma, phase, p, high_precision=False):
    """Independent real-axis values/jets obtained from the CW polynomial."""
    b = mp if high_precision else np
    pi, sqrt, exp, log = b.pi, b.sqrt, b.exp, b.log
    co, si = PHASE[phase]
    P, F, F2 = p['P'], p['F'], p['F_squared_eV2']
    U = p['Ui'] * exp(-2 * (chi - p['chi_i']))
    s = p['holomorphic_mass_i_eV'] ** 2 * (1 + p['xi'] / chi ** 2) / (1 + p['xi'] / p['chi_i'] ** 2)
    M = sqrt(s)
    r = -p['xi'] / (chi * (chi ** 2 + p['xi']))
    schi = 2 * r * s
    L = log(s / p['holomorphic_mass_i_eV'] ** 2)
    A = p['C5'] * kk ** 2 / P
    c = (12 * A + 3 * A ** 2) / (1 + A) ** 2
    beta = -6 * A / (1 + A)
    wm = F * sqrt(U / 2)
    rw = wm / (mg * P)
    ft_re, ft_im = mg * (co - rw), -mg * si
    X = ft_re ** 2 + ft_im ** 2
    Xchi = 2 * ft_re * mg * rw
    Xy = 2 * mg * wm * si / (3 * P)
    pp, ppchi, ppy = s * X, schi * X + s * Xchi, s * Xy
    cgchi = -c * (ppchi * (L - 1) + pp * 2 * r) / (8 * pi ** 2)
    cgy = -c * ppy * (L - 1) / (8 * pi ** 2)
    d = 2 * sigma * U / mv ** 2
    dchi = -2 * d
    dy = -12 * sqrt(2) * sigma * A * F * sqrt(U) * mg * si / mv ** 2
    dgchi = (schi * d * L + s * dchi * (L - 1)) / (8 * pi ** 2)
    dgy = s * dy * (L - 1) / (8 * pi ** 2)
    eta = p['rho'] / (6 * P * mg ** 2)
    hc = sqrt(2 * U) / F
    B0re, B0im = M * (hc * r + 2 * eta * ft_re), M * 2 * eta * ft_im
    dBy_im = 4 * sqrt(2) * sigma * M * F * sqrt(U) / mv ** 2
    bg_tree_y = B0im * dBy_im * L / (8 * pi ** 2)
    bg_5d_y = (-beta * M * ft_im) * dBy_im * L / (8 * pi ** 2)
    rr = p['rho'] / (3 * P * mg ** 2)
    baseline1 = 4 * (1 + rr) * U / (3 * P)
    baseline2 = -2 * rr * mg * wm * co / F2
    baseline3 = 2 * rr * mg ** 2 * (2 + rr) * ((co - rw) ** 2 + si ** 2)
    baseline = baseline1 + baseline2 + baseline3
    stabilization = 12 * U / mv ** 2
    correction = -24 * A * mg ** 2
    mass2 = baseline + stabilization + correction
    masschi = (-2 * baseline1 - baseline2 + 2 * rr * mg ** 2 * (2 + rr)
               * (2 * co * rw - 2 * rw ** 2) - 2 * stabilization)
    Jy = -6 * sqrt(2) * A * F * sqrt(U) * mg * si
    # The following are local quadratic estimates only, never a solved orbit.
    with np.errstate(divide='ignore', invalid='ignore'):
        yy = -Jy / (F2 * mass2)
        vv = -Jy ** 2 / (2 * F2 * mass2)
        # Expand 2m²+(m²)' first: its U/Λ² terms cancel identically.
        residual = baseline2 + 4 * rr * mg ** 2 * (2 + rr) * (1 - co * rw) + 2 * correction
        vgchi = 36 * A ** 2 * mg ** 2 * si ** 2 * U * residual / mass2 ** 2
    return dict(U=U, s=s, M=M, r=r, L=L, A=A, c=c, beta=beta,
                d=d, dchi=dchi, dy=dy, deltaBy_im=dBy_im,
                B0re=B0re, B0im=B0im, cgchi=cgchi, cgy=cgy,
                dgchi=dgchi, dgy=dgy, bg_tree_y=bg_tree_y, bg_5d_y=bg_5d_y,
                baseline=baseline, stabilization=stabilization,
                correction=correction, mass2=mass2, masschi=masschi,
                Jy=Jy, valley_y=yy, valley_v=vv, valley_gchi=vgchi,
                metric_fraction=sqrt(6) * F * abs(yy) / mv)


def summarize(mg, mv, kk, sigma, phase, p, chi):
    x = point(np.asarray(chi), mg, mv, kk, sigma, phase, p)
    cforce = np.max(np.hypot(x['cgchi'], x['cgy']) / (2 * x['U']))
    contact = np.max(np.hypot(x['dgchi'], x['dgy']) / (2 * x['U']))
    treeB = np.max(np.abs(x['bg_tree_y']) / (2 * x['U']))
    loopB = np.max(np.abs(x['bg_5d_y']) / (2 * x['U']))
    stable = bool(np.all(x['mass2'] > 0))
    valley = float(np.max(np.abs(x['valley_gchi']) / (2 * x['U']))) if stable else None
    metric = float(np.max(x['metric_fraction'])) if stable else None
    pointforce = (np.hypot(x['cgchi'], x['cgy']) + np.hypot(x['dgchi'], x['dgy'])
                  + np.abs(x['bg_tree_y']) + np.abs(x['bg_5d_y'])) / (2 * x['U'])
    axisforce = float(np.max(pointforce))
    sumforce = float(np.max(pointforce + np.abs(x['valley_gchi']) / (2 * x['U']))) if stable else axisforce
    ft = np.sqrt(3) * p['Mpl_eV'] * mg / mv ** 2
    fz = np.sqrt(np.max(x['U'])) / mv ** 2
    M5 = (p['P'] * kk / np.pi) ** (1 / 3)
    mass_hierarchy = max(p['holomorphic_mass_i_eV'] / mv, mv / kk, kk / M5,
                         np.sqrt(48) * mg / mv, mg / mv,
                         np.sqrt(np.max(np.abs(x['mass2']))) / mv, p['H_ref_eV'] / mv)
    perturb = 2 * max(1, sigma ** 2) / (16 * np.pi ** 2)
    criticalU = (24 * x['A'] * mg ** 2 - 4 * p['rho'] / (3 * p['P'])) / (12 / mv ** 2 + 4 / (3 * p['P']))
    criticalchi = float(p['chi_i'] + np.log(p['Ui'] / criticalU) / 2) if criticalU > 0 else None
    return dict(mg=mg, mv=mv, kk=kk, sigma=sigma, phase=phase, A=float(x['A']),
                cforce=float(cforce), contactforce=float(contact), treeBforce=float(treeB),
                loopBforce=float(loopB), valleyforce=valley, force_sum=float(sumforce), axisforce=axisforce,
                valley_energy=float(np.max(np.abs(x['valley_v']) / x['U'])) if stable else None,
                min_mass2=float(np.min(x['mass2'])), max_mass2=float(np.max(x['mass2'])),
                min_baseline=float(np.min(x['baseline'])), max_stabilization=float(np.max(x['stabilization'])),
                max_abs_Jy=float(np.max(np.abs(x['Jy']))),
                metric_fraction=metric, FT_ratio=float(ft), FZ_ratio=float(fz),
                mass_hierarchy=float(mass_hierarchy), perturbativity=float(perturb),
                stable=stable, force_ok=bool(sumforce <= .1),
                auxiliary_ok=bool(max(ft, fz) <= .1), hierarchy_ok=bool(mass_hierarchy <= .1),
                metric_ok=bool(stable and metric <= .1), perturbative_ok=bool(perturb <= .1),
                necessary_intersection=bool(stable and sumforce <= .1 and max(ft, fz) <= .1
                     and mass_hierarchy <= .1 and metric <= .1 and perturb <= .1),
                asymptotic_mass2=float(4 * p['rho'] / (3 * p['P']) - 24 * x['A'] * mg ** 2),
                critical_U=float(criticalU), critical_chi=criticalchi)


def cw(s, d, bre, bim):
    """Two complex scalars and two Weyl fermions, unexpanded determinant."""
    bb = mp.sqrt(bre ** 2 + bim ** 2)
    def f(x):
        return x ** 2 * (mp.log(x / mp.mpf('1e22')) - mp.mpf('1.5'))
    return (f(s + d + bb) + f(s + d - bb) - 2 * f(s)) / (32 * mp.pi ** 2)


def direct_potential_mass(chi, mg, mv, phase, p, quartic):
    """Dimensionless complete 2x2 no-scale inverse at fixed finite-rho jet.

    This independent check sets A=0; the A-dependent correction is audited in
    the separate geometry report. It does not integrate out T before Hessian.
    """
    eps = p['F_squared_eV2'] / p['P']
    rr = p['rho'] / (3 * p['P'] * mg ** 2)
    q = rr / (1 + rr)
    co, si = PHASE[phase]
    U = p['Ui'] * mp.exp(-2 * (chi - p['chi_i']))
    rw = mp.sqrt(p['F_squared_eV2'] * U / 2) / (mg * p['P'])
    def potential(y):
        k = eps * y ** 2
        lam = mv ** 2 / p['P']
        kcorr = k / lam if quartic else mp.mpf(0)
        Om = 1 - k * (1 - kcorr) / 3
        oi = mp.matrix([1, 1j * mp.sqrt(2 * eps) * y * (1 - 2 * kcorr) / 3])
        oij = mp.matrix([[q, 0], [0, -(1 - 6 * kcorr) / 3]])
        metric = -3 * (oij / Om - oi * oi.transpose_conj() / Om ** 2)
        kval = -3 * oi / Om
        wclock = -rw * mp.exp(-1j * y)
        w = mp.mpc(co, si) + wclock
        dw = kval * w + mp.matrix([0, -mp.sqrt(2 / eps) * wclock])
        return mp.re(((dw.transpose_conj() * metric ** -1 * dw)[0] - 3 * abs(w) ** 2) / Om ** 3)
    return mg ** 2 * mp.diff(potential, 0, 2) / eps


def derive():
    checks, rows, precision_results = [], [], {}
    def check(label, value, limit):
        checks.append(dict(name=label, error=mp.nstr(value, 30), tolerance=str(limit), passed=bool(value <= limit)))
    for digits in (180, 240):
        mp.mp.dps = digits
        p, chis = parameters(True)
        results = []
        # Each case has independent direct determinant derivatives and a 2x2
        # metric inverse, not an algebraic copy of a main output.
        for mgtext in ('1e-12', '1e-6', '1e-5', '1'):
            mg = mp.mpf(mgtext)
            for phase in PHASE:
                for index in (0, 512, 1024):
                    chi, mv, kk = chis[index], mp.mpf('1e12'), mp.mpf('1e14')
                    r = point(chi, mg, mv, kk, mp.mpf(1), phase, p, True)
                    label = f'dps{digits}:{mgtext}:{phase}:{index}'
                    # Direct no-scale finite-rho Hessian, with/without quartic.
                    m0 = direct_potential_mass(chi, mg, mv, phase, p, False)
                    m1 = direct_potential_mass(chi, mg, mv, phase, p, True)
                    err0 = abs(m0 - r['baseline']) / (abs(m0) + abs(r['baseline']))
                    err1 = abs((m1 - m0) - r['stabilization']) / r['stabilization']
                    check(label + ':finite_rho_mass', err0, mp.mpf('1e-70'))
                    check(label + ':quartic_mass', err1, mp.mpf('1e-70'))
                    # Direct determinant responds to the new contact and B jet.
                    # d0=0 isolates this selected operator; its effect on the
                    # expansion is O(d0/s), below every stated comparison scale.
                    def delta_potential_x(xx):
                        v = point(xx, mg, mv, kk, mp.mpf(1), phase, p, True)
                        return cw(v['s'], v['d'], v['B0re'], v['B0im']) - cw(v['s'], 0, v['B0re'], v['B0im'])
                    gradx = mp.diff(delta_potential_x, chi)
                    def delta_potential_y(y):
                        return (cw(r['s'], r['d'] + r['dy'] * y,
                                   r['B0re'], r['B0im'] + r['deltaBy_im'] * y)
                                - cw(r['s'], 0, r['B0re'], r['B0im']))
                    grady = mp.diff(delta_potential_y, 0)
                    sc = abs(r['dgchi']) + abs(r['dgy']) + abs(r['bg_tree_y'])
                    errx = abs(gradx - r['dgchi']) / sc
                    erry = abs(grady - r['dgy'] - r['bg_tree_y']) / sc
                    check(label + ':direct_CW_x', errx, mp.mpf('1e-35'))
                    check(label + ':direct_CW_y', erry, mp.mpf('1e-35'))
                    co, si = PHASE[phase]
                    rw = p['F'] * mp.sqrt(r['U'] / 2) / (mg * p['P'])
                    b5re = -r['beta'] * r['M'] * mg * (co - rw)
                    b5im = r['beta'] * r['M'] * mg * si
                    def delta_potential_y5(y):
                        return (cw(r['s'], r['d'] + r['dy'] * y, r['B0re'] + b5re,
                                   r['B0im'] + b5im + r['deltaBy_im'] * y)
                                - cw(r['s'], 0, r['B0re'] + b5re, r['B0im'] + b5im))
                    grady5 = mp.diff(delta_potential_y5, 0)
                    sc5 = sc + abs(r['bg_5d_y'])
                    check(label + ':direct_CW_y_with_5d_B', abs(grady5 - r['dgy'] - r['bg_tree_y'] - r['bg_5d_y']) / sc5, mp.mpf('1e-35'))
                    # Exact derivative of the quadratic valley formula verifies
                    # its force without asserting that this is a solved valley.
                    if r['mass2'] > 0:
                        vx = mp.diff(lambda xx: point(xx, mg, mv, kk, 1, phase, p, True)['valley_v'], chi)
                        vscale = abs(r['valley_gchi']) + mp.mpf('1e-200')
                        check(label + ':quadratic_valley_derivative', abs(vx - r['valley_gchi']) / vscale, mp.mpf('1e-60'))
                    row = dict(digits=digits, mg_eV=mgtext, MV_eV='1e12', KK_eV='1e14',
                               phase=phase, chi_index=index, chi=mp.nstr(chi, 30))
                    for key in ('U', 'd', 'dchi', 'dy', 'deltaBy_im', 'dgchi', 'dgy', 'bg_tree_y',
                                'bg_5d_y', 'baseline', 'stabilization', 'mass2', 'valley_y', 'valley_gchi'):
                        row[key] = mp.nstr(r[key], 90)
                    row.update(direct_mass0=mp.nstr(m0, 90), direct_mass1=mp.nstr(m1, 90),
                               direct_CW_x=mp.nstr(gradx, 90), direct_CW_y=mp.nstr(grady, 90),
                               direct_CW_y_with_5d_B=mp.nstr(grady5, 90))
                    rows.append(row)
                    results.append((m0, m1, gradx, grady))
        precision_results[digits] = results
        print('precision', digits, 'points', len(results), flush=True)
    # Compare each quantity in its own units; do not add masses to gradients.
    # Exact phase zeros remain exact zeros and are not divided by a tiny signal.
    for i, (lo, hi) in enumerate(zip(precision_results[180], precision_results[240])):
        errors = [abs(a - b) / (abs(b) if b else 1) for a, b in zip(lo, hi)]
        check(f'precision_repeat:{i}', max(errors), mp.mpf('1e-70'))
    p, chis = parameters(True)
    rr = point(chis[-1], mp.mpf('1e-6'), mp.mpf('1e12'), mp.mpf('1e14'), 1, 'quadrature', p, True)
    # Shared-current matching: -J^2/(2MV²) produces k² and kSQ with one ratio.
    for sigma in (-3, -1, 0, 1, 3):
        r = point(chis[-1], mp.mpf('1e-6'), mp.mpf('1e12'), mp.mpf('1e14'), sigma, 'quadrature', p, True)
        for key in ('d', 'dy', 'deltaBy_im', 'dgchi', 'dgy', 'bg_tree_y', 'bg_5d_y'):
            check(f'sigma_linearity:{sigma}:{key}', abs(r[key] - sigma * rr[key]) / (abs(rr[key]) + mp.mpf('1e-200')), mp.mpf('1e-100'))
        check(f'sigma_mass_independent:{sigma}', abs(r['mass2'] - rr['mass2']), mp.mpf('1e-100'))
    summary = dict(status='pass' if all(c['passed'] for c in checks) else 'fail',
                   checks_passed=sum(c['passed'] for c in checks), checks_total=len(checks),
                   high_precision_points_per_precision=36, precisions=[180, 240],
                   chi_coordinate_count=len(chis), checks=checks,
                   scope='Independent local analytic numeric audit; not full UV spectrum, global stability or observation.',
                   independence='No main implementation read/import/execute; shared protocol and historic input; same AI system.',
                   hashes={str(p.relative_to(REPO)): sha(p) for p in
                           (Path(__file__), ROOT / 'protocol.md', PARENT, OLD / 'trajectories.csv')},
                   runtime=dict(python=platform.python_version(), numpy=np.__version__, mpmath=mp.__version__),
                   utc=datetime.now(timezone.utc).isoformat())
    write_csv('independent_high_precision.csv', rows)
    dump('independent_summary.json', summary)
    p, chi = parameters()
    representatives = []
    for mg in (1e-12, 1e-6, 1e-5, 1):
        for mv in (1e11, 1e12, 1e13):
            for sigma in (-1, 0, 1):
                for phase in PHASE:
                    representatives.append(summarize(mg, mv, 1e14, sigma, phase, p, chi))
    write_csv('independent_representatives.csv', representatives)
    print('checks', summary['checks_passed'], '/', summary['checks_total'], flush=True)
    return summary['status'] == 'pass'


def compare():
    p, chi = parameters()
    fields = {
        'A_5d': 'A', 'min_my2_eV2': 'min_mass2', 'max_my2_eV2': 'max_mass2',
        'clockoff_my2_leading_eV2': 'asymptotic_mass2',
        'RSS_c_2loop_force_ratio': 'cforce', 'contact_d_1loop_force_ratio': 'contactforce',
        'contact_Btree_1loop_force_ratio': 'treeBforce', 'contact_B5d_2loop_force_ratio': 'loopBforce',
        'quadratic_valley_force_ratio': 'valleyforce', 'quadratic_valley_energy_over_U': 'valley_energy',
        'valley_metric_boundary_fraction': 'metric_fraction', 'selected_feedback_sum_ratio': 'force_sum',
        'known_axis_feedback_without_valley': 'axisforce', 'max_hierarchy_ratio': 'mass_hierarchy',
        'FTcanonical_over_vector2': 'FT_ratio', 'max_coupling_loop_parameter': 'perturbativity',
        'Ucrit_leading_eV4': 'critical_U', 'chi_critical_leading': 'critical_chi'}
    boolfields = {'local_positive_curvature': 'stable', 'hierarchy_pass': 'hierarchy_ok',
                  'auxiliary_pass': 'auxiliary_ok', 'perturbative_coupling_pass': 'perturbative_ok',
                  'necessary_conditions_pass': 'necessary_intersection'}
    maxima = {k: 0. for k in fields}
    worst = {}
    mismatches, counts, hashes = [], {}, {}
    for name in ('scan.csv', 'representatives.csv', 'sigma_sensitivity.csv'):
        path = OUT / name
        hashes[str(path.relative_to(REPO))] = sha(path)
        with path.open() as f:
            source = list(csv.DictReader(f))
        counts[name] = len(source)
        for index, row in enumerate(source):
            mg, mv, kk, sigma = [float(row[k]) for k in ('mG_eV', 'M_vector_eV', 'mKK_eV', 'sigma')]
            v = summarize(mg, mv, kk, sigma, row['phase'], p, chi)
            for mainkey, key in fields.items():
                expected = v[key]
                if mainkey == 'selected_feedback_sum_ratio' and not v['stable']:
                    expected = None
                if expected is None:
                    if row[mainkey] != '':
                        mismatches.append(dict(file=name, row=index, field=mainkey, expected=None, got=row[mainkey]))
                    continue
                if row[mainkey] == '':
                    mismatches.append(dict(file=name, row=index, field=mainkey, expected=expected, got=''))
                    continue
                got = float(row[mainkey])
                if mainkey in ('min_my2_eV2', 'max_my2_eV2'):
                    # Use component scales at cancellations; never hide them
                    # behind a denominator near the zero of the total mass.
                    scale = (v['max_stabilization'] + abs(24 * v['A'] * mg ** 2) + abs(v['min_baseline']))
                else:
                    scale = abs(expected) if expected != 0 else 1.
                error = abs(got - expected) / scale
                if error > maxima[mainkey]:
                    maxima[mainkey] = error
                    worst[mainkey] = dict(file=name, row=index, expected=expected, got=got)
            for mainkey, key in boolfields.items():
                if (row[mainkey] == 'True') != v[key]:
                    mismatches.append(dict(file=name, row=index, field=mainkey, expected=v[key], got=row[mainkey]))
        print('compared', name, len(source), flush=True)
    front = OUT / 'conditional_frontiers.csv'
    hashes[str(front.relative_to(REPO))] = sha(front)
    with front.open() as f:
        source = list(csv.DictReader(f))
    counts[front.name] = len(source)
    front_mass, front_aux = 0., 0.
    Umin = p['Ui'] * np.exp(-2 * (chi[-1] - p['chi_i']))
    for index, row in enumerate(source):
        kk, mv = float(row['mKK_eV']), float(row['M_vector_eV'])
        A = p['C5'] * kk ** 2 / p['P']
        mcurv = np.sqrt((12 * Umin / mv ** 2 + 4 * (Umin + p['rho']) / (3 * p['P'])) / (24 * A))
        maux = .1 * mv ** 2 / (np.sqrt(3) * p['Mpl_eV'])
        front_mass = max(front_mass, abs(float(row['local_curvature_mG_max_eV']) / mcurv - 1))
        front_aux = max(front_aux, abs(float(row['auxiliary_mG_max_eV']) / maux - 1))
        gaps = p['holomorphic_mass_i_eV'] / mv <= .1 and mv / kk <= .1
        if (row['two_mass_gaps_pass'] == 'True') != gaps:
            mismatches.append(dict(file=front.name, row=index, field='two_mass_gaps_pass', expected=gaps))
    maxima.update(frontier_curvature=front_mass, frontier_auxiliary=front_aux)
    tolerance = 1e-11
    checks = [dict(name=key, max_normalized_error=error, tolerance=tolerance,
                   passed=bool(error <= tolerance)) for key, error in maxima.items()]
    checks.append(dict(name='classification_and_missing_value_agreement', differences=len(mismatches), passed=not mismatches))
    hashes.update({str(path.relative_to(REPO)): sha(path) for path in (Path(__file__), ROOT / 'protocol.md', PARENT, OLD / 'trajectories.csv')})
    result = dict(status='pass' if all(c['passed'] for c in checks) else 'fail',
                  checks_passed=sum(c['passed'] for c in checks), checks_total=len(checks), checks=checks,
                  row_counts=counts, worst=worst, differences=mismatches[:50], hashes=hashes,
                  maximum_normalized_error=max(maxima.values()),
                  scope='Internal independent local implementation; no main code read, import or execution.',
                  utc=datetime.now(timezone.utc).isoformat())
    dump('independent_comparison.json', result)
    print('comparison', result['checks_passed'], '/', result['checks_total'], 'max', result['maximum_normalized_error'], flush=True)
    return result['status'] == 'pass'


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--derive', action='store_true')
    ap.add_argument('--compare', action='store_true')
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    okay = True
    if args.derive:
        okay = derive()
    if args.compare:
        okay = compare() and okay
    raise SystemExit(0 if okay else 1)
