#!/usr/bin/env python3
"""Independent selected-operator verification; never imports the production code.

Derives kernels from the two charged scalar eigenvalues and one Dirac pair.
Only the specified local Q^2 sector is checked, not complete SUGRA loops.
"""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import platform
import time

import mpmath as mp
import numpy as np

HERE = Path(__file__).resolve().parents[1]
REPO = HERE.parents[1]
OUT = HERE / 'results'
INPUT = REPO / 'experiments/sugra_shift_v01/results/inputs.json'
TRAJECTORY = REPO / 'experiments/sugra_shift_v01/results/trajectories.csv'
RAW = json.loads(INPUT.read_text())
with TRAJECTORY.open() as stream:
    COORDS = [r for r in csv.DictReader(stream) if r['case'] == 'shift_axis']
PHASES = {'zero': (1, 0), 'quadrature': (0, 1), 'pi': (-1, 0)}
MASSES = ('1e-6', '1', '1e5', '1e7')
SCALES = ('0.5', '1', '2', '10')
INDICES = (0, 173, 701, 1024)
CHECKS = []


def m(x):
    return mp.mpf(str(x))


def params():
    return dict(P=m(RAW['Mpl_eV']) ** 2, F=mp.sqrt(m(RAW['F_squared_eV2'])),
                Ui=m(RAW['potential_unit_eV4']) * m(RAW['W_i']),
                ci=m(RAW['chi_i']), xi=m(RAW['xi']), Mi=m('1e11'))


def record(name, residual, tolerance):
    CHECKS.append(dict(name=name, residual=str(abs(residual)), tolerance=str(tolerance),
                       passed=bool(abs(residual) <= m(tolerance))))


def norm(vec):
    return mp.sqrt(sum(v * v for v in vec))


def local_data(chi, mg, phase):
    p = params()
    co, si = PHASES[phase]
    U = p['Ui'] * mp.exp(-2 * (chi - p['ci']))
    r = -p['F'] * mp.sqrt(U / 2) / p['P']
    X = mg ** 2 + 2 * mg * co * r + r ** 2
    Xd = (-2 * r * (mg * co + r), -2 * mg * r * si / 3)
    den = 1 + p['xi'] / p['ci'] ** 2
    s = p['Mi'] ** 2 * (1 + p['xi'] / chi ** 2) / den
    sd = (-2 * p['Mi'] ** 2 * p['xi'] / (den * chi ** 3), mp.mpf(0))
    product = s * X
    pd = tuple(sd[j] * X + s * Xd[j] for j in range(2))
    return dict(U=U, r=r, X=X, Xd=Xd, s=s, sd=sd, p=product, pd=pd)


def kernels(d, scale=mp.mpf(1)):
    L = mp.log(d['s'] / (params()['Mi'] * scale) ** 2)
    # Factor analytically before evaluation: at L=0 the two large sd*X terms
    # cancel, while the quadrature X_chi term can be more than 80 orders smaller.
    cc = tuple((-d['sd'][j] * d['X'] * L + d['s'] * d['Xd'][j] * (1 - L)) /
               (8 * mp.pi ** 2) for j in range(2))
    bb = tuple((d['sd'][j] * d['X'] * (L + 1) + d['s'] * d['Xd'][j] * L) /
               (16 * mp.pi ** 2) for j in range(2))
    nn = tuple(v / (8 * mp.pi ** 2) for v in d['pd'])
    return cc, bb, nn


def series_force(d, c, b2, nu, scale):
    kk = kernels(d, scale)
    return tuple(c * kk[0][j] + b2 * kk[1][j] + nu * kk[2][j] for j in range(2))


def direct_spectrum(d, c, b2, nu, scale):
    """Direct logarithms of physical eigenvalues; no small-split expansion."""
    s, X = d['s'], d['X']
    delta = -c * X
    h = mp.sqrt(b2 * s * X)
    sp, sm = s + delta + h, s + delta - h
    if min(sp, sm) <= 0:
        raise ValueError('Tachyon in independent direct spectrum')
    mu2 = (params()['Mi'] * scale) ** 2

    def f(x):
        return x ** 2 * (mp.log(x / mu2) - m('1.5'))

    def df(x):
        return 2 * x * (mp.log(x / mu2) - 1)

    potential = (f(sp) + f(sm) - 2 * f(s)) / (32 * mp.pi ** 2)
    potential += nu * s * X / (8 * mp.pi ** 2)
    derivatives = []
    for j in range(2):
        ds = d['sd'][j]
        dd = -c * d['Xd'][j]
        dh = b2 * d['pd'][j] / (2 * h) if h else mp.mpf(0)
        value = (df(sp) * (ds + dd + dh) + df(sm) * (ds + dd - dh)
                 - 2 * df(s) * ds) / (32 * mp.pi ** 2)
        derivatives.append(value + nu * d['pd'][j] / (8 * mp.pi ** 2))
    return potential, tuple(derivatives), abs(delta / s), abs(h / s)


def write_csv(path, rows):
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def string(x):
    return mp.nstr(x, 50)


def derive():
    started = time.time()
    OUT.mkdir(exist_ok=True)
    mp.mp.dps = 80
    ip = params()
    record('1025 coordinates', len(COORDS) - 1025, 0)
    record('first chi equals fixed initial chi', m(COORDS[0]['chi']) - ip['ci'], 0)
    record('exact initial charged mass', local_data(ip['ci'], m(1), 'zero')['s'] /
           ip['Mi'] ** 2 - 1, '1e-75')
    basis_rows, budget_rows, all_budgets = [], [], {}
    max_null = mp.mpf(0)
    max_rg = mp.mpf(0)
    max_yjet = mp.mpf(0)
    singular_rows = []

    for mg_text in MASSES:
        mg = m(mg_text)
        for phase, (co, si) in PHASES.items():
            maxima = [[mp.mpf(0)] * 3 for _ in range(3)]
            indices = [[0] * 3 for _ in range(3)]
            matrix = []
            for idx, row in enumerate(COORDS):
                chi = m(row['chi'])
                d = local_data(chi, mg, phase)
                kk = kernels(d)
                budget = 2 * d['U']
                scale = max(norm(k) for k in kk)
                err = norm(tuple(kk[0][j] + 2 * kk[1][j] - kk[2][j]
                                 for j in range(2))) / scale
                max_null = max(max_null, err)
                out = dict(mG_eV=mg_text, phase=phase, index=idx, chi=string(chi),
                           U_eV4=string(d['U']))
                for n, name in enumerate(('c', 'beta2', 'nu')):
                    ratios = (abs(kk[n][0]) / budget, abs(kk[n][1]) / budget,
                              norm(kk[n]) / budget)
                    for j, direction in enumerate(('chi', 'y', 'norm')):
                        out[f'K_{name}_{direction}_over_2U'] = string(ratios[j])
                        if ratios[j] > maxima[n][j]:
                            maxima[n][j], indices[n][j] = ratios[j], idx
                    out[f'K_{name}_chi_signed'] = string(kk[n][0])
                    out[f'K_{name}_y_signed'] = string(kk[n][1])
                basis_rows.append(out)
                matrix.extend([[float(kk[n][j] / budget) for n in range(3)]
                               for j in range(2)])
            # Column normalization is stated explicitly; no singular value has units.
            matrix = np.asarray(matrix)
            col = np.linalg.norm(matrix, axis=0)
            sv = np.linalg.svd(matrix / col, compute_uv=False)
            singular_rows.append(dict(mG_eV=mg_text, phase=phase,
                                      column_l2_normalized_singular_values=sv.tolist(),
                                      smallest_over_largest=float(sv[-1] / sv[0])))
            record(f'normalized response rank two {mg_text} {phase}', sv[-1] / sv[0], '1e-13')
            limits = [m('.1') / x[2] for x in maxima]
            all_budgets[(mg_text, phase)] = limits
            for n, name in enumerate(('c', 'beta2', 'nu')):
                budget_rows.append(dict(mG_eV=mg_text, phase=phase, coefficient=name,
                                        limit_10pct=string(limits[n]),
                                        beta_limit_10pct=string(mp.sqrt(limits[n])) if n == 1 else '',
                                        max_chi=string(maxima[n][0]), max_y=string(maxima[n][1]),
                                        max_norm=string(maxima[n][2]), index_chi=indices[n][0],
                                        index_y=indices[n][1], index_norm=indices[n][2]))
            # Differentiate the full auxiliary combination to verify the y factor 1/3.
            for idx in INDICES:
                chi = m(COORDS[idx]['chi'])
                d = local_data(chi, mg, phase)

                def full_x(y):
                    S = 3 * mg * mp.mpc(co, si) + d['r'] * (3 + 2j * y) * mp.exp(-1j * y)
                    return abs(S) ** 2 / 9

                # This full norm is dominated by mG^2. The original 80-digit
                # differentiation lost about 33 digits; its failed audit is saved.
                with mp.workdps(160):
                    recovered = mp.diff(full_x, mp.mpf(0))
                yscale = max(abs(mg * d['r']), abs(d['r'] ** 2))
                max_yjet = max(max_yjet, abs(recovered - d['Xd'][1]) / yscale)
                # RG closure uses arbitrary independent finite coefficients, not a fit.
                cc, bb, nn = m('.37') * limits[0], m('.43') * limits[1], -m('.19') * limits[2]
                baseline = series_force(d, cc, bb, nn, mp.mpf(1))
                for scale_text in SCALES:
                    mu = m(scale_text)
                    running_nu = nn + 2 * (-cc + bb / 2) * mp.log(mu)
                    shifted = series_force(d, cc, bb, running_nu, mu)
                    err = norm(tuple(shifted[j] - baseline[j] for j in range(2))) / (2 * d['U'])
                    max_rg = max(max_rg, err)
    record('kernel exact null identity all grid points', max_null, '1e-75')
    record('RG running closure independent coefficients', max_rg, '1e-75')
    record('transverse full S derivative factor one third', max_yjet, '1e-70')
    write_csv(OUT / 'independent_kernels.csv', basis_rows)
    write_csv(OUT / 'independent_budgets.csv', budget_rows)

    # Direct spectra use independent sampling positions, all requested matching scales,
    # and both signs of the common split. Coefficients are frozen at reference scale.
    exact_rows = []
    maxima = dict(precision_force_over_2U=mp.mpf(0), neglected_force_over_2U=mp.mpf(0),
                  common_split_over_s=mp.mpf(0), B_over_s=mp.mpf(0))
    for mg_text in MASSES:
        for phase in PHASES:
            cc0, bb0, nn0 = all_budgets[(mg_text, phase)]
            scenarios = {
                'contact_positive': (m('.73') * cc0, mp.mpf(0), mp.mpf(0)),
                'contact_negative': (-m('.61') * cc0, mp.mpf(0), mp.mpf(0)),
                'B_squared': (mp.mpf(0), m('.67') * bb0, mp.mpf(0)),
                'local_only': (mp.mpf(0), mp.mpf(0), m('.41') * nn0),
                'mixed': (m('.31') * cc0, m('.29') * bb0, -m('.23') * nn0),
                'contact_finite_shift': (m('.73') * cc0, mp.mpf(0), -m('.73') * cc0),
                'B_finite_shift': (mp.mpf(0), m('.67') * bb0, -m('.67') * bb0 / 2),
            }
            for idx in INDICES:
                for scale_text in SCALES:
                    for name, coefficients in scenarios.items():
                        results = {}
                        for precision in (160, 220):
                            with mp.workdps(precision):
                                d = local_data(m(COORDS[idx]['chi']), m(mg_text), phase)
                                scale = m(scale_text)
                                cc, bb, nn = coefficients
                                nu = nn + 2 * (-cc + bb / 2) * mp.log(scale)
                                value, force, split, h = direct_spectrum(d, cc, bb, nu, scale)
                                approximation = series_force(d, cc, bb, nu, scale)
                                results[precision] = (value, force, approximation, split, h, d['U'])
                        v, f, ap, ds, hs, U = results[220]
                        prec = norm(tuple(f[j] - results[160][1][j] for j in range(2))) / (2 * U)
                        trunc = norm(tuple(f[j] - ap[j] for j in range(2))) / (2 * U)
                        for key, value in [('precision_force_over_2U', prec),
                                           ('neglected_force_over_2U', trunc),
                                           ('common_split_over_s', ds), ('B_over_s', hs)]:
                            maxima[key] = max(maxima[key], value)
                        exact_rows.append(dict(mG_eV=mg_text, phase=phase, index=idx,
                                               scale=scale_text, scenario=name,
                                               c=string(coefficients[0]), beta2=string(coefficients[1]),
                                               nu0=string(coefficients[2]), direct_potential_eV4=string(v),
                                               direct_chi=string(f[0]), direct_y=string(f[1]),
                                               leading_chi=string(ap[0]), leading_y=string(ap[1]),
                                               precision_force_over_2U=string(prec),
                                               neglected_force_over_2U=string(trunc),
                                               common_split_over_s=string(ds), B_over_s=string(hs)))
    record('160 to 220 digit direct force convergence', maxima['precision_force_over_2U'], '1e-85')
    record('direct CW omitted orders below budget', maxima['neglected_force_over_2U'], '1e-45')
    record('common mass split small', maxima['common_split_over_s'], '1e-45')
    record('B split small', maxima['B_over_s'], '1e-24')
    write_csv(OUT / 'independent_direct_cw.csv', exact_rows)
    out = dict(scope='Specified selected Q^2 determinant plus one finite local coefficient; not full SUGRA matching',
               direct_spectrum_points=len(exact_rows), kernel_points=len(basis_rows), budget_rows=len(budget_rows),
               precisions=[160, 220], independent_indices=list(INDICES), scales=list(SCALES),
               max_residuals={k: string(v) for k, v in maxima.items()},
               max_null_residual=string(max_null), max_rg_residual=string(max_rg),
               max_yjet_residual=string(max_yjet), singular_values=singular_rows,
               checks=CHECKS, passed=sum(c['passed'] for c in CHECKS), total=len(CHECKS),
               elapsed_seconds=time.time() - started,
               hashes={str(p.relative_to(REPO)): hashlib.sha256(p.read_bytes()).hexdigest()
                       for p in (Path(__file__), INPUT, TRAJECTORY, HERE / 'protocol.md')},
               python=platform.python_version(), mpmath=mp.__version__, numpy=np.__version__,
               production_code_imported_or_read=False,
               retained_diagnostics=[
                   dict(file='independent_initial_diagnostic.json', passed=21, total=22,
                        issue='80-digit differentiation of full |S|^2 lost significant digits',
                        residual='1.5165949892991368e-47', original_tolerance='1e-70',
                        correction='Use 160 digits for this independent differentiation; tolerance unchanged'),
                   dict(issue='Potential cancellation in initial quadrature Kc identified by code review',
                        correction='Algebraically factor sd*X cancellation before numerical evaluation')],
               tolerance_convention='Direct errors use common reference force 2U, including near-zero components; no relative division by a vanishing component.',
               scope_limits=['Loop-order counting is external to this determinant check.',
                             'No mixed gravity/fermion/ghost determinants or unknown Wilson coefficients are computed.',
                             'RG closure retains only terms linear in c and quadratic in beta.',
                             'The local-only check verifies normalization, not a prediction of nu.',
                             'Column-normalized singular values diagnose a two-shape degeneracy, not a protective symmetry.'])
    path = OUT / 'independent_summary.json'
    path.write_text(json.dumps(out, indent=2) + '\n')
    print(json.dumps({k: out[k] for k in ('passed', 'total', 'kernel_points', 'direct_spectrum_points',
                                         'max_residuals', 'max_null_residual', 'max_rg_residual', 'elapsed_seconds')}, indent=2))
    return out


def compare():
    """Compare archived production outputs after the independent derivation exists."""
    mp.mp.dps = 80
    checks = []

    def audit(name, residual, tolerance='1e-11'):
        checks.append(dict(name=name, residual=str(abs(residual)), tolerance=tolerance,
                           passed=bool(abs(residual) <= m(tolerance))))

    def read(name):
        with (OUT / name).open() as stream:
            return list(csv.DictReader(stream))

    def key(row):
        return float(row['mG_eV']), row['phase']

    independent = {(key(r), r['coefficient']): r for r in read('independent_budgets.csv')}
    primary_budgets = read('coefficient_budgets.csv')
    max_budget_error = mp.mpf(0)
    for row in primary_budgets:
        for coef, primary in [('c', 'c_10percent'), ('beta2', 'beta2_10percent'),
                              ('nu', 'abs_nu_10percent')]:
            own = m(independent[(key(row), coef)]['limit_10pct'])
            error = abs(m(row[primary]) / own - 1)
            max_budget_error = max(max_budget_error, error)
            audit(f'budget {key(row)} {coef}', error)

    own_kernels = {}
    grouped = {}
    for row in read('independent_kernels.csv'):
        group = key(row)
        idx = int(row['index'])
        kk = tuple(tuple(m(row[f'K_{coef}_{dr}_signed']) for dr in ('chi', 'y'))
                   for coef in ('c', 'beta2', 'nu'))
        own_kernels[(group, idx)] = (kk, m(row['U_eV4']))
        grouped.setdefault(group, []).append((idx, kk, m(row['U_eV4'])))

    primary_kernels = read('response_kernels.csv')
    primary_by_phase = {}
    for row in primary_kernels:
        primary_by_phase.setdefault(key(row), []).append(row)
    max_kernel_error = mp.mpf(0)
    for group, rows in primary_by_phase.items():
        for n, coef in enumerate(('c', 'beta2', 'nu')):
            unit_scale = max(norm(own_kernels[(group, int(r['index']))][0][n]) for r in rows)
            error = max(norm(tuple(m(r[f'{coef}_force_{dr}_eV4']) -
                                   own_kernels[(group, int(r['index']))][0][n][j]
                                   for j, dr in enumerate(('chi', 'y')))) / unit_scale
                        for r in rows)
            max_kernel_error = max(max_kernel_error, error)
            audit(f'production kernel global vector norm {group} {coef}', error)

    rg_rows = read('rg_closure.csv')
    rg_max_error = mp.mpf(0)
    for row in rg_rows:
        cc, bb, nn = (m(row[k]) for k in ('c', 'beta2', 'nu0'))
        log_mu = mp.log(m(row['mu_over_Mi']))
        # Under a scale change the retained explicit log shifts by -2 log(mu/Mi).
        # This gives the finite frozen-nu error without using production formulas.
        frozen_shift = 2 * (cc - bb / 2) * log_mu
        reference_max = mp.mpf(0)
        frozen_max = mp.mpf(0)
        for _, kk, U in grouped[key(row)]:
            reference = tuple(cc * kk[0][j] + bb * kk[1][j] + nn * kk[2][j]
                              for j in range(2))
            frozen = tuple(reference[j] + frozen_shift * kk[2][j] for j in range(2))
            reference_max = max(reference_max, norm(reference) / (2 * U))
            frozen_max = max(frozen_max, norm(frozen) / (2 * U))
        scale = max(reference_max, frozen_max, m('1e-100'))
        error = max(abs(m(row['reference_max_force_ratio']) - reference_max),
                    abs(m(row['matched_max_force_ratio']) - reference_max),
                    abs(m(row['incorrectly_frozen_nu_max_force_ratio']) - frozen_max)) / scale
        rg_max_error = max(rg_max_error, error)
        audit(f"RG archived output {key(row)} {row['case']} {row['mu_over_Mi']}", error)

    # Changed finite boundary conditions are model choices. Independently reproduce
    # their one-coefficient limits rather than interpreting them as scale variation.
    max_finite_error = mp.mpf(0)
    for row in primary_budgets:
        cmax = mp.mpf(0)
        bmax = mp.mpf(0)
        for _, kk, U in grouped[key(row)]:
            cmax = max(cmax, norm(tuple(kk[0][j] - kk[2][j] for j in range(2))) / (2 * U))
            bmax = max(bmax, norm(tuple(kk[1][j] - kk[2][j] / 2 for j in range(2))) / (2 * U))
        c_limit = m('.1') / cmax
        b_limit = mp.sqrt(m('.1') / bmax)
        for label, value, own in [('c finite boundary', 'c_limit_if_nu0_minus_c', c_limit),
                                  ('B finite boundary', 'beta_limit_if_nu0_minus_beta2_over2', b_limit)]:
            error = abs(m(row[value]) / own - 1)
            max_finite_error = max(max_finite_error, error)
            audit(f'{label} {key(row)}', error)
    comparison = dict(passed=sum(c['passed'] for c in checks), total=len(checks), checks=checks,
                      production_kernel_rows=len(primary_kernels), production_rg_rows=len(rg_rows),
                      max_budget_relative_error=string(max_budget_error),
                      max_kernel_global_norm_error=string(max_kernel_error),
                      max_rg_output_relative_error=string(rg_max_error),
                      max_finite_boundary_relative_error=string(max_finite_error),
                      production_code_imported_or_read=False,
                      production_output_hashes={name: hashlib.sha256((OUT / name).read_bytes()).hexdigest()
                                               for name in ('coefficient_budgets.csv', 'response_kernels.csv',
                                                            'rg_closure.csv')})
    (OUT / 'independent_comparison.json').write_text(json.dumps(comparison, indent=2) + '\n')
    print(json.dumps({k: v for k, v in comparison.items() if k != 'checks'}, indent=2))
    return comparison


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--derive', action='store_true')
    parser.add_argument('--compare', action='store_true')
    args = parser.parse_args()
    summaries = []
    if args.derive or not args.compare:
        summaries.append(derive())
    if args.compare:
        summaries.append(compare())
    raise SystemExit(0 if all(s['passed'] == s['total'] for s in summaries) else 1)
