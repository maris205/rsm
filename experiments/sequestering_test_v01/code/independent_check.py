#!/usr/bin/env python3
"""Independent source conversion, K/W Hessians, and selected threshold checks.

No import of the primary implementation.  A retained finite 5D coefficient does
not constitute full matching of the hybrid 4D stabilization/clock theory.
"""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import platform
from datetime import datetime, timezone

import mpmath as mp

ROOT = Path(__file__).resolve().parents[1]
REPO = ROOT.parents[1]
OUT = ROOT / 'results'
PARENT = REPO / 'experiments/sugra_shift_v01/results'


def dump(name, obj):
    (OUT / name).write_text(json.dumps(obj, ensure_ascii=False, indent=2) + '\n')


def write_csv(name, rows):
    with (OUT / name).open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def text(x):
    return mp.nstr(x, 90)


def parameters():
    raw = json.loads((PARENT / 'inputs.json').read_text())
    par = {k: mp.mpf(str(v)) for k, v in raw.items() if isinstance(v, (int, float))}
    with (PARENT / 'trajectories.csv').open() as stream:
        chi = [mp.mpf(r['chi']) for r in csv.DictReader(stream) if r['case'] == 'shift_axis']
    par['P'] = par['Mpl_eV'] ** 2
    par['epsilon_exact'] = par['F_squared_eV2'] / par['P']
    par['rho'] = 3 * par['omega_lambda'] * par['P'] * par['H_ref_eV'] ** 2
    par['Ui'] = par['potential_unit_eV4'] * par['W_i']
    return par, chi


def coeff_A(mkk, par):
    return mp.zeta(3) * mkk ** 2 / (48 * mp.pi ** 4 * par['P'])


def gfun(t, A, eta, v=0):
    return t - A / t ** 2 + eta * (t - 1) ** 2 + (t - 1) ** 4 + v ** 4


def hfun(t, A, kappa):
    return 1 + kappa * A / t ** 3


def vacuum(A, mg, par):
    """Once fix eta and solve the clock-off stationary radius at V=rho."""
    target = par['rho'] / (mg ** 2 * par['P'])

    def potential(t, eta):
        g = gfun(t, A, eta)
        p = 1 + 2 * A / t ** 3 + 2 * eta * (t - 1) + 4 * (t - 1) ** 3
        q = -6 * A / t ** 4 + 2 * eta + 12 * (t - 1) ** 2
        return 3 * q / (g ** 2 * (p ** 2 - g * q))

    eta0 = 3 * A + target / 6
    t0 = 1 - A + target / 36
    tstar, eta = mp.findroot(
        (lambda t, e: potential(t, e) * gfun(t, A, e) ** 3 / target - 1,
         lambda t, e: mp.diff(lambda x: potential(x, e), t) / max(A, target)),
        (t0, eta0), tol=mp.power(10, -mp.mp.dps + 30), maxsteps=40)
    gv = gfun(tstar, A, eta)
    return tstar, eta, lambda t, e: gv ** 3 * potential(t, e)


def normalized_parameters(par, tstar, A, eta, kappa):
    """Fix holomorphic inputs once using the physical clock-off vacuum."""
    local = dict(par)
    gv = gfun(tstar, A, eta)
    hv = hfun(tstar, A, kappa)
    local['epsilon_hol'] = gv * par['epsilon_exact'] / hv
    local['F2_hol'] = gv * par['F_squared_eV2'] / hv
    local['Ui_hol'] = gv ** 2 * hv * par['Ui']
    local['Wc_normalization'] = gv ** mp.mpf('1.5')
    return local


def potential_matrix(t, y, chi, mg, phase, A, eta, kappa, par, clock=True, v=0):
    """Direct two-field Kahler inverse; scaled potential V/(mg² P)."""
    eps = par['epsilon_hol']
    g = gfun(t, A, eta, v)
    p = 1 + 2 * A / t ** 3 + 2 * eta * (t - 1) + 4 * (t - 1) ** 3
    q = -6 * A / t ** 4 + 2 * eta + 12 * (t - 1) ** 2
    h = hfun(t, A, kappa)
    hp = -3 * kappa * A / t ** 4
    hpp = 12 * kappa * A / t ** 5
    Om = g - h * eps * y ** 2 / 3
    ot = p - 4j * v ** 3 - hp * eps * y ** 2 / 3
    oz = 1j * h * mp.sqrt(2 * eps) * y / 3
    ott = q + 12 * v ** 2 - hpp * eps * y ** 2 / 3
    otzb = -1j * hp * mp.sqrt(2 * eps) * y / 3
    ozz = -h / 3
    oi = mp.matrix([ot, oz])
    oij = mp.matrix([[ott, otzb], [mp.conj(otzb), ozz]])
    metric = -3 * (oij / Om - oi * oi.transpose_conj() / Om ** 2)
    kval = -3 * oi / Om
    co, si = {'zero': (1, 0), 'quadrature': (0, 1), 'pi': (-1, 0)}[phase]
    wc = mp.mpc(co, si) * par['Wc_normalization']
    if clock:
        U = par['Ui_hol'] * mp.exp(-2 * (chi - par['chi_i']))
        wr = -mp.sqrt(par['F2_hol'] * U / 2) / (mg * par['P'])
        wclock = wr * mp.exp(-1j * y)
        wz = -mp.sqrt(2 / eps) * wclock
    else:
        wclock, wz = mp.mpf('0'), mp.mpf('0')
    w = wc + wclock
    dw = kval * w + mp.matrix([0, wz])
    value = (dw.transpose_conj() * metric ** -1 * dw)[0] - 3 * abs(w) ** 2
    return mp.re(value / Om ** 3)


def clockoff_mass_formula(t, mg, A, eta, kappa):
    g = gfun(t, A, eta)
    p = 1 + 2 * A / t ** 3 + 2 * eta * (t - 1) + 4 * (t - 1) ** 3
    q = -6 * A / t ** 4 + 2 * eta + 12 * (t - 1) ** 2
    h = hfun(t, A, kappa)
    j = -3 * kappa * A / t ** 4
    k = 12 * kappa * A / t ** 5
    D = p ** 2 - g * q
    numerator = (-g * h ** 2 * q ** 2 - 2 * g * h * j * p * q
                 - g * h * k * p ** 2 + 2 * g * j ** 2 * p ** 2
                 + 2 * h ** 2 * p ** 2 * q)
    return 2 * mg ** 2 * g ** 3 * numerator / (g ** 2 * h ** 2 * D ** 2)


def direct_mass(t, chi, mg, phase, A, eta, kappa, par, clock):
    vyy = mp.diff(lambda y: potential_matrix(t, y, chi, mg, phase, A, eta,
                                             kappa, par, clock=clock), 0, 2)
    normalization = par['epsilon_hol'] * hfun(t, A, kappa) / gfun(t, A, eta)
    return mg ** 2 * vyy / normalization


def kernels(chi, mg, phase, par):
    """Response from analytic derivatives of the CW determinant expansion."""
    co, si = {'zero': (1, 0), 'quadrature': (0, 1), 'pi': (-1, 0)}[phase]
    U = par['Ui'] * mp.exp(-2 * (chi - par['chi_i']))
    wr = -mp.sqrt(par['F_squared_eV2'] * U / 2) / par['P']
    X = (mg * co + wr) ** 2 + (mg * si) ** 2
    Xc = -2 * wr * (mg * co + wr)
    Xy = -2 * mg * wr * si / 3
    s = mp.mpf('1e22') * (1 + par['xi'] / chi ** 2) / (1 + par['xi'] / par['chi_i'] ** 2)
    L = mp.log(s / mp.mpf('1e22'))
    sc = -2 * s * par['xi'] / (chi * (chi ** 2 + par['xi']))
    kc = [(-sc * X * L + s * Xc * (1 - L)) / (8 * mp.pi ** 2),
          s * Xy * (1 - L) / (8 * mp.pi ** 2)]
    kn = [(sc * X + s * Xc) / (8 * mp.pi ** 2), s * Xy / (8 * mp.pi ** 2)]
    kb = [(kn[i] - kc[i]) / 2 for i in range(2)]
    return U, kc, kn, kb


def mass_refined(U, mg, phase, A, kap, par):
    co, si = {'zero': (1, 0), 'quadrature': (0, 1), 'pi': (-1, 0)}[phase]
    rr = par['rho'] / (3 * par['P'] * mg ** 2)
    wm = mp.sqrt(par['F_squared_eV2'] * U / 2)
    rw = wm / (mg * par['P'])
    shape = (co - rw) ** 2 + si ** 2
    return (4 * (1 + rr) * U / (3 * par['P'])
            - 2 * rr * mg * wm * co / par['F_squared_eV2']
            + 2 * rr * mg ** 2 * (2 + rr) * shape - 24 * kap * A * mg ** 2)


def run_derive():
    if not (ROOT / 'protocol.md').is_file():
        raise RuntimeError('Production requires the experiment protocol first.')
    OUT.mkdir(parents=True, exist_ok=True)
    mp.mp.dps = 100
    par, chis = parameters()
    checks = []

    def check(name, value, tolerance):
        checks.append(dict(name=name, value=text(value), tolerance=text(tolerance),
                           passed=bool(mp.isfinite(value) and abs(value) < tolerance)))

    # Derive normalization from RSS F-frame coefficients, including M5_R³=2P/L0.
    mkk = mp.mpf('1e12')
    L0 = 2 * mp.pi / mkk
    m5cube = 2 * par['P'] / L0
    from_pure = (mp.zeta(3) / (4 * mp.pi ** 2 * L0 ** 2)) / (3 * par['P'])
    from_matter = mp.zeta(3) / (6 * mp.pi ** 2 * m5cube * L0 ** 3)
    expected = coeff_A(mkk, par)
    check('RSS_A_from_pure', abs(from_pure / expected - 1), mp.mpf('1e-90'))
    check('RSS_same_A_from_matter', abs(from_matter / expected - 1), mp.mpf('1e-90'))
    for kap in (-1, 0, 1):
        ell = lambda t: mp.log(hfun(t, expected, kap))
        d1 = mp.diff(ell, 1)
        d2 = mp.diff(ell, 1, 2)
        exact1 = -3 * kap * expected / (1 + kap * expected)
        exact2 = (12 * kap * expected + 3 * kap ** 2 * expected ** 2) / (1 + kap * expected) ** 2
        scale = max(abs(expected), mp.mpf('1e-100'))
        check(f'logmetric_first_{kap}', abs(d1 - exact1) / scale, mp.mpf('1e-65'))
        check(f'logmetric_second_{kap}', abs(d2 - exact2) / scale, mp.mpf('1e-65'))

    budgets = []
    for mg_s in ('1e-24', '1e-16', '1'):
        mg = mp.mpf(mg_s)
        for phase in ('zero', 'quadrature', 'pi'):
            ratios = []
            for chi in chis:
                U, kc, _, _ = kernels(chi, mg, phase, par)
                ratios.append(mp.sqrt(sum(v ** 2 for v in kc)) / (2 * U))
            peak = max(ratios)
            limit = mp.mpf('.1') / peak
            kk_limit = mp.sqrt(limit * par['P'] * 48 * mp.pi ** 4 / (12 * mp.zeta(3)))
            budgets.append(dict(mG_eV=mg_s, phase=phase, c_unit_max=text(peak),
                                c_10percent=text(limit), leading_mKK_10percent_eV=text(kk_limit),
                                leading_product_eV2=text(mg * kk_limit),
                                maximizing_index=ratios.index(peak)))

    rows_by_precision = {}
    for dps in (140, 220):
        mp.mp.dps = dps
        par, chis = parameters()
        rows = []
        for mg_s in ('1e-24', '1e-16', '1'):
            mg = mp.mpf(mg_s)
            for kk_s in ('1e12', '1e14'):
                A = coeff_A(mp.mpf(kk_s), par)
                tstar, eta, v0 = vacuum(A, mg, par)
                for kap in (-1, 0, 1):
                    local = normalized_parameters(par, tstar, A, eta, kap)
                    exact = direct_mass(tstar, chis[0], mg, 'zero', A, eta, kap, local, False)
                    analytic = clockoff_mass_formula(tstar, mg, A, eta, kap)
                    leading = 4 * par['rho'] / (3 * par['P']) - 24 * kap * A * mg ** 2
                    scale = abs(4 * par['rho'] / (3 * par['P'])) + abs(24 * kap * A * mg ** 2)
                    check(f'full_matrix_vs_formula_{dps}_{mg_s}_{kk_s}_{kap}',
                          abs(exact - analytic) / scale, mp.mpf('1e-65'))
                    rows.append(dict(dps=dps, case='clock_off', mG_eV=mg_s, mKK_eV=kk_s,
                                     kappa=kap, phase='zero', index=-1,
                                     A=text(A), delta=text(tstar - 1), eta=text(eta),
                                     mass_direct_eV2=text(exact), mass_formula_eV2=text(analytic),
                                     mass_leading_eV2=text(leading), normalization_eV2=text(scale),
                                     leading_normalized_difference=text(abs(exact - leading) / scale),
                                     mass_refined_eV2='', refined_normalized_difference='',
                                     tachyon_time_seconds=text(par['hbar_eV_s'] / mp.sqrt(-exact)) if exact < 0 else '',
                                     vacuum_relative_error=text(abs(v0(tstar, eta) * mg ** 2 * par['P'] / par['rho'] - 1)),
                                     finite_clock_y_slope_eV4='0'))
        # Independent selected finite-clock slices, keeping the once-fixed eta.
        for mg_s in ('1e-24', '1e-16', '1'):
            mg, kk_s = mp.mpf(mg_s), '1e12'
            A = coeff_A(mp.mpf(kk_s), par)
            tstar, eta, _ = vacuum(A, mg, par)
            for kap in (-1, 0, 1):
                local = normalized_parameters(par, tstar, A, eta, kap)
                for phase in ('zero', 'quadrature', 'pi'):
                    for index in (0, len(chis) - 1):
                        chi = chis[index]
                        vx = lambda t: potential_matrix(t, 0, chi, mg, phase, A, eta, kap, local)
                        tv = mp.findroot(lambda t: mp.diff(vx, t), tstar,
                                         tol=mp.power(10, -dps + 30), maxsteps=40)
                        exact = direct_mass(tv, chi, mg, phase, A, eta, kap, local, True)
                        U = par['Ui'] * mp.exp(-2 * (chi - par['chi_i']))
                        leading = 4 * (U + par['rho']) / (3 * par['P']) - 24 * kap * A * mg ** 2
                        scale = 4 * (U + par['rho']) / (3 * par['P']) + abs(24 * kap * A * mg ** 2)
                        refined = mass_refined(U, mg, phase, A, kap, par)
                        check(f'refined_mass_{dps}_{mg_s}_{kap}_{phase}_{index}',
                              abs(exact - refined) / scale, mp.mpf('1e-12'))
                        slope = mg ** 2 * par['P'] * mp.diff(
                            lambda y: potential_matrix(tv, y, chi, mg, phase, A, eta, kap, local), 0)
                        rows.append(dict(dps=dps, case='finite_clock_y0_slice', mG_eV=mg_s,
                                         mKK_eV=kk_s, kappa=kap, phase=phase, index=index,
                                         A=text(A), delta=text(tv - 1), eta=text(eta),
                                         mass_direct_eV2=text(exact), mass_formula_eV2='',
                                         mass_leading_eV2=text(leading), normalization_eV2=text(scale),
                                         leading_normalized_difference=text(abs(exact - leading) / scale),
                                         mass_refined_eV2=text(refined), refined_normalized_difference=text(abs(exact - refined) / scale),
                                         tachyon_time_seconds=text(par['hbar_eV_s'] / mp.sqrt(-exact)) if exact < 0 else '',
                                         vacuum_relative_error='', finite_clock_y_slope_eV4=text(slope)))
        rows_by_precision[dps] = rows
        print(f'Direct K/W batch {dps} digits: {len(rows)} points', flush=True)
    low, high = rows_by_precision[140], rows_by_precision[220]
    precision_errors = [abs(mp.mpf(l['mass_direct_eV2']) - mp.mpf(h['mass_direct_eV2'])) /
                        mp.mpf(h['normalization_eV2']) for l, h in zip(low, high)]
    check('140_to_220_mass_max_normalized_difference', max(precision_errors), mp.mpf('1e-65'))
    write_csv('independent_budgets.csv', budgets)
    write_csv('independent_kahler_140.csv', low)
    write_csv('independent_kahler_220.csv', high)
    summary = dict(created_utc=datetime.now(timezone.utc).isoformat(),
                   python_version=platform.python_version(), mpmath_version=mp.__version__,
                   code_independence='No import or execution of the primary implementation.',
                   scope='Direct 4D two-chiral-field K/W Hessians and selected coefficient susceptibility; not full loop matching.',
                   checks=checks, check_count=len(checks), checks_passed=sum(c['passed'] for c in checks),
                   direct_points_per_precision=len(high), precision_digits=[140, 220],
                   max_precision_normalized_difference=text(max(precision_errors)),
                   max_clockoff_leading_normalized_difference=max(float(r['leading_normalized_difference']) for r in high if r['case'] == 'clock_off'),
                   max_finiteclock_leading_normalized_difference=max(float(r['leading_normalized_difference']) for r in high if r['case'] != 'clock_off'),
                   max_finiteclock_refined_normalized_difference=max(float(r['refined_normalized_difference']) for r in high if r['case'] != 'clock_off'),
                   source_hashes={str(p.relative_to(REPO)): sha(p) for p in
                                  [Path(__file__), ROOT / 'protocol.md', PARENT / 'inputs.json', PARENT / 'trajectories.csv']})
    dump('independent_summary.json', summary)
    print(json.dumps({k: v for k, v in summary.items() if k not in ('checks', 'source_hashes')}, indent=2), flush=True)
    if not all(c['passed'] for c in checks):
        raise SystemExit(1)


def run_compare():
    """Compare output tables only, using independently implemented equations."""
    mp.mp.dps = 100
    par, chis = parameters()
    units = {}
    for phase in ('zero', 'quadrature', 'pi'):
        ratios = [[], [], []]
        for chi in chis:
            U, kc, kn, kb = kernels(chi, mp.mpf(1), phase, par)
            for storage, vector in zip(ratios, (kc, [kc[i] - kn[i] for i in range(2)], kb)):
                storage.append(mp.sqrt(sum(x ** 2 for x in vector)) / (2 * U))
        units[phase] = [max(v) for v in ratios]
    representative_units = {}

    def actual_units(mg, phase):
        key = (str(mg), phase)
        if key not in representative_units:
            ratios = [[], [], []]
            for chi in chis:
                U, kc, kn, kb = kernels(chi, mg, phase, par)
                for storage, vector in zip(ratios, (kc, [kc[i] - kn[i] for i in range(2)], kb)):
                    storage.append(mp.sqrt(sum(x ** 2 for x in vector)) / (2 * U))
            representative_units[key] = [max(v) for v in ratios]
        return representative_units[key]
    Uvalues = [par['Ui'] * mp.exp(-2 * (c - par['chi_i'])) for c in chis]
    Umin, Umax = min(Uvalues), max(Uvalues)
    errors = {}
    bool_failures = []
    checked_rows = 0

    def error(field, observed, expected, scale=None):
        residual = abs(mp.mpf(observed) - expected) / (scale if scale is not None else max(abs(expected), mp.mpf('1e-300')))
        errors[field] = max(errors.get(field, mp.mpf(0)), residual)

    for filename in ('scan.csv', 'representatives.csv'):
        with (OUT / filename).open() as stream:
            for n, row in enumerate(csv.DictReader(stream)):
                checked_rows += 1
                mg, kk = mp.mpf(row['mG_eV']), mp.mpf(row['mKK_eV'])
                kap, phase = int(row['kappa']), row['phase']
                A = coeff_A(kk, par)
                c, beta = (12 * A + 3 * A ** 2) / (1 + A) ** 2, -6 * A / (1 + A)
                response = actual_units(mg, phase)
                force = c * response[0]
                changed = c * response[1]
                term = -24 * kap * A * mg ** 2
                mymin = mass_refined(Umin, mg, phase, A, kap, par)
                mymax = mass_refined(Umax, mg, phase, A, kap, par)
                m5 = (par['P'] * kk / mp.pi) ** (mp.mpf(1) / 3)
                hierarchy = max(mp.mpf('1e11') / kk, mp.sqrt(48) * mg / kk, mg / kk, kk / m5)
                expected = {'A_5d': A, 'c_Q': c, 'beta_Q': beta,
                            'selected_c_2loop_force_ratio': force,
                            'changed_finite_boundary_force_ratio': changed,
                            'beta_squared_3loop_force_ratio': beta ** 2 * response[2],
                            'mass_product_eV2': mg * kk, 'eta_leading': 3 * A + par['rho'] / (6 * mg ** 2 * par['P']),
                            'delta_t_leading': -A + par['rho'] / (36 * mg ** 2 * par['P']),
                            'M5_eV': m5, 'modulus_mass_eV': mp.sqrt(48) * mg,
                            'charged_over_KK': mp.mpf('1e11') / kk, 'KK_over_M5': kk / m5,
                            'max_hierarchy_ratio': hierarchy}
                for key, value in expected.items():
                    error(key, row[key], value)
                mass_scale = 4 * (Umax + par['rho']) / (3 * par['P']) + abs(term)
                error('min_my2_eV2', row['min_my2_eV2'], mymin, mass_scale)
                error('max_my2_eV2', row['max_my2_eV2'], mymax, mass_scale)
                if mymin < 0:
                    error('local_growth_efold_seconds', row['local_growth_efold_seconds'],
                          par['hbar_eV_s'] / mp.sqrt(-mymin))
                booleans = {'hierarchy_pass': hierarchy <= mp.mpf('.1'),
                            'selected_force_pass': force <= mp.mpf('.1'),
                            'changed_boundary_force_pass': changed <= mp.mpf('.1'),
                            'local_y_stable': mymin > 0,
                            'necessary_tests_intersection': hierarchy <= mp.mpf('.1') and force <= mp.mpf('.1') and mymin > 0}
                for key, expected_bool in booleans.items():
                    if (row[key] == 'True') != bool(expected_bool):
                        bool_failures.append(dict(file=filename, row=n, field=key, actual=row[key], expected=bool(expected_bool)))
    with (OUT / 'conditional_limits.csv').open() as stream:
        for row in csv.DictReader(stream):
            kk = mp.mpf(row['mKK_eV'])
            A = coeff_A(kk, par)
            expected = {'c_force_mG_10percent_eV': mp.sqrt(mp.mpf('.1') / (12 * A * units['zero'][0])),
                        'changed_boundary_mG_10percent_eV': mp.sqrt(mp.mpf('.1') / (12 * A * units['zero'][1])),
                        'kappa1_y_stability_mG_eV': mp.sqrt((Umin + par['rho']) / (18 * par['P'] * A)),
                        'stability_rho_only_mG_eV': mp.sqrt(par['rho'] / (18 * par['P'] * A))}
            for key, value in expected.items():
                error(key, row[key], value)
    tolerance = mp.mpf('1e-11')
    checks = [dict(name=key, value=text(value), tolerance=text(tolerance), passed=bool(value < tolerance))
              for key, value in errors.items()]
    checks.append(dict(name='classification_mismatches', value=len(bool_failures), tolerance=0, passed=not bool_failures))
    summary = dict(created_utc=datetime.now(timezone.utc).isoformat(), checked_rows=checked_rows,
                   python_version=platform.python_version(), mpmath_version=mp.__version__,
                   checks=checks, checks_passed=sum(c['passed'] for c in checks), check_count=len(checks),
                   boolean_mismatches=bool_failures, max_normalized_difference=text(max(errors.values())),
                   scope='Main leading-order table comparison; no reading or importing main implementation.',
                   source_hashes={str(p.relative_to(REPO)): sha(p) for p in [Path(__file__), ROOT / 'protocol.md',
                                  OUT / 'scan.csv', OUT / 'representatives.csv', OUT / 'conditional_limits.csv']})
    dump('independent_comparison.json', summary)
    print(json.dumps({k: v for k, v in summary.items() if k not in ('checks', 'source_hashes')}, indent=2), flush=True)
    if not all(c['passed'] for c in checks):
        raise SystemExit(1)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--derive', action='store_true')
    parser.add_argument('--compare', action='store_true')
    args = parser.parse_args()
    if args.derive:
        run_derive()
    if args.compare:
        run_compare()
