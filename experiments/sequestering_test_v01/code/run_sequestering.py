#!/usr/bin/env python3
"""Selected finite 5D threshold embedded conditionally in a local 4D clock model.

No full UV embedding, complete two-loop matching, or new cosmic trajectory.
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
MG = np.logspace(-24, 2, 105)
KK = np.logspace(11, 14, 61)
KAPPAS = (-1, 0, 1)
PHASES = {'zero': (1., 0.), 'quadrature': (0., 1.), 'pi': (-1., 0.)}
C5 = float(zeta(3, 1) / (48 * np.pi ** 4))
MI = 1e11


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
    parent = REPO / 'experiments/sugra_shift_v01/results'
    par = json.loads((parent / 'inputs.json').read_text())
    with (parent / 'trajectories.csv').open() as stream:
        rows = [r for r in csv.DictReader(stream) if r['case'] == 'shift_axis']
    chi = np.array([float(r['chi']) for r in rows])
    P = par['Mpl_eV'] ** 2
    par['rho_lambda_eV4'] = 3 * par['omega_lambda'] * P * par['H_ref_eV'] ** 2
    U = par['potential_unit_eV4'] * par['W_i'] * np.exp(-2 * (chi - par['chi_i']))
    return par, chi, U


def kernels(mg, phase, par, chi, U):
    P = par['Mpl_eV'] ** 2
    F = np.sqrt(par['F_squared_eV2'])
    co, si = PHASES[phase]
    xi = par['xi']
    rw = F * np.sqrt(U / 2) / (mg * P)
    X = mg ** 2 * ((co - rw) ** 2 + si ** 2)
    dX = mg ** 2 * np.stack((2 * rw * (co - rw), 2 * rw * si / 3))
    relative_s = xi * (chi ** -2 - par['chi_i'] ** -2) / (1 + xi / par['chi_i'] ** 2)
    s = MI ** 2 * (1 + relative_s)
    L = np.log1p(relative_s)
    r = -xi / (chi * (chi ** 2 + xi))
    ds = np.stack((2 * r * s, np.zeros_like(s)))
    dp = X * ds + s * dX
    kc = -(X * L * ds + s * (L - 1) * dX) / (8 * np.pi ** 2)
    kb = (L * dp + X * ds) / (16 * np.pi ** 2)
    kn = dp / (8 * np.pi ** 2)
    baseline_tree = par['rho_lambda_eV4'] * np.stack((2 * rw * (co - rw), 2 * rw * si / 3))
    return dict(kc=kc, kb=kb, kn=kn, baseline_tree=baseline_tree,
                s=s, X=X, dX=dX, ds=ds, r=r, rw=rw, L=L)


def ratio(gradient, U):
    return np.hypot(gradient[0], gradient[1]) / (2 * U)


def record(mg, kk, kappa, phase, par, chi, U, f):
    P = par['Mpl_eV'] ** 2
    F = np.sqrt(par['F_squared_eV2'])
    rho = par['rho_lambda_eV4']
    A = C5 * kk ** 2 / P
    c = (12 * A + 3 * A * A) / (1 + A) ** 2
    beta = -6 * A / (1 + A)
    qforce = c * f['kc']
    changed_boundary = c * (f['kc'] - f['kn'])
    # One-loop metric times pre-existing tree-clock B: selected underlying two-loop interference.
    M = np.sqrt(f['s'])
    co, si = PHASES[phase]
    ft = mg * (co - f['rw'] - 1j * si)
    ft_x = mg * f['rw']
    ft_y = -1j * mg * f['rw'] / 3
    r = f['r']
    xi = par['xi']
    rp = xi * (3 * chi ** 2 + xi) / (chi ** 2 * (chi ** 2 + xi) ** 2)
    hc = np.sqrt(2 * U) / F
    # The baseline is tree vacuum retuning only. Scalar-loop eta from the parent is higher order here.
    eta_rho = rho / (6 * mg ** 2 * P)
    B0 = M * (hc * r + 2 * eta_rho * ft)
    B0x = M * (hc * (r * r + rp - r) + 2 * eta_rho * (r * ft + ft_x))
    B0y = 1j * M * (hc * (r * r + rp + r) + 2 * eta_rho * (-mg * f['rw'] / 3 - r * ft))
    dB = beta * M * ft
    dBx = beta * M * (r * ft + ft_x)
    dBy = beta * M * (1j * r * ft + ft_y)
    cross = 2 * np.real(np.conj(B0) * dB)
    crossjet = 2 * np.real(np.stack((np.conj(B0x) * dB + np.conj(B0) * dBx,
                                   np.conj(B0y) * dB + np.conj(B0) * dBy)))
    crossforce = (f['L'] * crossjet + cross * f['ds'] / f['s']) / (16 * np.pi ** 2)
    # Additional transverse clock auxiliary jet at q=0 and leading in the new metric.
    # This has zero axis value, so it is distinct from the hQ-induced delta B above.
    hZ_logprime = -3 * kappa * A / (1 + kappa * A)
    extra_B0_y = 2j * M * (r * ft - mg * f['rw'] / 3) * hZ_logprime
    hZ_Bjet_yforce = 2 * np.real(np.conj(B0) * extra_B0_y) * f['L'] / (16 * np.pi ** 2)
    extra_y = -6 * np.sqrt(2) * kappa * A * F * np.sqrt(U) * mg * si
    # Finite-rho, finite-clock local Hessian at h=1, retaining the 1/epsilon-enhanced
    # interference omitted by 4(U+rho)/(3P). Corrections from the tiny radial shift
    # and from A*rho are independently checked with the full K/W potential.
    rrho = rho / (3 * P * mg ** 2)
    wm = F * np.sqrt(U / 2)
    wshape = (co - f['rw']) ** 2 + si ** 2
    base_mass2 = (4 * (1 + rrho) * U / (3 * P)
                  - 2 * rrho * mg * wm * co / F ** 2
                  + 2 * rrho * mg ** 2 * (2 + rrho) * wshape)
    simpler_mass2 = 4 * (U + rho) / (3 * P)
    delta_mass2 = -24 * kappa * A * mg ** 2
    my2 = base_mass2 + delta_mass2
    min_my2 = float(np.min(my2))
    M5 = (P * kk / np.pi) ** (1 / 3)
    modulus = np.sqrt(48) * mg
    hierarchy = max(MI / kk, modulus / kk, mg / kk, kk / M5)
    qnorm = float(np.max(ratio(qforce, U)))
    changednorm = float(np.max(ratio(changed_boundary, U)))
    bnorm = float(np.max(ratio(crossforce, U)))
    # Growth time applies to a frozen local linear system, not a nonlinear cosmic trajectory.
    tau = par['hbar_eV_s'] / np.sqrt(-min_my2) if min_my2 < 0 else None
    if min_my2 < 0:
        H = par['H_ref_eV']
        gamma = 2 * (-min_my2) / (np.sqrt(9 * H * H + 4 * (-min_my2)) + 3 * H)
        damped_tau = par['hbar_eV_s'] / gamma
    else:
        damped_tau = None
    return dict(mG_eV=float(mg), mKK_eV=float(kk), kappa=int(kappa), phase=phase,
                A_5d=float(A), c_Q=float(c), beta_Q=float(beta),
                selected_c_2loop_force_ratio=qnorm,
                changed_finite_boundary_force_ratio=changednorm,
                treeB_cross_2loop_force_ratio=bnorm,
                hclock_Bjet_2loop_force_ratio=float(np.max(np.abs(hZ_Bjet_yforce) / (2 * U))),
                beta_squared_3loop_force_ratio=float(np.max(ratio(beta ** 2 * f['kb'], U))),
                hclock_center_force_ratio=float(np.max(np.abs(extra_y) / (2 * U))),
                rho_center_force_ratio=float(np.max(ratio(f['baseline_tree'], U))),
                min_my2_eV2=min_my2, max_my2_eV2=float(np.max(my2)),
                min_my2_over_Href2=min_my2 / par['H_ref_eV'] ** 2,
                finite_clock_mass_correction_over_simple=float(np.max(np.abs(base_mass2 - simpler_mass2) / simpler_mass2)),
                local_growth_efold_seconds=float(tau) if tau is not None else None,
                constant_Href_growth_efold_seconds=float(damped_tau) if damped_tau is not None else None,
                mass_product_eV2=float(mg * kk),
                eta_leading=float(eta_rho + 3 * A),
                delta_t_leading=float(eta_rho / 6 - A),
                M5_eV=float(M5), modulus_mass_eV=float(modulus),
                charged_over_KK=float(MI / kk), KK_over_M5=float(kk / M5),
                max_hierarchy_ratio=float(hierarchy),
                hierarchy_pass=bool(hierarchy <= .1 * (1 + 1e-14)),
                selected_force_pass=bool(qnorm <= .1),
                changed_boundary_force_pass=bool(changednorm <= .1),
                local_y_stable=bool(min_my2 > 0),
                necessary_tests_intersection=bool(qnorm <= .1 and min_my2 > 0 and hierarchy <= .1 * (1 + 1e-14)))


def main():
    OUT.mkdir(exist_ok=True, parents=True)
    par, chi, U = inputs()
    checks, grid, reps = [], [], []

    def check(name, value, tolerance):
        checks.append(dict(name=name, value=float(value), tolerance=float(tolerance),
                           passed=bool(np.isfinite(value) and value < tolerance)))

    for mg in MG:
        f = kernels(float(mg), 'zero', par, chi, U)
        for kk in KK:
            for kappa in KAPPAS:
                grid.append(record(float(mg), float(kk), kappa, 'zero', par, chi, U, f))
    for mg in (1e-24, 1e-18, 1e-16, 1.):
        for phase in PHASES:
            f = kernels(mg, phase, par, chi, U)
            for kk in (1e11, 1e12, 1e13, 1e14):
                for kappa in KAPPAS:
                    reps.append(record(mg, kk, kappa, phase, par, chi, U, f))
    P = par['Mpl_eV'] ** 2
    r1 = kernels(1., 'zero', par, chi, U)
    unit_c = float(np.max(ratio(r1['kc'], U)))
    unit_shifted = float(np.max(ratio(r1['kc'] - r1['kn'], U)))
    product_force = np.sqrt(.1 * P / (12 * C5 * unit_c))
    product_changed = np.sqrt(.1 * P / (12 * C5 * unit_shifted))
    product_stable = np.sqrt((np.min(U) + par['rho_lambda_eV4']) / (18 * C5))
    product_rho_only = np.sqrt(par['rho_lambda_eV4'] / (18 * C5))
    limits = [dict(mKK_eV=kk, c_force_mG_10percent_eV=product_force / kk,
                   changed_boundary_mG_10percent_eV=product_changed / kk,
                   kappa1_y_stability_mG_eV=product_stable / kk,
                   stability_rho_only_mG_eV=product_rho_only / kk,
                   charged_threshold_hierarchy_pass=bool(MI / kk <= .1))
              for kk in (1e11, 1e12, 1e13, 1e14)]

    check('grid_count', abs(len(grid) - 105 * 61 * 3), .5)
    check('representative_count', abs(len(reps) - 4 * 3 * 4 * 3), .5)
    check('A_small', max(r['A_5d'] for r in grid), 1e-20)
    check('local_vacuum_shift_small', max(abs(r['delta_t_leading']) for r in grid), 1e-8)
    check('eta_small', max(abs(r['eta_leading']) for r in grid), 1e-8)
    check('canonical_C5', abs(C5 / .0002570894757699417 - 1), 1e-14)
    check('prior_c_kernel', abs(unit_c / 8.116894967570333e30 - 1), 1e-13)
    check('source_radion_beta_c_relation', max(abs(r['beta_Q'] / r['c_Q'] + .5) for r in grid), 1e-13)
    check('kappa0_positive_y', 0 if all(r['local_y_stable'] for r in grid if r['kappa'] == 0) else 1, .5)
    check('negative_kappa_control_positive', 0 if all(r['local_y_stable'] for r in grid if r['kappa'] == -1) else 1, .5)
    check('kappa1_both_outcomes', 0 if {r['local_y_stable'] for r in grid if r['kappa'] == 1} == {True, False} else 1, .5)
    for kk in (1e11, 1e12, 1e13, 1e14):
        for factor in (.5, 2.):
            mg = factor * product_stable / kk
            rec = record(mg, kk, 1, 'zero', par, chi, U, kernels(mg, 'zero', par, chi, U))
            check(f'stability_boundary_{kk}_{factor}', 0 if rec['local_y_stable'] == (factor < 1) else 1, .5)
        for phase in PHASES:
            sample = next(r for r in reps if r['mG_eV'] == 1. and r['mKK_eV'] == kk and r['kappa'] == 1 and r['phase'] == phase)
            check(f'phase_does_not_cure_tachyon_{kk}_{phase}', 0 if not sample['local_y_stable'] else 1, .5)
    for r in reps:
        prefix = f'{r["mG_eV"]}_{r["mKK_eV"]}_{r["kappa"]}_{r["phase"]}'
        # First-order algebra has a known tiny clock-superpotential correction.
        expected = .1 * (r['mass_product_eV2'] / product_force) ** 2
        check(prefix + '_force_product_scaling', abs(r['selected_c_2loop_force_ratio'] / expected - 1), 1e-8)
    sample1 = next(r for r in reps if r['mG_eV'] == 1. and r['mKK_eV'] == 1e12 and r['kappa'] == 1 and r['phase'] == 'zero')
    samplelow = next(r for r in reps if r['mG_eV'] == 1e-18 and r['mKK_eV'] == 1e12 and r['kappa'] == 1 and r['phase'] == 'zero')
    check('growth_time_1eV_1TeV', abs(sample1['local_growth_efold_seconds'] / 20.404077416553772 - 1), 1e-12)
    check('low_mass_necessary_tests', 0 if samplelow['necessary_tests_intersection'] else 1, .5)
    csvout(OUT / 'scan.csv', grid)
    csvout(OUT / 'representatives.csv', reps)
    csvout(OUT / 'conditional_limits.csv', limits)
    dump(OUT / 'inputs.json', dict(parent_parameters=par, C5=C5, masses_eV=MG.tolist(), KK_eV=KK.tolist(),
                                  kappas=KAPPAS, phases=PHASES, chi_points=len(chi), initial_mass_eV=MI,
                                  hierarchy_ratio_limit=.1, force_budget=.1,
                                  main_finite_potential_boundary='nu0=0; selected underlying two-loop contribution only',
                                  alternative_finite_potential_boundary='nu0=-c_Q; different matching assumption',
                                  status='Conditional 4D graft of a finite 5D sector, not a complete UV embedding.'))
    sources = [Path(__file__), ROOT / 'protocol.md',
               REPO / 'experiments/sugra_shift_v01/results/inputs.json',
               REPO / 'experiments/sugra_shift_v01/results/trajectories.csv']
    summary = dict(created_utc=datetime.now(timezone.utc).isoformat(), python=platform.python_version(),
                   numpy=np.__version__, scipy=scipy.__version__, source_hashes={str(p.relative_to(REPO)): sha(p) for p in sources},
                   scan_rows=len(grid), representative_rows=len(reps),
                   checks_passed=sum(c['passed'] for c in checks), check_count=len(checks), checks=checks,
                   force_product_bound_eV2=float(product_force), changed_boundary_product_bound_eV2=float(product_changed),
                   kappa1_stability_product_bound_eV2=float(product_stable), rho_only_stability_product_bound_eV2=float(product_rho_only),
                   sample_1eV_1TeV=sample1, sample_1e18small_1TeV=samplelow,
                   pass_counts={str(k): sum(r['necessary_tests_intersection'] for r in grid if r['kappa'] == k) for k in KAPPAS},
                   interpretation='local_y_stable means positive Vyy/Gyy on the stated local slice; pass is only intersection of necessary local tests, not full covariant rolling-background stability, full protection or a new cosmic solution.',
                   unresolved=['shift-clock brane coefficient is an added assumption', 'radion axion quartic stabilization needs new sector',
                               'co-located clock-matter local operators', 'complete two-loop finite matching',
                               'real Standard Model and cosmological solution', 'Riemann origin of interactions and inverse-log exponent'])
    dump(OUT / 'summary.json', summary)
    print('Checks', summary['checks_passed'], '/', summary['check_count'])
    print('Rows', len(grid), len(reps))
    print('Product bounds', product_force, product_stable)
    print('1 eV sample', sample1)
    print('Small-mass sample', samplelow)
    if not all(c['passed'] for c in checks):
        raise SystemExit(1)


if __name__ == '__main__':
    main()
