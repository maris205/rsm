#!/usr/bin/env python3
"""Conditional Wilson-coefficient sensitivities and a selected one-loop RG closure.

No new cosmology, no full SUGRA loop matching, and no observational limits.
"""
from pathlib import Path
from datetime import datetime, timezone
import csv
import hashlib
import json
import platform

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
REPO = ROOT.parents[1]
OUT = ROOT / 'results'
MASSES = (1e-6, 1., 1e5, 1e7)
PHASES = {'zero': (1., 0.), 'quadrature': (0., 1.), 'pi': (-1., 0.)}
MU_FACTORS = (.5, 1., 2., 10.)
MI = 1e11
BUDGET = .1


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def dump(path, value):
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + '\n')


def csvout(path, rows):
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def inputs():
    source = REPO / 'experiments/sugra_shift_v01/results'
    par = json.loads((source / 'inputs.json').read_text())
    with (source / 'trajectories.csv').open() as stream:
        old = [row for row in csv.DictReader(stream) if row['case'] == 'shift_axis']
    chi = np.array([float(row['chi']) for row in old])
    return par, chi


def fields(mg, phase, par, chi):
    P = par['Mpl_eV'] ** 2
    F = np.sqrt(par['F_squared_eV2'])
    xi = par['xi']
    U = par['potential_unit_eV4'] * par['W_i'] * np.exp(-2 * (chi - par['chi_i']))
    wm = F * np.sqrt(U / 2)
    rw = wm / (mg * P)
    co, si = PHASES[phase]
    X = mg ** 2 * ((co - rw) ** 2 + si ** 2)
    gradX = mg ** 2 * np.stack((2 * rw * (co - rw), 2 * rw * si / 3))
    relative_s = xi * (chi ** -2 - par['chi_i'] ** -2) / (1 + xi / par['chi_i'] ** 2)
    s = MI ** 2 * (1 + relative_s)
    L = np.log1p(relative_s)
    r = -xi / (chi * (chi ** 2 + xi))
    grads = np.stack((2 * r * s, np.zeros_like(s)))
    p = s * X
    gradp = grads * X + s * gradX
    gradL = grads / s
    # Stable c-kernel, including the nonzero X jet when L(initial) is exactly zero.
    kc = -(X * L * grads + s * (L - 1) * gradX) / (8 * np.pi ** 2)
    kb = (L * gradp + p * gradL) / (16 * np.pi ** 2)
    kn = gradp / (8 * np.pi ** 2)
    return dict(U=U, wm=wm, rw=rw, s=s, L=L, X=X, gradX=gradX,
                grads=grads, p=p, gradp=gradp, gradL=gradL, r=r, kc=kc, kb=kb, kn=kn)


def force(c, beta2, nu0, factor, f, run=True):
    A = -c + beta2 / 2
    nu = nu0 + (2 * A * np.log(factor) if run else 0.)
    L = f['L'] - 2 * np.log(factor)
    return (f['gradp'] * (A * L + c + nu) + f['p'] * A * f['gradL']) / (8 * np.pi ** 2)


def norm_ratio(gradient, U):
    return np.hypot(gradient[0], gradient[1]) / (2 * U)


def baseline_b_interference(mg, phase, beta, par, chi, f):
    """Leading-center B from the parent's scalar-retuned eta, plus its 2D jet.

    This bounds interference at the one-coefficient beta budget, not at arbitrary beta.
    It is a local expansion check, not a completed new valley solution.
    """
    P = par['Mpl_eV'] ** 2
    F = np.sqrt(par['F_squared_eV2'])
    rho = 3 * par['omega_lambda'] * P * par['H_ref_eV'] ** 2
    vT0 = -3 * 48 ** 2 * mg ** 4 / (64 * np.pi ** 2)
    eta = (rho - vT0) / (6 * mg ** 2 * P)
    co, si = PHASES[phase]
    ft = mg * ((co - f['rw']) - 1j * si)
    M = np.sqrt(f['s'])
    hclock = np.sqrt(2 * f['U']) / F
    r = f['r']
    xi = par['xi']
    rp = xi * (3 * chi ** 2 + xi) / (chi ** 2 * (chi ** 2 + xi) ** 2)
    B0 = M * (hclock * r + 2 * eta * ft)
    B0x = M * (hclock * (r * r + rp - r) + 2 * eta * (r * ft + mg * f['rw']))
    B0y = 1j * M * (hclock * (r * r + rp + r) + 2 * eta * (-mg * f['rw'] / 3 - r * ft))
    delta = beta * M * ft
    dx = beta * M * (r * ft + mg * f['rw'])
    dy = 1j * beta * M * (r * ft - mg * f['rw'] / 3)
    cross = 2 * np.real(np.conj(B0) * delta)
    gradcross = 2 * np.real(np.stack((np.conj(B0x) * delta + np.conj(B0) * dx,
                                    np.conj(B0y) * delta + np.conj(B0) * dy)))
    gradient = (f['L'] * gradcross + cross * f['gradL']) / (16 * np.pi ** 2)
    return float(np.max(norm_ratio(gradient, f['U']))), float(np.max(np.abs(B0) / np.abs(delta)))


def main():
    OUT.mkdir(exist_ok=True, parents=True)
    par, chi = inputs()
    rows, kernels, rgrows, degeneracy, checks = [], [], [], [], []

    def check(name, value, tolerance):
        checks.append(dict(name=name, value=float(value), tolerance=float(tolerance),
                           passed=bool(np.isfinite(value) and value < tolerance)))

    for mg in MASSES:
        for phase in PHASES:
            f = fields(mg, phase, par, chi)
            scaled = [k / (2 * f['U']) for k in (f['kc'], f['kb'], f['kn'])]
            unit = np.array([np.max(np.hypot(k[0], k[1])) for k in scaled])
            c_limit, beta2_limit, nu_limit = BUDGET / unit
            b_limit = np.sqrt(beta2_limit)
            prefix = f'{mg:g}_{phase}'
            vectors = np.stack([k.flatten() for k in scaled], axis=1)
            vectors /= np.linalg.norm(vectors, axis=0)
            singular = np.linalg.svd(vectors, compute_uv=False)
            null_residual = np.max(np.abs(f['kc'] + 2 * f['kb'] - f['kn'])) / np.max(np.abs(f['kn']))
            interference, relative_b = baseline_b_interference(mg, phase, b_limit, par, chi, f)
            row = dict(mG_eV=mg, phase=phase, coordinates=len(chi),
                       c_unit_max_force_ratio=unit[0], beta2_unit_max_force_ratio=unit[1],
                       nu_unit_max_force_ratio=unit[2], c_10percent=c_limit,
                       abs_beta_10percent=b_limit, beta2_10percent=beta2_limit,
                       abs_nu_10percent=nu_limit,
                       c_limit_if_nu0_minus_c=BUDGET / np.max(norm_ratio(f['kc'] - f['kn'], f['U'])),
                       beta_limit_if_nu0_minus_beta2_over2=np.sqrt(BUDGET / np.max(norm_ratio(f['kb'] - f['kn'] / 2, f['U']))),
                       max_abs_d_over_s_at_c_limit=float(np.max(c_limit * f['X'] / f['s'])),
                       max_abs_B_over_s_at_beta_limit=float(np.max(b_limit * np.sqrt(f['X'] / f['s']))),
                       baseline_B_cross_max_force_ratio=interference,
                       max_abs_baseline_B_over_added_B=relative_b,
                       max_modulus_mass_over_charged_mass=float(np.max(np.sqrt(48 * f['X'] / f['s']))),
                       null_relative_residual=null_residual)
            rows.append(row)
            degeneracy.append(dict(mG_eV=mg, phase=phase, normalized_singular_values=singular.tolist(),
                                   null_vector_c_beta2_nu=[1, 2, -1], null_relative_residual=null_residual,
                                   interpretation='Only A=-c+beta²/2 and D=c+nu0 enter this selected leading potential.'))
            check(prefix + '_mass_normalization', abs(f['s'][0] / MI ** 2 - 1), 1e-14)
            check(prefix + '_kernel_null', null_residual, 1e-14)
            check(prefix + '_small_d_split', row['max_abs_d_over_s_at_c_limit'], 1e-8)
            check(prefix + '_small_B_split', row['max_abs_B_over_s_at_beta_limit'], 1e-8)
            check(prefix + '_B_interference_over_budget', interference / BUDGET, 1e-12)
            check(prefix + '_threshold_order', row['max_modulus_mass_over_charged_mass'], .01)
            check(prefix + '_numerical_rank2', singular[2] / singular[0], 1e-13)
            for label, c, b2, nu in [('c', c_limit, 0., 0.), ('beta2', 0., beta2_limit, 0.),
                                     ('nu', 0., 0., nu_limit),
                                     ('mixed', .2 * c_limit, .3 * beta2_limit, .1 * nu_limit)]:
                # The stable kernel form is the reference at mu0.
                reference = c * f['kc'] + b2 * f['kb'] + nu * f['kn']
                reference_budget = float(np.max(norm_ratio(reference, f['U'])))
                for factor in MU_FACTORS:
                    running = force(c, b2, nu, factor, f)
                    frozen = force(c, b2, nu, factor, f, run=False)
                    error = float(np.max(norm_ratio(running - reference, f['U'])))
                    rgrows.append(dict(mG_eV=mg, phase=phase, case=label, mu_over_Mi=factor,
                                       c=c, beta2=b2, nu0=nu,
                                       reference_max_force_ratio=reference_budget,
                                       matched_max_force_ratio=float(np.max(norm_ratio(running, f['U']))),
                                       incorrectly_frozen_nu_max_force_ratio=float(np.max(norm_ratio(frozen, f['U']))),
                                       matched_max_absolute_error_over_tree_force=error))
                    check(prefix + f'_RG_{label}_{factor}', error / BUDGET, 1e-11)
            if mg == 1.:
                for i in range(len(chi)):
                    item = dict(mG_eV=mg, phase=phase, index=i, chi=chi[i], U_eV4=f['U'][i],
                                s_eV2=f['s'][i], X_eV2=f['X'][i], L_at_Mi=f['L'][i])
                    for label, k in zip(('c', 'beta2', 'nu'), (f['kc'], f['kb'], f['kn'])):
                        item[label + '_force_chi_eV4'] = k[0, i]
                        item[label + '_force_y_eV4'] = k[1, i]
                    kernels.append(item)

    # Cross-mass scaling is a diagnostic of the stated large-Wc approximation.
    for phase in PHASES:
        selected = [r for r in rows if r['phase'] == phase]
        base = next(r for r in selected if r['mG_eV'] == 1.)
        for r in selected:
            check(f'scaling_c_{phase}_{r["mG_eV"]}', abs(r['c_10percent'] * r['mG_eV'] ** 2 / base['c_10percent'] - 1), 1e-12)
            check(f'scaling_beta_{phase}_{r["mG_eV"]}', abs(r['abs_beta_10percent'] * r['mG_eV'] / base['abs_beta_10percent'] - 1), 1e-12)

    csvout(OUT / 'coefficient_budgets.csv', rows)
    csvout(OUT / 'response_kernels.csv', kernels)
    csvout(OUT / 'rg_closure.csv', rgrows)
    dump(OUT / 'identifiability.json', degeneracy)
    dump(OUT / 'inputs.json', dict(parent_parameters=par, masses_eV=MASSES, phases=PHASES,
                                  matching_scale_factors=MU_FACTORS, initial_charged_mass_eV=MI,
                                  budget=BUDGET, chi_points=len(chi),
                                  default_finite_boundary='nu0=0 at mu0=Mi; one coefficient varied at a time',
                                  scope='Local-center selected leading EFT susceptibility, not full loop matching or a new cosmic trajectory.'))
    sources = [Path(__file__), ROOT / 'protocol.md',
               REPO / 'experiments/sugra_shift_v01/results/inputs.json',
               REPO / 'experiments/sugra_shift_v01/results/trajectories.csv']
    summary = dict(created_utc=datetime.now(timezone.utc).isoformat(), python=platform.python_version(),
                   numpy=np.__version__, source_hashes={str(p.relative_to(REPO)): sha(p) for p in sources},
                   budget_rows=len(rows), kernel_rows=len(kernels), rg_rows=len(rgrows),
                   checks_passed=sum(c['passed'] for c in checks), check_count=len(checks), checks=checks,
                   example_1eV_zero=next(r for r in rows if r['mG_eV'] == 1. and r['phase'] == 'zero'),
                   max_RG_error_over_tree_force=max(r['matched_max_absolute_error_over_tree_force'] for r in rgrows),
                   max_baseline_B_cross_force_ratio=max(r['baseline_B_cross_max_force_ratio'] for r in rows),
                   actual_model_assumptions='Finite local Wilson coefficients not fixed by previous vacuum retuning.',
                   unfinished=['full O(Q²) matter matching', 'wavefunction and auxiliary 1PI terms',
                               'full SUGRA loops and mixed two-loop contributions', 'UV coefficient protection',
                               'new valley/cosmic solution', 'real Standard Model observable'])
    dump(OUT / 'summary.json', summary)
    print('Checks', summary['checks_passed'], '/', summary['check_count'])
    print(json.dumps(summary['example_1eV_zero'], indent=2))
    if not all(c['passed'] for c in checks):
        raise SystemExit(1)


if __name__ == '__main__':
    main()
