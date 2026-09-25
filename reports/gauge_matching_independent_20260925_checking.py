#!/usr/bin/env python3
"""Independent KL threshold audit: reads only the archived SUGRA experiment.

The primary/new gauge matching implementation is deliberately not imported.
The algebra follows KL hep-th/9402005v2 (3.4),(3.7),(3.13),(3.19),(3.21)--(3.29).
This is a local one-loop matching audit, not a new background integration.
"""
from pathlib import Path
from datetime import datetime, timezone
import csv
import hashlib
import json
import platform

import mpmath as mp
import sympy as sy

REPO = Path(__file__).resolve().parents[1]
PREFIX = Path(__file__).with_suffix('')
ARCHIVE = REPO / 'experiments/sugra_shift_v01/results'
INPUTS = ARCHIVE / 'inputs.json'
TRAJECTORIES = ARCHIVE / 'trajectories.csv'
SOURCE = REPO / 'build/sugra_shift_theory_sources/kaplunovsky_louis_9402005v2.pdf'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    mp.mp.dps = 180
    checks = []

    def symbolic(name, expression):
        residual = sy.simplify(expression)
        checks.append({'name': name, 'passed': residual == 0, 'residual': str(residual)})

    def numerical(name, residual, tolerance):
        checks.append({'name': name, 'passed': bool(abs(residual) < tolerance),
                       'residual': mp.nstr(residual, 20),
                       'tolerance': mp.nstr(tolerance, 10)})

    # k=K/Mp^2; l=ln|Mhol|; z=ln(Z+ Z-). Coefficients use physical ln mass.
    k, ell, z, f, lam = sy.symbols('k ell z f lam', real=True)
    pi = sy.pi
    mlog = ell + k / 2 - z / 2
    uv = f + 2 * k / (16 * pi**2) - z / (8 * pi**2)
    threshold = 2 * (lam - mlog) / (8 * pi**2)
    low = f + 2 * (lam - ell) / (8 * pi**2)
    symbolic('KL_UV_plus_physical_threshold_equals_holomorphic_low', uv + threshold - low)
    symbolic('K_anomaly_cancels_mass_K', sy.diff(uv + threshold, k))
    symbolic('Konishi_cancels_mass_wavefunction', sy.diff(uv + threshold, z))
    symbolic('holomorphic_mass_log_coefficient', sy.diff(low, ell) + 1 / (4 * pi**2))
    symbolic('complex_scalar_Dirac_beta_sum', sy.Rational(1, 3) * 2 + sy.Rational(4, 3) - 2)
    # Real part j of a holomorphic Kahler transformation in Planck units.
    j = sy.symbols('j', real=True)
    transformed_uv = uv.subs({k: k + 2*j, f: f - 2*j/(8*pi**2)}, simultaneous=True)
    transformed_mlog = mlog.subs({k: k + 2*j, ell: ell-j}, simultaneous=True)
    transformed_low = low.subs({f: f - 2*j/(8*pi**2), ell: ell-j}, simultaneous=True)
    symbolic('Kahler_invariant_UV_physical_g', transformed_uv-uv)
    symbolic('Kahler_invariant_physical_mass', transformed_mlog-mlog)
    symbolic('Kahler_invariant_low_energy_g', transformed_low-low)
    # Reparameterize charged fields Q'_+-=exp(u_+-)Q_+-; u is Re(sum u_i).
    u = sy.symbols('u', real=True)
    symbolic('charged_redefinition_physical_mass',
             mlog.subs({ell: ell-u, z: z-2*u}, simultaneous=True)-mlog)
    symbolic('charged_redefinition_UV_g_with_Konishi_anomaly',
             uv.subs({f: f-u/(4*pi**2), z:z-2*u}, simultaneous=True)-uv)
    symbolic('charged_redefinition_low_g',
             low.subs({f:f-u/(4*pi**2), ell:ell-u}, simultaneous=True)-low)
    # An imposed constant physical UV g needs a compensating, generally
    # nonharmonic real f variation; it is not the fixed-holomorphic-f theory.
    f_phys = sy.symbols('f_phys', real=True)
    fixed_phys_low = (uv+threshold).subs(f, f_phys-k/(8*pi**2)+z/(8*pi**2))
    symbolic('fixed_physical_UV_recovers_mass_only_threshold',
             fixed_phys_low-(f_phys+threshold))
    # Full mass determinant leaves a nonholomorphic soft remainder.
    aa, bb, ls = sy.symbols('a b ls', real=True)
    weighted = sy.Rational(1, 3)*(ls+sy.log(1+aa+bb)) + sy.Rational(1, 3)*(ls+sy.log(1+aa-bb)) + sy.Rational(4, 3)*ls
    symbolic('split_threshold_log_separation',
             sy.expand(weighted)-2*ls-(sy.log(1+aa+bb)+sy.log(1+aa-bb))/3)
    delta_series = sy.series(sy.log((1+aa)**2-bb**2), bb, 0, 3).removeO()
    symbolic('split_threshold_even_in_B',
             sy.log((1+aa)**2-bb**2)-sy.log((1+aa)**2-(-bb)**2))
    symbolic('split_zero_soft_limit', (sy.log((1+aa)**2-bb**2)).subs({aa:0,bb:0}))
    symbolic('leading_common_soft_coefficient', sy.diff(sy.log((1+aa)**2-bb**2), aa).subs({aa:0,bb:0})-2)
    symbolic('leading_B2_coefficient', sy.diff(sy.log((1+aa)**2-bb**2),bb,2).subs({aa:0,bb:0})/2+1)

    raw = json.loads(INPUTS.read_text())
    inp = {key: mp.mpf(str(value)) for key, value in raw.items()
           if isinstance(value, (int, float))}
    eps, ci = inp['epsilon'], inp['chi_i']
    alpha = inp['alpha_reference']
    mi = inp['holomorphic_mass_i_eV']
    fs, mp2 = inp['F_squared_eV2'], inp['Mpl_eV']**2
    rows = list(csv.DictReader(TRAJECTORIES.open()))

    def spectrum(model, chi):
        g = (1+(ci/chi)**2)/2
        mh2 = mi**2*g
        rr = -ci**2/(2*chi**3*g)
        U = inp['potential_unit_eV4']*inp['W_i']*mp.exp(-2*(chi-ci))
        if model == 'canonical':
            kval = eps*chi**2/2
            aval = 1-eps*chi/2
            dj = aval**2-eps
            bracket = aval*(rr+eps*chi/2)+eps/2
        else:
            kval = mp.mpf(0)
            dj = 1-eps
            bracket = rr+(eps/2 if model == 'shift' else 0)
        s = mp.exp(kval)*mh2
        d = mp.exp(kval)*U*dj/mp2 if model != 'global' else mp.mpf(0)
        h = mp.exp(kval)*mp.sqrt(2*U*mh2/fs)*abs(bracket)
        remainder = mp.log1p(d/s+h/s)+mp.log1p(d/s-h/s)
        return dict(k=kval, mh2=mh2, s=s, d=d, h=h, remainder=remainder)

    samples, histories = [], []
    for model in ('global', 'canonical', 'shift'):
        axis = [r for r in rows if r['case'] == model+'_axis']
        final = spectrum(model, mp.mpf(axis[-1]['chi']))
        initial = spectrum(model, mp.mpf(axis[0]['chi']))
        dlhol = mp.log(initial['mh2']/final['mh2'])
        dr = initial['remainder']-final['remainder']
        q = alpha/(2*mp.pi)*dlhol+alpha/(12*mp.pi)*dr
        qlocal = alpha/(2*mp.pi)*mp.log(initial['s']/final['s'])+alpha/(12*mp.pi)*dr
        histories.append({'model':model, 'z_initial':axis[0]['z'],
                          'alpha_ppm_fixed_holomorphic_f':mp.nstr(1e6*q/(1-q),25),
                          'alpha_ppm_mass_only_fixed_physical_UV':mp.nstr(1e6*qlocal/(1-qlocal),25),
                          'soft_delta_g_inverse2_initial_vs_today':mp.nstr(-dr/(48*mp.pi**2),25),
                          'soft_delta_alpha_fraction_at_first_order':mp.nstr(alpha*dr/(12*mp.pi),25)})
        for index in (0,1,256,512,768,1024):
            row = axis[index]
            p = spectrum(model,mp.mpf(row['chi']))
            ln_s = mp.log(p['s']/final['s'])
            dk = p['k']-final['k']
            # Direct three-particle weighted determinant, no truncated series.
            direct = ((mp.log((p['s']+p['d']+p['h'])/(final['s']+final['d']+final['h']))
                       +mp.log((p['s']+p['d']-p['h'])/(final['s']+final['d']-final['h'])))/3
                      +mp.mpf(4)/3*ln_s)
            split_delta = p['remainder']-final['remainder']
            stable = 2*ln_s+split_delta/3
            numerical(f'{model}_{index}_weighted_direct_vs_stable', direct-stable, mp.mpf('1e-175'))
            direct_low = dk/(8*mp.pi**2)-direct/(16*mp.pi**2)
            hol_low = -mp.log(p['mh2']/final['mh2'])/(8*mp.pi**2)-split_delta/(48*mp.pi**2)
            numerical(f'{model}_{index}_KL_vs_holomorphic', direct_low-hol_low,mp.mpf('1e-175'))
            checks.append({'name':f'{model}_{index}_positive_scalar_spectrum',
                           'passed':bool(p['s']+p['d']-p['h']>0)})
            # Reconstruct source mass ratios without the primary module.
            for key,val in [('d_over_s',p['d']/p['s']),('h_over_s',p['h']/p['s'])]:
                ref = mp.mpf(row[key])
                err = abs(ref-val)/max(abs(val),mp.mpf('1e-150'))
                numerical(f'{model}_{index}_{key}_archive_consistency',err,mp.mpf('1e-11'))
            leading = 2*p['d']/p['s']-(p['h']/p['s'])**2
            err = abs(p['remainder']-leading)/max(abs(leading),mp.mpf('1e-150'))
            numerical(f'{model}_{index}_tiny_soft_series',err,mp.mpf('1e-84'))
            samples.append({'model':model,'index':index,'chi':row['chi'],
                            'd_over_s':mp.nstr(p['d']/p['s'],25),
                            'h_over_s':mp.nstr(p['h']/p['s'],25),
                            'soft_log_remainder':mp.nstr(p['remainder'],25),
                            'delta_g_inverse2_KL':mp.nstr(hol_low,25)})

    result = {'created_utc':datetime.now(timezone.utc).isoformat(),
              'method':'Independent KL algebra; 180 decimal digit direct three-mass logarithms; archived tree trajectories only',
              'scope':'Minimal unit-charge U(1) pair, fixed holomorphic UV f and fixed physical cutoff, no surviving charged matter; local one-loop adiabatic matching',
              'precision_digits':mp.mp.dps,'python':platform.python_version(),
              'mpmath':mp.__version__,'sympy':sy.__version__,
              'materials':{str(p.relative_to(REPO)):sha(p) for p in (INPUTS,TRAJECTORIES,Path(__file__).resolve(),SOURCE)},
              'primary_source':'https://arxiv.org/abs/hep-th/9402005v2',
              'primary_equations':['(2.21)','(3.4)','(3.7)','(3.13)','(3.19)','(3.21)--(3.29)'],
              'read_new_gauge_matching_primary_code':False,
              'checks':checks,'passed':sum(c['passed'] for c in checks),'total':len(checks),
              'histories':histories,'samples':samples}
    output = Path(__file__).with_suffix('.json')
    output.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({'passed':result['passed'],'total':result['total'],'histories':histories},indent=2))
    if result['passed'] != result['total']:
        raise SystemExit('Independent KL matching checks failed')


if __name__ == '__main__':
    main()
