#!/usr/bin/env python3
"""Independent WLS cross-check; does not import the recovery implementation.

Uses rational integer-power transfer formulae and covariance-based GLS normal
equations, rather than the production log1p transfer and nuisance QR projection.
All synthetic responses are created here; no observed alpha mean is read.
"""
from pathlib import Path
import csv
import hashlib
import json
import numpy as np
from scipy.linalg import cho_factor, cho_solve
from scipy.stats import chi2, norm

ROOT = Path(__file__).resolve().parents[1]
summary = json.loads((ROOT/'results/summary.json').read_text())
rows = list(csv.DictReader((ROOT/'results/observation_design.csv').open()))
recovery = list(csv.DictReader((ROOT/'results/recovery.csv').open()))
shapes = list(csv.DictReader((ROOT/'results/shape_geometry.csv').open()))
YEAR = 365.25*86400
T0 = summary['background']['t0_years']
L0 = summary['background']['L0']
checks = []
details = []


def check(name, ok, **info):
    checks.append(dict(name=name, passed=bool(ok), **info))


def trans(z, s, l0=L0):
    # The common prefactor in age(z)/age(0) cancels.
    z = np.asarray(z, dtype=np.longdouble)
    om = np.longdouble('.315')
    a = np.sqrt((1-om)/om)
    ell = np.log(np.arcsinh(a/(1+z)**np.longdouble('1.5'))/np.arcsinh(a))
    ll = np.longdouble(l0)
    if s == 1:
        v = T0*ll*ell/(ll+ell)
    elif s == 2:
        v = T0*ll*ell*(2*ll+ell)/(2*(ll+ell)**2)
    elif s == 3:
        v = T0*ll*ell*(3*ll**2+3*ll*ell+ell**2)/(3*(ll+ell)**3)
    return np.asarray(v, dtype=float)


def construct(spec):
    rr = [r for r in rows if r['dataset'] == ('ESPRESSO2022' if spec['name'].startswith('espresso') else 'King2012')]
    z = np.array([float(r['z']) for r in rr])
    sig = np.array([float(r['sigma_ppm']) for r in rr])
    group = np.array([r['instrument'] for r in rr])
    qso = np.array([r['qso'] for r in rr])
    nq = len(z)
    clk = spec['n_total'] > nq
    xx = np.column_stack([1e-13*trans(z, s) for s in (1,2,3)])
    nuisance = np.column_stack([group == g for g in spec['nuisance']]).astype(float) if spec['nuisance'] else np.empty((nq,0))
    if clk:
        sig = np.r_[sig, summary['current_clock_sigma_d']]
        xx = np.vstack((xx,np.ones(3)))
        nuisance = np.vstack((nuisance,np.zeros(nuisance.shape[1])))
    corr = np.eye(len(sig))
    rho = spec['assumed_same_qso_rho']
    corr[:nq,:nq] = (1-rho)*np.eye(nq)+rho*(qso[:,None] == qso[None,:])
    covariance = sig[:,None]*sig[None,:]*corr
    factor = cho_factor(covariance,lower=True)
    ci = lambda a: cho_solve(factor,a)
    wls = []
    for j in range(3):
        design = np.column_stack((xx[:,j],nuisance))
        scale = np.linalg.norm(design,axis=0)
        scaled = design/scale
        cov_scaled = np.linalg.inv(scaled.T@ci(scaled))
        cov = cov_scaled/scale[:,None]/scale[None,:]
        operator = (cov_scaled@ci(scaled).T)/scale[:,None]
        wls.append(dict(sigma=np.sqrt(cov[0,0]),op=operator[0],design=design,operator=operator))
    return xx, nuisance, covariance, ci, wls


rng = np.random.default_rng(2026092271)
for spec in summary['designs']:
    name = spec['name']
    xx,b,cov,ci,wls = construct(spec)
    sigs = np.array([w['sigma'] for w in wls])
    expected = np.array([spec['sigma_d'][str(s)] for s in (1,2,3)])
    err = np.max(abs(sigs/expected-1))
    check(name+'_GLS_sigma',err < 2e-10,max_relative_error=float(err))
    operators = np.stack([w['op']/w['sigma'] for w in wls])
    corr = operators@cov@operators.T
    check(name+'_score_covariance',np.max(abs(corr-np.array(spec['score_correlation'])))<2e-12)
    if spec['n_total'] > spec['n_astronomical']:
        frac = sigs[1]**2/summary['current_clock_sigma_d']**2
        check(name+'_clock_information_fraction',abs(frac-spec['clock_information_fraction'])<1e-12)
    for row in [r for r in recovery if r['design']==name]:
        truth = 2 if row['true_s']=='0' else int(row['true_s'])
        amplitude = float(row['d_injected'])
        mu = (wls[1]['op']@(amplitude*xx[:,truth-1]))/sigs[1]
        cutoff = np.sqrt(spec['fixed_s2_5pct_threshold'])
        power = norm.sf(cutoff-mu)+norm.cdf(-cutoff-mu)
        check(name+f"_power_s{truth}_d{amplitude}",abs(power-float(row['fixed_s2_analytic_power']))<1e-11)
    for row in [r for r in shapes if r['design']==name]:
        truth,fit = int(row['true_s'])-1,int(row['fit_s'])-1
        signal = 2.5*xx[:,truth]
        beta = wls[fit]['operator']@signal
        residual = signal-wls[fit]['design']@beta
        lam = float(residual@ci(residual))
        target = float(row['delta_chi2_at_d2p5'])
        check(name+f'_shape_{truth+1}_to_{fit+1}',abs(lam/target-1)<2e-7,
              relative_error=float(abs(lam/target-1)))
    # Full observation simulations use a fresh seed; compare to exact Gaussian laws.
    nmc = 16000
    scores = []
    chol = np.linalg.cholesky(cov)
    for k in range(0,nmc,1000):
        noise = chol@rng.standard_normal((len(cov), min(1000,nmc-k)))
        scores.append(operators@noise)
    scores = np.hstack(scores)
    fixed = float(np.mean(scores[1]**2 > spec['fixed_s2_5pct_threshold']))
    scan = float(np.mean(np.max(scores**2,axis=0)>spec['scan_5pct_threshold']))
    se = np.sqrt(.05*.95/nmc)
    combined = np.sqrt(.05*.95*(1/nmc+1/summary['n_calibration']))
    check(name+'_full_observation_null',abs(fixed-.05)<5*se,
          false_positive=fixed,mc_se=float(se))
    check(name+'_full_observation_scan',abs(scan-.05)<5*combined,
          false_positive=scan,calibration_and_evaluation_mc_se=float(combined))
    mainzero = next(r for r in recovery if r['design']==name and r['true_s']=='0')
    published_scan = float(mainzero['scan_detection'])
    published_scan_se = np.sqrt(published_scan*(1-published_scan)/summary['n_evaluation'])
    details.append(dict(design=name,independent_sigma_d=sigs.tolist(),
         independent_full_observation_null=fixed,independent_full_observation_scan=scan,
         stored_scan_false_positive=published_scan,stored_scan_conditional_mc_se=float(published_scan_se),
         combined_scan_mc_se_at_nominal=float(np.sqrt(.05*.95*(1/summary['n_evaluation']+1/summary['n_calibration'])))))

result=dict(status='PASS' if all(c['passed'] for c in checks) else 'FAIL',
    method='Independent rational transfer, covariance GLS normal equations, full-observation noise',
    reviewed_main_script_sha256=summary['script_sha256'],
    check_script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    checks=checks,design_details=details)
(ROOT/'reports/independent_statistical_checks.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(dict(status=result['status'],checks=len(checks),failed=[c for c in checks if not c['passed']],design_details=details),indent=2))
if result['status'] != 'PASS':
    raise SystemExit(1)
