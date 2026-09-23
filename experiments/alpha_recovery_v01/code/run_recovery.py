#!/usr/bin/env python3
"""Conditional Gaussian injection/recovery; no measured responses are fitted.

The three sufficient scores are simulated by a QR factor, preserving their
near-singular correlation without clipping eigenvalues. All writes stay here.
"""
from pathlib import Path
import csv
import hashlib
import json
import platform

import numpy as np
import scipy
from scipy.linalg import solve_triangular
from scipy.optimize import brentq, minimize_scalar
from scipy.stats import norm, chi2

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / 'data/source'
YEAR = 365.25 * 86400.0
TP = 5.391247e-44
H0, OM = 67.4, .315
MPC_KM = 3.0856775814913673e19
SEED, NCAL, NMC = 20260922, 30000, 20000
EXPONENTS = (1, 2, 3)


def age(z):
    return 2 / (3 * H0 / MPC_KM * YEAR * np.sqrt(1 - OM)) * np.arcsinh(
        np.sqrt((1 - OM) / OM) / (1 + np.asarray(z)) ** 1.5)


T0 = float(age(0.0))
L0 = float(np.log(T0 * YEAR / TP))


def transfer(z, s, l0=L0):
    ell = np.log(age(z) / T0)
    if s <= 0 or np.any(l0 + ell <= 0):
        raise ValueError('Positive exponent and positive logarithm required')
    return -T0 * l0 / s * np.expm1(-s * np.log1p(ell / l0))


def read_csv(path):
    with path.open() as f:
        return list(csv.DictReader(f))


def inputs():
    rows = [x.split() for x in (SOURCE/'raw/King2012_tablea1.dat').read_text().splitlines() if x.strip()]
    retained = [r for r in rows if int(r[8]) == 0]
    # Deliberately do not read column 4 (the measured alpha response).
    k = dict(z=np.array([float(r[3]) for r in retained]),
             sigma=np.array([np.hypot(10*float(r[5]), {1:0, 2:17.43, 3:9.05}[int(r[7])]) for r in retained]),
             group=np.array([r[6] for r in retained]),
             qso=np.array([r[1] for r in retained]),
             ids=np.array([r[0] for r in retained]))
    mr = read_csv(SOURCE/'modern_measurements.csv')
    e = dict(z=np.array([float(r['absorber_z']) for r in mr]),
             sigma=np.array([np.hypot(float(r['sigma_stat_ppm']), float(r['sigma_sys_ppm'])) for r in mr]),
             group=np.array(['ESPRESSO']*len(mr)), qso=np.array([r['target'] for r in mr]),
             ids=np.array([r['measurement_id'] for r in mr]))
    cr = next(r for r in read_csv(SOURCE/'clock_constraints.csv') if r['id']=='Filzinger2023')
    cs = float(cr['sigma_per_year'])/1e-19
    assert len(k['z']) == 293 and np.isclose(cs, 2.5, rtol=1e-14)
    return k, e, cs


class Design:
    def __init__(self, name, data, clock_sigma=None, offsets=False, rho=0., error_reduction=1.):
        self.name, self.data, self.clock_sigma = name, data, clock_sigma
        self.offsets, self.rho = offsets, rho
        self.nq = len(data['z'])
        self.sigma = data['sigma']/error_reduction
        if clock_sigma is not None:
            self.sigma = np.r_[self.sigma, clock_sigma]
        self.n = len(self.sigma)
        self.corr = np.eye(self.n)
        same = data['qso'][:, None] == data['qso'][None, :]
        self.corr[:self.nq, :self.nq] = (1-rho)*np.eye(self.nq) + rho*same
        self.chol = np.linalg.cholesky(self.corr) if rho else None
        self.labels = sorted(set(data['group'])) if offsets else []
        b = np.column_stack([data['group']==g for g in self.labels]).astype(float) if offsets else np.empty((self.nq, 0))
        if clock_sigma is not None:
            b = np.vstack([b, np.zeros((1, b.shape[1]))])
        self.b = b
        bw = self.whiten(b)
        self.qb = np.linalg.qr(bw, mode='reduced')[0] if offsets else np.empty((self.n, 0))
        self.x = np.column_stack([self.raw_column(s) for s in EXPONENTS])
        self.r = self.residualize(self.whiten(self.x))
        self.length = np.linalg.norm(self.r, axis=0)
        self.unit = self.r/self.length
        # QR is stable even when the candidate curves almost coincide.
        self.q, self.rfactor = np.linalg.qr(self.unit, mode='reduced')

    def whiten(self, x):
        z = x/self.sigma if x.ndim == 1 else x/self.sigma[:, None]
        return solve_triangular(self.chol, z, lower=True) if self.chol is not None else z

    def residualize(self, x):
        return x-self.qb@(self.qb.T@x)

    def raw_column(self, s, l0=L0):
        x = 1e-13*transfer(self.data['z'], s, l0)
        return np.r_[x, 1.] if self.clock_sigma is not None else x

    def projected_column(self, s, l0=L0):
        return self.residualize(self.whiten(self.raw_column(s, l0)))

    def separation(self, true_s, fit_s, d, fit_l0=L0):
        v = d*self.projected_column(true_s)
        w = self.projected_column(fit_s, fit_l0)
        a = np.dot(w, v)/np.dot(w, w)
        res = v-a*w
        return float(np.dot(res, res)), float(a)

    def noise_scores(self, rng, count):
        return self.rfactor.T @ rng.standard_normal((self.rfactor.shape[0], count))


def json_write(path, value):
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2,
        default=lambda x: x.item() if isinstance(x, np.generic) else x.tolist())+'\n')


def csv_write(path, rows):
    with path.open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def normal_power(mu, threshold):
    cut = np.sqrt(threshold)
    return float(norm.sf(cut-mu)+norm.cdf(-cut-mu))


def run():
    out = ROOT/'results'
    out.mkdir(exist_ok=True, parents=True)
    king, espresso, clock_sigma = inputs()
    specs = [('king_only', king, None, False, 0.),
             ('king_only_offsets', king, None, True, 0.),
             ('king_clock', king, clock_sigma, False, 0.),
             ('king_clock_offsets', king, clock_sigma, True, 0.),
             ('espresso_clock', espresso, clock_sigma, False, 0.),
             ('king_only_offsets_rho025', king, None, True, .25),
             ('king_clock_offsets_rho025', king, clock_sigma, True, .25)]
    designs = [Design(*s) for s in specs]
    checks=[]
    def check(name, passed, **kw):
        checks.append(dict(name=name, passed=bool(passed), **kw))
    summary, recoveries, shapes, mismatches = [], [], [], []
    threshold = float(chi2.ppf(.95, 1))
    children = np.random.SeedSequence(SEED).spawn(len(designs)*3)
    for j, design in enumerate(designs):
        calibration = design.noise_scores(np.random.default_rng(children[3*j]), NCAL)
        scan_cut = float(np.quantile(np.max(calibration**2, axis=0), .95))
        noise = design.noise_scores(np.random.default_rng(children[3*j+1]), NMC)
        item=dict(name=design.name, n_astronomical=design.nq, n_total=design.n,
                  nuisance=design.labels, assumed_same_qso_rho=design.rho,
                  sigma_d={str(s):float(1/design.length[i]) for i,s in enumerate(EXPONENTS)},
                  sigma_D0_per_year=float(1e-19/design.length[1]),
                  clock_information_fraction=float((1/clock_sigma**2)/design.length[1]**2) if design.clock_sigma else 0.,
                  fixed_s2_5pct_threshold=threshold, scan_5pct_threshold=scan_cut,
                  score_correlation=design.unit.T@design.unit)
        summary.append(item)
        # Independent full-observation WLS: QR-score calculations must agree.
        rng = np.random.default_rng(children[3*j+2])
        for fit_idx, s in enumerate(EXPONENTS):
            xw = design.whiten(np.column_stack([design.x[:,fit_idx], design.b]))
            yw = rng.standard_normal((design.n,8))+2.5*design.whiten(design.x[:,1])[:,None]
            if design.labels:
                yw += design.whiten(design.b)@np.full((len(design.labels),8),3.0)
            beta = np.linalg.lstsq(xw,yw,rcond=None)[0][0]
            score_beta = design.unit[:,fit_idx]@yw/design.length[fit_idx]
            err = float(np.max(np.abs(beta-score_beta)))
            check(design.name+f'_full_WLS_s{s}',err < 1e-8, max_abs_d_error=err)
        check(design.name+'_QR_covariance', np.max(abs(design.rfactor.T@design.rfactor-design.unit.T@design.unit))<2e-14)
        # Zero truth once; nonzero signed amplitudes for each of the three laws.
        truths=[(0,0.)]+[(s,d) for s in EXPONENTS for d in (-12.5,-5.,-2.5,2.5,5.,12.5)]
        for true_s,d in truths:
            true_idx = 1 if true_s==0 else EXPONENTS.index(true_s)
            mean=d*(design.unit.T@design.r[:,true_idx])
            scores=noise+mean[:,None]
            gains=scores**2
            zerr=scores[true_idx]-d*design.length[true_idx]
            coverage=float(np.mean(abs(zerr)<=norm.ppf(.975)))
            fixed_rate=float(np.mean(gains[1]>threshold))
            fixed_expected=normal_power(mean[1], threshold)
            scan_rate=float(np.mean(np.max(gains,axis=0)>scan_cut))
            scores_aic=np.vstack([np.zeros(NMC),gains-2])
            winners=np.argmax(scores_aic,axis=0)
            forced=np.argmax(gains,axis=0)
            sorted_gain=np.sort(gains,axis=0)
            separated=sorted_gain[-1]-sorted_gain[-2]>=2
            row=dict(design=design.name,true_s=true_s,d_injected=d,
                     D0_injected_per_year=d*1e-19,
                     bias_d=float(np.mean(zerr)/design.length[true_idx]),
                     bias_in_standard_errors=float(np.mean(zerr)),
                     empirical_standardized_sd=float(np.std(zerr,ddof=1)),
                     coverage95=coverage,coverage_mc_se=float(np.sqrt(coverage*(1-coverage)/NMC)),
                     fixed_s2_detection=fixed_rate,fixed_s2_analytic_power=fixed_expected,
                     detection_mc_se=float(np.sqrt(fixed_rate*(1-fixed_rate)/NMC)),
                     scan_detection=scan_rate,
                     aic_constant=float(np.mean(winners==0)),
                     aic_s1=float(np.mean(winners==1)),aic_s2=float(np.mean(winners==2)),aic_s3=float(np.mean(winners==3)),
                     forced_s1=float(np.mean(forced==0)),forced_s2=float(np.mean(forced==1)),forced_s3=float(np.mean(forced==2)),
                     shape_gap_at_least_2=float(np.mean(separated)),
                     true_shape_wins_with_gap_2=float(np.mean(separated & (forced==true_idx))) if true_s else 0.)
            recoveries.append(row)
            se=np.sqrt(max(fixed_expected*(1-fixed_expected),1/NMC)/NMC)
            check(design.name+f'_power_s{true_s}_d{d}',abs(fixed_rate-fixed_expected)<6*se+1/NMC)
        # Coverage check once per distinct true direction, not once per amplitude.
        for i,s in enumerate(EXPONENTS):
            coverage=float(np.mean(abs(noise[i])<=norm.ppf(.975)))
            check(design.name+f'_coverage_s{s}',abs(coverage-.95)<6*np.sqrt(.95*.05/NMC))
        for true_s in EXPONENTS:
            v=design.projected_column(true_s)
            for fit_s in EXPONENTS:
                if fit_s==true_s:continue
                lam, amp=design.separation(true_s,fit_s,2.5)
                shapes.append(dict(design=design.name,true_s=true_s,fit_s=fit_s,
                    residual_angle_sine=float(np.sqrt(lam)/(2.5*np.linalg.norm(v))),
                    delta_chi2_at_d2p5=lam,delta_chi2_at_d12p5=25*lam,
                    best_fitted_d_at_d2p5=amp))
    # A deliberately misspecified covariance, under the same-qso rho=0.25 noise.
    for j, name in enumerate(['king_only_offsets','king_clock_offsets']):
        des=next(d for d in designs if d.name==name)
        corr=np.eye(des.n)
        same=king['qso'][:,None]==king['qso'][None,:]
        corr[:des.nq,:des.nq]=.75*np.eye(des.nq)+.25*same
        u=des.unit[:,1]
        inflation=float(np.sqrt(u@corr@u))
        noise=np.random.default_rng(SEED+900+j).normal(0,inflation,NMC)
        mismatches.append(dict(fitted_design=name,true_same_qso_rho=.25,
            fitted_same_qso_rho=0,standard_error_inflation=inflation,
            null_false_positive=float(np.mean(noise**2>threshold)),
            null_false_positive_analytic=normal_power(0,threshold/inflation**2),
            coverage95=float(np.mean(abs(noise)<norm.ppf(.975))),
            coverage95_analytic=float(2*norm.cdf(norm.ppf(.975)/inflation)-1)))
    # One absorber plus its own arbitrary offset carries no historical shape info.
    unident=Design('espresso_clock_free_offset',espresso,clock_sigma,True)
    null_history=unident.r[:unident.nq,:]
    check('single_absorber_free_offset_erases_history',np.max(abs(null_history))<1e-12)
    # Scale diagnostics use no measured means and make no discovery claim.
    main=next(d for d in designs if d.name=='king_clock_offsets')
    sensitivity=[]
    for l0 in (10.,30.,70.,L0,280.):
        for s in (1,3):
            truth=2.5*main.projected_column(2,l0)
            candidate=main.projected_column(s,l0)
            best=float(candidate@truth/(candidate@candidate))
            res=truth-best*candidate
            sensitivity.append(dict(true_s=2,fit_s=s,common_L0=l0,
                common_tstar_seconds=float(T0*YEAR*np.exp(-l0)),
                delta_chi2_at_d2p5=float(res@res),best_d=best))
    profiles=[]
    for s in (1,3):
        objective=lambda l0:main.separation(2,s,2.5,l0)[0]
        opt=minimize_scalar(objective,bounds=(6.,400.),method='bounded',options={'xatol':1e-9})
        profiles.append(dict(true_s=2,fit_s=s,true_L0=L0,bounds_L0=[6.,400.],
            best_fit_L0=float(opt.x),fixed_scale_delta_chi2=objective(L0),
            profiled_scale_delta_chi2=float(opt.fun),optimization_success=bool(opt.success)))
    # Ideal conditional precision requirements: telescope offsets retained.
    def ideal_distance(log10factor,d):
        des=Design('ideal',king,clock_sigma,True,error_reduction=10.**log10factor)
        return min(des.separation(2,s,d)[0] for s in (1,3))
    budget=[]
    for d in (2.5,12.5):
        for target in (1.,9.):
            root=brentq(lambda r:ideal_distance(r,d)-target,0,10,xtol=1e-8)
            factor=10.**root
            budget.append(dict(d_injected=d,target_expected_delta_chi2=target,
                astronomical_error_reduction=factor,
                median_king_error_ppm_after=float(np.median(king['sigma'])/factor),
                clock_sigma_d_unchanged=clock_sigma,
                achieved_min_delta_chi2=ideal_distance(root,d)))
    budget_curve=[dict(error_reduction=10.**x,min_delta_chi2_d2p5=ideal_distance(x,2.5),
                       min_delta_chi2_d12p5=ideal_distance(x,12.5)) for x in np.linspace(0,8,65)]
    design_rows=[]
    for label,data in [('King2012',king),('ESPRESSO2022',espresso)]:
        for i in range(len(data['z'])):
            design_rows.append(dict(dataset=label,id=data['ids'][i],qso=data['qso'][i],
                instrument=data['group'][i],z=data['z'][i],sigma_ppm=data['sigma'][i]))
    sources=[SOURCE/'raw/King2012_tablea1.dat',SOURCE/'raw/King2012_provenance.json',
             SOURCE/'clock_constraints.csv',SOURCE/'modern_measurements.csv',ROOT/'protocol.md']
    hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}
    result=dict(scope='Synthetic responses at existing observational designs; no measured alpha or clock means fitted.',
        protocol_sha256=hashes['protocol.md'],
        seed=SEED,n_calibration=NCAL,n_evaluation=NMC,background=dict(H0=H0,Omega_m=OM,tstar_s=TP,t0_years=T0,L0=L0),
        current_clock_sigma_d=clock_sigma,king_rows=len(king['z']),king_unique_qsos=len(set(king['qso'])),
        z_range=[float(min(king['z'])),float(max(king['z']))],
        designs=summary,covariance_misspecification=mismatches,
        reference_scale_sensitivity=sensitivity,reference_scale_profile=profiles,
        ideal_precision_budget=budget,ideal_precision_curve=budget_curve,
        single_absorber_free_offset=dict(sigma_d=(1/unident.length).tolist(),shape_columns_identical=bool(np.max(np.ptp(unident.r,axis=1))<1e-12)),
        versions=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__),source_hashes=hashes,
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        checks=checks,all_checks_passed=all(c['passed'] for c in checks))
    csv_write(out/'observation_design.csv',design_rows)
    csv_write(out/'recovery.csv',recoveries)
    csv_write(out/'shape_geometry.csv',shapes)
    json_write(out/'summary.json',result)
    print(json.dumps(dict(status='PASS' if result['all_checks_passed'] else 'FAIL',
        scenarios=len(recoveries),checks=len(checks),failed=[c for c in checks if not c['passed']],
        primary_design=summary[3],ideal_precision_budget=budget),indent=2,default=lambda x:x.tolist()))
    if not result['all_checks_passed']:
        raise SystemExit('Validation failed; results retained for diagnosis.')


if __name__=='__main__':
    run()
