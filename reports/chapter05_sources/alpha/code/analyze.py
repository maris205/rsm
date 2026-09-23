#!/usr/bin/env python3
"""Reproduce all summary-data fits. No network calls, stochastic fits or clipping.

Internally d is the current fractional alpha drift in 1e-19 / Julian year.
Quasar response and errors are ppm. Gamma is always dimensionless.
"""
from pathlib import Path
import csv
import hashlib
import json
import platform
import sys

import numpy as np
import scipy
from scipy.integrate import quad
from scipy.stats import chi2 as chi2_dist
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
YEAR = 365.25 * 86400.0
MPC_KM = 3.0856775814913673e19
TP = 5.391247e-44  # fixed reference scale in seconds
H0, OM = 67.4, 0.315
MODELS = ['log2', 'log1', 'log3', 'log_time', 'linear_time', 'exponential', 'inverse_time']
LABELS = {'log2': r'$1/\ln^2(t/t_*)$', 'log1': r'$1/\ln(t/t_*)$',
          'log3': r'$1/\ln^3(t/t_*)$', 'log_time': r'$\ln(t/t_0)$',
          'linear_time': r'$t/t_0$', 'exponential': r'$\exp(-t/t_0)$',
          'inverse_time': r'$t_0/t$'}


def age(z, h0=H0, om=OM):
    """Exact cosmic age (years) in the declared flat matter+Lambda background."""
    z = np.asarray(z, dtype=float)
    return 2 / (3 * (h0 / MPC_KM * YEAR) * np.sqrt(1-om)) * np.arcsinh(
        np.sqrt((1-om)/om) / (1+z)**1.5)


def transfer(z, model='log2', h0=H0, om=OM, tstar=TP):
    """y(z)/D0 in years, evaluated stably even near z=0."""
    t0 = age(0, h0, om)
    ratio = age(z, h0, om)/t0
    u = np.log(ratio)
    if model in ('log1', 'log2', 'log3'):
        p = int(model[-1])
        l0 = np.log(t0 * YEAR / tstar)
        if np.any(l0+u <= 0):
            raise ValueError('The logarithmic model requires t > tstar.')
        return -t0 * l0/p * np.expm1(-p*np.log1p(u/l0))
    if model == 'log_time':
        return t0*u
    if model == 'linear_time':
        return t0*(ratio-1)
    if model == 'exponential':
        return -t0*np.expm1(1-ratio)
    if model == 'inverse_time':
        return -t0*np.expm1(-u)
    raise ValueError(model)


def gamma_per_d(h0=H0, om=OM, tstar=TP):
    t0 = age(0, h0, om)
    return -0.5e-19*t0*np.log(t0*YEAR/tstar)**3


def read_king(include_outliers=False, budget='published'):
    path = ROOT/'data/raw/King2012_tablea1.dat'
    rows = [line.split() for line in path.read_text().splitlines() if line.strip()]
    if any(len(row) != 9 for row in rows):
        raise ValueError('Expected 9 columns in the CDS King table.')
    rows = [r for r in rows if include_outliers or int(r[8]) == 0]
    z = np.array([float(r[3]) for r in rows])
    y = np.array([float(r[4]) for r in rows])*10  # 1e-5 -> ppm
    stat = np.array([float(r[5]) for r in rows])*10
    extra = np.array([{1:0, 2:17.43, 3:9.05}[int(r[7])] for r in rows])
    if budget == 'statistical':
        extra *= 0
    elif budget == 'dipole':
        extra = np.array([{1:0, 2:16.30, 3:9.05}[int(r[7])] for r in rows])
    elif budget != 'published':
        raise ValueError(budget)
    return dict(z=z, y=y, sigma=np.hypot(stat, extra),
                sample=np.array([r[6] for r in rows]),
                id=np.array([int(r[0]) for r in rows]),
                qso=np.array([r[1] for r in rows]))


def subset(data, mask):
    return {key: value[mask] for key,value in data.items()}


def solve(design, y, sigma):
    if design.shape[1] == 0:
        return np.array([]), np.empty((0,0)), float(np.sum((y/sigma)**2))
    a = design / sigma[:,None]
    b = y/sigma
    u,s,vt = np.linalg.svd(a, full_matrices=False)
    if s[-1] <= s[0]*1e-12:
        raise ValueError('Rank-deficient design matrix.')
    beta = vt.T @ ((u.T @ b)/s)
    cov = (vt.T/s**2) @ vt
    return beta, cov, float(np.sum(((y-design@beta)/sigma)**2))


def fit(data, model='log2', clock=None, offsets=False, h0=H0, om=OM, tstar=TP):
    z,y,sigma = data['z'],data['y'],data['sigma']
    nuisance = np.column_stack([data['sample']==x for x in sorted(set(data['sample']))]).astype(float) if offsets else np.empty((len(z),0))
    design = np.column_stack([1e-13*transfer(z,model,h0,om,tstar),nuisance])
    if clock is not None:
        design = np.vstack([design, np.r_[1., np.zeros(nuisance.shape[1])]])
        nuisance = np.vstack([nuisance, np.zeros((1,nuisance.shape[1]))])
        y = np.r_[y,clock[0]]
        sigma = np.r_[sigma,clock[1]]
    beta,cov,c2 = solve(design,y,sigma)
    b0,cov0,c20 = solve(nuisance,y,sigma)
    k = design.shape[1]
    gain = c20-c2
    r = dict(model=model,n=len(y),k=k,offsets=offsets,d=beta[0],sigma_d=np.sqrt(cov[0,0]),
             chi2=c2,chi2_null=c20,delta_chi2=gain,delta_aic=2-gain,
             delta_bic=np.log(len(y))-gain,dof=len(y)-k,
             nuisance_ppm=beta[1:].tolist(),nuisance_null_ppm=b0.tolist(),
             covariance=cov.tolist(),h0=h0,omega_m=om,tstar_s=tstar)
    if len(y)>k:
        r['nominal_goodness_p'] = chi2_dist.sf(c2,len(y)-k)
    if model == 'log2':
        r['gamma'] = beta[0]*gamma_per_d(h0,om,tstar)
        r['sigma_gamma'] = np.sqrt(cov[0,0])*abs(gamma_per_d(h0,om,tstar))
    return r


def load_csv(path):
    with path.open() as f:
        return list(csv.DictReader(f))


def verified_clock():
    rows = load_csv(ROOT/'data/clock_constraints.csv')
    r = next(row for row in rows if row['id']=='Filzinger2023')
    return (float(r['mean_per_year'])/1e-19, float(r['sigma_per_year'])/1e-19)


def modern_optical():
    rows = load_csv(ROOT/'data/modern_measurements.csv')
    return dict(z=np.array([float(r['absorber_z']) for r in rows]),
                y=np.array([float(r['delta_alpha_ppm']) for r in rows]),
                sigma=np.array([np.hypot(float(r['sigma_stat_ppm']),float(r['sigma_sys_ppm'])) for r in rows]),
                sample=np.array([r['measurement_id'] for r in rows]))


def verification():
    t0=age(0)
    age_diffs=[]
    for z in [0.,.1,1.,2.,4.2]:
        numeric=quad(lambda a: np.sqrt(a)/np.sqrt(OM+(1-OM)*a**3),0,1/(1+z),epsabs=1e-12)[0]/(H0/MPC_KM*YEAR)
        age_diffs.append(abs(numeric-age(z))/age(z))
    # Independent finite difference in physical cosmic time for Gamma=1.
    dt=1e5
    f=lambda t:np.log(t*YEAR/TP)**-2
    fd=(f(t0+dt)-f(t0-dt))/(2*dt)
    analytic=-2/(t0*np.log(t0*YEAR/TP)**3)
    z=np.array([.2,1.,4.])
    direct=(np.log(age(z)*YEAR/TP)**-2-f(t0))*1e6
    normalized=1e-13*transfer(z)/gamma_per_d()
    data=read_king()
    r=fit(data)
    x=1e-13*transfer(data['z'])
    w=1/data['sigma']**2
    exact_d=np.sum(w*x*data['y'])/np.sum(w*x*x)
    checks=dict(catalog_rows=295,retained_rows=len(data['z']),
                retained_keck=int(np.sum(data['sample']=='Keck')),
                retained_vlt=int(np.sum(data['sample']=='VLT')),
                age_quadrature_max_relative_error=max(age_diffs),
                derivative_relative_error=abs(fd/analytic-1),
                gamma_ppm_normalization_max_relative_error=float(np.max(abs(direct/normalized-1))),
                svd_vs_scalar_d_difference=float(abs(exact_d-r['d'])),
                model_zero_at_z0=all(float(transfer(0,m))==0 for m in MODELS))
    assert len(data['z'])==293
    assert max(age_diffs)<1e-9
    assert checks['derivative_relative_error']<1e-6
    assert checks['gamma_ppm_normalization_max_relative_error']<1e-10
    assert checks['svd_vs_scalar_d_difference']<1e-8
    assert checks['model_zero_at_z0']
    return checks


def dump_json(path,obj):
    path.write_text(json.dumps(obj,indent=2,default=lambda v: v.item() if isinstance(v,np.generic) else str(v))+'\n')


def save_tables(results):
    fields=['dataset','model','n','k','d','sigma_d','gamma','sigma_gamma','chi2','chi2_null','delta_chi2','delta_aic','delta_bic']
    with (ROOT/'results/fits.csv').open('w') as f:
        writer=csv.DictWriter(f,fields,extrasaction='ignore')
        writer.writeheader()
        for dataset,fits in results['fits'].items():
            for r in fits:writer.writerow(dict(dataset=dataset,**r))
    rows=[]
    names={'king':'King293','king_clock':'King293 + clock','king_offsets':'King293, offsets','king_offsets_clock':'King293 + clock, offsets',
           'keck':'Keck140','vlt':'VLT153','optical':'ESPRESSO','optical_clock':'ESPRESSO + clock'}
    for key,label in names.items():
        r=results['fits'][key][0]
        aic = f"{r['delta_aic']:+.3f}" if r['n'] > 10 else '---'
        rows.append(f"{label} & {r['n']} & ${r['gamma']:+.5f}\\pm{r['sigma_gamma']:.5f}$ & {r['chi2']:.3f} & {r['delta_chi2']:.3f} & {aic} \\\\")
    (ROOT/'paper/results_primary.tex').write_text(
        '\\begingroup\n\\setlength{\\tabcolsep}{5pt}\n\\begin{tabular}{lrrrrr}\n\\toprule\nData & $N$ & $\\Gamma$ & $\\chi^2$ & $\\Delta\\chi^2$ & $\\Delta\\mathrm{AIC}$\\\\\n\\midrule\n'
        +'\n'.join(rows)+'\n\\bottomrule\n\\end{tabular}\n\\endgroup\n')
    rows=[]
    for r,rj in zip(results['fits']['king'],results['fits']['king_clock']):
        rows.append(f"{LABELS[r['model']]} & {r['chi2']:.3f} & {r['delta_aic']:+.3f} & {rj['chi2']:.3f} & {rj['delta_aic']:+.3f} \\\\")
    (ROOT/'paper/results_shapes.tex').write_text(
        '\\begingroup\n\\setlength{\\tabcolsep}{6pt}\n\\begin{tabular}{lrrrr}\n\\toprule\nHistory & $\\chi^2_{\\rm QSO}$ & $\\Delta\\mathrm{AIC}_{\\rm QSO}$ & $\\chi^2_{\\rm QSO+clk}$ & $\\Delta\\mathrm{AIC}_{\\rm QSO+clk}$\\\\\n\\midrule\n'
        +'\n'.join(rows)+'\n\\bottomrule\n\\end{tabular}\n\\endgroup\n')


def figures(results):
    plt.rcParams.update({'font.family':'serif','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.bbox':'tight'})
    data=read_king()
    optical=modern_optical()
    dclock,sclock=results['clock']
    z=np.linspace(0,4.2,401)
    colors={'Keck':'#0072B2','VLT':'#D55E00'}
    fig,axs=plt.subplots(1,2,figsize=(11,4.2))
    for telescope,c in colors.items():
        m=data['sample']==telescope
        axs[0].errorbar(data['z'][m],data['y'][m],yerr=data['sigma'][m],fmt='.',ms=3,lw=.55,alpha=.4,color=c,label=telescope)
    r=results['fits']['king'][0]
    axs[0].plot(z,1e-13*transfer(z)*r['d'],color='black',label='Joint-telescope log-squared fit')
    axs[0].axhline(0,color='gray',lw=.8,ls='--')
    axs[0].set(xlabel='Absorber redshift',ylabel=r'$\Delta\alpha/\alpha_0$ (ppm)',title='(a) Legacy King sample (293 rows)')
    axs[0].legend(fontsize=8)
    prediction=1e-13*transfer(z)*dclock
    band=1.96e-13*np.abs(transfer(z))*sclock
    axs[1].fill_between(z,prediction-band,prediction+band,color='#009E73',alpha=.25,label='Clock transfer: 95% pointwise interval')
    axs[1].plot(z,prediction,color='#009E73',label='Clock central value (not a detection)')
    axs[1].axhline(0,color='gray',lw=.8,ls='--')
    axs[1].set(xlabel='Redshift',ylabel=r'$\Delta\alpha/\alpha_0$ (ppm)',title='(b) Homogeneous log-squared prediction')
    axs[1].legend(fontsize=8,loc='lower left')
    fig.tight_layout();fig.savefig(ROOT/'figures/alpha_history.pdf');fig.savefig(ROOT/'figures/alpha_history.png',dpi=160);plt.close(fig)

    fig,axs=plt.subplots(1,2,figsize=(11,4.2))
    keys=['king','keck','vlt','king_offsets','optical']
    labs=['King293','Keck140','VLT153','King293 + offsets','ESPRESSO (1 absorber)']
    for i,key in enumerate(keys):
        r=results['fits'][key][0]
        axs[0].errorbar(r['gamma'],i,xerr=r['sigma_gamma'],fmt='o',color='#0072B2',capsize=3)
    gc=dclock*gamma_per_d();gs=sclock*abs(gamma_per_d())
    axs[0].axvspan(gc-1.96*gs,gc+1.96*gs,color='#009E73',alpha=.3,label='Clock 95% interval')
    axs[0].axvline(0,color='gray',ls='--',lw=.8)
    axs[0].set(yticks=range(len(labs)),yticklabels=labs,xlabel=r'Dimensionless coupling $\Gamma$',title='(a) Separate summary-data fits')
    axs[0].invert_yaxis();axs[0].legend(fontsize=8)
    for i,key in enumerate(['king_clock','king_offsets_clock','optical_clock']):
        r=results['fits'][key][0]
        axs[1].errorbar(r['gamma'],i,xerr=1.96*r['sigma_gamma'],fmt='o',capsize=4,color='#009E73')
    axs[1].axvline(0,color='gray',ls='--',lw=.8)
    axs[1].set(yticks=range(3),yticklabels=['King + clock','King + offsets + clock','ESPRESSO + clock'],xlabel=r'$\Gamma$ (95% interval)',title='(b) Clock-dominated joint fits')
    axs[1].invert_yaxis();fig.tight_layout();fig.savefig(ROOT/'figures/coupling_constraints.pdf');fig.savefig(ROOT/'figures/coupling_constraints.png',dpi=160);plt.close(fig)

    fig,axs=plt.subplots(1,2,figsize=(11,4.2))
    zz=z[1:]
    for model in ['log1','log2','log3','log_time','linear_time','exponential','inverse_time']:
        axs[0].plot(z,-transfer(z,model)/1e9,label=LABELS[model],lw=2 if model=='log2' else 1)
    axs[0].set(xlabel='Redshift',ylabel=r'$-\mathcal{T}(z)$ (Gyr)',title='(a) Equal present-drift normalization')
    axs[0].legend(fontsize=8,ncol=2)
    for model in ['log1','log3','log_time']:
        axs[1].plot(zz,100*(transfer(zz,model)/transfer(zz,'log2')-1),label=LABELS[model])
    axs[1].axhline(0,color='gray',lw=.8,ls='--')
    axs[1].set(xlabel='Redshift',ylabel='Relative difference from p = 2 (%)',title='(b) Inverse-log exponent degeneracy')
    axs[1].legend(fontsize=9);fig.tight_layout();fig.savefig(ROOT/'figures/shape_comparison.pdf');fig.savefig(ROOT/'figures/shape_comparison.png',dpi=160);plt.close(fig)


def main():
    for folder in ['results','figures','paper']:(ROOT/folder).mkdir(exist_ok=True)
    clock=verified_clock()
    king=read_king()
    optical=modern_optical()
    results=dict(background={'H0_km_s_Mpc':H0,'Omega_m':OM,'t0_years':float(age(0)),
                            'tstar_s':TP,'L0':float(np.log(age(0)*YEAR/TP)),
                            'gamma_per_d':float(gamma_per_d()),'year_s':YEAR},clock=clock,fits={})
    for key,data,c,offsets in [('king',king,None,False),('king_clock',king,clock,False),
          ('king_offsets',king,None,True),('king_offsets_clock',king,clock,True),
          ('keck',subset(king,king['sample']=='Keck'),None,False),
          ('vlt',subset(king,king['sample']=='VLT'),None,False),
          ('optical',optical,None,False),('optical_clock',optical,clock,False),
          ('king_statistical',read_king(budget='statistical'),None,False),
          ('king_all295',read_king(include_outliers=True),None,False),
          ('king_dipole_budget',read_king(budget='dipole'),None,False)]:
        results['fits'][key]=[fit(data,model,c,offsets) for model in MODELS]
    results['background_sensitivity']=[fit(king,'log2',clock,h0=h,om=o) for h,o in [(67.4,.315),(73.,.3),(67.4,.30),(67.4,.33)]]
    results['scale_sensitivity']=[dict(tstar_s=ts,L0=float(np.log(age(0)*YEAR/ts)),T_z4_yr=float(transfer(4,tstar=ts))) for ts in [TP,1.,YEAR]]
    results['clock_transfer']=[dict(z=z,T_years=float(transfer(z)),central_ppm=float(1e-13*transfer(z)*clock[0]),
                                  sigma_ppm=float(1e-13*abs(transfer(z))*clock[1]),
                                  p1_relative=float(transfer(z,'log1')/transfer(z)-1),
                                  p3_relative=float(transfer(z,'log3')/transfer(z)-1)) for z in [1,2,4]]
    q=results['fits']['king'][0];k=results['fits']['keck'][0];v=results['fits']['vlt'][0]
    results['diagnostics']={'king_clock_parameter_difference_sigma':abs(q['d']-clock[0])/np.hypot(q['sigma_d'],clock[1]),
        'keck_vlt_parameter_difference_sigma':abs(k['d']-v['d'])/np.hypot(k['sigma_d'],v['sigma_d']),
        'clock_fraction_amplitude_information':(1/clock[1]**2)/(1/clock[1]**2+1/q['sigma_d']**2)}
    results['verification']=verification()
    dump_json(ROOT/'results/results.json',results)
    save_tables(results)
    figures(results)
    files=[p for folder in ['data','code'] for p in (ROOT/folder).rglob('*') if p.is_file() and '__pycache__' not in p.parts]
    dump_json(ROOT/'results/run_manifest.json',{'python':sys.version,'platform':platform.platform(),
             'numpy':np.__version__,'scipy':scipy.__version__,'matplotlib':matplotlib.__version__,
             'command':'python code/analyze.py','inputs':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(files)}})
    print(json.dumps({'verification':results['verification'],'primary':results['fits']['king'][0],
                      'joint':results['fits']['king_clock'][0],'diagnostics':results['diagnostics']},indent=2))


if __name__=='__main__':
    main()
