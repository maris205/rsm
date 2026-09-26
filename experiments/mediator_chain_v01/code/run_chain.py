#!/usr/bin/env python3
"""Gapped-chain matching and selected finite-F static response.

No observational inference or claim of a fully rematched cosmological model.
"""
from pathlib import Path
from datetime import datetime, timezone
import csv
import hashlib
import json
import numpy as np
import scipy
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parents[1]
REPO = ROOT.parents[1]
OUT = ROOT / 'results'


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def dump(p, data):
    p.write_text(json.dumps(data, ensure_ascii=False, indent=2, allow_nan=False) + '\n')


def csvout(p, rows):
    with p.open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator='\n')
        w.writeheader(); w.writerows(rows)


def chain(q, n, lam):
    g00 = q**-2 * (-np.expm1(-2*n*np.log(q))) / (-np.expm1(-2*(n+1)*np.log(q)))
    tau = (q-q**-1)*np.exp(-n*np.log(q))/(-np.expm1(-2*n*np.log(q)))
    i = np.arange(n)
    col = np.exp((i-n-1)*np.log(q)) * (-np.expm1(-2*(i+1)*np.log(q))) / (-np.expm1(-2*(n+1)*np.log(q)))
    h2 = float(col @ col)
    eigen = q*q+1-2*q*np.cos(np.arange(1,n+1)*np.pi/(n+1))
    M = lam*np.sqrt(g00)
    return dict(q=q, n=n, Lambda_eV=lam, G00=g00, tau=tau, zeta=1., h2=h2,
                M_eV=M, min_mass_eV=M*np.sqrt(eigen[0]),
                max_mass_eV=M*np.sqrt(eigen[-1]), condition_H=eigen[-1]/eigen[0])


def cubic(kappa):
    """Positive root, stable also when x is far below machine epsilon."""
    k = np.asarray(kappa, dtype=float)
    x = np.minimum(k, np.cbrt(k))
    for _ in range(24):
        x -= (x*(1+x)**2-k)/((1+x)*(1+3*x))
    return x


def rematch(k0):
    if k0 >= 27/8:
        return None
    # Use log(x) so the solver does not round a tiny positive root to zero.
    def equation(v):
        x=np.exp(v)
        return v+3*np.log1p(1.5*x)-4*np.log1p(x)-np.log(k0)
    logx=brentq(equation, np.log(k0)-4, max(5.,np.log(k0)+4), xtol=1e-12)
    x=np.exp(logx)
    delta=x*(.5+x)/(1+x)**2
    return dict(x0=x, alpha=1-delta, one_minus_alpha=delta)


def fields(mg, p, chi, phase):
    co, si = (1.,0.) if phase=='zero' else (0.,1.)
    U=p['Ui']*np.exp(-2*(chi-p['chi_i']))
    rw=np.sqrt(p['F_squared_eV2']*U/2)/(p['P']*mg)
    S=3*p['P']*mg**2*((co-rw)**2+si**2)
    # The y jet comes from the full auxiliary-field geometry, not from
    # treating F_T as a holomorphic W alone; its factor is one third.
    dSx=6*p['P']*mg**2*rw*(co-rw)
    dSy=2*p['P']*mg**2*rw*si
    return U,S,dSx,dSy


def response(c, mg, p, chi, phase='zero'):
    U,S,dx,dy=fields(mg,p,chi,phase)
    k=4*c['h2']*S/c['M_eV']**4
    x=cubic(k)
    coeff=-x/(1+x)
    R=np.abs(coeff)*np.hypot(dx,dy)/(2*U)
    k0=12*c['h2']*p['P']*mg**2/c['M_eV']**4
    match=rematch(k0)
    if match is None:
        Rm=None
    else:
        xm=cubic(k/match['alpha']**3)
        # Keep 1-alpha separately: at low mG it cannot be obtained by
        # subtracting two ordinary double-precision numbers near one.
        cm=(match['one_minus_alpha']-match['alpha']*xm)/(match['alpha']*(1+xm))
        Rm=np.abs(cm)*np.hypot(dx,dy)/(2*U)
    return dict(U=U,S=S,kappa=k,x=x,ratio=R,ratio_rematched=Rm,match=match,
                dVchi=coeff*dx,dVy=coeff*dy)


def record(c, mg, p, chi, phase='zero'):
    z=response(c,mg,p,chi,phase)
    a=c['tau']*p['P']/(p['A']*c['Lambda_eV']**2)
    f0=np.sqrt(3*p['P'])*mg
    aux=f0/c['min_mass_eV']**2
    light=1e11/c['min_mass_eV']
    heavy=c['max_mass_eV']/1e14
    myoff=4*c['tau']*f0**2/c['Lambda_eV']**2-24*p['A']*mg**2+4*p['rho']/(3*p['P'])
    r=float(np.max(z['ratio']))
    rm=None if z['ratio_rematched'] is None else float(np.max(z['ratio_rematched']))
    gaps=light<=.1*(1+1e-13) and heavy<=.1 and aux<=.1
    return dict(**{k:float(v) if isinstance(v,np.floating) else v for k,v in c.items()},
                mG_eV=float(mg),phase=phase,rH=float(a),
                hidden_self_mass_leading_eV=2*f0/c['Lambda_eV'],
                clockoff_my2_leading_eV2=float(myoff),
                max_x=float(np.max(z['x'])),
                max_static_force_ratio=r,max_once_rematched_force_ratio=rm,
                alpha_delta=None if z['match'] is None else float(z['match']['one_minus_alpha']),
                max_abs_static_dVchi_eV4=float(np.max(np.abs(z['dVchi']))),
                max_abs_static_dVy_eV4=float(np.max(np.abs(z['dVy']))),
                FX_over_min_mass2=float(aux),MQ_over_min_mass=float(light),
                max_mass_over_KK=float(heavy),gap_conditions_pass=bool(gaps),
                static_response_and_gap_conditions_pass=bool(gaps and r<=.1 and myoff>0),
                once_rematched_response_and_gap_conditions_pass=bool(gaps and rm is not None and rm<=.1 and myoff>0))


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    parent=REPO/'experiments/persistent_restoration_v01/results/inputs.json'
    p=json.loads(parent.read_text())['parent_parameters']
    p['A']=float(json.loads(parent.read_text())['C5']*1e28/p['P'])
    trajectory=REPO/'experiments/sugra_shift_v01/results/trajectories.csv'
    with trajectory.open() as f:
        chi=np.array([float(r['chi']) for r in csv.DictReader(f) if r['case']=='shift_axis'])
    checks=[]
    def check(name, ok, detail=''):
        checks.append(dict(name=name,passed=bool(ok),detail=detail))
    matches=[]
    for q in (2.,3.,5.):
        for n in range(2,261):
            c=chain(q,n,1e12)
            matches.append(dict(**c,rH=c['tau']*p['P']/(p['A']*1e24)))
        for n in (2,3,5,8,12):
            H=np.diag(np.full(n,q*q+1))+np.diag(np.full(n-1,-q),1)+np.diag(np.full(n-1,-q),-1)
            G=np.linalg.inv(H); c=chain(q,n,1e12)
            check(f'inverse_{q}_{n}',np.isclose(c['tau'],G[0,-1]/G[0,0],rtol=2e-13,atol=0)
                  and np.isclose(c['G00'],G[0,0],rtol=2e-13))
            check(f'spectrum_{q}_{n}',np.allclose(np.linalg.eigvalsh(H),
                  (q*q+1-2*q*np.cos(np.arange(1,n+1)*np.pi/(n+1))),rtol=2e-13))
            check(f'self_response_{q}_{n}',np.isclose(c['h2'],(G@G)[-1,-1],rtol=2e-13))
    chosen=[]
    for lam in (1e12,2e12):
        for q in (2.,3.,5.):
            candidates=[chain(q,n,lam) for n in range(2,261)]
            valid=[c for c in candidates if c['tau']*p['P']/(p['A']*lam**2)>2]
            c=min(valid,key=lambda c:abs(np.log(c['tau']*p['P']/(p['A']*lam**2)/4)))
            chosen.append(c)
    scan=[record(c,mg,p,chi) for c in chosen for mg in np.logspace(-20,-4,65)]
    reps=[record(c,mg,p,chi,phase) for c in chosen for mg in (1e-6,1e-14,1e-15)
          for phase in ('zero','quadrature')]
    limits=[]
    for c in chosen:
        def fun(lmg,key):
            z=response(c,10.**lmg,p,chi)
            v=z[key]
            return np.log10(np.max(v)) + 1 if v is not None else 100.
        limit=10.**brentq(lambda m:fun(m,'ratio'),-20,-4,xtol=1e-12)
        limitm=10.**brentq(lambda m:fun(m,'ratio_rematched'),-20,-4,xtol=1e-12)
        limits.append(dict(q=c['q'],n=c['n'],Lambda_eV=c['Lambda_eV'],
                           mG_static_10percent_eV=limit,mG_once_rematched_10percent_eV=limitm,
                           MQ_over_min_mass=1e11/c['min_mass_eV'],
                           interpretation='Selected static force only; not a fully matched cosmological bound.'))
    for r in reps:
        c=next(v for v in chosen if v['q']==r['q'] and v['Lambda_eV']==r['Lambda_eV'])
        z=response(c,r['mG_eV'],p,chi,r['phase'])
        check(f'cubic_{r["q"]}_{r["Lambda_eV"]}_{r["mG_eV"]}_{r["phase"]}',
              np.max(np.abs(z['x']*(1+z['x'])**2/z['kappa']-1))<1e-13)
        match=z['match']
        if match:
            x=match['x0']; k0=12*c['h2']*p['P']*r['mG_eV']**2/c['M_eV']**4
            check(f'rematch_{r["q"]}_{r["Lambda_eV"]}_{r["mG_eV"]}_{r["phase"]}',
                  abs(x*(1+1.5*x)**3/(1+x)**4/k0-1)<1e-11)
    check('matching_count',len(matches)==777)
    check('scan_count',len(scan)==390)
    check('representative_count',len(reps)==36)
    check('old_scale_response_failure',all(r['max_static_force_ratio']>1e20 for r in reps if r['mG_eV']==1e-6))
    check('oneTeV_gap_failure',all(not r['gap_conditions_pass'] for r in reps if r['Lambda_eV']==1e12))
    check('twoTeV_low_mass_selected_pass',all(r['static_response_and_gap_conditions_pass'] for r in reps
          if r['Lambda_eV']==2e12 and r['mG_eV']==1e-15))
    check('rematch_branch_finite_bound',rematch(27/8) is None)
    for name,rows in [('chain_matching.csv',matches),('response_scan.csv',scan),
                      ('representatives.csv',reps),('response_limits.csv',limits)]:
        csvout(OUT/name,rows)
    dump(OUT/'inputs.json',dict(parameters=p,chi_first=float(chi[0]),chi_last=float(chi[-1]),
         chi_points=len(chi),chosen_chains=chosen,mG_grid=np.logspace(-20,-4,65).tolist(),
         phase_scan='aligned',transverse_representatives='full-auxiliary leading y jet; not a holomorphic-W substitute'))
    result=dict(created_utc=datetime.now(timezone.utc).isoformat(),numpy=np.__version__,scipy=scipy.__version__,
         checks_passed=sum(c['passed'] for c in checks),check_count=len(checks),checks=checks,
         representatives=[r for r in reps if r['q']==3 and r['phase']=='zero'],limits=limits,
         source_hashes={str(path.relative_to(REPO)):sha(path) for path in
              (Path(__file__),ROOT/'protocol.md',parent,trajectory)},
         scope=['Only declared endpoint matching, local static response and selected gap conditions.',
                'One-time leading vacuum rematch is not full radial recentering or complete matching.',
                'Rigid Gaussian spectrum and full SUGRA response are separate independently scoped reports.',
                'No Riemann origin, observational fit, locality protection, or full cosmic trajectory established.'])
    dump(OUT/'summary.json',result)
    print(f'Main checks {result["checks_passed"]}/{result["check_count"]}')
    for c in checks:
        if not c['passed']: print(c)
    for r in result['representatives']:
        print({k:r[k] for k in ('n','Lambda_eV','mG_eV','tau','rH','max_static_force_ratio',
              'max_once_rematched_force_ratio','MQ_over_min_mass','static_response_and_gap_conditions_pass')})
    if result['checks_passed']!=result['check_count']: raise SystemExit(1)


if __name__=='__main__':
    main()
