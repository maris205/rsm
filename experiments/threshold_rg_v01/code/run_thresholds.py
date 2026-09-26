#!/usr/bin/env python3
"""One-loop Gaussian threshold matching at leading source order S^2.

S is an external source here; no all-order interacting RG or complete SUGRA
matching is implied. Finite boundaries are held fixed when scales change.
"""
from pathlib import Path
from datetime import datetime, timezone
import csv
import hashlib
import json
import numpy as np
from scipy.optimize import brentq

ROOT=Path(__file__).resolve().parents[1]
REPO=ROOT.parents[1]
OUT=ROOT/'results'

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def dump(p,v): p.write_text(json.dumps(v,indent=2,ensure_ascii=False,allow_nan=False)+'\n')
def csvout(p,rows):
    with p.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)

def chain(n,lam,q=3.):
    t=np.arange(1,n+1)*np.pi/(n+1)
    ev=q*q+1-2*q*np.cos(t)
    w=2*np.sin(t)**2/((n+1)*ev)
    g00=float(np.sum(w));h2=float(np.sum(w/ev));M2=lam*lam*g00
    z=(ev[:,None]-ev[None,:])/ev[None,:]
    safe=np.where(z==0,1,z)
    correction=(1+z)*np.log1p(z)/safe-1
    small=abs(z)<1e-5
    correction[small]=z[small]/2-z[small]**2/6+z[small]**3/12-z[small]**4/20
    dd=2*(np.log(g00*ev[None,:])+correction)
    LH=float(np.sum(w[:,None]*w[None,:]*dd))
    # Hats multiply coefficients by Lambda^4, avoiding huge dimensionful
    # cancellations and simplifying the one-loop scale audit.
    a=2*h2/g00**2
    k=1/(2*np.pi**2)
    h=3*LH/(8*np.pi**2*g00**2)
    tau=(q-q**-1)*np.exp(-n*np.log(q))/(-np.expm1(-2*n*np.log(q)))
    return dict(n=n,q=q,Lambda_eV=lam,G00=g00,h2=h2,M2_eV2=M2,
                tau=tau,a_hat=a,k_hat=k,h_hat=h,b_eVminus2=4/lam**2,
                min_mass_eV=np.sqrt(M2*ev[0]),max_mass_eV=np.sqrt(M2*ev[-1]),
                endpoint_weight_sum=float(np.sum(w)),LH_at_Lambda=LH)

def fields(mg,p,chi,phase='zero'):
    co,si=(1.,0.) if phase=='zero' else (0.,1.)
    U=p['Ui']*np.exp(-2*(chi-p['chi_i']))
    rw=np.sqrt(p['F_squared_eV2']*U/2)/(p['P']*mg)
    S=3*p['P']*mg**2*((co-rw)**2+si**2)
    dx=6*p['P']*mg**2*rw*(co-rw)
    dy=2*p['P']*mg**2*rw*si
    return U,S,dx,dy

def coefficients(c,S,mu_h,mu_l):
    la=c['Lambda_eV'];a=c['a_hat'];k=c['k_hat'];h=c['h_hat']
    # All c and h coefficients below are hats, i.e. multiplied by Lambda^4.
    uv=-a+5*k*np.log(mu_h/la)
    hh=h-3*k*np.log(mu_h/la)
    low=uv+hh+2*k*np.log(mu_l/mu_h)
    ell=np.log(c['b_eVminus2']*S/mu_l**2)
    vhat=low+k*(ell-1.5)
    dhat=low+k*(ell-1)
    return dict(cUV_hat=uv,heavy_threshold_hat=hh,cL_hat=low,ell=ell,
                potential_hat=vhat,derivative_hat=dhat)

def physics(c,mg,p,chi,phase='zero'):
    U,S,dx,dy=fields(mg,p,chi,phase)
    ell=np.log(c['b_eVminus2']*S/c['Lambda_eV']**2)
    tree=-2*c['a_hat']*S/c['Lambda_eV']**4
    heavy=2*c['h_hat']*S/c['Lambda_eV']**4
    light=2*c['k_hat']*S*(ell-1)/c['Lambda_eV']**4
    force=tree+heavy+light
    factor=np.hypot(dx,dy)/(2*U)
    delta=c['b_eVminus2']*S/c['Lambda_eV']**2
    scalar_loop=24*delta*np.abs(ell)/(16*np.pi**2)
    goldstino_loop=4*delta*np.abs(ell)/(16*np.pi**2)
    return dict(U=U,S=S,dx=dx,dy=dy,ell=ell,tree=tree,heavy=heavy,light=light,
                total=force,ratio=np.abs(force)*factor,tree_ratio=np.abs(tree)*factor,
                loop_ratio=np.abs(heavy+light)*factor,
                loop_over_new_tree=np.abs((heavy+light)/tree),delta=delta,
                scalar_loop=scalar_loop,goldstino_loop=goldstino_loop,
                mX=np.sqrt(c['b_eVminus2']*S))

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    parent=REPO/'experiments/mediator_chain_v01/results/inputs.json'
    p=json.loads(parent.read_text())['parameters']
    traj=REPO/'experiments/sugra_shift_v01/results/trajectories.csv'
    with traj.open() as f:
        chi=np.array([float(r['chi']) for r in csv.DictReader(f) if r['case']=='shift_axis'])
    cases=[chain(128,1e12),chain(127,2e12)]
    checks=[]
    def check(name,ok,detail=''):checks.append(dict(name=name,passed=bool(ok),detail=detail))
    scale=[];bench=[]
    for c in cases:
        for mg in (1e-17,1e-15,1e-14,1e-12,1e-6):
            f=physics(c,mg,p,chi)
            S=f['S'][-1];mx=f['mX'][-1];lam=c['Lambda_eV']
            ref=coefficients(c,S,lam,lam)
            lows=[.1*mx,.3*mx,mx,3*mx,10*mx,.5*lam,lam]
            for mh_factor in (.5,.75,1.,1.25,2.):
                for ml in lows:
                    z=coefficients(c,S,mh_factor*lam,ml)
                    err=max(abs(z['potential_hat']/ref['potential_hat']-1),
                            abs(z['derivative_hat']/ref['derivative_hat']-1))
                    check(f'path_n{c["n"]}_mG{mg}_h{mh_factor}_l{ml}',err<1e-12)
                    scale.append(dict(n=c['n'],Lambda_eV=lam,mG_eV=mg,muH_eV=mh_factor*lam,
                                      muL_eV=float(ml),**{k:float(v) for k,v in z.items()},
                                      normalized_error=float(err)))
            atmass=coefficients(c,S,lam,mx)
            chainrule_d=atmass['cL_hat']-c['k_hat']
            wrong_norun=-c['a_hat']+c['h_hat']-c['k_hat']
            wrong_nochain=atmass['cL_hat']-1.5*c['k_hat']
            check(f'field_scale_chainrule_n{c["n"]}_{mg}',abs(chainrule_d/ref['derivative_hat']-1)<1e-12)
            check(f'wrong_reset_detected_n{c["n"]}_{mg}',abs(wrong_norun/ref['derivative_hat']-1)>.05)
            bench.append(dict(n=c['n'],Lambda_eV=lam,mG_eV=mg,mX_eV=float(mx),
                  log_mX2_over_Lambda2=float(f['ell'][-1]),
                  cL_hat_at_Lambda=-c['a_hat']+c['h_hat'],cL_hat_at_mX=float(atmass['cL_hat']),
                  force_hat_at_fixed_Lambda=float(ref['derivative_hat']),
                  force_hat_at_mX_with_chainrule=float(chainrule_d),
                  wrong_force_hat_with_reset=float(wrong_norun),
                  wrong_force_hat_no_chainrule=float(wrong_nochain),
                  max_tree_force_ratio=float(np.max(f['tree_ratio'])),
                  max_selected_loop_force_ratio=float(np.max(f['loop_ratio'])),
                  max_selected_total_force_ratio=float(np.max(f['ratio'])),
                  selected_loop_over_new_tree=float(f['loop_over_new_tree'][-1]),
                  low_scalar_log_parameter=float(np.max(f['scalar_loop'])),
                  low_goldstino_log_parameter=float(np.max(f['goldstino_loop'])),
                  leading_S2_scope='old 1e-6 control is an expansion, not the archived exact spectrum'))
    c=cases[1]
    mass=[]
    for mg in np.logspace(-20,-12,65):
        f=physics(c,mg,p,chi)
        mass.append(dict(mG_eV=float(mg),max_selected_total_force_ratio=float(np.max(f['ratio'])),
                         max_tree_force_ratio=float(np.max(f['tree_ratio'])),
                         max_selected_loop_force_ratio=float(np.max(f['loop_ratio'])),
                         selected_loop_over_new_tree=float(f['loop_over_new_tree'][-1]),
                         max_scalar_log_parameter=float(np.max(f['scalar_loop'])),
                         max_goldstino_log_parameter=float(np.max(f['goldstino_loop'])),
                         mX_eV=float(f['mX'][-1]),
                         selected_force_under_10percent=bool(np.max(f['ratio'])<=.1)))
    profiles=[]
    for phase in ('zero','quadrature'):
        f=physics(c,1e-15,p,chi,phase)
        for j,x in enumerate(chi):
            profiles.append(dict(phase=phase,chi=float(x),U_eV4=float(f['U'][j]),
                                 tree_ratio=float(f['tree_ratio'][j]),loop_ratio=float(f['loop_ratio'][j]),
                                 total_ratio=float(f['ratio'][j]),
                                 scalar_log_parameter=float(f['scalar_loop'][j])))
    scale_curve=[]
    f=physics(c,1e-15,p,chi);S=f['S'][-1]
    for ratio in np.logspace(-14,0,141):
        z=coefficients(c,S,c['Lambda_eV'],ratio*c['Lambda_eV'])
        scale_curve.append(dict(mu_over_Lambda=float(ratio),cL_hat=float(z['cL_hat']),
                                explicit_loop_force_hat=float(c['k_hat']*(z['ell']-1)),
                                total_force_hat=float(z['derivative_hat'])))
    bound=10.**brentq(lambda logm:np.log10(np.max(physics(c,10.**logm,p,chi)['ratio']))+1,-20,-12)
    primary=next(r for r in bench if r['n']==127 and r['mG_eV']==1e-15)
    check('primary_total_small',primary['max_selected_total_force_ratio']<.1)
    check('primary_loop_exceeds_new_tree',primary['selected_loop_over_new_tree']>1)
    check('primary_light_loop_parameters_small',primary['low_scalar_log_parameter']<1e-20)
    check('scale_rows',len(scale)==350);check('mass_rows',len(mass)==65);check('profile_rows',len(profiles)==2050)
    for name,rows in [('scale_audit.csv',scale),('benchmarks.csv',bench),('mass_scan.csv',mass),
                      ('clock_profiles.csv',profiles),('scale_curve.csv',scale_curve)]:csvout(OUT/name,rows)
    dump(OUT/'inputs.json',dict(parameters=p,chains=cases,chi_first=float(chi[0]),chi_last=float(chi[-1]),
         chi_count=len(chi),boundary='cUV(Lambda)=-a; extra finite source-squared boundary zero, unchanged under scale variations.',
         approximation='Leading S^2 and the selected one-loop rigid Gaussian sector; S externally prescribed.'))
    result=dict(created_utc=datetime.now(timezone.utc).isoformat(),checks_passed=sum(v['passed'] for v in checks),
                check_count=len(checks),checks=checks,primary=primary,selected_total_10percent_mG_eV=float(bound),
                chains=cases,max_path_error=max(r['normalized_error'] for r in scale),
                source_hashes={str(v.relative_to(REPO)):sha(v) for v in (Path(__file__),ROOT/'protocol.md',parent,traj)},
                scope=['one-loop leading source-squared matching only','no deletion of finite-boundary logarithm',
                       'small low-energy repeated-loop parameters do not bound all UV matching',
                       'no protected chain locality, complete SUGRA, observations, or cosmic-time integration'])
    dump(OUT/'summary.json',result)
    print(f'Main {result["checks_passed"]}/{result["check_count"]}')
    print(json.dumps(primary,indent=2))
    print('Selected 10% mG:',bound)
    for v in checks:
        if not v['passed']:print(v)
    if result['checks_passed']!=result['check_count']:raise SystemExit(1)

if __name__=='__main__':main()
