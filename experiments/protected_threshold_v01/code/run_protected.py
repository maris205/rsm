#!/usr/bin/env python3
"""Protected charged spectrum: cancellation-safe one-loop diagnostic, not a fit."""
from pathlib import Path
import csv
import hashlib
import json
import math
import platform
import time
from datetime import datetime, timezone
import numpy as np
import scipy
from scipy.integrate import solve_ivp

ROOT=Path(__file__).resolve().parents[1]
REPO=ROOT.parents[1]
OLD=REPO/'experiments/charged_threshold_v01'
OUT=ROOT/'results'
ETA=(0.,1.,100.,1e4,1e6)
MI=1e11
NGRID=np.linspace(-np.log(5.2),0,1025)


def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(p,obj): p.write_text(json.dumps(obj,ensure_ascii=False,indent=2,allow_nan=False)+'\n')
def label(eta): return f'soft_{eta:g}'


def mass_terms(chi,mi,inputs):
    ci=inputs['chi_i']; us=inputs['potential_unit_eV4']; fs=inputs['F_squared_eV2']
    u=chi-ci
    gm1=.5*np.expm1(-2*np.log1p(u/ci))
    g=1+gm1; gp=-ci**2/chi**3; gpp=3*ci**2/chi**4
    logg=np.log1p(gm1)
    W=inputs['W_i']*np.exp(-2*u)
    y2=mi**2*gp**2/(2*fs*g)
    ly=2*gpp/gp-gp/g
    h2=us*W*y2
    lh=-2+ly
    return dict(g=g,gm1=gm1,gp=gp,gpp=gpp,logg=logg,W=W,y2=y2,ly=ly,
                h2=h2,lh=lh,s=mi**2*g,sp=mi**2*gp)


def _raw_loop(chi,mi,eta,mu_ratio,inputs):
    r=mass_terms(chi,mi,inputs)
    d=8*np.pi**2*inputs['potential_unit_eV4']*eta/mi**2
    s,sp,h2=r['s'],r['sp'],r['h2']
    L=r['logg']-2*np.log(mu_ratio)
    t=s+d; Lt=L+np.log1p(d/s)
    common=4*s*d*(L-1)+2*d*d*L+(2/3)*d**3/s
    commonp=4*d*L*sp+2*d*d*sp/s-(2/3)*d**3*sp/s**2
    split=2*h2*Lt-h2*h2/(6*t*t)
    splitp=2*h2*(r['lh']*Lt+sp/t)-h2*h2/(6*t*t)*(2*r['lh']-2*sp/t)
    r.update(d=d,L=L,U_eV4=(common+split)/(32*np.pi**2),
             Up_eV4=(commonp+splitp)/(32*np.pi**2),split_raw=split,
             common_raw=common,force_F_eV4=splitp/(32*np.pi**2))
    return r


def evaluate_loop(chi,mi,eta,mu_ratio,inputs):
    """Return raw, derivative, and initial-constant-subtracted potential in eV^4.

    Keep d^3 and h^4 terms. Domain and independent direct-supertrace errors are
    archived. No addition of another Utree/(1+deltaK) potential.
    """
    r=_raw_loop(chi,mi,eta,mu_ratio,inputs)
    first=_raw_loop(inputs['chi_i'],mi,eta,mu_ratio,inputs)
    x=r['gm1']; d=r['d']; L0=-2*np.log(mu_ratio)
    # (1+x)log(1+x)-x, stable at the matching point.
    rem=sum((-1)**n*x**n/(n*(n-1)) for n in range(2,31))
    common_delta=4*d*mi**2*(rem+x*L0)+2*d*d*r['logg']+(2/3)*d**3/mi**2*(-x/(1+x))
    r['deltaU_eV4']=(common_delta+r['split_raw']-first['split_raw'])/(32*np.pi**2)
    r['d_over_s']=d/r['s']; r['h_over_s']=np.sqrt(r['h2'])/r['s']
    r['deltaK']=-r['y2']*r['L']/(16*np.pi**2)
    r['deltaKp']=-r['y2']*(r['ly']*r['L']+r['gp']/r['g'])/(16*np.pi**2)
    return r


class Model:
    def __init__(self,inputs,eta): self.i,self.eta=inputs,eta
    def terms(self,n,u,q):
        chi=self.i['chi_i']+u
        r=evaluate_loop(chi,MI,self.eta,1.,self.i)
        A=1+r['deltaK']; Ap=r['deltaKp']
        V=r['W']+r['deltaU_eV4']/self.i['potential_unit_eV4']
        Vp=-2*r['W']+r['Up_eV4']/self.i['potential_unit_eV4']
        eps=self.i['epsilon']
        rad=self.i['omega_r']*np.exp(-4*n); mat=self.i['omega_m']*np.exp(-3*n)
        den=1-eps*A*q*q/6
        num=rad+mat+self.i['omega_lambda']+eps*V/3
        e2=num/den
        if np.any(e2<=0): raise ValueError('Nonpositive Friedmann E2')
        f=-(4*rad+3*mat)/(2*e2)-eps*A*q*q/2
        r.update(A=A,Ap=Ap,V=V,Vp=Vp,E=np.sqrt(e2),E2=e2,den=den,num=num,f=f,chi=chi,
                 rad=rad,mat=mat)
        return r
    def rhs(self,n,y):
        u,q,tau,work=y; r=self.terms(n,u,q)
        return [q,-(3+r['f'])*q-r['Ap']*q*q/(2*r['A'])-r['Vp']/(r['A']*r['E2']),
                1/r['E'],3*r['A']*r['E2']*q*q]
    def solve(self,tight=False):
        def domain(n,y):
            r=self.terms(n,y[0],y[1])
            return min(20-abs(y[0]),.25-abs(r['logg']),r['den']-1e-4,r['num']-1e-10,
                       .01-r['d_over_s'],.01-r['h_over_s'],r['A'])
        domain.terminal=True; domain.direction=-1
        return solve_ivp(self.rhs,(NGRID[0],0),[0,self.i['q_i'],self.i['tau_i'],0],
                         method='DOP853',rtol=2e-12 if tight else 2e-10,
                         atol=2e-14 if tight else 2e-12,max_step=.0025 if tight else .005,
                         dense_output=True,events=domain)
    def outputs(self,sol,n):
        u,q,tau,work=sol.sol(n); r=self.terms(n,u,q)
        rho=.5*r['A']*r['E2']*q*q+r['V']
        ri=self.terms(NGRID[0],0,self.i['q_i'])
        rhoi=.5*ri['A']*ri['E2']*self.i['q_i']**2+ri['V']
        # Threshold retains nonzero common log even when tiny splittings round off.
        logplus=r['logg']+np.log1p(r['d_over_s']+r['h_over_s'])
        logminus=r['logg']+np.log1p(r['d_over_s']-r['h_over_s'])
        weighted=(logplus+logminus)/3+4*r['logg']/3
        r.update(N=n,z=np.expm1(-n),u=u,q=q,tau=tau,work=work,rho=rho,
                 t_Gyr=tau/self.i['H_ref_s_inv']/(365.25*86400*1e9),
                 energy_residual=(rho+work-rhoi)/rhoi,
                 Omega_chi=self.i['epsilon']*rho/(3*r['E2']),
                 loop_force_over_tree=abs(r['Up_eV4'])/(2*self.i['potential_unit_eV4']*r['W']))
        if abs(sol.t[-1])<1e-12:
            c=self.i['alpha_reference']/(4*np.pi)
            B=1-c*(weighted-weighted[-1])
            delta=(1-B)/B
            # stable direct small difference instead of subtracting B from one
            d=c*(weighted-weighted[-1]); delta=d/(1-d)
            s_clock=np.expm1(-2*np.log(r['chi']/r['chi'][-1]))/r['chi'][-1]**2
            gamma=self.i['alpha_reference']/(2*np.pi)*self.i['chi_i']**2/(1+self.i['chi_i']**2/r['chi'][-1]**2)
            lnt=np.log(tau/(self.i['H_ref_s_inv']*self.i['tstar_s']))
            age=np.expm1(-2*np.log(lnt/lnt[-1]))/lnt[-1]**2
            r.update(delta_alpha=delta,alpha_linear=gamma*s_clock,B=B,
                     clock_shape=s_clock/s_clock[0],age_shape=age/age[0])
        return r


def write_csv(p,rows):
    with p.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n')
        w.writeheader(); w.writerows(rows)


def budget(inputs,old):
    chi=old['chi']; q=old['q']; E=old['E']; W=old['tree']; us=inputs['potential_unit_eV4']
    records=[]
    for mi in (1.,1e6,1e11):
        for mu in (.5,1.,2.):
            f=evaluate_loop(chi,mi,0,mu,inputs)
            s=evaluate_loop(chi,mi,1,mu,inputs)
            coefficient=np.max(abs(s['Up_eV4']-f['Up_eV4'])/(2*us*W))
            eta_max=.1/coefficient
            softmax=np.sqrt(8*np.pi**2*us*eta_max/mi**2)
            # Static transverse potential curvature + coordinate metric correction.
            xi=inputs['chi_i']**2
            D=(3*chi**4+3*chi**2*xi+2*xi**2)/(chi**2*(chi**2+xi)**2)
            EE=-xi*(3*chi**2+xi)/(chi**2*(chi**2+xi)**2)
            hyy=f['y2']*(EE-D*f['L'])/(8*np.pi**2)
            transverse=hyy*(W/E**2+q*q/2)/(1+f['deltaK'])
            records.append(dict(mass_eV=mi,mu_over_mi=mu,
                max_F_force_over_tree=float(np.max(abs(f['Up_eV4'])/(2*us*W))),
                max_abs_deltaK=float(np.max(abs(f['deltaK']))),
                max_h_over_s=float(np.max(f['h_over_s'])),
                max_soft_force_per_eta=float(coefficient),eta_for_force_point1=float(eta_max),
                msoft_for_force_point1_eV=float(softmax),
                min_transverse_coord_mass2_over_H2=float(np.min(transverse)),
                max_abs_transverse_coord_mass2_over_H2=float(np.max(abs(transverse)))))
    # Large-soft example: exact functions are numerically safe at d/s~1e-4.
    mi=1e11; d=1e18; r=mass_terms(chi,mi,inputs); ss=r['s']; L=r['logg']
    lp=L+np.log1p(d/ss)
    Up=((ss+d)*(lp-1)-ss*(L-1))*r['sp']/(8*np.pi**2)
    sug=np.exp(inputs['epsilon']*chi**2/2)*((1-inputs['epsilon']*chi/2)**2-1.5*inputs['epsilon'])
    return dict(records=records,large_soft_example=dict(mass_eV=mi,msoft_eV=1e9,
        max_force_over_tree=float(np.max(abs(Up)/(2*us*W))),
        scope='Fixed common soft term, no claim of actual visible-sector transmission'),
        minimal_canonical_SUGRA_over_global_range=[float(min(sug)),float(max(sug))],
        interpretation='Force<0.1 is a diagnostic on the old tree trajectory; not an experimental bound. Varying mu at fixed finite inputs compares matching choices, not an RG-invariant prediction.')


def main():
    started=time.time(); OUT.mkdir(parents=True,exist_ok=True)
    inputs=json.loads((OLD/'results/inputs.json').read_text())
    inputs.update(protected_mass_eV=MI,eta_soft=list(ETA),mu_over_mi=1.,
                  case_definitions=[{'case':label(eta),'eta_soft':eta,'mass_eV':MI,
                                     'mu_over_mi':1.,'kinetic':'canonical_chi'} for eta in ETA],
                  parent_inputs_sha256=sha(OLD/'results/inputs.json'),
                  protocol_sha256=sha(ROOT/'protocol.md'))
    dump(OUT/'inputs.json',inputs)
    old=np.genfromtxt(OLD/'results/trajectories.csv',names=True,delimiter=',',dtype=None,encoding='utf8')
    old=old[old['case']=='chi_probe']
    dump(OUT/'budgets.json',budget(inputs,old))
    rows=[]; records=[]; checks=[]
    def check(name,value,tol=1e-7):
        checks.append(dict(name=name,value=float(value),tolerance=tol,passed=bool(np.isfinite(value) and value<tol)))
    for eta in ETA:
        try:
            m=Model(inputs,eta); sol=m.solve(); n=NGRID[NGRID<=sol.t[-1]+1e-13]
            r=m.outputs(sol,n)
            rec=dict(case=label(eta),eta_soft=eta,complete=bool(sol.success and abs(sol.t[-1])<1e-12),
                status=sol.message,N_last=float(sol.t[-1]),msoft_eV=float(np.sqrt(r['d'])),
                chi_change=float(r['u'][-1]),q_min=float(np.min(r['q'])),
                clock_reverses=bool(np.any(r['q']<0)),max_energy_residual=float(np.max(abs(r['energy_residual']))),
                max_force_ratio=float(np.max(r['loop_force_over_tree'])),
                max_abs_deltaK=float(np.max(abs(r['deltaK']))),max_d_over_s=float(np.max(r['d_over_s'])),
                max_h_over_s=float(np.max(r['h_over_s'])),min_A=float(np.min(r['A'])),
                min_scalar_mass_ratio=float(np.min(1+r['d_over_s']-r['h_over_s'])),
                max_abs_Omega_chi=float(np.max(abs(r['Omega_chi']))))
            check(label(eta)+'_energy',rec['max_energy_residual'])
            if rec['complete']:
                tight=m.solve(tight=True); rt=m.outputs(tight,n)
                if not tight.success or abs(tight.t[-1])>1e-12: raise RuntimeError('Tight integration incomplete')
                for key in ('u','q','E','tau','rho','delta_alpha'):
                    scale=max(float(np.max(abs(rt[key]))),1e-12 if key=='delta_alpha' else 1.)
                    check(label(eta)+'_convergence_'+key,np.max(abs(r[key]-rt[key]))/scale)
                rec.update(delta_alpha_zi_ppm=float(r['delta_alpha'][0]*1e6),
                    linear_error_over_peak=float(np.max(abs(r['delta_alpha']-r['alpha_linear']))/np.max(abs(r['delta_alpha']))),
                    clock_age_shape_error=float(np.max(abs(r['clock_shape']-r['age_shape']))))
                if eta==0:
                    for key in ('u','q','E','tau'):
                        check('parent_tree_'+key,np.max(abs(r[key]-old[key]))/max(float(np.max(abs(old[key]))),1.))
            records.append(rec)
            for j in range(len(n)):
                row={'case':label(eta)}
                for key,val in r.items():
                    ar=np.asarray(val); row[key]=float(ar if ar.ndim==0 else ar[j])
                rows.append(row)
            print(json.dumps(rec),flush=True)
        except Exception as error:
            records.append(dict(case=label(eta),complete=False,exception=repr(error)))
            checks.append(dict(name=label(eta)+'_execution',passed=False,exception=repr(error)))
    write_csv(OUT/'trajectories.csv',rows)
    result=dict(run_utc=datetime.now(timezone.utc).isoformat(),python=platform.python_version(),
        numpy=np.__version__,scipy=scipy.__version__,elapsed_seconds=time.time()-started,
        protocol_sha256=sha(ROOT/'protocol.md'),code_sha256=sha(__file__),
        summaries=records,checks=checks,checks_passed=sum(x['passed'] for x in checks),check_count=len(checks))
    dump(OUT/'summary.json',result)
    print('Checks',result['checks_passed'],'/',len(checks),flush=True)
    if not all(x['passed'] for x in checks): raise SystemExit(1)


if __name__=='__main__': main()
