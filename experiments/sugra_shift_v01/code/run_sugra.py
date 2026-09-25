#!/usr/bin/env python3
"""Tree SUGRA bosonic-sector test with spectator quantum budgets; no data fit."""
from pathlib import Path
from datetime import datetime, timezone
import csv
import hashlib
import json
import time
import platform
import numpy as np
import scipy
from scipy.integrate import solve_ivp

ROOT=Path(__file__).resolve().parents[1]
REPO=ROOT.parents[1]
PARENT=REPO/'experiments/charged_threshold_v01'
OUT=ROOT/'results'
MODELS=('global','canonical','shift')
KAPPAS=(0.,1.,10.,100.)
NGRID=np.linspace(-np.log(5.2),0,1025)
MI=1e11


def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(p,obj): p.write_text(json.dumps(obj,ensure_ascii=False,indent=2,allow_nan=False)+'\n')
def write_csv(path,rows):
    keys=list(dict.fromkeys(k for row in rows for k in row))
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=keys,lineterminator='\n'); w.writeheader(); w.writerows(rows)


def potential(model,chi,y,inputs):
    eps=inputs['epsilon']; U=inputs['W_i']*np.exp(-2*(chi-inputs['chi_i']))
    zero=np.zeros_like(np.asarray(chi)+np.asarray(y))
    if model=='global':
        return dict(U=U,V=U,Vx=-2*U,Vy=zero,Vxx=4*U,Vxy=zero,Vyy=zero)
    if model=='shift':
        P=1-1.5*eps+eps**2*y*y; base=U*np.exp(eps*y*y)
        V=base*P; Vy=base*2*eps*y*(P+eps)
        Vyy=base*(2*eps*(P+eps)+4*eps**2*y*y*P+8*eps**3*y*y)
        return dict(U=U,V=V,Vx=-2*V,Vy=Vy,Vxx=4*V,Vxy=-2*Vy,Vyy=Vyy)
    if model!='canonical': raise ValueError(model)
    a=1-eps*chi/2; h=-2+eps*chi
    P=a*a+eps**2*y*y/4-1.5*eps
    base=U*np.exp(eps*(chi*chi+y*y)/2)
    V=base*P; Vx=base*(h*P-eps*a)
    Vy=base*eps*y*(P+eps/2)
    Vxx=base*((h*h+eps)*P-2*eps*a*h+eps**2/2)
    Vxy=base*eps*y*(h*(P+eps/2)-eps*a)
    Vyy=base*(eps*(P+eps/2)+eps**2*y*y*(P+eps))
    return dict(U=U,V=V,Vx=Vx,Vy=Vy,Vxx=Vxx,Vxy=Vxy,Vyy=Vyy)


def particle_budget(model,chi,inputs):
    """Real-axis physical spectrum and small-splitting local flat CW only."""
    eps=inputs['epsilon']; ci=inputs['chi_i']; us=inputs['potential_unit_eV4']
    fs=inputs['F_squared_eV2']; mp2=inputs['Mpl_eV']**2
    x=.5*np.expm1(-2*np.log1p((chi-ci)/ci)); g=1+x
    gp=-ci**2/chi**3; gpp=3*ci**2/chi**4
    r=gp/(2*g); rp=(gpp/g-(gp/g)**2)/2
    U=us*inputs['W_i']*np.exp(-2*(chi-ci))
    if model=='canonical':
        k=eps*chi*chi/2; kp=eps*chi; ki=eps*ci*ci/2
        a=1-eps*chi/2; ap=-eps/2; J=a*a-eps; Jp=-eps*a
        T=a*(r+eps*chi/2)+eps/2
        Tp=ap*(r+eps*chi/2)+a*(rp+eps/2)
    else:
        k=0.; kp=0.; ki=0.; J=1-eps; Jp=0.
        T=r+(eps/2 if model=='shift' else 0.)
        Tp=rp
    s=np.exp(k)*MI**2*g; sp=s*(kp+2*r)
    d=np.exp(k)*U*J/mp2 if model!='global' else np.zeros_like(s)
    dp=d*(kp-2+Jp/J)
    h2=2*U*np.exp(2*k)*MI**2*g/fs*T*T
    h2p=h2*(2*kp-2+2*r+2*Tp/T)
    L=(k-ki)+np.log1p(x)
    t=s+d; Lt=L+np.log1p(d/s)
    common=4*s*d*(L-1)+2*d*d*L+(2/3)*d**3/s
    commonp=4*(s*dp*(L-1)+d*sp*L)+4*d*dp*L+2*d*d*sp/s+2*d*d*dp/s-(2/3)*d**3*sp/s**2
    split=2*h2*Lt-h2*h2/(6*t*t)
    splitp=2*h2p*Lt+2*h2*(sp+dp)/t-h2*h2/(6*t*t)*(2*h2p/h2-2*(sp+dp)/t)
    U1=(common+split)/(32*np.pi**2); U1p=(commonp+splitp)/(32*np.pi**2)
    Vx=potential(model,chi,0.,inputs)['Vx']*us
    hratio=np.sqrt(h2)/s; dratio=d/s
    weighted=2*L+(np.log1p(dratio+hratio)+np.log1p(dratio-hratio))/3
    return dict(s=s,sp=sp,d=d,dp=dp,h2=h2,h2p=h2p,U1=U1,U1p=U1p,
                loop_force_ratio=abs(U1p/Vx),physical_mass_eV=np.sqrt(s),
                d_over_s=dratio,h_over_s=hratio,weighted_threshold=weighted,
                log_s_over_initial=L,log_hol_mass2_ratio=np.log1p(x))


class Background:
    def __init__(self,model,yi,inputs):
        self.model,self.yi,self.i=model,yi,inputs
        self.label=f'{model}_axis' if yi==0 else f'{model}_y_{yi:g}'
        V=potential(model,inputs['chi_i'],yi,inputs)['V']
        n=inputs['N_i']; w=inputs['physical_wchi_i']
        ei2=inputs['omega_r']*np.exp(-4*n)+inputs['omega_m']*np.exp(-3*n)+inputs['omega_lambda']+inputs['epsilon']*(.5*w*w+V)/3
        self.Ei=np.sqrt(ei2); self.qi=w/self.Ei
        self.initial=np.array([0.,yi,self.qi,0.,inputs['tau_i'],0.])
    def terms(self,n,state):
        u,y,qx,qy,tau,work=state; chi=self.i['chi_i']+u
        r=potential(self.model,chi,y,self.i)
        rad=self.i['omega_r']*np.exp(-4*n); mat=self.i['omega_m']*np.exp(-3*n)
        q2=qx*qx+qy*qy; den=1-self.i['epsilon']*q2/6
        num=rad+mat+self.i['omega_lambda']+self.i['epsilon']*r['V']/3
        e2=num/den
        if np.any(e2<=0): raise ValueError('Nonpositive Friedmann E2')
        f=-(4*rad+3*mat)/(2*e2)-self.i['epsilon']*q2/2
        r.update(chi=chi,E2=e2,E=np.sqrt(e2),f=f,den=den,num=num,q2=q2,
                 rho=.5*e2*q2+r['V'],rad=rad,mat=mat)
        return r
    def rhs(self,n,state):
        r=self.terms(n,state); qx,qy=state[2:4]
        return [qx,qy,-(3+r['f'])*qx-r['Vx']/r['E2'],
                -(3+r['f'])*qy-r['Vy']/r['E2'],1/r['E'],3*r['E2']*r['q2']]
    def solve(self,tight=False):
        def domain(n,state):
            r=self.terms(n,state)
            return min(20-abs(state[0]),1-abs(state[1]),r['den']-1e-4,r['num']-1e-10)
        domain.terminal=True; domain.direction=-1
        sol=solve_ivp(self.rhs,(NGRID[0],0),self.initial,method='DOP853',
                      rtol=2e-12 if tight else 2e-10,atol=2e-14 if tight else 2e-12,
                      max_step=.0025 if tight else .005,dense_output=True,events=domain)
        return sol
    def outputs(self,sol,n):
        states=sol.sol(n); u,y,qx,qy,tau,work=states; r=self.terms(n,states)
        initial_rho=self.terms(NGRID[0],self.initial)['rho']
        r.update(N=n,z=np.expm1(-n),u=u,y=y,qx=qx,qy=qy,tau=tau,work=work,
                 wx=r['E']*qx,wy=r['E']*qy,
                 t_Gyr=tau/self.i['H_ref_s_inv']/(365.25*86400*1e9),
                 Omega_chi=self.i['epsilon']*r['rho']/(3*r['E2']),
                 energy_residual=(r['rho']+work-initial_rho)/initial_rho)
        if self.yi==0:
            p=particle_budget(self.model,r['chi'],self.i); r.update(p)
            r['transverse_mass2_over_H2']=r['Vyy']/r['E2']
            r['curvature_dimensional_marker']=self.i['H_ref_eV']**2*r['E2']*p['s']/(16*np.pi**2*self.i['potential_unit_eV4']*r['U'])
            r['H_over_mass']=self.i['H_ref_eV']*r['E']/p['physical_mass_eV']
            r['mass_adiabatic']=abs(.5*p['sp']/p['s']*r['wx']*self.i['H_ref_eV']/p['physical_mass_eV'])
            if abs(sol.t[-1])<1e-12:
                d=self.i['alpha_reference']/(4*np.pi)*(p['weighted_threshold']-p['weighted_threshold'][-1])
                r['delta_alpha']=d/(1-d)
        return r


def linear_mode(bg,sol,kappa,tight=False):
    def terms(n):
        state=sol.sol(n); r=bg.terms(n,state)
        mass=r['Vyy']; eps=bg.i['epsilon']; chi=r['chi']
        if bg.model=='global': massp=0.
        elif bg.model=='shift': massp=-2*mass
        else:
            a=1-eps*chi/2
            massp=eps*r['U']*np.exp(eps*chi*chi/2)*((-2+eps*chi)*(a*a-eps)-eps*a)
        omega=kappa*kappa*np.exp(-2*n)+mass
        omegap=-2*kappa*kappa*np.exp(-2*n)+massp*state[2]
        return r,omega,omegap
    def rhs(n,v):
        r,om,omp=terms(n); out=[]
        for j in (0,3):
            D,Dp,_=v[j:j+3]
            out.extend([Dp,-(3+r['f'])*Dp-om/r['E2']*D,
                        3*r['E2']*Dp*Dp-.5*omp*D*D])
        return out
    result=solve_ivp(rhs,(NGRID[0],0),[1,0,0,0,1,0],method='DOP853',
                     rtol=2e-12 if tight else 2e-10,atol=2e-14 if tight else 2e-12,
                     max_step=.001 if tight else .002,dense_output=True)
    if not result.success: raise RuntimeError(result.message)
    v=result.sol(NGRID); r,om,omp=terms(NGRID)
    output=dict(N=NGRID,z=np.expm1(-NGRID),omega2=om,omega2_N=omp)
    for i,j in enumerate((0,3)):
        D,Dp,work=v[j:j+3]
        energy=.5*r['E2']*Dp*Dp+.5*om*D*D
        initial_energy=energy[0]
        output.update({f'D{i}':D,f'D{i}_N':Dp,f'energy{i}':energy,
                       f'energy_residual{i}':(energy+work-initial_energy)/max(initial_energy,1e-30)})
    return output


def main():
    start=time.time(); OUT.mkdir(parents=True,exist_ok=True)
    p=PARENT/'results/inputs.json'; inputs=json.loads(p.read_text())
    old=np.genfromtxt(PARENT/'results/trajectories.csv',names=True,delimiter=',',dtype=None,encoding='utf8')
    old=old[old['case']=='chi_probe']
    inputs.update(physical_wchi_i=float(old['E'][0])*inputs['q_i'],
                  holomorphic_mass_i_eV=MI,parent_inputs_sha256=sha(p),
                  protocol_sha256=sha(ROOT/'protocol.md'),
                  case_definitions=[{'model':m,'y_initial':y} for m in MODELS for y in (0.,.001,.1)])
    dump(OUT/'inputs.json',inputs)
    rows=[]; checks=[]; summaries=[]; axes={}; mode_rows=[]; mode_summaries=[]; outputs={}
    def check(name,value,tol=1e-7):
        checks.append(dict(name=name,value=float(value),tolerance=tol,passed=bool(np.isfinite(value) and value<tol)))
    for model in MODELS:
        for yi in (0.,.001,.1):
            bg=Background(model,yi,inputs)
            try:
                sol=bg.solve(); complete=sol.success and abs(sol.t[-1])<1e-12
                n=NGRID[NGRID<=sol.t[-1]+1e-13]; r=bg.outputs(sol,n); outputs[bg.label]=r
                rec=dict(case=bg.label,model=model,y_initial=yi,complete=bool(complete),message=sol.message,
                         N_last=float(sol.t[-1]),chi_change=float(r['u'][-1]),y_final=float(r['y'][-1]),
                         H_today_over_Href=float(r['E'][-1]),age_today_Gyr=float(r['t_Gyr'][-1]),
                         min_qchi=float(np.min(r['qx'])),max_abs_Omega_chi=float(np.max(abs(r['Omega_chi']))),
                         max_energy_residual=float(np.max(abs(r['energy_residual']))))
                check(bg.label+'_energy',rec['max_energy_residual'])
                if complete:
                    st=bg.solve(tight=True)
                    if not st.success or abs(st.t[-1])>1e-12: raise RuntimeError('Tight run incomplete')
                    rt=bg.outputs(st,n)
                    for key in ('u','y','qx','qy','E','tau'):
                        scale=max(float(np.max(abs(rt[key]))),1e-6 if key in ('y','qy') else 1.)
                        check(bg.label+'_convergence_'+key,np.max(abs(r[key]-rt[key]))/scale)
                if yi==0 and complete:
                    axes[model]=(bg,sol,r)
                    rec.update(delta_alpha_zi_ppm=float(r['delta_alpha'][0]*1e6),
                        physical_mass_initial_GeV=float(r['physical_mass_eV'][0]/1e9),
                        mass_y2_H2_min=float(np.min(r['transverse_mass2_over_H2'])),
                        mass_y2_H2_max=float(np.max(r['transverse_mass2_over_H2'])),
                        max_loop_force_ratio=float(np.max(r['loop_force_ratio'])),
                        max_curvature_dimensional_marker=float(np.max(r['curvature_dimensional_marker'])),
                        max_h_over_s=float(np.max(r['h_over_s'])),max_d_over_s=float(np.max(abs(r['d_over_s']))),
                        min_scalar_mass_ratio=float(np.min(1+r['d_over_s']-r['h_over_s'])),
                        max_H_over_mass=float(np.max(r['H_over_mass'])),
                        max_mass_adiabatic=float(np.max(r['mass_adiabatic'])))
                    check(bg.label+'_spectator_loop_budget',rec['max_loop_force_ratio'])
                    check(bg.label+'_splitting_domain',max(rec['max_h_over_s'],rec['max_d_over_s']),.01)
                    if rec['min_scalar_mass_ratio']<=0: raise RuntimeError('Nonpositive charged scalar mass')
                    if model=='global':
                        for key,newkey in (('u','u'),('E','E'),('tau','tau'),('q','qx')):
                            check('old_reference_'+key,np.max(abs(old[key]-r[newkey]))/max(float(np.max(abs(old[key]))),1.))
                summaries.append(rec)
                for j in range(len(n)):
                    row=dict(case=bg.label,model=model,y_initial=yi)
                    for key,value in r.items():
                        a=np.asarray(value); row[key]=float(a if a.ndim==0 else a[j])
                    rows.append(row)
                print(json.dumps(rec),flush=True)
            except Exception as err:
                summaries.append(dict(case=bg.label,complete=False,exception=repr(err)))
                checks.append(dict(name=bg.label+'_execution',passed=False,exception=repr(err)))
    for model,(bg,sol,r) in axes.items():
        for kappa in KAPPAS:
            try:
                lm=linear_mode(bg,sol,kappa); tight=linear_mode(bg,sol,kappa,tight=True)
                for key in ('D0','D0_N','D1','D1_N'):
                    check(f'{model}_k{kappa:g}_'+key,np.max(abs(lm[key]-tight[key]))/max(float(np.max(abs(tight[key]))),1.))
                for j in (0,1): check(f'{model}_k{kappa:g}_energy{j}',np.max(abs(lm[f'energy_residual{j}'])),1e-6)
                sm=dict(model=model,kappa=kappa,D0_final=float(lm['D0'][-1]),D1_final=float(lm['D1'][-1]),
                        max_abs_D0=float(np.max(abs(lm['D0']))),max_abs_D1=float(np.max(abs(lm['D1']))),
                        omega2_min=float(np.min(lm['omega2'])),omega2_N_max=float(np.max(lm['omega2_N'])),
                        max_energy_residual=float(max(np.max(abs(lm['energy_residual0'])),np.max(abs(lm['energy_residual1'])))))
                if kappa==0:
                    nl=outputs[f'{model}_y_0.001']
                    difference=float(np.max(abs(nl['y']/.001-lm['D0'])))
                    sm['small_nonlinear_vs_linear']=difference
                    check(model+'_small_nonlinear',difference,1e-5)
                    if model=='global': check('massless_homogeneous_constant',np.max(abs(lm['D0']-1)),1e-10)
                mode_summaries.append(sm)
                for j in range(len(NGRID)):
                    mode_rows.append(dict(model=model,kappa=kappa,**{key:float(val[j]) for key,val in lm.items()}))
                print('mode',json.dumps(sm),flush=True)
            except Exception as err:
                checks.append(dict(name=f'{model}_k{kappa:g}_execution',passed=False,exception=repr(err)))
    write_csv(OUT/'trajectories.csv',rows); write_csv(OUT/'linear_modes.csv',mode_rows)
    result=dict(run_utc=datetime.now(timezone.utc).isoformat(),elapsed_seconds=time.time()-start,
                python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,
                protocol_sha256=sha(ROOT/'protocol.md'),code_sha256=sha(__file__),
                summaries=summaries,mode_summaries=mode_summaries,checks=checks,
                checks_passed=sum(x['passed'] for x in checks),check_count=len(checks))
    dump(OUT/'summary.json',result)
    print('Checks',result['checks_passed'],'/',len(checks),flush=True)
    if not all(x['passed'] for x in checks): raise SystemExit(1)


if __name__=='__main__': main()
