#!/usr/bin/env python3
"""Retain the uneaten Higgs chiral: geometry and conditional finite-F valley."""
from pathlib import Path
import csv,hashlib,json,math
import mpmath as mp
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp,quad

ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1]
mp.mp.dps=110
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def fmt(x):return mp.nstr(x,60)
def write(p,x):p.write_text(json.dumps(x,ensure_ascii=False,indent=2)+'\n')
def csvwrite(name,rows):
    with (ROOT/'results'/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)
def null(q,n):
    w=[mp.mpf(1)]+[-mp.mpf(q)**a for a in range(1,n+1)]
    norm=mp.sqrt(mp.fsum(a*a for a in w));return [a/norm for a in w]
def solve_h(q,n,rhs):
    d=[mp.mpf(q*q+1)]*n;b=list(rhs)
    for j in range(1,n):
        f=-mp.mpf(q)/d[j-1];d[j]-=f*(-q);b[j]-=f*b[j-1]
    x=[mp.mpf(0)]*n;x[-1]=b[-1]/d[-1]
    for j in range(n-2,-1,-1):x[j]=(b[j]+q*x[j+1])/d[j]
    return x
def source_response(q,n):
    rhs=[mp.mpf(0)]*n;rhs[-1]=1
    g=solve_h(q,n,rhs);x=solve_h(q,n,g)
    a=[q*x[0]]+[x[j-1]-q*x[j] for j in range(1,n)]+[x[-1]]
    return mp.fsum(y*y for y in g),a
def troot(B,kappa):
    if not kappa:return mp.mpf(0)
    t=min(kappa/(B*B),mp.root(kappa,3))
    for _ in range(80):
        new=t-(t*(B+t)**2-kappa)/((B+t)*(B+3*t))
        if abs(new-t)<mp.mpf('1e-100')*max(abs(new),mp.mpf('1e-100')):return new
        t=new
    raise RuntimeError('positive t root did not converge')

def main():
    (ROOT/'results').mkdir(parents=True,exist_ok=True)
    checks=[]
    def check(name,ok,detail=''):checks.append(dict(name=name,passed=bool(ok),detail=detail))
    inputs=json.loads((REPO/'experiments/threshold_rg_v01/results/inputs.json').read_text())
    P=mp.mpf(str(inputs['parameters']['P']));hb=mp.mpf(str(inputs['parameters']['hbar_eV_s']))
    q,n=3,127;Lambda=mp.mpf('2e12');g=mp.mpf(1);v=Lambda/(2*g*q);mg=mp.mpf('1e-15')
    S=3*P*mg*mg;u=null(q,n);P00=1-u[0]**2
    cI=1/(v*Lambda*mp.sqrt(P00));cX=q*cI
    betaI=2*v*cI*u[0];betaX=2*v*cX*u[-1];b=-2*v*v*cX*u[-1]
    tau=-q*u[0]*u[-1]/P00;h,avec=source_response(q,n);kappa=S*cX*cX*h/g**2
    check('normalized_null',abs(mp.fsum(x*x for x in u)-1)<mp.mpf('1e-100'))
    check('heavy_source_orthogonal',abs(mp.fsum(u[j]*avec[j] for j in range(n+1)))<mp.mpf('1e-100'))
    check('source_pseudoinverse_diagonal',abs(avec[-1]/h-1)<mp.mpf('1e-100'))
    check('source_product_cross_identity',abs(betaI*betaX*Lambda**2/(-4*tau)-1)<mp.mpf('1e-100'))
    endpoint=[]
    for qq in (2,3,5):
        for nn in (2,4,8,16,32,64,127,128):
            uu=null(qq,nn);pp=1-uu[0]**2;tt=-qq*uu[0]*uu[-1]/pp
            endpoint.append(dict(q=qq,n=nn,u_left=fmt(uu[0]),u_right=fmt(uu[-1]),
                tau=fmt(tt),right_light_weight=fmt(uu[-1]**2),right_heavy_weight=fmt(1-uu[-1]**2),
                betaI_times_Lambda=fmt(2*uu[0]/mp.sqrt(pp)),
                betaX_times_Lambda=fmt(2*qq*uu[-1]/mp.sqrt(pp))))
            check(f'endpoint_tau_{qq}_{nn}',abs(tt/((qq-1/mp.mpf(qq))/(mp.mpf(qq)**nn-mp.mpf(qq)**(-nn)))-1)<mp.mpf('1e-95'))
    csvwrite('endpoint_scaling.csv',endpoint)
    csvwrite('mode_profile.csv',[dict(link=j,u=fmt(x),light_weight=fmt(x*x),heavy_diagonal=fmt(1-x*x)) for j,x in enumerate(u)])

    # Exact rational short chains: no inverse of a singular Gram matrix is assumed.
    for qq in (2,3):
        for nn in (1,2,4):
            Q=sp.zeros(nn+1,nn);Q[0,0]=qq
            for j in range(1,nn):Q[j,j-1]=1;Q[j,j]=-qq
            Q[nn,nn-1]=1;ww=sp.Matrix([1]+[-qq**j for j in range(1,nn+1)])
            A=Q*Q.T;Ap=Q*(Q.T*Q).inv()**2*Q.T
            check(f'rational_projector_{qq}_{nn}',A*Ap==sp.eye(nn+1)-ww*ww.T/(ww.T*ww)[0])
            check(f'finite_source_no_stationarity_{qq}_{nn}',ww.T*A==sp.zeros(1,nn+1) and ww[-1]!=0)
    physical=[]
    for j in range(101):
        r=mp.mpf(j)/100;B=1+b*r;t=troot(B,kappa);D=B+t
        kin0=mp.fsum(x*x/mp.sqrt(1+x*x*r*r) for x in u)
        alpha=S*cX/(g*g*D*D)
        dralpha=-2*S*cX*b/(g*g*D*D*(D+2*t))
        dif=[-2*v*v*r*u[a]+alpha*avec[a] for a in range(n+1)]
        dd=[-2*v*v*u[a]+dralpha*avec[a] for a in range(n+1)]
        kinfinite=mp.fsum(dd[a]**2/(2*mp.sqrt(4*v**4+dif[a]**2)*v*v) for a in range(n+1))
        energy=1/D+t/(2*D*D);slope=-b/(D*D)
        physical.append(dict(r=fmt(r),Dflat_metric_X=fmt(B),finiteF_metric_X=fmt(D),
            heavy_metric_shift=fmt(t),V_over_S=fmt(energy),dV_dr_over_S=fmt(slope),
            d2V_dr2_over_S=fmt(2*b*b/(D*D*(D+2*t))),
            radial_kinetic_B=fmt(kin0),finiteF_radial_kinetic_B=fmt(kinfinite),
            kahler_metric_C=fmt(1/kin0),clock_metric=fmt(1+cI*dif[0]),
            force_relative_Dflat_error=fmt(abs((B/D)**2-1)),
            kinetic_relative_Dflat_error=fmt(abs(kinfinite/kin0-1))))
        check('physical_monotonic_'+str(j),slope<0 and D>0 and 1+cI*dif[0]>0)
        check('physical_heavy_shift_'+str(j),abs(t*D**2/kappa-1)<mp.mpf('1e-95'))
    csvwrite('physical_valley.csv',physical)
    stress=[]
    for kap in ('0','1e-6','1e-3','0.1','1'):
        kk=mp.mpf(kap)
        for j in range(61):
            r=mp.mpf(j)/20;B=1+b*r;t=troot(B,kk);D=B+t
            stress.append(dict(kappa=kap,r=fmt(r),hidden_metric=fmt(D),t=fmt(t),
                V_hat=fmt(kk/D+t*t/2),dV_hat_dr=fmt(-kk*b/D**2),
                V_over_S=fmt(1/D+t/(2*D*D)) if kk else 'undefined_zero_source'))
            check('stress_root_'+kap+'_'+str(j),abs(t*D**2-kk)<mp.mpf('1e-90'))
    csvwrite('source_stress.csv',stress)

    # A limited homogeneous Minkowski transient using the leading small-F geometry.
    uf=np.array([float(x) for x in u]);bf=float(b);target=.1/bf
    def metric(r):return np.sum(uf*uf/np.sqrt(1+uf*uf*r*r))
    def metricprime(r):return -r*np.sum(uf**4/(1+uf*uf*r*r)**1.5)
    def rhs(s,y):
        r,rd=y;B=metric(r)
        return [rd,bf/(B*(1+bf*r)**2)-metricprime(r)*rd*rd/(2*B)]
    def event(s,y):return y[0]-target
    event.terminal=True;event.direction=1
    sol=solve_ivp(rhs,(0,2),[0,0],method='DOP853',rtol=2e-12,atol=2e-14,
                  events=event,dense_output=True,max_step=.025)
    if not sol.success or not len(sol.t_events[0]):raise RuntimeError('local transient event not reached')
    st=float(sol.t_events[0][0]);timeunit=float(hb*v/mp.sqrt(S))
    # r=y² removes the integrable inverse-square-root endpoint in energy quadrature.
    sq=quad(lambda y:np.sqrt(2*metric(y*y)*(1+bf*y*y)/bf),0,math.sqrt(target),epsabs=1e-13,epsrel=1e-13)[0]
    check('transient_energy_quadrature',abs(st/sq-1)<1e-10,abs(st/sq-1))
    transient=[];maxenergy=0
    for s in np.linspace(0,st,201):
        r,rd=sol.sol(s);energy=.5*metric(r)*rd*rd+1/(1+bf*r)
        maxenergy=max(maxenergy,abs(energy-1))
        transient.append(dict(s=s,local_elapsed_seconds=s*timeunit,r=r,dr_ds=rd,
            hidden_metric_change=bf*r,energy_over_S=energy,kinetic_B=metric(r)))
    check('transient_energy_conservation',maxenergy<1e-10,maxenergy)
    csvwrite('local_transient.csv',transient)
    benchmark=dict(q=q,n=n,g=fmt(g),Lambda_eV=fmt(Lambda),v_eV=fmt(v),mg_eV=fmt(mg),S_eV4=fmt(S),
        cI_eVminus2=fmt(cI),cX_eVminus2=fmt(cX),u_left=fmt(u[0]),u_right=fmt(u[-1]),
        betaI_eVminus1=fmt(betaI),betaX_eVminus1=fmt(betaX),tau=fmt(tau),h=fmt(h),kappa=fmt(kappa),b=fmt(b),
        leading_sigma_force_at_origin_eV3=fmt(-S*betaX),
        leading_sigma_Hessian_eV2=fmt(2*S*betaX**2),
        leading_clock_partial_Hessian_eV2=fmt(4*tau*S/Lambda**2),
        leading_clock_covariant_Hessian_eV2=fmt(2*tau*S/Lambda**2),
        leading_X_partial_Hessian_eV2=fmt(S*(4/Lambda**2+betaX**2/2)),
        hessian_scope='Leading local entries at a nonstationary origin; not vacuum particle masses or full rolling stability.',
        transient_time_unit_seconds=timeunit,transient_event_s=st,transient_energy_quadrature_s=sq,
        transient_event_seconds=st*timeunit,transient_event_r=target,max_transient_energy_error=maxenergy,
        max_finiteF_force_relative_error=fmt(max(mp.mpf(x['force_relative_Dflat_error']) for x in physical)),
        max_finiteF_kinetic_relative_error=fmt(max(mp.mpf(x['kinetic_relative_Dflat_error']) for x in physical)))
    deps=['experiments/threshold_rg_v01/results/inputs.json','experiments/light_modulus_v01/protocol.md',
          'experiments/light_modulus_v01/code/run_light_modulus.py']
    out=dict(scope='Conditional rigid light-modulus geometry and finite-F stationary obstruction; no new cosmic solution.',
        source_hashes={p:sha(REPO/p) for p in deps},precision_digits=mp.mp.dps,benchmark=benchmark,
        row_counts={'endpoint_scaling':len(endpoint),'mode_profile':len(u),'physical_valley':len(physical),
                    'source_stress':len(stress),'local_transient':len(transient)},
        check_count=len(checks),checks_passed=sum(c['passed'] for c in checks),checks=checks,
        retained_physical_failure='No finite interior stationary point for S>0,cX!=0 on specified X=Za=0 branch.',
        exclusions=['all field branches','SUGRA stabilization','loop effects','large-field UV completion',
                    'cosmic lifetime prediction','Riemann derivation of physical interactions or inverse-log exponent'])
    write(ROOT/'results/summary.json',out)
    print(json.dumps(benchmark,indent=2));print('Main consistency',out['checks_passed'],'/',len(checks))
    for c in checks:
        if not c['passed']:print(c)
    if out['checks_passed']!=len(checks):raise SystemExit(1)

if __name__=='__main__':main()
