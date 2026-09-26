#!/usr/bin/env python3
"""Two specified Higgs stabilizations; deterministic algebra, not data fitting."""
from pathlib import Path
import csv
import hashlib
import json
import math
import mpmath as mp
import numpy as np
import sympy as sp
from scipy.linalg import eigh_tridiagonal

ROOT=Path(__file__).resolve().parents[1]
REPO=ROOT.parents[1]
mp.mp.dps=120

def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()
def write_json(path,obj): path.write_text(json.dumps(obj,indent=2,ensure_ascii=False)+'\n')
def table(name,rows):
    with (ROOT/'results'/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)
def text(x): return mp.nstr(x,65)

def charge(q,n,eta=None):
    Q=sp.zeros(n+1,n if eta is None else n+1)
    Q[0,0]=q
    for a in range(1,n): Q[a,a-1]=1;Q[a,a]=-q
    Q[n,n-1]=1
    if eta is not None:Q[n,n]=eta
    return Q

def determinant(q,n,eta,s):
    previous,current=mp.mpf(1),mp.mpf(q*q+1)+s
    for _ in range(1,n): previous,current=current,(q*q+1+s)*current-q*q*previous
    return (eta*eta+s)*current-eta*eta*previous

def main():
    (ROOT/'results').mkdir(parents=True,exist_ok=True)
    checks=[]
    def check(name,passed,detail=''): checks.append(dict(name=name,passed=bool(passed),detail=detail))
    rows=[]
    for q in (2,3,5):
        for n in (2,4,8,16,32,64,127,128,160):
            D=(q**(n+1)-1)//(q-1); W=(q**(2*(n+1))-1)//(q*q-1)
            for r in ('0.1','0.5','0.9','0.99'):
                logmu=mp.log10(mp.mpf(W)/2)/2+(D-2)*mp.log10(mp.mpf(r))
                rows.append(dict(q=q,n=n,degree_exact=str(D),v_over_cutoff=r,
                    log10_degree=text(mp.log10(D)),log10_mu_over_v=text(logmu)))
            check(f'degree_sum_{q}_{n}',D==sum(q**a for a in range(n+1)))
    table('phase_lifting.csv',rows)

    for q in (2,3):
        for n in (2,3,5):
            Q=charge(q,n);G=(Q.T*Q).inv();P=Q*G*Q.T
            w=sp.Matrix([1]+[-q**a for a in range(1,n+1)])
            check(f'old_null_{q}_{n}',Q.T*w==sp.zeros(n,1))
            check(f'old_projector_{q}_{n}',P==sp.eye(n+1)-w*w.T/(w.T*w)[0])
            for eta in (sp.Rational(1,2),sp.Integer(1),sp.Integer(2)):
                Qp=charge(q,n,eta);H=Qp.T*Qp;Gp=H.inv(); a=sp.eye(n+1)[:,0];b=Qp[n,:].T
                label=f'{q}_{n}_{eta}'
                check('full_row_kernel_'+label,Qp*Gp*Qp.T==sp.eye(n+1))
                check('local_cross_'+label,(a.T*Gp*b)[0]==0)
                check('local_self_'+label,(b.T*Gp*b)[0]==1 and Gp[0,0]==sp.Rational(1,q*q))
                check('site_cross_'+label,Gp[0,n-1]==sp.Rational(1,q**(n+1)))
                check('physical_h2_'+label,(b.T*Gp*Gp*b)[0]==1/eta**2)
                c=sp.Matrix([sp.Rational(1,q**n)]+[-sp.Rational(1,q**(n-a)) for a in range(1,n)]+[0])
                check('distributed_source_'+label,Qp.T*c==sp.eye(n+1)[:,n-1])
                check('determinant_'+label,H.det()==q**(2*n)*eta**2)
                s=sp.Rational(1,10)
                explicit=(Qp*(H+s*sp.eye(n+1)).inv()*Qp.T)[0,n]/q
                check('momentum_resolvent_'+label,explicit==s*q**(n-1)/(H+s*sp.eye(n+1)).det())

    q,n,Lambda=3,127,mp.mpf('2e12');M=Lambda/q
    f=(mp.mpf(q)**(2*n)-1)/(mp.mpf(q)**(2*(n+1))-1)
    oldtau=(q-1/mp.mpf(q))/(mp.mpf(q)**n-mp.mpf(q)**(-n))
    tau=mp.mpf(q)**(1-n)
    etas=sorted(set([round(float(x),8) for x in np.linspace(.1,3,59)]+[.5,1.,2.]))
    scan=[]
    for eta in etas:
        eig=eigh_tridiagonal(np.array([q*q+1.]*n+[eta*eta]),
                             np.array([-float(q)]*(n-1)+[eta]),eigvals_only=True)
        mvec=float(M)*np.sqrt(eig)
        radial1=float(M)/math.sqrt(2);radial3=3*radial1
        scan.append(dict(q=q,n=n,eta=eta,M_eV=float(M),Lambda_eV=float(Lambda),
            lightest_vector_eV=mvec[0],heaviest_vector_eV=mvec[-1],
            radial_mass_lambda_over_g_1_eV=radial1,radial_mass_lambda_over_g_3_eV=radial3,
            all_heavy_gap_lambda_over_g_3_eV=min(mvec[0],radial3),
            probe_over_gap_lambda_over_g_3=1e11/min(mvec[0],radial3),
            gap_criterion_lambda_over_g_3=bool(1e11/min(mvec[0],radial3)<=.1),
            abstract_tau=text(tau),physical_local_tau='0',
            physical_local_zeta=q*q,abstract_zeta=text(q*q*f/(1-f)),
            physical_source_h2=1/eta**2))
        check(f'positive_spectrum_eta_{eta}',eig[0]>0)
    table('gauge_scan.csv',scan)

    momentum=[]
    for eta in (mp.mpf('.5'),mp.mpf('1'),mp.mpf('2')):
        for exponent in np.linspace(-12,1,131):
            s=mp.power(10,str(float(exponent)))
            det=determinant(q,n,eta,s)
            localcross=s*mp.mpf(q)**(n-1)/det
            taulocal=q*q*localcross
            momentum.append(dict(eta=str(eta),s=text(s),p_over_M=text(mp.sqrt(s)),
                                 local_cross=text(localcross),local_tau=text(taulocal),
                                 ratio_to_abstract_static=text(taulocal/tau)))
        s=mp.mpf('1e-12')
        approx=s/(eta*eta*mp.mpf(q)**(n+1))
        exact=s*mp.mpf(q)**(n-1)/determinant(q,n,eta,s)
        check('small_s_limit_eta_'+str(eta),abs(exact/approx-1)<mp.mpf('1e-9'))
    table('momentum_profile.csv',momentum)

    D=(q**(n+1)-1)//(q-1);W=(q**(2*(n+1))-1)//(q*q-1)
    logr_for_mu_v=-mp.log(mp.mpf(W)/2)/(2*(D-2))
    benchmark=dict(q=q,n=n,Lambda_eV=text(Lambda),M_eV=text(M),
        original_tau=text(oldtau),abstract_tau_after_gauging=text(tau),
        abstract_enhancement=text(tau/oldtau),physical_local_tau='0',
        physical_local_zeta=q*q,abstract_zeta=text(q*q*f/(1-f)),
        degree_exact=str(D),degree_scientific=text(mp.mpf(D)),
        log10_mu_over_v_at_ratio_0p9=text(mp.log10(mp.mpf(W)/2)/2+(D-2)*mp.log10(mp.mpf('.9'))),
        one_minus_ratio_for_leading_mu_equal_v=text(-mp.expm1(logr_for_mu_v)),
        shared_left_row_coefficient_for_old_site_source=text(mp.mpf(q)**(-n)),
        radial_mass_lambda_over_g_3_eV=text(3*M/mp.sqrt(2)),
        eta_benchmarks=[r for r in scan if r['eta'] in (.5,1.,2.)])
    inputfiles=['experiments/threshold_rg_v01/results/inputs.json',
                'experiments/higgs_stabilization_v01/protocol.md',
                'experiments/higgs_stabilization_v01/code/run_stabilization.py']
    result=dict(scope='Deterministic conditional Higgs spectrum/current matching; no full SUGRA or observational fit.',
        input_hashes={p:sha(REPO/p) for p in inputfiles},precision=mp.mp.dps,
        benchmark=benchmark,row_counts={'phase_lifting':len(rows),'gauge_scan':len(scan),'momentum_profile':len(momentum)},
        check_count=len(checks),checks_passed=sum(c['passed'] for c in checks),checks=checks,
        physical_limitations=['Original-field controlled polynomial lift is exceedingly weak.',
          'Extra gauge route removes flat chiral but also the leading local endpoint cross contact.',
          'A nonzero abstract site cross uses a distributed Higgs source with explicitly tiny shared overlap.',
          'F=D=0 heavy spectrum is not the full broken-SUSY or SUGRA spectrum.',
          'p=0 alone does not eliminate higher-superderivative auxiliary-field effects.'])
    write_json(ROOT/'results/summary.json',result)
    print(json.dumps(benchmark,indent=2));print('Main consistency',result['checks_passed'],'/',len(checks))
    for c in checks:
        if not c['passed']:print(c)
    if result['checks_passed']!=len(checks):raise SystemExit(1)

if __name__=='__main__':main()
