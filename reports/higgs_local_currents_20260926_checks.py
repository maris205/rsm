#!/usr/bin/env python3
"""Independent polynomial-Higgs current matching; no main implementation import."""
from pathlib import Path
import hashlib
import json
import platform
import sys
import sympy as sp
import mpmath as mp

ROOT = Path(__file__).resolve().parents[1]
OUT = Path(__file__).with_suffix('.json')
checks = []

def check(name, passed, detail=''):
    checks.append({'name': name, 'passed': bool(passed), 'detail': str(detail)})

def eq(name, a, b):
    check(name, sp.simplify(a-b) == 0, sp.simplify(a-b))

def Qmatrix(n, q, eta=None):
    Q = sp.zeros(n+1, n+(eta is not None))
    Q[0,0] = q
    for a in range(1,n):
        Q[a,a-1], Q[a,a] = 1,-q
    Q[n,n-1] = 1
    if eta is not None:
        Q[n,n] = eta
    return Q

def matrixeq(name, a, b):
    check(name, all(sp.simplify(x) == 0 for x in a-b))

for n in (1,2,3,5):
    for q in (sp.Rational(3),sp.Rational(5,3)):
        tag=f'n{n}_q{q}'
        Q=Qmatrix(n,q)
        H=Q.T*Q
        w=sp.Matrix([1]+[-q**a for a in range(1,n+1)])
        W=(w.T*w)[0]
        P=Q*H.inv()*Q.T
        matrixeq(tag+'_null_vector',Q.T*w,sp.zeros(n,1))
        matrixeq(tag+'_projector',P,sp.eye(n+1)-w*w.T/W)
        matrixeq(tag+'_idempotent',P*P,P)
        eq(tag+'_end_cross',P[0,n],q**n/W)
        eq(tag+'_self_ratio',P[n,n]/P[0,0],1/q**2)
        tausite=(H.inv())[0,n-1]/(H.inv())[0,0]
        eq(tag+'_cross_normalization',P[0,n]/P[0,0],tausite/q)
        eq(tag+'_rescaled_right_tau',q*P[0,n]/P[0,0],tausite)
        eq(tag+'_rescaled_right_self',q*q*P[n,n]/P[0,0],1)
        for eta in (sp.Rational(1),sp.Rational(1,3)):
            tag2=tag+f'_eta{eta}'
            R=Qmatrix(n,q,eta)
            Hp=R.T*R
            invH=Hp.inv()
            matrixeq(tag2+'_square_identity',R*invH*R.T,sp.eye(n+1))
            eq(tag2+'_determinant',R.det()**2,eta**2*q**(2*n))
            D=R*invH**2*R.T
            matrixeq(tag2+'_derivative_identity',D,(R*R.T).inv())
            eq(tag2+'_derivative_end',D[0,n],-1/(eta**2*q**n))
            eq(tag2+'_site_end_cross',invH[0,n-1],q**(-n-1))
            eq(tag2+'_site_left_self',invH[0,0],q**(-2))
            eq(tag2+'_site_right_self',invH[n-1,n-1],(1-q**(-2*n))/(q*q-1))
            s=sp.Rational(1,1000)
            C=R*(Hp+s*sp.eye(n+1)).inv()*R.T
            matrixeq(tag2+'_finite_momentum_identity',C,
                     sp.eye(n+1)-s*(R*R.T+s*sp.eye(n+1)).inv())
            check(tag2+'_finite_end_positive',C[0,n] > 0)

# Exact separable quotient, independent of the chain inverse.
j=sp.symbols('j', real=True)
x=-sp.atanh(j)
eq('quotient_stationary',sp.sinh(x)+j*sp.cosh(x),0)
eq('quotient_value',sp.cosh(x)+j*sp.sinh(x),sp.sqrt(1-j*j))
eq('quotient_second_jet',sp.diff(sp.sqrt(1-j*j),j,2).subs(j,0),-1)
I,S,a,b=sp.symbols('I S a b',real=True)
eff=sp.sqrt(1-a*a*I*I)+sp.sqrt(1-b*b*S*S)
eq('quotient_no_cross_exact',sp.diff(eff,I,S),0)

# Direct numerical resolvents at two tiny momenta, not an expansion-only test.
mp.mp.dps=100
finite_results=[]
for n in (2,5,9):
    for eta in (mp.mpf('1'),mp.mpf('0.1')):
        R=mp.matrix(Qmatrix(n,sp.Rational(3),sp.Rational(str(eta))).tolist())
        Hp=R.T*R
        for ss in ('1e-20','1e-40'):
            s=mp.mpf(ss)
            C=R*(Hp+s*mp.eye(n+1))**-1*R.T
            leading=s/(eta**2*3**n)
            rel=abs(C[0,n]/leading-1)
            check(f'mp_n{n}_eta{eta}_s{ss}_leading',rel < 100*s/eta**2,mp.nstr(rel,20))
            finite_results.append({'n':n,'eta':str(eta),'s':ss,
                                   'exact_cross':mp.nstr(C[0,n],45),
                                   'leading_cross':mp.nstr(leading,45),
                                   'relative_error':mp.nstr(rel,20)})

# Original low-energy benchmark: 127 sites, but independent 128x128 solve.
mp.mp.dps=180
n=127
q=mp.mpf(3)
eta=mp.mpf(1)
R=mp.zeros(n+1)
R[0,0]=q
for z in range(1,n):
    R[z,z-1],R[z,z]=1,-q
R[n,n-1],R[n,n]=1,eta
Hp=R.T*R
eR=mp.matrix([0]*n+[1])
right_source=R.T*eR
sol=mp.lu_solve(Hp,right_source)
cross=(R*sol)[0]
selfR=(R*sol)[n]
check('n127_direct_square_endpoint_cross',abs(cross)<mp.mpf('1e-150'),mp.nstr(cross,40))
check('n127_direct_square_endpoint_self',abs(selfR-1)<mp.mpf('1e-150'),mp.nstr(selfR-1,40))
W=sum(q**(2*z) for z in range(n+1))
P0n=q**n/W
tau_site=(q-q**-1)/(q**n-q**-n)
tau_local=P0n/(1-1/W)
check('n127_local_tau_normalization',abs(q*tau_local/tau_site-1)<mp.mpf('1e-150'))
benchmark={'n':n,'q':'3','eta':'1','digits':mp.mp.dps,
           'W':mp.nstr(W,45),'before_gauging_P0n':mp.nstr(P0n,45),
           'before_gauging_P00':mp.nstr(1-1/W,45),
           'before_gauging_Pnn':mp.nstr(1-q**(2*n)/W,45),
           'tau_equal_local_couplings':mp.nstr(tau_local,45),
           'tau_right_coupling_times_q':mp.nstr(tau_site,45),
           'after_gauging_cross_exact':'0',
           'after_gauging_endpoint_derivative_coefficient':mp.nstr(q**-n/eta**2,45),
           'derivative_to_old_contact_ratio_over_s':mp.nstr(W/(eta**2*q**(2*n)),45),
           'dimensionless_s_illustration':'1e-90',
           'derivative_cross_at_illustration':mp.nstr(mp.mpf('1e-90')*q**-n/eta**2,45),
           'illustration_is_observational_bound':False}

inputs=['reports/chain_locality_breaking_20260926.md',
        'experiments/threshold_rg_v01/protocol.md']
result={'status':'pass' if all(x['passed'] for x in checks) else 'fail',
        'passed':sum(x['passed'] for x in checks),'total':len(checks),
        'scope':'Own tree-level rigid-SUSY polynomial local Higgs-current matching; not full SUGRA, loop matching, or observations.',
        'precision':{'exact':'sympy rational','finite_resolvents_digits':100,'large_direct_solve_digits':180},
        'tolerances_changed_after_run':False,'previous_failed_runs':[],
        'python':platform.python_version(),'sympy':sp.__version__,'mpmath':mp.__version__,
        'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'input_sha256':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in inputs},
        'benchmark':benchmark,'finite_resolvents':finite_results,'checks':checks}
OUT.write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n')
print(json.dumps({k:result[k] for k in ('status','passed','total')},ensure_ascii=False))
if result['status']!='pass':
    print(json.dumps([x for x in checks if not x['passed']],indent=2))
    sys.exit(1)
