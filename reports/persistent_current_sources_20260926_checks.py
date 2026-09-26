#!/usr/bin/env python3
"""Exact, bounded algebra for one invariant massive-vector current.
Own derivations only; no loop integration or observational fit.
"""
import json
from pathlib import Path
import sympy as s

checks=[]
def check(name, expr):
    remainder=s.simplify(expr)
    checks.append({'name':name,'passed':bool(remainder == 0),'remainder':str(remainder)})

z,zb,x,xb,q,qb=s.symbols('z zb x xb q qb')
sigma,tau=s.symbols('sigma tau', real=True)
lam,fx2,P,mg2,A=s.symbols('Lambda FX2 P mG2 A',positive=True)
k=-(z-zb)**2/2
j=k+sigma*q*qb+tau*x*xb
DK=-j**2/lam**2
zero={z:0,zb:0,x:0,xb:0,q:0,qb:0}
gxx=s.diff(DK,x,xb)
check('hidden_metric_jet',gxx+(2*tau*k+2*sigma*tau*q*qb+4*tau**2*x*xb)/lam**2)
DV=-fx2*gxx
chi,y,xr,xi,qr,qi=s.symbols('chi y xr xi qr qi',real=True)
realmap={z:(chi+s.I*y)/s.sqrt(2),zb:(chi-s.I*y)/s.sqrt(2),x:(xr+s.I*xi)/s.sqrt(2),xb:(xr-s.I*xi)/s.sqrt(2),q:(qr+s.I*qi)/s.sqrt(2),qb:(qr-s.I*qi)/s.sqrt(2)}
VR=s.expand(DV.subs(realmap))
expected={'chi':0,'y':4*tau*fx2/lam**2,'xr':4*tau**2*fx2/lam**2,'xi':4*tau**2*fx2/lam**2,'qr':2*sigma*tau*fx2/lam**2,'qi':2*sigma*tau*fx2/lam**2}
for field in (chi,y,xr,xi,qr,qi): check('mass_'+str(field),s.diff(VR,field,2)-expected[str(field)])
check('clock_Hermitian',s.diff(DV,z,zb)-2*tau*fx2/lam**2)
check('clock_holomorphic',s.diff(DV,z,2)+2*tau*fx2/lam**2)
check('charged_clock_correlation',expected['qr']-sigma*expected['y']/2)
check('hidden_clock_correlation',expected['xr']-tau*expected['y'])
u=s.Matrix([1,sigma,tau]); C=u*u.T/lam**2
for i in range(3):
    for j2 in range(i+1,3):check(f'rank_one_minor_{i}_{j2}',C[i,i]*C[j2,j2]-C[i,j2]**2)
check('rank_one_cross_relation',C[0,1]*C[0,2]-C[0,0]*C[1,2])
check('rank_one_det',C.det())
check('nonzero_eigenvalue_trace',s.trace(C)-(1+sigma**2+tau**2)/lam**2)
check('matrix_current_square',DK+(s.Matrix([k,q*qb,x*xb]).T*C*s.Matrix([k,q*qb,x*xb]))[0])
for n in range(4):
    for v in (z,zb,x,xb,q,qb): check(f'lower_jet_{n}_{v}',s.diff(DK,v,n).subs(zero))
# Canonical radion displacement X = sqrt(3P) (T-T*) at g=p=1.
delta,v=s.symbols('delta v',real=True)
DOmega=3*tau**2*P*(delta**2+v**2)**2/(16*lam**2)
qgeom=s.diff(DOmega,delta,2)+s.diff(DOmega,v,2)
check('radion_mixed_quartic',qgeom-3*tau**2*P*(delta**2+v**2)/lam**2)
Vgeom=3*P*mg2*qgeom
Vcanon=Vgeom.subs({delta:s.sqrt(2/(3*P))*xr,v:s.sqrt(2/(3*P))*xi})
check('radion_full_quartic_mass_real',s.diff(Vcanon,xr,2)-12*tau**2*P*mg2/lam**2)
check('radion_full_quartic_mass_imag',s.diff(Vcanon,xi,2)-12*tau**2*P*mg2/lam**2)
check('RSS_bookkeeping_sign',expected['y'].subs({fx2:3*P*mg2,tau:-2*A*lam**2/P})+24*A*mg2)
check('leading_positive_tau_threshold',expected['y'].subs({fx2:3*P*mg2,tau:2*A*lam**2/P})-24*A*mg2)
# Numerics are examples of dimensions/hierarchy only, not bounds or fits.
Mp=2.435e27; Lambda=1e12; kk=1e14
C5=float(s.zeta(3)/(48*s.pi**4)); Aval=C5*kk**2/Mp**2
samples={
    'Mp_eV':Mp,'Lambda_eV':Lambda,'mKK_eV':kk,'C5':C5,'A':Aval,
    'tau_leading_stability_threshold_ignoring_rho':2*Aval*Lambda**2/Mp**2,
    'tau_RSS_mass_reparameterization_not_new_match':-2*Aval*Lambda**2/Mp**2,
    'FX_over_MV2_at_mG_1e_minus_6_eV':3**.5*Mp*1e-6/Lambda**2,
}
result={'scope':'Own local current-square algebra; no naturalness proof, complete spectrum or full loop matching', 'symbolic_engine':'sympy '+s.__version__,'passed':sum(c['passed'] for c in checks),'total':len(checks),'checks':checks,'illustrative_numbers':samples,'independent_review':'A separate internal agent independently derived masses, rank-one identities, holomorphic Hessian and local SUGRA scope without importing this script.'}
Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:result[k] for k in ('passed','total','illustrative_numbers')},indent=2))
if result['passed']!=result['total']:raise SystemExit(1)
