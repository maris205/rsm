#!/usr/bin/env python3
"""Independent finite-auxiliary-field audit; imports no experiment implementation.

The tested rigid branch has X=Z_a=0, W_clock=0, real positive Higgs
expectation values, and positive Kahler metric. No SUGRA claim is made.
"""
from pathlib import Path
import hashlib
import json
import platform
import sys
import numpy as np
import scipy
from scipy.optimize import least_squares
import sympy as sp

HERE = Path(__file__).resolve()
OUT = HERE.with_suffix('.json')
checks = []

def check(name, ok, detail=None):
    checks.append({'name': name, 'passed': bool(ok), 'detail': detail})

def eq(name, a, b):
    diff = sp.simplify(a-b)
    check(name, diff == 0, str(diff))

def charge(n, q, eta, symbolic=False):
    R = sp.zeros(n+1) if symbolic else np.zeros((n+1,n+1))
    R[0,0] = q
    for a in range(1,n):
        R[a,a-1],R[a,a] = 1,-q
    R[n,n-1],R[n,n] = 1,eta
    return R

# 1. Symbolic finite-S Gaussian shift, including the auxiliary denominator.
for n in (1,2,3,4):
    for eta in (sp.Rational(1,2), sp.Rational(2)):
        R = charge(n,sp.Rational(3),eta,True)
        H = R.T*R
        a = sp.zeros(n+1,1); a[0] = 1
        b = R[n,:].T
        inv = H.inv()
        eq(f'gaussian_n{n}_eta{eta}_overlap', (a.T*inv*b)[0], 0)
        kval,M = sp.symbols('k M', real=True, positive=True)
        u = sp.Matrix(sp.symbols(f'u0:{n+1}',real=True))
        j = u-sp.sqrt(2)*inv*a*kval/M
        eq(f'gaussian_n{n}_eta{eta}_F_denominator',
           1+sp.sqrt(2)*(b.T*j)[0]/M,
           1+sp.sqrt(2)*(b.T*u)[0]/M)
        check(f'gaussian_n{n}_eta{eta}_D_shift',
              all(sp.simplify(x)==0 for x in H*j+sp.sqrt(2)*a*kval/M-H*u))
        eq(f'gaussian_n{n}_eta{eta}_hidden_susceptibility',
           (b.T*inv*inv*b)[0],eta**-2)

# 2. Explicit complex-coordinate metric block check for the canonical Higgs model.
T,Tb,X,Xb = sp.symbols('T Tb X Xb')
z = sp.symbols('z0:2'); zb = sp.symbols('zb0:2')
p = sp.symbols('p0:2'); pb = sp.symbols('pb0:2')
m = sp.symbols('m0:2'); mb = sp.symbols('mb0:2')
cI,cX = sp.symbols('cI cX',real=True)
k = -(T-Tb)**2/2
d = [p[i]*pb[i]-m[i]*mb[i] for i in range(2)]
K = k+X*Xb+sum(p[i]*pb[i]+m[i]*mb[i]+z[i]*zb[i] for i in range(2))
K += cI*k*d[0]+cX*X*Xb*d[1]
fields = [T,X,*p,*m,*z]
bars = [Tb,Xb,*pb,*mb,*zb]
branch = {X:0,Xb:0,**{a:0 for a in z+zb}}
metric = sp.Matrix([[sp.diff(K,a,b).subs(branch) for b in bars] for a in fields])
for j in range(len(fields)):
    eq(f'canonical_X_metric_column_{j}',metric[1,j],1+cX*d[1] if j==1 else 0)
for i in (6,7):
    for j in range(len(fields)):
        eq(f'canonical_Z_metric_{i}_{j}',metric[i,j],int(i==j))
f,lam,v = sp.symbols('f lam v',real=True)
W = f*X+lam*sum(z[i]*(p[i]*m[i]-v*v) for i in range(2))
for i,field in enumerate(fields):
    expected = f if i==1 else lam*(p[i-6]*m[i-6]-v*v) if i in (6,7) else 0
    eq(f'canonical_W_derivative_{i}',sp.diff(W,field).subs(branch),expected)

# 3. Direct heavy-field stationary solves, using log amplitudes rather than d.
#    phi_a=exp(r_a+xi_a), tilde_phi_a=exp(r_a-xi_a), v=1.
#    Starting from zero every time avoids seeding the analytic stationary answer.
rows = []
max_gradient = 0.
max_energy_span = 0.
min_metric_schur = np.inf
for n in (1,2,4):
    for eta in (.5,1.,2.):
        R = charge(n,3.,eta)
        gram = R@R.T
        count = n+1
        g,lam,c_i,c_x = .7,2.,.2,.3
        for S in (.01,.5,4.):
            energies=[]
            for kval in (0.,.01,.1,.5):
                def parts(z):
                    rr,xi = z[:count],z[count:]
                    prod = np.exp(2*rr)
                    dif = 2*prod*np.sinh(2*xi)
                    sums = 2*prod*np.cosh(2*xi)
                    source = dif.copy()
                    source[0] += c_i*kval*sums[0]
                    force = g*g*(gram@source)
                    denom = 1+c_x*dif[-1]
                    energy = S/denom+lam*lam*np.sum((prod-1)**2)+.5*g*g*np.dot(source,gram@source)
                    # Independent analytic derivatives of the full nonlinear potential.
                    gd = force.copy()
                    gd[-1] -= S*c_x/(denom*denom)
                    grad_r = 2*dif*gd+4*lam*lam*prod*(prod-1)
                    grad_xi = 2*sums*gd
                    grad_r[0] += 2*c_i*kval*sums[0]*force[0]
                    grad_xi[0] += 2*c_i*kval*dif[0]*force[0]
                    return energy,np.r_[grad_r,grad_xi],prod,dif,sums,force,denom
                solve = least_squares(lambda z:parts(z)[1],np.zeros(2*count),
                                      ftol=1e-13,xtol=1e-13,gtol=1e-13,max_nfev=3000)
                energy,grad,prod,dif,sums,force,denom = parts(solve.x)
                gradmax = float(np.max(np.abs(grad)))
                direct_k_force = float(c_i*sums[0]*force[0])
                schur = float(1+c_i*dif[0]-2*kval*c_i*c_i*(
                    (sums[0]+dif[0])/2/(1+c_i*kval)+
                    (sums[0]-dif[0])/2/(1-c_i*kval)))
                tag=f'nonlinear_n{n}_eta{eta}_S{S}_k{kval}'
                check(tag+'_stationarity',solve.success and gradmax<1e-10,gradmax)
                check(tag+'_left_pair_product',abs(prod[0]-1)<1e-10,float(prod[0]-1))
                check(tag+'_clock_force',abs(direct_k_force)<1e-11,direct_k_force)
                check(tag+'_metric_positive',denom>0 and schur>0 and 1-c_i*kval>0,
                      {'hidden_block':float(denom),'clock_schur':schur})
                # Auxiliary equation follows from the solved field values, not input.
                t=c_x*dif[-1]
                rhs=S*c_x*c_x/(g*g*eta*eta)
                check(tag+'_auxiliary_branch_equation',abs(t*(1+t)**2-rhs)<1e-10,
                      float(t*(1+t)**2-rhs))
                energies.append(float(energy))
                max_gradient=max(max_gradient,gradmax)
                min_metric_schur=min(min_metric_schur,schur)
                rows.append({'n':n,'eta':eta,'S':S,'k':kval,'energy':float(energy),
                             'gradient_max':gradmax,'dV_dk':direct_k_force,
                             'clock_metric_schur':schur,'hidden_metric':float(denom),
                             'left_difference':float(dif[0]),'right_difference':float(dif[-1]),
                             'n_function_evaluations':int(solve.nfev)})
            span=max(energies)-min(energies)
            max_energy_span=max(max_energy_span,span)
            check(f'nonlinear_n{n}_eta{eta}_S{S}_finite_k_energy_flat',span<1e-11,span)

# 4. Exact nonlinear coordinate transformation and stationarity identities.
d0,p0,h,target = sp.symbols('d0 p0 h target',real=True)
s0=sp.sqrt(d0*d0+4*p0*p0)
xi=sp.symbols('xi',real=True)
dxi=2*p0*sp.sinh(2*xi)
sxi=2*p0*sp.cosh(2*xi)
eq('higgs_scaling_derivative',sp.diff(dxi+h*sxi,xi),2*(sxi+h*dxi))
eq('higgs_scaling_product_invariant',sp.diff(p0*sp.exp(xi)*sp.exp(-xi),xi),0)
# The positive branch of d+h sqrt(d²+4p²)=target requires |h|<1.
for pv in (.3,1.,2.):
    for hv in (-.9,-.1,0.,.1,.9):
        for tv in (-2.,0.,3.):
            dv=(tv-hv*np.sqrt(tv*tv+4*pv*pv*(1-hv*hv)))/(1-hv*hv)
            residual=dv+hv*np.sqrt(dv*dv+4*pv*pv)-tv
            check(f'exact_nonlinear_map_p{pv}_h{hv}_target{tv}',abs(residual)<1e-11,float(residual))

payload={'scope':'Rigid W=fX branch, X=Z_a=0, W_clock=0, positive Kahler metric. No full-SUGRA or general higher-superderivative null result.',
         'status':'PASS' if all(x['passed'] for x in checks) else 'FAIL',
         'total_checks':len(checks),'passed_checks':sum(x['passed'] for x in checks),
         'numerical_stationary_points':len(rows),'max_stationarity_residual':max_gradient,
         'max_fixed_source_energy_span':max_energy_span,'minimum_clock_metric_schur':min_metric_schur,
         'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,'sympy':sp.__version__,
         'script_sha256':hashlib.sha256(HERE.read_bytes()).hexdigest(),'checks':checks,'stationary_results':rows}
OUT.write_text(json.dumps(payload,ensure_ascii=False,indent=2)+'\n')
print(json.dumps({k:v for k,v in payload.items() if k not in ('checks','stationary_results')},indent=2))
if payload['status']!='PASS':
    print(json.dumps([x for x in checks if not x['passed']],indent=2))
    sys.exit(1)
