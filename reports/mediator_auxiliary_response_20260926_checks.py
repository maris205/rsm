#!/usr/bin/env python3
"""Independent rigid Stückelberg finite-auxiliary calculation.

Only historical inputs are read. No main-stage implementation is imported.
The full local rigid spectrum is derived from K and W in the companion report.
Selected Gaussian matching uses a fixed supersymmetric subtraction at mu=Lambda.
This is not a complete supergravity calculation.
"""
from pathlib import Path
import hashlib
import json
import math

import mpmath as mp
import numpy as np
from scipy.linalg import eigh

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
INPUT = ROOT / 'experiments/persistent_restoration_v01/results/inputs.json'
OUTPUT = HERE / 'mediator_auxiliary_response_20260926_checks.json'
mp.mp.dps = 100
checks = []


def check(name, condition, detail=None):
    checks.append({'name': name, 'passed': bool(condition), 'detail': detail})


def close(name, x, y, tol='1e-75'):
    err = abs(x-y) / max(mp.mpf(1), abs(x), abs(y))
    check(name, err < mp.mpf(tol), {'normalized_error': mp.nstr(err, 12)})


def tridiagonal(n, q):
    return mp.matrix([[q*q+1 if i == j else (-q if abs(i-j) == 1 else 0)
                       for j in range(n)] for i in range(n)])


def determinant(n, q):
    return (q**(2*n+2)-1)/(q*q-1)


def shift(kappa):
    if not kappa:
        return mp.mpf(0)
    return mp.findroot(lambda t: t*(1+t)**2-kappa, (0, kappa), solver='anderson')


def trace(m):
    return sum(m[i,i] for i in range(m.rows))


derivation_cases = []
for n in (1, 2, 4, 8):
    for qs in ('1.01', '2', '3'):
        q = mp.mpf(qs)
        H = tridiagonal(n,q)
        inv = H**-1
        ev, O = mp.eigsy(H)
        ih = O * mp.diag([1/mp.sqrt(v) for v in ev]) * O.T
        e = mp.matrix([0]*(n-1)+[1])
        u = ih*e
        h1 = inv[n-1,n-1]
        h2 = (inv*inv)[n-1,n-1]
        close(f'n{n}_q{qs}_endpoint_inverse', h1, determinant(n-1,q)/determinant(n,q))
        close(f'n{n}_q{qs}_cross_inverse', inv[0,n-1], q**(n-1)/determinant(n,q))
        close(f'n{n}_q{qs}_h2', h2, sum((q**(n-1-j)*determinant(j,q)/determinant(n,q))**2 for j in range(n)))
        for ss in ('0.000001', '0.1', '1'):
            S = mp.mpf(ss)
            M = mp.mpf('2')
            gx = mp.sqrt(2)
            a = gx/M
            kap = 2*gx*gx*h2*S/M**4
            t = shift(kap)
            z = 1+t
            j = 2*S*a/M**2/z**2 * (inv*inv*e)
            V = S/z + M**2/4*(j.T*H*H*j)[0]
            label = f'n{n}_q{qs}_S{ss}'
            close(label+'_stationary_vector', mp.norm(M**2/2*H*H*j-S*a/z**2*e), 0)
            close(label+'_shift', a*j[n-1], t)
            close(label+'_potential', V, S*(1+mp.mpf('1.5')*t)/z**2)
            fun = lambda s: s*(1+mp.mpf('1.5')*shift(2*gx*gx*h2*s/M**4))/(1+shift(2*gx*gx*h2*s/M**4))**2
            close(label+'_envelope', mp.diff(fun,S), 1/z)
            C = a*a*S/z**3
            B = M*M*H+C*u*u.T
            R = M*M*H+4*C*u*u.T
            close(label+'_real_hessian', mp.norm(2*ih*(M*M/2*H*H+2*S*a*a/z**3*e*e.T)*ih-R), 0)
            # Use the full X-dependent inverse metric and D potential.
            def potential_x(xr):
                x2=xr*xr/(2*z)
                dv=H*j+a*x2*e
                return S/(z-a*a*h1*x2)+M*M/4*(dv.T*dv)[0]
            close(label+'_x_scalar_hessian',mp.diff(potential_x,0,2),2*C*h1)
            # Real f can be chosen by a rigid field rephasing. Fermion mass
            # connects n chiral Stückelberg fermions to n gauginos plus psi_X.
            hh=O*mp.diag([mp.sqrt(v) for v in ev])*O.T
            D=M*hh
            b=-a*mp.sqrt(S)/z**mp.mpf('1.5')*u
            FM=mp.zeros(2*n+1)
            for i in range(n):
                for k in range(n):
                    FM[i,n+k]=D[i,k]
                    FM[n+k,i]=D[k,i]
                FM[n+i,2*n]=b[i]
                FM[2*n,n+i]=b[i]
            ef=mp.eigsy(FM,eigvals_only=True)
            expected=mp.eigsy(B,eigvals_only=True)
            close(label+'_fermion_massless_goldstino',ef[n],0)
            for k in range(n):
                close(label+f'_fermion_dirac_{k}',ef[n+1+k]**2,expected[k])
            check(label+'_spectrum_positive',min(ev)>0 and min(expected)>0 and min(mp.eigsy(R,eigvals_only=True))>0 and C>0)
            # Gaussian heavy O(C) term cancels for arbitrary H and rank-one u.
            close(label+'_heavy_linear_supertrace',trace(4*u*u.T)-4*trace(u*u.T),0)
            derivation_cases.append({'n':n,'q':qs,'S':ss,'t':mp.nstr(t,30)})

# Exact F^4 coefficient without subtracting two almost equal energies.
for eps in ('1e-4','1e-8','1e-16'):
    S=mp.mpf(eps); lam=mp.mpf('0.137')
    t=shift(2*lam*S)
    correction=S*(-mp.mpf('0.5')*t-t*t)/(1+t)**2
    residual=correction+lam*S*S
    check('F4_asymptotics_'+eps,abs(residual)/(lam*S*S)<5*lam*S)

parent=json.loads(INPUT.read_text())
p=parent['parent_parameters']
P=mp.mpf(str(p['P']))
Fc=mp.sqrt(mp.mpf(str(p['F_squared_eV2'])))
mg=mp.mpf('1e-6')
chi_i=mp.mpf(str(parent['chi_first']))
chi_f=mp.mpf(str(parent['chi_last']))
Ui=mp.mpf(str(p['Ui']))
q=mp.mpf(3); n=128; gx=mp.sqrt(2); Lam=mp.mpf('1e12')
col=[q**(n-1-j)*determinant(j,q)/determinant(n,q) for j in range(n)]
h1=col[-1]; h2=sum(v*v for v in col)
M=Lam*mp.sqrt(h1)
tau=col[0]/h1
reference=[]
for chi,point in ((chi_i,'old_first'),(chi_f,'old_last')):
    U=Ui*mp.exp(-2*(chi-chi_i))
    rw=Fc*mp.sqrt(U/2)/(mg*P)
    for phase,c,s in (('aligned',mp.mpf(1),mp.mpf(0)),('quadrature',mp.mpf(0),mp.mpf(1))):
        S=3*P*mg**2*(1-2*c*rw+rw*rw)
        Schi=6*P*mg**2*rw*(c-rw)
        # Inherited on-axis auxiliary jet. It is not the naive holomorphic
        # |exp(i theta)-rw exp(-iy)|^2 continuation, whose y derivative is 3x.
        Sy=2*P*mg**2*rw*s
        t=shift(2*gx**2*h2*S/M**4)
        d=-t/(1+t)
        r={'point':point,'phase':phase,'chi':chi,'U_eV4':U,'rw':rw,'S_eV4':S,
           'S_chi_eV4':Schi,'S_y_inherited_eV4':Sy,'t':t,
           'delta_tree_dV_dS':d,'delta_tree_Vchi_eV4':d*Schi,
           'delta_tree_Vy_inherited_eV4':d*Sy,
           'axial_force_ratio':abs(d*Schi)/(2*U),
           'transverse_force_ratio_inherited':abs(d*Sy)/(2*U),
           'gradient_force_ratio_inherited':abs(d)*mp.sqrt(Schi*Schi+Sy*Sy)/(2*U),
           'delta_V_eV4':-S*t*(mp.mpf('.5')+t)/(1+t)**2,
           'C_eV2':gx*gx*S/M**2/(1+t)**3,
           'mX_eV':mp.sqrt(2*gx*gx*S/M**2/(1+t)**3*h1)}
        reference.append({k:(mp.nstr(v,35) if isinstance(v,mp.mpf) else v) for k,v in r.items()})

def gaussian_spectrum(n,ng=24,gravitino=1e-6):
    # Fixed subtraction mu=Lambda, including all rigid physical modes.
    H=np.diag(np.full(n,10.0))+np.diag(np.full(n-1,-3.0),1)+np.diag(np.full(n-1,-3.0),-1)
    ev,O=eigh(H); e=np.zeros(n);e[-1]=1
    u=O@((O.T@e)/np.sqrt(ev)); h1=float(u@u)
    h2=float(np.sum((O.T@e)**2/ev**2))
    la=1e12; mass=la*np.sqrt(h1); S=3*float(P)*gravitino**2
    kap=4*h2*S/mass**4
    t=kap
    for _ in range(30): t=kap/(1+t)**2
    C=2*S/mass**2/(1+t)**3; A=mass*mass*H; uu=np.outer(u,u)
    fp=lambda z:2*z*(np.log(z/la**2)-1)
    fun=lambda z:z*z*(np.log(z/la**2)-1.5)
    def d_heavy(c):
        lb,Qb=eigh(A+4*c*uu)
        lf,Qf=eigh(A+c*uu)
        return 4*((Qb.T@u)**2@fp(lb)-(Qf.T@u)**2@fp(lf))/(64*np.pi**2)
    # Integrate the derivative from C=0, avoiding subtracting O(10^50)
    # vacuum energies to extract O(10^36) residuals.
    x,w=np.polynomial.legendre.leggauss(ng)
    hv=C*sum(ww*d_heavy(C*zz) for zz,ww in zip((x+1)/2,w/2))
    mx2=2*C*h1
    lv=2*fun(mx2)/(64*np.pi**2)
    dc=C/(S*(1+3*t))
    hd=d_heavy(C)*dc
    ld=4*h1*fp(mx2)/(64*np.pi**2)*dc
    return {'n':n,'quadrature_order':ng,'mG_eV':gravitino,'fixed_mu_eV':la,'t':t,'C_eV2':C,
            'V1_heavy_eV4':hv,'V1_X_eV4':lv,'V1_total_eV4':hv+lv,
            'dV1_heavy_dS':hd,'dV1_X_dS':ld,'dV1_total_dS':hd+ld,
            'loop_to_finiteF_tree_slope_ratio':(hd+ld)/(-t/(1+t)),
            'vector_min_mass_eV':mass*np.sqrt(ev.min()),
            'vector_max_mass_eV':mass*np.sqrt(ev.max()),
            'real_heavy_min_mass_eV':float(np.sqrt(eigh(A+4*C*uu,eigvals_only=True).min())),
            'dirac_min_mass_eV':float(np.sqrt(eigh(A+C*uu,eigvals_only=True).min())),
            'X_real_mass_eV':float(np.sqrt(mx2))}

gaussian=[gaussian_spectrum(128,16),gaussian_spectrum(128,32),gaussian_spectrum(96,32)]
for k in ('V1_heavy_eV4','V1_X_eV4','dV1_heavy_dS','dV1_X_dS'):
    for other in gaussian[1:]:
        err=abs(other[k]/gaussian[0][k]-1)
        check('gaussian_convergence_'+k+'_n'+str(other['n'])+'_Q'+str(other['quadrature_order']),err<2e-5,{'relative_error':err})
for r in reference:
    d=gaussian[1]['dV1_total_dS']
    r['selected_rigid_gaussian_axial_force_ratio']=mp.nstr(abs(mp.mpf(d)*mp.mpf(r['S_chi_eV4']))/(2*mp.mpf(r['U_eV4'])),35)
    r['selected_rigid_gaussian_transverse_force_ratio_inherited']=mp.nstr(abs(mp.mpf(d)*mp.mpf(r['S_y_inherited_eV4']))/(2*mp.mpf(r['U_eV4'])),35)


def small_splitting_benchmark(n,la,gravitino):
    """Stable heavy C² coefficient; light X term remains unexpanded."""
    q=mp.mpf(3); la=mp.mpf(la); gravitino=mp.mpf(gravitino)
    col=[q**(n-1-j)*determinant(j,q)/determinant(n,q) for j in range(n)]
    h1=col[-1]; h2=sum(v*v for v in col); mass=la*mp.sqrt(h1)
    S=3*P*gravitino*gravitino;t=shift(4*h2*S/mass**4)
    C=2*S/mass**2/(1+t)**3;cs=C/(S*(1+3*t))
    H=np.diag(np.full(n,10.0))+np.diag(np.full(n-1,-3.0),1)+np.diag(np.full(n-1,-3.0),-1)
    ev,O=eigh(H); w=O[-1,:]**2/ev
    # Dimensionless x=m_vector² / mu². Divided difference of F'(s)
    # is scale invariant and finite on coincident eigenvalues.
    x=float(h1)*ev
    dd=np.empty((n,n))
    for i in range(n):
        for j in range(n):
            dd[i,j]=(2*np.log(x[i]) if i==j else
                     2*((x[i]*np.log(x[i])-x[j]*np.log(x[j]))/(x[i]-x[j])-1))
    L=float(w@dd@w)
    hv=6*C*C*mp.mpf(L)/(64*mp.pi**2)
    hd=12*C*mp.mpf(L)*cs/(64*mp.pi**2)
    mx2=2*C*h1
    lv=2*mx2**2*(mp.log(mx2/la**2)-mp.mpf('1.5'))/(64*mp.pi**2)
    ld=8*h1*mx2*(mp.log(mx2/la**2)-1)*cs/(64*mp.pi**2)
    d=-t/(1+t)
    fields={'n':n,'Lambda_eV':la,'mG_eV':gravitino,'M_eV':mass,
            'tau_exchange':col[0]/h1,'t':t,'C_eV2':C,'X_mass_eV':mp.sqrt(mx2),
            'X_log_m2_over_mu2':mp.log(mx2/la**2),
            'heavy_C2_divided_difference_coefficient':mp.mpf(L),
            'heavy_rank_one_shift_to_min_mass2':4*C*h1/(mass**2*mp.mpf(float(ev.min()))),
            'V1_heavy_C2_eV4':hv,'V1_X_eV4':lv,
            'dV1_heavy_C2_dS':hd,'dV1_X_dS':ld,'dV1_total_dS':hd+ld,
            'loop_to_finiteF_tree_slope_ratio':(hd+ld)/d,'delta_tree_dV_dS':d}
    response=[]
    for chi,point in ((chi_i,'old_first'),(chi_f,'old_last')):
        U=Ui*mp.exp(-2*(chi-chi_i));rw=Fc*mp.sqrt(U/2)/(gravitino*P)
        schi=6*P*gravitino**2*rw*(1-rw);sy=2*P*gravitino**2*rw
        response.append({'point':point,'aligned_tree_force_ratio':mp.nstr(abs(d*schi)/(2*U),35),
                         'quadrature_transverse_tree_force_ratio':mp.nstr(abs(d*sy)/(2*U),35),
                         'aligned_selected_rigid_loop_force_ratio':mp.nstr(abs((hd+ld)*schi)/(2*U),35),
                         'quadrature_transverse_selected_rigid_loop_force_ratio':mp.nstr(abs((hd+ld)*sy)/(2*U),35)})
    return {**{k:(mp.nstr(v,35) if isinstance(v,mp.mpf) else v) for k,v in fields.items()},
            'force_response':response,
            'scope':'Leading C² heavy trace plus unexpanded light-X Gaussian; fixed mu=Lambda; full SUGRA spectrum not supplied.'}


small_splitting=[small_splitting_benchmark(128,'1e12','1e-6'),small_splitting_benchmark(127,'2e12','1e-15')]
for key,exact in (('V1_heavy_C2_eV4',gaussian[1]['V1_heavy_eV4']),('dV1_heavy_C2_dS',gaussian[1]['dV1_heavy_dS'])):
    err=abs(float(small_splitting[0][key])/exact-1)
    check('C2_asymptotic_at_resolved_splitting_'+key,err<.002,{'relative_error':err,'meaning':'Small-C expansion error, not solver tolerance.'})
asymptotic_convergence={'exact':gaussian_spectrum(128,32,1e-7),
                        'C2':small_splitting_benchmark(128,'1e12','1e-7')}
for key,exact in (('V1_heavy_C2_eV4',asymptotic_convergence['exact']['V1_heavy_eV4']),
                  ('dV1_heavy_C2_dS',asymptotic_convergence['exact']['dV1_heavy_dS'])):
    err=abs(float(asymptotic_convergence['C2'][key])/exact-1)
    check('C2_asymptotic_smaller_splitting_'+key,err<.002,
          {'relative_error':err,'meaning':'Same 0.2 percent target at mG=1e-7 eV; original failure retained.'})

# A scalar n=1 spectrum permits direct high precision CW differentiation.
for ss in ('0.0001','0.01','0.1'):
    S=mp.mpf(ss);hh=mp.mpf(10);mass=mp.mpf(2);g=mp.sqrt(2);h=1/hh;h2s=h*h
    def cw(s):
        t=shift(2*g*g*h2s*s/mass**4);c=g*g*s/mass**2/(1+t)**3
        A=mass**2*hh; mu2=mp.mpf(9)
        fun=lambda z:z*z*(mp.log(z/mu2)-mp.mpf('1.5'))
        return (fun(A+4*c*h)-4*fun(A+c*h)+3*fun(A)+2*fun(2*c*h))/(64*mp.pi**2)
    t=shift(2*g*g*h2s*S/mass**4);C=g*g*S/mass**2/(1+t)**3
    fp=lambda z:2*z*(mp.log(z/9)-1)
    derivative=(4*h*(fp(mass**2*hh+4*C*h)-fp(mass**2*hh+C*h))+4*h*fp(2*C*h))*C/(S*(1+3*t))/(64*mp.pi**2)
    close('gaussian_highprecision_derivative_S'+ss,mp.diff(cw,S),derivative)

out={'status':'complete','scope':'Independent local GLOBAL SUSY calculation; inherited auxiliary jets are external response probes, not a complete SUGRA potential.',
     'execution_history':[{'attempt':1,'outcome':'implementation error before any physical calculation',
                           'detail':'mp.matrix(n,n,lambda) produced a zero matrix; replaced by explicit nested-list construction. No threshold or scientific condition changed.'},
                          {'attempt':2,'outcome':'509/509 original exact checks passed'},
                          {'attempt':3,'outcome':'510/511, added two asymptotic accuracy tests; one retained 0.247784 percent versus 0.2 percent failure'},
                          {'attempt':4,'outcome':'Additional bounded smaller-splitting convergence checks; see final counts below'}],
     'arithmetic_digits':100,'input_sha256':hashlib.sha256(INPUT.read_bytes()).hexdigest(),
     'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
     'model':{'J':'S+Sbar+sqrt(2) M V','D_potential':'M^2 j^T H^2 j / 4',
              'q':3,'n':128,'gX':'sqrt(2)','Lambda_eV':'1e12','M_eV':mp.nstr(M,35),
              'endpoint_Hinv':mp.nstr(h1,35),'endpoint_Hinv2':mp.nstr(h2,35),'tau_exchange':mp.nstr(tau,35)},
     'reference_response':reference,'gaussian_spectrum':gaussian,
     'small_splitting_benchmarks':small_splitting,
     'asymptotic_convergence_benchmark':asymptotic_convergence,
     'derivation_cases':derivation_cases,'checks':checks,
     'passed':sum(c['passed'] for c in checks),'total':len(checks),
     'failures':[c for c in checks if not c['passed']]}
OUTPUT.write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({'passed':out['passed'],'total':out['total'],'failures':out['failures'],'model':out['model'],'gaussian':gaussian[1]},indent=2))
raise SystemExit(0 if not out['failures'] else 1)
