#!/usr/bin/env python3
"""Independent algebra and high-precision local checks of specified stabilization.

Does not import this round's other implementations or fit observational data.
Dimensionless numerical examples verify asymptotic structure, not cosmology.
"""
from pathlib import Path
import hashlib
import json
import platform
import sympy as s
import mpmath as mp

OUT = Path(__file__).with_suffix('.json')
checks = []

def check(name, lhs, rhs=0):
    residual = s.factor(s.simplify(s.expand(lhs-rhs)))
    checks.append(dict(name=name, passed=residual == 0, residual=str(residual)))

def numeric(name, value, target, tol):
    err = abs(value-target)
    checks.append(dict(name=name, passed=bool(err <= tol),
                       absolute_error=mp.nstr(err, 12), tolerance=str(tol)))

P,Y=s.symbols('P Y', positive=True, real=True)
kz,kzb,J,Jb,W,Wb,wz,wzb,ff=s.symbols('kz kzb J Jb W Wb wz wzb ff')
G=s.Matrix([[3*P/Y**2,-kzb/Y**2],[-kz/Y**2,1/Y+kz*kzb/(3*P*Y**2)]])
Gi=s.Matrix([[Y**2/(3*P)+Y*kz*kzb/(9*P**2),Y*kzb/(3*P)],
             [Y*kz/(3*P),Y]])
for row in range(2):
    for col in range(2):
        check(f'exact_Kahler_metric_inverse_{row}_{col}',(G*Gi-s.eye(2))[row,col])
D=s.Matrix([J-3*W/Y,wz+kz*W/(P*Y)])
Db=s.Matrix([Jb-3*Wb/Y,wzb+kzb*Wb/(P*Y)])
A=wz+kz*J/(3*P)
Ab=wzb+kzb*Jb/(3*P)
N=A*Ab+ff+Y*J*Jb/(3*P)-(J*Wb+Jb*W)/P
raw=((Db.T*Gi*D)[0]+Y*ff-3*W*Wb/P)/Y**3
check('direct_DW_potential_with_nilpotent_F',raw,N/Y**2)
check('unstabilized_J_zero_exact_cancellation',(N/Y**2).subs({J:0,Jb:0}),(wz*wzb+ff)/Y**2)

# Expand the entire bosonic potential, then canonically normalize Q by sqrt(Y).
qp,qm,qpb,qmb,M,Mb,Mz,Mzb=s.symbols('qp qm qpb qmb M Mb Mz Mzb')
norm=qp*qpb+qm*qmb
Yq=Y-norm/(3*P)
Wq=W+M*qp*qm
Wqb=Wb+Mb*qpb*qmb
Aq=wz+Mz*qp*qm+kz*J/(3*P)
Aqb=wzb+Mzb*qpb*qmb+kzb*Jb/(3*P)
Qp=M*qm+qpb*J/(3*P)
Qpb=Mb*qmb+qp*Jb/(3*P)
Qm=M*qp+qmb*J/(3*P)
Qmb=Mb*qpb+qm*Jb/(3*P)
Vq=(Aq*Aqb+ff+Qp*Qpb+Qm*Qmb+Yq*J*Jb/(3*P)
    -(J*Wqb+Jb*Wq)/P)/Yq**2
zero={qp:0,qm:0,qpb:0,qmb:0}
V=N/Y**2
diag=M*Mb/Y+2*V/(3*P)
B=(Ab*Mz-Jb*M/(3*P))/Y
Bb=(A*Mzb-J*Mb/(3*P))/Y
check('charged_diagonal_Qplus',Y*s.diff(Vq,qp,qpb).subs(zero),diag)
check('charged_diagonal_Qminus',Y*s.diff(Vq,qm,qmb).subs(zero),diag)
check('charged_holomorphic_B',Y*s.diff(Vq,qp,qm).subs(zero),B)
check('charged_antiholomorphic_B',Y*s.diff(Vq,qpb,qmb).subs(zero),Bb)
check('charged_no_hermitian_mixing',s.diff(Vq,qp,qmb).subs(zero))
check('charged_no_same_charge_holomorphic',s.diff(Vq,qp,qp).subs(zero))
check('canonical_Dirac_mass_squared',M*Mb/Y**3/(1/Y)**2,M*Mb/Y)
check('common_soft_term_exact',diag-M*Mb/Y,2*V/(3*P))

# One constant uplift: exact stationarity and prescribed rho in clock-off vacuum.
a,b,cR,cI,kap,rho=s.symbols('a b cR cI kap rho',real=True)
r=a+s.I*b
Yr=1+2*a
Wr=kap*(cR+s.I*cI+r*r/2)
Jr=kap*r
fr=s.symbols('fr',real=True)
Nr=fr+Yr*Jr*s.conjugate(Jr)/(3*P)-(Jr*s.conjugate(Wr)+s.conjugate(Jr)*Wr)/P
Vref=s.simplify(s.expand_complex(Nr/Yr**2))
eta=P*rho/kap**2
cRstar=(a-s.Rational(3,2)*a*a-b*b/2-6*eta*Yr)/3
cIstar=b*(1-a)/3
fstar=(1-2*a)*(kap**2*(a*a+b*b)/(3*P)+rho*Yr)
sub={cR:cRstar,cI:cIstar,fr:fstar}
check('reference_uplift_energy',Vref.subs(sub),rho)
check('reference_uplift_real_stationarity',s.diff(Vref,a).subs(sub))
check('reference_uplift_imag_stationarity',s.diff(Vref,b).subs(sub))
check('modulus_mass_in_SUSY_Minkowski',s.diff(Vref,a,2).subs({a:0,b:0,cR:0,cI:0,fr:0})/(6*P),kap**2/(9*P**2))
check('axion_mass_in_SUSY_Minkowski',s.diff(Vref,b,2).subs({a:0,b:0,cR:0,cI:0,fr:0})/(6*P),kap**2/(9*P**2))
u=s.symbols('u',real=True)
Vwithw=(fr+Yr*Jr*s.conjugate(Jr)/(3*P)-(Jr*(s.conjugate(Wr)+u)+s.conjugate(Jr)*(Wr+u))/P)/Yr**2
check('envelope_linear_real_w_coefficient',s.diff(Vwithw,u),-2*kap*a/(P*Yr**2))
cv=s.symbols('cv',real=True)
check('quadrature_finite_kappa_real_shift_leading',
      s.series((a-s.Rational(3,2)*a*a-b*b/2).subs({a:s.Rational(9,2)*cv**2,b:3*cv}),cv,0,3).removeO())
check('quadrature_finite_kappa_imag_shift_leading',
      s.series((b*(1-a)-3*cv).subs({a:s.Rational(9,2)*cv**2,b:3*cv}),cv,0,3).removeO())

# Strong-stabilization limit: eliminate the auxiliary-sized J before freezing T.
Dgeom=Y+kz*kzb/(3*P)
Jstar=(3*W-kzb*wz)/Dgeom
Jbstar=(3*Wb-kz*wzb)/Dgeom
check('J_stationary_in_strong_limit',s.diff(N,J).subs({Jb:Jbstar}))
Neff=wz*wzb+ff-(3*W-kzb*wz)*(3*Wb-kz*wzb)/(3*P*Dgeom)
check('strong_limit_completed_square',N.subs({J:Jstar,Jb:Jbstar}),Neff)
# G rows are holomorphic and columns antiholomorphic; G^{-1} has the
# opposite index ordering, so the contravariant F vector is Gi.T * Db.
check('strong_limit_F_T_zero',((Gi.T*Db)[0]).subs({Jb:Jbstar}))

F,rootU,m=s.symbols('F rootU m',positive=True,real=True)
y,theta=s.symbols('y theta',real=True)
ep=F**2/P
Y0=1-ep*y*y/3
D0=1+ep*y*y/3
C=s.sqrt(2)*F*m*rootU
w=-F*rootU/s.sqrt(2)*s.exp(-s.I*y)+m*P*s.exp(s.I*theta)
wprime=rootU*s.exp(-s.I*y)
kprime=-s.I*s.sqrt(2)*F*y
f2=3*m*m*P+rho
Vstrong=(rootU**2+f2-(3*w-s.conjugate(kprime)*wprime)*s.conjugate(3*w-s.conjugate(kprime)*wprime)/(3*P*D0))/Y0**2
Vform=(f2-3*m*m*P/D0+rootU**2*(1-ep*(9+4*y*y)/(6*D0))
       +C*(3*s.cos(y+theta)+2*y*s.sin(y+theta))/D0)/Y0**2
check('strong_limit_full_y_potential',s.expand_complex(Vstrong),Vform)
Vaxis=(1-3*ep/2)*rootU**2+rho+3*C*s.cos(theta)
check('strong_limit_axis_potential',Vform.subs(y,0),Vaxis)
check('strong_limit_axis_transverse_tadpole',s.diff(Vform,y).subs(y,0),-C*s.sin(theta))
check('strong_limit_hidden_transverse_mass',s.diff(Vform,y,2).subs({y:0,rootU:0,rho:0})/F**2,2*m*m)
t,aa=s.symbols('t aa',real=True)
asym=s.series((Vform-rho-(1-3*ep/2)*rootU**2-3*C*s.cos(theta)).subs({y:aa*t,m:1/t}),t,0,1).removeO()
q=F*F*aa*aa-s.sqrt(2)*F*rootU*aa*s.sin(theta)
check('heavy_y_asymptotic_quadratic',asym,q)
astar=rootU*s.sin(theta)/(s.sqrt(2)*F)
check('heavy_y_local_minimum',s.diff(q,aa).subs(aa,astar))
check('heavy_y_tree_potential_shift',q.subs(aa,astar),-rootU**2*s.sin(theta)**2/2)

# Bilinear alone remains when d is sequestered; finite local kernel only.
ss,mu2,h=s.symbols('ss mu2 h',positive=True,real=True)
func=lambda x:x*x*(s.log(x/mu2)-s.Rational(3,2))
cw=(func(ss+h*s.sqrt(ss))+func(ss-h*s.sqrt(ss))-2*func(ss))/(32*s.pi**2)
coef=s.diff(cw,h,2).subs(h,0)/2
check('pure_B_CW_coefficient',coef,ss*s.log(ss/mu2)/(16*s.pi**2))
check('pure_B_CW_force_derivative',s.diff(coef,ss),(s.log(ss/mu2)+1)/(16*s.pi**2))

# High precision examples: solve both T components and tune f once, then add w.
mp.mp.dps=80
def refvac(kappa,phase,cc=mp.mpf('1'),rh=mp.mpf('1e-12')):
    cv=cc*mp.e**(1j*phase)/kappa
    et=rh/kappa**2
    aa,bb=mp.findroot(lambda aa,bb:(aa-mp.mpf('1.5')*aa*aa-bb*bb/2-3*cv.real-6*et*(1+2*aa),
                                   bb*(1-aa)-3*cv.imag),(3*cv.real,3*cv.imag),tol=mp.mpf('1e-65'))
    yy=1+2*aa
    fval=(1-2*aa)*(kappa*kappa*(aa*aa+bb*bb)/3+rh*yy)
    return aa,bb,fval

def potential(aa,bb,chi,yy,kappa,phase,fval,Fv=mp.mpf('.01'),uv=mp.mpf('1e-6')):
    Tv=mp.mpf('.5')+aa+1j*bb
    z=Fv*(chi+1j*yy)/mp.sqrt(2)
    Yv=2*Tv.real-Fv*Fv*yy*yy/3
    wzv=uv*mp.e**(-chi-1j*yy)
    wv=-Fv*wzv/mp.sqrt(2)
    jv=kappa*(Tv-mp.mpf('.5'))
    wtot=wv+mp.e**(1j*phase)+kappa*(Tv-mp.mpf('.5'))**2/2
    az=wzv-1j*mp.sqrt(2)*Fv*yy*jv/3
    return (abs(az)**2+fval+Yv*abs(jv)**2/3-2*(jv*mp.conj(wtot)).real)/Yv**2

examples=[]
for kval in (mp.mpf('1e3'),mp.mpf('1e5'),mp.mpf('1e7')):
    for phlabel,ph in (('aligned',mp.mpf('0')),('quadrature',mp.pi/2),('opposed',mp.pi)):
        av,bv,fv=refvac(kval,ph)
        vfun=lambda x,z:potential(x,z,mp.mpf('0'),mp.mpf('0'),kval,ph,fv,uv=mp.mpf('0'))
        numeric(f'vacuum_energy_{kval}_{phlabel}',vfun(av,bv),mp.mpf('1e-12'),mp.mpf('1e-60'))
        numeric(f'vacuum_grad_real_{kval}_{phlabel}',mp.diff(lambda x:vfun(x,bv),av),0,mp.mpf('1e-56'))
        numeric(f'vacuum_grad_imag_{kval}_{phlabel}',mp.diff(lambda x:vfun(av,x),bv),0,mp.mpf('1e-56'))
        # Local T response at a finite, small clock amplitude with fixed f.
        vclock=lambda x,z:potential(x,z,mp.mpf('0'),mp.mpf('0'),kval,ph,fv)
        at,bt=mp.findroot(lambda x,z:(mp.diff(lambda v:vclock(v,z),x),mp.diff(lambda v:vclock(x,v),z)),
                         (av,bv),tol=mp.mpf('1e-55'))
        Jt=kval*(at+1j*bt)
        Wt=mp.e**(1j*ph)-mp.mpf('.01e-6')/mp.sqrt(2)+kval*(at+1j*bt)**2/2
        H=mp.matrix([[mp.diff(lambda x:vfun(x,bv),av,2),mp.diff(lambda x:mp.diff(lambda z:vfun(x,z),bv),av)],
                     [mp.diff(lambda x:mp.diff(lambda z:vfun(x,z),bv),av),mp.diff(lambda z:vfun(av,z),bv,2)]])
        eig=mp.eigsy(H,eigvals_only=True)
        checks.append(dict(name=f'vacuum_Hessian_positive_{kval}_{phlabel}',passed=bool(min(eig)>0),
                           min_eigenvalue=mp.nstr(min(eig),12)))
        rec=dict(kappa=mp.nstr(kval,12),phase=phlabel,a=mp.nstr(av,20),b=mp.nstr(bv,20),
                 f_squared=mp.nstr(fv,20),ReJ_reference=mp.nstr(kval*av,20),
                 J_to_3W_at_clock_relative_error=mp.nstr(abs(Jt-3*Wt/(1+2*at))/abs(Jt),12),
                 physical_modulus_min_mass_to_kappa_over3=mp.nstr(mp.sqrt(min(eig)*(1+2*av)**2/6)/(kval/3),12))
        if phlabel=='quadrature':
            rec['ReJ_over_leading_4p5_over_kappa']=mp.nstr(kval*av/(mp.mpf('4.5')/kval),16)
        examples.append(rec)

out=dict(scope='Specified sequestered no-scale K with quadratic W stabilization; symbolic identities and local dimensionless checks; no cosmological trajectories or observations',
         python=platform.python_version(),sympy=s.__version__,mpmath=mp.__version__,mpmath_dps=80,
         script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
         checks=checks,passed=sum(x['passed'] for x in checks),total=len(checks),numeric_examples=examples,
         sources=[dict(url='https://arxiv.org/pdf/1012.1858',scope='Abstract, introduction, sections 2.1–2.2; separable Kahler geometry does not automatically suppress B for explicit visible supersymmetric masses; accessed 2026-09-25'),
                  dict(url='https://arxiv.org/html/1209.0499v2',scope='Section 2.1, equations 9–12 and section 2.2 introductory caveat; strongly stabilized modulus can have small D_T W while W_T remains nonzero; different uplift K from our model'),
                  dict(url='https://arxiv.org/html/hep-ph/9707209v3',scope='Opened for general soft-term conventions; no formula or numerical result imported into this calculation')],
         anomalies=[dict(kind='source_transport',url='https://arxiv.org/html/1012.1858v2',status='HTTP 404',resolution='Read authoritative arXiv PDF instead; no claimed read of unavailable version'),
                    dict(kind='check_implementation',attempt=1,passed=64,total=65,
                         failed_check='strong_limit_F_T_zero',
                         cause='Contravariant auxiliary vector used inverse metric without the transpose needed by the stored holomorphic/antiholomorphic matrix convention',
                         resolution='Changed this diagnostic contraction to Gi.T * Db. Potential, spectrum, stationary equations, and numerical examples were unchanged; reran full bounded script')])
OUT.write_text(json.dumps(out,indent=2,ensure_ascii=False)+'\n')
print(json.dumps({k:out[k] for k in ('passed','total')}))
if not all(x['passed'] for x in checks):
    print(json.dumps([x for x in checks if not x['passed']],indent=2))
    raise SystemExit(1)
