#!/usr/bin/env python3
"""Independent algebra audit of the specified nilpotent hidden-sector model.

No integration, observational fit, or test of physical truth. The charged
quadratic action is expanded independently of the cosmological implementation.
"""
from pathlib import Path
import hashlib
import json
import platform
import sympy as s

OUT = Path(__file__).with_suffix('.json')
checks = []

def check(name, lhs, rhs=0):
    residual = s.trigsimp(s.simplify(s.expand(lhs-rhs)))
    checks.append(dict(name=name, passed=residual == 0, residual=str(residual)))

F, P, u, m, rho = s.symbols('F P u m rho', positive=True, real=True)
y, theta = s.symbols('y theta', real=True)
ep = F**2/P
b = 1-3*ep/2
C = s.sqrt(2)*F*m*u
phase = theta+y
w = -F*u/s.sqrt(2)*s.exp(-s.I*y)+m*P*s.exp(s.I*theta)
A = u*s.exp(-s.I*y)-s.I*s.sqrt(2)*F*y*w/P
f2 = 3*m*m*P+rho
raw = s.expand_complex(A*s.conjugate(A)+f2-3*w*s.conjugate(w)/P)
inside = u*u*(b+ep*ep*y*y)+rho+2*F*F*m*m*y*y+C*((3-2*ep*y*y)*s.cos(phase)+2*y*s.sin(phase))
check('full_F_potential_direct_DW', raw, inside)
V = s.exp(ep*y*y)*inside
# d/dchi acts as -u*d/du because sqrt(U) is exponential.
Dchi = lambda expr: -u*s.diff(expr,u)
Vaxis = b*u*u+rho+3*C*s.cos(theta)
check('axis_potential',V.subs(y,0),Vaxis)
check('axis_longitudinal_force',Dchi(V).subs(y,0),-2*b*u*u-3*C*s.cos(theta))
check('axis_transverse_tadpole',s.diff(V,y).subs(y,0),-C*s.sin(theta))
check('axis_longitudinal_hessian',Dchi(Dchi(V)).subs(y,0),4*b*u*u+3*C*s.cos(theta))
check('axis_mixed_hessian',Dchi(s.diff(V,y)).subs(y,0),C*s.sin(theta))
Vyy = 2*ep*(1-ep/2)*u*u+2*ep*rho+4*F*F*m*m+(1+2*ep)*C*s.cos(theta)
check('axis_transverse_hessian',s.diff(V,y,2).subs(y,0),Vyy)
for label,angle,co,si in [('aligned',0,1,0),('orthogonal',s.pi/2,0,1),('opposed',s.pi,-1,0)]:
    check(label+'_longitudinal_force',Dchi(V).subs({y:0,theta:angle}),-2*b*u*u-3*C*co)
    check(label+'_transverse_force',s.diff(V,y).subs({y:0,theta:angle}),-C*si)
check('zero_Wc_still_has_hidden_Lambda_y_dependence',V.subs(m,0),s.exp(ep*y*y)*(u*u*(b+ep*ep*y*y)+rho))
check('cross_gradient_norm_phase_bound_identity',9*C*C*s.cos(theta)**2+C*C*s.sin(theta)**2,C*C*(1+8*s.cos(theta)**2))

# Direct charged-field expansion includes |D_X W|^2=f^2 BEFORE setting X=0.
qp,qm,qpb,qmb,t=s.symbols('qp qm qpb qmb t')
ww,wwb,mm,mmb,Az,Azb,Cz,Czb,ff=s.symbols('w wb M Mb A Ab Cz Czb ff')
norm=qp*qpb+qm*qmb
v0=Az*Azb+ff-3*ww*wwb/P
Dz=Az+t*t*Cz*qp*qm
Dzb=Azb+t*t*Czb*qpb*qmb
Dp=t*(mm*qm+qpb*ww/P)
Dpb=t*(mmb*qmb+qp*wwb/P)
Dm=t*(mm*qp+qmb*ww/P)
Dmb=t*(mmb*qpb+qm*wwb/P)
WV=ww+t*t*mm*qp*qm
WbV=wwb+t*t*mmb*qpb*qmb
poly=s.expand((1+t*t*norm/P)*(Dz*Dzb+ff+Dp*Dpb+Dm*Dmb-3*WV*WbV/P))
quad=s.expand(poly).coeff(t,2)
diag=mm*mmb+ww*wwb/P**2+v0/P
B=Azb*Cz-wwb*mm/P
Bb=Az*Czb-ww*mmb/P
check('charged_action_direct_expansion_with_hidden_F',quad,diag*norm+B*qp*qm+Bb*qpb*qmb)
check('hidden_F_enters_diagonal',s.diff(diag,ff),1/P)
check('constant_hidden_F_no_direct_B',s.diff(B,ff),0)
M,r=s.symbols('M r',real=True)
wa=w.subs(y,0)
m32=s.expand_complex(wa*s.conjugate(wa)/P**2)
dexpected=m*m+(1-ep)*u*u/P+rho/P+2*C*s.cos(theta)/P
check('physical_gravitino_mass',m32,m*m+F*F*u*u/(2*P*P)-C*s.cos(theta)/P)
check('charged_axis_diagonal',m32+Vaxis/P,dexpected)
Baxis=u*s.sqrt(2)*M*r/F-s.conjugate(wa)*M/P
Bexpected=M*(s.sqrt(2)*u/F*(r+ep/2)-m*s.exp(-s.I*theta))
check('charged_axis_complex_B',Baxis,Bexpected)
check('charged_B_modulus_squared',s.expand_complex(Baxis*s.conjugate(Baxis)),M*M*((s.sqrt(2)*u/F*(r+ep/2))**2+m*m-2*m*s.sqrt(2)*u/F*(r+ep/2)*s.cos(theta)))
check('nilpotent_F_positive_modulus',s.exp(ep*y*y)*f2,s.exp(ep*y*y)*(3*m*m*P+rho))

# Controlled asymptotic tree argument: y=a/m at m >> sqrt(U)/F.
# This is not a cosmological integration of the fast heavy mode.
a,delta=s.symbols('a delta',real=True)
asym=s.series((V-rho-3*C*s.cos(theta)-b*u*u).subs({y:a*delta,m:1/delta}),delta,0,1).removeO()
q=2*F*F*a*a-s.sqrt(2)*F*u*a*s.sin(theta)
check('heavy_y_asymptotic_quadratic',asym,q)
astar=u*s.sin(theta)/(2*s.sqrt(2)*F)
check('heavy_y_stationary_shift',s.diff(q,a).subs(a,astar))
check('heavy_y_finite_tree_correction',q.subs(a,astar),-u*u*s.sin(theta)**2/4)

# Leading finite local CW kernel: d=m^2 and |B|^2=m^2 s.
masssq,mu2,hid=s.symbols('masssq mu2 hid',positive=True,real=True)
func=lambda x: x*x*(s.log(x/mu2)-s.Rational(3,2))
cw=(func(masssq+hid*hid+hid*s.sqrt(masssq))+func(masssq+hid*hid-hid*s.sqrt(masssq))-2*func(masssq))/(32*s.pi**2)
coef=s.diff(cw,hid,2).subs(hid,0)/2
check('CW_hidden_m2_coefficient',coef,masssq*(3*s.log(masssq/mu2)-2)/(16*s.pi**2))
check('CW_hidden_force_masssq_derivative',s.diff(coef,masssq),(3*s.log(masssq/mu2)+1)/(16*s.pi**2))
check('CW_hidden_force_at_initial_scale_is_nonzero',s.diff(coef,masssq).subs(masssq,mu2),1/(16*s.pi**2))

out=dict(scope='specified conditional model: algebra identities only',
         python=platform.python_version(),sympy=s.__version__,
         script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
         checks=checks,passed=sum(x['passed'] for x in checks),total=len(checks))
OUT.write_text(json.dumps(out,indent=2,ensure_ascii=False)+'\n')
print(json.dumps({k:out[k] for k in ('passed','total')}))
if not all(x['passed'] for x in checks):
    raise SystemExit(1)
