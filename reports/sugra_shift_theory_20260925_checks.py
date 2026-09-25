#!/usr/bin/env python3
"""Independent algebra identities for the specified one-clock SUGRA candidate.

No cosmological integration and no observational fit. Uses formal conjugate
symbols when expanding the charged quadratic action.
"""
from pathlib import Path
import json
import sympy as s

OUT = Path(__file__).with_suffix('.json')
checks = []

def check(name, lhs, rhs=0):
    residual = s.simplify(s.expand(lhs-rhs))
    checks.append(dict(name=name, passed=residual == 0, residual=str(residual)))

Z, Zb, Mp, F, chi, y, ci, Ui, xi, mi = s.symbols(
    'Z Zb Mp F chi y ci Ui xi mi', positive=True, real=True)
ep = F**2/Mp**2
Kcan=Z*Zb
Ksh=-(Z-Zb)**2/2
check('shift_metric_is_one',s.diff(Ksh,Z,Zb),1)
check('canonical_metric_is_one',s.diff(Kcan,Z,Zb),1)
check('K_difference_is_holomorphic_plus_conjugate',Ksh-Kcan,-(Z**2+Zb**2)/2)
sub={Z:F*(chi+s.I*y)/s.sqrt(2),Zb:F*(chi-s.I*y)/s.sqrt(2)}
check('shift_K_in_real_coordinates',Ksh.subs(sub),F**2*y**2)
check('canonical_K_in_real_coordinates',Kcan.subs(sub),F**2*(chi**2+y**2)/2)

w=-F*s.sqrt(Ui)/s.sqrt(2)*s.exp(-s.sqrt(2)*(Z-F*ci/s.sqrt(2))/F)
wb=w.xreplace({Z:Zb})
U=Ui*s.exp(-2*(chi-ci))
pot={}
for name,K in [('shift',Ksh),('canonical',Kcan)]:
    A=s.diff(w,Z)+s.diff(K,Z)*w/Mp**2
    Ab=s.diff(wb,Zb)+s.diff(K,Zb)*wb/Mp**2
    V=s.simplify((s.exp(K/Mp**2)*(A*Ab-3*w*wb/Mp**2)).subs(sub))
    pot[name]=V

Vsh=U*s.exp(ep*y**2)*(1-3*ep/2+ep**2*y**2)
Vca=U*s.exp(ep*(chi**2+y**2)/2)*((1-ep*chi/2)**2+ep**2*y**2/4-3*ep/2)
check('shift_full_potential_from_DW',pot['shift'],Vsh)
check('canonical_full_potential_from_DW',pot['canonical'],Vca)
check('shift_longitudinal_force',s.diff(Vsh,chi),-2*Vsh)
check('shift_transverse_gradient',s.diff(Vsh,y),
      2*ep*y*U*s.exp(ep*y**2)*(1-ep/2+ep**2*y**2))
check('shift_axis_cross_hessian_zero',s.diff(Vsh,chi,y).subs(y,0))
check('shift_axis_transverse_mass',s.diff(Vsh,y,2).subs(y,0)/F**2,
      2*U/Mp**2*(1-ep/2))
check('canonical_axis_transverse_mass',s.diff(Vca,y,2).subs(y,0)/F**2,
      U*s.exp(ep*chi**2/2)/Mp**2*((1-ep*chi/2)**2-ep))

# Direct quadratic expansion in Q_+, Q_- and independent conjugates.
qp,qm,qpb,qmb,t=s.symbols('qp qm qpb qmb t')
ww,wwb,mm,mmb,Az,Azb,C,Cb,kk=s.symbols('w wb M Mb A Ab C Cb k')
norm=qp*qpb+qm*qmb
v0=Az*Azb-3*ww*wwb/Mp**2
Dz=Az+t**2*C*qp*qm
Dzb=Azb+t**2*Cb*qpb*qmb
Dp=t*(mm*qm+qpb*ww/Mp**2)
Dpb=t*(mmb*qmb+qp*wwb/Mp**2)
Dm=t*(mm*qp+qmb*ww/Mp**2)
Dmb=t*(mmb*qpb+qm*wwb/Mp**2)
WV=ww+t**2*mm*qp*qm
WbV=wwb+t**2*mmb*qpb*qmb
poly=s.expand((1+t**2*norm/Mp**2)*(Dz*Dzb+Dp*Dpb+Dm*Dmb-3*WV*WbV/Mp**2))
quad=s.expand(poly).coeff(t,2)
diag=mm*mmb+ww*wwb/Mp**4+v0/Mp**2
B=Azb*C-wwb*mm/Mp**2
Bb=Az*Cb-ww*mmb/Mp**2
check('charged_direct_quadratic_expansion',quad,diag*norm+B*qp*qm+Bb*qpb*qmb)
check('charged_diagonal_V_plus_m32',s.diff(quad,qp,qpb),diag)
check('charged_B_term_sign',s.diff(quad,qp,qm),B)

M=mi*s.sqrt(1+xi/(s.sqrt(2)*Z/F)**2)
Mchi=M.subs(Z,F*chi/s.sqrt(2))
r=-xi/(chi*(chi**2+xi))
check('log_holomorphic_mass_derivative',s.diff(Mchi,chi)/Mchi,r)
for name,K in [('shift',Ksh),('canonical',Kcan)]:
    A=s.diff(w,Z)+s.diff(K,Z)*w/Mp**2
    Ab=s.diff(wb,Zb)+s.diff(K,Zb)*wb/Mp**2
    Bfull=s.exp(K/Mp**2)*(Ab*(s.diff(M,Z)+s.diff(K,Z)*M/Mp**2)-wb*M/Mp**2)
    Baxis=s.simplify(Bfull.subs(sub).subs(y,0))
    m32=(s.exp(K/Mp**2)*w*wb/Mp**4).subs(sub).subs(y,0)
    d=s.simplify(m32+pot[name].subs(y,0)/Mp**2)
    if name=='shift':
        Bexpected=s.sqrt(2*U)*Mchi/F*(r+ep/2)
        dexpected=(1-ep)*U/Mp**2
    else:
        a=1-ep*chi/2
        Bexpected=s.exp(ep*chi**2/2)*s.sqrt(2*U)*Mchi/F*(a*(r+ep*chi/2)+ep/2)
        dexpected=s.exp(ep*chi**2/2)*U/Mp**2*(a*a-ep)
    check(name+'_charged_B',Baxis,Bexpected)
    check(name+'_charged_gravity_diagonal',d,dexpected)

check('shift_true_Kahler_transform_preserves_G',
      (Ksh-Kcan)/Mp**2+(Z**2+Zb**2)/(2*Mp**2))
check('canonical_physical_log_mass_squared_derivative',
      s.diff(s.log(s.exp(ep*chi**2/2)*Mchi**2),chi),ep*chi+2*r)
check('shift_physical_log_mass_squared_derivative',s.diff(s.log(Mchi**2),chi),2*r)

out=dict(scope='algebra identities, not a test of physical validity',
         checks=checks, passed=sum(x['passed'] for x in checks), total=len(checks))
OUT.write_text(json.dumps(out,indent=2,ensure_ascii=False)+'\n')
print(json.dumps({k:out[k] for k in ('passed','total')}))
if not all(x['passed'] for x in checks):
    raise SystemExit(1)
