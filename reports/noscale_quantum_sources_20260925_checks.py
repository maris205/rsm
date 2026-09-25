#!/usr/bin/env python3
"""Independent symbolic checks of no-scale soft/gauge-interface identities.
No production imports, no observations, no loop-completeness assertion.
M_P=1, F^i=-exp(K/2) K^{i jbar} D_jbar Wbar.
"""
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
import sympy as s

checks=[]
def check(name, difference):
    v=s.factor(s.simplify(difference))
    checks.append({'name':name,'residual':str(v),'passed':v==0})

t,tb,z,zb=s.symbols('t tb z zb')
w,wb,wt,wbt,wz,wbz=s.symbols('w wb wt wbt wz wbz')
Y=t+tb+(z-zb)**2/6
K=-3*s.log(Y)
qs=[t,z];qbs=[tb,zb]
g=s.Matrix([[s.diff(K,i,j) for j in qbs] for i in qs])
ki=s.Matrix([s.diff(K,i) for i in qs]);kb=s.Matrix([s.diff(K,j) for j in qbs])
inv=g.inv().T
for n, item in enumerate((ki.T*inv+s.Matrix([[Y,0]]))):
    check('ideal_K_inverse_contraction_'+str(n),item)
check('ideal_no_scale_norm',(ki.T*inv*kb)[0]-3)
F=-Y**s.Rational(-3,2)*inv*(s.Matrix([wbt,wbz])+kb*wb)
fc=Y**s.Rational(-3,2)*wb+(ki.T*F)[0]/3
check('ideal_Fcomp_general_WT',fc-Y**s.Rational(-1,2)*wbt/3)
check('ideal_Fcomp_WT_zero',fc.subs(wbt,0))
check('ideal_FZ_real_axis',F[1].subs(zb,z)+wbz/s.sqrt(t+tb))
# Universal untwisted metrics: product Z+ Z- = exp(2K/3).
fiK,fiM,mbar,V,m2=s.symbols('fiK fiM mbar V m2')
b_over_M=-fiM-fiK/3-mbar
check('B_is_minus_Fcomp_minus_clock_derivative',b_over_M+fiM+(mbar+fiK/3))
check('ideal_B_has_no_m32',b_over_M.subs(fiK,-3*mbar)+fiM)
check('strong_superpotential_FT_zero_B_residual',b_over_M.subs({fiK:0,fiM:0})+mbar)
check('common_soft_all_F_included',m2+V-(V+3*m2)/3-2*V/3)
# General gauge group, D'Eramo/Thaler/Thomas Eq.30--32 combined.
TG,TR,g2=s.symbols('TG TR g2')
mg=-g2/(16*s.pi**2)*((3*TG-TR)*mbar+(TG-TR)*fiK+2*TR*fiK/3)
check('gauge_total_reduces_to_Fcomp',mg+g2/(16*s.pi**2)*(3*TG-TR)*(mbar+fiK/3))
check('gauge_noscale_total_zero',mg.subs(fiK,-3*mbar))
check('gauge_U1_pair_coefficient',mg.subs({TG:0,TR:2})-2*g2/(16*s.pi**2)*(mbar+fiK/3))
check('gauge_empty_IR_beta_zero',mg.subs({TG:0,TR:0}))
# A general curved G(T,Tbar), on the shift-clock real axis.
y,gt,gb,gtt=s.symbols('Y G_T G_barT G_TbarT',positive=True)
A=gt*gb-y*gtt
KT=-3*gt/y; KTbar=-3*gb/y;KTT=3*A/y**2
FT=-y**s.Rational(-3,2)/KTT*(wbt+KTbar*wb)
FC=y**s.Rational(-3,2)*wb+KT*FT/3
check('curved_G_Fcomp_formula',FC-y**s.Rational(-1,2)/A*(gt*wbt/3-gtt*wb))
check('curved_G_WT_zero_formula',FC.subs(wbt,0)+y**s.Rational(-1,2)*gtt*wb/A)
check('curved_G_center_zero',FC.subs({wbt:0,gtt:0}))
# Quartic derivative jets at T=Tbar=1/2, not a solution claim.
a,b=s.symbols('a b');G=t+tb+a*(t+tb-1)**4+b*(t-tb)**4
c={t:s.Rational(1,2),tb:s.Rational(1,2)}
check('quartic_center_G',G.subs(c)-1)
check('quartic_center_G_T',s.diff(G,t).subs(c)-1)
check('quartic_center_G_TbarT',s.diff(G,t,tb).subs(c))
check('quartic_center_third_mixed',s.diff(G,t,t,tb).subs(c))
r,eta=s.symbols('r eta',real=True)
check('quartic_mixed_curvature_real_axis',s.diff(G,t,tb).subs({t:r/2,tb:r/2})-12*a*(r-1)**2)
check('quartic_mixed_curvature_imag_axis',s.diff(G,t,tb).subs({t:s.Rational(1,2)+s.I*eta,tb:s.Rational(1,2)-s.I*eta})-48*b*eta**2)
# Positive rolling energy exerts modulus force even at the cancellation point.
U=s.symbols('U',positive=True)
check('positive_clock_energy_modulus_force',s.diff(U/G**2,t).subs(c)+2*U)

here=Path(__file__).resolve()
out={'run_utc':datetime.now(timezone.utc).isoformat(),'scope':__doc__,
     'sympy_version':s.__version__,'script_sha256':hashlib.sha256(here.read_bytes()).hexdigest(),
     'total':len(checks),'passed':sum(c['passed'] for c in checks),'checks':checks,
     'limitations':['Symbolic local identities only; no cosmological trajectory is solved.',
       'No UV matching, complete SUGRA loops, or observational claim.',
       'Quartic central derivative jets do not prove a stationary point with U>0.']}
here.with_suffix('.json').write_text(json.dumps(out,indent=2,ensure_ascii=False)+'\n')
print(json.dumps({'total':out['total'],'passed':out['passed']}))
if out['passed']!=out['total']: raise SystemExit(1)
