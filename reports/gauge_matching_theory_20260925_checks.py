#!/usr/bin/env python3
"""One-loop U(1) heavy-pair matching identities, not an observational test.

KL denotes Kaplunovsky--Louis hep-th/9402005v2, equations 2.21,
3.3--3.7, 3.13, 3.19 and 3.26--3.29. No previous matching code is imported.
The tiny-splitting identity is a component one-loop threshold specialization,
not an assertion of exact supersymmetric holomorphy at finite breaking.
"""
from datetime import datetime, timezone
from pathlib import Path
import hashlib
import json
import sympy as S

OUT = Path(__file__).with_suffix('.json')
checks = []


def check(name, lhs, rhs=0):
    residual = S.simplify(S.expand(lhs-rhs))
    checks.append(dict(name=name, passed=bool(residual == 0), residual=str(residual)))


p2=S.pi**2
b,c,tg=S.symbols('b c tg', real=True)
f,k,z,l,q=S.symbols('f k log_Zprod log_Mhol_squared log_Lambda_over_mu', real=True)
g,rg=S.symbols('g_inverse_squared log_g_inverse_squared', real=True)
# KL F=g^-2-TG/(8pi^2)log(g^-2)+sum(T logZ)/(8pi^2)+b/(8pi^2)log(mu).
check('KL_rearrangement_signs',
      g-tg*rg/(8*p2)+z/(8*p2)-(f+c*k/(16*p2)+b*q/(8*p2)),
      g-(f+c*k/(16*p2)+b*q/(8*p2)+tg*rg/(8*p2)-z/(8*p2)))
charges=[S.Integer(1),S.Integer(-1)]
bb=sum(x*x for x in charges)-3*0
cc=sum(x*x for x in charges)-0
check('U1_pair_b_coefficient',bb,2)
check('U1_pair_c_coefficient',cc,2)
check('component_beta_sum',S.Rational(1,3)*2+S.Rational(4,3),bb)

log_s=k+l-z
uv=f+cc*k/(16*p2)-z/(8*p2)  # field-independent RG logarithm drops in differences
local=-bb*log_s/(16*p2)
ir=uv+local
check('K_and_matter_metric_cancel_in_IR',ir,f-l/(8*p2))
check('KL_c_jump_equals_nonHiggs_index_sum',cc-0,bb)
check('holomorphic_threshold_real_part',-2*l/2/(8*p2),-l/(8*p2))

j,u=S.symbols('Re_J_over_Mp_squared log_abs_Uplus_Uminus',real=True)
change={k:k+2*j,l:l-2*j,f:f-cc*j/(8*p2)}
check('Kahler_physical_mass_invariant',log_s.subs(change,simultaneous=True),log_s)
check('Kahler_UV_g_invariant',uv.subs(change,simultaneous=True),uv)
check('Kahler_IR_g_invariant',ir.subs(change,simultaneous=True),ir)
konishi={z:z-2*u,l:l-2*u,f:f-u/(4*p2)}
check('Konishi_physical_mass_invariant',log_s.subs(konishi,simultaneous=True),log_s)
check('Konishi_UV_g_invariant',uv.subs(konishi,simultaneous=True),uv)
check('Konishi_IR_g_invariant',ir.subs(konishi,simultaneous=True),ir)

# B fixes physical UV g at each real-axis point. It is a different model boundary.
physical_boundary=ir.subs(f, -k/(8*p2)+z/(8*p2))
check('fixed_physical_UV_leaves_local_threshold',physical_boundary,local)
x,y,mp=S.symbols('x y Mp', real=True, nonzero=True)
kcan=(x*x+y*y)/(mp*mp)
ref_compensator=-(x*x-y*y)/(8*p2*mp*mp)
residual=ref_compensator+kcan/(8*p2)
check('canonical_holomorphic_f_compensates_only_axis',residual,y*y/(4*p2*mp*mp))
check('canonical_K_nonharmonic_off_axis',S.diff(kcan,x,2)+S.diff(kcan,y,2),4/mp**2)

ss,dd,hh=S.symbols('s d h', positive=True)
r,v,t=S.symbols('r v t',real=True)
det=(ss+dd+hh)*(ss+dd-hh)
check('split_scalar_determinant',det,ss**2*((1+dd/ss)**2-(hh/ss)**2))
check('split_threshold_log_s_coefficient',2*S.Rational(1,3)+S.Rational(4,3),2)
R=(1+r)**2-v*v
series=S.series(S.log(R.subs({r:t*r,v:t*v})),t,0,4).removeO()
check('soft_log_determinant_small_split_series',series,
      2*r*t-(r*r+v*v)*t*t+(S.Rational(2,3)*r**3+2*r*v*v)*t**3)
rp,vp=S.symbols('rprime vprime', real=True)
check('soft_drift_derivative_includes_both_splittings',
      S.diff(S.log(R),r)*rp+S.diff(S.log(R),v)*vp,
      (2*(1+r)*rp-2*v*vp)/R)

alpha0,D=S.symbols('alpha0 delta_inverse_g_squared',real=True)
relative=1/(1+4*S.pi*alpha0*D)-1
check('alpha_inverse_conversion_and_today_normalization',
      relative,-4*S.pi*alpha0*D/(1+4*S.pi*alpha0*D))
xi,ch,ch0=S.symbols('xi chi chi0',positive=True)
Mratio=(1+xi/ch**2)/(1+xi/ch0**2)
check('inverse_square_mass_gives_log_response_derivative',
      S.diff(S.log(Mratio),ch),-2*xi/(ch*(ch**2+xi)))

out=dict(
    created_utc=datetime.now(timezone.utc).isoformat(),
    scope='Algebra identities of the specified one-loop minimal U(1); no observed alpha, no full SM or SUGRA completion.',
    source='https://arxiv.org/pdf/hep-th/9402005v2',
    code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    sympy_version=S.__version__,checks=checks,
    passed=sum(x['passed'] for x in checks),total=len(checks))
OUT.write_text(json.dumps(out,indent=2,ensure_ascii=False)+'\n')
print(json.dumps({k:out[k] for k in ('passed','total')}))
if not all(x['passed'] for x in checks):
    raise SystemExit(1)
