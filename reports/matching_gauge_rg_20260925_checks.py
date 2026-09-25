#!/usr/bin/env python3
"""Independent threshold/RG algebra; no main-experiment imports or fits."""
from __future__ import annotations
import hashlib
import json
import math
from pathlib import Path
import sys
import mpmath as mp
import sympy as sp

ROOT=Path(__file__).resolve().parents[1]
OUT=Path(__file__).with_suffix('.json')
checks=[]
def eq(name, expression):
    residual=sp.simplify(sp.expand(expression))
    checks.append({'name':name,'pass':residual==0,'residual':str(residual)})
def cond(name, passed, detail):
    checks.append({'name':name,'pass':bool(passed),'detail':detail})

s,d,h,mu,mu0=sp.symbols('s d h mu mu0',positive=True)
x,csub=sp.symbols('x csub',real=True)
F=lambda x: x**2*(sp.log(x/mu**2)-sp.Rational(3,2))
V=(F(s+d+h)+F(s+d-h)-2*F(s))/(32*sp.pi**2)
J=(2*s*d+d*d+h*h)/(8*sp.pi**2)
eq('exact_CW_scale_derivative',mu*sp.diff(V,mu)+J)
# Independent soft-matrix contraction (Martin Eq7.11), omitting gaugino.
soft=sp.diag(d,d)
mass=sp.Matrix([[0,sp.sqrt(s)],[sp.sqrt(s),0]])
bmat=sp.Matrix([[0,h],[h,0]])
Jmatrix=(sp.trace(soft*soft)+2*sp.trace(soft*mass*mass)+sp.trace(bmat*bmat))/(16*sp.pi**2)
eq('vacuum_beta_matrix_vs_mass_eigenvalues',Jmatrix-J)
eq('RG_completed_potential',mu*sp.diff(V+J*sp.log(mu/mu0),mu))
for p in (s,d,h):
    eq('RG_commutes_with_'+str(p)+'_derivative',mu*sp.diff(sp.diff(V+J*sp.log(mu/mu0),p),mu))
# Finite subtraction: replacing 3/2 by 3/2+delta_c.
delta_c=sp.symbols('delta_c')
finite_shift=-delta_c*((s+d+h)**2+(s+d-h)**2-2*s*s)/(32*sp.pi**2)
eq('finite_subtraction_local_polynomial',finite_shift+delta_c*J/2)
q=sp.symbols('q')
series=sp.series(V.subs({d:q*d,h:q*h}),q,0,4).removeO().expand()
expected=q*s*d*(sp.log(s/mu**2)-1)/(8*sp.pi**2)+q*q*(d*d+h*h)*sp.log(s/mu**2)/(16*sp.pi**2)+q**3*(d**3/(48*sp.pi**2*s)+d*h*h/(16*sp.pi**2*s))
eq('CW_small_splitting_cubic_series',sp.expand_log(series-expected,force=True))
mt,mv=sp.symbols('mt mv',positive=True)
VT=(F(mt**2)+F(mv**2))/(64*sp.pi**2)
JT=(mt**4+mv**4)/(32*sp.pi**2)
eq('two_real_scalar_scale_derivative',mu*sp.diff(VT,mu)+JT)
# Convention-fixed anomaly + superpotential threshold.
k,Fc,AM,Ap,Am,UV=sp.symbols('k Fcomp A_M A_plus A_minus UV')
AS=Ap+Am
S=Fc+AM-AS
above=UV+2*k*(Fc-AS)
threshold=-2*k*S
below=UV-2*k*AM
eq('above_plus_threshold_equals_IR',above+threshold-below)
eq('constant_holomorphic_mass_decoupling',(above+threshold-UV).subs(AM,0))
eq('arbitrary_linear_metric_Aplus_cancels',sp.diff(above+threshold,Ap))
eq('arbitrary_linear_metric_Aminus_cancels',sp.diff(above+threshold,Am))
eq('compensator_cancels_below_empty_U1',sp.diff(above+threshold,Fc))
cond('variable_mass_not_zero_when_betaIR_zero',sp.simplify((above+threshold-UV).subs({Fc:0,Ap:0,Am:0}))!=0,'Residual is -2 k A_M for A_M != 0.')
cond('missing_threshold_negative_control',sp.simplify((above-UV).subs({AM:0,Ap:0,Am:0}))!=0,'Above-threshold 2 k Fcomp would falsely survive constant-mass decoupling.')
# Holomorphic f and field redefinition, derivative represented abstractly.
f,m,j,u=sp.symbols('f m j u')
coef=sp.Rational(1,4)/sp.pi**2
fIR=f-coef*m # m = log M
for name, fp, mp_ in [('Kahler',f-coef*j,m-j),('Konishi_pair',f-coef*u,m-u)]:
    eq(name+'_holomorphic_matching_invariance',fp-coef*mp_-fIR)
# EFT stage's selected d and B basis, keeping sd and B^2 only.
p,X,c,b2,L,nu,ell=sp.symbols('p X c beta_squared L nu ell')
A=-c+b2/2
selected=p*(A*L+c+nu)/(8*sp.pi**2)
eq('selected_basis_scale_completion',sp.diff(selected.subs({L:L-2*ell,nu:nu+2*A*ell}),ell))
Kc=sp.diff(selected,c);Kb=sp.diff(selected,b2);Kn=sp.diff(selected,nu)
eq('selected_basis_rank_two_identity',Kc+2*Kb-Kn)
eq('selected_basis_finite_degeneracy',selected.subs({c:c+ell,b2:b2+2*ell,nu:nu-ell},simultaneous=True)-selected)
# At frozen background the omitted d^2 term is a separate order in X/s.
leading= s*d*(L-1)/(8*sp.pi**2)+h*h*L/(16*sp.pi**2)
eq('selected_basis_from_general_CW',leading.subs({d:-c*X,h*h:b2*s*X}).subs(s*X,p)-selected.subs(nu,0))
eq('selected_basis_RG_from_general_J',J.subs({d:-c*X,h*h:b2*s*X})-c*c*X*X/(8*sp.pi**2)-2*A*s*X/(8*sp.pi**2))

# Direct arbitrary-precision determinant, avoiding shared Taylor implementations.
mp.mp.dps=100
samples=[('3.2','.07','.11'),('1e22','1e-54','1e-21'),('9','.4','0')]
max_scale=mp.mpf(0)
for idx,row in enumerate(samples):
    ss,dd,hh=map(mp.mpf,row)
    def v(mm):
        fn=lambda xx:xx*xx*(mp.log(xx/(mm*mm))-mp.mpf('1.5'))
        return (fn(ss+dd+hh)+fn(ss+dd-hh)-2*fn(ss))/(32*mp.pi**2)
    jj=(2*ss*dd+dd*dd+hh*hh)/(8*mp.pi**2)
    refmu=mp.sqrt(ss)
    ref=v(refmu)
    for factor in ('0.5','2','7'):
        mm=refmu*mp.mpf(factor)
        corrected=v(mm)+jj*mp.log(mm/refmu)
        err=abs(corrected-ref)/max(abs(ref),abs(jj),mp.mpf('1e-90'))
        max_scale=max(max_scale,err)
        cond(f'direct_high_precision_RG_sample_{idx}_factor_{factor}',err<mp.mpf('1e-18'),str(err))

inp=json.loads((ROOT/'experiments/noscale_protection_v01/results/inputs.json').read_text())['parent_inputs']
Mq=inp['holomorphic_mass_i_eV'];MT=math.sqrt(48.)
ratio=Mq/MT
cond('threshold_hierarchy_Q_above_T',ratio>1e10,{'M_Q_eV':Mq,'m_T_eV_for_mG1_a1':MT,'ratio':ratio})
Fclock=math.sqrt(inp['F_squared_eV2']);U=inp['potential_unit_eV4']*inp['W_i'];chi=inp['chi_i'];xi=inp['xi']
AMnum=math.sqrt(2*U)/Fclock*xi/(chi*(chi*chi+xi))
g2=4*math.pi*inp['alpha_reference']
MIR=-2*g2/(16*math.pi**2)*AMnum
cond('clock_mass_derivative_nonzero',AMnum>0,{'A_M_eV_at_initial_center':AMnum,'conditional_IR_gaugino_eV':MIR})
result={'scope':'Independent algebraic and high-precision checks of a specified leading threshold and selected CW RG identities; not full SUGRA matching.','python':sys.version,'sympy':sp.__version__,'mpmath':mp.__version__,'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'parent_input_sha256':hashlib.sha256((ROOT/'experiments/noscale_protection_v01/results/inputs.json').read_bytes()).hexdigest(),'checks':checks,'passed':sum(c['pass'] for c in checks),'total':len(checks),'max_direct_RG_relative_residual':str(max_scale),'hierarchy_ratio_MQ_over_MT':ratio,'initial_center_A_M_eV':AMnum,'initial_center_conditional_Mlambda_IR_eV':MIR}
OUT.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
print(json.dumps({k:result[k] for k in ['passed','total','max_direct_RG_relative_residual','hierarchy_ratio_MQ_over_MT','initial_center_conditional_Mlambda_IR_eV']},indent=2))
if result['passed'] != result['total']:
    raise SystemExit(1)
