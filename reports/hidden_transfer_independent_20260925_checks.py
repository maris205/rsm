#!/usr/bin/env python3
"""Independent algebra checks; does not import the new experiment implementation."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import numpy as np
import sympy as s

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
checks = []


def check(name, expression):
    reduced = s.trigsimp(s.expand(expression))
    if reduced != 0:
        reduced = s.simplify(reduced)
    checks.append({"name": name, "passed": reduced == 0,
                   "residual": str(reduced)})


# Independent real/imaginary construction from D_Z W and W, in F_chi^2 H^2 units.
ep, y, beta, lam, amp, x, th = s.symbols(
    "epsilon y beta lambda amp x theta", real=True)
q = amp*s.exp(-x)
dr = q*(s.cos(y)+ep*y*s.sin(y)) + s.sqrt(2)*beta*y*s.sin(th)
di = q*(ep*y*s.cos(y)-s.sin(y)) - s.sqrt(2)*beta*y*s.cos(th)
wr = -q*s.cos(y)/s.sqrt(2)+beta*s.cos(th)/ep
wi = q*s.sin(y)/s.sqrt(2)+beta*s.sin(th)/ep
direct_bracket = dr**2+di**2+3*beta**2/ep+lam-3*ep*(wr**2+wi**2)
cross = (3-2*ep*y*y)*s.cos(y+th)+2*y*s.sin(y+th)
compact_bracket = lam + q*q*(1-3*ep/2+ep**2*y*y) + 2*beta**2*y*y + s.sqrt(2)*beta*q*cross
check("complete_two_field_potential", direct_bracket-compact_bracket)
V = s.exp(ep*y*y)*compact_bracket
c = 1-3*ep/2
C = s.sqrt(2)*beta*q
axis_expectations = {
    "V": lam+c*q*q+3*C*s.cos(th),
    "Vx": -2*c*q*q-3*C*s.cos(th),
    "Vy": -C*s.sin(th),
    "Vxx": 4*c*q*q+3*C*s.cos(th),
    "Vxy": C*s.sin(th),
    "Vyy": 2*ep*lam+2*ep*(1-ep/2)*q*q+4*beta*beta+(1+2*ep)*C*s.cos(th),
}
expressions = {"V": V, "Vx": s.diff(V,x), "Vy": s.diff(V,y),
               "Vxx": s.diff(V,x,2), "Vxy": s.diff(V,x,y), "Vyy": s.diff(V,y,2)}
for key, expression in expressions.items():
    check("axis_"+key, expression.subs(y,0)-axis_expectations[key])
for phase, label in [(0,"aligned"),(s.pi/2,"quadrature"),(s.pi,"anti_aligned")]:
    for key in ("Vx","Vy"):
        check(label+"_"+key, expressions[key].subs({y:0,th:phase})-axis_expectations[key].subs(th,phase))

# Large beta, theta=pi/2, y=v/beta: retain the heavy direction before minimizing.
v = s.symbols("v", real=True)
heavy = s.limit(V.subs({th:s.pi/2,y:v/beta}),beta,s.oo)
heavy_expected = lam+c*q*q+2*v*v-s.sqrt(2)*q*v
check("heavy_quadrature_rescaled_potential", heavy-heavy_expected)
vstar = q/(2*s.sqrt(2))
check("heavy_quadrature_stationary",s.diff(heavy,v).subs(v,vstar))
check("heavy_quadrature_effective_potential",heavy.subs(v,vstar)-(lam+(c-s.Rational(1,4))*q*q))

# Charged scalar Hessian from the full potential, at y=0, independently expanded.
# W=w_r+i w_i, D_Z W=D, M and M_Z are real on this axis.
p, D, M, Mz, f, ar, ai, a, b = s.symbols("Mp D M Mz f wr wi a b", real=True)
W = ar+s.I*ai

def charged_direct(phase):
    qp = a
    qm = b*phase
    Wq = W+M*qp*qm
    Dz = D+Mz*qp*qm
    Dp = M*qm+s.conjugate(qp)*Wq/p**2
    Dm = M*qp+s.conjugate(qm)*Wq/p**2
    # exp(K/Mp^2) expanded to quadratic; sufficient for mass Hessian.
    prefactor = 1+(a*a+b*b)/p**2
    return s.expand(prefactor*(Dz*s.conjugate(Dz)+f*f+Dp*s.conjugate(Dp)+Dm*s.conjugate(Dm)-3*Wq*s.conjugate(Wq)/p**2))

d_general = (D*D+f*f)/p**2-2*(ar*ar+ai*ai)/p**4
B_general = D*Mz-M*s.conjugate(W)/p**2
direct_real = charged_direct(1)
direct_imag = charged_direct(s.I)
for variable in (a,b):
    diagonal = s.diff(direct_real,variable,2).subs({a:0,b:0})/2
    check("charged_diagonal_"+str(variable),diagonal-(M*M+d_general))
check("charged_B_real",s.diff(direct_real,a,b).subs({a:0,b:0})-2*s.re(B_general))
check("charged_B_imag",s.diff(direct_imag,a,b).subs({a:0,b:0})+2*s.im(B_general))
F, rootU, mg, rho, r = s.symbols("F sqrtU mg rho r", real=True)
sub = {ar:-F*rootU/s.sqrt(2)+mg*p*p*s.cos(th),
       ai:mg*p*p*s.sin(th),D:rootU,f*f:3*mg*mg*p*p+rho,Mz:s.sqrt(2)*M*r/F}
d_target = mg*mg+rho/p**2+(1-F*F/p**2)*rootU**2/p**2+2*s.sqrt(2)*F*mg*rootU*s.cos(th)/p**2
B_target = s.sqrt(2)*rootU*M/F*(r+F*F/(2*p*p))-mg*M*(s.cos(th)-s.I*s.sin(th))
check("charged_hidden_diagonal",d_general.subs(sub)-d_target)
check("charged_hidden_holomorphic_B",B_general.subs(sub)-B_target)

# Heavy-hidden but m_G << M charged threshold, with mu fixed.
ss, ds, hs, mu, mm, sp = s.symbols("s d h2 mu m sp", positive=True)
L = s.log(ss/mu**2)
small_split_cw = (4*ss*ds*(L-1)+2*hs*L)/(32*s.pi**2)
heavy_cw = small_split_cw.subs({ds:mm**2,hs:mm**2*ss})
check("CW_hidden_leading_potential",heavy_cw-mm**2*ss*(6*L-4)/(32*s.pi**2))
check("CW_hidden_leading_derivative",s.diff(heavy_cw,ss)*sp-mm**2*sp*(3*L+1)/(16*s.pi**2))

inputs_path=REPO/'experiments/sugra_shift_v01/results/inputs.json'
trajectory_path=REPO/'experiments/sugra_shift_v01/results/trajectories.csv'
inputs=json.loads(inputs_path.read_text())
rows=np.genfromtxt(trajectory_path,names=True,delimiter=',',dtype=None,encoding='utf8')
rows=rows[rows['case']=='shift_axis']
chi=rows['chi']; eps=inputs['epsilon']; H=inputs['H_ref_eV']; Mp=inputs['Mpl_eV']
u=inputs['W_i']*np.exp(-2*(chi-inputs['chi_i']))
coefficient=1-1.5*eps
rhoL=3*inputs['omega_lambda']*H*H*Mp*Mp
bounds={}
for phase,label in [(0.,'aligned'),(np.pi/2,'quadrature'),(np.pi,'anti_aligned')]:
    bound=.1*np.sqrt(2)*coefficient*np.sqrt(u)*H/np.sqrt(9*np.cos(phase)**2+np.sin(phase)**2)
    bounds[label+'_frozen_axis_force_0p1_mG_eV']=float(np.min(bound))
bounds.update(rho_lambda_eV4=rhoL,nilpotent_sqrtf_floor_eV=rhoL**.25,
              sqrtf_at_least_100GeV_order_one_mG_min_eV=np.sqrt(1e44-rhoL)/(np.sqrt(3)*Mp),
              U_dimless_initial=float(u[0]),U_dimless_final=float(u[-1]))
mass2=1e22*.5*(1+(inputs['chi_i']/chi)**2)
mass2_chi=-1e22*inputs['chi_i']**2/chi**3
logmass=np.log(mass2/1e22)
loop_coefficient=np.abs(mass2_chi*(3*logmass+1))/(32*np.pi**2*coefficient*u*inputs['potential_unit_eV4'])
max_loop_coefficient=float(np.max(loop_coefficient))
bounds.update(leading_CW_force_ratio_coefficient_per_eV2=max_loop_coefficient,
              leading_CW_force_0p1_formal_mG_eV=float(np.sqrt(.1/max_loop_coefficient)),
              leading_CW_force_ratio_at_nominal_100GeV_cutoff=max_loop_coefficient*bounds['sqrtf_at_least_100GeV_order_one_mG_min_eV']**2)
out={"created_utc":datetime.now(timezone.utc).isoformat(),"scope":"independent symbolic and frozen-parent-background budget; no cosmological trajectory integration",
     "script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
     "inputs_sha256":hashlib.sha256(inputs_path.read_bytes()).hexdigest(),
     "parent_trajectory_sha256":hashlib.sha256(trajectory_path.read_bytes()).hexdigest(),
     "checks":checks,"passed":sum(c['passed'] for c in checks),"total":len(checks),
     "frozen_axis_budgets":bounds}
output=HERE/'hidden_transfer_independent_20260925_checks.json'
output.write_text(json.dumps(out,indent=2,ensure_ascii=False)+'\n')
print(f"{out['passed']}/{out['total']} symbolic checks passed")
print(json.dumps(bounds,indent=2))
if out['passed']!=out['total']:
    raise SystemExit(1)
