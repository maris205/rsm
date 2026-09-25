#!/usr/bin/env python3
"""Independent local operator algebra and conditional CW transfer checks.

This does not import the production EFT matching code and does not compute a
complete supergravity loop action. No observational data or statistical tests.
"""
from pathlib import Path
from datetime import datetime, timezone
import csv
import hashlib
import json
import platform
import numpy as np
import sympy as sp
import mpmath as mp

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
checks = []


def check(name, expression):
    residual = sp.simplify(expression)
    checks.append({"name": name, "passed": residual == 0,
                   "residual": str(residual)})


P, V, m2, curvature = sp.symbols("P V m2 curvature", positive=True)
check("metric_deformation_common_soft_mass",
      m2 + V/P - (V+3*m2*P)/(3*P) - curvature
      - (2*V/(3*P)-curvature))
Ki, Fi, logMi, elli = sp.symbols("Ki Fi logMi elli")
check("bilinear_deformation", -Fi*(Ki/P+logMi-2*Ki/(3*P)-elli)
      - (-Fi*Ki/(3*P)-Fi*logMi+Fi*elli))
ai = sp.symbols("ai")
check("holomorphic_Q_rescaling_b_invariance",
      -Fi*(logMi-ai)+Fi*(elli-ai)-(-Fi*logMi+Fi*elli))
tau, taub, c = sp.symbols("tau taub c", real=True)
ell = c*tau*taub
check("contact_value", ell.subs({tau: 0, taub: 0}))
check("contact_first_derivative", sp.diff(ell,tau).subs({tau: 0, taub: 0}))
check("contact_irreducible_curvature", sp.diff(ell,tau,taub)-c)

# Complete auxiliary equations of the undeformed hidden no-scale metric.
O = sp.symbols("Omega", positive=True)
kz, kzb, W, Wb, wz, wzb = sp.symbols("kz kzb W Wb wz wzb")
metric = sp.Matrix([[3*P/O**2, -kzb/O**2],
                    [-kz/O**2, 1/O+kz*kzb/(3*P*O**2)]])
Dbar = sp.Matrix([-3*Wb/O, wzb+kzb*Wb/(P*O)])
aux = -O**sp.Rational(-3,2)*metric.T.inv()*Dbar
check("complete_auxiliary_FT", aux[0]-(3*Wb-kz*wzb)/(3*P*sp.sqrt(O)))
check("complete_auxiliary_FZ", aux[1]+wzb/sp.sqrt(O))

# Independent direct scalar-potential expansion at the clock-off center.
Q, Qb, mg = sp.symbols("Q Qb mg", real=True)
metric_Q = sp.Matrix([[3+(2+c)*Q*Qb, -Q], [-Qb, 1]])
Ds = sp.Matrix([(-3-Q*Qb)*mg, Qb*mg])
Dbs = sp.Matrix([(-3-Q*Qb)*mg, Q*mg])
Vq = sp.exp(Q*Qb)*((Dbs.T*metric_Q.inv()*Ds)[0]-3*mg**2)
check("direct_full_potential_contact_mass",
      sp.diff(Vq,Q,Qb).subs({Q: 0,Qb: 0})+c*mg**2)

y, wr, co, si = sp.symbols("y wr co si", real=True)
S = 3*mg*P*(co+sp.I*si)+wr*(3+2*sp.I*y)*sp.exp(-sp.I*y)
Sb = 3*mg*P*(co-sp.I*si)+wr*(3-2*sp.I*y)*sp.exp(sp.I*y)
X = S*Sb/(9*P**2)
check("contact_X_chi_derivative", (-wr*sp.diff(X,wr)).subs(y,0)
      +2*wr*(mg*P*co+wr)/P**2)
check("contact_X_y_one_third", sp.diff(X,y).subs(y,0)+2*mg*wr*si/(3*P))

# Kähler-bilinear examples, in a local gauge W>0, F^T=mg.
# Jets are H, d_T H, d_Tbar H, d_T d_Tbar H; d_T ln(Z+Z-)=-2.
h = sp.symbols("h")
def h_terms(h0, ht, hb, htb):
    mu = mg*h0+mg*hb
    b = 2*mg**2*h0+mg**2*hb-mg**2*(ht+2*h0)-mg**2*(htb+2*hb)
    return mu, b

for name, jets, expected in (
    ("holomorphic_constant", (h,0,0,0), (h*mg,0)),
    ("antiholomorphic_linear", (0,0,h,0), (h*mg,-h*mg**2)),
    ("inverse_t_no_scale_bilinear", (h,-h,-h,2*h), (0,0)),
):
    got = h_terms(*jets)
    check(name+"_fermion_mass",got[0]-expected[0])
    check(name+"_scalar_bilinear",got[1]-expected[1])

s, mu, d, h2 = sp.symbols("s mu d h2", positive=True)
x = sp.symbols("x", positive=True)
F = x**2*(sp.log(x/mu**2)-sp.Rational(3,2))
check("CW_first_derivative",sp.diff(F,x)-2*x*(sp.log(x/mu**2)-1))
check("CW_second_derivative",sp.diff(F,x,2)-2*sp.log(x/mu**2))
X0, Xa, sa, A, nu, b2 = sp.symbols("X Xa sa A nu b2")
L = sp.log(s/mu**2)
Vc = -c*s*X0*(L-1)/(8*sp.pi**2)
grad = sp.diff(Vc,s)*sa+sp.diff(Vc,X0)*Xa
check("contact_force_kernel",grad+c*(X0*L*sa+s*(L-1)*Xa)/(8*sp.pi**2))
check("contact_scale_derivative",mu*sp.diff(Vc,mu)-c*s*X0/(4*sp.pi**2))
Vsel = s*X0*((-c+b2/2)*L+c+nu)/(8*sp.pi**2)
check("three_coefficients_two_shapes",sp.diff(Vsel,c)+2*sp.diff(Vsel,b2)-sp.diff(Vsel,nu))
check("selected_RG_cancellation",mu*sp.diff(Vsel,mu)+sp.diff(Vsel,nu)*2*(-c+b2/2))

# Numeric material passport: historical chi values are coordinates, not a new
# background solution of this deformed action.
parpath = REPO/"experiments/sugra_shift_v01/results/inputs.json"
trajpath = REPO/"experiments/sugra_shift_v01/results/trajectories.csv"
par = json.loads(parpath.read_text())
with trajpath.open() as stream:
    rows = [row for row in csv.DictReader(stream) if row["case"] == "shift_axis"]
chi = np.array([float(row["chi"]) for row in rows])
PP = par["Mpl_eV"]**2
UU = par["potential_unit_eV4"]*par["W_i"]*np.exp(-2*(chi-par["chi_i"]))
ww = -np.sqrt(par["F_squared_eV2"]*UU/2)
mi = par["holomorphic_mass_i_eV"]
den = 1+par["xi"]/par["chi_i"]**2
ss = mi**2*(1+par["xi"]/chi**2)/den
sx = -2*mi**2*par["xi"]/chi**3/den
LL = np.log(ss/mi**2)
tolerances = []
for mass in (1e-6,1.,1e5):
    for name, cosine, sine in (("zero",1.,0.),("quadrature",0.,1.),("pi",-1.,0.)):
        XX = mass**2+2*mass*cosine*ww/PP+(ww/PP)**2
        XXx = -2*mass*cosine*ww/PP-2*(ww/PP)**2
        XXy = -2*mass*sine*ww/(3*PP)
        kx = (ss*(LL-1)*XXx+XX*LL*sx)/(8*np.pi**2)
        ky = ss*(LL-1)*XXy/(8*np.pi**2)
        rat = np.hypot(kx,ky)/(2*UU)
        idx = int(np.argmax(rat))
        tolerances.append({"mG_eV":mass,"phase":name,"max_R_per_abs_c":float(rat[idx]),
                           "conditional_c_10percent":float(.1/rat[idx]),
                           "max_chi":float(chi[idx]),"initial_R_per_abs_c":float(rat[0])})

# Direct 160-digit CW differentiation, with c the inserted soft coefficient.
# The first 100-digit run lost ~38 digits at L=0 and failed the unchanged
# 1e-65 threshold. Preserve that attempt in the machine-readable history.
mp.mp.dps = 160
numeric_checks = []
for index in (0,512,1024):
    sm = mp.mpf(str(ss[index])); xm = mp.mpf("1")
    um = mp.mpf(str(UU[index])); smx = mp.mpf(str(sx[index]))
    xm_x = -2*mp.mpf(str(ww[index]))/mp.mpf(str(PP))
    mum = mp.mpf(str(mi))
    def ff(xx):
        return xx**2*(mp.log(xx/mum**2)-mp.mpf("1.5"))
    def vv(cc, shift):
        sv = sm+smx*shift; xv = xm+xm_x*shift
        return (2*ff(sv-cc*xv)-2*ff(sv))/(32*mp.pi**2)
    raw = mp.diff(lambda sh: mp.diff(lambda cc: vv(cc,sh),mp.mpf(0)),mp.mpf(0))
    ln = mp.log(sm/mum**2)
    analytic = -(xm*ln*smx+sm*(ln-1)*xm_x)/(8*mp.pi**2)
    rel = abs(raw-analytic)/max(abs(analytic),mp.mpf("1e-90"))
    numeric_checks.append({"name":f"direct_CW_mixed_derivative_{index}","passed":bool(rel<mp.mpf("1e-65")),
                           "relative_error":str(rel)})

result = {
    "date_utc":datetime.now(timezone.utc).isoformat(),
    "scope":"Independent local operator algebra; conditional selected-CW tolerance, not complete SUGRA or observation",
    "versions":{"python":platform.python_version(),"sympy":sp.__version__,"numpy":np.__version__,"mpmath":mp.__version__},
    "script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    "input_sha256":{str(p.relative_to(REPO)):hashlib.sha256(p.read_bytes()).hexdigest() for p in (parpath,trajpath)},
    "symbolic_checks":checks,"numeric_checks":numeric_checks,
    "checks_total":len(checks)+len(numeric_checks),
    "checks_passed":sum(z["passed"] for z in checks+numeric_checks),
    "conditional_boundary":"mu_Q=100GeV, local potential nu(mu_Q)=0, no other operators, leading d/s and |B|/s; c varied singly",
    "direct_CW_decimal_precision":mp.mp.dps,
    "old_coordinate_count":len(chi),"final_mass_squared_ratio":float(ss[-1]/ss[0]),
    "final_log_mass_squared_ratio":float(LL[-1]),"conditional_tolerances":tolerances,
}
out = HERE/"matching_operator_audit_20260925_checks.json"
if out.exists():
    previous = json.loads(out.read_text())
    history = previous.get("execution_history", [])
    if previous["checks_passed"] != previous["checks_total"]:
        history.append({"date_utc":previous["date_utc"],"checks_passed":previous["checks_passed"],
                        "checks_total":previous["checks_total"],"script_sha256":previous["script_sha256"],
                        "decimal_precision":previous.get("direct_CW_decimal_precision",100),
                        "failed_checks":[z for z in previous["symbolic_checks"]+previous["numeric_checks"] if not z["passed"]],
                        "resolution":"Increase direct CW precision to 160 digits; retain 1e-65 criterion."})
    result["execution_history"] = history
out.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n")
print(json.dumps({"checks_passed":result["checks_passed"],"checks_total":result["checks_total"],
                  "c10_mG1eV":tolerances[3]["conditional_c_10percent"]},indent=2))
raise SystemExit(0 if result["checks_passed"] == result["checks_total"] else 1)
