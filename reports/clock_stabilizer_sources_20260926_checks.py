"""Local current matching and convention checks, not a full UV calculation."""
from pathlib import Path
import hashlib
import json

import mpmath as mp
import sympy as sp

rows = []


def identity(name, expr):
    rows.append({"name": name, "passed": sp.simplify(expr) == 0})


def check(name, value):
    rows.append({"name": name, "passed": bool(value)})


J, k, S, gz, gq, M, a = sp.symbols("J k S gz gq M a", nonzero=True, real=True)
current = gz * k + gq * S
boundary = k + S + J**2 / 2 + J * current / M
sol = -current / M
effective = boundary.subs(J, sol)
identity("vector_stationarity", sp.diff(boundary, J).subs(J, sol))
identity("negative_current_square", effective - k - S + current**2 / (2 * M**2))
identity("clock_quartic", sp.expand(effective).coeff(k, 2) + gz**2 / (2 * M**2))
identity("mixed_contact", sp.diff(effective, k, S) + gz * gq / M**2)
identity("matter_quartic", sp.expand(effective).coeff(S, 2) + gq**2 / (2 * M**2))
identity("Lambda_equals_M_at_gsqrt2", (2*M**2/gz**2).subs(gz, sp.sqrt(2)) - M**2)
identity("general_positive_weight", (a * J**2 / 2 + J * current / M).subs(J, -current/(a*M)) + current**2/(2*a*M**2))

z, zb, y, chi, L, fz2 = sp.symbols("z zb y chi L fz2", real=True)
kz = -(z - zb)**2 / 2
keff = kz - kz**2 / L**2
identity("k_y_normalization", kz.subs({z:(chi+sp.I*y)/sp.sqrt(2), zb:(chi-sp.I*y)/sp.sqrt(2)})-y**2)
identity("local_clock_metric", sp.diff(keff,z,zb) - (1-6*kz/L**2))
for n in range(4):
    identity(f"quartic_y_jet_{n}_zero_on_axis", sp.diff(-y**4/L**2,y,n).subs(y,0))
identity("mixed_metric_curvature", sp.diff(sp.log(1-gz*gq*kz/M**2),z,zb).subs(z,zb)+gz*gq/M**2)
identity("mixed_soft_sign", -fz2 * sp.diff(sp.log(1-gz*gq*kz/M**2),z,zb).subs(z,zb)-gz*gq*fz2/M**2)

V, s, sb, xi, xib, v2 = sp.symbols("V s sb xi xib v2", real=True)
j = s+sb+sp.sqrt(2)*M*V
identity("Stueckelberg_gauge_invariance", j.subs({V:V+xi+xib,s:s-sp.sqrt(2)*M*xi,sb:sb-sp.sqrt(2)*M*xib}, simultaneous=True)-j)
identity("canonical_Stueckelberg_metric", sp.diff(j**2/2,s,sb)-1)
# Wess--Zumino V^2|theta^4=-v_mu v^mu/2; canonical Maxwell coefficient -1/4.
identity("physical_vector_mass", sp.expand(j**2/2).coeff(V,2)*(-v2/2)+M**2*v2/2)
old = sp.symbols("S_old", real=True)
identity("Kors_Nath_redefinition", (j**2/2).subs({s:sp.sqrt(2)*old,sb:sp.sqrt(2)*old})-(M*V+2*old)**2)

mp.mp.dps=60
Mp=mp.mpf("2.435e27")
MV=mp.mpf("1e12")
benchmarks=[]
for mg_s in ["1", "1e-5", "1e-6", "1e-18"]:
    mg=mp.mpf(mg_s)
    ratio=mp.sqrt(3)*Mp*mg/MV**2
    benchmarks.append({"mG_eV":mg_s,"F_T_can_over_MV2":str(ratio),"under_declared_0p1":bool(ratio<=mp.mpf("0.1"))})
check("one_eV_fails_sufficient_hierarchy", not benchmarks[0]["under_declared_0p1"])
check("new_examples_pass_sufficient_hierarchy", all(b["under_declared_0p1"] for b in benchmarks[1:]))
check("charged_over_mediator_is_0p1", mp.mpf("1e11")/MV == mp.mpf("0.1"))
check("mediator_over_KK_is_0p01", MV/mp.mpf("1e14") == mp.mpf("0.01"))
check("local_metric_positive_small_patch", 1-6*mp.mpf("0.1")**2>0)
check("quartic_truncation_not_global_metric", 1-6*mp.mpf("0.5")**2<0)

out={
 "scope":"Independent algebra and reference conventions; not full heavy-multiplet spectrum, loop calculation or UV completion",
 "passed":sum(r["passed"] for r in rows),"total":len(rows),"checks":rows,
 "benchmarks":benchmarks,
 "mG_upper_eV_for_F_ratio_0p1":str(mp.mpf("0.1")*MV**2/(mp.sqrt(3)*Mp)),
 "versions":{"sympy":sp.__version__,"mpmath":mp.__version__},
 "script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
}
dest=Path(__file__).with_suffix(".json")
dest.write_text(json.dumps(out,ensure_ascii=False,indent=2)+"\n")
print(json.dumps({"passed":out["passed"],"total":out["total"],"destination":str(dest)}))
if out["passed"]!=out["total"]:
    raise SystemExit(1)
