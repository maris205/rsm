#!/usr/bin/env python3
"""All-cell comparison, independent formulas; never imports the main experiment."""
from pathlib import Path
import csv
import hashlib
import json
import math
import platform

import mpmath as mp
import numpy as np
import scipy
from scipy.integrate import quad
from scipy.optimize import brentq

mp.mp.dps = 200
HERE = Path(__file__).resolve()
REPO = HERE.parent.parent
EXP = REPO / "experiments/light_modulus_v01"
RESULTS = EXP / "results"
source_path = REPO / "experiments/threshold_rg_v01/results/inputs.json"
source = json.loads(source_path.read_text())
summary = json.loads((RESULTS / "summary.json").read_text())
groups = []
comparisons = {}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def check(name, ok, detail=None):
    groups.append({"name": name, "passed": bool(ok), "detail": detail})


def normed_residual(actual, expected, transient=False):
    a, e = mp.mpf(str(actual)), mp.mpf(expected)
    if transient:
        allowance = mp.mpf("5e-12") * abs(e)
    else:
        allowance = mp.mpf("2e-54") * max(abs(e), mp.mpf("1e-200"))
    if e == 0:
        allowance = mp.mpf("1e-13") if transient else mp.mpf("1e-90")
    return abs(a-e)/allowance


def compare_rows(filename, expected, transient=False):
    with (RESULTS / filename).open() as f:
        rows = list(csv.DictReader(f))
    failed = []
    max_residual = mp.mpf(0)
    max_location = None
    cells, numeric, textcells = 0, 0, 0
    if len(rows) != len(expected):
        failed.append({"reason": "row count", "actual": len(rows), "expected": len(expected)})
    for i, (row, exp) in enumerate(zip(rows, expected)):
        if set(row) != set(exp):
            failed.append({"row": i, "reason": "column schema", "actual": sorted(row), "expected": sorted(exp)})
        for key, value in exp.items():
            cells += 1
            if isinstance(value, str):
                textcells += 1
                ok = row.get(key) == value
                resid = mp.mpf(0) if ok else mp.inf
            else:
                numeric += 1
                resid = normed_residual(row.get(key), value, transient)
                ok = resid <= 1
            if resid > max_residual:
                max_residual, max_location = resid, [i, key]
            if not ok:
                failed.append({"row": i, "column": key, "actual": row.get(key),
                               "expected": str(value), "residual_over_tolerance": str(resid)})
    item = {"rows": len(rows), "cells_compared": cells, "numeric_cells": numeric,
            "text_cells": textcells, "failures": failed,
            "max_residual_over_tolerance": mp.nstr(max_residual, 12), "worst_cell": max_location}
    comparisons[filename] = item
    check("all cells: " + filename, not failed, item)
    return rows


def geometric_mode(q, n):
    q = mp.mpf(q)
    norm = mp.sqrt((q**(2*n+2)-1)/(q*q-1))
    return [1/norm] + [-q**a/norm for a in range(1, n+1)]


q, n = 3, 127
L = mp.mpf("2e12")
g = mp.mpf(1)
v = L/(2*g*q)
mg = mp.mpf("1e-15")
Pgrav = mp.mpf(str(source["parameters"]["P"]))
hbar = mp.mpf(str(source["parameters"]["hbar_eV_s"]))
S = 3*Pgrav*mg**2
u = geometric_mode(q, n)
p00 = 1-u[0]**2
ci, cx = 1/(v*L*mp.sqrt(p00)), q/(v*L*mp.sqrt(p00))
bi, bx = 2*v*ci*u[0], 2*v*cx*u[-1]
tau = (q-1/mp.mpf(q))/(mp.mpf(q)**n-mp.mpf(q)**(-n))
b = -2*v*v*cx*u[-1]

# H = d I - q adjacency, d = 2q cosh(a). For the endpoint Green function,
# H^-1[j,n-1] = sinh((j+1)a)/(q sinh((n+1)a)).
# -d/dd of this Green function is H^-2[j,n-1]. No matrix solve is reused.
a = mp.log(q)
hinv2_last = []
for j in range(n):
    green = mp.sinh((j+1)*a)/(q*mp.sinh((n+1)*a))
    derivfactor = ((n+1)*mp.coth((n+1)*a)-(j+1)*mp.coth((j+1)*a))/(2*q*mp.sinh(a))
    hinv2_last.append(green*derivfactor)
avec = [q*hinv2_last[0]] + [hinv2_last[j-1]-q*hinv2_last[j] for j in range(1, n)] + [hinv2_last[-1]]
h = avec[-1]
kappa = S*cx*cx*h/(g*g)
check("closed Green derivative is orthogonal to null mode", abs(mp.fsum(u[j]*avec[j] for j in range(n+1))) < mp.mpf("1e-185"))


def cubic_shift(B, kap):
    if kap == 0:
        return mp.mpf(0)
    # Cardano for D³-BD²-kap=0, then t=D-B; unlike main's Newton root.
    mid = B**3/27+kap/2
    disc = mp.sqrt(kap*B**3/27+kap**2/4)
    D = B/3+mp.root(mid+disc, 3)+mp.root(mid-disc, 3)
    return D-B


endpoint = []
for qq in (2, 3, 5):
    for nn in (2, 4, 8, 16, 32, 64, 127, 128):
        uu = geometric_mode(qq, nn)
        pp = 1-uu[0]**2
        endpoint.append(dict(q=qq, n=nn, u_left=uu[0], u_right=uu[-1],
            tau=(qq-1/mp.mpf(qq))/(mp.mpf(qq)**nn-mp.mpf(qq)**(-nn)),
            right_light_weight=uu[-1]**2, right_heavy_weight=1-uu[-1]**2,
            betaI_times_Lambda=2*uu[0]/mp.sqrt(pp), betaX_times_Lambda=2*qq*uu[-1]/mp.sqrt(pp)))
compare_rows("endpoint_scaling.csv", endpoint)
compare_rows("mode_profile.csv", [dict(link=j, u=x, light_weight=x*x, heavy_diagonal=1-x*x) for j,x in enumerate(u)])

physical = []
for j in range(101):
    r = mp.mpf(j)/100
    baseline = 1+b*r
    t = cubic_shift(baseline, kappa)
    D = baseline+t
    # Dimensionless currents da=d_a/(2v²), differentiated by the implicit cubic.
    shift = S*cx/(2*g*g*v*v*D*D)
    dt_dr = -2*b*t/(D+2*t)
    dshift_dr = -2*shift*(b+dt_dr)/D
    da = [-r*u[j]+shift*avec[j] for j in range(n+1)]
    dar = [-u[j]+dshift_dr*avec[j] for j in range(n+1)]
    Bzero = mp.fsum(x*x/mp.sqrt(1+r*r*x*x) for x in u)
    Bfinite = mp.fsum(dar[j]**2/mp.sqrt(1+da[j]**2) for j in range(n+1))
    physical.append(dict(r=r, Dflat_metric_X=baseline, finiteF_metric_X=D,
        heavy_metric_shift=t, V_over_S=(baseline+mp.mpf("1.5")*t)/D**2,
        dV_dr_over_S=-b/D**2, d2V_dr2_over_S=2*b*(b+dt_dr)/D**3,
        radial_kinetic_B=Bzero, finiteF_radial_kinetic_B=Bfinite,
        kahler_metric_C=1/Bzero, clock_metric=1+2*v*v*ci*da[0],
        force_relative_Dflat_error=abs((baseline/D)**2-1),
        kinetic_relative_Dflat_error=abs(Bfinite/Bzero-1)))
compare_rows("physical_valley.csv", physical)

stress = []
for label in ("0", "1e-6", "1e-3", "0.1", "1"):
    kap = mp.mpf(label)
    for j in range(61):
        r = mp.mpf(j)/20
        baseline = 1+b*r
        t = cubic_shift(baseline, kap)
        D = baseline+t
        stress.append(dict(kappa=kap, r=r, hidden_metric=D, t=t,
            V_hat=kap/D+t*t/2, dV_hat_dr=-kap*b/D**2,
            V_over_S=(baseline+mp.mpf("1.5")*t)/D**2 if kap else "undefined_zero_source"))
compare_rows("source_stress.csv", stress)

# All transient rows reconstructed by energy quadrature and scalar inversion.
# This never integrates the equations of motion or invokes solve_ivp.
uf = np.array([float(x) for x in u])
weights = uf*uf
bf = float(b)
target = mp.mpf("0.1")/b
def Bfloat(r):
    return float(np.sum(weights/np.sqrt(1+weights*r*r)))

def energy_time(r):
    if r == 0:
        return 0.0
    return quad(lambda z: math.sqrt(2*Bfloat(z*z)*(1+bf*z*z)/bf),
                0, math.sqrt(r), epsabs=2e-14, epsrel=2e-14)[0]

event_s = energy_time(float(target))
timeunit = hbar*v/mp.sqrt(S)
transient = []
for j in range(201):
    s = j*event_s/200
    if j == 0:
        r = rd = 0.0
    elif j == 200:
        r = float(target)
        rd = math.sqrt(2*bf*r/((1+bf*r)*Bfloat(r)))
    else:
        r = brentq(lambda rr: energy_time(rr)-s, 0, float(target)*1.001, xtol=1e-20, rtol=2e-14)
        rd = math.sqrt(2*bf*r/((1+bf*r)*Bfloat(r)))
    transient.append(dict(s=s, local_elapsed_seconds=mp.mpf(str(s))*timeunit,
        r=r, dr_ds=rd, hidden_metric_change=bf*r, energy_over_S=1, kinetic_B=Bfloat(r)))
actual_transient = compare_rows("local_transient.csv", transient, transient=True)

benchmark_expected = dict(q=q, n=n, g=g, Lambda_eV=L, v_eV=v, mg_eV=mg, S_eV4=S,
    cI_eVminus2=ci, cX_eVminus2=cx, u_left=u[0], u_right=u[-1], betaI_eVminus1=bi,
    betaX_eVminus1=bx, tau=tau, h=h, kappa=kappa, b=b,
    leading_sigma_force_at_origin_eV3=-S*bx, leading_sigma_Hessian_eV2=2*S*bx**2,
    leading_clock_partial_Hessian_eV2=4*tau*S/L**2,
    leading_clock_covariant_Hessian_eV2=2*tau*S/L**2,
    leading_X_partial_Hessian_eV2=S*(4/L**2+bx**2/2),
    transient_time_unit_seconds=timeunit, transient_event_s=event_s,
    transient_energy_quadrature_s=event_s, transient_event_seconds=mp.mpf(str(event_s))*timeunit,
    transient_event_r=target,
    max_transient_energy_error=max(abs(float(row["energy_over_S"])-1) for row in actual_transient),
    max_finiteF_force_relative_error=max(row["force_relative_Dflat_error"] for row in physical),
    max_finiteF_kinetic_relative_error=max(row["kinetic_relative_Dflat_error"] for row in physical))
benchmark_failures = []
for key, value in benchmark_expected.items():
    ratio = normed_residual(summary["benchmark"][key], value, "transient" in key)
    if ratio > 1:
        benchmark_failures.append(dict(key=key, actual=summary["benchmark"][key], expected=str(value), residual_over_tolerance=str(ratio)))
check("all numeric summary benchmark values", not benchmark_failures,
      {"cells_compared": len(benchmark_expected), "failures": benchmark_failures})
check("summary benchmark numeric schema complete", set(summary["benchmark"])-{"hessian_scope"} == set(benchmark_expected))
check("summary row counts", summary["row_counts"] == {Path(k).stem:v["rows"] for k,v in comparisons.items()})
check("main check record aggregation", summary["check_count"] == len(summary["checks"]) and summary["checks_passed"] == sum(c["passed"] for c in summary["checks"]) == 549)
check("summary source hashes current", all(sha(REPO/k) == v for k,v in summary["source_hashes"].items()))
check("scope statements retain nonstationarity and cosmology exclusion",
      "nonstationary" in summary["benchmark"]["hessian_scope"] and "No finite interior stationary" in summary["retained_physical_failure"] and "cosmic lifetime prediction" in summary["exclusions"])

paths = [source_path, EXP/"protocol.md", EXP/"code/run_light_modulus.py", RESULTS/"summary.json"] + sorted(RESULTS.glob("*.csv"))
out = {
    "scope": "Internal independent implementation after reading protocol/schema/main code; same-model-family review, not blind or external peer review.",
    "independent_methods": ["closed geometric null mode", "hyperbolic Green-function derivative for H^-2 and (QQ^T)^+", "Cardano cubic branch", "dimensionless Higgs-current kinetic reconstruction", "energy quadrature plus scalar inversion for every transient row, no ODE integration"],
    "versions": {"python": platform.python_version(), "mpmath": mp.__version__, "numpy": np.__version__, "scipy": scipy.__version__},
    "precision_decimal_digits": mp.mp.dps,
    "tolerances": {"analytic_numeric": "2e-54 * max(abs(expected),1e-200); exact-zero absolute 1e-90", "transient": "5e-12 * abs(expected); exact-zero absolute 1e-13", "undefined": "exact string comparison"},
    "input_hashes": {str(p.relative_to(REPO)): sha(p) for p in paths}, "script_sha256": sha(HERE),
    "csv_rows_compared": sum(v["rows"] for v in comparisons.values()),
    "csv_cells_compared": sum(v["cells_compared"] for v in comparisons.values()),
    "csv_numeric_cells": sum(v["numeric_cells"] for v in comparisons.values()),
    "csv_text_cells": sum(v["text_cells"] for v in comparisons.values()),
    "summary_numeric_cells_compared": len(benchmark_expected), "passed": sum(c["passed"] for c in groups),
    "total": len(groups), "checks": groups,
    "event_seconds_independent": str(mp.mpf(str(event_s))*timeunit),
}
destination = HERE.with_suffix(".json")
if destination.exists():
    previous = json.loads(destination.read_text())
    if previous.get("passed") != previous.get("total"):
        out["previous_failed_runs"] = previous.get("previous_failed_runs", []) + [{k:v for k,v in previous.items() if k != "previous_failed_runs"}]
destination.write_text(json.dumps(out, ensure_ascii=False, indent=2)+"\n")
print(json.dumps({k:out[k] for k in ("passed", "total", "csv_rows_compared", "csv_cells_compared", "csv_numeric_cells", "csv_text_cells", "summary_numeric_cells_compared")}, indent=2))
for item in groups:
    if not item["passed"]:
        print(json.dumps(item, indent=2, ensure_ascii=False))
if out["passed"] != out["total"]:
    raise SystemExit(1)
