#!/usr/bin/env python3
"""Independent table recomputation. Does not read/import run_chain.py."""
from pathlib import Path
from functools import lru_cache
import csv
import hashlib
import json
import math

import mpmath as mp
import numpy as np
from scipy.optimize import brentq

BASE = Path(__file__).resolve().parents[1]
ROOT = BASE.parents[1]
mp.mp.dps = 90
inputs_path = ROOT/"experiments/persistent_restoration_v01/results/inputs.json"
trajectory_path = ROOT/"experiments/sugra_shift_v01/results/trajectories.csv"
parent = json.loads(inputs_path.read_text())
parameters = parent["parent_parameters"]
P = parameters["P"]
F = math.sqrt(parameters["F_squared_eV2"])
A = float(mp.zeta(3)/(48*mp.pi**4)*mp.mpf("1e14")**2/mp.mpf(str(P)))
rho, Ui = parameters["rho"], parameters["Ui"]
with trajectory_path.open() as source:
    old_rows = [r for r in csv.DictReader(source) if r["case"] == "shift_axis"]
chi = np.array([float(r["chi"]) for r in old_rows])
U = Ui*np.exp(-2*(chi-chi[0]))
MQ, KK = 1e11, 1e14
checks = []
aggregate = {}
failures = []
execution_notes = []


def compare(label, field, expected, obtained, relative_tolerance=1e-10):
    """No absolute floor that could silently accept an incorrect tiny coupling."""
    if isinstance(expected, (bool, np.bool_)):
        error = 0.0 if bool(expected) == (str(obtained).lower() == "true") else 1.0
        passed = error == 0
    elif isinstance(expected, str):
        error = 0.0 if str(obtained) == expected else 1.0
        passed = error == 0
    else:
        expected, actual = float(expected), float(obtained)
        error = abs(actual-expected)/abs(expected) if expected != 0 else abs(actual)
        passed = math.isfinite(error) and (error == 0 if expected == 0 else error <= relative_tolerance)
    key = f"{label.split(':')[0]}:{field}"
    item = aggregate.setdefault(key, {"count":0,"passed":0,"max_relative_error":0.0})
    item["count"] += 1
    item["passed"] += int(passed)
    item["max_relative_error"] = max(item["max_relative_error"], error)
    checks.append(bool(passed))
    if not passed:
        failures.append({"row":label,"field":field,"independent":str(expected),
                         "main":str(obtained),"relative_error":error,
                         "tolerance":relative_tolerance})


@lru_cache(None)
def chain(q, n, cutoff):
    qq = mp.mpf(str(q))
    f = (qq**(2*n)-1)/(qq**(2*n+2)-1)
    endpoint = qq**(n-1)*(qq*qq-1)/(qq**(2*n+2)-1)
    # Independent source-column norm, not the analytic derivative used in main.
    column = [qq**(-j-2)*(1-qq**(-2*(n-j)))/(1-qq**(-2*n-2)) for j in range(n)]
    h2 = mp.fsum(g*g for g in column)
    eigen_low = qq*qq+1-2*qq*mp.cos(mp.pi/(n+1))
    eigen_high = qq*qq+1+2*qq*mp.cos(mp.pi/(n+1))
    M = mp.mpf(str(cutoff))*mp.sqrt(f)
    return dict(q=q,n=n,Lambda_eV=cutoff,G00=float(f),tau=float(endpoint/f),
                zeta=1.0,h2=float(h2),M_eV=float(M),
                min_mass_eV=float(M*mp.sqrt(eigen_low)),
                max_mass_eV=float(M*mp.sqrt(eigen_high)),
                condition_H=float(eigen_high/eigen_low),
                rH=float(endpoint/f)*P/(A*cutoff**2))


def cubic(kappa):
    """Monotone Newton solve, retaining relative precision as kappa tends to zero."""
    kappa = np.asarray(kappa, dtype=float)
    x = np.minimum(kappa,np.cbrt(kappa))
    for _ in range(48):
        residual = x*(1+x)**2-kappa
        updated = x-residual/((1+x)*(1+3*x))
        if np.all(np.abs(updated-x) <= 4e-16*np.maximum(np.abs(updated),1e-300)):
            x = updated
            break
        x = updated
    return x


def rematch(kappa0):
    if not 0 <= kappa0 < 27/8:
        raise ValueError("Outside declared one-time vacuum branch.")
    if kappa0 == 0:
        return 0.0, 1.0, 0.0
    lk = math.log(kappa0)
    root = brentq(lambda y:y+3*np.log1p(1.5*np.exp(y))-4*np.log1p(np.exp(y))-lk,
                  -750.0, 40.0, xtol=5e-14, rtol=1e-15)
    x0 = math.exp(root)
    delta = x0*(0.5+x0)/(1+x0)**2
    return x0,1-delta,delta


def response(q, n, cutoff, mg, phase):
    data = chain(q,n,cutoff).copy()
    cosine,sine = (1.,0.) if phase == "zero" else (0.,1.)
    rw = F*np.sqrt(U/2)/(P*mg)
    S0 = 3*P*mg**2
    S = S0*((cosine-rw)**2+sine*sine)
    Schi = 2*S0*rw*(cosine-rw)
    # Full auxiliary-field y jet is one third of a naive phase rotation.
    Sy = (2/3)*S0*rw*sine
    kappa = 4*data["h2"]*S/data["M_eV"]**4
    x = cubic(kappa)
    static_slope = -x/(1+x)
    x0,alpha,delta = rematch(4*data["h2"]*S0/data["M_eV"]**4)
    rematched_x = cubic(kappa/alpha**3)
    rematched_slope = (delta-alpha*rematched_x)/(alpha*(1+rematched_x))
    gradient = np.hypot(Schi,Sy)/(2*U)
    ratio = np.abs(static_slope)*gradient
    ratio_matched = np.abs(rematched_slope)*gradient
    massgap = data["min_mass_eV"]
    Fgap = math.sqrt(S0)/massgap**2
    mqgap = MQ/massgap
    maxgap = data["max_mass_eV"]/KK
    gapok = (Fgap<=0.1 and mqgap<=0.1 and maxgap<=0.1)
    data.update(mG_eV=mg,phase=phase,
        hidden_self_mass_leading_eV=2*math.sqrt(S0)/cutoff,
        clockoff_my2_leading_eV2=4*rho/(3*P)+12*A*mg**2*(data["rH"]-2),
        max_x=float(np.max(x)), max_static_force_ratio=float(np.max(ratio)),
        max_once_rematched_force_ratio=float(np.max(ratio_matched)),
        alpha_delta=delta,
        max_abs_static_dVchi_eV4=float(np.max(np.abs(static_slope*Schi))),
        max_abs_static_dVy_eV4=float(np.max(np.abs(static_slope*Sy))),
        FX_over_min_mass2=Fgap,MQ_over_min_mass=mqgap,max_mass_over_KK=maxgap,
        gap_conditions_pass=gapok,
        static_response_and_gap_conditions_pass=gapok and np.max(ratio)<=0.1,
        once_rematched_response_and_gap_conditions_pass=gapok and np.max(ratio_matched)<=0.1)
    return data


def read_table(name):
    with (BASE/"results"/name).open() as source:
        return list(csv.DictReader(source))


tables = {name:read_table(name) for name in ["chain_matching.csv","response_scan.csv",
                                          "representatives.csv","response_limits.csv"]}
for name,count in [("chain_matching.csv",777),("response_scan.csv",390),
                   ("representatives.csv",36),("response_limits.csv",6)]:
    compare("count:"+name,"rows",count,len(tables[name]),0)
compare("input:trajectory","rows",1025,len(chi),0)
compare("input:trajectory","first_chi",parent["chi_first"],chi[0],0)
compare("input:trajectory","last_chi",parent["chi_last"],chi[-1],0)
matching_keys = {(float(r["q"]),int(r["n"]),float(r["Lambda_eV"])) for r in tables["chain_matching.csv"]}
compare("coverage:matching","exact_grid",
        matching_keys == {(q,n,1e12) for q in [2.,3.,5.] for n in range(2,261)}, True,0)
scan_keys = {(float(r["q"]),float(r["Lambda_eV"]),round(math.log10(float(r["mG_eV"])),8),r["phase"])
             for r in tables["response_scan.csv"]}
compare("coverage:scan","exact_grid",
        scan_keys == {(q,c,-20+i/4,"zero") for q in [2.,3.,5.] for c in [1e12,2e12] for i in range(65)},True,0)
rep_keys = {(float(r["q"]),float(r["Lambda_eV"]),round(math.log10(float(r["mG_eV"]))),r["phase"])
            for r in tables["representatives.csv"]}
compare("coverage:representatives","exact_grid",
        rep_keys == {(q,c,lm,phase) for q in [2.,3.,5.] for c in [1e12,2e12]
                     for lm in [-6,-14,-15] for phase in ["zero","quadrature"]},True,0)

for index,row in enumerate(tables["chain_matching.csv"]):
    expected = chain(float(row["q"]),int(row["n"]),float(row["Lambda_eV"]))
    for field, value in expected.items():
        compare(f"matching:{index}",field,value,row[field])

selected = {}
for cutoff in [1e12,2e12]:
    for q in [2.,3.,5.]:
        candidates = [chain(q,n,cutoff) for n in range(2,261)]
        valid = [c for c in candidates if c["rH"]>2]
        selected[q,cutoff] = min(valid,key=lambda c:abs(c["rH"]-4))["n"]

classification_counts = {}
for filename in ["response_scan.csv","representatives.csv"]:
    group = filename.removesuffix(".csv")
    classification_counts[group] = {"gap":0,"static_and_gap":0,"once_rematched_and_gap":0}
    for index,row in enumerate(tables[filename]):
        q,n,cutoff,mg = float(row["q"]),int(row["n"]),float(row["Lambda_eV"]),float(row["mG_eV"])
        compare(f"{group}:{index}","selected_n",selected[q,cutoff],n,0)
        expected = response(q,n,cutoff,mg,row["phase"])
        for field,value in expected.items():
            compare(f"{group}:{index}",field,value,row[field])
        for label,field in [("gap","gap_conditions_pass"),("static_and_gap","static_response_and_gap_conditions_pass"),
                            ("once_rematched_and_gap","once_rematched_response_and_gap_conditions_pass")]:
            classification_counts[group][label] += int(expected[field])

independent_limits = []
for index,row in enumerate(tables["response_limits.csv"]):
    q,n,cutoff = float(row["q"]),int(row["n"]),float(row["Lambda_eV"])
    compare(f"limits:{index}","selected_n",selected[q,cutoff],n,0)
    item = {"q":q,"n":n,"Lambda_eV":cutoff}
    for result_field,response_field in [
        ("mG_static_10percent_eV","max_static_force_ratio"),
        ("mG_once_rematched_10percent_eV","max_once_rematched_force_ratio")]:
        root = brentq(lambda lm:math.log10(response(q,n,cutoff,10**lm,"zero")[response_field])+1,
                      -20.,-6.,xtol=2e-13,rtol=1e-15)
        item[result_field] = 10**root
        compare(f"limits:{index}",result_field,item[result_field],row[result_field])
    item["MQ_over_min_mass"] = MQ/chain(q,n,cutoff)["min_mass_eV"]
    compare(f"limits:{index}","MQ_over_min_mass",item["MQ_over_min_mass"],row["MQ_over_min_mass"])
    independent_limits.append(item)

# High-precision independent roots at fixed S/S0, including otherwise invisible alpha shifts.
high_precision = []
for q,cutoff,mg in [(3.,1e12,1e-6),(3.,1e12,1e-14),(3.,2e12,1e-14),(2.,1e12,1e-20),(2.,1e12,1e-4)]:
    n = selected[q,cutoff]
    c = chain(q,n,cutoff)
    S0 = mp.mpf(3)*mp.mpf(str(P))*mp.mpf(str(mg))**2
    kp = 4*mp.mpf(str(c["h2"]))*S0/mp.mpf(str(c["M_eV"]))**4
    xs = mp.findroot(lambda x:x*(1+x)**2-kp,(kp/2,kp))
    xr = mp.findroot(lambda x:x*(1+mp.mpf("1.5")*x)**3/(1+x)**4-kp,(kp/2,kp*2))
    ar = (1+mp.mpf("1.5")*xr)/(1+xr)**2
    xd,ad,dd = rematch(float(kp))
    compare(f"precision:q{q},cutoff{cutoff},mg{mg}","static_x",float(xs),float(cubic(float(kp))))
    compare(f"precision:q{q},cutoff{cutoff},mg{mg}","vacuum_x",float(xr),xd)
    compare(f"precision:q{q},cutoff{cutoff},mg{mg}","alpha_delta",float(1-ar),dd)
    leading_slope = -xr/(2+3*xr)
    compare(f"precision:q{q},cutoff{cutoff},mg{mg}","slope_identity",float(leading_slope),
            (dd-ad*xd)/(ad*(1+xd)))
    high_precision.append({"q":q,"n":n,"Lambda_eV":cutoff,"mG_eV":mg,
                           "static_x":mp.nstr(xs,45),"vacuum_x":mp.nstr(xr,45),
                           "alpha_delta":mp.nstr(1-ar,45),"remaining_slope":mp.nstr(leading_slope,45)})

input_files = [inputs_path,trajectory_path,BASE/"protocol.md",
               BASE/"results/inputs.json"]+[BASE/"results"/name for name in tables]
hashes = {str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest() for path in input_files}
payload = {
    "scope":"Independent endpoint, local static response, once-only leading vacuum rematch and classification comparison; no main code read or imported.",
    "code_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    "passed":sum(checks),"total":len(checks),"failures":failures,
    "relative_tolerance":1e-10,"zero_expected_tolerance":0,
    "high_precision_digits":mp.mp.dps,
    "inputs":{"P":P,"F":F,"A":A,"rho":rho,"Ui":Ui,"MQ_eV":MQ,"KK_eV":KK,
              "chi_count":len(chi),"chi_first":float(chi[0]),"chi_last":float(chi[-1])},
    "source_hashes":hashes,"aggregate":aggregate,
    "selected_chains":[{"q":q,"Lambda_eV":c,"n":n} for (q,c),n in selected.items()],
    "classification_counts":classification_counts,
    "independent_limits":independent_limits,
    "high_precision_checks":high_precision,
    "execution_notes":["First comparison passed 21729/21729. Added exact-zero checks and three explicit grid-coverage checks before final archive; thresholds were tightened, not loosened."]+execution_notes,
}
destination = BASE/"results/independent_comparison.json"
destination.write_text(json.dumps(payload,ensure_ascii=False,indent=2)+"\n")
print(json.dumps({"passed":payload["passed"],"total":payload["total"],
                  "failures":failures[:20],"classification_counts":classification_counts},indent=2))
raise SystemExit(0 if not failures else 1)
