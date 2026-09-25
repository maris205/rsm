#!/usr/bin/env python3
"""Independent arbitrary-precision supertrace; does not integrate a cosmology.

The potential is evaluated by direct subtraction of the three determinants.
No small-splitting series from the production implementation is imported.
Input masses and background samples are interpreted as exact decimal inputs;
high arithmetic precision does not improve their physical/input precision.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
from pathlib import Path
import platform
from datetime import datetime, timezone

import mpmath as mp

ROOT = Path(__file__).resolve().parents[1]
REPO = ROOT.parents[1]
PREVIOUS = REPO / "experiments/charged_threshold_v01/results"


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_inputs():
    return json.loads((PREVIOUS / "inputs.json").read_text(), parse_float=str)


def mpval(x):
    return mp.mpf(str(x))


def field(inputs, u, mass_eV, eta_soft, f_split=True, mu_ratio="1"):
    """Evaluate masses and chi derivatives, with g_i exactly equal to one."""
    ci = mpval(inputs["chi_i"])
    u = mpval(u)
    chi = ci + u
    mi = mpval(mass_eV)
    eta = mpval(eta_soft)
    fsq = mpval(inputs["F_squared_eV2"])
    unit = mpval(inputs["potential_unit_eV4"])
    wi = mpval(inputs["W_i"])
    g = (1 + (ci / chi) ** 2) / 2
    gp = -ci**2 / chi**3
    gpp = 3 * ci**2 / chi**4
    mass = mi * mp.sqrt(g)
    massp = mi * gp / (2 * mp.sqrt(g))
    tree = unit * wi * mp.exp(-2 * u)
    if f_split:
        h = mp.sqrt(tree) * mp.sqrt(2 / fsq) * abs(massp)
        hp = h * (-1 + gpp / gp - gp / (2 * g))
    else:
        h = hp = mp.mpf("0")
    s = mass**2
    sp = mi**2 * gp
    d = 8 * mp.pi**2 * unit * eta / mi**2
    plus, minus = s + d + h, s + d - h
    if minus <= 0:
        raise ValueError("Tachyonic scalar mass squared in independent benchmark")
    return dict(chi=chi, u=u, g=g, gp=gp, gpp=gpp, s=s, sp=sp,
                d=d, h=h, hp=hp, plus=plus, minus=minus,
                mu2=(mi * mpval(mu_ratio))**2,
                tree=tree, tree_prime=-2 * tree)


def raw_supertrace(p):
    """Direct MSbar determinants and their analytic first derivatives."""
    def f(x):
        return x**2 * (mp.log(x / p["mu2"]) - mp.mpf("1.5"))

    def fp(x):
        return 2 * x * (mp.log(x / p["mu2"]) - 1)

    potential = (f(p["plus"]) + f(p["minus"]) - 2 * f(p["s"])) / (32 * mp.pi**2)
    derivative = (fp(p["plus"]) * (p["sp"] + p["hp"])
                  + fp(p["minus"]) * (p["sp"] - p["hp"])
                  - 2 * fp(p["s"]) * p["sp"]) / (32 * mp.pi**2)
    return potential, derivative


def evaluate(inputs, u, mass_eV, eta_soft="0", f_split=True, today_u=None, mu_ratio="1"):
    p = field(inputs, u, mass_eV, eta_soft, f_split, mu_ratio)
    initial = field(inputs, "0", mass_eV, eta_soft, f_split, mu_ratio)
    v, vp = raw_supertrace(p)
    vi, _ = raw_supertrace(initial)
    p.update(potential_raw=v, potential_prime=vp, potential_initial=vi,
             potential_subtracted=v - vi,
             force_ratio=abs(vp / p["tree_prime"]),
             splitting_ratio=p["h"] / p["s"],
             soft_ratio=p["d"] / p["s"])
    if today_u is not None:
        now = field(inputs, today_u, mass_eV, eta_soft, f_split, mu_ratio)
        alpha0 = mpval(inputs["alpha_reference"])
        weighted_log = (mp.log(p["plus"] / now["plus"])
                        + mp.log(p["minus"] / now["minus"])
                        + 4 * mp.log(p["s"] / now["s"])) / 3
        b = 1 - alpha0 * weighted_log / (4 * mp.pi)
        p.update(gauge_denominator=b, delta_alpha=1 / b - 1)
    return p


def relative_error(a, b):
    return abs(a - b) / max(abs(a), abs(b), mp.mpf("1e-1000"))


def decimal(x):
    return mp.nstr(x, 110, min_fixed=-5, max_fixed=8)


def canonical_cases():
    """Predeclared physical grid in protocol.md; no fitting or selection."""
    return [{"case": f"m{mass}_eta{eta}_mu{mu}", "mass_eV": mass,
             "eta_soft": eta, "mu_ratio": mu, "f_split": True}
            for mass in ("1", "1e6", "1e11")
            for eta in ("0", "1", "100", "1e4", "1e6")
            for mu in ("0.5", "1", "2")]


def compare_production(report, inputs):
    """Compare an imported public function, without sharing its expansion."""
    primary_path = ROOT / "code/run_protected.py"
    if not primary_path.exists():
        return
    spec = importlib.util.spec_from_file_location("protected_primary", primary_path)
    primary = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(primary)
    floats = json.loads((PREVIOUS / "inputs.json").read_text())
    pairs = {"U_eV4": "potential_raw", "Up_eV4": "potential_prime",
             "deltaU_eV4": "potential_subtracted", "s": "s", "sp": "sp", "d": "d"}
    maxima = {key: mp.mpf("0") for key in pairs}
    comparisons = []
    for record in report["records"]:
        ref = record["values"]
        got = primary.evaluate_loop(floats["chi_i"] + float(ref["u"]),
                                    float(record["mass_eV"]), float(record["eta_soft"]),
                                    float(record["mu_ratio"]), floats)
        errors = {}
        for name, reference in pairs.items():
            target = mpval(ref[reference])
            error = relative_error(mpval(float(got[name])), target)
            errors[name] = decimal(error)
            maxima[name] = max(maxima[name], error)
        h2_error = relative_error(mpval(float(got["h2"])), mpval(ref["h"])**2)
        errors["h2"] = decimal(h2_error)
        maxima["h2"] = max(maxima.get("h2", mp.mpf("0")), h2_error)
        comparisons.append({"case": record["case"], "index": record["index"],
                            "relative_errors": errors})
    # Read the actually archived new trajectories to check the full threshold
    # output; these samples are additional to the frozen potential grid.
    trajectory_path = ROOT / "results/trajectories.csv"
    threshold_records = []
    threshold_max = mp.mpf("0")
    if trajectory_path.exists():
        with trajectory_path.open() as handle:
            trajectories = list(csv.DictReader(handle))
        labels = sorted(set(row["case"] for row in trajectories))
        for label in labels:
            subset = [row for row in trajectories if row["case"] == label]
            if len(subset) != 1025 or "delta_alpha" not in subset[0]:
                continue
            eta = label.split("_", 1)[1]
            today = subset[-1]["u"]
            for index in (0, 1, 256, 512, 768, 1024):
                ref = evaluate(inputs, subset[index]["u"], "1e11", eta,
                               True, today, "1")
                got = mpval(subset[index]["delta_alpha"])
                error = relative_error(got, ref["delta_alpha"])
                threshold_max = max(threshold_max, error)
                threshold_records.append({"case": label, "index": index,
                                          "reference_delta_alpha": decimal(ref["delta_alpha"]),
                                          "production_delta_alpha": decimal(got),
                                          "relative_error": decimal(error)})
    passed = all(x < mp.mpf("1e-8") for x in maxima.values())
    if threshold_records:
        passed = passed and threshold_max < mp.mpf("1e-8")
    result = {
        "scope": "Compare direct determinants with the public production series; relative errors use each actual nonzero physical quantity, with exact zero required at a shared subtraction reference.",
        "source_sha256": {str(primary_path.relative_to(REPO)): sha(primary_path),
                          "independent_precision.py": sha(Path(__file__))},
        "potential_points": len(comparisons),
        "max_relative_errors": {key: decimal(value) for key, value in maxima.items()},
        "threshold_points": len(threshold_records),
        "max_threshold_relative_error": decimal(threshold_max),
        "tolerance": "1e-8", "all_checks_pass": bool(passed),
        "comparisons": comparisons, "threshold_comparisons": threshold_records}
    if threshold_records:
        result["source_sha256"][str(trajectory_path.relative_to(REPO))] = sha(trajectory_path)
    (ROOT / "results/independent_comparison.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False) + "\n")
    print(json.dumps({key: result[key] for key in ("potential_points", "max_relative_errors",
                                                 "threshold_points", "max_threshold_relative_error",
                                                 "all_checks_pass")}, indent=2))
    if not passed:
        raise SystemExit("Independent comparison against production failed")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dps", type=int, default=320)
    parser.add_argument("--confirm-dps", type=int, default=400)
    args = parser.parse_args()
    if args.dps < 120 or args.confirm_dps <= args.dps:
        raise ValueError("Require at least 120 digits and a higher precision repeat")
    mp.mp.dps = args.confirm_dps
    inputs = read_inputs()
    cases = canonical_cases()
    with (PREVIOUS / "trajectories.csv").open() as handle:
        rows = [r for r in csv.DictReader(handle) if r["case"] == "chi_probe"]
    indices = sorted(set([0, 1, len(rows) // 4, len(rows) // 2,
                          3 * len(rows) // 4, len(rows) - 1]))
    today_u = rows[-1]["u"]
    records = []
    max_precision_error = mp.mpf("0")
    max_derivative_error = mp.mpf("0")
    compared = ("potential_raw", "potential_prime", "potential_initial",
                "potential_subtracted", "force_ratio", "delta_alpha")
    checks = []
    for case in cases:
        for index in indices:
            u = rows[index]["u"]
            values = []
            for precision in (args.dps, args.confirm_dps):
                with mp.workdps(precision):
                    values.append(evaluate(inputs, u, case["mass_eV"],
                                           case["eta_soft"], case["f_split"], today_u,
                                           case["mu_ratio"]))
            errors = {k: relative_error(values[0][k], values[1][k]) for k in compared}
            max_precision_error = max(max_precision_error, *errors.values())
            p = values[1]
            numeric_derivative = mp.diff(
                lambda offset: raw_supertrace(field(inputs, offset, case["mass_eV"],
                                                     case["eta_soft"], case["f_split"],
                                                     case["mu_ratio"]))[0],
                mpval(u))
            derivative_error = relative_error(numeric_derivative, p["potential_prime"])
            max_derivative_error = max(max_derivative_error, derivative_error)
            records.append(dict(case=case["case"], index=index, N=rows[index]["N"],
                                mass_eV=str(case["mass_eV"]), eta_soft=str(case["eta_soft"]),
                                mu_ratio=case["mu_ratio"],
                                f_split=case["f_split"],
                                values={k: decimal(v) for k, v in p.items()},
                                precision_relative_errors={k: decimal(v) for k, v in errors.items()},
                                differentiated_potential_prime=decimal(numeric_derivative),
                                derivative_relative_error=decimal(derivative_error)))
    unbroken = evaluate(inputs, rows[0]["u"], "1e11", "0", False, today_u)
    checks.append({"name": "exact_susy_potential_and_force_cancel",
                   "pass": unbroken["potential_raw"] == unbroken["potential_prime"] == 0})
    checks.append({"name": "exact_susy_charged_threshold_does_not_cancel",
                   "pass": unbroken["delta_alpha"] != 0})
    h4_diagnostics = []
    for mass in ("1", "1e6", "1e11"):
        p = evaluate(inputs, "0", mass)
        h4 = -p["h"]**4 / (192 * mp.pi**2 * p["s"]**2)
        ratio = p["potential_raw"] / h4
        leading_force = p["h"]**2 * p["sp"] / p["s"] / (16 * mp.pi**2)
        checks.append({"name": f"initial_h4_residual_m{mass}",
                       "pass": abs(ratio - 1) < mp.mpf("1e-60")})
        checks.append({"name": f"initial_F_force_nonzero_m{mass}",
                       "pass": p["potential_prime"] != 0 and
                               relative_error(p["potential_prime"], leading_force) < mp.mpf("1e-60")})
        h4_diagnostics.append({"mass_eV": mass, "raw_initial_potential": decimal(p["potential_raw"]),
                               "leading_h4_potential": decimal(h4),
                               "ratio_to_h4": decimal(ratio),
                               "initial_force": decimal(p["potential_prime"]),
                               "leading_initial_force": decimal(leading_force),
                               "h_over_s": decimal(p["h"] / p["s"])})
    checks.extend([
        {"name": "precision_repeat_relative_error_below_1e-100",
         "pass": bool(max_precision_error < mp.mpf("1e-100"))},
        {"name": "direct_differentiation_relative_error_below_1e-100",
         "pass": bool(max_derivative_error < mp.mpf("1e-100"))}])
    source_paths = [Path(__file__), PREVIOUS / "inputs.json",
                    PREVIOUS / "trajectories.csv"]
    if (ROOT / "protocol.md").exists():
        source_paths.append(ROOT / "protocol.md")
    report = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "scope": "Independent direct high-precision flat-space CW determinants and thresholds on saved tree chi trajectory; no background ODE or observational test.",
        "arithmetic_digits": [args.dps, args.confirm_dps],
        "input_precision_note": "Decimal inputs and sampled u are fixed; arithmetic digits are not physical precision.",
        "renormalization": "Fixed mu/m_i = 0.5, 1, 2 per case; remove only initial constant after evaluating full supertrace.",
        "g_definition": "g=(1+(chi_i/chi)^2)/2, so g(chi_i)=1 exactly.",
        "python": platform.python_version(), "mpmath": mp.__version__,
        "source_sha256": {str(p.relative_to(REPO)) if p.is_relative_to(REPO) else str(p): sha(p)
                          for p in source_paths},
        "max_precision_relative_error": decimal(max_precision_error),
        "max_derivative_relative_error": decimal(max_derivative_error),
        "checks": checks, "all_checks_pass": all(c["pass"] for c in checks),
        "unbroken_control": {k: decimal(unbroken[k]) for k in
                             ("potential_raw", "potential_prime", "delta_alpha")},
        "initial_h4_diagnostics": h4_diagnostics,
        "records": records}
    target = ROOT / "results/independent_precision.json"
    target.write_text(json.dumps(report, ensure_ascii=False, indent=2, allow_nan=False) + "\n")
    print(json.dumps({k: report[k] for k in ("arithmetic_digits", "max_precision_relative_error",
                                            "max_derivative_relative_error", "all_checks_pass")}, indent=2))
    print(f"Saved {len(records)} benchmark samples to {target}")
    if not report["all_checks_pass"]:
        raise SystemExit(1)
    compare_production(report, inputs)


if __name__ == "__main__":
    main()
