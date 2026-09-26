#!/usr/bin/env python3
"""Verify archive identity and claim scope; not physical or empirical validation."""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
from decimal import Decimal, InvalidOperation
import hashlib
import json
from pathlib import Path
import platform
import re
import subprocess

import matplotlib
import mpmath
import numpy
import scipy
import sympy


ROOT = Path(__file__).resolve().parents[1]
REPO = ROOT.parents[1]
REPORTS = REPO / "reports"
MANIFEST = ROOT / "materials_manifest.json"
AUDIT = ROOT / "results/archive_audit.json"
BASELINE = "fbfae9851ce0c38452d18d2ec8d70b1d840e6c6d"
ROOT_DOCS = ("README.md", "WORK_PLAN.md", "RESEARCH_DIRECTION.md")
ROW_COUNTS = {
    "endpoint_scaling.csv": 24,
    "mode_profile.csv": 128,
    "physical_valley.csv": 101,
    "source_stress.csv": 305,
    "local_transient.csv": 201,
}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path):
    return json.loads(path.read_text())


def write(path, value):
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False) + "\n")


def inventory():
    paths = [p for p in ROOT.rglob("*") if p.is_file()
             and p not in (MANIFEST, AUDIT) and "__pycache__" not in p.parts]
    paths += [p for p in REPORTS.glob("light_modulus*_20260926*") if p.is_file()]
    paths += [REPO / p for p in ROOT_DOCS]
    paths += [REPO / "experiments/threshold_rg_v01/results/inputs.json"]
    return sorted(set(paths))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-manifest", action="store_true")
    args = parser.parse_args()
    if args.write_manifest:
        write(MANIFEST, {
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "stage": "light_modulus_v01",
            "baseline_git_commit": BASELINE,
            "material_type": "Conditional theoretical derivation and deterministic numerical consistency audit.",
            "question": "Does retaining the original Higgs light modulus preserve a finite controlled stationary background?",
            "input": "Original rank-n gauge chain, canonical paired Higgs fields, local endpoint norm currents and declared finite hidden source.",
            "output": "Quotient geometry, physical light-endpoint couplings, conditional finite-F heavy valley, stationary obstruction, local Hessian and bounded local transient.",
            "authorship": "AI-assisted internal derivations and independent implementations; no external peer review or new observations.",
            "source_read_record": "reports/light_modulus_geometry_20260926_sources.json",
            "versions": {
                "python": platform.python_version(), "numpy": numpy.__version__,
                "scipy": scipy.__version__, "sympy": sympy.__version__,
                "mpmath": mpmath.__version__, "matplotlib": matplotlib.__version__,
            },
            "limitations": [
                "The no-finite-interior-stationary-point result applies only on the stated rigid X=Za=0, no-clock-superpotential, positive-metric branch with nonzero Higgs pairs.",
                "The algebraic large-field continuation can meet kinetic or EFT boundaries; it is not an established physical infinite runaway.",
                "Zero-source quotient geometry and leading local Hessians do not replace all-order finite-F geometry or vacuum particle masses.",
                "The fixed-source homogeneous Minkowski local transient is not cosmic time, an observational forecast, or a SUGRA lifetime.",
                "No complete SUGRA embedding, new light-modulus stabilization, quantum correction, or other-X/Za-branch exclusion.",
                "No derivation of 1/ln²t driving or real physical interactions from Riemann structure.",
            ],
            "retained_failures": [
                "Initial stationarity audit 1423/1435: 12 second-order finite-difference truncation errors; original code and JSON retained.",
                "Fourth-order centered difference fixed the numerical check with the same step, samples and tolerances.",
                "Absence of a finite internal stationary point is a retained physical failure, not a failed implementation test.",
            ],
            "files": [{"path": str(p.relative_to(REPO)), "bytes": p.stat().st_size,
                       "sha256": sha(p)} for p in inventory()],
            "exclusions": ["manifest/audit self-reference", "Python caches", "ignored fetched-page caches"],
        })

    manifest = read(MANIFEST)
    checks = []

    def check(name, passed, detail=""):
        checks.append({"name": name, "passed": bool(passed), "detail": detail})

    def verify_hashes(values, label, base=REPO):
        for rel, expected in values.items():
            path = base / rel
            check(label + ":" + rel, path.is_file() and sha(path) == expected)

    check("baseline_record", manifest["baseline_git_commit"] == BASELINE)
    for item in manifest["files"]:
        path = REPO / item["path"]
        check("identity:" + item["path"], path.is_file()
              and path.stat().st_size == item["bytes"] and sha(path) == item["sha256"])
        if path.suffix == ".md":
            broken = []
            for link in re.findall(r'!?\[[^\]]*\]\(([^)\n]+)\)', path.read_text()):
                if "://" in link or link.startswith(("#", "mailto:")):
                    continue
                dest = (path.parent / link.split("#", 1)[0].strip("<>")).resolve()
                if not dest.exists() and dest != AUDIT.resolve():
                    broken.append(link)
            check("local_links:" + item["path"], not broken, broken)
    check("inventory_complete", {str(p.relative_to(REPO)) for p in inventory()}
          == {item["path"] for item in manifest["files"]})

    summary = read(ROOT / "results/summary.json")
    check("main_checks", summary["checks_passed"] == summary["check_count"] == 549
          and len(summary["checks"]) == 549 and all(c["passed"] for c in summary["checks"]))
    verify_hashes(summary["source_hashes"], "main_source")
    check("main_log_completion", "Main consistency 549 / 549" in (ROOT / "logs/main_run.log").read_text())
    check("main_stationary_failure_retained", "No finite interior stationary point" in summary["retained_physical_failure"]
          and "X=Za=0" in summary["retained_physical_failure"])
    check("nonstationary_hessian_scope", "nonstationary" in summary["benchmark"]["hessian_scope"]
          and "not vacuum particle masses" in summary["benchmark"]["hessian_scope"])

    specifications = [
        ("geometry", 507, "checks_total", "checks_passed", "source_hashes"),
        ("stationarity", 1435, "total", "passed", "input_sha256"),
        ("response", 110, "total", "passed", None),
    ]
    results = {}
    for stem, expected, total_key, passed_key, sources_key in specifications:
        path = REPORTS / f"light_modulus_{stem}_20260926_checks.json"
        value = read(path)
        results[stem] = value
        check(stem + "_checks", value[total_key] == value[passed_key] == expected
              and len(value["checks"]) == expected and all(c["passed"] for c in value["checks"]))
        check(stem + "_code_identity", sha(path.with_suffix(".py")) == value["script_sha256"])
        if sources_key:
            verify_hashes(value[sources_key], stem + "_input", REPORTS)

    initial = read(REPORTS / "light_modulus_stationarity_20260926_initial_checks.json")
    check("stationarity_initial_failure_retained", initial["total"] == 1435 and initial["passed"] == 1423
          and len(initial["failed"]) == 12 and all("_envelope" in c["name"] for c in initial["failed"]))
    check("stationarity_initial_code_identity", sha(REPORTS / "light_modulus_stationarity_20260926_initial_checks.py")
          == initial["script_sha256"])
    old_code = (REPORTS / "light_modulus_stationarity_20260926_initial_checks.py").read_text()
    new_code = (REPORTS / "light_modulus_stationarity_20260926_checks.py").read_text()
    check("stationarity_thresholds_and_step_retained", all(s in old_code and s in new_code
          for s in ("width=2e-4", "derr<2e-8", "res<2e-8", "rerr<2e-8 and verr<2e-9")))
    check("stationarity_amplitude_solver_coverage", results["stationarity"]["nonlinear_solve_count"] == 216
          and len(results["stationarity"]["nonlinear_rows"]) == 216
          and all(r["solver_success"] for r in results["stationarity"]["nonlinear_rows"]))
    check("stationarity_sample_metrics_positive", results["stationarity"]["min_hidden_metric"] > 0
          and results["stationarity"]["min_clock_schur"] > 0)

    total_rows, total_cells, numeric_cells, undefined_cells = 0, 0, 0, 0
    for name, expected in ROW_COUNTS.items():
        path = ROOT / "results" / name
        with path.open() as stream:
            reader = csv.DictReader(stream)
            names = reader.fieldnames
            rows = list(reader)
        total_rows += len(rows)
        total_cells += len(rows) * len(names)
        valid = all(set(row) == set(names) for row in rows)
        try:
            for row in rows:
                for column, value in row.items():
                    if value == "undefined_zero_source":
                        undefined_cells += 1
                        valid = valid and name == "source_stress.csv" and column == "V_over_S" and Decimal(row["kappa"]) == 0
                    else:
                        numeric_cells += 1
                        valid = valid and Decimal(value).is_finite()
        except (InvalidOperation, TypeError, ValueError):
            valid = False
        check("table:" + name, len(rows) == expected and valid and b"\r" not in path.read_bytes(),
              {"rows": len(rows), "columns": len(names), "finite_numeric_or_explicit_undefined_cells": valid})
        check("row_count_record:" + name, summary["row_counts"][path.stem] == len(rows))
    check("table_total_coverage", total_rows == 759 and total_cells == 5583
          and numeric_cells == 5522 and undefined_cells == 61,
          {"rows": total_rows, "cells": total_cells, "numeric_cells": numeric_cells,
           "explicit_undefined_cells": undefined_cells})

    comparison_path = REPORTS / "light_modulus_output_comparison_20260926.json"
    comparison = read(comparison_path)
    # This schema is fixed to the independent comparison record at archive freeze.
    verify_comparison(comparison, comparison_path, check, verify_hashes)

    for ext, magic in (("png", b"\x89PNG"), ("pdf", b"%PDF")):
        path = ROOT / "figures" / ("light_modulus_diagnostics." + ext)
        check("figure:" + ext, path.is_file() and path.stat().st_size > 1000
              and path.read_bytes().startswith(magic))

    sources = read(REPORTS / "light_modulus_geometry_20260926_sources.json")
    check("source_scope_accurate", sources["arxiv"] == "0904.0370v2"
          and "Abstract only" in sources["read_scope"]
          and sources["human_read_attestation"] is False
          and "chain-specific" in sources["does_not_support"])
    response_scope = results["response"]["conditional_transient"]["scope"]
    protocol = (ROOT / "protocol.md").read_text()
    check("transient_scope_restricted", "Minkowski" in response_scope
          and "not a cosmological solution or SUGRA lifetime" in response_scope
          and "不是实验界或宇宙寿命" in protocol)
    check("rigid_branch_scope_recorded", "X=Za=0" in protocol and "刚性分支" in protocol
          and "不排除" in protocol and "1/ln²t仍为基础响应设想" in protocol)

    def permitted(path):
        return path in ROOT_DOCS or path.startswith("experiments/light_modulus_v01/") or (
            path.startswith("reports/light_modulus") and "_20260926" in path)

    changed = subprocess.check_output(["git", "diff", "--name-only", "-z", BASELINE, "--"], cwd=REPO).decode().split("\0")
    disallowed = [p for p in changed if p and not permitted(p)]
    check("historical_files_and_manuscript_unchanged", not disallowed, disallowed)
    untracked = subprocess.check_output(["git", "ls-files", "--others", "--exclude-standard", "-z"], cwd=REPO).decode().split("\0")
    disallowed_untracked = [p for p in untracked if p and not permitted(p)]
    check("untracked_files_within_authorized_stage", not disallowed_untracked, disallowed_untracked)

    result = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "manifest_sha256": sha(MANIFEST), "audit_script_sha256": sha(Path(__file__)),
        "baseline_git_commit": BASELINE,
        "passed": sum(c["passed"] for c in checks), "total": len(checks), "checks": checks,
        "scope": "Archive identity, source scope, completeness and retained checks; not physical or empirical validation.",
    }
    write(AUDIT, result)
    print("Archive consistency", result["passed"], "/", result["total"])
    for item in checks:
        if not item["passed"]:
            print(item)
    raise SystemExit(0 if result["passed"] == result["total"] else 1)


def verify_comparison(value, path, check, verify_hashes):
    """Bind the complete comparison coverage to its code and input identities."""
    check("independent_comparison_checks", value["passed"] == value["total"] == 12
          and len(value["checks"]) == 12 and all(c["passed"] for c in value["checks"]))
    check("independent_comparison_full_coverage", value["csv_rows_compared"] == 759
          and value["csv_cells_compared"] == 5583 and value["csv_numeric_cells"] == 5522
          and value["csv_text_cells"] == 61 and value["summary_numeric_cells_compared"] == 30)
    check("independent_comparison_code_identity", sha(path.with_suffix(".py")) == value["script_sha256"])
    verify_hashes(value["input_hashes"], "independent_comparison_input")
    check("independent_comparison_scope", "after reading protocol/schema/main code" in value["scope"]
          and "not blind or external peer review" in value["scope"])


if __name__ == "__main__":
    main()
