#!/usr/bin/env python3
"""Small independent Fourier/CIC diagnostics; never run the production matrix.

The expected values are closed-form single-mode fields, CIC's linear response,
and the six-mode first-shell power. No time integration or spectrum fit occurs.
"""

from __future__ import annotations

import hashlib
import json
import platform
import time
import traceback
import types
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import scipy


ROOT = Path(__file__).resolve().parents[1]
PM_PATH = ROOT / "code/pm.py"
REPORT_PATH = ROOT / "reports/independent_pm_checks.json"


def sha256(data):
    return hashlib.sha256(data).hexdigest()


def run_checks(pm, checks, details):
    def upper(name, value, limit):
        checks.append({"name": name, "value": float(value),
                       "strict_upper_threshold": float(limit),
                       "passed": bool(np.isfinite(value) and value < limit)})

    # An oblique cosine fixes the signs of all three half-cell phase factors.
    mode = np.array([1, 2, -1])
    amplitude = 1e-3
    coefficients = np.zeros((32, 32, 32), dtype=complex)
    coefficients[tuple(mode % 32)] = amplitude / 2
    coefficients[tuple((-mode) % 32)] = amplitude / 2
    modes = {"coefficients": coefficients, "metadata": {}}
    wave = 2 * np.pi * mode
    momentum_factor = .02**2 * 200.0 * .99
    initial_rows = {}
    for n in (32, 64):
        initial = pm.initial_conditions(modes, n, .02, 200.0, .99)
        selected = np.array([0, 1, n + 2, n*n + 3, n**3 // 3, n**3 - 1])
        points = initial["lagrangian"][selected]
        # Direct closed-form evaluation, with no FFT or PM helper.
        expected = -amplitude * np.sin(points @ wave)[:, None] * wave / (wave @ wave)
        error = float(np.max(np.abs(initial["displacement"][selected] - expected)))
        imaginary = initial["metadata"]["max_ifft_imaginary_displacement"]
        momentum_error = float(np.max(np.abs(initial["momenta"] -
                                            momentum_factor * initial["displacement"])))
        initial_rows[str(n)] = {
            "selected_flat_indices": selected.tolist(), "positions": points.tolist(),
            "expected_displacements": expected.tolist(),
            "actual_displacements": initial["displacement"][selected].tolist(),
            "direct_evaluation_max_abs_error": error,
            "ifft_max_imaginary": float(imaginary),
            "common_momentum_relation_max_abs_error": momentum_error,
        }
        upper(f"cell_center_direct_mode:n={n}", error, 1e-16)
        upper(f"cell_center_imaginary:n={n}", imaginary, 1e-16)
        upper(f"common_reference_momentum:n={n}", momentum_error, 1e-17)
    details["shared_continuous_field"] = {
        "source_grid": 32, "mode": mode.tolist(), "density_cosine_amplitude": amplitude,
        "momentum_factor": momentum_factor,
        "expected_displacement": "-A sin(2*pi*m.q) (2*pi*m)/|2*pi*m|^2, in box units",
        "note": "32 and 64 cell centers are not nested; each is checked at its own physical positions.",
        "grids": initial_rows,
    }

    cic_rows = {}
    density_amplitude = 2e-6
    for n in (16, 32):
        centers = (np.arange(n) + .5) / n
        q = np.stack(np.meshgrid(centers, centers, centers, indexing="ij"), axis=-1).reshape(-1, 3)
        displacement = np.zeros_like(q)
        displacement[:, 0] = -density_amplitude / (2*np.pi) * np.sin(2*np.pi*q[:, 0])
        force, diagnostics = pm.PMGrid(n).force(q + displacement)
        projected = float(np.sum(force[:, 0] * displacement[:, 0]) /
                          np.sum(displacement[:, 0]**2))
        kh = 2*np.pi/n
        expected = float(np.sin(kh) / kh)
        force_residual = float(np.max(np.abs(force - expected*displacement)) /
                               np.max(np.abs(displacement)))
        transverse = float(np.max(np.abs(force[:, 1:])))
        cic_rows[str(n)] = {
            "analytic_F": expected, "projected_F": projected,
            "projected_F_abs_error": abs(projected - expected),
            "full_force_residual_over_displacement_amplitude": force_residual,
            "transverse_force_max_abs": transverse, "force_diagnostics": diagnostics,
        }
        upper(f"CIC_projected_linear_response:n={n}", abs(projected - expected), 1e-8)
        upper(f"CIC_finite_amplitude_remainder:n={n}", force_residual, 3*density_amplitude)
        upper(f"CIC_raw_mass:n={n}", diagnostics["mass_relative_error"], 1e-12)
        upper(f"CIC_transverse_force:n={n}", transverse, 1e-14)
    details["CIC_linear_response"] = {
        "initial_density_amplitude": density_amplitude,
        "geometry": "mesh nodes j/N, initial particle positions (j+1/2)/N",
        "analytic_response": "F = sin(kh)/(kh), kh = 2*pi/N",
        "remainder_note": "The full finite-amplitude force has higher-order terms; the projected fundamental checks F.",
        "grids": cic_rows,
    }

    n, box, amplitude = 32, 50.0, .02
    delta = np.broadcast_to((amplitude*np.cos(2*np.pi*np.arange(n)/n))[:, None, None], (n, n, n)).copy()
    fundamental = 2*np.pi/box
    result = pm.power_spectrum(delta, box, [.99*fundamental, 1.01*fundamental])
    # This shell has six modes; only the +x/-x pair has power L^3 A^2/4 each.
    expected_power = box**3 * amplitude**2 / 12
    relative_error = float(abs(result["pk"][0]/expected_power - 1))
    upper("single_axis_cosine_shell_power", relative_error, 1e-12)
    checks.append({"name": "first_shell_mode_count", "value": int(result["count"][0]),
                   "expected": 6, "passed": bool(result["count"][0] == 6)})
    details["power_normalization"] = {
        "grid": n, "box_mpc_h": box, "density_cosine_amplitude": amplitude,
        "actual_shell_power": float(result["pk"][0]), "expected_shell_power": expected_power,
        "relative_error": relative_error, "full_fft_mode_count": int(result["count"][0]),
        "formula": "L^3 A^2/12 for the six-mode first shell of a single-axis cosine",
    }


def main():
    started = time.perf_counter()
    # Execute exactly these source bytes, so the recorded hash identifies the
    # implementation actually checked even if another worker edits the file.
    source = PM_PATH.read_bytes()
    source_hash = sha256(source)
    pm = types.ModuleType("independent_pm_snapshot")
    pm.__file__ = str(PM_PATH)
    checks, details = [], {}
    try:
        exec(compile(source, str(PM_PATH), "exec"), pm.__dict__)
        run_checks(pm, checks, details)
    except Exception:
        checks.append({"name": "execution_exception", "passed": False,
                       "traceback": traceback.format_exc()})
    final_source_hash = sha256(PM_PATH.read_bytes())
    checks.append({"name": "source_unchanged_during_checks",
                   "passed": source_hash == final_source_hash})
    failures = [row for row in checks if not row["passed"]]
    report = {
        "status": "FAIL" if failures else "PASS",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "scope": "Independent small Fourier/CIC diagnostics only; no production matrix, time integration, observed-data fit or statistical detection.",
        "code_sha256": {"pm.py": source_hash, "independent_pm_checks.py": sha256(Path(__file__).read_bytes())},
        "pm_disk_sha256_after_checks": final_source_hash,
        "protocol_sha256": sha256((ROOT / "protocol.md").read_bytes()),
        "environment": {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__},
        "checks_passed": len(checks) - len(failures), "checks_total": len(checks),
        "checks": checks, "failures": failures, "detail": details,
        "wall_seconds": time.perf_counter() - started,
    }
    if REPORT_PATH.exists():
        previous = json.loads(REPORT_PATH.read_text())
        history = previous.pop("previous_runs", [])
        report["previous_runs"] = history + [previous]
    REPORT_PATH.parent.mkdir(parents=True, exist_ok=True)
    REPORT_PATH.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    print(json.dumps({"status": report["status"], "checks_passed": report["checks_passed"],
                      "checks_total": report["checks_total"], "pm_sha256": source_hash,
                      "report": str(REPORT_PATH)}, indent=2))
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
