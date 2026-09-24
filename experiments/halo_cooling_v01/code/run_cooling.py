#!/usr/bin/env python3
"""Run the frozen halo-to-atomic-cooling diagnostic and its source checks.

Production cooling uses cooling.py exclusively. The archived Grackle file is
imported only inside validate_against_official(). No hydro or stars are evolved.
Do not execute until protocol.md is frozen. Optional --events reads the
explicit event interface described in attach_event_readouts().
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import platform
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import scipy

import cooling as c


ROOT = Path(__file__).resolve().parents[1]
OLD_BACKGROUND = ROOT.parent / "cosmic_bridge_v01" / "code" / "background.py"
BACKGROUND_CSV = {"reference": ROOT.parent / "cosmic_bridge_v01" / "results" / "background_eps0.csv",
                  "clock": ROOT.parent / "cosmic_bridge_v01" / "results" / "background_eps1e-4.csv"}
OFFICIAL_PATH = ROOT / "sources" / "grackle_primordial_equilibrium.py"
SOURCE_METADATA = ROOT / "sources" / "grackle_source_metadata.json"
COMMIT = "af7939494ce65007887ada7b98d1813df6843346"
PROTOCOL_SHA256 = "f3a0cc5abf21a4374439db756cbf8c4a3103bcabfa2bee6f9bd1921b28e608f2"
RELATIVE_TARGET = 1e-10
DENSITY_TARGET = 1e-12
ROOT_TARGET = 1e-8


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot import {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def jsonable(obj):
    if isinstance(obj, dict):
        return {key: jsonable(value) for key, value in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [jsonable(value) for value in obj]
    if isinstance(obj, np.ndarray):
        return jsonable(obj.tolist())
    if isinstance(obj, np.generic):
        return jsonable(obj.item())
    if isinstance(obj, float) and not np.isfinite(obj):
        return None
    return obj


def write_json(path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(jsonable(obj), ensure_ascii=False, indent=2,
                               allow_nan=False) + "\n", encoding="utf-8")


def write_csv(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: json.dumps(value, ensure_ascii=False)
                             if isinstance(value, (dict, list)) else value
                             for key, value in jsonable(row).items()})


class Checks:
    def __init__(self):
        self.items = []

    def condition(self, name, condition, **details):
        self.items.append({"name": name, "passed": bool(condition), **details})

    def close(self, name, actual, expected, rtol=RELATIVE_TARGET, atol=0.0,
              units="dimensionless"):
        a, b = np.broadcast_arrays(np.asarray(actual, dtype=float), np.asarray(expected, dtype=float))
        finite = bool(np.all(np.isfinite(a)) and np.all(np.isfinite(b)))
        absolute = np.abs(a - b)
        rel = np.divide(absolute, np.abs(b), out=np.full_like(absolute, np.nan), where=b != 0.0)
        finite_rel = rel[np.isfinite(rel)]
        max_rel = float(np.max(finite_rel)) if finite_rel.size else None
        passed = finite and bool(np.all(absolute <= atol + rtol * np.abs(b)))
        abs_fallback = (absolute > rtol * np.abs(b)) & (absolute <= atol + rtol * np.abs(b))
        self.condition(name, passed, max_absolute_error=float(np.max(absolute)),
                       max_relative_error_nonzero_expected=max_rel, relative_tolerance=rtol,
                       absolute_tolerance=atol, absolute_tolerance_units=units,
                       n_values=int(a.size), n_zero_expected=int(np.count_nonzero(b == 0.0)),
                       n_values_needing_absolute_tolerance=int(np.count_nonzero(abs_fallback)),
                       all_inputs_finite=finite)

    @property
    def passed(self):
        return all(item["passed"] for item in self.items)


def validate_against_official(checks):
    """Only validation calls the pinned third-party source; production never does."""
    official = load_module("archived_grackle_cie_validation_only", OFFICIAL_PATH)
    t = np.logspace(4.0, 8.0, 161)
    rate_map = {
        "C_HI": "GammaeHI", "C_HeI": "GammaeHeI", "C_HeII": "GammaeHeII",
        "A_HII": "alphaHII", "A_HeII_rad": "alphaHeII",
        "A_HeII_die": "alphad", "A_HeIII": "alphaHeIII",
    }
    for name, function in rate_map.items():
        checks.close(f"official_rate_{name}", c.rates_cen(t)[name],
                     getattr(official, function)(t, rates="cen"), units="cm^3 s^-1")
    process_map = dict(zip(c.PROCESS_NAMES, (
        "ceHI", "ceHeII", "ciHI", "ciHeI", "ciHeII", "rHII", "rHeII",
        "rHeIII", "drHeII", "freefree")))
    species_map = {
        "nHI_cm3": "nHI", "nHII_cm3": "nHII", "nHeI_cm3": "nHeI",
        "nHeII_cm3": "nHeII", "nHeIII_cm3": "nHeIII", "ne_cm3": "ne",
    }
    curve_rows = []
    for nh in (1e-6, 1e-4, 1e-2):
        own = c.atomic_cooling(t, nh)
        for name, function in species_map.items():
            checks.close(f"official_species_{name}_nH{nh:g}", own[name],
                         getattr(official, function)(t, nh, rates="cen"), units="cm^-3")
        for name, function in process_map.items():
            checks.close(f"official_process_{name}_nH{nh:g}", own["losses_erg_cm3_s"][name],
                         getattr(official, function)(t, nh, rates="cen"), units="erg cm^-3 s^-1")
        checks.close(f"official_total_nH{nh:g}", own["Lvol_erg_cm3_s"],
                     official.total_cooling(t, nh, rates="cen"), units="erg cm^-3 s^-1")
        checks.close(f"H_nuclei_conserved_nH{nh:g}", own["nHI_cm3"] + own["nHII_cm3"], nh)
        checks.close(f"He_nuclei_conserved_nH{nh:g}",
                     own["nHeI_cm3"] + own["nHeII_cm3"] + own["nHeIII_cm3"], own["nHe_cm3"])
        checks.close(f"charge_conserved_nH{nh:g}", own["ne_cm3"],
                     own["nHII_cm3"] + own["nHeII_cm3"] + 2.0 * own["nHeIII_cm3"])
        checks.close(f"particles_include_neutrals_nH{nh:g}", own["n_tot_cm3"],
                     own["ne_cm3"] + own["n_ions_cm3"] + own["n_neutral_cm3"])
        rate = c.rates_cen(t)
        for label, left, right in (
            ("H", rate["C_HI"] * own["nHI_cm3"], rate["A_HII"] * own["nHII_cm3"]),
            ("HeI", rate["C_HeI"] * own["nHeI_cm3"],
             (rate["A_HeII_rad"] + rate["A_HeII_die"]) * own["nHeII_cm3"]),
            ("HeII", rate["C_HeII"] * own["nHeII_cm3"], rate["A_HeIII"] * own["nHeIII_cm3"]),
        ):
            checks.close(f"CIE_pair_balance_{label}_nH{nh:g}", left, right)
        checks.condition(f"positive_losses_and_cooling_nH{nh:g}",
                         all(np.all(values >= 0.0) for values in own["losses_erg_cm3_s"].values())
                         and np.all(np.isfinite(own["tcool_s"])) and np.all(own["tcool_s"] > 0.0))
        checks.condition(f"atomic_cools_at_least_as_fast_as_ff_nH{nh:g}",
                         np.all(own["tcool_s"] <= own["tcool_freefree_s"] * (1.0 + 1e-14)))
        relax = c.chemical_relaxation(t, nh)
        for index, temperature in enumerate(t):
            row = {"T_K": float(temperature), "nH_cm3": nh,
                   "Lvol_erg_cm3_s": float(own["Lvol_erg_cm3_s"][index]),
                   "Lambda_per_nH_squared_erg_cm3_s": float(own["Lambda_per_nH_squared_erg_cm3_s"][index]),
                   "tcool_s": float(own["tcool_s"][index]),
                   "mu": float(own["mu"][index]),
                   "tau_chem_s": float(relax["tau_chem_s"][index])}
            row.update({f"L_{name}_erg_cm3_s": float(value[index])
                        for name, value in own["losses_erg_cm3_s"].items()})
            curve_rows.append(row)
    low, high = c.atomic_cooling(t, 1e-6), c.atomic_cooling(t, 1e-2)
    factor = 1e4
    for process in c.PROCESS_NAMES:
        checks.close(f"two_body_density_scaling_{process}",
                     high["losses_erg_cm3_s"][process] / low["losses_erg_cm3_s"][process],
                     factor ** 2, rtol=DENSITY_TARGET)
    checks.close("cooling_time_density_scaling", high["tcool_s"] / low["tcool_s"],
                 1.0 / factor, rtol=DENSITY_TARGET)
    checks.close("mu_density_independence", high["mu"], low["mu"], rtol=DENSITY_TARGET)
    checks.close("chemical_time_density_scaling",
                 c.chemical_relaxation(t, 1e-2)["tau_chem_s"] /
                 c.chemical_relaxation(t, 1e-6)["tau_chem_s"], 1.0 / factor, rtol=DENSITY_TARGET)
    # Limiting states are algebraic diagnostics outside the production T range.
    neutral, ionized = c.cie_state(100.0, 1.0), c.cie_state(1e9, 1.0)
    checks.close("neutral_mu_limit_T100K", neutral["mu"], c.MU_NEUTRAL, rtol=1e-12)
    checks.close("ionized_mu_limit_T1e9K", ionized["mu"], c.MU_IONIZED, rtol=1e-6)
    checks.close("ionized_electron_limit_T1e9K", ionized["ne_cm3"],
                 (c.X_H + c.Y_HE / 2.0) / c.X_H, rtol=1e-6)
    checks.close("ionized_ion_limit_T1e9K", ionized["n_ions_cm3"],
                 (c.X_H + c.Y_HE / 4.0) / c.X_H, rtol=1e-6)
    # Independent direct eigenvalues check the stable helium quadratic times.
    for temp in np.logspace(4.0, 8.0, 9):
        rate = c.rates_cen(temp)
        a1 = rate["A_HeII_rad"] + rate["A_HeII_die"]
        a2, c0, c1 = rate["A_HeIII"], rate["C_HeI"], rate["C_HeII"]
        generator = np.array([[-c0, a1, 0.0], [c0, -a1-c1, a2], [0.0, c1, -a2]])
        eigenvalues = np.sort(np.linalg.eigvals(generator).real)
        slow_decay = -eigenvalues[1] * float(c.cie_state(temp, 1e-4)["ne_cm3"])
        checks.close(f"He_relaxation_eigenvalue_T{temp:g}",
                     c.chemical_relaxation(temp, 1e-4)["tau_He_s"], 1.0 / slow_decay,
                     rtol=1e-8, units="s")
    return curve_rows


def background_context(bg, z):
    state = bg.evaluate_a(1.0 / (1.0 + z))
    h0_s = c.H0_KM_S_MPC * 1e5 / c.MPC_CM
    return {"background_E": float(state["E"]), "background_age_Gyr": float(state["t_Gyr"]),
            "background_H_s_inv": h0_s * float(state["E"]),
            "background_Hubble_time_Gyr": 1.0 / (h0_s * float(state["E"]) * c.GYR_S)}


def attach_event_readouts(path, backgrounds):
    """Read an explicit event JSON contract, never infer event meaning.

    Primary contract: spherical_collapse.py's sphere_summary.json event_rows;
    select exactly one event='Delta200' for each background and use its
    M_solar_physical. An explicit {"events": [...]} interface with
    M200m_Msun_phys is also accepted. Additional fields are retained. This
    interface consumes a Delta=200m crossing, not a virialization assertion.
    """
    if path is None:
        return {"status": "pending_spherical_event_input", "readouts": []}
    data = json.loads(path.read_text(encoding="utf-8"))
    if "event_rows" in data:
        events = [event for event in data["event_rows"] if event.get("event") == "Delta200"]
    else:
        events = data["events"]
    if not isinstance(events, list):
        raise ValueError("events must be a list under the explicit event contract")
    outputs = []
    labels = [event["label"] for event in events]
    if sorted(labels) != sorted(backgrounds):
        raise ValueError("Require exactly one Delta200 event per background, without duplicates")
    for event in events:
        label = event["label"]
        if label not in backgrounds:
            raise ValueError(f"Unknown background label {label!r}")
        if "z_event" in event:
            z = float(event["z_event"])
        else:
            z = 1.0 / float(event["a_event"]) - 1.0
        if "a_event" in event and abs((1.0 + z) * float(event["a_event"]) - 1.0) > 1e-8:
            raise ValueError("Inconsistent event a and z")
        mass = event["M_solar_physical"] if "M_solar_physical" in event else event["M200m_Msun_phys"]
        if "Delta" in event and not np.isclose(float(event["Delta"]), 200.0, rtol=1e-8):
            raise ValueError("Cooling formation readout requires a Delta200 event")
        if "epsilon" in event and float(event["epsilon"]) != backgrounds[label].epsilon:
            raise ValueError("Event epsilon does not match its background label")
        micro = c.halo_cooling(float(mass), z)
        context = background_context(backgrounds[label], z)
        outputs.append({"label": label, "source_event": event, "microphysics": micro,
                        "background": context,
                        "interpretation": "Prescribed uniform gas at Tchar at a Delta=200m event; "
                        "not a measured gas temperature, virialization proof, or star-formation prediction."})
    return {"status": "event_readouts_computed", "input_sha256": sha256(path),
            "input_path": str(path.resolve()), "readouts": outputs}


def main():
    started = time.monotonic()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--events", type=Path)
    args = parser.parse_args()
    protocol = ROOT / "protocol.md"
    if not protocol.is_file():
        raise RuntimeError("The frozen experiment protocol must exist before execution")
    if sha256(protocol) != PROTOCOL_SHA256:
        raise RuntimeError("Experiment protocol differs from its pre-production frozen hash")
    for path in (OFFICIAL_PATH, SOURCE_METADATA, OLD_BACKGROUND):
        if not path.is_file():
            raise RuntimeError(f"Required provenance/input file missing: {path}")
    metadata = json.loads(SOURCE_METADATA.read_text(encoding="utf-8"))
    if metadata.get("commit") != COMMIT:
        raise RuntimeError("Archived Grackle metadata does not identify the pinned commit")
    for source in metadata["files"]:
        path = ROOT / "sources" / source["local_path"]
        if sha256(path) != source["sha256"]:
            raise RuntimeError(f"Archived source/license hash mismatch: {path}")
    checks = Checks()
    curves = validate_against_official(checks)
    bg_module = load_module("cooling_context_background", OLD_BACKGROUND)
    backgrounds = {"reference": bg_module.Background(epsilon=0.0),
                   "clock": bg_module.Background(epsilon=1e-4)}
    # This run re-evaluates the same background ODE, checked against every saved
    # CSV sample; it does not silently substitute a newly chosen cosmology.
    for label, bg in backgrounds.items():
        saved = np.genfromtxt(BACKGROUND_CSV[label], delimiter=",", names=True)
        reevaluated = bg.evaluate_a(saved["a"])
        for field in ("E", "t_Gyr"):
            checks.close(f"saved_background_{label}_{field}", reevaluated[field], saved[field],
                         rtol=1e-10, units="Gyr" if field == "t_Gyr" else "dimensionless")
    checks.close("background_H0_matches_microphysics", bg_module.H0_KM_S_MPC, c.H0_KM_S_MPC,
                 rtol=0.0)
    checks.close("background_Omega_m_matches_microphysics", bg_module.OMEGA_M0, c.OMEGA_M0,
                 rtol=0.0)
    grid_rows, max_pair_difference = [], 0.0
    for z in c.REDSHIFT_GRID:
        for logm in c.MASS_LOG10_GRID:
            paired = {}
            for label, bg in backgrounds.items():
                micro = c.halo_cooling(10.0 ** logm, float(z), omega_m0=bg.omega_m0)
                paired[label] = micro
                ctx = background_context(bg, z)
                grid_rows.append({"background": label, "epsilon": bg.epsilon, **micro, **ctx,
                                  "tcool_over_Hubble_time": micro["tcool_s"] * ctx["background_H_s_inv"]})
            # Same M,z and fixed physical densities imply exact microphysical identity.
            a, b = paired["reference"], paired["clock"]
            checks.condition(f"fixed_Mz_exact_null_logM{logm:g}_z{z:g}", a == b)
            for key, value in a.items():
                if isinstance(value, (float, int)) and np.isfinite(value):
                    max_pair_difference = max(max_pair_difference, abs(float(value) - float(b[key])))
            checks.condition(f"temperature_root_logM{logm:g}_z{z:g}",
                             a["Tchar_root_relative_residual"] <= ROOT_TARGET,
                             relative_residual=a["Tchar_root_relative_residual"], threshold=ROOT_TARGET)
            checks.close(f"tff_tdyn_ratio_logM{logm:g}_z{z:g}", a["tff_s"] / a["tdyn_s"],
                         np.pi / (2.0 * np.sqrt(2.0)), rtol=1e-12)
            radius_cm = a["R200m_kpc_phys"] * c.KPC_CM
            checks.close(f"mass_radius_density_logM{logm:g}_z{z:g}",
                         4.0 * np.pi / 3.0 * radius_cm ** 3 * a["rho_total_mean_g_cm3"],
                         a["M200m_Msun_phys"] * c.M_SUN_G, rtol=1e-12, units="g")
    events = attach_event_readouts(args.events, backgrounds)
    manifest = {
        "created_utc": datetime.now(timezone.utc).isoformat(), "python": platform.python_version(),
        "argv": sys.argv, "wall_seconds_before_output": time.monotonic() - started,
        "numpy": np.__version__, "scipy": scipy.__version__, "protocol_sha256": sha256(protocol),
        "code_sha256": {path.name: sha256(path) for path in (Path(__file__), Path(c.__file__), OLD_BACKGROUND)},
        "official_source_sha256": sha256(OFFICIAL_PATH), "source_metadata_sha256": sha256(SOURCE_METADATA),
        "official_commit": COMMIT,
        "background_csv_sha256": {label: sha256(path) for label, path in BACKGROUND_CSV.items()},
        "background_method": "Same Background code re-evaluated and checked against all saved CSV E/age samples.",
        "source_and_license_sha256": {source["local_path"]: source["sha256"] for source in metadata["files"]},
        "constants": {name: getattr(c, name) for name in ("G_CGS", "K_B_CGS", "M_P_G", "M_SUN_G",
                      "MPC_CM", "GYR_S", "X_H", "Y_HE", "H0_KM_S_MPC", "OMEGA_M0", "OMEGA_B0", "DELTA_M")},
        "mass_log10_Msun_phys": c.MASS_LOG10_GRID, "redshifts": c.REDSHIFT_GRID,
        "temperature_validation_K": {"log10_min": 4, "log10_max": 8, "count": 161},
        "density_validation_nH_cm3": [1e-6, 1e-4, 1e-2],
        "scope": "Historical KWH96 electron-collision CIE; X=.76,Y=.24; no UV, metals, molecules, "
                 "Compton, hydro, star formation, or clock-dependent atomic/G multipliers.",
        "temperature_scope": "10^4.5 to 10^8 K is a declared display regime, not a universal accuracy guarantee.",
        "time_scope": "tcool is instantaneous 3/2 n_tot kT / Lvol; no integrated cooling trajectory. "
                      "Chemical times are local linear CIE relaxation times at fixed T,density.",
        "density_scope": "Uniform gas at fb*200*rho_m; optical thinness is assumed, not checked from density alone.",
        "nonfinite_serialization": "Nonfinite diagnostic times serialize to null, never to a finite surrogate.",
    }
    unique_grid = [row for row in grid_rows if row["background"] == "reference"]
    domain_counts = {}
    for regime in sorted({row["temperature_regime"] for row in unique_grid}):
        selected = [row for row in unique_grid if row["temperature_regime"] == regime]
        domain_counts[regime] = {
            "n_unique_Mz_states": len(selected),
            "n_CIE_slower_than_cooling": sum(row["CIE_relaxation_slower_than_cooling"] for row in selected),
            "n_CIE_slower_than_dynamics": sum(row["CIE_relaxation_slower_than_dynamics"] for row in selected),
        }
    unique_counts = {"n_unique_Mz_states": len(unique_grid),
                     "n_CIE_slower_than_cooling": sum(row["CIE_relaxation_slower_than_cooling"] for row in unique_grid),
                     "n_CIE_slower_than_dynamics": sum(row["CIE_relaxation_slower_than_dynamics"] for row in unique_grid),
                     "by_temperature_regime": domain_counts}
    report = {"manifest": manifest, "passed": checks.passed, "n_checks": len(checks.items),
              "checks": checks.items, "fixed_Mz_microphysics_max_absolute_difference": max_pair_difference,
              "physical_warning_is_not_code_failure": True,
              "n_grid_rows_with_CIE_slower_than_cooling": sum(row["CIE_relaxation_slower_than_cooling"] for row in grid_rows),
              "n_grid_rows_with_CIE_slower_than_dynamics": sum(row["CIE_relaxation_slower_than_dynamics"] for row in grid_rows),
              "unique_grid_warning_counts": unique_counts,
              "event_status": events["status"]}
    reports = ROOT / "reports"
    reports.mkdir(parents=True, exist_ok=True)
    attempt = 1
    while (reports / f"cooling_validation_attempt_{attempt:03d}.json").exists():
        attempt += 1
    write_json(reports / f"cooling_validation_attempt_{attempt:03d}.json", report)
    if not checks.passed:
        failed = [item["name"] for item in checks.items if not item["passed"]]
        raise RuntimeError(f"Cooling validation failed; retained attempt {attempt}: {failed}")
    write_json(reports / "cooling_validation.json", report)
    write_csv(ROOT / "results" / "cooling_curve.csv", curves)
    write_csv(ROOT / "results" / "cooling_grid.csv", grid_rows)
    write_json(ROOT / "results" / "cooling_events.json", events)
    write_json(ROOT / "results" / "cooling_summary.json", {
        "manifest": manifest, "n_grid_rows": len(grid_rows), "n_curve_rows": len(curves),
        "validation_passed": checks.passed, "validation_n_checks": len(checks.items),
        "fixed_Mz_microphysics_max_absolute_difference": max_pair_difference,
        "unique_grid_warning_counts": unique_counts,
        "event_status": events["status"],
        "warnings": "tcool<tdyn is a conditional thermal diagnostic, not evidence that stars formed. "
                    "A CIE time warning reports approximation limits, not failed code.",
        "output_sha256": {name: sha256(ROOT / "results" / name) for name in
                          ("cooling_curve.csv", "cooling_grid.csv", "cooling_events.json")},
    })
    print(json.dumps({"passed": True, "checks": len(checks.items), "grid_rows": len(grid_rows),
                      "event_status": events["status"], "attempt": attempt}))


if __name__ == "__main__":
    main()
