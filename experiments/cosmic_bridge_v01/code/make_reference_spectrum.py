#!/usr/bin/env python3
"""Generate a fixed external CAMB linear matter-spectrum input.

This is a standard transfer input to the late-time experiment, not a spectrum
derived from the Riemann model. No observational fit or nonlinear correction
is performed. A frozen protocol hash is required before CAMB is called.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import platform
import sys
import time

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
H0, OB, OC, OR_BACKGROUND = 67.4, .049, .266, 9.2e-5
H = H0 / 100
PARAMETERS = {
    "H0_km_s_Mpc": H0, "h": H, "Omega_b0": OB, "Omega_c0": OC,
    "Omega_m0": OB + OC, "ombh2": OB * H**2, "omch2": OC * H**2,
    "mnu_eV": 0., "num_massive_neutrinos": 0, "N_eff": 3.046,
    "TCMB_K": 2.7255, "As": 2.1e-9, "ns": .965, "tau_reion": .054,
    "Omega_K0": 0., "w": -1., "wa": 0., "pivot_scalar_Mpc_inv": .05,
    "nrun": 0., "r_tensor": 0., "nonlinear": False,
}
DOCS = {
    "parameter_model": "https://camb.readthedocs.io/en/latest/model.html",
    "matter_power_and_sigma8": "https://camb.readthedocs.io/en/latest/results.html",
    "density_variable_and_gauge": "https://camb.readthedocs.io/en/latest/transfer_variables.html",
}


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def row_at(redshifts, values, target):
    redshifts = np.asarray(redshifts, dtype=float)
    selected = np.flatnonzero(np.isclose(redshifts, target, rtol=0, atol=1e-8))
    if len(selected) != 1:
        raise ValueError(f"Expected exactly one explicitly calculated row at z={target}")
    return values[selected[0]]


def generate(args):
    protocol = args.protocol.resolve()
    protocol_hash = sha256(protocol)
    if protocol_hash != args.protocol_sha256:
        raise ValueError("Protocol bytes do not match the explicitly frozen SHA256.")
    if args.npoints != 2049 or args.kmin != .001 or args.kmax != 10.:
        raise ValueError("This experiment fixes 2049 output points over [.001,10] h/Mpc.")
    if args.growth_diagnostic:
        raise ValueError("The frozen v01 protocol explicitly excludes the optional k-dependent growth diagnostic.")
    try:
        import camb
        from camb import model
    except ImportError as exc:
        raise SystemExit("CAMB is unavailable. Use the separately prepared isolated CAMB environment.") from exc

    started = time.monotonic()
    # kmax supplied to CAMB is physical Mpc^-1, not h/Mpc.
    physical_kmax = 1.2 * H * args.kmax
    pars = camb.CAMBparams()
    pars.set_cosmology(H0=H0, ombh2=PARAMETERS["ombh2"], omch2=PARAMETERS["omch2"],
                       omk=0., mnu=0., num_massive_neutrinos=0, nnu=3.046,
                       TCMB=2.7255, tau=.054)
    pars.InitPower.set_params(As=2.1e-9, ns=.965, nrun=0., r=0., pivot_scalar=.05)
    pars.set_dark_energy(w=-1., wa=0.)
    pars.set_accuracy(AccuracyBoost=2., lAccuracyBoost=2., lSampleBoost=1.)
    pars.WantCls = False
    pars.DoLensing = False
    pars.WantTensors = False
    pars.set_matter_power(redshifts=[49., 0.], kmax=physical_kmax,
                          k_per_logint=50, nonlinear=False, silent=True)
    pars.NonLinear = model.NonLinear_none
    results = camb.get_results(pars)

    native_kh, native_z, native_pk = results.get_linear_matter_power_spectrum(
        var1="delta_tot", var2="delta_tot", hubble_units=True,
        k_hunit=True, have_power_spectra=True, nonlinear=False)
    physical_k, physical_z, physical_pk = results.get_linear_matter_power_spectrum(
        var1="delta_tot", var2="delta_tot", hubble_units=False,
        k_hunit=False, have_power_spectra=True, nonlinear=False)
    kh, power_z, pk = results.get_matter_power_spectrum(
        minkh=args.kmin, maxkh=args.kmax, npoints=args.npoints,
        var1="delta_tot", var2="delta_tot", have_power_spectra=True)
    p49 = np.asarray(row_at(power_z, pk, 49.))
    p0 = np.asarray(row_at(power_z, pk, 0.))
    transfer_z = np.asarray(results.Params.Transfer.PK_redshifts, dtype=float)
    sigma8 = np.asarray(results.get_sigma8(), dtype=float)
    fsigma8 = np.asarray(results.get_fsigma8(), dtype=float)
    sigma_reference = {f"z{int(z)}": {
        "sigma8": float(row_at(transfer_z, sigma8, z)),
        "fsigma8": float(row_at(transfer_z, fsigma8, z)),
        "fsigma8_over_sigma8": float(row_at(transfer_z, fsigma8, z) / row_at(transfer_z, sigma8, z)),
    } for z in (49., 0.)}

    omegas = {name: float(results.get_Omega(name, 0.))
              for name in ("baryon", "cdm", "photon", "neutrino", "nu", "de", "K")}
    omega_r = omegas["photon"] + omegas["neutrino"]
    omega_m = omegas["baryon"] + omegas["cdm"] + omegas["nu"]
    background_comparison = {}
    for z in (0., 49.):
        e2 = OR_BACKGROUND * (1+z)**4 + (OB+OC)*(1+z)**3 + 1-OR_BACKGROUND-OB-OC
        h_smooth = H0 * np.sqrt(e2)
        h_camb = float(results.hubble_parameter(z))
        background_comparison[f"z{int(z)}"] = {
            "H_camb_km_s_Mpc": h_camb, "H_smooth_reference_km_s_Mpc": float(h_smooth),
            "relative_H_difference_smooth_minus_camb": float(h_smooth/h_camb-1),
        }

    checks = []
    def check(name, passed, value=None, threshold=None):
        checks.append({"name": name, "passed": bool(passed), "value": value, "threshold": threshold})
    check("frozen_protocol_hash", protocol_hash == args.protocol_sha256, protocol_hash)
    check("nonlinear_corrections_disabled", results.Params.NonLinear == model.NonLinear_none,
          str(results.Params.NonLinear))
    check("native_k_support_covers_output", native_kh.min() <= args.kmin and native_kh.max() >= args.kmax,
          [float(native_kh.min()), float(native_kh.max())])
    check("native_unit_redshifts_identical", np.array_equal(native_z, physical_z))
    k_unit_error = float(np.max(np.abs(physical_k / (H*native_kh)-1)))
    p_unit_error = float(np.max(np.abs(physical_pk / (native_pk/H**3)-1)))
    check("k_units_physical_equals_h_times_tabulated", k_unit_error < 1e-10, k_unit_error, 1e-10)
    check("power_units_physical_equals_tabulated_over_h_cubed", p_unit_error < 1e-10, p_unit_error, 1e-10)
    check("output_shape", kh.shape == p49.shape == p0.shape == (2049,))
    check("finite_positive_linear_spectra", np.isfinite(pk).all() and np.all(pk > 0))
    check("strictly_increasing_output_k", np.all(np.diff(kh) > 0))
    check("fixed_output_endpoints", np.isclose(kh[0], .001, rtol=1e-13)
          and np.isclose(kh[-1], 10., rtol=1e-13), [float(kh[0]), float(kh[-1])])
    spacing_error = float(np.max(np.abs(np.diff(np.log(kh)) - np.log(1e4)/2048)))
    check("uniform_log_sampling", spacing_error < 1e-12, spacing_error, 1e-12)
    check("specified_matter_density", abs(omega_m-.315) < 1e-10, omega_m, 1e-10)
    check("massive_neutrino_density_zero", abs(omegas["nu"]) < 1e-15, omegas["nu"], 1e-15)
    check("flat_density_closure", abs(sum(omegas.values())-1) < 1e-10, sum(omegas.values()), 1e-10)
    check("positive_sigma8", np.isfinite(sigma8).all() and np.all(sigma8 > 0))

    scope = ("External standard CAMB linear total-matter transfer input; not generated by the Riemann model, "
             "not SMICA phase reconstruction, and not an observational fit or self-consistent scalar-era transfer.")
    numerical = {"output_kmin_h_mpc": args.kmin, "output_kmax_h_mpc": args.kmax,
                 "output_npoints": args.npoints, "camb_kmax_mpc_inv": physical_kmax,
                 "kmax_headroom_factor": 1.2, "k_per_logint": 50,
                 "AccuracyBoost": 2., "lAccuracyBoost": 2., "lSampleBoost": 1.,
                 "WantCls": False, "DoLensing": False, "NonLinear": "NonLinear_none",
                 "requested_redshifts": [49., 0.], "growth_diagnostic_k_dependent": False}
    metadata = {
        "format": "rsm-cosmic-bridge-linear-reference-v1", "scope": scope,
        "status": "PASS" if all(c["passed"] for c in checks) else "FAIL",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "protocol_sha256": protocol_hash, "generator_sha256": sha256(__file__),
        "input_note_sha256": sha256(ROOT / "REFERENCE_SPECTRUM_INPUT.md"),
        "software": {"camb_version": camb.__version__, "numpy_version": np.__version__,
                     "python_version": platform.python_version(), "platform": platform.platform(),
                     "OMP_NUM_THREADS": os.environ.get("OMP_NUM_THREADS")},
        "parameters_requested": PARAMETERS, "parameters_actual": {
            "YHe": float(results.Params.YHe), "N_eff": float(results.Params.N_eff),
            "H0": float(results.Params.H0), "ombh2": float(results.Params.ombh2),
            "omch2": float(results.Params.omch2), "omnuh2": float(results.Params.omnuh2),
            "note": "All underlying parameter values and version-dependent defaults are in camb_parameters.txt."},
        "numerical_settings": numerical,
        "density_variable": "delta_tot = density-weighted CDM+baryons; no massive neutrinos in this input",
        "gauge": "CAMB synchronous-gauge density contrast",
        "units": {"k_h_mpc": "h/Mpc", "pk_z49_mpc_over_h3": "(Mpc/h)^3",
                  "pk_z0_mpc_over_h3": "(Mpc/h)^3", "sigma8": "dimensionless, radius 8 Mpc/h",
                  "As_pivot": "physical k=0.05 Mpc^-1"},
        "sigma8_reference": sigma_reference,
        "growth_velocity_convention": "CAMB fsigma8/sigma8 is diagnostic only. Particle initial momenta use the shared smooth-radiation regular-growth background convention.",
        "background_density_comparison": {"camb_Omega0_components": omegas,
            "camb_Omega_r0": omega_r, "smooth_background_Omega_r0": OR_BACKGROUND,
            "delta_Omega_r0_camb_minus_background": omega_r-OR_BACKGROUND,
            "relative_delta_Omega_r0": omega_r/OR_BACKGROUND-1,
            "camb_Omega_m0": omega_m, "H_comparison": background_comparison,
            "note": "CAMB includes standard linear radiation/baryon transfer; the later background uses a separately prescribed smooth-radiation approximation."},
        "array_order": {"power_redshifts_returned_by_camb": np.asarray(power_z).tolist(),
                        "sigma_redshifts_returned_by_camb": transfer_z.tolist()},
        "checks": checks, "checks_passed": sum(c["passed"] for c in checks),
        "checks_total": len(checks), "official_documentation": DOCS,
        "elapsed_seconds": time.monotonic()-started,
    }
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    embedded = {"scope": scope, "protocol_sha256": protocol_hash,
                "parameters": PARAMETERS, "units": metadata["units"], "camb_version": camb.__version__}
    np.savez_compressed(output / "reference_linear_spectrum.npz", k_h_mpc=kh,
                        pk_z49_mpc_over_h3=p49, pk_z0_mpc_over_h3=p0,
                        redshifts=np.array([49., 0.]),
                        sigma8=np.array([sigma_reference["z49"]["sigma8"], sigma_reference["z0"]["sigma8"]]),
                        metadata_json=np.array(json.dumps(embedded, sort_keys=True)))
    np.savetxt(output / "reference_linear_spectrum.csv", np.column_stack((kh, p49, p0)),
               delimiter=",", header="k_h_mpc,pk_z49_mpc_over_h3,pk_z0_mpc_over_h3", comments="", fmt="%.17e")
    (output / "camb_parameters.txt").write_text(str(results.Params) + "\n")
    metadata["output_files"] = {name: {"sha256": sha256(output / name), "bytes": (output / name).stat().st_size}
                                for name in ("reference_linear_spectrum.npz", "reference_linear_spectrum.csv", "camb_parameters.txt")}
    (output / "reference_spectrum_metadata.json").write_text(json.dumps(metadata, indent=2, allow_nan=False) + "\n")
    print(json.dumps({"status": metadata["status"], "checks_passed": metadata["checks_passed"],
                      "checks_total": metadata["checks_total"], "camb_version": camb.__version__,
                      "protocol_sha256": protocol_hash, "sigma8": sigma_reference,
                      "Omega_r0_CAMB": omega_r, "relative_delta_Omega_r0": omega_r/OR_BACKGROUND-1,
                      "output_directory": str(output)}, indent=2))
    return 0 if metadata["status"] == "PASS" else 1


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--protocol", type=Path, default=ROOT / "protocol.md")
    parser.add_argument("--protocol-sha256", required=True)
    parser.add_argument("--output-dir", type=Path, default=ROOT / "inputs")
    parser.add_argument("--kmin", type=float, default=.001)
    parser.add_argument("--kmax", type=float, default=10.)
    parser.add_argument("--npoints", type=int, default=2049)
    parser.add_argument("--growth-diagnostic", action="store_true", help="Reserved option; rejected by frozen v01 protocol")
    sys.exit(generate(parser.parse_args()))
