#!/usr/bin/env python3
"""Independent audit of saved sphere and atomic-cooling artifacts.

No production sphere/cooling function is imported. The sphere is solved in
cosmic time with physical radius; atomic species are reconstructed from the
archived upstream reaction coefficients. This is an implementation/unit audit,
not evidence for the clock model or a validation of CIE as an astrophysical
approximation. Does not rerun the large particle-mesh experiment.
"""
from __future__ import annotations

import csv
import hashlib
import importlib.util
import json
import platform
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import scipy
from scipy.integrate import quad, solve_ivp
from scipy.interpolate import CubicSpline

ROOT = Path(__file__).resolve().parents[1]
BRIDGE = ROOT.parent / "cosmic_bridge_v01"
PROTOCOL_SHA = "f3a0cc5abf21a4374439db756cbf8c4a3103bcabfa2bee6f9bd1921b28e608f2"
G = 6.67430e-8
KB = 1.380649e-16
MP = 1.67262192369e-24
MSUN = 1.988409870698051e33
MPC = 3.0856775814913673e24
GYR = 365.25 * 86400 * 1e9
H = 0.674
H0 = 67.4e5 / MPC
OM, OB, OR = 0.315, 0.049, 9.2e-5
X, Y = 0.76, 0.24
QL, QC, AI = 2.5, 4.0, 0.02


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_csv(path):
    # CSV cells contain quoted lists and commas; parse the CSV grammar before
    # numeric conversion rather than relying on NumPy delimiter splitting.
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    columns, dtype = {}, []
    for name in rows[0]:
        values = [row[name] for row in rows]
        if all(value.lower() in ("true", "false") for value in values):
            array = np.array([value.lower() == "true" for value in values], dtype=bool)
        else:
            try:
                array = np.array(values, dtype=float)
            except ValueError:
                array = np.array(values, dtype=str)
        columns[name] = array
        dtype.append((name, array.dtype))
    result = np.empty(len(rows), dtype=dtype)
    for name, values in columns.items():
        result[name] = values
    return result


class Audit:
    def __init__(self):
        self.checks = []
        self.inputs = {}
        self.details = {}

    def track(self, path):
        self.inputs[str(path.relative_to(ROOT.parent.parent))] = sha(path)
        return path

    def check(self, name, passed, **detail):
        self.checks.append({"name": name, "passed": bool(passed), **detail})

    def close(self, name, actual, expected, rtol=1e-10, atol=0.0):
        a, b = np.broadcast_arrays(np.asarray(actual, dtype=float),
                                   np.asarray(expected, dtype=float))
        delta = abs(a - b)
        nonzero = abs(b) > 0
        relative = np.divide(delta, abs(b), out=np.zeros_like(delta), where=nonzero)
        self.check(name, np.all(np.isfinite(a)) and np.all(np.isfinite(b))
                   and np.all(delta <= atol + rtol * abs(b)),
                   max_absolute_error=float(np.max(delta)),
                   max_relative_error_nonzero_expected=float(np.max(relative)),
                   n_values=int(a.size), relative_tolerance=rtol,
                   absolute_tolerance=atol)


def audit_sphere(audit):
    summary = json.loads(audit.track(ROOT / "results/sphere_summary.json").read_text())
    bg_summary = json.loads(audit.track(BRIDGE / "results/background_summary.json").read_text())
    audit.track(ROOT / "code/spherical_collapse.py")
    audit.check("sphere_protocol_hash", summary["protocol_sha256"] == PROTOCOL_SHA)
    audit.check("sphere_code_hash", summary["code_sha256"] == sha(ROOT / "code/spherical_collapse.py"))
    # Rebuild the stated finite-start prescription from delta, f, and E, rather
    # than treating the saved y, y_N, and p as an oracle.
    delta = summary["delta_i"]
    yi = (1 + delta) ** (-1 / 3)
    fref = bg_summary["models"]["eps0"]["initial_pm"]["f_regular"]
    eref = summary["initial"]["reference"]["E_i"]
    yn_ref = -fref * delta * yi / (3 * (1 + delta))
    common_p = AI**2 * eref * yn_ref
    q = np.linspace(0, 5, 20001)
    s = np.clip((q - QL) / (QC - QL), 0, 1)
    w = 1 - 10*s**3 + 15*s**4 - 6*s**5
    dw = (-30*s**2 + 60*s**3 - 30*s**4) / (QC - QL)
    y = (1 + delta*w) ** (-1/3)
    dy = -delta*dw / (3*(1 + delta*w)**(4/3))
    jac = y*y*(y + q*dy)
    local_delta = 1/jac - 1
    audit.check("compensated_map_positive_jacobian", np.all(jac > 0),
                minimum_jacobian=float(jac.min()))
    audit.close("compensated_exterior_identity", y[q >= QC], 1, rtol=0, atol=0)
    audit.check("compensation_has_underdense_shell", np.any(local_delta < 0),
                minimum_local_delta=float(local_delta.min()))
    # The coefficient in production is rounded to nine significant figures;
    # this 3e-8 tolerance explicitly permits that known constant rounding only.
    rhom = 3*H0**2*OM/(8*np.pi*G)
    physical_mass = 4*np.pi/3 * rhom * (QL/H*MPC)**3/MSUN
    audit.close("continuous_core_mass_CGS_to_Msun_h",
                summary["inputs"]["core_mass_Msun_h"], physical_mass*H, rtol=3e-8)
    models = {}
    for label, key in (("reference", "eps0"), ("clock", "eps1e-4")):
        table = read_csv(audit.track(BRIDGE / f"results/background_{key}.csv"))
        n = np.log(table["a"])
        loge = CubicSpline(n, np.log(table["E"]), extrapolate=False)
        pressure = (np.zeros(len(table)) if label == "reference" else
                    table["Omega_chi"]*(1 + 3*table["w_chi"]))
        ps = CubicSpline(n, pressure, extrapolate=False)
        lam = bg_summary["models"][key]["omega_lambda"]
        initial = summary["initial"][label]
        ei = np.exp(loge(np.log(AI)))
        yn = common_p / (AI**2 * ei)
        audit.close(f"{label}_finite_start_y", initial["y_i"], yi, rtol=2e-14)
        audit.close(f"{label}_matched_momentum_yN", initial["y_N_i"], yn, rtol=2e-13)
        audit.close(f"{label}_matched_actual_p", initial["common_peculiar_momentum_per_RL"],
                    common_p, rtol=2e-13)
        # Independent time variable and coordinates: dN/dtau=E, ds/dtau=v,
        # dv/dtau=-Om/(2s^2)+(Lambda-Or/a^4-E^2 Omega_chi(1+3w)/2)s.
        # Here tau=H0(t-ti), s=r/RL, v=dr/(H0 RL dt).
        def rhs(tau, state):
            nn, radius, vel = state
            e = np.exp(loge(nn))
            smooth = lam - OR*np.exp(-4*nn) - 0.5*e*e*ps(nn)
            return (e, vel, -OM/(2*radius**2) + smooth*radius)

        events = []
        for threshold in (100, 200, 500):
            def cross(tau, state, threshold=threshold):
                return state[1]*np.exp(-state[0]) - threshold**(-1/3)
            cross.direction = -1
            cross.terminal = threshold == 500
            events.append(cross)
        def turn(tau, state):
            return state[2]
        turn.direction, turn.terminal = -1, False
        events.append(turn)
        state0 = (np.log(AI), AI*yi, AI*ei*(yi+yn))
        taumax = quad(lambda nn: np.exp(-loge(nn)), np.log(AI), np.log(.6),
                     epsabs=1e-12, epsrel=2e-11)[0]
        sol = solve_ivp(rhs, (0, taumax), state0, method="Radau", rtol=2e-10,
                        atol=2e-12, max_step=.001, events=events, dense_output=True)
        audit.check(f"{label}_cosmic_time_solver", sol.success)
        event_details = {}
        for j, threshold in enumerate((100, 200, 500)):
            audit.check(f"{label}_event_{threshold}_exists", len(sol.t_events[j]) == 1)
            tau = sol.t_events[j][0]
            nn, radius, vel = sol.y_events[j][0]
            target = summary["models"][label]["events"][str(threshold)]
            age = table["t_Gyr"][0] + tau/(H0*GYR)
            audit.close(f"{label}_cosmic_time_a{threshold}", np.exp(nn), target["a"], rtol=3e-7)
            audit.close(f"{label}_integrated_age{threshold}", age, target["t_Gyr"], rtol=3e-7)
            audit.close(f"{label}_physical_v{threshold}", 100*QL*vel,
                        target["physical_radial_velocity_km_s"], rtol=3e-7)
            audit.close(f"{label}_threshold{threshold}", (np.exp(nn)/radius)**3,
                        threshold, rtol=1e-10)
            event_details[str(threshold)] = {"a":float(np.exp(nn)), "t_Gyr":float(age),
                                            "velocity_km_s":float(100*QL*vel)}
        turn_n = sol.y_events[3][0, 0]
        audit.close(f"{label}_physical_turnaround_a", np.exp(turn_n),
                    summary["models"][label]["turnaround"]["a"], rtol=3e-7)
        audit.check(f"{label}_Delta200_still_infall", event_details["200"]["velocity_km_s"] < 0)
        trajectory = read_csv(audit.track(ROOT / f"results/sphere_{label}.csv"))
        audit.check(f"{label}_trajectory_has_4097_rows", len(trajectory) == 4097)
        audit.close(f"{label}_trajectory_Delta_definition", trajectory["Delta"],
                    trajectory["y"]**(-3), rtol=2e-13)
        audit.close(f"{label}_trajectory_physical_radius", trajectory["physical_radius_Mpc_h"],
                    trajectory["a"]*QL*trajectory["y"], rtol=2e-13)
        audit.close(f"{label}_trajectory_velocity_h_conversion",
                    trajectory["physical_radial_velocity_km_s"],
                    67.4*(QL/H)*trajectory["a"]*trajectory["E"]*
                    (trajectory["y"]+trajectory["y_N"]), rtol=3e-11, atol=1e-11)
        # Integrate d(H0 t)/dN directly for independent trajectory samples.
        samples = trajectory[np.arange(0, len(trajectory), 64)]
        taus = np.array([quad(lambda nn: np.exp(-loge(nn)), np.log(AI), nn,
                             epsabs=1e-12, epsrel=2e-11)[0] for nn in samples["N"]])
        states = sol.sol(taus)
        audit.close(f"{label}_physical_radius_trajectory", states[1]*np.exp(-states[0]),
                    samples["y"], rtol=3e-7)
        ni = np.linspace(np.log(.0201), np.log(.55), 2001)
        e2 = np.exp(2*loge(ni))
        acceleration_pressure = -OM*np.exp(-3*ni)/2 + lam - OR*np.exp(-4*ni) - e2*ps(ni)/2
        audit.close(f"{label}_background_pressure_acceleration", acceleration_pressure,
                    e2*(1+loge(ni, 1)), rtol=3e-7)
        models[label] = event_details
    for event in summary["event_rows"]:
        audit.close(f"event_{event['label']}_{event['event']}_mass_h_conversion",
                    event["M_solar_physical"]*H, event["M_Msun_h"], rtol=2e-14)
        audit.close(f"event_{event['label']}_{event['event']}_radius_h_conversion",
                    event["R_physical_Mpc"]*H, event["R_physical_Mpc_h"], rtol=2e-14)
    shift = models["clock"]["200"]["t_Gyr"]-models["reference"]["200"]["t_Gyr"]
    audit.close("independent_clock_Delta200_age_difference_Myr", shift*1000,
                summary["paired_event_shifts"]["200"]["delta_t_Myr"], rtol=2e-4)
    audit.details["sphere"] = {"models":models, "continuous_mass_Msun_phys_CGS":physical_mass,
                                "delta_t200_Myr":shift*1000,
                                "core_map_minimum_jacobian":float(jac.min())}


def source_state(upstream, temperature, nh):
    """CIE balance reconstructed via stable normalized reaction weights."""
    t, nh = np.broadcast_arrays(np.asarray(temperature), np.asarray(nh))
    args = {"rates":"cen"}
    hC, hA = upstream.GammaeHI(t, **args), upstream.alphaHII(t, **args)
    heC0, heC1 = upstream.GammaeHeI(t, **args), upstream.GammaeHeII(t, **args)
    heA1 = upstream.alphaHeII(t, **args)+upstream.alphad(t, **args)
    heA2 = upstream.alphaHeIII(t, **args)
    weights = np.array([heA1*heA2, heC0*heA2, heC0*heC1])
    helium = Y/(4*X)*nh*weights/weights.sum(axis=0)
    hi, hii = nh*hA/(hA+hC), nh*hC/(hA+hC)
    hei, heii, heiii = helium
    ne = hii+heii+2*heiii
    ntot = nh+Y/(4*X)*nh+ne
    return {"nHI_cm3":hi, "nHII_cm3":hii, "nHeI_cm3":hei, "nHeII_cm3":heii,
            "nHeIII_cm3":heiii, "ne_cm3":ne, "n_tot_cm3":ntot,
            "mu":nh/X/ntot, "h_decay":ne*(hA+hC),
            "He_coefficients":(heC0, heC1, heA1, heA2)}


def audit_cooling(audit):
    summary = json.loads(audit.track(ROOT / "results/cooling_summary.json").read_text())
    events = json.loads(audit.track(ROOT / "results/cooling_events.json").read_text())
    for filename in ("cooling.py", "run_cooling.py"):
        audit.track(ROOT/"code"/filename)
    meta = json.loads(audit.track(ROOT/"sources/grackle_source_metadata.json").read_text())
    for item in meta["files"]:
        path = audit.track(ROOT/"sources"/item["local_path"])
        audit.check("source_hash_"+item["local_path"], sha(path) == item["sha256"])
    spec = importlib.util.spec_from_file_location("upstream_cie_audit", ROOT/"sources/grackle_primordial_equilibrium.py")
    upstream = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(upstream)
    audit.check("cooling_frozen_protocol", summary["manifest"]["protocol_sha256"] == PROTOCOL_SHA)
    for name, expected in summary["manifest"]["code_sha256"].items():
        path = (BRIDGE/"code/background.py" if name == "background.py" else ROOT/"code"/name)
        audit.track(path)
        audit.check("cooling_code_hash_"+name, sha(path) == expected)
    curve = read_csv(audit.track(ROOT/"results/cooling_curve.csv"))
    grid = read_csv(audit.track(ROOT/"results/cooling_grid.csv"))
    audit.check("cooling_curve_grid_size", len(curve) == 483 and len(grid) == 350)
    process = {"excitation_HI":"ceHI", "excitation_HeII":"ceHeII", "ionization_HI":"ciHI",
               "ionization_HeI":"ciHeI", "ionization_HeII":"ciHeII", "recombination_HII":"rHII",
               "recombination_HeII":"rHeII", "recombination_HeIII":"rHeIII",
               "dielectronic_HeII":"drHeII", "freefree":"freefree"}
    t, nh = curve["T_K"], curve["nH_cm3"]
    for name, fn in process.items():
        audit.close("saved_curve_upstream_"+name, curve["L_"+name+"_erg_cm3_s"],
                    getattr(upstream, fn)(t, nh, rates="cen"), rtol=1e-10)
    source = source_state(upstream, t, nh)
    loss = upstream.total_cooling(t, nh, rates="cen")
    audit.close("saved_curve_upstream_total", curve["Lvol_erg_cm3_s"], loss)
    audit.close("saved_curve_energy_over_loss", curve["tcool_s"], 1.5*source["n_tot_cm3"]*KB*t/loss)
    audit.close("saved_curve_explicit_nH_squared_normalization", curve["Lambda_per_nH_squared_erg_cm3_s"], loss/nh**2)
    audit.close("saved_curve_mean_molecular_weight", curve["mu"], source["mu"])
    # All 350 saved states are audited, including below-range temperatures.
    t, nh = grid["Tchar_K"], grid["nH_cm3"]
    state = source_state(upstream, t, nh)
    for name in ("nHI_cm3", "nHII_cm3", "nHeI_cm3", "nHeII_cm3", "nHeIII_cm3", "ne_cm3", "n_tot_cm3", "mu"):
        audit.close("saved_grid_species_"+name, grid[name], state[name], rtol=1e-10)
    rhom = 3*H0**2*OM/(8*np.pi*G)*(1+grid["z"])**3
    radius = (3*grid["M200m_Msun_phys"]*MSUN/(4*np.pi*200*rhom))**(1/3)
    v2 = G*grid["M200m_Msun_phys"]*MSUN/radius
    tdyn = (3/(4*np.pi*G*200*rhom))**.5
    tff = (3*np.pi/(32*G*200*rhom))**.5
    audit.close("saved_grid_radius_physical_CGS", grid["R200m_kpc_phys"], radius/(MPC/1000))
    audit.close("saved_grid_gas_density_fb", grid["rho_g_g_cm3"], OB/OM*200*rhom)
    audit.close("saved_grid_hydrogen_nuclei_density", nh, X*OB/OM*200*rhom/MP)
    audit.close("saved_grid_temperature_self_consistent_mu", t/state["mu"], MP*v2/(2*KB), rtol=1e-8)
    audit.close("saved_grid_dynamical_time_total_density", grid["tdyn_s"], tdyn)
    audit.close("saved_grid_freefall_total_density", grid["tff_s"], tff)
    audit.close("saved_grid_tff_tdyn_ratio", grid["tff_s"]/grid["tdyn_s"], np.pi/(2*np.sqrt(2)))
    audit.close("saved_grid_thermal_energy_including_neutrals", grid["u_thermal_erg_cm3"],
                1.5*state["n_tot_cm3"]*KB*t)
    # The official utility has ill-conditioned low-temperature He abundance
    # divisions; its full cooling oracle is used only over the declared atomic
    # range. Below-range states are retained and explicitly flagged.
    inside = t >= 1e4
    atom_loss = upstream.total_cooling(t[inside], nh[inside], rates="cen")
    audit.close("saved_grid_upstream_total_atomic_range", grid["Lvol_erg_cm3_s"][inside], atom_loss)
    audit.close("saved_grid_upstream_tcool_atomic_range", grid["tcool_s"][inside],
                1.5*state["n_tot_cm3"][inside]*KB*t[inside]/atom_loss)
    # An independent charge-squared free-free identity catches missing HeIII.
    gff = 1.1+.34*np.exp(-(5.5-np.log10(t))**2/3)
    ff = 1.42e-27*gff*np.sqrt(t)*state["ne_cm3"]*(state["nHII_cm3"]+state["nHeII_cm3"]+4*state["nHeIII_cm3"])
    audit.close("saved_grid_freefree_HeIII_Z_squared", grid["L_freefree_erg_cm3_s"], ff)
    audit.close("saved_grid_local_H_relaxation", grid["tau_H_s"], 1/state["h_decay"])
    c0, c1, a1, a2 = state["He_coefficients"]
    slow = []
    for cc0, cc1, aa1, aa2, ne in zip(c0, c1, a1, a2, state["ne_cm3"]):
        generator = np.array([[-cc0, aa1, 0], [cc0, -aa1-cc1, aa2], [0, cc1, -aa2]])
        eig = np.sort(np.linalg.eigvals(generator).real)
        slow.append(1/(-eig[1]*ne))
    audit.close("saved_grid_local_He_relaxation_eigenvalues", grid["tau_He_s"], slow, rtol=1e-7)
    audit.check("saved_grid_CIE_warning_vs_cooling",
                np.all(grid["CIE_relaxation_slower_than_cooling"] == (grid["tau_chem_s"] > grid["tcool_s"])))
    audit.check("saved_grid_CIE_warning_vs_dynamics",
                np.all(grid["CIE_relaxation_slower_than_dynamics"] == (grid["tau_chem_s"] > tdyn)))
    audit.check("below_atomic_range_retained_and_flagged",
                np.all(grid["temperature_regime"][~inside] == "below_atomic_diagnostic_range"),
                n_rows_below_atomic_range=int(np.count_nonzero(~inside)))
    ref = grid[grid["background"] == "reference"]
    clock = grid[grid["background"] == "clock"]
    audit.check("same_Mz_grid_order", np.array_equal(ref["z"], clock["z"])
                and np.array_equal(ref["M200m_Msun_phys"], clock["M200m_Msun_phys"]))
    null_fields = [name for name in grid.dtype.names if name not in
                   ("background", "epsilon", "background_E", "background_age_Gyr",
                    "background_H_s_inv", "background_Hubble_time_Gyr", "tcool_over_Hubble_time")]
    audit.check("same_Mz_all_microphysics_exact_null", all(np.array_equal(ref[name], clock[name]) for name in null_fields),
                fields_checked=null_fields, n_pairs=len(ref))
    audit.check("same_Mz_age_context_differs", np.any(ref["background_age_Gyr"] != clock["background_age_Gyr"]))
    audit.check("sphere_event_cooling_completed", events["status"] == "event_readouts_computed"
                and len(events["readouts"]) == 2)
    sphere = json.loads((ROOT/"results/sphere_summary.json").read_text())
    bylabel = {row["label"]:row for row in sphere["event_rows"] if row["event"] == "Delta200"}
    for readout in events["readouts"]:
        label = readout["label"]
        actual, source = readout["microphysics"], bylabel[label]
        audit.close(f"{label}_cooling_event_physical_mass", actual["M200m_Msun_phys"], source["M_solar_physical"])
        audit.close(f"{label}_cooling_event_redshift", actual["z"], source["z_event"])
        audit.close(f"{label}_cooling_event_radius_matches_sphere",
                    actual["R200m_kpc_phys"] / 1000, source["R_physical_Mpc"], rtol=3e-8)
        audit.close(f"{label}_cooling_event_age", readout["background"]["background_age_Gyr"], source["t_Gyr"], rtol=3e-7)
    main = ref[ref["temperature_regime"] == "declared_main_temperature_range"]
    audit.details["cooling"] = {"n_curve_rows":len(curve), "n_grid_rows":len(grid),
                               "n_unique_states":len(ref), "n_unique_below_atomic_range":int(np.count_nonzero(ref["Tchar_K"] < 1e4)),
                               "n_unique_CIE_slower_than_cooling":int(np.count_nonzero(ref["CIE_relaxation_slower_than_cooling"])),
                               "n_unique_CIE_slower_than_dynamics":int(np.count_nonzero(ref["CIE_relaxation_slower_than_dynamics"])),
                               "n_main_temperature_states":len(main),
                               "n_main_CIE_slower_than_cooling":int(np.count_nonzero(main["CIE_relaxation_slower_than_cooling"])),
                               "n_main_CIE_slower_than_dynamics":int(np.count_nonzero(main["CIE_relaxation_slower_than_dynamics"]))}


def main():
    audit = Audit()
    audit.track(Path(__file__))
    audit.track(ROOT/"protocol.md")
    audit.check("protocol_bytes_frozen", sha(ROOT/"protocol.md") == PROTOCOL_SHA)
    for name, fn in (("sphere", audit_sphere), ("cooling", audit_cooling)):
        try:
            fn(audit)
        except Exception as exc:
            audit.check(name+"_audit_completed", False, exception=repr(exc))
    changed = []
    for relative, oldhash in audit.inputs.items():
        path = ROOT.parent.parent / relative
        if sha(path) != oldhash:
            changed.append(relative)
    audit.check("inputs_unchanged_during_audit", not changed, changed_paths=changed)
    result = {"created_utc":datetime.now(timezone.utc).isoformat(),
              "python":platform.python_version(), "numpy":np.__version__, "scipy":scipy.__version__,
              "scope":"Independent formulas and saved-artifact audit; no production physics module imported; no large PM rerun.",
              "protocol_sha256":PROTOCOL_SHA, "input_sha256":audit.inputs,
              "n_checks":len(audit.checks), "n_passed":sum(x["passed"] for x in audit.checks),
              "passed":all(x["passed"] for x in audit.checks), "checks":audit.checks,
              "details":audit.details,
              "limitations":["Shared input backgrounds and physical constants are assumed; this does not establish their observational validity.",
                             "Independent sphere interpolation is cubic, unlike production PCHIP; integration uses physical radius and cosmic time.",
                             "CIE oracle is the archived historical analytic Grackle utility, not modern full chemistry or an observation.",
                             "Full upstream cooling comparison excludes T<1e4 K; such states remain in saved outputs with scope flags.",
                             "No inference of virial equilibrium, gas thermalization, stars, galaxies, or resolved PM clock signal follows."]}
    path = ROOT/"results/physics_independent_validation.json"
    # Preserve prior failed results inside the owned JSON instead of silently
    # overwriting their diagnostics after corrections.
    if path.exists():
        previous = json.loads(path.read_text())
        history = previous.get("prior_failed_attempts", [])
        if not previous.get("passed", False):
            history.append({key:value for key,value in previous.items() if key != "prior_failed_attempts"})
        if history:
            result["prior_failed_attempts"] = history
    path.write_text(json.dumps(result, indent=2, ensure_ascii=False, allow_nan=False)+"\n")
    print(json.dumps({"passed":result["passed"], "checks":result["n_checks"],
                      "n_passed":result["n_passed"],
                      "failed":[x["name"] for x in result["checks"] if not x["passed"]]}))
    if not result["passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
