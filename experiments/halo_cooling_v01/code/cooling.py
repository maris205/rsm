"""Primordial atomic CIE diagnostic, independently transcribed from KWH96.

Source: Katz, Weinberg & Hernquist 1996, ApJS 105, 19, Tables 1 and 2
(https://arxiv.org/pdf/astro-ph/9509107). CGS throughout unless named otherwise.
This is the historical electron-collision equilibrium subset, not the full
modern Grackle chemical network. X=.76, Y=.24; no UV, metals, molecules, dust,
Compton exchange, optical-depth correction, or arbitrary clock multiplier.

Cooling times are instantaneous thermal energy / loss rate. Halo temperature
and uniform gas density are prescribed characteristic states, not hydro results.
There is deliberately no epsilon, chi, or H(z) argument in the microphysics.
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import brentq


G_CGS = 6.67430e-8
K_B_CGS = 1.380649e-16
M_P_G = 1.67262192369e-24
M_SUN_G = 1.988409870698051e33
MPC_CM = 3.0856775814913673e24
KPC_CM = MPC_CM / 1000.0
GYR_S = 365.25 * 86400.0 * 1e9
X_H = 0.76
Y_HE = 0.24
H0_KM_S_MPC = 67.4
OMEGA_M0 = 0.315
OMEGA_B0 = 0.049
DELTA_M = 200.0
MU_NEUTRAL = 1.0 / (X_H + Y_HE / 4.0)
MU_IONIZED = 1.0 / (2.0 * X_H + 3.0 * Y_HE / 4.0)
MASS_LOG10_GRID = np.arange(8.0, 14.0 + 0.125, 0.25)
REDSHIFT_GRID = np.array([0.0, 0.5, 1.0, 2.0, 4.0, 6.0, 9.0])
PROCESS_NAMES = (
    "excitation_HI", "excitation_HeII", "ionization_HI",
    "ionization_HeI", "ionization_HeII", "recombination_HII",
    "recombination_HeII", "recombination_HeIII", "dielectronic_HeII",
    "freefree",
)


def _positive(value, name):
    value = np.asarray(value, dtype=float)
    if np.any(~np.isfinite(value)) or np.any(value <= 0.0):
        raise ValueError(f"{name} must be finite and positive")
    return value


def rates_cen(temperature_K):
    """Seven KWH96 Table 2 reaction coefficients, all in cm^3 s^-1."""
    t = _positive(temperature_K, "temperature_K")
    s = np.sqrt(t)
    d5 = 1.0 + np.sqrt(t / 1e5)
    d6 = 1.0 + (t / 1e6) ** 0.7
    common = t ** -0.5 * (t / 1e3) ** -0.2 / d6
    return {
        "C_HI": 5.85e-11 * s * np.exp(-157809.1 / t) / d5,
        "C_HeI": 2.38e-11 * s * np.exp(-285335.4 / t) / d5,
        "C_HeII": 5.68e-12 * s * np.exp(-631515.0 / t) / d5,
        "A_HII": 8.4e-11 * common,
        "A_HeII_rad": 1.5e-10 * t ** -0.6353,
        "A_HeII_die": 1.9e-3 * t ** -1.5 * np.exp(-470000.0 / t)
        * (1.0 + 0.3 * np.exp(-94000.0 / t)),
        "A_HeIII": 3.36e-10 * common,
    }


def cie_state(temperature_K, nH_cm3):
    """Nonzero-electron CIE branch; abundances are number densities in cm^-3.

    Hydrogen and helium fractions are normalized separately by their nuclei.
    Product-form helium weights remain finite when low-temperature ionization
    rates underflow. ne=0 is possible in the formal cold limit; this does not
    assert that equilibrium can be attained from exactly neutral initial gas.
    """
    t, nh = np.broadcast_arrays(_positive(temperature_K, "temperature_K"),
                               _positive(nH_cm3, "nH_cm3"))
    r = rates_cen(t)
    a1 = r["A_HeII_rad"] + r["A_HeII_die"]
    a2, c0, c1 = r["A_HeIII"], r["C_HeI"], r["C_HeII"]
    hden = r["A_HII"] + r["C_HI"]
    fh0, fh1 = r["A_HII"] / hden, r["C_HI"] / hden
    w0, w1, w2 = a1 * a2, c0 * a2, c0 * c1
    wsum = w0 + w1 + w2
    f0, f1, f2 = w0 / wsum, w1 / wsum, w2 / wsum
    nhe = Y_HE * nh / (4.0 * X_H)
    h0, h1, he0, he1, he2 = nh * fh0, nh * fh1, nhe * f0, nhe * f1, nhe * f2
    ne = h1 + he1 + 2.0 * he2
    nions = h1 + he1 + he2
    nneutral = h0 + he0
    ntot = nh + nhe + ne
    rho = nh * M_P_G / X_H
    return {
        "T_K": t, "nH_cm3": nh, "nHe_cm3": nhe, "rho_g_g_cm3": rho,
        "nHI_cm3": h0, "nHII_cm3": h1, "nHeI_cm3": he0,
        "nHeII_cm3": he1, "nHeIII_cm3": he2, "ne_cm3": ne,
        "n_ions_cm3": nions, "n_neutral_cm3": nneutral, "n_tot_cm3": ntot,
        "f_HI": fh0, "f_HII": fh1, "f_HeI": f0, "f_HeII": f1,
        "f_HeIII": f2, "mu": rho / (M_P_G * ntot),
    }


def gaunt_factor(temperature_K):
    t = _positive(temperature_K, "temperature_K")
    return 1.1 + 0.34 * np.exp(-((5.5 - np.log10(t)) ** 2) / 3.0)


def atomic_cooling(temperature_K, nH_cm3):
    """Table 1 losses in erg cm^-3 s^-1, plus explicitly normalized Lambda.

    Every named process is positive cooling. The returned coefficient uses
    n_H^2 normalization, explicitly named Lambda_per_nH_squared_erg_cm3_s.
    """
    state = cie_state(temperature_K, nH_cm3)
    t = state["T_K"]
    ne = state["ne_cm3"]
    h0, h1 = state["nHI_cm3"], state["nHII_cm3"]
    he0, he1, he2 = state["nHeI_cm3"], state["nHeII_cm3"], state["nHeIII_cm3"]
    s, d5 = np.sqrt(t), 1.0 + np.sqrt(t / 1e5)
    rc = s * (t / 1e3) ** -0.2 / (1.0 + (t / 1e6) ** 0.7)
    losses = {
        "excitation_HI": 7.50e-19 * np.exp(-118348.0 / t) / d5 * ne * h0,
        "excitation_HeII": 5.54e-17 * t ** -0.397 * np.exp(-473638.0 / t) / d5 * ne * he1,
        "ionization_HI": 1.27e-21 * s * np.exp(-157809.1 / t) / d5 * ne * h0,
        "ionization_HeI": 9.38e-22 * s * np.exp(-285335.4 / t) / d5 * ne * he0,
        "ionization_HeII": 4.95e-22 * s * np.exp(-631515.0 / t) / d5 * ne * he1,
        "recombination_HII": 8.70e-27 * rc * ne * h1,
        "recombination_HeII": 1.55e-26 * t ** 0.3647 * ne * he1,
        "recombination_HeIII": 3.48e-26 * rc * ne * he2,
        "dielectronic_HeII": 1.24e-13 * t ** -1.5 * np.exp(-470000.0 / t)
        * (1.0 + 0.3 * np.exp(-94000.0 / t)) * ne * he1,
        "freefree": 1.42e-27 * gaunt_factor(t) * s * ne * (h1 + he1 + 4.0 * he2),
    }
    total = sum(losses.values())
    energy = 1.5 * state["n_tot_cm3"] * K_B_CGS * t
    tcool = np.divide(energy, total, out=np.full_like(total, np.inf), where=total > 0.0)
    ff_time = np.divide(energy, losses["freefree"], out=np.full_like(total, np.inf),
                        where=losses["freefree"] > 0.0)
    return {
        **state, "losses_erg_cm3_s": losses, "Lvol_erg_cm3_s": total,
        "Lambda_per_nH_squared_erg_cm3_s": total / state["nH_cm3"] ** 2,
        "u_thermal_erg_cm3": energy, "tcool_s": tcool, "tcool_freefree_s": ff_time,
    }


def chemical_relaxation(temperature_K, nH_cm3):
    """Local linear CIE relaxation times at fixed temperature and density.

    At equilibrium the derivative of ne multiplying each balanced reaction
    vanishes. Hydrogen and helium linearized fraction blocks consequently use
    equilibrium ne. For helium, the slow eigenvalue of the three-state chain
    is evaluated by the stable quadratic product rather than subtraction.
    These times do not describe the finite journey from arbitrary neutral gas.
    """
    state = cie_state(temperature_K, nH_cm3)
    r = rates_cen(state["T_K"])
    ne = state["ne_cm3"]
    a1 = r["A_HeII_rad"] + r["A_HeII_die"]
    a2, c0, c1 = r["A_HeIII"], r["C_HeI"], r["C_HeII"]
    s = a1 + a2 + c0 + c1
    prod = a1 * a2 + c0 * a2 + c0 * c1
    discriminant = np.maximum(s * s - 4.0 * prod, 0.0)
    he_slow = 2.0 * prod / (s + np.sqrt(discriminant))
    h_decay = ne * (r["A_HII"] + r["C_HI"])
    he_decay = ne * he_slow
    th = np.divide(1.0, h_decay, out=np.full_like(ne, np.inf), where=h_decay > 0.0)
    the = np.divide(1.0, he_decay, out=np.full_like(ne, np.inf), where=he_decay > 0.0)
    return {"tau_H_s": th, "tau_He_s": the, "tau_chem_s": np.maximum(th, the)}


def characteristic_temperature(vc_squared_cgs):
    """Solve T/mu(T)=mp Vc^2/(2 kB); CIE mu is independent of density."""
    v2 = float(_positive(vc_squared_cgs, "vc_squared_cgs"))
    target = M_P_G * v2 / (2.0 * K_B_CGS)
    # Endpoint bounds follow directly from neutral and fully ionized mu.
    low = np.log(target * MU_IONIZED * (1.0 - 1e-9))
    high = np.log(target * MU_NEUTRAL * (1.0 + 1e-9))

    def residual(logt):
        temperature = np.exp(logt)
        mu = float(cie_state(temperature, 1.0)["mu"])
        return np.log(temperature / mu) - np.log(target)

    temperature = float(np.exp(brentq(residual, low, high, xtol=1e-13, rtol=1e-13)))
    mu = float(cie_state(temperature, 1.0)["mu"])
    return temperature, float(abs(temperature / mu / target - 1.0))


def halo_cooling(mass_Msun_phys, z, h0_km_s_mpc=H0_KM_S_MPC,
                 omega_m0=OMEGA_M0, omega_b0=OMEGA_B0):
    """One prescribed M200m state; returned numerical fields are scalar CGS.

    The result has no dependency on a clock parameter. Different backgrounds
    with the same H0, Omega_m0 and Omega_b0 give identical results at fixed M,z.
    A later caller may attach background H(z), age, or distinct formation z.
    """
    mass = float(_positive(mass_Msun_phys, "mass_Msun_phys"))
    z = float(z)
    if not np.isfinite(z) or z < 0.0:
        raise ValueError("This experiment requires finite z >= 0")
    if not (0.0 < omega_b0 < omega_m0 and h0_km_s_mpc > 0.0):
        raise ValueError("Require 0 < Omega_b0 < Omega_m0 and H0 > 0")
    h0 = h0_km_s_mpc * 1e5 / MPC_CM
    rhom = 3.0 * h0 * h0 * omega_m0 / (8.0 * np.pi * G_CGS) * (1.0 + z) ** 3
    rho_total = DELTA_M * rhom
    rho_g = omega_b0 / omega_m0 * rho_total
    mass_g = mass * M_SUN_G
    radius = (3.0 * mass_g / (4.0 * np.pi * rho_total)) ** (1.0 / 3.0)
    v2 = G_CGS * mass_g / radius
    tchar, root_residual = characteristic_temperature(v2)
    nh = X_H * rho_g / M_P_G
    atom = atomic_cooling(tchar, nh)
    relax = chemical_relaxation(tchar, nh)
    tdyn = radius / np.sqrt(v2)
    tff = np.sqrt(3.0 * np.pi / (32.0 * G_CGS * rho_total))
    tcool, tauchem = float(atom["tcool_s"]), float(relax["tau_chem_s"])
    if tchar < 1e4:
        temperature_regime = "below_atomic_diagnostic_range"
    elif tchar < 10.0 ** 4.5:
        temperature_regime = "low_temperature_sensitive"
    elif tchar <= 1e8:
        temperature_regime = "declared_main_temperature_range"
    else:
        temperature_regime = "above_declared_main_temperature_range"
    warnings = []
    if temperature_regime != "declared_main_temperature_range":
        warnings.append(temperature_regime)
    if tauchem > tcool:
        warnings.append("local_CIE_relaxation_slower_than_instantaneous_cooling")
    if tauchem > tdyn:
        warnings.append("local_CIE_relaxation_slower_than_dynamical_time")
    out = {key: float(value) for key, value in atom.items() if key != "losses_erg_cm3_s"}
    out.update({f"L_{key}_erg_cm3_s": float(value)
                for key, value in atom["losses_erg_cm3_s"].items()})
    out.update({key: float(value) for key, value in relax.items()})
    out.update({
        "M200m_Msun_phys": mass, "z": z, "a": 1.0 / (1.0 + z),
        "rho_m_mean_g_cm3": rhom, "rho_total_mean_g_cm3": rho_total,
        "R200m_kpc_phys": radius / KPC_CM, "V200_km_s": np.sqrt(v2) / 1e5,
        "Tchar_K": tchar, "Tchar_root_relative_residual": root_residual,
        "tdyn_s": tdyn, "tff_s": tff, "tcool_Gyr": tcool / GYR_S,
        "tdyn_Gyr": tdyn / GYR_S, "tff_Gyr": tff / GYR_S,
        "tcool_over_tdyn": tcool / tdyn, "tcool_over_tff": tcool / tff,
        "tauchem_over_tcool": tauchem / tcool if np.isfinite(tcool) else None,
        "tauchem_over_tdyn": tauchem / tdyn,
        "temperature_regime": temperature_regime,
        "CIE_relaxation_slower_than_cooling": bool(tauchem > tcool),
        "CIE_relaxation_slower_than_dynamics": bool(tauchem > tdyn),
        "warnings": warnings,
    })
    return out
