"""Independent local algebra and scale checks; does not run cosmological ODEs."""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

import sympy as sp


def main() -> None:
    chi, mu, m_inf, xi, ref, alpha0 = sp.symbols(
        "chi mu m_inf xi ref alpha0", positive=True
    )
    pi = sp.pi
    s = m_inf**2 * (1 + xi / chi**2)
    ds = sp.diff(s, chi)
    dds = sp.diff(ds, chi)
    potential = s**2 * (sp.log(s / mu**2) - sp.Rational(3, 2)) / (32 * pi**2)
    checks: list[dict] = []

    def equal(name: str, actual: sp.Expr, expected: sp.Expr) -> None:
        difference = sp.simplify(actual - expected)
        checks.append({"name": name, "passed": difference == 0, "residual": str(difference)})

    equal("CW first derivative", sp.diff(potential, chi), s * ds * (sp.log(s / mu**2) - 1) / (16 * pi**2))
    equal("CW second derivative", sp.diff(potential, chi, 2), (ds**2 * sp.log(s / mu**2) + s * dds * (sp.log(s / mu**2) - 1)) / (16 * pi**2))
    x = sp.symbols("x")
    equal("bubble Feynman parameter integral", sp.integrate(x * (1 - x), (x, 0, 1)), sp.Rational(1, 6))
    equal("complex scalar kinetic normalization", 2 * sp.Rational(1, 2) * sp.Rational(1, 6) * ds**2 / (16 * pi**2 * s), ds**2 / (96 * pi**2 * s))
    e = sp.symbols("e", positive=True)
    equal("beta e to beta alpha conversion", sp.diff(e**2 / (4 * pi), e) * e**3 / (48 * pi**2), (e**2 / (4 * pi))**2 / (6 * pi))
    b = 1 - alpha0 / (12 * pi) * sp.log(s / s.subs(chi, ref))
    a = alpha0 / b
    equal("threshold logarithmic derivative", sp.diff(a, chi) / a, a * ds / (12 * pi * s))
    g, gp, v, acc, vp, h, source = sp.symbols("g gp v acc vp h source")
    rho_dot = sp.Rational(1, 2) * gp * v**3 + g * v * acc + vp * v
    equal("scalar continuity and source sign", (rho_dot + 3 * h * g * v**2).subs(acc, -3 * h * v - gp * v**2 / (2 * g) - (vp + source) / g), -source * v)
    r, r_ref, f, velocity = sp.symbols("r r_ref f velocity", positive=True)
    equal("exponential metric preserves R kinetic energy", f**2 * (r / r_ref)**2 * (velocity / r)**2 / 2, f**2 * velocity**2 / (2 * r_ref**2))
    gv = sp.Function("G")(chi)
    vv = sp.Function("V")(chi)
    equal("canonical potential curvature", sp.diff(sp.diff(vv, chi) / sp.sqrt(gv), chi) / sp.sqrt(gv), sp.diff(vv, chi, 2) / gv - sp.diff(gv, chi) * sp.diff(vv, chi) / (2 * gv**2))
    s0 = sp.symbols("s0", positive=True)
    u_s = s0**2 * (sp.log(s0 / mu**2) - sp.Rational(3, 2)) / (32 * pi**2)
    equal("CW reference slope with respect to log mass", (2 * s0 * sp.diff(u_s, s0)).subs(mu**2, s0), -s0**2 / (8 * pi**2))

    alpha = 1 / 137.035999084
    mp = 2.435e27
    hbar_evs = 6.582119569e-16
    mpc_km = 3.085677581491367e19
    h_ref = 67.4 / mpc_km * hbar_evs
    rho_ref = 3 * mp**2 * h_ref**2
    delta_alpha = 1e-6
    scales = []
    for mass_ev in (1.0, 1e6, 1e11):
        delta_v = 3 / (4 * math.pi * alpha) * mass_ev**4 * delta_alpha
        scales.append({"mass_eV": mass_ev, "delta_U1_eV4": delta_v, "delta_U1_over_reference_critical_density": delta_v / rho_ref})
    c0, xi0, mass0 = 5.0, 1.0, 1e11
    mi2 = mass0**2 / (1 + xi0 / c0**2)
    derivative = -2 * mi2 * xi0 / c0**3
    second = 6 * mi2 * xi0 / c0**4
    vp0 = -mass0**2 * derivative / (16 * math.pi**2)
    vpp0 = -mass0**2 * second / (16 * math.pi**2)
    result = {
        "scope": "Independent algebra and dimensional-scale checks, not a cosmological or observational validation",
        "passed": all(c["passed"] for c in checks),
        "check_count": len(checks),
        "checks": checks,
        "references": [{"url": "https://arxiv.org/pdf/1412.1837v2", "locators": ["Eq. (2.54), PDF p. 28", "Eq. (2.55) and following coefficient list, PDF p. 30"]}],
        "reference_scales": {"H_eV": h_ref, "critical_density_eV4": rho_ref, "fractional_alpha_change_assumed": delta_alpha, "delta_log_mass_linear": 6 * math.pi * delta_alpha / alpha},
        "unsubtracted_loop_scale_examples": scales,
        "conditional_chi5_xi1_100GeV_FMpl_example": {
            "potential_first_derivative_eV4": vp0,
            "potential_second_derivative_eV4": vpp0,
            "loop_metric_over_tree_metric": derivative**2 / (96 * math.pi**2 * mass0**2 * mp**2),
            "sqrt_abs_potential_curvature_over_H": math.sqrt(abs(vpp0) / mp**2) / h_ref,
            "negative_curvature_not_stable_oscillator_mass": vpp0 < 0,
        },
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    out = Path(__file__).with_suffix(".json")
    out.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({"passed": result["passed"], "check_count": len(checks), "output": str(out)}, ensure_ascii=False))
    if not result["passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
