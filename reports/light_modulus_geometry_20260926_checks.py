#!/usr/bin/env python3
"""Independent retained-Higgs-modulus quotient and current checks.

No main experiment implementation is imported.  Exact rational linear algebra
checks the quotient projection; high precision component differentiation checks
the scalar metric.  This is the constrained, rigid, two-derivative tree theory,
not a finite-auxiliary-field relaxation or a supergravity computation.
"""
from pathlib import Path
import hashlib
import json
import math
import platform
import sys

import mpmath as mp
import sympy as sp

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
CHECKS = []


def check(name, passed, detail=None):
    row = {"name": name, "passed": bool(passed)}
    if detail is not None:
        row["detail"] = detail
    CHECKS.append(row)


def near(a, b, tol=mp.mpf("1e-70")):
    return abs(a-b) <= tol*max(mp.mpf(1), abs(a), abs(b))


def charge(q, n):
    Q = sp.zeros(n+1, n)
    Q[0, 0] = q
    for a in range(1, n):
        Q[a, a-1], Q[a, a] = 1, -q
    Q[n, n-1] = 1
    return Q


for q in (2, 3, 5):
    for n in (1, 2, 4, 6):
        tag = f"exact_q{q}_n{n}"
        Q = charge(q, n)
        w = sp.Matrix([1]+[-q**a for a in range(1, n+1)])
        W = (w.T*w)[0]
        L = w*w.T/W
        P = Q*(Q.T*Q).inv()*Q.T
        check(tag+"_null", Q.T*w == sp.zeros(n, 1))
        check(tag+"_primitive", math.gcd(*map(int, w)) == 1)
        check(tag+"_one_complex_modulus", 2*(n+1)-(n+1)-Q.rank() == 1)
        check(tag+"_projector", P == sp.eye(n+1)-L)
        check(tag+"_heavy_light_completeness", P+L == sp.eye(n+1))
        check(tag+"_opposite_endpoint_residues", P[0, n] == -L[0, n])
        check(tag+"_right_weights_ratio", sp.cancel(L[n,n]/P[n,n]-(q*q-1)/(1-sp.Rational(1,q**(2*n)))) == 0)
        # General positive weights, independently of the explicit Higgs solution.
        S = sp.diag(*[sp.Rational(2*a+3, a+2) for a in range(n+1)])
        weighted = S*Q*(Q.T*S*Q).inv()*Q.T*S
        h_w = (w.T*S.inv()*w)[0]
        check(tag+"_weighted_source_kernel", weighted == S-w*w.T/h_w)
        check(tag+"_source_kernel_has_one_null", weighted*S.inv()*w == sp.zeros(n+1,1))
        # Residual global charge and local gauge/gravitational anomalies.
        check(tag+"_global_gravity_anomaly", sum(w)+sum(-w) == 0)
        check(tag+"_global_gauge_squared_anomaly", all(
            sum(w[a]*Q[a,i]*Q[a,j]+(-w[a])*(-Q[a,i])*(-Q[a,j])
                for a in range(n+1)) == 0
            for i in range(n) for j in range(n)))
        check(tag+"_global_cubed_anomaly", sum(x**3+(-x)**3 for x in w) == 0)


mp.mp.dps = 100
numeric_examples = []
for q, n, v_text in ((2,1,"0.7"), (3,2,"1.0"), (3,4,"2.3"), (5,8,"1.3")):
    v = mp.mpf(v_text)
    w = mp.matrix([1]+[-mp.mpf(q)**a for a in range(1,n+1)])
    W = sum(x*x for x in w)
    u = w/mp.sqrt(W)
    Q = mp.matrix(charge(q,n).tolist())

    def radii(zz):
        return [mp.sqrt(4*v**4+ua**2*zz**2) for ua in u]

    def real_logs(zz):
        return [mp.asinh(ua*zz/(2*v*v))/2 for ua in u]

    def real_invariant(zz):
        return sum(ua*rr for ua,rr in zip(u,real_logs(zz)))

    def quotient_K(zz):
        return sum(radii(zz))

    for scale_text in ("-100", "-1", "0", "0.1", "1", "10", "100"):
        z = mp.mpf(scale_text)*v*v
        tag = f"metric_q{q}_n{n}_v{v_text}_zoverv2_{scale_text}"
        s = radii(z)
        h = sum(ua*ua/sa for ua,sa in zip(u,s))
        r = real_logs(z)
        # Direct fields, before using the proposed d/s kinetic formulas.
        A = [v*v*mp.exp(2*rr) for rr in r]
        B = [v*v*mp.exp(-2*rr) for rr in r]
        d = mp.matrix([aa-bb for aa,bb in zip(A,B)])
        check(tag+"_products", all(near(aa*bb,v**4) for aa,bb in zip(A,B)))
        check(tag+"_Dflat", all(near(x,0) for x in Q.T*d))
        check(tag+"_difference_coordinate", all(near(d[a],u[a]*z) for a in range(n+1)))
        direct_radial = sum(
            mp.diff(lambda zz: v*mp.exp(real_logs(zz)[a]), z)**2
            +mp.diff(lambda zz: v*mp.exp(-real_logs(zz)[a]), z)**2
            for a in range(n+1))
        check(tag+"_component_radial_metric", near(direct_radial,h/4))
        check(tag+"_holomorphic_R_derivative", near(mp.diff(real_invariant,z),h/2))
        # Differentiate K(R(z)) twice by the chain rule, not its claimed metric.
        R1 = mp.diff(real_invariant,z)
        R2 = mp.diff(real_invariant,z,2)
        K1 = mp.diff(quotient_K,z)
        K2 = mp.diff(quotient_K,z,2)
        KRR = (K2*R1-K1*R2)/(R1**3)
        G = KRR/(8*v*v)
        check(tag+"_Kahler_metric", near(G,1/(2*v*v*h)))
        check(tag+"_positive_metric", h>0 and G>0)
        # Eliminate gauge velocities for a generic initial phase velocity.
        S = mp.diag(s)
        phidot = mp.matrix([mp.mpf(a+1)/(a+2) for a in range(n+1)])
        gauge_velocity = mp.lu_solve(Q.T*S*Q,Q.T*S*phidot)
        horizontal = phidot-Q*gauge_velocity
        theta_dot = (u.T*phidot)[0]
        angular = (horizontal.T*S*horizontal)[0]
        check(tag+"_gauge_projected_angular_metric", near(angular,theta_dot**2/h))
        check(tag+"_horizontal_velocity", all(near(horizontal[a],u[a]*theta_dot/(s[a]*h)) for a in range(n+1)))
        check(tag+"_complex_radial_agreement", near(G*2*v*v*R1**2,h/4))
        check(tag+"_complex_angular_agreement", near(G*2*v*v,1/h))
        numeric_examples.append({"q":q,"n":n,"v":v_text,"z_over_v2":scale_text,
                                 "h_v2":mp.nstr(h*v*v,35),"K_CbarC":mp.nstr(G,35)})
    h0 = sum(ua*ua/(2*v*v) for ua in u)
    check(f"canonical_q{q}_n{n}_radial", near(h0/4,1/(8*v*v)))
    check(f"canonical_q{q}_n{n}_angle", near(1/h0,2*v*v))
    check(f"canonical_q{q}_n{n}_complex", near(1/(2*v*v*h0),1))
    # Norm-only source terms and product superpotential are exactly phase blind.
    angles = [mp.mpf("0.17")*(a+1) for a in range(n+1)]
    alpha = mp.mpf("0.43")
    for a in range(n+1):
        phi = v*mp.exp(real_logs(v*v)[a]+1j*(angles[a]+u[a]*alpha))
        tilde = v*mp.exp(-real_logs(v*v)[a]-1j*(angles[a]+u[a]*alpha))
        check(f"phase_q{q}_n{n}_a{a}_product_invariant", near(phi*tilde,v*v))
        check(f"phase_q{q}_n{n}_a{a}_norm_invariant", near(abs(phi)**2-abs(tilde)**2,u[a]*v*v))


mp.mp.dps = 180
q, n = 3, 127
W_int = sum(q**(2*a) for a in range(n+1))
W = mp.mpf(W_int)
u0, un = 1/mp.sqrt(W), -mp.mpf(q)**n/mp.sqrt(W)
tau_site = (q-1/mp.mpf(q))/(mp.mpf(q)**n-mp.mpf(q)**(-n))
tau_local_equal_c = -u0*un/(1-u0*u0)
check("benchmark_nonzero_u0", u0>0 and u0<mp.mpf("1e-60"))
check("benchmark_site_normalization", near(q*tau_local_equal_c,tau_site,mp.mpf("1e-165")))
check("benchmark_right_un_suppressed", abs(un)>mp.mpf("0.94"))
check("benchmark_right_light_heavy_ratio", near(un*un/(1-un*un),mp.mpf(8),mp.mpf("1e-120")))
check("benchmark_equal_source_light_ratio", near(u0/un,-mp.mpf(q)**(-n),mp.mpf("1e-165")))

benchmark = {
    "q":q,"n":n,"W_exact":str(W_int),
    "u0":mp.nstr(u0,65),"un":mp.nstr(un,65),
    "left_over_right_light_coupling_equal_c":mp.nstr(u0/un,65),
    "left_over_right_light_coupling_cX_eq_qcI":mp.nstr(u0/(q*un),65),
    "tau_local_equal_c":mp.nstr(tau_local_equal_c,65),
    "tau_old_site_cX_eq_qcI":mp.nstr(tau_site,65),
    "right_light_to_heavy_squared_weight":mp.nstr(un*un/(1-un*un),65),
    "left_light_squared_weight":mp.nstr(u0*u0,65),
    "right_light_squared_weight":mp.nstr(un*un,65),
    "phase_f_over_v":mp.nstr(2*u0,65),
    "phase_period_note":"For primitive w and integer gauge charges, theta period is 2pi/sqrt(W); a=2v theta at the origin. This alone is not a scattering cutoff."
}

inputs = ["higgs_phase_operators_20260926.md", "higgs_local_currents_20260926.md"]
result = {
    "stage":"light_modulus_v01_independent_geometry",
    "date":"2026-09-26", "baseline_commit":"fbfae9851ce0c38452d18d2ec8d70b1d840e6c6d",
    "scope":"Rigid classical canonical K, exact product constraint and D-flat quotient; leading two-derivative/auxiliary expansion. No finite-F heavy relaxation, supergravity, loop matching, cosmological time, or observational inference.",
    "implementation_independent":True,
    "environment":{"python":platform.python_version(),"mpmath":mp.__version__,"sympy":sp.__version__},
    "source_hashes":{name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in inputs},
    "script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    "checks_total":len(CHECKS),"checks_passed":sum(x["passed"] for x in CHECKS),
    "checks_failed":sum(not x["passed"] for x in CHECKS),
    "checks":CHECKS,"benchmark":benchmark,"numeric_examples":numeric_examples,
    "initial_failures":[],
    "verification_limit":"These are internal algebra/implementation checks, not empirical support or external peer review."
}
out = HERE/"light_modulus_geometry_20260926_checks.json"
out.write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n")
print(json.dumps({k:result[k] for k in ("checks_total","checks_passed","checks_failed")},ensure_ascii=False))
print(json.dumps(benchmark,indent=2,ensure_ascii=False))
sys.exit(0 if result["checks_failed"] == 0 else 1)
