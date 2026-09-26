#!/usr/bin/env python3
"""Independent exact algebra audit; no experiment implementation is imported."""
from pathlib import Path
import hashlib
import json

import mpmath as mp
import numpy as np
import sympy as sy

mp.mp.dps = 220
HERE = Path(__file__).resolve().parent
checks = []


def check(name, obtained, expected, tol="1e-100"):
    obtained, expected = mp.mpf(obtained), mp.mpf(expected)
    error = abs(obtained - expected) / max(abs(expected), mp.mpf("1e-170"))
    checks.append({"name": name, "passed": bool(error <= mp.mpf(tol)),
                   "relative_error": mp.nstr(error, 12), "tolerance": tol})


def flag(name, value):
    checks.append({"name": name, "passed": bool(value)})


def solve(q, rhs):
    """Independent positive tridiagonal elimination."""
    n = len(rhs)
    diagonal = [q*q + 1 for _ in range(n)]
    right = list(rhs)
    for i in range(1, n):
        factor = -q / diagonal[i-1]
        diagonal[i] -= factor * -q
        right[i] -= factor * right[i-1]
    result = [mp.mpf(0) for _ in range(n)]
    result[-1] = right[-1] / diagonal[-1]
    for i in reversed(range(n-1)):
        result[i] = (right[i] + q*result[i+1]) / diagonal[i]
    return result


def endpoint_closed(q, n):
    f = (q**(2*n) - 1) / (q**(2*n+2) - 1)
    e = q**(n-1)*(q*q - 1) / (q**(2*n+2) - 1)
    return f, e


def derivative_moments(q, n):
    eta = mp.log(q)
    f, e = endpoint_closed(q, n)
    D = n*mp.coth(n*eta) - (n+1)*mp.coth((n+1)*eta)
    Dp = -n*n/mp.sinh(n*eta)**2 + (n+1)**2/mp.sinh((n+1)*eta)**2
    E = mp.coth(eta) - (n+1)*mp.coth((n+1)*eta)
    Ep = -1/mp.sinh(eta)**2 + (n+1)**2/mp.sinh((n+1)*eta)**2
    den2 = 2*q*mp.sinh(eta)
    den3 = 8*q*q*mp.sinh(eta)**2
    return (-f*D/den2, -e*E/den2,
            f*(D*D+Dp-D*mp.coth(eta))/den3,
            e*(E*E+Ep-E*mp.coth(eta))/den3)


cases = [(str(q), n) for q in ["1.01", "1.2", "2", "3", "5"]
         for n in [1, 2, 8, 31]]
cases += [("2", 201), ("2", 202), ("3", 127), ("3", 128),
          ("5", 87), ("5", 88)]
matrix_rows = []
for q_string, n in cases:
    q = mp.mpf(q_string)
    f, e = endpoint_closed(q, n)
    tau = e/f
    rhs = [mp.mpf(0) for _ in range(n)]
    rhs[0] = 1
    g0 = solve(q, rhs)
    ge = list(reversed(g0))
    g20 = solve(q, g0)
    s2, c2, s3, c3 = derivative_moments(q, n)
    prefix = f"q={q_string},n={n}"
    check(prefix+": inverse diagonal", g0[0], f, "1e-150")
    check(prefix+": inverse endpoint", g0[-1], e, "1e-150")
    check(prefix+": tau hyperbolic", tau, mp.sinh(mp.log(q))/mp.sinh(n*mp.log(q)), "1e-150")
    check(prefix+": S2 Thomas norm", mp.fdot(g0, g0), s2, "1e-140")
    check(prefix+": C2 Thomas overlap", mp.fdot(g0, ge), c2, "1e-140")
    check(prefix+": S3 Thomas norm", mp.fdot(g0, g20), s3, "1e-140")
    check(prefix+": C3 Thomas overlap", mp.fdot(ge, g20), c3, "1e-140")
    eta_a = lambda a: mp.acosh(a/(2*q))
    f_a = lambda a: mp.sinh(n*eta_a(a))/(q*mp.sinh((n+1)*eta_a(a)))
    e_a = lambda a: mp.sinh(eta_a(a))/(q*mp.sinh((n+1)*eta_a(a)))
    check(prefix+": S2 diagonal derivative", -mp.diff(f_a, q*q+1), s2, "1e-140")
    check(prefix+": C2 diagonal derivative", -mp.diff(e_a, q*q+1), c2, "1e-140")
    check(prefix+": S3 diagonal derivative", mp.diff(f_a, q*q+1, 2)/2, s3, "1e-140")
    check(prefix+": C3 diagonal derivative", mp.diff(e_a, q*q+1, 2)/2, c3, "1e-140")
    lambdas = [q*q+1-2*q*mp.cos(mp.pi*k/(n+1)) for k in range(1, n+1)]
    weights = [mp.mpf(2)/(n+1)*mp.sin(mp.pi*k/(n+1))**2 for k in range(1, n+1)]
    for power, diag, cross in [(1, f, e), (2, s2, c2), (3, s3, c3)]:
        check(prefix+f": S{power} spectral sum", mp.fsum(w/l**power for w,l in zip(weights,lambdas)), diag)
        check(prefix+f": C{power} spectral sum", mp.fsum((-1)**k*w/l**power for k,(w,l) in enumerate(zip(weights,lambdas))), cross)
    flag(prefix+": strictly positive finite gap", lambdas[0] > (q-1)**2)
    flag(prefix+": upper bound", lambdas[-1] < (q+1)**2)
    flag(prefix+": endpoint determinant", (f*f-e*e > 0) if n > 1 else f == e)
    flag(prefix+": inverse positivity", all(x > 0 for x in g0))
    qf = float(q)
    numeric_H = np.diag(np.full(n, qf*qf+1))
    if n > 1:
        numeric_H += np.diag(np.full(n-1, -qf), 1)+np.diag(np.full(n-1, -qf), -1)
    numeric_eigs = np.linalg.eigvalsh(numeric_H)
    flag(prefix+": independent numeric eigenspectrum",
         np.allclose(numeric_eigs, np.array([float(l) for l in lambdas]), rtol=2e-12, atol=2e-14))
    momentum = mp.mpf("0.017")
    diag = q*q+1+momentum
    eta = mp.acosh(diag/(2*q))
    resolved_f = mp.sinh(n*eta)/(q*mp.sinh((n+1)*eta))
    resolved_e = mp.sinh(eta)/(q*mp.sinh((n+1)*eta))
    check(prefix+": diagonal finite momentum", mp.fsum(w/(l+momentum) for w,l in zip(weights,lambdas)), resolved_f)
    check(prefix+": cross finite momentum", mp.fsum((-1)**k*w/(l+momentum) for k,(w,l) in enumerate(zip(weights,lambdas))), resolved_e)
    matrix_rows.append({"q": q_string, "n": n, "G00": mp.nstr(f, 35),
                        "tau": mp.nstr(tau, 35), "S2": mp.nstr(s2, 35),
                        "C2": mp.nstr(c2, 35), "S3": mp.nstr(s3, 35),
                        "C3": mp.nstr(c3, 35)})

# Quartic curvature is evaluated from independent complex K derivatives.
X, XB, Z, ZB, Q, QB = sy.symbols("X XB Z ZB Q QB")
ta, ze, sg, LL = sy.symbols("tau zeta sigma Lambda", real=True)
kk = -(Z-ZB)**2/2
source = kk+sg*Q*QB
K = X*XB+kk+Q*QB-(source**2+2*ta*source*X*XB+ze*(X*XB)**2)/LL**2
metric = sy.Matrix([[sy.diff(K, a, b) for b in [XB,ZB,QB]] for a in [X,Z,Q]])
metric_fn = sy.lambdify((X,XB,Z,ZB,Q,QB,ta,ze,sg,LL), metric, "mpmath")

curvature_rows = []
for tau_string, zeta_string, sigma_string in [("0.3", "1", "1"), ("1e-61", "1", "1"),
                                               ("0.4", "0.16", "-1"), ("0.7", "2", "0.2")]:
    tau, zeta, sigma = map(mp.mpf, (tau_string, zeta_string, sigma_string))
    for phase in [mp.mpf(0), mp.pi/2, mp.mpf("0.47")]:
        for us in ["0", "0.2", "1.3"]:
            U, f, cutoff = mp.mpf(us), mp.mpf("1.7"), mp.mpf("3.1")
            FX, FZ = mp.sqrt(f)*mp.exp(1j*phase), -mp.sqrt(U)
            W = mp.matrix([-mp.conj(FX), -mp.conj(FZ), 0])
            def V(x, a, b, qr=mp.mpf(0)):
                xx = (x+1j*a)/mp.sqrt(2)
                zz = 1j*b/mp.sqrt(2)
                qq = qr/mp.sqrt(2)
                G = metric_fn(xx,mp.conj(xx),zz,mp.conj(zz),qq,qq,tau,zeta,sigma,cutoff)
                return mp.re((W.conjugate().T*(G**-1)*W)[0])
            C = FX*mp.conj(FZ)
            expected = mp.matrix([[4*zeta*f+2*tau*U, 0, -4*tau*mp.im(C)],
                                  [0,4*zeta*f+2*tau*U,4*tau*mp.re(C)],
                                  [-4*tau*mp.im(C),4*tau*mp.re(C),4*tau*f+12*U]])/cutoff**2
            label = f"curvature tau={tau_string},zeta={zeta_string},sigma={sigma_string},phase={mp.nstr(phase,4)},U={us}"
            for i in range(3):
                for j in range(i,3):
                    orders = [0,0,0]
                    orders[i] += 1; orders[j] += 1
                    actual = mp.diff(V, (0,0,0), tuple(orders))
                    error = abs(actual-expected[i,j])/max(f, U, mp.mpf(1))
                    flag(label+f": H{i}{j}", error < mp.mpf("1e-145"))
            actual_soft = mp.diff(lambda qr: V(0,0,0,qr), 0, 2)
            check(label+": charged mass", actual_soft, 2*sigma*(U+tau*f)/cutoff**2, "1e-140")
            determinant = (4*zeta*f+2*tau*U)*(4*tau*f+12*U)-16*tau*tau*f*U
            check(label+": determinant identity", determinant,
                  16*zeta*tau*f*f+(48*zeta-8*tau*tau)*f*U+24*tau*U*U, "1e-145")
            flag(label+": positive added Hessian", determinant > 0)
            curvature_rows.append({"tau": tau_string, "zeta": zeta_string, "sigma": sigma_string,
                                   "phase": mp.nstr(phase, 15), "U": us,
                                   "clock_mass2": mp.nstr(expected[2,2], 35),
                                   "hidden_mass2": mp.nstr(expected[0,0], 35)})

P = mp.mpf("2.435e27")**2
Lambda = mp.mpf("1e12")
mG = mp.mpf("1e-6")
C5 = mp.zeta(3)/(48*mp.pi**4)
A = C5*mp.mpf("1e14")**2/P
target = 4*A*Lambda**2/P
FX = mp.sqrt(3*P)*mG
hidden_mass = 2*FX/Lambda
physical = []
for qi in [2,3,5]:
    q = mp.mpf(qi)
    nstar = mp.asinh(mp.sinh(mp.log(q))/target)/mp.log(q)
    for n in [int(mp.floor(nstar)), int(mp.ceil(nstar))]:
        f, e = endpoint_closed(q,n)
        tau = e/f
        M = Lambda*mp.sqrt(f)
        gap = M*mp.sqrt(q*q+1-2*q*mp.cos(mp.pi/(n+1)))
        upper = M*mp.sqrt(q*q+1+2*q*mp.cos(mp.pi/(n+1)))
        s2,c2,s3,c3 = derivative_moments(q,n)
        row = dict(q=qi,n=n)
        quantities = {
            "n_continuous_target": nstar, "tau": tau, "rH": tau*P/(A*Lambda**2),
            "M_eV": M,"mass_gap_eV": gap,"maximum_mass_eV": upper,
            "hidden_mass_eV": hidden_mass,"hidden_mass_to_gap": hidden_mass/gap,
            "raw_FX_over_M2": FX/M**2,
            "source_weighted_FX_response": FX*mp.sqrt(s2)/M**2,
            "diagonal_momentum_fraction_at_hidden_mass": hidden_mass**2/M**2*s2/f,
            "cross_momentum_fraction_at_hidden_mass": hidden_mass**2/M**2*c2/e,
            "rank_one_hidden_mass_eV": hidden_mass*tau,
        }
        row.update({key: mp.nstr(value,40) for key,value in quantities.items()})
        physical.append(row)
        flag(f"q{qi},n{n}: restoring threshold", tau*P/(A*Lambda**2)>2)
        flag(f"q{qi},n{n}: new hidden mass below gap", hidden_mass/gap < mp.mpf("0.1"))
        flag(f"q{qi},n{n}: upper mass below 100 TeV KK", upper < mp.mpf("1e14"))
    nmax = int(mp.floor(mp.asinh(mp.sinh(mp.log(q))/(2*A*Lambda**2/P))/mp.log(q)))
    flag(f"q{qi}: upper allowed integer agrees", nmax == int(mp.ceil(nstar)))
    fn, en = endpoint_closed(q,nmax+1)
    flag(f"q{qi}: next site fails restoring", en/fn*P/(A*Lambda**2)<2)

for qi,n in [(2,202),(3,128),(5,88)]:
    q = mp.mpf(qi)
    weights = [mp.mpf(2)/(n+1)*mp.sin(mp.pi*k/(n+1))**2 for k in range(1,n+1)]
    check(f"q{qi}: endpoint mode completeness", mp.fsum(weights), 1, "1e-150")
    f_at_one = mp.mpf(n)/(n+1)
    gap_at_one = mp.sqrt(f_at_one*(2-2*mp.cos(mp.pi/(n+1))))
    flag(f"n{n}: q1 small gap counterexample", gap_at_one < mp.mpf("0.04"))

payload = {
    "date": "2026-09-26", "working_precision_digits": mp.mp.dps,
    "scope": "Independent algebra, supersymmetric zero-auxiliary spectrum, and local quartic curvature; no broken-SUSY heavy spectrum or quantum matching.",
    "implementation_imports": ["mpmath", "numpy", "sympy"],
    "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    "inputs": {"P_eV2": mp.nstr(P,40),"Lambda_eV":"1e12","mG_eV":"1e-6", "KK_eV":"1e14", "A":mp.nstr(A,40),"target_tau":mp.nstr(target,40)},
    "passed": sum(c["passed"] for c in checks), "total": len(checks),
    "failures": [c for c in checks if not c["passed"]],
    "matrix_rows": matrix_rows, "curvature_rows": curvature_rows,
    "physical_rows": physical, "checks": checks,
}
destination = HERE/"mediator_chain_algebra_20260926_checks.json"
destination.write_text(json.dumps(payload,ensure_ascii=False,indent=2)+"\n")
print(json.dumps({"passed":payload["passed"],"total":payload["total"],"failures":payload["failures"],"output":str(destination)},indent=2))
raise SystemExit(0 if payload["passed"] == payload["total"] else 1)
