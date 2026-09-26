#!/usr/bin/env python3
"""Independent finite-F, original-charge-chain stationary audit.

No import of experiment implementation; exact matrices and direct nonlinear
Higgs-amplitude KKT solves. Output is internal consistency, not observations.
"""
from __future__ import annotations

import hashlib
import json
import platform
from pathlib import Path

import numpy as np
import scipy
from scipy.optimize import brentq, root
import sympy as sp


HERE = Path(__file__).resolve().parent
OUT = HERE / "light_modulus_stationarity_20260926_checks.json"
checks = []


def check(name, condition, **details):
    checks.append({"name": name, "passed": bool(condition), **details})


def charge(n, q, exact=False):
    Q = sp.zeros(n + 1, n) if exact else np.zeros((n + 1, n))
    Q[0, 0] = q
    for a in range(1, n):
        Q[a, a - 1], Q[a, a] = 1, -q
    Q[n, n - 1] = 1
    return Q


def exact_checks():
    for n in range(1, 6):
        for q in (2, 3, 5):
            Q = charge(n, q, True)
            w = sp.Matrix([1] + [-q**a for a in range(1, n + 1)])
            W = (w.T * w)[0]
            P = sp.eye(n + 1) - w * w.T / W
            H = Q.T * Q
            A = Q * Q.T
            Ap = Q * H.inv()**2 * Q.T
            tag = f"n{n}_q{q}"
            check(tag + "_left_null", Q.T * w == sp.zeros(n, 1))
            check(tag + "_rank", Q.rank() == n)
            check(tag + "_projector", Q * H.inv() * Q.T == P)
            check(tag + "_pseudoinverse_left", A * Ap == P)
            check(tag + "_pseudoinverse_right", Ap * A == P)
            check(tag + "_pseudoinverse_kernel", Ap * w == sp.zeros(n + 1, 1))
            check(tag + "_pseudoinverse_end", Ap[n, n] == (H.inv()**2)[n - 1, n - 1])
            # The force is -S cX e_n / D^2 + g^2 A r.
            # Dotting into w eliminates D force exactly for arbitrary r.
            rr = sp.Matrix(sp.symbols("r0:" + str(n + 1)))
            check(tag + "_D_force_null", (w.T * A * rr)[0] == 0)

    r, p, hh = sp.symbols("r p hh", real=True, positive=True)
    # Inverse real coordinate map for |hh|<1.
    disc = sp.sqrt(r*r + 4*p*p*(1-hh*hh))
    d = (r-hh*disc)/(1-hh*hh)
    s = (disc-hh*r)/(1-hh*hh)
    check("map_inverse_r", sp.simplify(d+hh*s-r) == 0)
    check("map_inverse_product", sp.simplify(s*s-d*d-4*p*p) == 0)

    S, c, un, B, t, D, g, h = sp.symbols("S c un B t D g h", nonzero=True, real=True)
    # At fixed z: t D^2=S c^2 h/g^2 and D=B+t.
    Dprime = c*un*D/(D+2*t)
    tprime = Dprime-c*un
    Vprime = -S*Dprime/D**2 + g*g*t*tprime/(c*c*h)
    Vprime = sp.simplify(Vprime.subs(g*g/(c*c*h), S/(t*D**2)))
    check("envelope_derivative", sp.simplify(Vprime+S*c*un/D**2) == 0)
    Vsecond = sp.diff(-S*c*un/D**2, D)*Dprime
    check("envelope_convexity", sp.simplify(Vsecond-2*S*(c*un)**2/(D**2*(D+2*t))) == 0)


def numeric_branch(n, q, S, cX, g, z):
    Q = charge(n, q)
    w = np.array([1.] + [-float(q)**a for a in range(1, n+1)])
    u = w / np.linalg.norm(w)
    A = Q @ Q.T
    Hinv = np.linalg.inv(Q.T @ Q)
    Ap = Q @ Hinv @ Hinv @ Q.T
    h = Ap[-1, -1]
    B = 1 + cX*u[-1]*z
    C = S*cX*cX*h/(g*g)
    # Solve for D rather than t; D^2(D-B)=C, D>max(B,0).
    low = max(B, 0.)
    high = max(low+1., 1.)
    while high*high*(high-B) < C:
        high *= 2
    D = brentq(lambda dd: dd*dd*(dd-B)-C, low, high,
               xtol=5e-15, rtol=5e-15)
    t = C/(D*D)
    alpha = S*cX/(g*g*D*D)
    rr = z*u + alpha*Ap[:, -1]
    V = S/D + g*g*t*t/(2*cX*cX*h)
    return Q, u, A, Ap, h, D, t, rr, V


def amplitude_kkt(vec, *, n, Q, u, S, cX, g, lam, cI, k, z, v):
    m = n+1
    xx, yy, mult = vec[:m], vec[m:2*m], vec[-1]
    AA, BB = np.exp(2*xx), np.exp(2*yy)
    pp = np.exp(xx+yy)
    dd, ss = AA-BB, AA+BB
    rr = dd.copy()
    rr[0] += cI*k*ss[0]
    D = 1+cX*dd[-1]
    ar = Q @ (Q.T @ rr)
    leftplus, leftminus = np.ones(m), np.ones(m)
    leftplus[0] += cI*k
    leftminus[0] -= cI*k
    px = 2*AA*leftplus
    py = -2*BB*leftminus
    radial = 2*lam*lam*(pp-v*v)*pp
    gradx = radial + g*g*ar*px
    grady = radial + g*g*ar*py
    gradx[-1] -= 2*S*cX*AA[-1]/D**2
    grady[-1] += 2*S*cX*BB[-1]/D**2
    Cx, Cy = u*px, u*py
    residual = np.r_[gradx+mult*Cx, grady+mult*Cy, u@rr-z]
    potential = S/D + lam*lam*np.sum((pp-v*v)**2) + g*g/2*np.dot(Q.T@rr, Q.T@rr)
    return residual, potential, rr, dd, ss, pp, D


def nonlinear_checks():
    rows = []
    for n in (1, 2, 4):
        for q in (2, 3):
            for S in (0.01, 0.5, 4.):
                for cX in (-0.3, 0.3):
                    for z in (-2., 0., 2.):
                        for k in (0., 0.4):
                            g, lam, cI, v = 0.7, 2., 0.2, 1.
                            Q, u, A, Ap, h, D, t, rr, V = numeric_branch(n, q, S, cX, g, z)
                            kw = dict(n=n,Q=Q,u=u,S=S,cX=cX,g=g,lam=lam,cI=cI,k=k,z=z,v=v)
                            # Same fixed zero log-amplitude initial point for all cases.
                            sol = root(lambda vv: amplitude_kkt(vv, **kw)[0],
                                       np.zeros(2*(n+1)+1), method="hybr",
                                       options={"xtol": 1e-10, "maxfev": 4000})
                            residual, Vnum, rnum, dnum, snum, pnum, Dnum = amplitude_kkt(sol.x, **kw)
                            res = float(np.max(np.abs(residual)))
                            rerr = float(np.max(np.abs(rnum-rr)))
                            verr = float(abs(Vnum-V))
                            # Full clock/Higgs metric positivity via Schur complement.
                            # K_TPhi = cI k_T Phi* ; |k_T|^2=2k.
                            schur = 1+cI*dnum[0] - 2*cI*cI*k*(
                                (snum[0]+dnum[0])/2/(1+cI*k)
                                +(snum[0]-dnum[0])/2/(1-cI*k))
                            deriv = -S*cX*u[-1]/D**2
                            width=2e-4
                            Vp=numeric_branch(n,q,S,cX,g,z+width)[-1]
                            Vm=numeric_branch(n,q,S,cX,g,z-width)[-1]
                            Vpp=numeric_branch(n,q,S,cX,g,z+2*width)[-1]
                            Vmm=numeric_branch(n,q,S,cX,g,z-2*width)[-1]
                            # Fourth-order centered derivative, same step and tolerance.
                            fd=(-Vpp+8*Vp-8*Vm+Vmm)/(12*width)
                            derr=float(abs(fd-deriv))
                            row = dict(n=n,q=q,S=S,cX=cX,z=z,k=k,solver_success=bool(sol.success),
                                       solver_message=str(sol.message),residual=res,r_error=rerr,
                                       energy_error=verr,derivative_fd_error=derr,
                                       hidden_metric=float(Dnum),clock_schur=float(schur),
                                       derivative=float(deriv),potential=float(V),
                                       kkt_multiplier=float(sol.x[-1]),
                                       product_error=float(np.max(np.abs(pnum-v*v))))
                            rows.append(row)
                            tag=f"n{n}_q{q}_S{S}_c{cX}_z{z}_k{k}"
                            check(tag+"_amplitude_KKT", sol.success and res<2e-8,
                                  residual=res,solver_success=bool(sol.success))
                            check(tag+"_constrained_match", rerr<2e-8 and verr<2e-9,
                                  r_error=rerr,energy_error=verr)
                            check(tag+"_product", np.max(np.abs(pnum-v*v))<2e-8)
                            check(tag+"_positive_metric", Dnum>0 and schur>0 and abs(cI*k)<1,
                                  hidden_metric=float(Dnum),clock_schur=float(schur))
                            check(tag+"_envelope", derr<2e-8 and abs(sol.x[-1]+deriv)<2e-8,
                                  derivative_fd_error=derr)
                            # Dot full r-gradient with u: exactly the nonzero modulus force.
                            force=g*g*A@rr
                            force[-1]-=S*cX/D**2
                            check(tag+"_no_full_stationarity",
                                  abs(float(u@force)-deriv)<2e-10 and abs(deriv)>1e-5)
    return rows


def extended_branch_checks():
    rows=[]
    for z in (-100.,-20.,-5.,0.,5.,20.,100.):
        Q,u,A,Ap,h,D,t,rr,V=numeric_branch(4,3,0.5,0.3,0.7,z)
        B=1+0.3*u[-1]*z
        grad=0.7**2*A@rr
        grad[-1]-=0.5*0.3/D**2
        perp=grad-u*(u@grad)
        check(f"extended_z{z}_unique_positive_branch",D>0 and t>0 and
              abs(t*D*D-0.5*0.3**2*h/0.7**2)<2e-10)
        check(f"extended_z{z}_perp_stationarity",np.max(np.abs(perp))<2e-9)
        rows.append(dict(z=z,B=float(B),D=float(D),t=float(t),potential=float(V),
                         clock_metric_k0=float(1+0.2*rr[0]),
                         derivative=float(-0.5*0.3*u[-1]/D**2)))
    # This diagnostic deliberately reports, rather than discards, the eventual
    # kinetic-domain boundary along the formal algebraic runaway.
    zboundary=brentq(lambda zz:1+0.2*numeric_branch(4,3,0.5,0.3,0.7,zz)[7][0],-1000.,0.)
    check("clock_metric_boundary_finite",-1000<zboundary<0)
    return rows,float(zboundary)


def main():
    exact_checks()
    rows=nonlinear_checks()
    extended,boundary=extended_branch_checks()
    source_names=["higgs_finite_auxiliary_20260926.md","higgs_extra_gauge_20260926.md"]
    payload={
        "scope":"Internal exact and numerical audit; no full SUGRA or nonzero-X/Z branch claim.",
        "date":"2026-09-26",
        "script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "input_sha256":{name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in source_names},
        "versions":{"python":platform.python_version(),"numpy":np.__version__,"scipy":scipy.__version__,"sympy":sp.__version__},
        "total":len(checks),"passed":sum(x["passed"] for x in checks),
        "failed":[x for x in checks if not x["passed"]],
        "max_KKT_residual":max(x["residual"] for x in rows),
        "max_energy_error":max(x["energy_error"] for x in rows),
        "max_r_error":max(x["r_error"] for x in rows),
        "max_derivative_fd_error":max(x["derivative_fd_error"] for x in rows),
        "min_hidden_metric":min(x["hidden_metric"] for x in rows),
        "min_clock_schur":min(x["clock_schur"] for x in rows),
        "nonlinear_solve_count":len(rows),"nonlinear_rows":rows,
        "extended_branch_diagnostics":extended,"clock_boundary_n4_q3_cI02":boundary,
        "checks":checks,
    }
    OUT.write_text(json.dumps(payload,ensure_ascii=False,indent=2)+"\n")
    print(json.dumps({k:v for k,v in payload.items() if k not in ("checks","nonlinear_rows","extended_branch_diagnostics")},indent=2))
    return 0 if payload["passed"]==payload["total"] else 1


if __name__=="__main__":
    raise SystemExit(main())
