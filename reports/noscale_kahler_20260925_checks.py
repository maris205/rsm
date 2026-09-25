#!/usr/bin/env python3
"""Independent symbolic audit of the quartic no-scale candidate.

All symbolic expressions use M_P=1.  This script does not import the main
experiment, and it does not integrate a cosmological background.
"""
import json
from pathlib import Path
import sympy as s

checks = []


def check(name, expr):
    residual = s.factor(s.cancel(expr))
    passed = residual == 0
    checks.append({"name": name, "passed": bool(passed), "residual": str(residual)})
    if not passed:
        raise AssertionError((name, residual))


O, p, pb, q, u, v, W, Wb, j, jb = s.symbols("O p pb q u v W Wb j jb")
E = p*pb-q*(O+u*v/3)
D = p*pb-O*q
g = s.Matrix([[3*D/O**2, -p*v/O**2],
              [-pb*u/O**2, 1/O+u*v/(3*O**2)]])
check("metric_determinant", g.det()-3*E/O**3)
DW = s.Matrix([-3*p*W/O, j+u*W/O])
DWb = s.Matrix([-3*pb*Wb/O, jb+v*Wb/O])
direct = ((DWb.T*g.inv()*DW)[0]-3*W*Wb)/O**3
closed = j*jb/O**2+q*(3*W-v*j)*(3*Wb-u*jb)/(3*O**2*E)
check("full_complex_T_and_Z_potential", direct-closed)
check("real_Z_axis_potential", closed.subs({u: 0, v: 0})-(j*jb/O**2+3*q*W*Wb/(O**2*D)))

d, n, a, b, eta, U, ww = s.symbols("delta nu a b eta U ww", real=True)
f = 1+d+eta*d**2+a*d**4+b*n**4
pt = 1+2*eta*d+4*a*d**3-4*s.I*b*n**3
ptb = s.conjugate(pt)
qt = 2*eta+12*a*d**2+12*b*n**2
Dt = pt*ptb-f*qt
V = U/f**2+3*qt*ww/(f**2*Dt)
V0 = V.subs(eta, 0)
origin = {d: 0, n: 0}
check("quartic_center_preserves_V_equals_U", V0.subs(origin)-U)
check("clock_tadpole_real_T", s.diff(V0, d).subs(origin)+2*U)
check("imaginary_T_tadpole_zero", s.diff(V0, n).subs(origin))
check("real_T_canonical_mass_clockoff", s.diff(V0.subs(U, 0), d, 2).subs(origin)*s.Rational(2, 3)-48*a*ww)
check("imaginary_T_canonical_mass_clockoff", s.diff(V0.subs(U, 0), n, 2).subs(origin)*s.Rational(2, 3)-48*b*ww)
check("finite_clock_real_coordinate_Hessian_at_center", s.diff(V0, d, 2).subs(origin)*s.Rational(2, 3)-(48*a*ww+4*U))
check("finite_clock_imag_coordinate_Hessian_at_center", s.diff(V0, n, 2).subs(origin)*s.Rational(2, 3)-48*b*ww)
check("quartic_metric_center", (3*Dt/f**2).subs({**origin, eta: 0})-3)

f0, pp, qq = s.symbols("f0 pp qq")
D0 = pp-f0*qq
mc = s.symbols("mbar")
check("compensator_cancellation", mc-mc*pp/D0+mc*f0*qq/D0)

# Expand the exact full matter potential, treating the two independent
# quadratic invariants R=|Q+|^2+|Q-|^2 and P=Q+Q- as formal variables.
R, P, Pb, M, Mb, Mz, Mzb = s.symbols("R P Pb M Mb Mz Mzb")
Om = f0-R/3
vmatter = (U+jb*Mz*P+j*Mzb*Pb+M*Mb*R)/Om**2
vmatter += qq*(3*W+M*P)*(3*Wb+Mb*Pb)/(3*Om**2*D0)
vbase = U/f0**2+3*qq*W*Wb/(f0**2*D0)
zeroq = {R: 0, P: 0, Pb: 0}
diag = s.diff(vmatter, R).subs(zeroq)*f0
hol = s.diff(vmatter, P).subs(zeroq)*f0
check("charged_canonical_common_soft", diag-M*Mb/f0-2*vbase/3)
check("charged_canonical_B", hol-(jb*Mz/f0+qq*M*Wb/(f0*D0)))
check("charged_B_at_exact_center", hol.subs(qq, 0)-jb*Mz/f0)
check("charged_fermion_mass", (f0**(-s.Rational(3,2))*M/(1/f0))**2-M**2/f0)

check("clock_shift_linear_equation", -2*U+72*a*ww*(U/(36*a*ww)))
check("shifted_q_leading", (12*a*d**2).subs(d, U/(36*a*ww))-U**2/(108*a*ww**2))
check("clock_valley_energy_leading", (-2*U*d+36*a*ww*d**2).subs(d, U/(36*a*ww))+U**2/(36*a*ww))

check("eta_center_potential", V.subs(origin)-(U+6*eta*ww/(1-2*eta)))
check("eta_center_real_T_tadpole", s.diff(V, d).subs(origin)-(-2*U-12*eta*ww*(1-eta)/(1-2*eta)**2))
check("eta_vacuum_shift_leading", -12*eta*ww+72*a*ww*(eta/(6*a)))
check("eta_combined_shift_leading", -2*U-12*eta*ww+72*a*ww*(eta/(6*a)+U/(36*a*ww)))

# At q=0, a shift-symmetric clock's imaginary direction is not made heavy
# by Wc: K keeps Omega=1-F^2*y^2/(3 M_P^2).
y, F = s.symbols("y F", real=True)
check("imaginary_clock_mass", s.diff(U/(1-F**2*y**2/3)**2, y, 2).subs(y, 0)/F**2-4*U/3)

A, B, wr = s.symbols("Wc_real Wc_imag w_real", real=True)
Sy = 3*(A+s.I*B)+wr*(3+2*s.I*y)*s.exp(-s.I*y)
Syb = s.conjugate(Sy)
S2 = s.expand(Sy*Syb)
check("uplift_clock_S_y", s.diff(S2, y).subs(y, 0)+6*B*wr)
check("uplift_clock_S_yy", s.diff(S2, y, 2).subs(y, 0)-(8*wr**2+6*A*wr))
Vy = U/(1-F**2*y**2/3)**2+q*S2/(3*(1-F**2*y**2/3)**2*(1-q-q*F**2*y**2/3))
check("uplift_y_tadpole", s.diff(Vy, y).subs(y, 0)+2*q*B*wr/(1-q))
vyy_expected = 4*F**2*U/3+q/(3*(1-q))*(8*wr**2+6*A*wr+6*F**2*(2+q/(1-q))*((A+wr)**2+B**2))
check("uplift_y_Hessian", s.diff(Vy, y, 2).subs(y, 0)-vyy_expected)

# The real-axis identity m_T^2=48*a*|W|^2 must not be differentiated
# transversely by replacing W only. Full geometry inserts S(y), and reduces
# the transverse linear Wc-w cross term by a factor of three.
f_quartic = f.subs(eta, 0)
q_quartic = qt.subs(eta, 0)
omega_y = f_quartic-F**2*y**2/3
E_y = (pt*ptb).subs(eta, 0)-q_quartic*(f_quartic+F**2*y**2/3)
V_complex = U/omega_y**2+q_quartic*S2/(3*omega_y**2*E_y)
V_y_axis = s.factor(s.diff(V_complex, y).subs(y, 0))
check("full_potential_real_T_mass_y_derivative", s.diff(V_y_axis, d, 2).subs(origin)*s.Rational(2,3)+32*a*B*wr)
check("full_potential_imag_T_mass_y_derivative", s.diff(V_y_axis, n, 2).subs(origin)*s.Rational(2,3)+32*b*B*wr)
W_y = A+s.I*B+wr*s.exp(-s.I*y)
check("heavy_mass_y_derivative_one_third_gravitino", s.diff(S2/9, y).subs(y, 0)-s.diff(W_y*s.conjugate(W_y), y).subs(y, 0)/3)

# Constant-W hidden-clock fermion at the center: D_Z W=0 and the
# connection cancels K_ZZ W, so the uneaten clockino has zero mass.
check("clockino_mass_clockoff_center", -W/f0+(pp/(3*D0))*(3*W/f0)-W*qq/D0)
check("clockino_mass_zero_q", (W*qq/D0).subs(qq, 0))

x, xp, mu2, c = s.symbols("x xp mu2 c", positive=True)
cw = (c*x)**2*(s.log(c*x/mu2)-s.Rational(3, 2))/(64*s.pi**2)
check("modulus_scalar_CW_chain_rule", s.diff(cw, x)*xp-c**2*x*xp*(s.log(c*x/mu2)-1)/(32*s.pi**2))
check("eta_vacuum_counterterm_slope", (s.diff(cw, x)-cw/x)*xp-c**2*x*xp*(s.log(c*x/mu2)-s.Rational(1,2))/(64*s.pi**2))

out = {"date": "2026-09-25", "method": "independent SymPy analytic checks",
       "units": "M_P=1 in symbolic expressions", "checks": checks,
       "passed": sum(c["passed"] for c in checks), "total": len(checks),
       "scope": "local classical geometry/spectrum and scalar CW chain rule; not a full SUGRA loop or cosmology calculation",
       "errors": []}
target = Path(__file__).with_suffix(".json")
target.write_text(json.dumps(out, indent=2)+"\n")
print(f'{out["passed"]}/{out["total"]} passed; {target.name}')
