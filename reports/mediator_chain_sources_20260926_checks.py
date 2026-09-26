#!/usr/bin/env python3
"""Small exact-algebra audit accompanying the bounded source review.

No experiment implementation is imported. This checks matching bookkeeping,
not loop protection, the physical SUGRA spectrum, or observational fit.
"""
import json
from pathlib import Path
import sympy as sp
import mpmath as mp

checks = []
def check(name, passed):
    checks.append({"name": name, "passed": bool(passed)})

q = sp.symbols("q", positive=True)
for n in (2, 3, 4, 5):
    h = sp.diag(*([q*q+1]*n))
    for i in range(n-1):
        h[i, i+1] = h[i+1, i] = -q
    inv = h.inv()
    dn = sum(q**(2*j) for j in range(n+1))
    dprev = sum(q**(2*j) for j in range(n))
    check(f"n{n}_determinant", sp.simplify(h.det()-dn) == 0)
    check(f"n{n}_endpoint_self", sp.simplify(inv[0,0]-dprev/dn) == 0)
    check(f"n{n}_endpoint_cross", sp.simplify(inv[0,n-1]-q**(n-1)/dn) == 0)
    check(f"n{n}_reflection", sp.simplify(inv[0,0]-inv[n-1,n-1]) == 0)
    # Boundary-pinned sum of local squares, unlike the zero-mode clockwork.
    v = sp.symbols(f"v0:{n}")
    local = q*q*v[0]**2+sum((v[i]-q*v[i+1])**2 for i in range(n-1))+v[-1]**2
    check(f"n{n}_local_square_decomposition", sp.expand((sp.Matrix(v).T*h*sp.Matrix(v))[0]-local) == 0)
    unpinned = h.copy()
    unpinned[0,0] -= q*q
    unpinned[-1,-1] -= 1
    null = sp.Matrix([q**(-i) for i in range(n)])
    check(f"n{n}_unpinned_zero_mode", all(sp.simplify(x) == 0 for x in unpinned*null))

j, sx, tau, zeta, a = sp.symbols("j sx tau zeta a", real=True)
check("complete_square_retains_self", sp.expand(j*j+2*tau*j*sx+zeta*sx*sx-(j+tau*sx)**2-(zeta-tau*tau)*sx*sx) == 0)
x, xb, y, f = sp.symbols("x xb y f", real=True)
contact = f*f*y*y*x*xb
check("bypass_has_mixed_nonholomorphic_curvature", sp.diff(contact,y,2,x,xb) == 2*f*f)

mp.mp.dps = 100
target = mp.mpf("2.9251519294813043e-61")
cases = []
def tratio(qv, n):
    return (qv*qv-1)*qv**(n-1)/(qv**(2*n)-1)
for qq in (2,3,4):
    qv = mp.mpf(qq)
    n = 2
    while tratio(qv,n) > target:
        n += 1
    t = tratio(qv,n)
    amplitude = target/t
    check(f"q{qq}_target_bracket", t <= target < tratio(qv,n-1))
    check(f"q{qq}_order_one_endpoint_ratio", 1 <= amplitude < qv)
    check(f"q{qq}_positive_spectral_gap", qv*qv+1-2*qv*mp.cos(mp.pi/(n+1)) > (qv-1)**2)
    check(f"q{qq}_self_exceeds_rank_one", 1-t*t > 0)
    cases.append({"q":qq,"n":n,"tau_equal_endpoint_couplings":mp.nstr(t,35),"target_gX_over_g0":mp.nstr(amplitude,35),"uniform_q_plus_one_percent_tau_ratio":mp.nstr(tratio(qv*mp.mpf('1.01'),n)/t,35)})

result={"scope":"Own exact algebra and illustrative hierarchy bookkeeping only; no all-order protection or full mediator matching computed", "precision_digits":100,"passed":sum(c['passed'] for c in checks),"total":len(checks),"failures":[c for c in checks if not c['passed']],"checks":checks,"target_tau":mp.nstr(target,35),"illustrative_cases":cases}
out=Path(__file__).with_suffix('.json')
out.write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n')
print(json.dumps({"passed":result['passed'],"total":result['total'],"cases":cases},ensure_ascii=False))
assert not result['failures']
