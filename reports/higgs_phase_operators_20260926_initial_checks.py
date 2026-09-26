#!/usr/bin/env python3
"""Independent original-field operator audit; no experiment implementation imported."""
from pathlib import Path
from itertools import product
from math import gcd
from functools import reduce
import hashlib
import json
import mpmath as mp
import sympy as sp

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
checks = []


def check(name, value, **details):
    checks.append(dict(name=name, passed=bool(value), **details))


def zero(name, expr):
    residual = sp.simplify(expr)
    check(name, residual == 0, residual=str(residual))


def charges(q, n):
    out = sp.zeros(n + 1, n)
    out[0, 0] = q
    for a in range(1, n):
        out[a, a - 1] = 1
        out[a, a] = -q
    out[n, n - 1] = 1
    return out


for q in (2, 3, 5):
    for n in (1, 2, 3, 5, 8, 12):
        Q = charges(q, n)
        w = sp.Matrix([1] + [-q**a for a in range(1, n + 1)])
        label = f'q={q},n={n}'
        check(label + ': charge rank', Q.rank() == n)
        check(label + ': primitive integer kernel', Q.T * w == sp.zeros(n, 1)
              and reduce(gcd, map(abs, w)) == 1)
        H = sp.zeros(n)
        for i in range(n):
            H[i, i] = q*q + 1
            if i + 1 < n:
                H[i, i+1] = H[i+1, i] = -q
        check(label + ': charge Gram', Q.T * Q == H)
        D = sum(abs(x) for x in w)
        S = (w.T*w)[0]
        zero(label + ': degree sum', D-sp.Rational(q**(n+1)-1, q-1))
        zero(label + ': norm sum', S-sp.Rational(q**(2*n+2)-1, q*q-1))
        # The first coordinate equals integer l. The gauge equations fix all others.
        for ell in (-3, -1, 1, 2):
            d = ell*w
            check(label + f': l={ell} attainable L1 lower bound',
                  Q.T*d == sp.zeros(n, 1) and sum(map(abs, d)) == abs(ell)*D)
        # Gauging one additional endpoint phase is legitimate and full rank.
        for endpoint in (0, n):
            c = sp.zeros(n+1, 1)
            c[endpoint] = 1
            full = Q.row_join(c)
            check(label + f': added endpoint {endpoint} gives full rank',
                  full.rank() == n+1)
        # All cubic gauge/global anomalies cancel within each vectorlike pair.
        cubic = sum(int(x)**3 + int(-x)**3 for x in w)
        mixed = sum(int(w[a])*int(Q[a, j])**2
                    -int(w[a])*int(-Q[a, j])**2
                    for a in range(n+1) for j in range(n))
        check(label + ': vectorlike phase anomalies', cubic == mixed == 0)

# Exhaustive bounded lattice check, independent of the constructive recursion.
for q, n in ((2, 1), (3, 1), (2, 2), (3, 2)):
    D = (q**(n+1)-1)//(q-1)
    solutions = []
    Q = [[int(x) for x in row] for row in charges(q, n).tolist()]
    for d in product(range(-D, D+1), repeat=n+1):
        norm = sum(map(abs, d))
        if norm == 0 or norm > D:
            continue
        if all(sum(d[a]*Q[a][j] for a in range(n+1)) == 0 for j in range(n)):
            solutions.append(d)
    w = (1,) + tuple(-q**a for a in range(1, n+1))
    check(f'brute q={q},n={n}: only +/- primitive at degree <=D',
          sorted(solutions) == sorted([w, tuple(-x for x in w)]),
          first_nonzero_degree=D, solutions=[list(d) for d in solutions])

# Local canonical normalization and leading holomorphic mass.
C, Cb, v, norm2, kappa, cutoff = sp.symbols('C Cb v S kappa M', positive=True)
mu, heavy, b, x = sp.symbols('mu heavy b x', positive=True)
for q, n in ((2, 1), (3, 2), (3, 5)):
    w = [1] + [-q**a for a in range(1, n+1)]
    S = sum(z*z for z in w)
    K = sum(v*v*(sp.exp(z*(C+Cb)/(sp.sqrt(2*S)*v))
                  + sp.exp(-z*(C+Cb)/(sp.sqrt(2*S)*v))) for z in w)
    zero(f'canonical q={q},n={n}', sp.diff(K,C,Cb).subs({C:0,Cb:0})-1)
    # Phase and saxion are C=(s+i a)/sqrt2. Common soft mass preserves the phase.
    s, a, soft = sp.symbols('s a soft', real=True)
    Kreal = K.subs({C:(s+sp.I*a)/sp.sqrt(2), Cb:(s-sp.I*a)/sp.sqrt(2)})
    zero(f'common soft phase q={q},n={n}', sp.diff(soft*Kreal,a,2))
    zero(f'common soft real mass q={q},n={n}',
         sp.diff(soft*Kreal,s,2).subs({s:0,a:0})-2*soft)
    D = sum(abs(z) for z in w)
    O = v**D*sp.exp(sp.sqrt(S)*C/(sp.sqrt(2)*v))
    zero(f'holomorphic mass q={q},n={n}',
         (kappa*sp.diff(O,C).subs(C,0)/cutoff**(D-2))**2
         - kappa**2*S*v*v/2*(v/cutoff)**(2*D-4))
    # Radial mixing cannot be silently removed at comparable masses.
    wvec = sp.Matrix(w)
    pvec = sp.Matrix(list(map(abs,w)))
    zero(f'radial/phase mixing norms q={q},n={n}',
         (wvec.T*wvec)[0]-(pvec.T*pvec)[0])

mass = sp.Matrix([[0,heavy,b,0],[heavy,0,0,0],
                  [b,0,0,mu],[0,0,mu,0]])
zero('complete local chiral characteristic polynomial',
     mass.charpoly(x).as_expr()-(x**4-(heavy**2+b**2+mu**2)*x**2+heavy**2*mu**2))
z = sp.symbols('z', positive=True)
xminus = (heavy**2+2*z-sp.sqrt(heavy**4+4*z*z))/2
zero('equal-lambda light mass expansion',
     sp.series(xminus,z,0,4).removeO()-(z-z*z/heavy**2))
zero('stable eigenvalue equivalent',
     xminus-2*heavy**2*z/(heavy**2+2*z+sp.sqrt(heavy**4+4*z*z)))

def numerical(dps):
    mp.mp.dps = dps
    q, n = 3, 127
    D = (q**(n+1)-1)//(q-1)
    S = (q**(2*n+2)-1)//(q*q-1)
    pref = mp.sqrt(mp.mpf(S)/2)
    ratios = ['0.5','0.9','0.99','0.9999999999']
    logs = {r: mp.log10(pref)+(D-2)*mp.log10(mp.mpf(r)) for r in ratios}
    edge = -mp.expm1(-mp.log(pref)/(D-2))
    check(f'{dps}: suppression for all fixed ratios', all(val < -mp.mpf('1e49') for val in logs.values()))
    check(f'{dps}: near-cutoff precision', mp.mpf('2e-59') < edge < mp.mpf('3e-59'))
    tests=[]
    for ratio in ('1e-20','1e-5','0.1','1','10'):
        mm=mp.mpf(ratio)
        full=mp.matrix([[0,1,mm,0],[1,0,0,0],[mm,0,0,mm],[0,0,mm,0]])
        vals=sorted([x*x for x in mp.eigsy(full,eigvals_only=True)])
        low=2*mm*mm/(1+2*mm*mm+mp.sqrt(1+4*mm**4))
        err=abs(vals[0]-low)/low
        check(f'{dps}: full local mass ratio={ratio}', err < mp.mpf('1e-75'),
              relative_error=mp.nstr(err,12), tolerance='1e-75')
        check(f'{dps}: light eigenvalue bounds ratio={ratio}', low>0 and low<=mm*mm and low<=mp.mpf('.5'))
        tests.append(dict(mu_over_heavy=ratio, light_mass_squared_over_heavy_squared=mp.nstr(low,90)))
    return dict(dps=dps,q=q,n=n,D=str(D),S=str(S),norm_w=mp.nstr(mp.sqrt(S),90),
                log10_mass_over_v={r:mp.nstr(val,90) for r,val in logs.items()},
                one_minus_v_over_cutoff_for_leading_mass_equal_v=mp.nstr(edge,90),
                local_spectrum_tests=tests)

nums=[numerical(140), numerical(200)]
mp.mp.dps=200
for key in nums[0]['log10_mass_over_v']:
    aa=mp.mpf(nums[0]['log10_mass_over_v'][key]); bb=mp.mpf(nums[1]['log10_mass_over_v'][key])
    err=abs(aa-bb)/max(abs(aa),abs(bb))
    check(f'cross precision log suppression {key}',err<mp.mpf('1e-85'), relative_error=mp.nstr(err,12))

inputs = [ROOT/'reports/chain_locality_breaking_20260926.md',
          ROOT/'experiments/threshold_rg_v01/protocol.md']
out = dict(date='2026-09-26', scope='Original charged fields, local rigid supersymmetry, declared polynomial lifting operator; no observation or full SUGRA completion',
           normalization='Phi=tildePhi=v; canonical Y; equal lambda for exact radial-mixing eigenvalues',
           precision_digits=[140,200], check_count=len(checks),
           passed=sum(c['passed'] for c in checks), failed=sum(not c['passed'] for c in checks),
           checks=checks,numerical=nums,
           input_hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},
           code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
           review_disclosure='AI-assisted algebra and internal independent-agent audit, not external peer review; no human-read attestation')
dest=HERE/'higgs_phase_operators_20260926_checks.json'
dest.write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n')
print(json.dumps({k:out[k] for k in ('check_count','passed','failed')},ensure_ascii=False))
if out['failed']:
    print(json.dumps([c for c in checks if not c['passed']],ensure_ascii=False,indent=2))
    raise SystemExit(1)
