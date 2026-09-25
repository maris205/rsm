#!/usr/bin/env python3
"""Independent two-real-scalar determinant diagnostic, not a SUGRA spectrum.

Only mpmath and the standard library are required. No main experiment code is
imported. Exact eigenvalues, divided differences, numerical parameter integrals,
and Schur complements independently evaluate the same specified scalar model.
"""
from pathlib import Path
import hashlib
import json
import platform

import mpmath as mp

HERE = Path(__file__).resolve()
OUT = HERE.with_suffix('.json')
checks = []


def m(value):
    return mp.mpf(str(value))


def check(name, actual, expected, tolerance='1e-40', floor='1e-100'):
    residual = abs(actual - expected) / max(abs(expected), m(floor))
    checks.append(dict(name=name, relative_residual=str(residual),
                       tolerance=tolerance, passed=bool(residual < m(tolerance))))


def f(x, mu2):
    return x*x*(mp.log(x/mu2) - m('1.5'))


def fp(x, mu2):
    return 2*x*(mp.log(x/mu2) - 1)


def divided_log(a, b, mu2):
    """[a(log(a/mu²)-1)-b(log(b/mu²)-1)]/(a-b)."""
    if a == b:
        return mp.log(a/mu2)
    # log1p avoids cancellation when a and b approach each other.
    e = (a-b)/b
    return mp.log(b/mu2) + (1+e)*mp.log1p(e)/e - 1


def q2_coefficient(a, b, u, v, w, mu2):
    return u*fp(a, mu2) + w*fp(b, mu2) + 2*v*v*divided_log(a, b, mu2)


def exact_trace(q, a, b, u, v, w, mu2):
    aa, bb = a+u*q*q, b+w*q*q
    root = mp.sqrt((aa-bb)**2+4*v*v*q*q)
    hi = (aa+bb+root)/2
    # determinant/high eigenvalue is stable for the hierarchical light eigenvalue.
    lo = (aa*bb-v*v*q*q)/hi
    return f(hi, mu2)+f(lo, mu2)


def extracted_coefficient(a, b, u, v, w, mu2, h):
    base = f(a, mu2)+f(b, mu2)
    c1 = (exact_trace(h, a, b, u, v, w, mu2)-base)/(h*h)
    c2 = (exact_trace(h/2, a, b, u, v, w, mu2)-base)/(h*h/4)
    return (4*c2-c1)/3


def numeric_integral(a, b, mu2):
    lo, hi = min(a, b), max(a, b)
    r = lo/hi
    points = [mp.mpf(0)]
    if r < m('.01'):
        points += [r, mp.sqrt(r)]
    points += [mp.mpf(1)]
    return mp.quad(lambda x: mp.log((lo+(hi-lo)*x)/mu2), points)


def run(precision):
    mp.mp.dps = precision
    pref = f'dps{precision}'
    cases = [
        ('separated', '2', '7', '.3', '.4', '.8', '5'),
        ('swapped', '7', '2', '.8', '.4', '.3', '5'),
        ('degenerate', '3', '3', '.3', '.4', '.8', '5'),
        ('near_degenerate', '3.000000000000000000000000000003', '3', '.3', '.4', '.8', '5'),
        ('hierarchical', '48', '1e22', '1e-34', '3e-6', '0', '48'),
    ]
    for name, *raw in cases:
        a, b, u, v, w, mu2 = map(m, raw)
        analytic = q2_coefficient(a, b, u, v, w, mu2)
        exact = extracted_coefficient(a, b, u, v, w, mu2, m('1e-12'))
        check(f'{pref}:{name}:exact_eigenvalues_vs_divided_difference', exact, analytic)
        check(f'{pref}:{name}:parameter_integral', numeric_integral(a, b, mu2),
              divided_log(a, b, mu2))
        check(f'{pref}:{name}:basis_swap', q2_coefficient(b, a, w, v, u, mu2), analytic)
        check(f'{pref}:{name}:explicit_mu_derivative',
              mp.diff(lambda z: q2_coefficient(a, b, u, v, w, mp.exp(2*z)), mp.log(mu2)/2),
              -4*(u*a+w*b+v*v))
        if a == b:
            check(f'{pref}:{name}:regular_degenerate_limit', analytic,
                  (u+w)*fp(a, mu2)+2*v*v*mp.log(a/mu2))
        for kval in ['0', '1', '1e22']:
            k2, q = m(kval), m('.01')
            x, y, z = k2+a+u*q*q, k2+b+w*q*q, v*q
            direct = mp.log(x*y-z*z)
            check(f'{pref}:{name}:schur_A:k2={kval}', mp.log(x)+mp.log(y-z*z/x), direct)
            check(f'{pref}:{name}:schur_B:k2={kval}', mp.log(y)+mp.log(x-z*z/y), direct)

    mpl, mg, s = m('2.435e27'), m(1), m('1e22')
    p, a = mpl*mpl, 48*mg*mg
    c = mp.sqrt(m(2)/3)/mpl
    u, v = s/(3*p), -c*s
    check(f'{pref}:mass_law_first_derivative', mp.diff(lambda sig: s*mp.exp(-c*sig), 0), v)
    check(f'{pref}:mass_law_second_derivative', mp.diff(lambda sig: s*mp.exp(-c*sig), 0, 2), 2*u)
    rows = []
    for label, mu2 in [('modulus_scale', a), ('charged_scale', s)]:
        la, lb = mp.log(a/mu2), mp.log(s/mu2)
        tad = u*a*(la-1)/(16*mp.pi**2)
        mix = v*v*divided_log(a, s, mu2)/(16*mp.pi**2)
        total = tad+mix
        exact = extracted_coefficient(a, s, u, v, m(0), mu2, m('1e-6'))/(32*mp.pi**2)
        check(f'{pref}:physical_hierarchy:{label}:exact_eigenvalues', exact, total)
        asymptotic = v*v*(lb-1)/(16*mp.pi**2)
        check(f'{pref}:physical_hierarchy:{label}:hierarchy_asymptotic', mix, asymptotic, '3e-19')
        rows.append(dict(scale=label, mu_eV=str(mp.sqrt(mu2)),
                         log_A_over_mu2=str(la), log_B_over_mu2=str(lb),
                         divided_log=str(divided_log(a, s, mu2)),
                         tadpole_delta_mass2_eV2=str(tad), mixed_delta_mass2_eV2=str(mix),
                         total_delta_mass2_eV2=str(total),
                         absolute_mixed_over_tadpole=str(abs(mix/tad)),
                         relative_total_over_tree_mass2=str(total/s)))
    beta = (u*a+v*v)/(8*mp.pi**2)
    # Frozen tree coefficients suffice at this one-loop order; their own running
    # inside the loop contributes at the next loop order.
    mu1, mu2 = mp.sqrt(a), mp.sqrt(s)
    c_run = beta*mp.log(mu2/mu1)
    delta1, delta2 = [m(row['total_delta_mass2_eV2']) for row in rows]
    check(f'{pref}:explicit_scale_cancellation_with_mass_counterterm', delta2+c_run, delta1)
    return dict(precision_digits=precision, Mpl_eV=str(mpl), mG_eV=str(mg),
                modulus_mass_eV=str(mp.sqrt(a)), charged_mass_eV=str(mp.sqrt(s)),
                charged_over_modulus_mass=str(mp.sqrt(s/a)),
                log_B_over_A=str(mp.log(s/a)), u=str(u), v_eV=str(v),
                reduced_scalar_beta_mass2_eV2=str(beta),
                counterterm_shift_between_scales_eV2=str(c_run), rows=rows)


if __name__ == '__main__':
    batches = [run(160), run(220)]
    for r1, r2 in zip(batches[0]['rows'], batches[1]['rows']):
        for key in ['tadpole_delta_mass2_eV2', 'mixed_delta_mass2_eV2', 'total_delta_mass2_eV2']:
            check(f'precision_repeat:{r1["scale"]}:{key}', m(r1[key]), m(r2[key]), '1e-140')
    result = dict(status='PASS' if all(c['passed'] for c in checks) else 'FAIL',
                  scope='Flat-kinetic two-real-scalar, constant-background 1PI diagnostic; not full SUGRA or a physical soft splitting.',
                  deterministic=True, statistical_tests='none',
                  python=platform.python_version(), mpmath=mp.__version__,
                  script_sha256=hashlib.sha256(HERE.read_bytes()).hexdigest(),
                  passed=sum(c['passed'] for c in checks), total=len(checks),
                  numeric_batches=batches, checks=checks)
    OUT.write_text(json.dumps(result, indent=2, ensure_ascii=False)+'\n')
    print(json.dumps(dict(status=result['status'], passed=result['passed'], total=result['total'], output=str(OUT))))
    raise SystemExit(0 if result['status']=='PASS' else 1)
