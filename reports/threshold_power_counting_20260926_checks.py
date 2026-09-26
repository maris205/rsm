#!/usr/bin/env python3
"""Independent low-EFT power counting; no main/parent implementation imported."""
from pathlib import Path
import hashlib
import json
import mpmath as mp
import sympy as sp

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
INPUT = ROOT / 'experiments/mediator_chain_v01/results/inputs.json'
checks = []


def symbolic(name, value):
    ok = sp.simplify(value) == 0
    checks.append({'name': name, 'passed': bool(ok), 'residual': str(sp.simplify(value))})


def close(name, a, b, tol='1e-65'):
    err = abs(a-b)/max(abs(a), abs(b), mp.mpf('1e-150'))
    checks.append({'name': name, 'passed': bool(err < mp.mpf(tol)),
                   'relative_error': mp.nstr(err, 12), 'tolerance': tol})


x, xb, La, f, fb, S, r, rho, mu = sp.symbols('x xb La f fb S r rho mu', positive=True)
K = x*xb - (x*xb)**2/La**2
g = sp.diff(K, x, xb)
symbolic('Kahler metric', g-(1-4*x*xb/La**2))
V = S/g
symbolic('complex scalar mass', sp.diff(V, x, xb).subs({x:0, xb:0})-4*S/La**2)
symbolic('origin quartic real radial potential', sp.expand(sp.series((S/(1-2*r*r/La**2)), r, 0, 6).removeO()).coeff(r,4)-4*S/La**4)
Gam = sp.diff(g,x)/g
symbolic('connection-generated Yukawa', -Gam*f-4*f*xb/(La**2*g))
R = sp.diff(g,x,xb)-sp.diff(g,x)*sp.diff(g,xb)/g
symbolic('four-fermion curvature', R+4/(La**2*g))
symbolic('zero-momentum four-Goldstino cancellation', (2*f/La**2)*(2*fb/La**2)/(4*f*fb/La**2)-1/La**2)
radial = r-r**3/(3*La**2)-r**5/(10*La**4)
symbolic('canonical radial derivative through quartic', sp.series(sp.diff(radial,r)**2-(1-2*r*r/La**2),r,0,6).removeO())
rinv = rho+rho**3/(3*La**2)
symbolic('radial inverse through cubic', sp.series(radial.subs(r,rinv)-rho,rho,0,4).removeO())
Vnormal = sp.series((S/(1-2*r*r/La**2)).subs(r,rinv),rho,0,6).removeO().expand()
symbolic('normal-coordinate quartic', Vnormal.coeff(rho,4)-sp.Rational(16,3)*S/La**4)
ell = sp.log(4*S/(La**2*mu**2))
V1 = S**2/(2*sp.pi**2*La**4)*(ell-sp.Rational(3,2))
symbolic('dimension-eight beta cancellation', mu*sp.diff(V1,mu)+S**2/(sp.pi**2*La**4))
symbolic('one-loop force', sp.diff(V1,S)-S/(sp.pi**2*La**4)*(ell-1))

# Wick factors are enumerated over real flavors, not copied from the final N=2 formula.
N = 2
wick_pot = sum(1+2*int(a==b) for a in range(N) for b in range(N))
wick_der = sum(1 for a in range(N) for b in range(N))
wick_projected = sum(int(a==b) for a in range(N) for b in range(N))
symbolic('potential Wick coefficient', sp.Integer(wick_pot)-8)
symbolic('derivative Wick coefficient', sp.Integer(wick_der)-4)
symbolic('projected derivative Wick coefficient', sp.Integer(wick_projected)-2)
symbolic('coordinate-invariant scalar diagram sum', sp.Rational(4,3)*wick_pot+sp.Rational(2,3)*(wick_der-wick_projected)-(wick_pot+wick_der))

par = json.loads(INPUT.read_text())['parameters']


def calculate(dps):
    mp.mp.dps = dps
    P = mp.mpf(str(par['P']))
    Lam = mp.mpf('2e12')
    mg = mp.mpf('1e-15')
    ss = 3*P*mg**2
    m2 = 4*ss/Lam**2
    log = mp.log(m2/Lam**2)
    C = 16*mp.pi**2
    yuk2 = 16*ss/Lam**4
    quartic = 96*ss/Lam**4
    tad = m2/C*(log-1)
    v1 = m2**2/(2*C)*(log-mp.mpf('1.5'))
    v2pot = 8*m2/Lam**2*tad**2
    v2der = 4*m2/Lam**2*tad**2
    v2 = v2pot+v2der
    v2closed = 3*ss**3/(mp.pi**4*Lam**8)*(log-1)**2
    close(f'{dps}: two-loop scalar equivalent formula', v2, v2closed)
    # Direct derivatives at fixed mu=Lambda. This is a force w.r.t. the source, not chi.
    fun1 = lambda z: z**2/(2*mp.pi**2*Lam**4)*(mp.log(4*z/Lam**4)-mp.mpf('1.5'))
    fun2 = lambda z: 3*z**3/(mp.pi**4*Lam**8)*(mp.log(4*z/Lam**4)-1)**2
    d1 = mp.diff(fun1,ss)
    d2 = mp.diff(fun2,ss)
    d1closed = ss/(mp.pi**2*Lam**4)*(log-1)
    d2closed = 3*ss**2/(mp.pi**4*Lam**8)*(3*(log-1)**2+2*(log-1))
    close(f'{dps}: CW source derivative', d1, d1closed)
    close(f'{dps}: scalar two-loop source derivative', d2, d2closed)
    q, n = mp.mpf(3), 127
    modes = [q*q+1-2*q*mp.cos(k*mp.pi/(n+1)) for k in range(1,n+1)]
    weights = [2*mp.sin(k*mp.pi/(n+1))**2/(n+1) for k in range(1,n+1)]
    h1 = mp.fsum(w/l for w,l in zip(weights,modes))
    h2 = mp.fsum(w/l**2 for w,l in zip(weights,modes))
    h1closed = (q**(2*n)-1)/(q**(2*n+2)-1)
    close(f'{dps}: endpoint completeness', mp.fsum(weights), mp.mpf(1))
    close(f'{dps}: inverse endpoint', h1, h1closed)
    M2 = Lam**2*h1
    c8tree = -2*h2/h1**2
    c8run = c8tree+mp.log(mp.sqrt(m2)/Lam)/mp.pi**2
    one_mu_low = -3*ss**2/(4*mp.pi**2*Lam**4)
    fixed_sum = c8tree*ss**2/Lam**4+v1
    run_sum = c8run*ss**2/Lam**4+one_mu_low
    close(f'{dps}: fixed-scale/running-coefficient potential identity', fixed_sum,run_sum)
    # Source-dependent mu requires differentiating the running coefficient too.
    runfun = lambda z: (c8tree+mp.log(mp.sqrt(4*z/Lam**2)/Lam)/mp.pi**2)*z**2/Lam**4-3*z**2/(4*mp.pi**2*Lam**4)
    close(f'{dps}: source-dependent mu total derivative',mp.diff(runfun,ss),2*c8tree*ss/Lam**4+d1)
    for frac in ['0.001','0.01','0.1','0.4']:
        rfrac = mp.mpf(frac)
        metric = 1-4*rfrac**2
        checks.append({'name':f'{dps}: metric at |X|/Lambda={frac}', 'passed':bool(metric>0), 'metric':mp.nstr(metric,20)})
    return {k: mp.nstr(v,60) for k,v in {
        'S_eV4':ss,'f_abs_eV2':mp.sqrt(ss),'sqrt_f_eV':ss**mp.mpf('0.25'),
        'mX_eV':mp.sqrt(m2),'mX2_eV2':m2,'mX2_over_Lambda2':m2/Lam**2,
        'ell':log,'lambda_O2':quartic,'y_abs_squared':yuk2,
        'epsilon_derivative':m2/(C*Lam**2),
        'epsilon_scalar':quartic/C,'epsilon_goldstino':yuk2/C,
        'scalar_log_parameter':quartic*abs(log)/C,
        'goldstino_log_parameter':yuk2*abs(log)/C,
        'VA_parameter_at_mX':m2**2/(C*ss),
        'V1_X_eV4':v1,'V2_potential_diagram_eV4':v2pot,
        'V2_derivative_diagram_eV4':v2der,'V2_scalar_diagrams_eV4':v2,
        'V2scalar_over_V1X':v2/v1,
        'dV2scalar_dS_over_dV1X_dS':d2/d1,
        'h1':h1,'h2':h2,'c8_tree':c8tree,'c8_low_mu':c8run,
        'CW_X_force_over_tree_force':d1/(2*c8tree*ss/Lam**4),
        'm_heavy_min_eV':mp.sqrt(M2*min(modes)),
        'm_heavy_max_eV':mp.sqrt(M2*max(modes)),
        'gX_loop_proxy':2/C,
        'max_endpoint_weight':max(weights),
    }.items()}


lo = calculate(90)
hi = calculate(140)
for key in lo:
    close('90/140 precision: '+key, mp.mpf(lo[key]),mp.mpf(hi[key]),tol='1e-55')
out = {'scope':'Independent origin-canonical low-EFT derivation and scalar-only two-loop diagnostic; not full SUSY or heavy matching.',
       'precision_digits':[90,140],
       'input_sha256':hashlib.sha256(INPUT.read_bytes()).hexdigest(),
       'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
       'input_path':str(INPUT.relative_to(ROOT)),
       'source_strength_convention':'clock-off S=3 P mg²; Lambda=2 TeV, mg=1e-15 eV, q=3, n=127',
       'results':hi,'checks':checks,
       'passed':sum(c['passed'] for c in checks),'total':len(checks),
       'failure_history':[],
       'external_sources':[
          {'url':'https://arxiv.org/pdf/0907.2441','read_scope':'sections 3.1-3.2, eqs 3.2-3.9 and 3.12-3.23; browser PDF extraction; no local page anchors',
           'use':'quartic Kahler sgoldstino example and low-energy nonlinear Goldstino matching only'},
          {'url':'https://arxiv.org/html/hep-ph/0111209','read_scope':'sections 1, 2.1-2.3 and 4, especially eqs 2.17-2.18 and 4.3,4.13',
           'use':'renormalized scalar tadpole convention and potential-only figure-eight normalization only'}]}
(HERE/'threshold_power_counting_20260926_checks.json').write_text(json.dumps(out,indent=2)+'\n')
print(f"{out['passed']}/{out['total']} passed")
print(json.dumps(hi,indent=2))
if out['passed'] != out['total']:
    raise SystemExit(1)
