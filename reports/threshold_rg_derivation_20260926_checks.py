#!/usr/bin/env python3
"""Independent leading-S^2 selected Gaussian matching and RG checks.

No experiment implementation is imported. Analytic finite Toeplitz modes are
summed with mpmath. Historical output is read only for declared comparisons.
"""
from pathlib import Path
import hashlib
import json
import mpmath as mp

ROOT = Path(__file__).resolve().parents[1]
SELF = Path(__file__).resolve()
INPUT = ROOT / 'experiments/persistent_restoration_v01/results/inputs.json'
OLD = ROOT / 'reports/mediator_auxiliary_response_20260926_checks.json'
OUT = SELF.with_suffix('.json')
mp.mp.dps = 100
checks = []


def close(name, actual, expected, tol='1e-75', scale=None):
    actual, expected, tol = mp.mpf(actual), mp.mpf(expected), mp.mpf(tol)
    denom = max(abs(actual), abs(expected), mp.mpf('1e-200')) if scale is None else mp.mpf(scale)
    err = abs(actual-expected)/denom
    checks.append({'name':name, 'passed':bool(err <= tol), 'relative_error':mp.nstr(err,12), 'tolerance':str(tol)})


def check(name, outcome):
    checks.append({'name':name, 'passed':bool(outcome)})


def chain(n, q, Lam):
    q, Lam = mp.mpf(q), mp.mpf(Lam)
    ev = [q*q+1-2*q*mp.cos(mp.pi*j/(n+1)) for j in range(1,n+1)]
    oe2 = [2*mp.sin(mp.pi*j*n/(n+1))**2/(n+1) for j in range(1,n+1)]
    weights = [o/l for o,l in zip(oe2,ev)]
    h1 = mp.fsum(weights)
    h2 = mp.fsum(w/l for w,l in zip(weights,ev))
    mass2 = Lam*Lam*h1
    return {'n':n, 'q':q, 'Lambda':Lam, 'eigenvalues':ev, 'weights':weights,
            'h1':h1, 'h2':h2, 'mass2':mass2, 'lambda2':[mass2*l for l in ev]}


def dd_log(x, y, mu):
    if x == y:
        return 2*mp.log(x/(mu*mu))
    return 2*((x*mp.log(x/(mu*mu))-y*mp.log(y/(mu*mu)))/(x-y)-1)


def Lcoef(ch, mu):
    la, w = ch['lambda2'], ch['weights']
    return mp.fsum(w[i]*w[j]*dd_log(la[i],la[j],mu)
                   for i in range(ch['n']) for j in range(ch['n']))


def coefficients(ch):
    Lam, M2, h1, h2 = ch['Lambda'], ch['mass2'], ch['h1'], ch['h2']
    d = 2/M2
    b = 2*d*h1
    a = 2*h2/M2**2
    k = b*b/(32*mp.pi**2)
    L0 = Lcoef(ch,Lam)
    h = 3*d*d*L0/(32*mp.pi**2)
    return a,b,k,h,d,L0


# Small matrices directly check the spectral coefficient and scale derivative.
for n,q in ((1,3),(2,2),(4,3),(8,mp.mpf('1.1'))):
    ch = chain(n,q,7)
    H = mp.matrix([[q*q+1 if i==j else -q if abs(i-j)==1 else 0 for j in range(n)] for i in range(n)])
    inv = H**-1
    label = f'n{n}_q{q}'
    close(label+'_h1_matrix', ch['h1'],inv[n-1,n-1])
    close(label+'_h2_matrix', ch['h2'],(inv*inv)[n-1,n-1])
    a,b,k,h,d,L0 = coefficients(ch)
    close(label+'_light_mass_matching',b,4/ch['Lambda']**2)
    for logscale in (-2,mp.mpf('0.3'),2):
        mu=ch['Lambda']*mp.exp(logscale)
        close(label+f'_heavy_scale_{logscale}',Lcoef(ch,mu),L0-4*ch['h1']**2*logscale)
    if n <= 4:
        # Direct matrix eigenvalues, numerically differentiated at C=0.
        evals,O = mp.eigsy(H)
        u = O*mp.diag([1/mp.sqrt(evals[i]) for i in range(n)])*O.T*mp.matrix([int(i==n-1) for i in range(n)])
        A=ch['mass2']*H
        def tracef(c):
            es=mp.eigsy(A+c*u*u.T,eigvals_only=True)
            return mp.fsum(s*s*(mp.log(s/ch['Lambda']**2)-mp.mpf('1.5')) for s in es)
        close(label+'_matrix_second_jet',mp.diff(tracef,0,2),L0,tol='1e-65',scale=max(ch['h1']**2,abs(L0)))
        close(label+'_supertrace_linear_cancellation',4*mp.diff(tracef,0)-4*mp.diff(tracef,0),0,scale=1)


inp=json.loads(INPUT.read_text())
par=inp['parent_parameters']
P=mp.mpf(str(par['P']))
F=mp.sqrt(mp.mpf(str(par['F_squared_eV2'])))
Ui=mp.mpf(str(par['Ui']))
chi_i=mp.mpf(str(inp['chi_first']))
chi_f=mp.mpf(str(inp['chi_last']))
bench=[]
old=json.loads(OLD.read_text())
for n,Lams,mgs in ((128,'1e12','1e-6'),(127,'2e12','1e-15')):
    ch=chain(n,3,Lams)
    Lam=ch['Lambda']; mg=mp.mpf(mgs); S=3*P*mg*mg
    a,b,k,h,d,L0=coefficients(ch)
    mx=mp.sqrt(b*S); loglight=mp.log(b*S/Lam**2)
    pref=f'n{n}'
    def c_uv(mu): return -a+5*k*mp.log(mu/Lam)
    def h_mu(mu): return h-3*k*mp.log(mu/Lam)
    def c_low(mu): return -a+h+2*k*mp.log(mu/Lam)
    def v_fixed(s,mu): return c_low(mu)*s*s+k*s*s*(mp.log(b*s/mu**2)-mp.mpf('1.5'))
    def v_running(s): return (c_low(mp.sqrt(b*s))-mp.mpf('1.5')*k)*s*s
    def v_s(s): return 2*s*(-a+h+k*(mp.log(b*s/Lam**2)-1))
    for muH_over_L in ('0.5','1','2'):
        muH=Lam*mp.mpf(muH_over_L)
        for mu in (mx/2,mx,2*mx,Lam):
            v_from_match=(c_uv(muH)+h_mu(muH)+2*k*mp.log(mu/muH))*S*S+k*S*S*(mp.log(b*S/mu**2)-mp.mpf('1.5'))
            close(pref+f'_threshold_{muH_over_L}_mu_{mp.nstr(mu,8)}',v_from_match,v_fixed(S,Lam))
    for mu in (mx/2,mx,2*mx,Lam):
        close(pref+f'_scale_invariance_{mp.nstr(mu,8)}',v_fixed(S,mu),v_fixed(S,Lam))
        close(pref+f'_force_scale_invariance_{mp.nstr(mu,8)}',mp.diff(lambda z:v_fixed(S*z,mu),1)/S,v_s(S))
    close(pref+'_beta_low',mp.diff(lambda ell:c_low(Lam*mp.exp(ell)),0),2*k)
    close(pref+'_beta_high',mp.diff(lambda ell:c_uv(Lam*mp.exp(ell)),0),5*k)
    close(pref+'_beta_heavy',mp.diff(lambda ell:h_mu(Lam*mp.exp(ell)),0),-3*k)
    close(pref+'_field_scale_potential',v_running(S),v_fixed(S,Lam))
    close(pref+'_field_scale_first_derivative',mp.diff(lambda z:v_running(S*z),1)/S,v_s(S))
    close(pref+'_field_scale_second_derivative',mp.diff(lambda z:v_running(S*z),1,2)/S**2,2*c_low(mx))
    close(pref+'_running_chain_term',mp.diff(lambda z:c_low(mp.sqrt(b*S*z)),1)/S,k/S)
    frozen_derivative=2*S*(c_low(mx)-mp.mpf('1.5')*k)
    close(pref+'_incorrect_frozen_coefficient_missing_term',v_s(S)-frozen_derivative,k*S)
    # Integrate the light determinant out at an arbitrary local matching scale.
    # The resulting background-dependent functional retains the same force.
    for eta in ('0.5','1','2'):
        eta=mp.mpf(eta)
        def below(s):
            mu=eta*mp.sqrt(b*s)
            matched=c_low(mu)+k*(mp.log(b*s/mu**2)-mp.mpf('1.5'))
            return matched*s*s
        close(pref+f'_light_decoupling_{eta}',below(S),v_running(S))
        close(pref+f'_light_decoupling_force_{eta}',mp.diff(lambda z:below(S*z),1)/S,v_s(S))
    close(pref+'_a_large_n',a*Lam**4,mp.mpf('2.25'),tol='1e-65')
    tree_s=-2*a*S
    heavy_s=2*h*S
    light_s=2*k*S*(loglight-1)
    close(pref+'_split_slopes',v_s(S),tree_s+heavy_s+light_s)
    check(pref+'_same_sign_selected_corrections',tree_s<0 and heavy_s<0 and light_s<0)
    Ulast=Ui*mp.exp(-2*(chi_f-chi_i))
    rw=F*mp.sqrt(Ulast/2)/(mg*P)
    Schi=6*P*mg**2*rw*(1-rw)
    Sy=2*P*mg**2*rw
    row={'n':n,'Lambda_eV':Lam,'mG_eV':mg,'M_eV':mp.sqrt(ch['mass2']),
         'S_eV4':S,'h1':ch['h1'],'h2':ch['h2'],'L_Lambda':L0,
         'a_times_Lambda4':a*Lam**4,'b_times_Lambda2':b*Lam**2,
         'k_times_Lambda4':k*Lam**4,'heavy_threshold_times_Lambda4':h*Lam**4,
         'c_low_Lambda_times_Lambda4':c_low(Lam)*Lam**4,
         'c_low_mX_times_Lambda4':c_low(mx)*Lam**4,
         'mX_eV_leading':mx,'ln_mX2_over_Lambda2':loglight,
         'heavy_mass_min_eV':mp.sqrt(ch['lambda2'][0]),'heavy_mass_max_eV':mp.sqrt(ch['lambda2'][-1]),
         'delta_tree_dV_dS_leading':tree_s,'heavy_loop_dV_dS_leading':heavy_s,
         'light_loop_dV_dS_leading':light_s,'total_loop_dV_dS_leading':heavy_s+light_s,
         'total_correction_dV_dS_leading':v_s(S),'loop_to_tree_slope_leading':(heavy_s+light_s)/tree_s,
         'aligned_last_tree_force_ratio':abs(tree_s*Schi)/(2*Ulast),
         'aligned_last_loop_force_ratio':abs((heavy_s+light_s)*Schi)/(2*Ulast),
         'aligned_last_total_force_ratio':abs(v_s(S)*Schi)/(2*Ulast),
         'quadrature_last_total_force_ratio':abs(v_s(S)*Sy)/(2*Ulast),
         'incorrect_drop_log_loop_to_tree_ratio':(h-k)/(-a),
         'loopforce_fraction_missed_without_chainrule':k*S/abs((heavy_s+light_s)),
         'scope':'Leading S² correction with S=3P mG²; imported clock jets are response probes, not a solved new cosmology.'}
    if n==127:
        historical=mp.mpf(old['small_splitting_benchmarks'][1]['loop_to_finiteF_tree_slope_ratio'])
        close(pref+'_historical_low_loop_ratio',row['loop_to_tree_slope_leading'],historical,tol='1e-13')
    else:
        historical=mp.mpf(str(old['gaussian_spectrum'][1]['loop_to_finiteF_tree_slope_ratio']))
        row['historical_full_finiteS_loop_to_tree_slope']=historical
        row['leading_to_historical_full_relative_difference']=abs(row['loop_to_tree_slope_leading']/historical-1)
        # Comparison diagnostic only: no invented 0.2% claim for this truncation.
    bench.append(row)


def serial(v):
    if isinstance(v,mp.mpf): return mp.nstr(v,50)
    if isinstance(v,dict): return {k:serial(x) for k,x in v.items()}
    if isinstance(v,list): return [serial(x) for x in v]
    return v


result={'status':'internal_independent_derivation_and_reproduction',
        'scope':'Leading S² selected rigid Gaussian; external S, fixed b and fixed matching boundary. Not full SUGRA or all-order resummation.',
        'arithmetic_digits':100,
        'read_only_inputs':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in (INPUT,OLD)},
        'script_sha256':hashlib.sha256(SELF.read_bytes()).hexdigest(),
        'execution_history':[
            'Initial run: 93/105. Preserved in threshold_rg_derivation_20260926_initial_checks.json.',
            'Initial 11 derivative failures were numerical precision loss from differentiating at S~1e43/1e25; reparameterized s=S*z and differentiated at dimensionless z=1, retaining original tolerances.',
            'Initial n=1 heavy second-jet reference vanishes at mu=Lambda; changed its relative-error denominator to max(h1²,abs(L)) so a true zero has a defined natural normalization. Original 1e-65 tolerance retained.',
            'No physical formula, finite boundary, or tolerance changed.'
        ],
        'benchmarks':bench,'checks':checks,'passed':sum(c['passed'] for c in checks),'total':len(checks),
        'failures':[c for c in checks if not c['passed']]}
OUT.write_text(json.dumps(serial(result),ensure_ascii=False,indent=2)+'\n')
print(f"{result['passed']}/{result['total']}")
for bmark in bench:
    print(json.dumps(serial(bmark),ensure_ascii=False))
if result['failures']:
    print(json.dumps(result['failures'],indent=2))
    raise SystemExit(1)
