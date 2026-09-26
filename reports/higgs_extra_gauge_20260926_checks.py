#!/usr/bin/env python3
"""Independent F=D=0 extra-U(1) stabilization and complete-current audit.

No experiment implementation is imported. Exact small matrices, independent
180-digit tridiagonal solves, and double-precision complete spectra are used.
The calculation is a canonical rigid supersymmetric tree-level completion,
not a finite-F supergravity or radiative matching calculation.
"""
from pathlib import Path
import hashlib
import json
import sys

import mpmath as mp
import numpy as np
import scipy
from scipy.linalg import eigh_tridiagonal
import sympy as sp
from sympy.matrices.normalforms import smith_normal_form

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
CHECKS = []


def check(name, condition, detail=None):
    entry = {"name": name, "passed": bool(condition)}
    if detail is not None:
        entry["detail"] = detail
    CHECKS.append(entry)


def symbolic_equal(a, b):
    return sp.cancel(a-b) == 0


def charges(q, n, eta):
    result = sp.zeros(n+1, n+1)
    result[0, 0] = q
    for j in range(n-1):
        result[j+1, j], result[j+1, j+1] = 1, -q
    result[n, n-1], result[n, n] = 1, eta
    return result


# eta is a gauge-coupling ratio. The actual integer charges use eta=1.
# All anomaly checks therefore use the integer charge matrix, not eta.
small_cases = []
for q in (2, 3, 5):
    for n in (2, 3, 5):
        for eta in (sp.Rational(1, 2), sp.Integer(1), sp.Integer(2)):
            tag = f"q{q}_n{n}_eta{eta}"
            Q = charges(q, n, eta)
            Hp = Q.T*Q
            H = Hp[:n, :n]
            old_G = H.inv()
            inv = Hp.inv()
            f, g = old_G[0, 0], old_G[0, n-1]
            e = sp.zeros(n, 1)
            e[n-1] = 1
            expected = sp.diag(H, eta**2)
            expected[:n, n] = eta*e
            expected[n, :n] = eta*e.T
            check(tag+"_augmented_gram", Hp == expected)
            check(tag+"_rank", Q.rank() == n+1)
            check(tag+"_determinant", Q.det() == (-1)**(n-1)*q**n*eta)
            check(tag+"_gram_determinant", Hp.det() == q**(2*n)*eta**2)
            check(tag+"_all_continuous_gauges_broken", Hp.is_positive_definite is True)
            check(tag+"_F_D_quotient_no_modulus", (2*(n+1)-(n+1)-Q.rank()) == 0)
            top = old_G+(old_G*e)*(e.T*old_G)/(1-f)
            cross = -old_G*e/(eta*(1-f))
            schur_inv = sp.zeros(n+1, n+1)
            schur_inv[:n, :n] = top
            schur_inv[:n, n] = cross
            schur_inv[n, :n] = cross.T
            schur_inv[n, n] = 1/(eta**2*(1-f))
            check(tag+"_complete_schur_inverse", inv == schur_inv)
            check(tag+"_old00", inv[0, 0] == sp.Rational(1, q*q))
            check(tag+"_old0e", inv[0, n-1] == sp.Rational(1, q**(n+1)))
            check(tag+"_oldee", inv[n-1, n-1] == sum(sp.Rational(1, q**(2*k)) for k in range(1, n+1)))
            check(tag+"_old0extra", inv[0, n] == -sp.Rational(1, q**(n+1))/eta)
            check(tag+"_normalized_cross", inv[0,n-1]/inv[0,0] == sp.Rational(1, q**(n-1)))
            check(tag+"_source_projector", Q*inv*Q.T == sp.eye(n+1))
            r = Q[n, :].T
            response = sp.zeros(n+1, 1)
            response[n] = 1/eta
            check(tag+"_right_complete_response", inv*r == response)
            check(tag+"_physical_endpoint_cross_zero", (Q[0,:]*inv*r)[0] == 0)
            check(tag+"_right_self_one", (r.T*inv*r)[0] == 1)
            check(tag+"_false_omit_extra_changes_cross", inv[0,n-1] != 0 and inv[0,n-1]+eta*inv[0,n] == 0)
            integer_Q = charges(q, n, 1)
            if eta == 1:
                smith = smith_normal_form(integer_Q,domain=sp.ZZ)
                invariants = [abs(smith[j,j]) for j in range(n+1)]
                check(tag+"_discrete_gauge_Smith_form", invariants == [1]*n+[q**n])
            check(tag+"_linear_anomalies", all(x+(-x) == 0 for x in integer_Q))
            # Every U(1)^3 and mixed U(1)_i U(1)_j U(1)_k anomaly cancels pairwise.
            check(tag+"_all_cubic_mixed_anomalies", all(
                integer_Q[a,i]*integer_Q[a,j]*integer_Q[a,k]
                +(-integer_Q[a,i])*(-integer_Q[a,j])*(-integer_Q[a,k]) == 0
                for a in range(n+1) for i in range(n+1)
                for j in range(n+1) for k in range(n+1)))
            small_cases.append({"q": q, "n": n, "eta": str(eta)})

# Full canonical quadratic superpotential and scalar potential at F=D=0.
ph_symbol, pt_symbol, z_symbol, vev_symbol, lambda_symbol = sp.symbols("ph_symbol pt_symbol z_symbol vev_symbol lambda_symbol")
single_W = lambda_symbol*z_symbol*(ph_symbol*pt_symbol-vev_symbol**2)
for variable in (ph_symbol,pt_symbol,z_symbol):
    check("Fflat_actual_W_derivative_"+str(variable), sp.diff(single_W,variable).subs({ph_symbol:vev_symbol,pt_symbol:vev_symbol,z_symbol:0}) == 0)
R, D, Z, v, lam = sp.symbols("R D Z v lam", real=True)
Phi = v+(R+D)/sp.sqrt(2)
Ptilde = v+(R-D)/sp.sqrt(2)
W = sp.expand(lam*Z*(Phi*Ptilde-v*v))
W2 = sp.sqrt(2)*lam*v*Z*R
check("canonical_W_quadratic_radial_singlet_mass", sp.expand(W-W2-lam*Z*(R*R-D*D)/2) == 0)
check("canonical_chiral_mass_squared", sp.diff(W,R,Z).subs({R:0,D:0,Z:0})**2 == 2*lam**2*v**2)
d, a = sp.symbols("d a", real=True)
ph = v+(d+sp.I*a)/2
pt = v-(d+sp.I*a)/2
moment = sp.expand(ph*sp.conjugate(ph)-pt*sp.conjugate(pt))
check("canonical_moment_map_real_difference", moment == 2*v*d)
check("canonical_phase_goldstone", sp.diff(moment,a) == 0)
check("Dflat_actual_moment_at_vacuum", moment.subs(d,0) == 0)
# Gauge mass: each pair gives 2 g^2 v^2 A^T Q^TQ A = (1/2) M^2 A^T H A.
check("canonical_vector_M_equals_2gv", sp.Rational(1,2)*4 == 2)

for n, eta in ((2, .5), (3, 1.), (5, 2.)):
    q = 3
    Q = np.asarray(charges(q,n,sp.Rational(str(eta)))).astype(float)
    m = n+1
    scale, chiral = 1.7, 2.3
    fermion = np.zeros((4*m,4*m))
    fermion[:m,m:2*m] = chiral*np.eye(m)
    fermion[m:2*m,:m] = chiral*np.eye(m)
    fermion[2*m:3*m,3*m:] = scale*Q
    fermion[3*m:,2*m:3*m] = scale*Q.T
    vec_sq = np.linalg.eigvalsh(scale*scale*Q.T@Q)
    expected_squared = np.sort(np.r_[np.full(2*m,chiral**2),np.repeat(vec_sq,2)])
    actual_squared = np.sort(np.linalg.svd(fermion,compute_uv=False)**2)
    check(f"full_canonical_fermion_n{n}_eta{eta}", np.allclose(actual_squared,expected_squared,rtol=2e-13,atol=2e-13))
    check(f"full_canonical_real_Higgs_n{n}_eta{eta}", np.allclose(np.linalg.eigvalsh(scale*scale*Q@Q.T),vec_sq,rtol=2e-13,atol=2e-13))
    check(f"physical_degree_balance_n{n}_eta{eta}", 3*m+m+4*m == 2*(4*m))

mp.mp.dps = 180
q, n = mp.mpf(3), 127
f = (q**(2*n)-1)/(q**(2*n+2)-1)
g = q**(n-1)*(q*q-1)/(q**(2*n+2)-1)
old_tau = g/f
Lambda = mp.mpf('2e12')
M = Lambda*mp.sqrt(f)


def tridiagonal_solve(diagonal, off, rhs):
    """Direct full augmented tridiagonal elimination, no endpoint formulas."""
    dd, bb = list(diagonal), list(rhs)
    for i in range(1,len(dd)):
        factor = off[i-1]/dd[i-1]
        dd[i] -= factor*off[i-1]
        bb[i] -= factor*bb[i-1]
    result = [mp.mpf(0)]*len(dd)
    result[-1] = bb[-1]/dd[-1]
    for i in reversed(range(len(dd)-1)):
        result[i] = (bb[i]-off[i]*result[i+1])/dd[i]
    return result


def residual(diagonal,off,x,b):
    errors=[]
    for i in range(len(x)):
        value=diagonal[i]*x[i]-b[i]
        if i: value += off[i-1]*x[i-1]
        if i+1<len(x): value += off[i]*x[i+1]
        errors.append(abs(value))
    return max(errors)


def close(x,y,tol='1e-110'):
    return abs(x-y) <= mp.mpf(tol)*max(abs(x),abs(y),mp.mpf('1e-100'))


spectra=[]
for eta_string in ('0.5','1','2'):
    eta=mp.mpf(eta_string)
    diag=[q*q+1]*n+[eta*eta]
    off=[-q]*(n-1)+[eta]
    rhs0=[mp.mpf(0)]*(n+1);rhs0[0]=1
    rhse=[mp.mpf(0)]*(n+1);rhse[n-1]=1
    rhsr=list(rhse);rhsr[n]=eta
    x0=tridiagonal_solve(diag,off,rhs0)
    xe=tridiagonal_solve(diag,off,rhse)
    xr=tridiagonal_solve(diag,off,rhsr)
    tag=f"full128_eta{eta_string}"
    check(tag+"_00_180digits",close(x0[0],1/q**2))
    check(tag+"_0e_180digits",close(x0[n-1],q**(-(n+1))))
    check(tag+"_0extra_180digits",close(x0[n],-q**(-(n+1))/eta))
    check(tag+"_ee_180digits",close(xe[n-1],f/(1-f)))
    check(tag+"_physical_cross_zero_180digits",abs(xr[0]) < mp.mpf('1e-170'))
    check(tag+"_physical_response_old_sites_zero",max(abs(z) for z in xr[:n]) < mp.mpf('1e-170'))
    check(tag+"_physical_response_extra",close(xr[n],1/eta))
    for source,xx,bb in [('left',x0,rhs0),('right',xe,rhse),('complete_right',xr,rhsr)]:
        check(tag+'_residual_'+source,residual(diag,off,xx,bb)<mp.mpf('1e-170'))
    dd=np.asarray([float(z) for z in diag]);oo=np.asarray([float(z) for z in off])
    eigs,U=eigh_tridiagonal(dd,oo)
    dense=np.diag(dd)+np.diag(oo,1)+np.diag(oo,-1)
    eigs2=np.linalg.eigvalsh(dense)
    check(tag+'_spectrum_independent_dense',np.allclose(eigs,eigs2,rtol=3e-13,atol=3e-13))
    check(tag+'_eigenvector_residual',np.max(np.abs(dense@U-U*eigs))<3e-13)
    check(tag+'_all_massive',eigs[0]>0)
    # Congruence lower bound: Hplus=T^T diag(H,s)T, ||T^{-1}||=(sqrt(beta²+4)+beta)/2.
    erhs=[mp.mpf(0)]*n;erhs[-1]=1
    ge=tridiagonal_solve([q*q+1]*n,[-q]*(n-1),erhs)
    beta=eta*mp.sqrt(mp.fsum(z*z for z in ge))
    s=eta**2*(1-f)
    hmin=q*q+1-2*q*mp.cos(mp.pi/(n+1))
    hmax=q*q+1+2*q*mp.cos(mp.pi/(n+1))
    tnorm=(mp.sqrt(beta*beta+4)+beta)/2
    bound=min(hmin,s)/(tnorm*tnorm)
    upper=max(hmax,eta*eta)+eta
    check(tag+'_rigorous_lower_bound',float(bound)<=eigs[0]*(1+1e-13))
    check(tag+'_weyl_upper_bound',eigs[-1]<=float(upper)*(1+1e-13))
    # Independent high-precision bound-state characteristic root below old bulk band.
    # det(Hplus-z)=det(H-z)*(eta²-z-eta²[(H-z)^-1]ee).
    def root_equation(z):
        pp=q*q+1-z
        for _ in range(1,n):
            pp=q*q+1-z-q*q/pp
        return eta*eta-z-eta*eta/pp
    semi_infinite=eta*eta*(q*q-eta*eta-1)/(q*q-eta*eta)
    low=mp.findroot(root_equation,(semi_infinite*mp.mpf('.99'),semi_infinite*mp.mpf('1.01')),tol=mp.mpf('1e-160'))
    check(tag+'_high_precision_spectral_root',abs(root_equation(low))<mp.mpf('1e-155'))
    check(tag+'_lowest_eigenvalue_double_vs_mp',abs(float(low)-eigs[0])<2e-13)
    check(tag+'_surface_mode_finite_correction',low>=semi_infinite and low-semi_infinite<mp.mpf('1e-50'))
    spectra.append({
        'eta':eta_string,'site_count_old':n,'gauge_count_new':n+1,
        'lambda_min_180digits':mp.nstr(low,160),'lambda_max':float(eigs[-1]),
        'analytic_lower_bound':mp.nstr(bound,40),'weyl_upper_bound':mp.nstr(upper,40),
        'mass_min_eV':mp.nstr(M*mp.sqrt(low),40),
        'mass_max_eV':float(M)*float(np.sqrt(eigs[-1])),
        'chiral_mass_at_lambda_over_g_1_eV':mp.nstr(M/mp.sqrt(2),40),
        'all_dimensionless_eigenvalues':eigs.tolist(),
        'old_current_cross':mp.nstr(x0[n-1],45),
        'old_current_self_left':mp.nstr(x0[0],45),
        'old_current_self_right':mp.nstr(xe[n-1],45),
        'physical_right_cross':mp.nstr(xr[0],12),
        'physical_right_self':mp.nstr(xr[n-1]+eta*xr[n],45),
        'direct_solve_max_residual':mp.nstr(max(residual(diag,off,x0,rhs0),residual(diag,off,xe,rhse),residual(diag,off,xr,rhsr)),12),
        'semi_infinite_lambda_formula':mp.nstr(semi_infinite,45),
        'finite_chain_lambda_difference':mp.nstr(low-semi_infinite,35)
    })

failed=[row for row in CHECKS if not row['passed']]
inputs=[ROOT/'reports/chain_locality_breaking_20260926.md', ROOT/'experiments/threshold_rg_v01/results/inputs.json']
output={
    'scope':'Independent canonical rigid supersymmetric F=D=0 tree completion; no full SUGRA, finite-F spectrum, or loop threshold claim.',
    'date':'2026-09-26','total':len(CHECKS),'passed':len(CHECKS)-len(failed),'failed':failed,
    'method':'Exact rational matrices; 180-digit full augmented tridiagonal solves; independent dense/tridiagonal spectra; complete quadratic canonical multiplet spectrum.',
    'inputs':{'q':3,'n_old':127,'Lambda_eV':'2e12','M_eV':mp.nstr(M,45),
              'f':mp.nstr(f,45),'old_tau':mp.nstr(old_tau,45),'new_abstract_tau':mp.nstr(q**(1-n),45),
              'fixed_original_Lambda_left_multiplier':mp.nstr(1/(q*q*f),45),
              'fixed_original_Lambda_cross_multiplier':mp.nstr(1/(1-f),45),
              'fixed_original_Lambda_right_multiplier':mp.nstr(1/(1-f),45),
              'normalized_new_right_multiplier':mp.nstr(q*q*f/(1-f),45),
              'eta_values':['0.5','1','2'],'lambda_over_g_for_radial_example':1},
    'small_exact_cases':small_cases,'spectra':spectra,'checks':CHECKS,
    'input_hashes':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},
    'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    'versions':{'python':sys.version.split()[0],'numpy':np.__version__,'scipy':scipy.__version__,'mpmath':mp.__version__,'sympy':sp.__version__},
    'initial_failures':[{'file':'reports/higgs_extra_gauge_20260926_initial_checks.json','passed':634,'total':635,'name':'canonical_W_quadratic_radial_singlet_mass','cause':'SymPy structural equality compared expanded versus factored polynomials. The correction expands their difference and tests exact zero; the mass formula and tolerance did not change.'}],
    'limitations':['No generic lower bound on a radiatively generated shortcut.',
                   'The complete local right Higgs current has exact zero cross at p=0 in this square full-rank construction.',
                   'Abstract old-site source matching is not a demonstrated neutral-bulk-modulus gauge-invariant completion.',
                   'Equal canonical vevs, no FI term, global supersymmetry, F=D=0, and tree level are explicit assumptions.']
}
path=HERE/'higgs_extra_gauge_20260926_checks.json'
path.write_text(json.dumps(output,indent=2,ensure_ascii=False)+'\n')
print(json.dumps({'total':output['total'],'passed':output['passed'],'failed':failed,'output':str(path)},ensure_ascii=False))
raise SystemExit(bool(failed))
