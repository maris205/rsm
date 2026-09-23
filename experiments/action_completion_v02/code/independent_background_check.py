#!/usr/bin/env python3
"""Independent q=dchi/dN + Radau verification; does not import primary solver.

Uses the published Lambda values, but derives all initial scalar conditions
independently. No observations fitted; this is the same dust-only toy model.
"""
from pathlib import Path
import hashlib
import json
import platform
import numpy as np
import scipy
from scipy.integrate import solve_ivp

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / 'results/background_summary.json'
MPC_KM = 3.0856775814913673e19


def main():
    source = json.loads(SOURCE.read_text())
    c = source['constants']
    om = c['omega_m']
    hi_s = c['Href_km_s_Mpc'] / MPC_KM
    hi_y = hi_s * c['year_seconds']
    ai, tstar = c['initial_a'], c['tstar_seconds']
    ni, nz = np.log(ai), -np.log1p(4.2)
    rows, checks = [], []
    for ref in source['baseline']:
        ep, ol = ref['epsilon'], ref['omega_lambda']
        # Initial dust+scalar scaling; Lambda enters the actual initial H.
        em2 = om * ai**-3 / (1.0 - .75*ep)
        taui = 2.0 / (3.0*np.sqrt(em2))
        wi = 1.0/taui
        Wi = .5*wi**2
        Ei2 = om*ai**-3 + ol + ep*(.5*wi**2+Wi)/3.0
        qi = wi/np.sqrt(Ei2)
        chii = np.log(taui/(hi_s*tstar))

        def eval_state(n, y):
            u, q, tau = y
            W = Wi*np.exp(-2*u)
            E2 = (om*np.exp(-3*n)+ol+ep*W/3)/(1-ep*q*q/6)
            return W, E2

        def rhs(n, y):
            u, q, tau = y
            W, E2 = eval_state(n, y)
            hprime_over_h = -1.5*om*np.exp(-3*n)/E2-.5*ep*q*q
            return (q, -(3+hprime_over_h)*q+2*W/E2, 1/np.sqrt(E2))

        sol = solve_ivp(rhs, (ni, 0.0), (0., qi, taui),
                        method='Radau', rtol=3e-12,
                        atol=(2e-14, 2e-14, 2e-16),
                        max_step=.03, dense_output=True)
        if not sol.success:
            raise RuntimeError(sol.message)
        u0, q0, tau0 = sol.y[:, -1]
        uz, qz, tauz = sol.sol(nz)
        W0, E02 = eval_state(0., sol.y[:, -1])
        E0 = np.sqrt(E02)
        chi0, chiz = chii+u0, chii+uz
        K0 = .5*E02*q0*q0
        ochi = ep*(K0+W0)/(3*E02)
        wchi = (K0-W0)/(K0+W0)
        # Factored differences avoid cancellation of two large chi values.
        d_inv2 = (u0-uz)*(chi0+chiz)/(chiz**2*chi0**2)
        tchi = -chi0**3/(2*q0*E0*hi_y)*d_inv2
        L0 = np.log(tau0/(hi_s*tstar))
        deltaL = np.log(tauz/tau0)
        Lz = L0+deltaL
        tlog = .5*(tau0/hi_y)*L0*deltaL*(Lz+L0)/Lz**2
        values = dict(E0=E0, chi0=chi0, w_chi0=wchi,
                      omega_chi0=ochi, age0_Gyr=tau0/hi_y/1e9,
                      T_chi_z4p2_year=tchi,
                      T_log_same_history_z4p2_year=tlog,
                      T_chi_over_same_history_log_z4p2=tchi/tlog)
        errors = {}
        for name, val in values.items():
            target = 1.0 if name == 'E0' else ref[name]
            error = abs(val-target) if name == 'chi0' else abs(val/target-1)
            bound = 1e-8 if name == 'chi0' else 1e-8
            errors[name] = float(error)
            checks.append(dict(name=f'epsilon_{ep:g}_{name}',
                               error=float(error), limit=bound,
                               metric='absolute' if name == 'chi0' else 'relative',
                               passed=bool(error < bound)))
        rows.append(dict(epsilon=ep, omega_lambda_from_primary=ol,
                         initial_tau=float(taui), initial_q=float(qi),
                         initial_E_including_Lambda=float(np.sqrt(Ei2)),
                         results={k:float(v) for k,v in values.items()},
                         errors=errors, nfev=sol.nfev, njev=sol.njev))
    result = dict(scope='Independent q=dchi/dN formulation and Radau integration; fixed primary Lambda values; same initial conditions, no radiation/EM source or observation fit.',
                  primary_summary_sha256=hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
                  independent_code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  software=dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__),
                  method='Radau', rtol=3e-12, atol=[2e-14,2e-14,2e-16],
                  cases=rows, checks=checks,
                  all_checks_passed=all(c['passed'] for c in checks))
    (ROOT/'results/independent_background_checks.json').write_text(json.dumps(result, ensure_ascii=False, indent=2)+'\n')
    print(json.dumps(dict(all_checks_passed=result['all_checks_passed'], check_count=len(checks),
                         cases=rows), indent=2))
    if not result['all_checks_passed']:
        raise SystemExit('Independent check failed.')

if __name__ == '__main__':
    main()
