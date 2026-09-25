#!/usr/bin/env python3
"""Frozen finite-window charged-threshold experiment; no observational fit."""
from pathlib import Path
import csv
import hashlib
import importlib.util
import json
import platform
import time
from datetime import datetime, timezone

import numpy as np
import scipy
from scipy.integrate import solve_ivp

ROOT = Path(__file__).resolve().parents[1]
REPO = ROOT.parents[1]
OUT = ROOT / 'results'
EPS = 1e-4
MPL = 2.435e27
HBAR = 6.582119569e-16
ALPHA = 1 / 137.035999084
C = ALPHA / (12 * np.pi)
HSEC = 67.4 / 3.0856775814913673e19
HEV = HSEC * HBAR
YEAR = 365.25 * 86400
FSQ = EPS * MPL**2
USCALE = FSQ * HEV**2
NGRID = np.linspace(-np.log(5.2), 0, 1025)
LABELS = [(k, ell) for k in ('chi', 'R') for ell in (0., 1., 100., 1e4)]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def dump(path, obj):
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2, allow_nan=False) + '\n')


def name(k, ell):
    return f'{k}_probe' if ell == 0 else f'{k}_loop_{ell:g}'


def get_inputs():
    path = REPO / 'experiments/cosmic_bridge_v01/code/background.py'
    spec = importlib.util.spec_from_file_location('prior_background', path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    reference = mod.Background(epsilon=EPS)
    first = reference.evaluate_n(NGRID[0])
    inputs = {
        'epsilon': EPS, 'Mpl_eV': MPL, 'hbar_eV_s': HBAR, 'H_ref_s_inv': HSEC,
        'H_ref_eV': HEV, 'F_squared_eV2': FSQ, 'potential_unit_eV4': USCALE,
        'alpha_reference': ALPHA, 'b_em': 1/3,
        'omega_m': mod.OMEGA_M0, 'omega_r': mod.OMEGA_R0,
        'omega_lambda': reference.omega_lambda, 'tstar_s': mod.TSTAR_RAD_S,
        'N_i': float(NGRID[0]), 'chi_i': first['chi'], 'q_i': first['q'],
        'tau_i': first['tau'], 'W_i': first['W'], 'xi': first['chi']**2,
        'reference_code_sha256': sha(path), 'protocol_sha256': sha(ROOT/'protocol.md'),
        'case_definitions': [{'case': name(k, ell), 'kinetic': k, 'ell': ell}
                             for k, ell in LABELS],
    }
    return inputs, reference


class Model:
    def __init__(self, inputs, kinetic, ell):
        self.i = inputs
        self.kinetic = kinetic
        self.ell = ell
        self.mass_i = (32*np.pi**2*USCALE*ell)**0.25 if ell else None

    def field(self, u):
        chi = self.i['chi_i'] + u
        logratio = np.log1p(u / self.i['chi_i'])
        gm1 = 0.5*np.expm1(-2*logratio)
        g = 1 + gm1
        logg = np.log1p(gm1)
        gp = -self.i['xi']/chi**3
        gpp = 3*self.i['xi']/chi**4
        K = np.ones_like(np.asarray(u)) if self.kinetic == 'chi' else np.exp(2*u)
        Kp = 0*K if self.kinetic == 'chi' else 2*K
        prefactor = self.mass_i**2 / (96*np.pi**2*FSQ) if self.mass_i else 0.
        dA = prefactor * gp**2 / g
        dAp = prefactor*(2*gp*gpp/g - gp**3/g**2)
        A, Ap = K + dA, Kp + dAp
        tree = self.i['W_i']*np.exp(-2*u)
        loop = self.ell*(g**2*logg - 1.5*np.expm1(2*logg))
        loopp = 2*self.ell*g*(logg-1)*gp
        V, Vp = tree + loop, -2*tree + loopp
        return {'chi':chi, 'g':g, 'gp':gp, 'logg':logg, 'K':K, 'dA':dA,
                'A':A, 'Ap':Ap, 'tree':tree, 'loop':loop, 'loopp':loopp,
                'V':V, 'Vp':Vp}

    def terms(self, n, u, q):
        r = self.field(u)
        rad = self.i['omega_r']*np.exp(-4*n)
        dust = self.i['omega_m']*np.exp(-3*n)
        denominator = 1-EPS*r['A']*q*q/6
        numerator = rad+dust+self.i['omega_lambda']+EPS*r['V']/3
        e2 = numerator/denominator
        if np.any(e2 <= 0):
            raise ValueError('Friedmann E^2 not positive in solver trial')
        f = -(4*rad+3*dust)/(2*e2)-EPS*r['A']*q*q/2
        r.update(rad=rad, dust=dust, denominator=denominator,
                 numerator=numerator, E2=e2, E=np.sqrt(e2), f=f)
        return r

    def rhs(self, n, state):
        u, q, tau, work = state
        r = self.terms(n, u, q)
        qp = -(3+r['f'])*q-0.5*r['Ap']/r['A']*q*q-r['Vp']/(r['A']*r['E2'])
        return q, qp, 1/r['E'], 3*r['A']*r['E2']*q*q

    def integrate(self, tight=False):
        def stop_offset(n,y): return 20-abs(y[0])
        def stop_log(n,y): return 0.25-abs(self.field(y[0])['logg'])
        def stop_denom(n,y): return self.terms(n,y[0],y[1])['denominator']-1e-4
        def stop_numerator(n,y): return self.terms(n,y[0],y[1])['numerator']-1e-10
        events = (stop_offset, stop_log, stop_denom, stop_numerator)
        for event in events:
            event.terminal, event.direction = True, -1
        solution = solve_ivp(self.rhs, (NGRID[0], 0),
                            (0, self.i['q_i'], self.i['tau_i'], 0),
                            method='DOP853', rtol=2e-12 if tight else 2e-10,
                            atol=2e-14 if tight else 2e-12,
                            max_step=0.005 if tight else 0.01,
                            dense_output=True, events=events)
        return solution

    def evaluate(self, solution, n):
        u,q,tau,work = solution.sol(n)
        r = self.terms(n,u,q)
        kinetic = 0.5*r['A']*r['E2']*q*q
        rho = kinetic+r['V']
        first = self.terms(NGRID[0], 0., self.i['q_i'])
        rho_i = 0.5*first['A']*first['E2']*self.i['q_i']**2+first['V']
        lnt = np.log(tau/(HSEC*self.i['tstar_s']))
        complete = abs(solution.t[-1])<1e-12
        r.update(u=u, q=q, tau=tau, N=n, a=np.exp(n), z=np.expm1(-n),
                 t_Gyr=tau/HSEC/(YEAR*1e9), kinetic=kinetic, rho=rho,
                 work=work, Omega_chi=EPS*rho/(3*r['E2']),
                 energy_residual=(rho+work-rho_i)/max(abs(rho_i),1.), lnt=lnt,
                 friedmann_residual=(r['E2']-r['rad']-r['dust']-self.i['omega_lambda']
                                     -EPS*rho/3)/r['E2'])
        if complete:
            today = self.field(solution.sol(0)[0])
            log_mass2_ratio = r['logg']-today['logg']
            d = C*log_mass2_ratio
            delta = d/(1-d)
            s = np.expm1(-2*np.log(r['chi']/today['chi']))/today['chi']**2
            st = np.expm1(-2*np.log(lnt/lnt[-1]))/lnt[-1]**2
            gamma = C*self.i['xi']/(1+self.i['xi']/today['chi']**2)
            r.update(delta_alpha=delta, delta_alpha_linear=gamma*s, s_chi=s, s_age=st,
                     gauge_denominator=1-d)
        if self.mass_i:
            mass = self.mass_i*np.sqrt(r['g'])
            r.update(mass_eV=mass, H_over_mass=HEV*r['E']/mass,
                     mass_adiabatic=np.abs(0.5*r['gp']/r['g']*q*r['E']*HEV/mass),
                     kinetic_loop_fraction=r['dA']/r['K'])
        return r


def metrics(model, solution, values):
    r = values
    result = {
        'case':name(model.kinetic,model.ell), 'kinetic':model.kinetic, 'ell':model.ell,
        'mass_initial_eV':model.mass_i, 'success':bool(solution.success),
        'complete':bool(abs(solution.t[-1])<1e-12), 'message':solution.message,
        'N_last':float(solution.t[-1]), 'rhs_evaluations':solution.nfev,
        'events':[list(x) for x in solution.t_events],
        'chi_change':float(r['u'][-1]), 'q_today_or_last':float(r['q'][-1]),
        'q_min':float(np.min(r['q'])), 'clock_reverses':bool(np.any(r['q']<0)),
        'H_today_over_Href_or_last':float(r['E'][-1]),
        'age_today_Gyr_or_last':float(r['t_Gyr'][-1]),
        'max_abs_Omega_chi':float(np.max(abs(r['Omega_chi']))),
        'max_energy_residual':float(np.max(abs(r['energy_residual']))),
        'max_friedmann_residual':float(np.max(abs(r['friedmann_residual']))),
        'min_kinetic_denominator':float(np.min(r['denominator'])),
        'min_A':float(np.min(r['A'])), 'min_mass_squared_ratio':float(np.min(r['g'])),
        'max_loop_force_over_tree':float(np.max(abs(r['loopp'])/(2*r['tree']))),
    }
    if result['complete']:
        amp = max(float(np.max(abs(r['delta_alpha']))),1e-30)
        result.update(delta_alpha_zi_ppm=float(r['delta_alpha'][0]*1e6),
                      delta_alpha_max_abs_ppm=amp*1e6,
                      threshold_linear_error_over_peak=float(np.max(abs(r['delta_alpha']-
                                                                 r['delta_alpha_linear']))/amp),
                      chi_age_shape_error=float(np.max(abs(r['s_chi']/r['s_chi'][0]-
                                                             r['s_age']/r['s_age'][0]))),
                      alpha_age_shape_error=float(np.max(abs(r['delta_alpha']/r['delta_alpha'][0]-
                                                               r['s_age']/r['s_age'][0]))),
                      alpha_drift_today_per_year=float(C*r['gp'][-1]/r['g'][-1]*
                                                       r['q'][-1]*r['E'][-1]*HSEC*YEAR),
                      min_gauge_denominator=float(np.min(r['gauge_denominator'])))
    if model.mass_i:
        for key in ('H_over_mass','mass_adiabatic','kinetic_loop_fraction'):
            result['max_'+key] = float(np.max(r[key]))
    return result


def budget(inputs, baseline):
    rho_crit = 3*MPL**2*HEV**2
    u = baseline['u']
    records=[]
    for mass in (1.,1e6,1e11):
        ell=mass**4/(32*np.pi**2*USCALE)
        model=Model(inputs,'chi',ell)
        r=model.field(u)
        records.append({
            'mass_initial_eV':mass, 'ell':ell,
            'max_abs_CW_change_over_rho_crit_ref':float(np.max(abs(r['loop']))*USCALE/rho_crit),
            'max_loop_force_over_tree':float(np.max(abs(r['loopp'])/(2*r['tree']))),
            'max_kinetic_loop_fraction':float(np.max(r['dA'])),
            'ppm_budget_over_rho_crit_ref':float(3*mass**4*1e-6/(4*np.pi*ALPHA*rho_crit)),
        })
    return {'scope':'Unprotected MSbar at fixed mu=m_i, initial constant removed; budgets on tree chi path only.',
            'rho_crit_ref_eV4':rho_crit, 'records':records,
            'mass_for_1ppm_CW_equals_rho_crit_eV':float((rho_crit*4*np.pi*ALPHA/(3e-6))**0.25)}


def write_csv(path, rows):
    keys=list(dict.fromkeys(k for row in rows for k in row))
    with path.open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=keys)
        writer.writeheader(); writer.writerows(rows)


def main():
    started=time.time()
    OUT.mkdir(parents=True,exist_ok=True)
    inputs, reference=get_inputs()
    dump(OUT/'inputs.json',inputs)
    checks=[]; summaries=[]; rows=[]; trajectories={}
    def check(label,value,tolerance=1e-7):
        checks.append({'name':label,'value':float(value),'tolerance':tolerance,
                       'passed':bool(np.isfinite(value) and value<tolerance)})
    for k,ell in LABELS:
        label=name(k,ell)
        try:
            model=Model(inputs,k,ell)
            solution=model.integrate()
            n=NGRID[NGRID<=solution.t[-1]+1e-13]
            r=model.evaluate(solution,n)
            trajectories[label]=r
            sm=metrics(model,solution,r)
            summaries.append(sm)
            for j in range(len(n)):
                row={'case':label}
                for key,value in r.items():
                    a=np.asarray(value)
                    row[key]=float(a if a.ndim==0 else a[j])
                rows.append(row)
            check(label+'_energy',sm['max_energy_residual'])
            check(label+'_friedmann',sm['max_friedmann_residual'])
            if sm['complete']:
                tight=model.integrate(tight=True)
                if not tight.success or abs(tight.t[-1])>1e-12:
                    raise RuntimeError('Tighter-tolerance run did not reach N=0')
                rt=model.evaluate(tight,n)
                diffs={}
                for key in ('u','q','E','tau','rho','delta_alpha'):
                    scale=max(float(np.max(abs(rt[key]))),1e-12 if key=='delta_alpha' else 1.)
                    diffs[key]=float(np.max(abs(r[key]-rt[key]))/scale)
                    check(label+'_convergence_'+key,diffs[key])
                sm['tight_comparison']=diffs
            print(label,json.dumps({x:sm.get(x) for x in ('complete','delta_alpha_zi_ppm',
                                                         'clock_reverses','max_energy_residual')}),flush=True)
        except Exception as error:
            summaries.append({'case':label,'exception':repr(error),'complete':False})
            checks.append({'name':label+'_execution','passed':False,'exception':repr(error)})
    if 'chi_probe' in trajectories:
        r=trajectories['chi_probe']; old=reference.evaluate_n(NGRID)
        for key,other in (('u',old['chi']-inputs['chi_i']),('E',old['E']),('tau',old['tau'])):
            check('prior_reference_'+key,np.max(abs(r[key]-other))/max(float(np.max(abs(other))),1.))
        dump(OUT/'radiative_budget.json',budget(inputs,r))
    write_csv(OUT/'trajectories.csv',rows)
    payload={'run_utc':datetime.now(timezone.utc).isoformat(), 'elapsed_seconds':time.time()-started,
             'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,
             'script_sha256':sha(__file__),'protocol_sha256':sha(ROOT/'protocol.md'),
             'inputs_sha256':sha(OUT/'inputs.json'),'summaries':summaries,'checks':checks,
             'check_pass_count':sum(x['passed'] for x in checks),'check_count':len(checks),
             'interpretation':'Implementation checks only; probe and light-mass examples are not viable-particle or observational claims.'}
    dump(OUT/'summary.json',payload)
    print('checks',payload['check_pass_count'],'/',len(checks),'seconds',payload['elapsed_seconds'],flush=True)
    if not all(x['passed'] for x in checks):
        raise SystemExit(1)


if __name__=='__main__':
    main()
