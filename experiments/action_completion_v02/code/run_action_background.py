#!/usr/bin/env python3
"""Autonomous homogeneous scalar background: a toy completion, no data fitting.

Natural units, reduced Planck mass. The matter + Lambda approximation omits
radiation even at the artificial initial a_i; it is not an early-Universe fit.
The scalar is evolved by its equation of motion, never reset to log(time).
"""
from pathlib import Path
import csv
import hashlib
import json
import platform

import numpy as np
import scipy
from scipy.integrate import solve_ivp, quad
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parents[1]
OM = 0.315
HREF_KM_S_MPC = 67.4
MPC_KM = 3.0856775814913673e19
YEAR_S = 365.25 * 86400.0
HREF_S = HREF_KM_S_MPC / MPC_KM
HREF_YR = HREF_S * YEAR_S
TSTAR_S = 5.391247e-44
D0_YR = 2.5e-19  # An illustrative normalization, not an observed detection.
EPSILONS = (1e-6, 1e-4, 1e-2)  # Not asserted to be observationally allowed.
Z = np.linspace(0.0, 4.2, 841)


def json_write(path, obj):
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2,
        default=lambda x: x.item() if isinstance(x, np.generic) else x.tolist()) + '\n')


def csv_write(path, rows):
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


class Background:
    def __init__(self, epsilon, ai=1e-3, omega_lambda=None,
                 rtol=2e-11, atol=2e-13, w_factor=1.0):
        self.epsilon, self.ai = epsilon, ai
        self.ni = float(np.log(ai))
        self.rtol, self.atol, self.w_factor = rtol, atol, w_factor
        self.ei_matter = np.sqrt(OM * ai**-3 / (1 - .75 * epsilon))
        self.tau_i = 2 / (3 * self.ei_matter)
        self.w_i = 1 / self.tau_i
        self.W_i = .5 / self.tau_i**2
        self.chi_i = float(np.log(self.tau_i / (HREF_S * TSTAR_S)))
        if omega_lambda is None:
            self.omega_lambda = brentq(self.present_closure, .0, .9,
                                       xtol=2e-13, rtol=2e-13)
        else:
            self.omega_lambda = omega_lambda
        self.solution = self.integrate(self.omega_lambda)
        self.calls = self.solution.nfev

    def e2(self, n, u, w, omega_lambda):
        return OM * np.exp(-3 * n) + omega_lambda + self.epsilon * (
            .5 * w*w + self.W_i * np.exp(-2*u)) / 3

    def integrate(self, omega_lambda):
        def rhs(n, y):
            u, w, _tau = y
            potential = self.W_i * np.exp(-2*u)
            e = np.sqrt(self.e2(n, u, w, omega_lambda))
            return (w/e, -3*w + 2*potential/e, 1/e)
        sol = solve_ivp(rhs, (self.ni, 0.),
            (0., self.w_i*self.w_factor, self.tau_i), method='DOP853',
            rtol=self.rtol, atol=(self.atol, self.atol, self.atol*.01),
            dense_output=True)
        if not sol.success:
            raise RuntimeError(sol.message)
        return sol

    def present_closure(self, omega_lambda):
        sol = self.integrate(omega_lambda)
        u, w, _ = sol.y[:, -1]
        return self.e2(0., u, w, omega_lambda) - 1.

    def evaluate_n(self, n):
        n = np.asarray(n)
        u, w, tau = self.solution.sol(n)
        W = self.W_i * np.exp(-2*u)
        energy = .5*w*w + W
        chi = self.chi_i + u
        e2 = self.e2(n, u, w, self.omega_lambda)
        return dict(n=n, u=u, w=w, tau=tau, W=W, energy=energy, chi=chi,
                    E=np.sqrt(e2), E2=e2,
                    omega_chi=self.epsilon*energy/(3*e2),
                    equation_of_state=(.5*w*w-W)/energy,
                    exact_log=np.log(tau/(HREF_S*TSTAR_S)),
                    delta_chi_from_exact_log=u-np.log(tau/self.tau_i),
                    age_year=tau/HREF_YR)

    def transfer(self, z):
        p = self.evaluate_n(0.)
        v = self.evaluate_n(-np.log1p(z))
        chi_dot_yr = p['w'] * HREF_YR
        # chi^-2 - chi0^-2 = chi0^-2 expm1[-2 log1p((chi-chi0)/chi0)].
        # Compute chi-chi0 from u-u0 so subtracting two large chi's is avoided.
        d = v['u']-p['u']
        tchi = -p['chi']/(2*chi_dot_yr) * np.expm1(-2*np.log1p(d/p['chi']))
        log_ratio_age = np.log(v['tau']/p['tau'])
        t2 = -.5*p['age_year']*p['exact_log'] * np.expm1(
            -2*np.log1p(log_ratio_age/p['exact_log']))
        beta = -D0_YR*p['chi']**3/(2*chi_dot_yr)
        qalpha = D0_YR/(chi_dot_yr*np.sqrt(self.epsilon))
        v.update(T_chi_year=tchi, T_log_same_history_year=t2,
                 alpha_ratio=1+D0_YR*tchi, beta=beta, qalpha0=qalpha)
        return v


def lcdm_transfer(z):
    tau = 2/(3*np.sqrt(1-OM)) * np.arcsinh(np.sqrt((1-OM)/OM)/(1+z)**1.5)
    tau0 = 2/(3*np.sqrt(1-OM)) * np.arcsinh(np.sqrt((1-OM)/OM))
    L0 = np.log(tau0/(HREF_S*TSTAR_S))
    return -.5*tau0/HREF_YR*L0 * np.expm1(-2*np.log1p(np.log(tau/tau0)/L0))


def relative_change(a, b):
    return float(np.max(np.abs(a/b-1)))


def run():
    for name in ('results', 'reports', 'figures'):
        (ROOT/name).mkdir(exist_ok=True)
    checks = []
    def check(name, passed, **evidence):
        checks.append(dict(name=name, passed=bool(passed), **evidence))

    result = dict(scope='Illustrative dust + Lambda + scalar background; no radiation, '
        'homogeneous electromagnetic source or data fitting; illustrative epsilon and D0 are not allowed bounds.',
        constants=dict(omega_m=OM, Href_km_s_Mpc=HREF_KM_S_MPC,
            year_seconds=YEAR_S, tstar_seconds=TSTAR_S, D0_example_per_year=D0_YR,
            initial_a=1e-3, epsilons=EPSILONS),
        software=dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__),
        code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        baseline=[], numerical_comparisons=[], matter_benchmarks=[], initial_condition_probes=[])
    rows = []
    models = []
    baseline_lcdm = lcdm_transfer(Z)

    for epsilon in EPSILONS:
        label = f'epsilon_{epsilon:g}'
        bg = Background(epsilon)
        models.append(bg)
        p = bg.evaluate_n(0.)
        v = bg.transfer(Z)
        tchi, tlog = v['T_chi_year'], v['T_log_same_history_year']
        mask = Z > 0
        identity_lhs = D0_YR/HREF_YR
        identity_rhs = float(v['qalpha0'] * np.sqrt(
            3*p['omega_chi']*(1+p['equation_of_state'])))
        item = dict(epsilon=epsilon, F_over_Mpl=np.sqrt(epsilon),
            omega_lambda=bg.omega_lambda, omega_chi0=p['omega_chi'],
            w_chi0=p['equation_of_state'], age0_Gyr=p['age_year']/1e9,
            chi0=p['chi'], exact_log0=p['exact_log'],
            delta_chi0_from_exact_log=p['delta_chi_from_exact_log'],
            chi_dot0_per_year=p['w']*HREF_YR, chi_dot0_times_age=p['w']*p['tau'],
            beta_for_example_D0=v['beta'], qalpha0_for_example_D0=v['qalpha0'],
            D0_over_H0=identity_lhs, qalpha_identity_rhs=identity_rhs,
            min_alpha_ratio=float(np.min(v['alpha_ratio'])),
            omega_chi_min=float(np.min(v['omega_chi'])),
            omega_chi_max=float(np.max(v['omega_chi'])),
            w_chi_min=float(np.min(v['equation_of_state'])),
            w_chi_max=float(np.max(v['equation_of_state'])),
            T_chi_z4p2_year=tchi[-1], T_log_same_history_z4p2_year=tlog[-1],
            T_chi_over_same_history_log_z4p2=tchi[-1]/tlog[-1],
            max_relative_transfer_deviation_same_history=relative_change(tchi[mask],tlog[mask]),
            T_chi_over_LCDM_log_z4p2=tchi[-1]/baseline_lcdm[-1],
            max_relative_transfer_deviation_LCDM=relative_change(tchi[mask],baseline_lcdm[mask]),
            max_relative_log_transfer_background_change=relative_change(tlog[mask],baseline_lcdm[mask]),
            max_abs_delta_chi_from_exact_log=float(np.max(np.abs(v['delta_chi_from_exact_log']))),
            solver_function_evaluations=bg.calls)
        result['baseline'].append(item)
        check(label+'_present_Friedmann_closure', abs(p['E']-1)<1e-10, E0=p['E'])
        check(label+'_positive_background_and_alpha',
            np.all(v['E']>0) and np.all(v['chi']>0) and np.all(v['alpha_ratio']>0),
            min_E=np.min(v['E']),min_chi=np.min(v['chi']), min_alpha_ratio=np.min(v['alpha_ratio']))
        check(label+'_canonical_energy_positivity',
            np.all(v['omega_chi']>0) and np.all(v['equation_of_state']>=-1)
            and np.all(v['equation_of_state']<=1))
        check(label+'_qalpha_identity', abs(identity_rhs/identity_lhs-1)<1e-12,
            relative_error=abs(identity_rhs/identity_lhs-1))

        # Independently integrate the continuity equation over disjoint intervals.
        edges = np.linspace(bg.ni, 0, 13)
        continuity_errors = []
        for a, b in zip(edges[:-1], edges[1:]):
            ea, eb = bg.evaluate_n(a)['energy'], bg.evaluate_n(b)['energy']
            loss = quad(lambda n: -3*float(bg.evaluate_n(n)['w'])**2,
                        a, b, epsabs=1e-10, epsrel=2e-10, limit=100)[0]
            continuity_errors.append(abs((eb-ea)-loss)/max(abs(eb-ea),1e-100))
        check(label+'_integrated_scalar_energy_balance',max(continuity_errors)<2e-8,
            max_relative_interval_residual=max(continuity_errors))

        # Five-point differences of E^2 use dense-output values, not the RHS.
        ngrid = np.linspace(bg.ni+.01,-.001,300)
        h = 1e-4
        fm2, fm1, fp1, fp2 = [bg.evaluate_n(ngrid+k*h)['E2'] for k in (-2,-1,1,2)]
        finite = (fm2-8*fm1+8*fp1-fp2)/(12*h)
        expected = -3*OM*np.exp(-3*ngrid)-epsilon*bg.evaluate_n(ngrid)['w']**2
        error = relative_change(finite,expected)
        check(label+'_finite_difference_Friedmann_derivative',error<2e-7,
              max_relative_error=error)

        tight = Background(epsilon,rtol=5e-13,atol=2e-15)
        early = Background(epsilon,ai=1e-4,rtol=5e-13,atol=2e-15)
        for kind, alternative, threshold in [('tighter_tolerance',tight,2e-8),
                                             ('earlier_initial_a',early,2e-6)]:
            vt = alternative.transfer(Z)
            err_t = relative_change(vt['T_chi_year'][mask],tchi[mask])
            err_e = relative_change(vt['E'],v['E'])
            err_age = relative_change(vt['age_year'],v['age_year'])
            comp = dict(epsilon=epsilon,change=kind,max_relative_Tchi_change=err_t,
                max_relative_E_change=err_e,max_relative_age_change=err_age,
                omega_lambda_change=alternative.omega_lambda-bg.omega_lambda,
                max_abs_chi_change=float(np.max(np.abs(vt['chi']-v['chi']))))
            result['numerical_comparisons'].append(comp)
            check(label+'_'+kind,max(err_t,err_e,err_age)<threshold,**comp)

        # Lambda=0 is an independently known exact matter scaling solution.
        matter = Background(epsilon,omega_lambda=0.)
        nv = np.linspace(matter.ni,0,1001)
        vm = matter.evaluate_n(nv)
        tau_exact = matter.tau_i*np.exp(1.5*(nv-matter.ni))
        matter_item = dict(epsilon=epsilon,
            max_relative_tau_error=relative_change(vm['tau'],tau_exact),
            max_relative_w_error=relative_change(vm['w'],1/tau_exact),
            max_abs_chi_minus_log_time=float(np.max(np.abs(vm['delta_chi_from_exact_log']))),
            max_relative_omega_chi_error=relative_change(vm['omega_chi'],np.full(len(nv),.75*epsilon)),
            max_abs_scalar_equation_of_state=float(np.max(np.abs(vm['equation_of_state']))))
        result['matter_benchmarks'].append(matter_item)
        check(label+'_pure_matter_analytic_scaling',
            max(value for key,value in matter_item.items() if key!='epsilon')<2e-8,**matter_item)
        for j,z in enumerate(Z):
            rows.append(dict(epsilon=epsilon,z=z,age_Gyr=v['age_year'][j]/1e9,
                E=v['E'][j],chi=v['chi'][j],exact_log_time=v['exact_log'][j],
                chi_minus_exact_log=v['delta_chi_from_exact_log'][j],
                omega_chi=v['omega_chi'][j],w_chi=v['equation_of_state'][j],
                T_chi_year=tchi[j],T_log_same_history_year=tlog[j],
                T_log_LCDM_year=baseline_lcdm[j],
                alpha_ratio_example=v['alpha_ratio'][j]))

    reference = models[1]
    ref = reference.transfer(Z)
    for w_factor in (.9,1.1):
        alt = Background(1e-4,w_factor=w_factor)
        val = alt.transfer(Z)
        pi = alt.evaluate_n(0.)
        result['initial_condition_probes'].append(dict(epsilon=1e-4,w_initial_factor=w_factor,
            omega_lambda=alt.omega_lambda,chi0_change=pi['chi']-reference.evaluate_n(0.)['chi'],
            max_relative_Tchi_change=relative_change(val['T_chi_year'][1:],ref['T_chi_year'][1:]),
            present_w_relative_change=pi['w']/reference.evaluate_n(0.)['w']-1,
            note='Finite-interval sensitivity only; not an attractor proof.'))
        check(f'initial_w_factor_{w_factor}_finite_positive',
            np.all(np.isfinite(val['T_chi_year'])) and np.all(val['E']>0))

    result['checks']=checks
    result['all_checks_passed']=all(c['passed'] for c in checks)
    json_write(ROOT/'results/background_summary.json',result)
    csv_write(ROOT/'results/background_grid.csv',rows)
    csv_write(ROOT/'results/background_summary.csv',result['baseline'])
    csv_write(ROOT/'results/initial_condition_probes.csv',result['initial_condition_probes'])
    make_plot(models)
    print(json.dumps(dict(all_checks_passed=result['all_checks_passed'],
        checks=len(checks),baseline=result['baseline'],
        failed_checks=[c for c in checks if not c['passed']]),indent=2))
    if not result['all_checks_passed']:
        raise SystemExit('Numerical checks failed; inspect background_summary.json.')


def make_plot(models):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(2,2,figsize=(10.5,7.6),constrained_layout=True)
    colors = ('#226c9a','#d18117','#ac3160')
    for bg,color in zip(models,colors):
        v=bg.transfer(Z)
        label=rf'$\epsilon={bg.epsilon:g}$'
        ax[0,0].plot(Z,v['delta_chi_from_exact_log'],color=color,label=label)
        ax[0,1].plot(Z,v['omega_chi']/bg.epsilon,color=color,label=label)
        ax[1,0].plot(Z[1:],v['T_chi_year'][1:]/v['T_log_same_history_year'][1:],color=color,label=label)
        ax[1,1].plot(Z,v['T_chi_year']/1e9,color=color,label=label)
    base=models[1].transfer(Z)
    ax[1,1].plot(Z,base['T_log_same_history_year']/1e9,color='#333333',ls='--',label='Exact-log ansatz, same history')
    ax[0,0].set_ylabel(r'$\chi-\ln(t/t_*)$')
    ax[0,1].set_ylabel(r'$\Omega_\chi/\epsilon$')
    ax[1,0].set_ylabel(r'$T_\chi/T_{\log,\mathrm{same\ history}}$')
    ax[1,1].set_ylabel('Drift-normalized transfer (Gyr)')
    for a in ax.flat:
        a.set_xlabel('Redshift z')
        a.grid(alpha=.2)
    ax[0,0].legend(frameon=False,fontsize=9)
    ax[1,1].legend(frameon=False,fontsize=8)
    ax[0,0].set_title('Autonomous field departs from exact log time')
    ax[0,1].set_title('Canonical scalar energy fraction')
    ax[1,0].set_title('Effect of the field history at fixed time history')
    ax[1,1].set_title('Illustrative background; no observations fitted')
    for ext in ('png','pdf'):
        fig.savefig(ROOT/f'figures/background_completion.{ext}',dpi=180)
    plt.close(fig)


if __name__=='__main__':
    run()
