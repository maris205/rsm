#!/usr/bin/env python3
"""Compare saved independent solutions; do not refit their phase or parameters."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('OMP_NUM_THREADS', '1')
from pathlib import Path
import csv
import hashlib
import json
import platform
import numpy as np
import scipy
from scipy.interpolate import CubicHermiteSpline
from scipy.optimize import brentq
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import NullLocator

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / 'results'
FIGURES = ROOT / 'figures'
used = {}


def load(name):
    path = RESULTS / name
    used[name] = hashlib.sha256(path.read_bytes()).hexdigest()
    return dict(np.load(path, allow_pickle=False))


def reconstruct(data, x, key='qhat'):
    values = data[key] @ np.exp(1j * np.outer(data['k'], x))
    if np.max(np.abs(values.imag)) > 1e-9:
        raise ValueError('Reality condition failed in Fourier reconstruction')
    return values.real


def relative_rms(value, reference):
    return float(np.linalg.norm(value-reference) / np.linalg.norm(reference))


def max_spatial_rms(value, reference):
    return float(np.sqrt(np.mean((value-reference)**2, axis=1)).max())


def modes(data, numbers=(1, 3, 5)):
    if 'x' in data:
        basis = np.cos(np.outer(data['x'], numbers))
        return (2/len(data['x']) * data['q'] @ basis,
                2/len(data['x']) * data['velocity'] @ basis)
    indices = [int(np.flatnonzero(data['k'] == k)[0]) for k in numbers]
    return 2*data['qhat'][:, indices].real, 2*data['velocity_hat'][:, indices].real


def crossing_frequency(t, q, velocity):
    spline = CubicHermiteSpline(t, q, velocity)
    intervals = np.flatnonzero((q[:-1] > 0) & (q[1:] <= 0))
    zeros = [brentq(spline, t[j], t[j+1], xtol=1e-13) for j in intervals]
    periods = np.diff(zeros)
    return dict(omega=float(2*np.pi / np.mean(periods)),
                descending_zero_times= zeros,
                period_min=float(min(periods)), period_max=float(max(periods)),
                definition='Mean descending-zero period over the finite window; not a fitted or asymptotic frequency.')


def tail_fraction(data, spacing, threshold):
    weights = (1+data['k']**2)[None, :]*np.abs(data['qhat'])**2 + np.abs(data['velocity_hat'])**2
    mask = np.abs(data['k'])*spacing > threshold
    return float(np.max(weights[:, mask].sum(axis=1)/weights.sum(axis=1)))


def csv_write(name, rows):
    with (RESULTS/name).open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)


def main():
    FIGURES.mkdir(exist_ok=True)
    reference = load('reference_continuum_K24_tight.npz')
    t = reference['t']
    cr, vr = modes(reference)
    rows = []
    checks = []
    def check(name, condition, **evidence):
        checks.append(dict(name=name, passed=bool(condition), **evidence))
    previous = None
    fields = {}
    for n in (64, 128, 256, 512):
        raw = load(f'lattice_spatial_N{n}.npz')
        exact = load(f'reference_lattice_N{n}.npz')
        corrected = load(f'reference_corrected_N{n}.npz')
        x = exact['x']; a = 2*np.pi/n
        qref = reconstruct(reference, x)
        vref = reconstruct(reference, x, 'velocity_hat')
        qcorr = reconstruct(corrected, x)
        vcorr = reconstruct(corrected, x, 'velocity_hat')
        cq, cv = modes(raw)
        row = dict(N=n, spacing=a,
            dt=json.loads(str(raw['parameters_json']))['actual_dt'],
            q_error_continuum=relative_rms(exact['q'], qref),
            q_error_corrected=relative_rms(exact['q'], qcorr),
            velocity_error_continuum=relative_rms(exact['velocity'], vref),
            velocity_error_corrected=relative_rms(exact['velocity'], vcorr),
            q_error_verlet_to_lattice_reference=relative_rms(raw['q'], exact['q']),
            q_error_verlet_to_continuum=relative_rms(raw['q'], qref),
            max_spatial_rms_continuum=max_spatial_rms(exact['q'], qref),
            max_spatial_rms_corrected=max_spatial_rms(exact['q'], qcorr),
            q_error_final_continuum=relative_rms(exact['q'][-1], qref[-1]),
            q_error_final_corrected=relative_rms(exact['q'][-1], qcorr[-1]),
            energy_relative_max=float(np.max(np.abs(raw['energy']/raw['energy'][0]-1))),
            cosine3_max=float(np.max(np.abs(cq[:, 1]))),
            cosine3_relative_error=relative_rms(cq[:, 1], cr[:, 1]),
            omega1_from_crossings=crossing_frequency(t, cq[:, 0], cv[:, 0])['omega'],
            raw_tail_norm_fraction_ka_gt_half=tail_fraction(exact, a, .5),
            corrected_tail_norm_fraction_ka_gt_half=tail_fraction(corrected, a, .5),
            raw_tail_norm_fraction_ka_gt_one=tail_fraction(exact, a, 1.),
            corrected_tail_norm_fraction_ka_gt_one=tail_fraction(corrected, a, 1.))
        row['correction_improvement_factor'] = row['q_error_continuum']/row['q_error_corrected']
        for kind in ('continuum', 'corrected'):
            row[f'order_{kind}_from_previous'] = (None if previous is None else
                float(np.log2(previous[f'q_error_{kind}']/row[f'q_error_{kind}'])))
        check(f'N{n}_times_align', np.array_equal(raw['t'], t) and np.array_equal(exact['t'], t)
              and np.array_equal(corrected['t'], t))
        check(f'N{n}_correction_improves', row['q_error_corrected'] < row['q_error_continuum'])
        if previous:
            check(f'N{n}_spatial_order', 1.8 < row['order_continuum_from_previous'] < 2.2,
                  order=row['order_continuum_from_previous'])
            check(f'N{n}_corrected_order', 3.5 < row['order_corrected_from_previous'] < 4.5,
                  order=row['order_corrected_from_previous'])
        rows.append(row); previous = row
        fields[n] = (raw, exact, corrected, qref, qcorr)
    temporal = []
    temporal_data = []
    ref128 = fields[128][1]
    for dt in (.01, .005, .0025, .00125):
        name = str(dt).replace('.', 'p')
        raw = load(f'lattice_time_N128_dt{name}.npz')
        row = dict(dt=dt, N=128,
            q_error=relative_rms(raw['q'], ref128['q']),
            velocity_error=relative_rms(raw['velocity'], ref128['velocity']),
            energy_relative_max=float(np.max(np.abs(raw['energy']/raw['energy'][0]-1))))
        row['order_from_previous'] = (None if not temporal else
            float(np.log2(temporal[-1]['q_error']/row['q_error'])))
        if temporal:
            check(f'dt{dt}_time_order', 1.8 < row['order_from_previous'] < 2.2,
                  order=row['order_from_previous'])
        temporal.append(row); temporal_data.append(raw)
    # An analytic whole-Brillouin-zone scan exposes the long-wave boundary.
    a = 2*np.pi/64
    dispersion = []
    for k in range(1, 33):
        lattice = np.sqrt(1+4/a**2*np.sin(k*a/2)**2)
        continuum = np.sqrt(1+k*k)
        corrected = np.sqrt(1+k*k-a*a*k**4/12)
        dispersion.append(dict(k=k, ka=k*a, omega_lattice=lattice,
            omega_continuum=continuum, omega_corrected=corrected,
            relative_continuum_error=float(continuum/lattice-1),
            relative_corrected_error=float(corrected/lattice-1)))
    linear = load('lattice_linear_N64_k1.npz')
    cl, _ = modes(linear)
    check('nonlinearity_generates_third_spatial_mode', np.max(np.abs(cr[:, 1])) > 1e-4)
    check('linear_control_has_no_third_mode', np.max(np.abs(cl[:, 1])) < 1e-10,
          maximum=float(np.max(np.abs(cl[:, 1]))))
    continuum_quality = []
    for name in ('reference_continuum_K16.npz', 'reference_continuum_K24.npz'):
        alternative = load(name)
        qr = reconstruct(reference, fields[128][1]['x'])
        qa = reconstruct(alternative, fields[128][1]['x'])
        continuum_quality.append(dict(file=name, relative_rms=relative_rms(qa, qr)))
    check('continuum_reference_below_smallest_corrected_error',
          max(r['relative_rms'] for r in continuum_quality) < rows[-1]['q_error_corrected']/20,
          note='Only checks the continuum portion; finest lattice tolerance has a separate audit.')
    result = dict(scope='Deterministic fixed-coefficient phi4 benchmark, finite interval [0,20]; no observations or fitted phases.',
        metric='Equal-weight global relative RMS over all shared times and sites; denominator is the reference array (second argument). Spatial comparisons use the respective continuum/corrected field; time comparisons use the same-lattice ODE reference.',
        tail_metric='Quadratic spectral norm fraction with weight (1+k^2)|qhat|^2+|velocity_hat|^2; excludes quartic interaction energy.',
        code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        software=dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__, matplotlib=matplotlib.__version__),
        spatial=rows, temporal=temporal, reference_quality=continuum_quality,
        continuum_modes=dict(cosine3_max=float(np.max(np.abs(cr[:, 1]))),
            cosine3_max_relative_to_initial_A=float(np.max(np.abs(cr[:, 1]))/.5),
            effective_frequency=crossing_frequency(t, cr[:, 0], vr[:, 0])),
        dispersion_selected=[row for row in dispersion if row['k'] in (1, 8, 24, 32)],
        checks=checks, all_checks_passed=all(c['passed'] for c in checks), input_sha256=used)
    (RESULTS/'comparison_summary.json').write_text(json.dumps(result, ensure_ascii=False, indent=2)+'\n')
    csv_write('comparison_spatial.csv', rows); csv_write('comparison_temporal.csv', temporal)
    csv_write('comparison_dispersion.csv', dispersion)
    plot(fields, reference, rows, temporal, temporal_data, dispersion)
    print(json.dumps({k:result[k] for k in ('spatial','temporal','continuum_modes','all_checks_passed')},indent=2))
    if not result['all_checks_passed']:
        raise SystemExit('Comparison check failed; report actual diagnostics.')


def save(fig, stem):
    for ext in ('png','pdf'):
        fig.savefig(FIGURES/f'{stem}.{ext}', dpi=190)
    plt.close(fig)


def plot(fields, reference, rows, temporal, temporal_data, dispersion):
    plt.rcParams.update({'font.size':10, 'axes.spines.top':False, 'axes.spines.right':False})
    t = reference['t']; raw, exact, corrected, qr, qc = fields[64]
    fig, axes = plt.subplots(2, 2, figsize=(11.5, 7.7), constrained_layout=True)
    for ax, data, title in zip(axes[0], (raw['q'], qr), ('Local nonlinear lattice, N = 64', 'Independent continuum field')):
        im = ax.pcolormesh(raw['x'], t, data, cmap='RdBu_r', shading='auto', vmin=-.5, vmax=.5, rasterized=True)
        ax.set(xlabel='Position x', ylabel='Time t', title=title)
    fig.colorbar(im, ax=axes[0].tolist(), label='Field amplitude', shrink=.9)
    cr, _ = modes(reference); cl, _ = modes(raw); ce, _ = modes(corrected)
    for ax, index, label in zip(axes[1], (0, 1), ('Fundamental spatial mode', 'Third spatial harmonic')):
        ax.plot(t, cr[:, index], color='#2070a2', label='Continuum', lw=1.6)
        ax.plot(t[::18], cl[::18, index], 'o', ms=3, mfc='none', color='#ce7924', label='Lattice samples')
        if index == 1:
            ax.plot(t, ce[:, index], '--', color='#994571', lw=1.2, label='With spatial correction')
        ax.set(xlabel='Time t', ylabel=f'Cosine coefficient C{1 if index==0 else 3}', title=label)
        ax.grid(alpha=.2); ax.legend(frameon=False, fontsize=9)
    fig.suptitle('A controlled lattice-to-field bridge: synthetic dynamics, no fitted observations', fontsize=13)
    save(fig, 'micro_macro_evolution')

    fig, axes = plt.subplots(2, 2, figsize=(11.5, 7.7), constrained_layout=True)
    n = np.array([r['N'] for r in rows]); e = np.array([r['q_error_continuum'] for r in rows])
    ec = np.array([r['q_error_corrected'] for r in rows])
    ax = axes[0,0]
    ax.loglog(n, e, 'o-', color='#2070a2', label='Continuum: expected second order')
    ax.loglog(n, ec, 's-', color='#994571', label='With correction: expected fourth order')
    ax.set(xlabel='Number of lattice sites N', ylabel='Relative field RMS error', title='Space: compare accurate ODE references')
    ax.set_xticks(n, labels=[str(v) for v in n]); ax.xaxis.set_minor_locator(NullLocator())
    ax.legend(frameon=False,fontsize=8.5); ax.grid(alpha=.2,which='both')
    ax = axes[0,1]
    dt=np.array([r['dt'] for r in temporal]); er=np.array([r['q_error'] for r in temporal])
    ax.loglog(dt, er, 'o-', color='#ce7924', label='Verlet versus same-lattice reference')
    ax.loglog(dt, er[-1]*(dt/dt[-1])**2, '--',color='#333333',label='Second-order guide')
    ax.set(xlabel='Time step', ylabel='Relative field RMS error', title='Time: fixed N = 128')
    ax.set_xticks(dt[::-1], labels=[f'{v:g}' for v in dt[::-1]]); ax.xaxis.set_minor_locator(NullLocator())
    ax.legend(frameon=False,fontsize=8.5); ax.grid(alpha=.2,which='both')
    ax = axes[1,0]
    ka=np.array([r['ka'] for r in dispersion])
    ax.semilogy(ka, np.abs([r['relative_continuum_error'] for r in dispersion]),color='#2070a2',label='Continuum frequency error')
    ax.semilogy(ka, np.abs([r['relative_corrected_error'] for r in dispersion]),color='#994571',label='Corrected frequency error')
    ax.axvline(.5,color='#777777',ls=':',lw=1)
    ax.set(xlabel='Dimensionless wave number ka',ylabel='Absolute relative frequency error',title='Long-wave expansion fails near grid scales')
    ax.legend(frameon=False,fontsize=8.5); ax.grid(alpha=.2,which='both')
    ax = axes[1,1]
    for entry, trajectory in zip(temporal, temporal_data):
        ax.plot(t, 1e6*(trajectory['energy']/trajectory['energy'][0]-1),label=f"dt = {entry['dt']:g}",lw=1)
    ax.set(xlabel='Time t',ylabel='Relative energy deviation (parts per million)',title='Symplectic evolution: energy error is finite')
    ax.legend(frameon=False,fontsize=8.5); ax.grid(alpha=.2)
    save(fig, 'micro_macro_convergence')


if __name__ == '__main__':
    main()
