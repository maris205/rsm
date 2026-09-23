#!/usr/bin/env python3
"""Compare independent trajectories; produce deterministic toy-model figures.

No observations, optimization, phase alignment or parameter fitting. The
comparison gates check implementation and separation of numerical errors.
"""
from pathlib import Path
import csv
import hashlib
import json
import platform

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FixedFormatter, NullFormatter

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results'
FIG = ROOT / 'figures'
INPUTS = {}
CHECKS = []


def read(name):
    path = OUT / (name + '.npz')
    INPUTS[path.name] = hashlib.sha256(path.read_bytes()).hexdigest()
    with np.load(path, allow_pickle=False) as data:
        return {key: data[key] for key in data.files}


def check(name, passed, **details):
    CHECKS.append(dict(name=name, passed=bool(passed), **details))


def rel_rms(q, ref):
    return float(np.sqrt(np.mean(np.abs(q-ref)**2) / np.mean(np.abs(ref)**2)))


def maxabs(q):
    return float(np.max(np.abs(q)))


def sample(data, x, key='qhat'):
    return (data[key] @ np.exp(1j*np.outer(data['k'], x))).real


def save_table(name, rows):
    with (OUT / name).open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def orders(rows, key):
    return [float(np.log2(a[key]/b[key])) for a,b in zip(rows,rows[1:])]


def main():
    FIG.mkdir(exist_ok=True)
    ext = read('reference_external_N128')
    frozen = read('reference_frozen_N128')
    ref = {i: read(f'reference_full_N128_I{i}') for i in (10,100,1000)}
    full = {i: read(f'lattice_main_I{i}') for i in (10,100,1000)}
    continuum = read('reference_continuum_I100_K24_tight')
    t = ext['t']
    drift_scale = maxabs(ext['mass_squared']-ext['mass_squared'][0])
    feedback = []
    for inertia in ref:
        d, kd = ref[inertia], full[inertia]
        check(f'I{inertia}_same_times_initial_field',
              np.array_equal(t, d['t']) and np.array_equal(t,kd['t'])
              and np.array_equal(d['q'][0],ext['q'][0])
              and np.array_equal(d['q'][0],kd['q'][0]))
        h0 = float(d['total_energy'][0])
        row = dict(I=inertia, initial_clock_energy=float(d['clock_energy'][0]),
                   final_R=float(d['R'][-1]),
                   final_R_fractional_departure=float(d['R'][-1]/ext['R'][-1]-1),
                   final_mass_squared=float(d['mass_squared'][-1]),
                   max_mass_difference_from_external=maxabs(d['mass_squared']-ext['mass_squared']),
                   mass_difference_over_external_drift=maxabs(d['mass_squared']-ext['mass_squared'])/drift_scale,
                   field_rms_vs_external=rel_rms(d['q'],ext['q']),
                   field_rms_vs_frozen=rel_rms(d['q'],frozen['q']),
                   field_energy_change=float(d['field_energy'][-1]-d['field_energy'][0]),
                   clock_energy_gain=float(d['clock_energy'][-1]-d['clock_energy'][0]),
                   independent_work=float(d['work'][-1]),
                   max_abs_energy_residual=maxabs(d['total_energy']-h0),
                   max_abs_field_work_residual=maxabs(d['field_energy']-d['field_energy'][0]-d['work']),
                   verlet_field_rms_vs_same_grid=rel_rms(kd['q'],d['q']),
                   verlet_max_R_error=maxabs(kd['R']-d['R']),
                   verlet_max_relative_energy_residual=maxabs(kd['Htotal']-kd['Htotal'][0])/float(kd['Htotal'][0]))
        feedback.append(row)
        check(f'I{inertia}_feedback_direction', np.min(d['R']-ext['R'])>=-1e-12
              and np.max(d['mass_squared']-ext['mass_squared'])<=1e-12)
        check(f'I{inertia}_energy_work_consistent',
              row['max_abs_energy_residual']<1e-9 and row['max_abs_field_work_residual']<1e-9)
        check(f'I{inertia}_independent_main_accuracy',
              row['verlet_field_rms_vs_same_grid']<1e-5 and row['verlet_max_R_error']<1e-6)
    for key in ('final_R_fractional_departure','max_mass_difference_from_external','field_rms_vs_external'):
        check(f'feedback_decreases_with_I_{key}', all(a[key]>b[key] for a,b in zip(feedback,feedback[1:])))
    control = dict(external_final_mass_squared=float(ext['mass_squared'][-1]),
                   external_field_rms_vs_frozen=rel_rms(ext['q'],frozen['q']),
                   external_drift_scale=drift_scale,
                   external_max_abs_field_work_residual=maxabs(ext['field_energy']-ext['field_energy'][0]-ext['work']),
                   frozen_max_abs_energy_residual=maxabs(frozen['field_energy']-frozen['field_energy'][0]))
    check('external_and_frozen_same_initial_field',np.array_equal(ext['q'][0],frozen['q'][0])
          and np.array_equal(t,frozen['t']))
    check('external_work_balance',control['external_max_abs_field_work_residual']<1e-9)
    spatial = []
    for n in (64,128,256):
        d = ref[100] if n==128 else read(f'reference_full_N{n}_I100')
        check(f'spatial_N{n}_time_and_initial_alignment',np.array_equal(t,d['t'])
              and np.array_equal(t,continuum['t'])
              and maxabs(d['q'][0]-sample(continuum,d['x'])[0])<1e-13)
        qref = sample(continuum, d['x'])
        spatial.append(dict(N=n, spacing=float(d['spacing']),
                            field_relative_rms=rel_rms(d['q'],qref),
                            max_field_rms=float(np.max(np.sqrt(np.mean((d['q']-qref)**2,axis=1)))),
                            max_R_error=maxabs(d['R']-continuum['R']),
                            max_mass_squared_error=maxabs(d['mass_squared']-continuum['mass_squared'])))
    spatial_orders = orders(spatial,'field_relative_rms')
    check('spatial_second_order',all(1.9 < o < 2.1 for o in spatial_orders),orders=spatial_orders)
    time_rows = []
    for dt,label in ((.01,'time_I100_dt0p01'),(.005,'time_I100_dt0p005'),
                     (.0025,'time_I100_dt0p0025'),(.00125,'main_I100')):
        kd = full[100] if label=='main_I100' else read('lattice_'+label)
        check(f'dt{dt}_alignment',np.array_equal(kd['t'],t)
              and np.array_equal(kd['q'][0],ref[100]['q'][0]))
        time_rows.append(dict(dt=dt,field_relative_rms=rel_rms(kd['q'],ref[100]['q']),
                              max_R_error=maxabs(kd['R']-ref[100]['R']),
                              max_abs_total_energy_residual=maxabs(kd['Htotal']-kd['Htotal'][0]),
                              max_abs_field_work_residual=maxabs(kd['field_work_residual'])))
    time_orders=orders(time_rows,'field_relative_rms')
    energy_orders=orders(time_rows,'max_abs_total_energy_residual')
    check('time_second_order',all(1.9 < o < 2.1 for o in time_orders),orders=time_orders)
    check('energy_second_order',all(1.9 < o < 2.1 for o in energy_orders),orders=energy_orders)
    c16=read('reference_continuum_I100_K16')
    c24=read('reference_continuum_I100_K24')
    x=2*np.pi*np.arange(256)/256
    errors = []
    for a,b in ((c16,c24),(c24,continuum)):
        dq=sample(a,x)-sample(b,x)
        errors.append(float(np.max(np.sqrt(np.mean(dq*dq,axis=1)))))
    precision=dict(max_continuum_cutoff_rms=errors[0],max_continuum_tolerance_rms=errors[1],
                   continuum_uncertainty_over_smallest_spatial_difference=max(errors)/spatial[-1]['max_field_rms'],
                   finest_time_error_over_N128_spatial_error=time_rows[-1]['field_relative_rms']/spatial[1]['field_relative_rms'])
    check('reference_separation',precision['continuum_uncertainty_over_smallest_spatial_difference']<.01)
    check('main_time_error_below_one_percent_spatial',precision['finest_time_error_over_N128_spatial_error']<.01)
    for path,rows in (('comparison_feedback.csv',feedback),('comparison_spatial.csv',spatial),('comparison_time.csv',time_rows)):
        save_table(path,rows)
    summary=dict(scope='Deterministic dimensionless trajectories, no fit or observations.',
                 normalization='Field RMS=sqrt(mean(|q-qref|^2)/mean(|qref|^2)) over all saved time/space samples. No phase or amplitude alignment. Feedback uses external denominator; frozen comparison uses frozen denominator. Separate R and coefficient differences.',
                 comparison_checks_note='Numerical consistency diagnostics; not preregistered statistical tests or physical confirmation.',
                 feedback=feedback,controls=control,spatial=spatial,spatial_orders=spatial_orders,
                 time=time_rows,time_orders=time_orders,energy_orders=energy_orders,precision=precision,
                 input_sha256=INPUTS,code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                 software=dict(python=platform.python_version(),numpy=np.__version__,matplotlib=matplotlib.__version__),
                 checks=CHECKS,all_checks_passed=all(c['passed'] for c in CHECKS))
    (OUT/'comparison_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    plots(t,ref,ext,frozen,spatial,time_rows)
    print(json.dumps({k:summary[k] for k in ('feedback','controls','spatial_orders','time_orders','precision','all_checks_passed')},indent=2))
    print(f'{len(CHECKS)} comparison checks; figures in {FIG}')
    if not summary['all_checks_passed']:
        raise SystemExit('A comparison check failed; inspect saved evidence.')


def plots(t,ref,ext,frozen,spatial,time_rows):
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,
                         'axes.titleweight':'semibold','legend.frameon':False,'savefig.dpi':180})
    colors={10:'#be4b39',100:'#2078a1',1000:'#5b9170'}
    fig,axs=plt.subplots(2,2,figsize=(11.4,7.5),layout='constrained')
    fig.suptitle('Inverse-log response with a dynamical clock',fontsize=15,weight='semibold')
    for inertia,d in ref.items():
        axs[0,0].plot(t,d['mass_squared'],color=colors[inertia],label=f'Full, I={inertia}')
        axs[0,1].plot(t,100*(d['R']/ext['R']-1),color=colors[inertia],label=f'I={inertia}')
    axs[0,0].plot(t,ext['mass_squared'],'k--',lw=1.4,label='Prescribed log law')
    axs[0,0].plot(t,frozen['mass_squared'],color='0.55',ls=':',label='Frozen coefficient')
    axs[0,0].set(title='(a) Restoring coefficient',ylabel=r'$M^2$',ylim=(.65,1.03))
    axs[0,0].legend(fontsize=8.5,ncol=2,loc='upper right')
    axs[0,1].set(title='(b) Clock backreaction',ylabel=r'$100(R/R_{\rm ext}-1)$ [%]')
    axs[0,1].legend()
    for d,color,style,label in ((ref[100],colors[100],'-','Full, I=100'),
                              (ext,'#d19532','--','Prescribed log law'),
                              (frozen,'0.5',':','Frozen coefficient')):
        axs[1,0].plot(t,d['cos_modes'][:,0],color=color,ls=style,label=label)
    axs[1,0].set(title='(c) Field mode, same initial state',ylabel=r'$C_1=2\,\mathrm{Re}\,\hat q_1$')
    axs[1,0].legend(fontsize=9,loc='lower left')
    d=ref[100]
    axs[1,1].plot(t,d['clock_energy']-d['clock_energy'][0],color=colors[100],label=r'$\Delta H_{\rm clock}$')
    axs[1,1].plot(t,d['field_energy']-d['field_energy'][0],color=colors[10],label=r'$\Delta H_{\rm field}$')
    axs[1,1].plot(t[::40],d['work'][::40],'o',ms=3.5,mfc='none',color='0.25',label=r'Independent $\int P_{\rm field}d\tau$')
    axs[1,1].axhline(0,color='0.7',lw=.7)
    axs[1,1].set(title='(d) Energy transfer, I=100',ylabel='Energy change [model units]')
    axs[1,1].legend(fontsize=9)
    for ax in axs.flat:
        ax.set_xlabel(r'Elapsed model time $\tau$')
        ax.grid(alpha=.15)
    fig.supxlabel('Chosen response; dimensionless benchmark, no observational fit.',fontsize=9)
    for suffix in ('png','pdf'):
        fig.savefig(FIG/f'log_clock_response.{suffix}')
    plt.close(fig)
    fig,axs=plt.subplots(2,2,figsize=(11.4,7.4),layout='constrained')
    fig.suptitle('Field evolution and independent numerical convergence',fontsize=15,weight='semibold')
    d=ref[100]
    delta=d['q']-frozen['q']
    for ax,arr,title,label in ((axs[0,0],d['q'],'(a) Full coupling, I=100',r'$q(x,\tau)$'),
                             (axs[0,1],delta,'(b) Departure from frozen coefficient',r'$q_{\rm full}-q_{\rm frozen}$')):
        bound=maxabs(arr)
        im=ax.imshow(arr,extent=(0,2*np.pi,t[-1],0),aspect='auto',cmap='RdBu_r',vmin=-bound,vmax=bound)
        ax.set(title=title,xlabel=r'Position $x$',ylabel=r'Elapsed time $\tau$')
        ax.set_xticks([0,np.pi,2*np.pi],['0',r'$\pi$',r'$2\pi$'])
        fig.colorbar(im,ax=ax,shrink=.9,label=label)
    for ax,rows,key,title,xlabel in ((axs[1,0],spatial,'spacing','(c) Space: coupled lattice vs continuum',r'Spacing $a$'),
                                  (axs[1,1],time_rows,'dt','(d) Time: Verlet vs independent ODE',r'Time step $\Delta\tau$')):
        xs=np.array([r[key] for r in rows]); ys=np.array([r['field_relative_rms'] for r in rows])
        ax.loglog(xs,ys,'o-',color=colors[100],label='Measured field error')
        ax.loglog(xs,ys[-1]*(xs/xs[-1])**2,'--',color='0.5',label='Second-order guide')
        ax.set(title=title,xlabel=xlabel,ylabel='Relative RMS over the full window')
        ax.xaxis.set_major_locator(FixedLocator(xs))
        ax.xaxis.set_major_formatter(FixedFormatter([f'{x:.4g}' for x in xs]))
        ax.xaxis.set_minor_formatter(NullFormatter())
        ax.grid(which='both',alpha=.15)
        ax.legend(fontsize=9)
    fig.supxlabel('Finite window 0 <= time <= 20; no phase alignment or parameter fitting.',fontsize=9)
    for suffix in ('png','pdf'):
        fig.savefig(FIG/f'log_clock_field_bridge.{suffix}')
    plt.close(fig)


if __name__=='__main__':
    main()
