#!/usr/bin/env python3
"""Execute the six frozen compensated-sphere PM cases and retain failures.

Run validation first. --preflight performs exactly three force evaluations;
it never integrates a production trajectory or changes scientific settings.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor,as_completed
from datetime import datetime,timezone
import gc
import json
from pathlib import Path
import platform
import resource
import time
import traceback

import numpy as np
import scipy
from scipy.optimize import brentq

import pm_refined as pm
import spherical_collapse as sphere

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results'
FROZEN_PROTOCOL='f3a0cc5abf21a4374439db756cbf8c4a3103bcabfa2bee6f9bd1921b28e608f2'
MATRIX=(
    ('ref_force128','reference',64,128,False,(0,0,0)),
    ('ref_main','reference',64,256,False,(0,0,0)),
    ('ref_halfstep','reference',64,256,True,(0,0,0)),
    ('ref_particles128','reference',128,256,False,(0,0,0)),
    ('ref_shift','reference',64,256,False,(.5,.25,.375)),
    ('clock_main','clock',64,256,False,(0,0,0)),
)


def write_json(path,obj):
    temp=path.with_suffix(path.suffix+'.tmp')
    temp.write_text(json.dumps(sphere.serializable(obj),indent=2,ensure_ascii=False,allow_nan=False)+'\n')
    temp.replace(path)


def check(name,value,threshold=None):
    passed=(False if value is None else bool(value) if threshold is None
            else bool(np.isfinite(value) and abs(value)<=threshold))
    result={'name':name,'value':value,'threshold':threshold,'passed':passed}
    if value is None:
        result['reason']='required numerical result or event is missing'
    return result


def inputs():
    if pm.sha(ROOT/'protocol.md')!=FROZEN_PROTOCOL:
        raise RuntimeError('Protocol differs from frozen bytes')
    summary=json.loads((OUT/'sphere_summary.json').read_text())
    if summary['protocol_sha256']!=FROZEN_PROTOCOL or summary['summary']['failures']:
        raise RuntimeError('A validated matching spherical ODE calibration is required')
    backgrounds=sphere.load_backgrounds()
    delta=float(summary['delta_i'])
    initials=sphere.initial_conditions(delta,backgrounds)
    ode={label:sphere.integrate_sphere(bg,initials[label]) for label,bg in backgrounds.items()}
    pm_backgrounds={label:pm.Background.from_csv(bg.path,.315,label) for label,bg in backgrounds.items()}
    return summary,backgrounds,pm_backgrounds,ode,delta


def preflight(args):
    _,backgrounds,_,_,delta=inputs()
    entries=[]
    for nparticle,nmesh in ((64,128),(64,256),(128,256)):
        initial=pm.compensated_initial(nparticle,delta,backgrounds['reference'],force_mesh=nmesh)
        mesh=pm.PMGrid(nmesh,args.fft_workers)
        start=time.perf_counter()
        _,diag=mesh.force(initial['positions'])
        entries.append({'nparticle':nparticle,'nmesh':nmesh,'force_seconds':time.perf_counter()-start,
                        'process_peak_rss_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                        'diagnostics':diag})
        print(json.dumps(entries[-1]),flush=True)
        del initial,mesh
        gc.collect()
    result={'protocol_sha256':FROZEN_PROTOCOL,'force_evaluations':3,'entries':entries,
            'timing_scope':'Sequential initial-state force calls, not total integration benchmark',
            'fft_workers':args.fft_workers,'code_sha256':pm.sha(__file__)}
    write_json(OUT/'pm_preflight.json',result)
    return result


def run_case(task,args,backgrounds,pm_backgrounds,ode,delta):
    name,label,nparticle,nmesh,halfstep,shift=task
    start=time.perf_counter()
    initial=pm.compensated_initial(nparticle,delta,backgrounds['reference'],shift,nmesh)
    bg=pm_backgrounds[label]
    x,p=initial['positions'].copy(),initial['momenta'].copy()
    mesh=pm.PMGrid(nmesh,args.fft_workers)
    force,first_diag=mesh.force(x)
    core_vectors=pm.minimum_vector(x[initial['core_mask']],initial['center'])
    core_r=np.linalg.norm(core_vectors,axis=1)
    radial_g=np.sum(force[initial['core_mask']]*core_vectors,axis=1)/core_r
    exact_g=-delta/3*core_vectors
    initial_force_diagnostic={
        'mean_radial_g':float(np.mean(radial_g)),
        'median_radial_g':float(np.median(radial_g)),
        'fraction_particles_with_inward_g':float(np.mean(radial_g<0)),
        'relative_L2_error_vs_uniform_core':float(np.linalg.norm(force[initial['core_mask']]-exact_g)/np.linalg.norm(exact_g)),
        'definition':'g=-grad(psi), laplacian(psi)=delta; radial_g=g dot (x-center)/|x-center|, negative means inward; box units',
        'reference':'continuum compensated spherical core g_exact=-(delta_i/3)*(x-center); fixed core labels only',
        'scope':'physical particle/force-grid discretization diagnostic, not an FFT implementation check',
    }
    mean_p0=np.mean(p,axis=0)
    a=.02
    history=[]
    event=None
    turnaround=None
    partial_event_calls=0
    nsteps=0
    ode_event=ode[label].events['200']

    def row_at(a,x,p,force_diag,step_diag=None):
        row=pm.core_diagnostics(x,p,a,bg,initial)
        row.update(force_diag)
        row.update(step_diag or {'max_drift_cells':0.,'actual_dt_over_tdyn':0.,'dt_H0':0.,'delta_ln_a':0.,'step_halvings':0})
        row['t_Gyr']=float(backgrounds[label].evaluate_a(a)['t_Gyr'])
        row['mean_momentum']=np.mean(p,axis=0).tolist()
        row['momentum_change_absolute']=float(np.linalg.norm(np.mean(p,axis=0)-mean_p0))
        if np.log(a)<=ode[label].nmax+1e-12:
            ref=ode[label].evaluate_a(a)
            row['ode_y']=float(ref['y'])
            row['relative_y_error_vs_ode']=row['y_median']/row['ode_y']-1
        else:
            row['ode_y']=None
            row['relative_y_error_vs_ode']=None
        return row

    history.append(row_at(a,x,p,first_diag))
    progress_path=OUT/f'pm_{name}_progress.json'
    print(json.dumps({'case':name,'status':'started','nparticle':nparticle,'nmesh':nmesh,
                      'initial_force_diagnostic':initial_force_diagnostic}),flush=True)
    while a<.55*(1-1e-14):
        right,step=pm.select_step(x,p,force,a,bg,initial,nmesh,.55,halfstep)
        next_x,next_p,next_force,diag=pm.advance_kdk(x,p,a,right,bg,mesh,force)
        next_y=pm.core_scale(next_x,initial)
        if next_y**-3>=200:
            # KDK positions depend only on the first kick and drift. Root finding
            # therefore needs no new force FFT until the final partial step.
            def residual(n):
                nonlocal partial_event_calls
                partial_event_calls+=1
                trial,_=pm.predict_positions(x,p,force,a,float(np.exp(n)),bg)
                return pm.core_scale(trial,initial)-200**(-1/3)
            event_n=brentq(residual,np.log(a),np.log(right),xtol=2e-12,rtol=2e-12)
            right=float(np.exp(event_n))
            next_x,next_p,next_force,diag=pm.advance_kdk(x,p,a,right,bg,mesh,force)
            next_y=pm.core_scale(next_x,initial)
            step['dt_H0']=pm.time_factor(bg,a,right)
            step['delta_ln_a']=float(np.log(right/a))
            step['terminal_partial_KDK']=True
        dt=step['dt_H0']
        tdyn_left=pm.dynamical_time(a,history[-1]['y_median'],bg.omega_m0)
        tdyn_right=pm.dynamical_time(right,next_y,bg.omega_m0)
        step['actual_dt_over_tdyn']=dt/min(tdyn_left,tdyn_right)
        step['max_drift_cells']=diag['max_drift_cells']
        row=row_at(right,next_x,next_p,diag,step)
        previous=history[-1]
        if turnaround is None and previous['physical_expansion_ratio_median']>0>=row['physical_expansion_ratio_median']:
            # Recorded bracket and interpolation, not an unreported exact root.
            v0,v1=previous['physical_expansion_ratio_median'],row['physical_expansion_ratio_median']
            fraction=v0/(v0-v1)
            ta=float(np.exp(np.log(a)+fraction*np.log(right/a)))
            turnaround={'a':ta,'t_Gyr':float(backgrounds[label].evaluate_a(ta)['t_Gyr']),
                        'bracket_a':[float(a),float(right)],'estimator':'linear interpolation of physical radial expansion ratio in ln(a)',
                        'bracket_relative_width':float(right/a-1)}
        history.append(row)
        x,p,force,a=next_x,next_p,next_force,right
        nsteps+=1
        if step.get('terminal_partial_KDK'):
            event={key:row[key] for key in ('a','z','t_Gyr','y_median','Delta_proxy','nonhomology_fraction',
                'axis_ratio_min_max','actual_enclosed_mean_density_ratio','actual_enclosed_particle_count_at_proxy_radius')}
            event['estimator']='first fixed-label median(r/q)^(-3)=200, temporary partial KDK root'
            break
        if nsteps%25==0:
            progress={'case':name,'status':'running','steps':nsteps,'a':a,'Delta_proxy':row['Delta_proxy'],
                      'wall_seconds':time.perf_counter()-start,'latest':row}
            write_json(progress_path,progress)
            print(json.dumps({key:progress[key] for key in ('case','steps','a','Delta_proxy','wall_seconds')}),flush=True)
        if nsteps>10000:
            raise RuntimeError('Frozen bounded run exceeded 10000 adaptive steps')

    stop_a=event['a'] if event is not None else a
    common_end=min(stop_a,float(ode_event['a']))
    common=[r for r in history if r['a']<=common_end*(1+1e-13) and r['relative_y_error_vs_ode'] is not None]
    diagnostics={
        'finite':bool(np.all(np.isfinite(x)) and np.all(np.isfinite(p))),
        'max_mass_relative_error':max(r['mass_relative_error'] for r in history),
        'max_total_force_normalized':max(r['total_force_normalized'] for r in history),
        'max_total_force_absolute':max(r['total_force_absolute'] for r in history),
        'max_momentum_change_absolute':max(r['momentum_change_absolute'] for r in history),
        'max_drift_cells':max(r['max_drift_cells'] for r in history),
        'max_actual_dt_over_tdyn':max(r['actual_dt_over_tdyn'] for r in history),
        'max_nonhomology_fraction':max(r['nonhomology_fraction'] for r in history),
        'max_axis_deviation':max(r['axis_deviation'] for r in history),
        'max_radius_relative_error_vs_ode_common_pre_event':max(abs(r['relative_y_error_vs_ode']) for r in common),
        'common_pre_event_end_a':common_end,'common_pre_event_samples':len(common),
        'event_relative_error_vs_ode':None if event is None else event['a']/ode_event['a']-1,
        'partial_event_root_position_evaluations':partial_event_calls,
        'steps':nsteps,'maximum_step_halvings':max(r['step_halvings'] for r in history),
    }
    checks=[check('finite',diagnostics['finite']),check('event200_present',event is not None),
            check('mass',diagnostics['max_mass_relative_error'],1e-12),
            check('total_force',diagnostics['max_total_force_normalized'],1e-11),
            check('actual_drift_cells',diagnostics['max_drift_cells'],.2),
            check('actual_dt_over_tdyn',diagnostics['max_actual_dt_over_tdyn'],.05),
            check('nonhomology',diagnostics['max_nonhomology_fraction'],.1),
            check('axis_deviation',diagnostics['max_axis_deviation'],.1)]
    if name=='ref_main':
        checks += [check('PM_vs_ODE_event',diagnostics['event_relative_error_vs_ode'],.02),
                   check('PM_vs_ODE_radius',diagnostics['max_radius_relative_error_vs_ode_common_pre_event'],.03)]
    arrays=pm.final_profiles(x,p,a,bg,initial,128)
    for key in ('a','t_Gyr','y_median','Delta_proxy','nonhomology_fraction','axis_deviation','ode_y','relative_y_error_vs_ode'):
        arrays['trajectory_'+key]=np.array([np.nan if r.get(key) is None else r[key] for r in history])
    arrays['final_a']=np.array(a)
    arrays['center_box']=initial['center']
    array_path=OUT/f'pm_{name}.npz'
    np.savez_compressed(array_path,**arrays)
    state_path=None
    if name in ('ref_main','clock_main'):
        state_path=OUT/f'pm_{name}_final_state.npz'
        np.savez_compressed(state_path,positions=x,momenta=p,center_box=initial['center'],a=np.array(a))
    metadata={
        'case':name,'background':label,'nparticle':nparticle,'nmesh':nmesh,'halfstep':halfstep,
        'fft_workers':args.fft_workers,'a_initial':.02,'a_max':.55,'a_final':a,'box_mpc_h':10.,
        'initial':initial['metadata'],'protocol_sha256':FROZEN_PROTOCOL,
        'legacy_pm_sha256':pm.LEGACY_SHA256,'pm_refined_sha256':pm.sha(pm.__file__),
        'runner_sha256':pm.sha(__file__),'background_source_sha256':pm.sha(backgrounds[label].path),
        'spherical_calibration_sha256':pm.sha(OUT/'sphere_summary.json'),
        'integrator':'KDK with state-dependent predeclared step bounds; no strict symplectic/time-reversal claim',
        'radius_estimator':'median(r/q), fixed original labels .25qL<=q<=.75qL; density proxy, not uniform-density proof',
        'profile_units':'r and q in box units; density ratios to mean; 128^2 full-box z projection',
        'audit_scope':'full float64 final x,p plus aggregates' if state_path else 'saved trajectories, diagnostics and shell/projected aggregates; no final particle archive',
        'clock_scope':'one clock pair only; no convergence claim for tiny paired PM response',
    }
    result={'status':'completed','metadata':metadata,'event200':event,'turnaround':turnaround,
            'initial_force_diagnostic':initial_force_diagnostic,
            'ode_event200':ode_event,'diagnostics':diagnostics,'history':history,'checks':checks,
            'failures':[c for c in checks if not c['passed']],
            'array_sha256':pm.sha(array_path),'state_sha256':None if state_path is None else pm.sha(state_path),
            'wall_seconds':time.perf_counter()-start,'process_peak_rss_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    write_json(OUT/f'pm_{name}.json',result)
    write_json(progress_path,{'case':name,'status':'completed','steps':nsteps,'a':a,'wall_seconds':result['wall_seconds']})
    print(json.dumps({'case':name,'status':'completed','steps':nsteps,'event200_a':None if event is None else event['a'],
                      'failures':len(result['failures']),'wall_seconds':result['wall_seconds']}),flush=True)
    return result


def main(args):
    OUT.mkdir(parents=True,exist_ok=True)
    if args.preflight:
        preflight(args)
        return 0
    validation_path=OUT/'refined_pm_validation.json'
    if not validation_path.exists():
        raise RuntimeError('Run refined PM validation before production')
    validation=json.loads(validation_path.read_text())
    if validation.get('protocol_sha256')!=FROZEN_PROTOCOL:
        raise RuntimeError('Validation protocol mismatch')
    if validation.get('implementation_failures',validation.get('failures')):
        raise RuntimeError('Implementation validation failures require explicit review before production')
    _,backgrounds,pm_backgrounds,ode,delta=inputs()
    started=time.perf_counter()
    runs={}
    selected=[task for task in MATRIX if not args.case or task[0] in args.case]
    def perform(task):
        name=task[0]
        try:
            path=OUT/f'pm_{name}.json'
            if path.exists() and not args.resume:
                raise RuntimeError('Existing result; use --resume to reuse exact completed result, never overwrite')
            if path.exists() and args.resume:
                previous=json.loads(path.read_text())
                if (previous.get('status')=='completed' and previous['metadata']['protocol_sha256']==FROZEN_PROTOCOL
                        and previous['metadata']['pm_refined_sha256']==pm.sha(pm.__file__)
                        and previous['metadata']['runner_sha256']==pm.sha(__file__)):
                    print(json.dumps({'case':name,'status':'reused_matching_completed_result'}),flush=True)
                    return name,previous
                raise RuntimeError('Existing result does not match current frozen code/protocol; retained without overwrite')
            result=run_case(task,args,backgrounds,pm_backgrounds,ode,delta)
        except Exception as error:
            failure={'status':'execution_failed','case':name,'error':repr(error),'traceback':traceback.format_exc(),
                     'protocol_sha256':FROZEN_PROTOCOL,'created_utc':datetime.now(timezone.utc).isoformat()}
            write_json(OUT/f'pm_{name}_execution_failure.json',failure)
            result=failure
            print(json.dumps({'case':name,'status':'execution_failed','error':repr(error)}),flush=True)
        gc.collect()
        return name,result
    # Each case owns its files; only this coordinating thread writes the summary.
    with ThreadPoolExecutor(max_workers=args.parallel) as pool:
        for future in as_completed([pool.submit(perform,task) for task in selected]):
            name,result=future.result()
            runs[name]=result
    comparisons={}
    global_checks=[]
    main_run=runs.get('ref_main',{})
    if main_run.get('status')=='completed' and main_run['event200'] is not None:
        a0=main_run['event200']['a']
        for name,threshold in (('ref_halfstep',.005),('ref_particles128',.01),('ref_shift',.01)):
            other=runs.get(name,{})
            if other.get('status')=='completed' and other['event200'] is not None:
                value=other['event200']['a']/a0-1
                comparisons[name+'_event_relative_to_main']=value
                global_checks.append(check(name+'_event_vs_main',value,threshold))
        clock=runs.get('clock_main',{})
        if clock.get('status')=='completed' and clock['event200'] is not None:
            comparisons['single_clock_pair_fractional_a200_shift']=clock['event200']['a']/a0-1
            comparisons['single_clock_pair_delta_t_Myr']=1000*(clock['event200']['t_Gyr']-main_run['event200']['t_Gyr'])
    global_checks.append(check('all_six_cases_completed',all(runs.get(t[0],{}).get('status')=='completed' for t in MATRIX)))
    compact={}
    for name,result in runs.items():
        if result.get('status')=='completed':
            compact[name]={key:result[key] for key in ('status','event200','turnaround','ode_event200','initial_force_diagnostic','diagnostics','checks','failures','wall_seconds','array_sha256','state_sha256')}
            compact[name]['json_sha256']=pm.sha(OUT/f'pm_{name}.json')
        else:
            compact[name]=result
    summary={'created_utc':datetime.now(timezone.utc).isoformat(),'protocol_sha256':FROZEN_PROTOCOL,
             'code_sha256':{'pm_refined.py':pm.sha(pm.__file__),'run_spherical_pm.py':pm.sha(__file__)},
             'input_sha256':{'sphere_summary.json':pm.sha(OUT/'sphere_summary.json'),
                             'refined_pm_validation.json':pm.sha(validation_path)},
             'environment':{'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__},
             'runs':compact,'comparisons':comparisons,'checks':global_checks,
             'failures':[c for c in global_checks if not c['passed']],
             'case_failure_count':sum(len(r.get('failures',[]))+(r.get('status')!='completed') for r in runs.values()),
             'wall_seconds':time.perf_counter()-started,'concurrent_production_jobs':args.parallel,
             'scope':'Controlled collisionless collapse, no observed-data fit, virialization claim, halo catalogue, gas evolution or resolved tiny clock response'}
    write_json(OUT/'spherical_pm_summary.json',summary)
    print(json.dumps({'status':'matrix_finished','completed':sum(r.get('status')=='completed' for r in runs.values()),
                      'case_failure_count':summary['case_failure_count'],'comparison_failure_count':len(summary['failures'])}),flush=True)
    return 1 if summary['case_failure_count'] or summary['failures'] else 0


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fft-workers',type=int,default=4,choices=(1,2,3,4))
    parser.add_argument('--parallel',type=int,default=2,choices=(1,2))
    parser.add_argument('--preflight',action='store_true')
    parser.add_argument('--resume',action='store_true')
    parser.add_argument('--case',action='append',choices=[t[0] for t in MATRIX])
    raise SystemExit(main(parser.parse_args()))
