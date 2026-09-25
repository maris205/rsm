#!/usr/bin/env python3
"""Matched particle/force-mesh sphere diagnostic; exactly two frozen cases.

Derived from the archived six-case runner without altering the force kernel,
initial calibration, adaptive KDK, or event definition. Old files are read-only.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor,as_completed
from datetime import datetime,timezone
import gc
import hashlib
import importlib.util
import sys
import json
from pathlib import Path
import platform
import resource
import time
import traceback

# Archived modules are imported read-only, including their cache directories.
sys.dont_write_bytecode=True

import numpy as np
import scipy
from scipy.optimize import brentq

ROOT=Path(__file__).resolve().parents[1]
EXPERIMENTS=ROOT.parent
ARCHIVE=EXPERIMENTS/'halo_cooling_v01'
OUT=ROOT/'results'
FROZEN_PROTOCOL='4409835c2fd93047d5df7bcb9675b0771e72fc0bfdfb0dd69e8a89c54cae9bf0'
ARCHIVED_PROTOCOL='f3a0cc5abf21a4374439db756cbf8c4a3103bcabfa2bee6f9bd1921b28e608f2'
SOURCE_SHA256={'halo_cooling_v01/code/pm_refined.py': 'f8e5bcbc32d726a0345282e228458894ac3b5e75baf0051f65cdbfb05a268602', 'halo_cooling_v01/code/spherical_collapse.py': '4e0afba9b783575ec02ad60ec2f6acfaf87da8a0824f6a978b115e89ba522e47', 'halo_cooling_v01/code/run_spherical_pm.py': '4c8f29acc2598935a30b027691ccf3f01ae34b6b80a27aac0f7cf5349fafa28b', 'halo_cooling_v01/results/sphere_summary.json': 'eff75c2d2a7144effad867560d54fafc0a3bec3e344f493031e236639b5bc275', 'halo_cooling_v01/results/refined_pm_validation.json': '8f97dd3a305e99ca2fbfd7e7f79c128812b5d57e6d35c690793784de2d7781f6', 'cosmic_bridge_v01/code/pm.py': '1b2607ac0b7ba35c5df9ceb7375910ba75b7bed65a1edec79c0583f067a5d982', 'cosmic_bridge_v01/results/background_summary.json': '2bd9d4d15cb60432543bdb485816767e335bda4b47783cfbaa06b4c34a808307', 'cosmic_bridge_v01/results/background_eps0.csv': '1efdc1e2ef28dccd06cf9d0d0aa83396ab4db2a6a7dc9bc51b48b7ef6a1091a0', 'cosmic_bridge_v01/results/background_eps1e-4.csv': '2f6fc8ea8e1b24cb87b41129d30005030488da0bb20974625187a8936942ef13'}
MATRIX=(
    ('ref_matched64','reference',64,64,True,(0,0,0)),
    ('ref_matched128','reference',128,128,True,(0,0,0)),
)
EVENT_THRESHOLD=.02
RADIUS_THRESHOLD=.03
PAIR_EVENT_THRESHOLD=.01
SHAPE_THRESHOLD=.1
MAX_TOTAL_WALL_SECONDS=4*3600
MAX_STEPS=10000
MAX_ARCHIVE_BYTES=100_000_000


def verify_archived_sources():
    for relative,expected in SOURCE_SHA256.items():
        actual=hashlib.sha256((EXPERIMENTS/relative).read_bytes()).hexdigest()
        if actual!=expected:
            raise RuntimeError(f'Archived source/input differs: {relative}')


def import_archived(name,relative):
    spec=importlib.util.spec_from_file_location(name,EXPERIMENTS/relative)
    module=importlib.util.module_from_spec(spec)
    sys.modules[name]=module
    spec.loader.exec_module(module)
    return module


verify_archived_sources()
pm=import_archived('matched_archived_pm_refined','halo_cooling_v01/code/pm_refined.py')
sphere=import_archived('matched_archived_spherical_collapse','halo_cooling_v01/code/spherical_collapse.py')


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
    if FROZEN_PROTOCOL=='PENDING_PARENT_FREEZE' or not (ROOT/'protocol.md').is_file():
        raise RuntimeError('Parent has not frozen the protocol; production is prohibited')
    if pm.sha(ROOT/'protocol.md')!=FROZEN_PROTOCOL:
        raise RuntimeError('Protocol differs from frozen bytes')
    verify_archived_sources()
    summary=json.loads((ARCHIVE/'results/sphere_summary.json').read_text())
    if summary['protocol_sha256']!=ARCHIVED_PROTOCOL or summary['summary']['failures']:
        raise RuntimeError('A validated archived spherical ODE calibration is required')
    backgrounds={'reference':sphere.load_backgrounds()['reference']}
    delta=float(summary['delta_i'])
    initials=sphere.initial_conditions(delta,backgrounds)
    ode={label:sphere.integrate_sphere(bg,initials[label]) for label,bg in backgrounds.items()}
    pm_backgrounds={label:pm.Background.from_csv(bg.path,.315,label) for label,bg in backgrounds.items()}
    return summary,backgrounds,pm_backgrounds,ode,delta


def force_shape_diagnostics(positions,force,initial,initial_q4=None):
    vectors=pm.minimum_vector(positions[initial['core_mask']],initial['center'])
    radius=np.linalg.norm(vectors,axis=1)
    if np.any(radius<=0):
        raise FloatingPointError('Angular direction undefined at zero core radius')
    directions=vectors/radius[:,None]
    core_force=force[initial['core_mask']]
    radial=np.sum(core_force*directions,axis=1)
    transverse=core_force-radial[:,None]*directions
    radial_l2=float(np.linalg.norm(radial))
    transverse_l2=float(np.linalg.norm(transverse))
    q4=float(np.mean(np.sum(directions**4,axis=1))-.6)
    return {'angular_cubic_Q4':q4,
            'angular_cubic_Q4_difference_from_initial':0. if initial_q4 is None else q4-initial_q4,
            'core_radial_force_L2':radial_l2,'core_transverse_force_L2':transverse_l2,
            'core_transverse_to_radial_force_L2_ratio':None if radial_l2==0 else transverse_l2/radial_l2,
            'scope':'initial/final descriptive diagnostics, not fitted or thresholded'}


def corrected_final_profiles(positions,momenta,a,background,initial,analysis_mesh=128):
    arrays=pm.final_profiles(positions,momenta,a,background,initial,analysis_mesh)
    arrays['projection_axis_x_box']=np.arange(analysis_mesh,dtype=float)/analysis_mesh
    arrays['projection_axis_y_box']=np.arange(analysis_mesh,dtype=float)/analysis_mesh
    return arrays


def save_final_state(directory,name,positions,momenta,center,a):
    directory=Path(directory)
    directory.mkdir(parents=True,exist_ok=True)
    outputs={}
    for quantity,value in (('positions',positions),('momenta',momenta)):
        path=directory/f'pm_{name}_final_{quantity}.npz'
        if path.exists():
            raise RuntimeError(f'Refusing to overwrite final state: {path}')
        if value.dtype!=np.float64 or value.ndim!=2 or value.shape[1]!=3 or not np.all(np.isfinite(value)):
            raise ValueError('Finite float64 (particles,3) terminal state required')
        np.savez_compressed(path,**{quantity:value,'center_box':np.asarray(center,dtype=float),'a':np.array(a)})
        size=path.stat().st_size
        if size>=MAX_ARCHIVE_BYTES:
            raise RuntimeError(f'Archive exceeds frozen byte ceiling: {path} ({size})')
        outputs[quantity]={'file':path.name,'sha256':pm.sha(path),'bytes':size,
                           'shape':list(value.shape),'dtype':str(value.dtype)}
    return outputs


def pair_comparison(runs):
    low,high=runs.get('ref_matched64',{}),runs.get('ref_matched128',{})
    available=all(r.get('status')=='completed' and r.get('event200') is not None for r in (low,high))
    difference=high['event200']['a']/low['event200']['a']-1 if available else None
    return {'matched128_event_relative_to_matched64':difference},[
        check('matched128_event_vs_matched64',difference,PAIR_EVENT_THRESHOLD),
        check('both_cases_completed',all(runs.get(t[0],{}).get('status')=='completed' for t in MATRIX))]


def locate_event_a(x,p,g,a_left,a_right,background,initial,target=200.):
    """Original partial-KDK position root; return scale factor and call count."""
    calls=0
    def residual(n):
        nonlocal calls
        calls+=1
        trial,_=pm.predict_positions(x,p,g,a_left,float(np.exp(n)),background)
        return pm.core_scale(trial,initial)-target**(-1/3)
    event_n=brentq(residual,np.log(a_left),np.log(a_right),xtol=2e-12,rtol=2e-12)
    return float(np.exp(event_n)),calls


def run_case(task,args,backgrounds,pm_backgrounds,ode,delta,deadline):
    name,label,nparticle,nmesh,halfstep,shift=task
    start=time.perf_counter()
    started_utc=datetime.now(timezone.utc).isoformat()
    initial=pm.compensated_initial(nparticle,delta,backgrounds['reference'],shift,nmesh)
    bg=pm_backgrounds[label]
    x,p=initial['positions'].copy(),initial['momenta'].copy()
    mesh=pm.PMGrid(nmesh,args.fft_workers)
    force_started=time.perf_counter()
    force,first_diag=mesh.force(x)
    initial_force_seconds=time.perf_counter()-force_started
    shape_initial=force_shape_diagnostics(x,force,initial)
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
    resource_stop_reason=None
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
        if time.perf_counter()>=deadline or nsteps>=MAX_STEPS:
            resource_stop_reason='total_wall_limit' if time.perf_counter()>=deadline else 'per_case_step_limit'
            break
        right,step=pm.select_step(x,p,force,a,bg,initial,nmesh,.55,halfstep)
        next_x,next_p,next_force,diag=pm.advance_kdk(x,p,a,right,bg,mesh,force)
        next_y=pm.core_scale(next_x,initial)
        if next_y**-3>=200:
            # KDK positions depend only on the first kick and drift. Root finding
            # therefore needs no new force FFT until the final partial step.
            right,calls=locate_event_a(x,p,force,a,right,bg,initial)
            partial_event_calls+=calls
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
        'resource_stop_reason':resource_stop_reason,
    }
    checks=[check('finite',diagnostics['finite']),check('event200_present',event is not None),
            check('mass',diagnostics['max_mass_relative_error'],1e-12),
            check('total_force',diagnostics['max_total_force_normalized'],1e-11),
            check('actual_drift_cells',diagnostics['max_drift_cells'],.2),
            check('actual_dt_over_tdyn',diagnostics['max_actual_dt_over_tdyn'],.05),
            check('nonhomology',diagnostics['max_nonhomology_fraction'],SHAPE_THRESHOLD),
            check('axis_deviation',diagnostics['max_axis_deviation'],SHAPE_THRESHOLD),
            check('resource_limits_not_reached',resource_stop_reason is None)]
    checks += [check('PM_vs_ODE_event',diagnostics['event_relative_error_vs_ode'],EVENT_THRESHOLD),
               check('PM_vs_ODE_radius',diagnostics['max_radius_relative_error_vs_ode_common_pre_event'],RADIUS_THRESHOLD)]
    shape_final=force_shape_diagnostics(x,force,initial,shape_initial['angular_cubic_Q4'])
    arrays=corrected_final_profiles(x,p,a,bg,initial,128)
    for key in ('a','t_Gyr','y_median','Delta_proxy','nonhomology_fraction','axis_deviation','ode_y','relative_y_error_vs_ode'):
        arrays['trajectory_'+key]=np.array([np.nan if r.get(key) is None else r[key] for r in history])
    arrays['final_a']=np.array(a)
    arrays['center_box']=initial['center']
    array_path=OUT/f'pm_{name}.npz'
    np.savez_compressed(array_path,**arrays)
    state_archives=save_final_state(OUT,name,x,p,initial['center'],a)
    metadata={
        'case':name,'background':label,'nparticle':nparticle,'nmesh':nmesh,'halfstep':halfstep,
        'fft_workers':args.fft_workers,'a_initial':.02,'a_max':.55,'a_final':a,'box_mpc_h':10.,
        'initial':initial['metadata'],'protocol_sha256':FROZEN_PROTOCOL,
        'legacy_pm_sha256':pm.LEGACY_SHA256,'pm_refined_sha256':pm.sha(pm.__file__),
        'runner_sha256':pm.sha(__file__),'background_source_sha256':pm.sha(backgrounds[label].path),
        'spherical_calibration_sha256':pm.sha(ARCHIVE/'results/sphere_summary.json'),
        'integrator':'KDK with state-dependent predeclared step bounds; no strict symplectic/time-reversal claim',
        'radius_estimator':'median(r/q), fixed original labels .25qL<=q<=.75qL; density proxy, not uniform-density proof',
        'profile_units':'r and q in box units; density ratios to mean; 128^2 full-box z projection',
        'audit_scope':'both cases: full float64 final x,p in separate archives plus all saved trajectories and aggregates',
        'projection_coordinate_convention':'CIC nodal samples x_j=j/128; archived old half-cell display coordinates corrected only in new products',
        'resolution_scope':'only two matched particle/force scales; independent phase/time convergence is not established',
        'clock_scope':'epsilon=0 only; no new clock response calculation',
        'angular_scope':'Q4 and transverse/radial force diagnostics only at initial/final state, descriptive without acceptance threshold',
        'initial_force_evaluation_scope':'one production initialization, reused for first KDK; independent initial-field diagnostic is separately counted and no runner preflight is performed',
        'source_sha256':SOURCE_SHA256,
    }
    status='completed' if resource_stop_reason is None else 'stopped_resource_limit'
    result={'status':status,'metadata':metadata,'event200':event,'turnaround':turnaround,
            'initial_force_diagnostic':initial_force_diagnostic,'initial_force_seconds':initial_force_seconds,
            'force_shape_diagnostics':{'initial':shape_initial,'final':shape_final},
            'ode_event200':ode_event,'diagnostics':diagnostics,'history':history,'checks':checks,
            'failures':[c for c in checks if not c['passed']],
            'array_sha256':pm.sha(array_path),'state_archives':state_archives,
            'started_utc':started_utc,'ended_utc':datetime.now(timezone.utc).isoformat(),
            'wall_seconds':time.perf_counter()-start,'process_peak_rss_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    write_json(OUT/f'pm_{name}.json',result)
    write_json(progress_path,{'case':name,'status':status,'steps':nsteps,'a':a,'wall_seconds':result['wall_seconds']})
    print(json.dumps({'case':name,'status':status,'steps':nsteps,'event200_a':None if event is None else event['a'],
                      'failures':len(result['failures']),'wall_seconds':result['wall_seconds']}),flush=True)
    return result


def main(args):
    OUT.mkdir(parents=True,exist_ok=True)
    _,backgrounds,pm_backgrounds,ode,delta=inputs()
    validation_path=OUT/'matched_runner_validation.json'
    if not validation_path.exists():
        raise RuntimeError('Run new matched-runner validation before production')
    validation=json.loads(validation_path.read_text())
    if not (validation.get('passed') is True and validation.get('runner_sha256')==pm.sha(__file__)
            and validation.get('protocol_sha256')==FROZEN_PROTOCOL):
        raise RuntimeError('Matching new runner/protocol implementation validation is required')
    # Prevent partial overwrite/retry. Both cases execute exactly once in a clean results directory.
    for task in MATRIX:
        if list(OUT.glob(f'pm_{task[0]}*')):
            raise RuntimeError(f'Existing production artifact for {task[0]}; no implicit retry or overwrite')
    started=time.perf_counter()
    started_utc=datetime.now(timezone.utc).isoformat()
    deadline=started+MAX_TOTAL_WALL_SECONDS
    runs={}
    def perform(task):
        name=task[0]
        try:
            result=run_case(task,args,backgrounds,pm_backgrounds,ode,delta,deadline)
        except Exception as error:
            failure={'status':'execution_failed','case':name,'error':repr(error),'traceback':traceback.format_exc(),
                     'protocol_sha256':FROZEN_PROTOCOL,'created_utc':datetime.now(timezone.utc).isoformat()}
            write_json(OUT/f'pm_{name}_execution_failure.json',failure)
            result=failure
            print(json.dumps({'case':name,'status':'execution_failed','error':repr(error)}),flush=True)
        gc.collect()
        return name,result
    with ThreadPoolExecutor(max_workers=args.parallel) as pool:
        for future in as_completed([pool.submit(perform,task) for task in MATRIX]):
            name,result=future.result()
            runs[name]=result
    comparisons,global_checks=pair_comparison(runs)
    compact={}
    for name,result in runs.items():
        if result.get('status')!='execution_failed':
            compact[name]={key:result[key] for key in ('status','event200','turnaround','ode_event200',
                'initial_force_diagnostic','initial_force_seconds','force_shape_diagnostics','diagnostics',
                'checks','failures','wall_seconds','array_sha256','state_archives')}
            compact[name]['json_sha256']=pm.sha(OUT/f'pm_{name}.json')
        else:
            compact[name]=result
    summary={'created_utc':datetime.now(timezone.utc).isoformat(),'protocol_sha256':FROZEN_PROTOCOL,
             'started_utc':started_utc,'command_argv':sys.argv,
             'code_sha256':{'run_matched_pm.py':pm.sha(__file__)},'archived_source_sha256':SOURCE_SHA256,
             'input_sha256':{'matched_runner_validation.json':pm.sha(validation_path)},
             'environment':{'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__},
             'runs':compact,'comparisons':comparisons,'checks':global_checks,
             'failures':[c for c in global_checks if not c['passed']],
             'case_failure_count':sum(len(r.get('failures',[]))+(r.get('status')=='execution_failed') for r in runs.values()),
             'wall_seconds':time.perf_counter()-started,'concurrent_production_jobs':args.parallel,
             'fft_workers_per_case':args.fft_workers,
             'resource_policy':{'maximum_wall_seconds':MAX_TOTAL_WALL_SECONDS,'maximum_steps_per_case':MAX_STEPS,
                 'enforcement':'checked between accepted steps; terminal archive may add a small bounded overrun'},
             'scope':'Two matched-scale epsilon=0 sphere diagnostics; no independent phase/time convergence, new clock response, halo structure or gas evolution claim'}
    write_json(OUT/'matched_pm_summary.json',summary)
    print(json.dumps({'status':'matrix_finished','completed':sum(r.get('status')=='completed' for r in runs.values()),
                      'case_failure_count':summary['case_failure_count'],'comparison_failure_count':len(summary['failures'])}),flush=True)
    return 1 if summary['case_failure_count'] or summary['failures'] else 0


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fft-workers',type=int,default=4,choices=(1,2,3,4))
    parser.add_argument('--parallel',type=int,default=2,choices=(1,2))
    raise SystemExit(main(parser.parse_args()))
