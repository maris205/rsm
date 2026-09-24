#!/usr/bin/env python3
"""Frozen matched-IC PM matrix; no fit to observed structure or CMB."""
from __future__ import annotations
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import hashlib
import json
import platform
import time
from pathlib import Path
import numpy as np
import scipy
from pm import Background, generate_initial_modes, initial_conditions, integrate

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results'

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def write_json(path, obj):
    def convert(x):
        if isinstance(x, np.ndarray): return x.tolist()
        if isinstance(x, np.generic): return x.item()
        if isinstance(x, Path): return str(x)
        raise TypeError(type(x).__name__)
    path.write_text(json.dumps(obj, indent=2, default=convert, allow_nan=False) + '\n')

def run(args):
    OUT.mkdir(parents=True, exist_ok=True)
    ref = Background.from_csv(args.reference, omega_m0=.315, label='reference')
    clock = Background.from_csv(args.clock, omega_m0=.315, label='clock_eps1e-4')
    modes = generate_initial_modes(args.spectrum, seed=20260924, base_grid=32,
                                   box_mpc_h=50., cutoff_fraction=.6)
    np.savez_compressed(OUT / 'initial_modes.npz', coefficients=modes['coefficients'])
    edges = np.arange(.5, 16.6, 1.) * (2*np.pi/50.)
    runs = {}
    start = time.perf_counter()
    def perform(task):
        nmesh, steps, label, bg, initial = task
        name = f'{label}_n{nmesh}_s{steps}'
        print(f'Start {name}', flush=True)
        run_started = time.perf_counter()
        result = integrate(initial, bg, nmesh, steps=steps, a_final=1.,
                           a_outputs=(.02,.1,.25,.5,1.), box_mpc_h=50.,
                           bin_edges=edges, retain_particles=False)
        snapshots = result['snapshots']
        arrays = {'a':np.array([s['a'] for s in snapshots]),
                  'delta':np.array([s['delta'] for s in snapshots],dtype=np.float32),
                  'projection':np.array([s['projection'] for s in snapshots],dtype=np.float32),
                  'k':snapshots[-1]['spectrum']['k'],
                  'pk':np.array([s['spectrum']['pk'] for s in snapshots]),
                  'mode_count':snapshots[-1]['spectrum']['count'], 'edges':edges}
        np.savez_compressed(OUT / f'{name}.npz', **arrays)
        info = {k:result[k] for k in ('metadata','diagnostics','history')}
        info['integration_wall_seconds'] = time.perf_counter()-run_started
        info.update(initial_metadata=initial['metadata'], snapshot_a=arrays['a'],
                    snapshot_step_kind=[s['output_step_kind'] for s in snapshots],
                    array_sha256=sha(OUT/f'{name}.npz'))
        write_json(OUT / f'{name}.json', info)
        record = {'pk_final':arrays['pk'][-1], 'k':arrays['k'],
                      'variance_final':float(np.var(snapshots[-1]['delta'])),
                      'diagnostics':result['diagnostics'],
                      'integration_wall_seconds':info['integration_wall_seconds'],
                      'npz_sha256':sha(OUT/f'{name}.npz'),
                      'json_sha256':sha(OUT/f'{name}.json')}
        print(f'Finished {name}', flush=True)
        return name, record
    initials = {n:initial_conditions(modes, n, a_initial=.02,
                        E_reference=float(ref.E(.02)), f_reference=args.f_initial)
                for n in (32,64)}
    tasks=[(n, steps, label, bg, initials[n])
           for n,steps in [(64,256),(64,128),(64,512),(32,256)]
           for label,bg in [('reference',ref),('clock',clock)]]
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        for future in as_completed([pool.submit(perform, task) for task in tasks]):
            name, record = future.result()
            runs[name] = record
    runs=dict(sorted(runs.items()))
    # Compare paired responses separately from the absolute spectrum.
    ratios = {}
    for nmesh,steps in [(64,256),(64,128),(64,512),(32,256)]:
        key=f'n{nmesh}_s{steps}'
        ratios[key]=runs[f'clock_{key}']['pk_final']/runs[f'reference_{key}']['pk_final']-1
    k=runs['reference_n64_s256']['k']
    low=(k<=.3)&np.isfinite(k)
    main=ratios['n64_s256']
    comparison = {
        'k_h_mpc':k, 'low_k_mask':low, 'paired_fractional_power_response':ratios,
        'absolute_reference_pk_fractional_change_256_to_512':
            runs['reference_n64_s256']['pk_final']/runs['reference_n64_s512']['pk_final']-1,
        'absolute_reference_pk_fractional_change_32_to_64':
            runs['reference_n32_s256']['pk_final']/runs['reference_n64_s256']['pk_final']-1,
        'paired_response_time_difference_256_minus_512':main-ratios['n64_s512'],
        'paired_response_space_difference_32_minus_64':ratios['n32_s256']-main,
        'low_k_response_exceeds_tested_time_and_space_differences':bool(np.all(
            (np.abs(main[low])>np.abs((main-ratios['n64_s512'])[low])) &
            (np.abs(main[low])>np.abs((ratios['n32_s256']-main)[low])))),
        'interpretation':'Differences at tested grids/steps are diagnostics, not certified error bounds or detection significance.'}
    # Units: rho_crit0 = 2.77536627e11 h^2 Msun/Mpc^3.
    resolutions = {str(n):{'particle_count':n**3,'mesh_mpc_over_h':50./n,
                          'particle_mass_msun_over_h':2.77536627e11*.315*50.**3/n**3,
                          '100_particle_mass_msun_over_h':100*2.77536627e11*.315*50.**3/n**3}
                   for n in (32,64)}
    checks=[]
    for name,result in runs.items():
        diag=result['diagnostics']
        checks += [dict(name=name+'_finite',value=diag['finite'],passed=bool(diag['finite'])),
                   dict(name=name+'_mass',value=diag['max_mass_relative_error'],threshold=1e-12,
                        passed=diag['max_mass_relative_error']<1e-12),
                   dict(name=name+'_total_force',value=diag['max_total_force_normalized'],threshold=1e-12,
                        passed=diag['max_total_force_normalized']<1e-12)]
    min_horizon_ratio = min(float(np.min(299792.458*(2*np.pi/50.)/(100*bg.a*bg.e)))
                            for bg in (ref,clock))
    manifest = {'protocol_sha256':sha(ROOT/'protocol.md'),
        'inputs':{str(Path(p).relative_to(ROOT)) if Path(p).is_relative_to(ROOT) else Path(p).name:sha(p)
                  for p in (args.reference,args.clock,args.spectrum)},
        'code_sha256':{p.name:sha(p) for p in (Path(__file__),ROOT/'code/pm.py')},
        'environment':{'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__},
        'box_mpc_over_h':50.,'f_reference_initial':args.f_initial,
        'common_p_initial_coefficient':.02**2*float(ref.E(.02))*args.f_initial,
        'seed':20260924,'initial_modes_metadata':modes['metadata'],
        'minimum_kfund_over_aHc':min_horizon_ratio,
        'initial_modes_sha256':sha(OUT/'initial_modes.npz'),
        'resolutions':resolutions,'runs':runs,'comparison':comparison,'checks':checks,
        'wall_seconds':time.perf_counter()-start,'concurrent_workers':args.workers,
        'timing_scope':'Concurrent execution timing, not a serial code-performance benchmark',
        'scope':'Matched late-time IC experiment; external CAMB spectrum; no scalar early transfer, gas, stars, halos or observed-data fit.'}
    write_json(OUT/'clustering_summary.json',manifest)
    print(json.dumps({'wall_seconds':manifest['wall_seconds'],
                      'low_k_paired_response':main[low].tolist(),
                      'low_k_resolution_diagnostic':comparison['low_k_response_exceeds_tested_time_and_space_differences']}),flush=True)

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--reference',type=Path,required=True)
    ap.add_argument('--clock',type=Path,required=True)
    ap.add_argument('--spectrum',type=Path,required=True)
    ap.add_argument('--f-initial',type=float,required=True)
    ap.add_argument('--workers',type=int,default=4,help='Independent runs in parallel; not a timing benchmark')
    run(ap.parse_args())
