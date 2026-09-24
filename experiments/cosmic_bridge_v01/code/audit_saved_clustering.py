#!/usr/bin/env python3
"""Audit saved PM products without importing any production PM functions.

Uses NumPy FFT, direct shell sums, and an independently written CIC deposit.
No particle trajectories are integrated. Partial records are explicitly marked.
"""
from __future__ import annotations
from datetime import datetime, timezone
from pathlib import Path
import hashlib
import json
import platform
import re
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / 'results'
REPORT = ROOT / 'reports/clustering_independent_audit.json'
EXPECTED = [f'{label}_n{n}_s{steps}' for n,steps in
            ((64,256),(64,128),(64,512),(32,256)) for label in ('reference','clock')]
BOX = 50.


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def raw_sha(array):
    return hashlib.sha256(np.ascontiguousarray(array).view(np.uint8)).hexdigest()


def independent_spectrum(delta, edges):
    n = delta.shape[0]
    integers = np.arange(n)
    integers = np.where(integers < (n+1)//2, integers, integers-n)
    k = (2*np.pi/BOX)*np.sqrt(integers[:,None,None]**2+
                              integers[None,:,None]**2+integers[None,None,:]**2)
    power = np.abs(np.fft.fftn(np.asarray(delta,dtype=float))/n**3)**2*BOX**3
    counts, means, centers = [], [], []
    for lo,hi in zip(edges[:-1],edges[1:]):
        mask=(k>=lo)&(k<hi)
        counts.append(int(mask.sum()))
        means.append(float(power[mask].mean()))
        centers.append(float(k[mask].mean()))
    return np.array(centers), np.array(means), np.array(counts)


def independent_initial(coefficients, n, factor):
    base=coefficients.shape[0]
    integers=np.arange(base)
    integers=np.where(integers<(base+1)//2,integers,integers-base)
    mode=np.stack(np.meshgrid(integers,integers,integers,indexing='ij'),axis=-1)
    use=np.any(mode!=0,axis=-1)&(np.abs(coefficients)>0)
    m=mode[use]
    k=2*np.pi*m
    phase=np.exp(2j*np.pi*np.sum(m,axis=1)*(.5/n))
    gradient=1j*k/np.sum(k*k,axis=1)[:,None]
    displacement=np.empty((n,n,n,3))
    for axis in range(3):
        modes=np.zeros((n,n,n),dtype=complex)
        modes[tuple((m%n).T)]=coefficients[use]*phase*gradient[:,axis]
        displacement[...,axis]=(np.fft.ifftn(modes)*n**3).real
    displacement=displacement.reshape(-1,3)
    centers=(np.arange(n)+.5)/n
    lag=np.stack(np.meshgrid(centers,centers,centers,indexing='ij'),axis=-1).reshape(-1,3)
    positions=(lag+displacement)%1
    scaled=positions*n
    lower=np.floor(scaled).astype(int)
    part=scaled-lower
    counts=np.zeros((n,n,n),dtype=float)
    for corner in np.ndindex(2,2,2):
        side=np.asarray(corner)
        cells=(lower+side)%n
        weights=np.prod(np.where(side,part,1-part),axis=1)
        np.add.at(counts,tuple(cells.T),weights)
    # Exactly one particle per cell: each count is rho/rho_bar.
    return counts-1, displacement, factor*displacement


def run():
    out={'created_utc':datetime.now(timezone.utc).isoformat(),
         'scope':'Read saved arrays and metadata; independent NumPy FFT/CIC reconstruction; '
                 'no PM integrations and no production pm.py imports.',
         'checks':[],'runs':{},'comparisons':{},'input_sha256':{},
         'software':{'python':platform.python_version(),'numpy':np.__version__},
         'reviewed_production_code_sha256':{p.name:sha(p) for p in
                (ROOT/'code/pm.py',ROOT/'code/run_clustering.py')},
         'audit_code_sha256':sha(__file__)}

    def check(name, value, threshold=0., **context):
        v=float(value)
        out['checks'].append({'name':name,'value':v,'threshold':threshold,
                              'passed':bool(np.isfinite(v) and abs(v)<=threshold),**context})

    def relative(a,b):
        return float(np.max(abs(np.asarray(a)/np.asarray(b)-1)))

    with np.load(RESULTS/'initial_modes.npz') as archive:
        coefficients=archive['coefficients']
    coefficient_sha=raw_sha(coefficients)
    bg=json.loads((RESULTS/'background_summary.json').read_text())
    ref=bg['models']['eps0']
    ef=ref['initial_pm']['E']
    ff=ref['regular_f_PM']
    factor=.02**2*ef*ff
    check('source_coefficients_finite',0 if np.isfinite(coefficients).all() else 1)
    check('source_coefficients_shape',0 if coefficients.shape==(32,32,32) else 1)
    check('source_DC_zero',abs(coefficients[0,0,0]))
    # Recreate the declared source realization independently with NumPy FFT.
    data=np.genfromtxt(ROOT/'inputs/reference_linear_spectrum.csv',delimiter=',',names=True)
    modes=np.fft.fftfreq(32)*32
    kmag=(2*np.pi/BOX)*np.sqrt(modes[:,None,None]**2+modes[None,:,None]**2+modes[None,None,:]**2)
    keep=(kmag>0)&(kmag<=.6*np.pi*32/BOX)
    expected_power=np.zeros((32,32,32))
    expected_power[keep]=np.exp(np.interp(np.log(kmag[keep]),np.log(data['k_h_mpc']),
                                         np.log(data['pk_z49_mpc_over_h3'])))
    white=np.random.default_rng(20260924).normal(size=(32,32,32))
    reconstructed_coeff=np.fft.fftn(white)*np.sqrt(expected_power/(BOX**3*32**3))
    check('source_coefficients_independent_relative_L2',
          np.linalg.norm(reconstructed_coeff-coefficients)/np.linalg.norm(coefficients),1e-13)
    initial_cache={}
    stored={}
    found=[]
    for name in EXPECTED:
        jp,ap=RESULTS/(name+'.json'),RESULTS/(name+'.npz')
        if not jp.exists() or not ap.exists():
            continue
        info=json.loads(jp.read_text())
        with np.load(ap) as archive:
            arrays={key:archive[key] for key in archive.files}
        found.append(name)
        stored[name]=(info,arrays)
        n,steps=map(int,re.search(r'_n(\d+)_s(\d+)$',name).groups())
        meta,ic,hist=info['metadata'],info['initial_metadata'],info['history']
        prefix=name+':'
        check(prefix+'array_file_hash',0 if sha(ap)==info['array_sha256'] else 1)
        check(prefix+'pm_code_hash',0 if meta['code_sha256']==out['reviewed_production_code_sha256']['pm.py'] else 1)
        check(prefix+'metadata_dimensions',0 if (meta['nmesh'],meta['steps'],meta['particles'])==(n,steps,n**3) else 1)
        check(prefix+'initial_metadata_copy',0 if meta['initial_conditions']==ic else 1)
        check(prefix+'source_coefficient_hash',0 if ic['coefficient_sha256']==coefficient_sha else 1)
        check(prefix+'source_spectrum_hash',0 if ic['source_sha256']==sha(ROOT/'inputs/reference_linear_spectrum.csv') else 1)
        bgfile=RESULTS/('background_eps0.csv' if name.startswith('reference_') else 'background_eps1e-4.csv')
        check(prefix+'background_hash',0 if meta['background']['source_sha256']==sha(bgfile) else 1)
        check(prefix+'reference_E_initial',relative(ic['E_reference_initial'],ef),1e-13)
        check(prefix+'reference_f_initial',relative(ic['f_reference_initial'],ff),1e-13)
        check(prefix+'common_momentum_factor',relative(ic['common_momentum_factor'],factor),1e-13)
        check(prefix+'snapshot_shapes',0 if arrays['delta'].shape==(5,n,n,n) and
              arrays['projection'].shape==(5,n,n) and arrays['pk'].shape==(5,16) else 1)
        check(prefix+'all_saved_arrays_finite',0 if all(np.isfinite(v).all() for v in arrays.values()) else 1)
        check(prefix+'snapshot_a',np.max(abs(arrays['a']-np.array([.02,.1,.25,.5,1.]))),1e-14)
        check(prefix+'saved_density_lower_bound',max(0.,-1.-float(arrays['delta'].min())),2e-7)
        check(prefix+'density_raw_mean_from_float32',np.max(abs(arrays['delta'].mean(axis=(1,2,3),dtype=np.float64))),1e-6)
        delta=np.asarray(arrays['delta'],dtype=float)
        reconstructed_projection=1+delta.mean(axis=3)
        projection_error=np.max(abs(reconstructed_projection-arrays['projection'])/
                                np.maximum(1.,abs(arrays['projection'])))
        check(prefix+'projection_from_saved_float32_density',projection_error,1e-6)
        k,power,count=independent_spectrum(delta[-1],arrays['edges'])
        power_error=relative(power,arrays['pk'][-1])
        check(prefix+'final_power_independent_FFT',power_error,1e-6,
              note='Tolerance includes float32 density quantization; saved power used float64 density.')
        check(prefix+'shell_count',np.max(abs(count-arrays['mode_count'])))
        check(prefix+'shell_centers',np.max(abs(k-arrays['k'])),1e-11)
        if n not in initial_cache:
            initial_cache[n]=independent_initial(coefficients,n,factor)
        initial_delta,s,p0=initial_cache[n]
        check(prefix+'initial_density_independent_CIC',
              np.max(abs(initial_delta-delta[0]))/max(1.,np.max(abs(initial_delta))),1e-7)
        check(prefix+'initial_displacement_rms',relative(np.sqrt(np.mean(s*s)),ic['displacement_rms_box_units']),1e-13)
        check(prefix+'initial_momentum_mean',np.max(abs(p0.mean(axis=0)-hist[0]['mean_momentum'])),1e-15)
        check(prefix+'initial_momentum_rms',relative(np.sqrt(np.mean(p0*p0)),hist[0]['momentum_rms']),1e-13)
        check(prefix+'history_length',len(hist)-(steps+1))
        check(prefix+'history_fixed_log_a',np.max(abs(np.array([row['a'] for row in hist])-
                                                    np.geomspace(.02,1,steps+1))),1e-13)
        pmean=np.asarray([row['mean_momentum'] for row in hist])
        pchange=np.linalg.norm(pmean-pmean[0],axis=1)
        check(prefix+'recorded_mean_momentum_change',
              np.max(abs(pchange-np.array([row['momentum_change_absolute'] for row in hist]))),1e-28)
        check(prefix+'recorded_force_normalization',
              np.max(abs(np.array([row['total_force_absolute']/max(row['mean_force_magnitude'],1e-300)-
                                    row['total_force_normalized'] for row in hist]))),1e-25)
        mapping={'max_mass_relative_error':'mass_relative_error',
                 'max_total_force_normalized':'total_force_normalized',
                 'max_total_force_absolute':'total_force_absolute',
                 'max_momentum_change_absolute':'momentum_change_absolute',
                 'max_drift_cells':'max_drift_cells',
                 'max_force_ifft_imaginary':'force_ifft_max_imaginary'}
        for key,field in mapping.items():
            check(prefix+'summary_'+key,info['diagnostics'][key]-max(row[field] for row in hist),1e-25)
        check(prefix+'mass_protocol_threshold',info['diagnostics']['max_mass_relative_error'],1e-12)
        check(prefix+'force_protocol_threshold',info['diagnostics']['max_total_force_normalized'],1e-12)
        check(prefix+'reported_final_finite',0 if info['diagnostics']['finite'] else 1)
        out['runs'][name]={'final_power_max_relative_quantization_difference':power_error,
                           'projection_max_scaled_difference':float(projection_error),
                           'independent_final_pk':power.tolist(),
                           'saved_final_pk':arrays['pk'][-1].tolist(),
                           'max_momentum_change_absolute':info['diagnostics']['max_momentum_change_absolute'],
                           'max_drift_cells':info['diagnostics']['max_drift_cells']}
        out['input_sha256'][str(jp.relative_to(ROOT))]=sha(jp)
        out['input_sha256'][str(ap.relative_to(ROOT))]=sha(ap)
    # Every same-grid run must share exactly the saved initial realization.
    for n in (32,64):
        names=[name for name in found if f'_n{n}_' in name]
        if not names:
            continue
        first=stored[names[0]][1]
        for name in names[1:]:
            later=stored[name][1]
            check(name+':identical_same_grid_initial_density',
                  0 if np.array_equal(first['delta'][0],later['delta'][0]) else 1)
            check(name+':identical_same_grid_initial_power',
                  0 if np.array_equal(first['pk'][0],later['pk'][0]) else 1)
    pair={}
    for n,steps in ((64,256),(64,128),(64,512),(32,256)):
        suffix=f'n{n}_s{steps}'
        if all(f'{label}_{suffix}' in stored for label in ('clock','reference')):
            pair[suffix]=stored['clock_'+suffix][1]['pk'][-1]/stored['reference_'+suffix][1]['pk'][-1]-1
    out['comparisons']['paired_fractional_power_response']={key:val.tolist() for key,val in pair.items()}
    summary_path=RESULTS/'clustering_summary.json'
    complete=all(name in found for name in EXPECTED) and summary_path.exists()
    if complete:
        summary=json.loads(summary_path.read_text())
        cmp=summary['comparison']
        for key,value in pair.items():
            check('summary_pair_'+key,np.max(abs(value-cmp['paired_fractional_power_response'][key])),1e-14)
        main=pair['n64_s256']
        time_difference=main-pair['n64_s512']
        space_difference=pair['n32_s256']-main
        absolute_time=stored['reference_n64_s256'][1]['pk'][-1]/stored['reference_n64_s512'][1]['pk'][-1]-1
        absolute_space=stored['reference_n32_s256'][1]['pk'][-1]/stored['reference_n64_s256'][1]['pk'][-1]-1
        for key,value in (('paired_response_time_difference_256_minus_512',time_difference),
                          ('paired_response_space_difference_32_minus_64',space_difference),
                          ('absolute_reference_pk_fractional_change_256_to_512',absolute_time),
                          ('absolute_reference_pk_fractional_change_32_to_64',absolute_space)):
            check('summary_'+key,np.max(abs(value-cmp[key])),1e-14)
        k=stored['reference_n64_s256'][1]['k']
        low=(k<=.3)&np.isfinite(k)
        resolved=bool(np.all((abs(main[low])>abs(time_difference[low])) &
                             (abs(main[low])>abs(space_difference[low]))))
        check('summary_low_k_mask',0 if np.array_equal(low,cmp['low_k_mask']) else 1)
        check('summary_resolution_flag',0 if resolved==cmp['low_k_response_exceeds_tested_time_and_space_differences'] else 1)
        out['comparisons'].update(low_k_shell_mean_h_mpc=k[low].tolist(),
              low_k_response=main[low].tolist(),low_k_time_difference=time_difference[low].tolist(),
              low_k_space_difference=space_difference[low].tolist(),
              low_k_absolute_time_difference=absolute_time[low].tolist(),
              low_k_absolute_space_difference=absolute_space[low].tolist(),
              low_k_resolution_flag=resolved,
              bin_note='Low-k selection uses shell means <=0.3 h/Mpc; the second shell extends to 0.3078 h/Mpc.')
        for name in EXPECTED:
            info,arrays=stored[name]
            srun=summary['runs'][name]
            for ext in ('json','npz'):
                check('summary_hash_'+name+'.'+ext,0 if srun[ext+'_sha256']==sha(RESULTS/(name+'.'+ext)) else 1)
            check('summary_final_power_'+name,np.max(abs(np.asarray(srun['pk_final'])-arrays['pk'][-1])),1e-14)
        for name,digest in summary['code_sha256'].items():
            check('summary_production_code_hash_'+name,0 if digest==sha(ROOT/'code'/name) else 1)
        for filename,digest in summary['inputs'].items():
            candidates=[ROOT/filename,ROOT/'inputs'/filename,ROOT/'results'/filename]
            matches=[p for p in candidates if p.exists()]
            check('summary_input_hash_'+filename,0 if any(sha(p)==digest for p in matches) else 1)
        check('summary_protocol_hash',0 if summary['protocol_sha256']==sha(ROOT/'protocol.md') else 1)
        check('summary_initial_modes_hash',0 if summary['initial_modes_sha256']==sha(RESULTS/'initial_modes.npz') else 1)
        check('summary_momentum_factor',relative(summary['common_p_initial_coefficient'],factor),1e-13)
        ratios=[]
        for name in ('eps0','eps1e-4'):
            table=np.genfromtxt(RESULTS/f'background_{name}.csv',delimiter=',',names=True)
            ratios.append(np.min(299792.458*(2*np.pi/BOX)/(100*table['a']*table['E'])))
        check('summary_minimum_horizon_ratio',relative(min(ratios),summary['minimum_kfund_over_aHc']),1e-13)
        out['input_sha256'][str(summary_path.relative_to(ROOT))]=sha(summary_path)
    for path in (ROOT/'protocol.md',RESULTS/'initial_modes.npz',RESULTS/'background_summary.json',
                 ROOT/'inputs/reference_linear_spectrum.csv'):
        out['input_sha256'][str(path.relative_to(ROOT))]=sha(path)
    out['complete']=complete
    out['completed_runs']=len(found)
    out['expected_runs']=len(EXPECTED)
    out['missing_runs']=[name for name in EXPECTED if name not in found]
    out['summary_available']=summary_path.exists()
    out['checks_passed']=sum(row['passed'] for row in out['checks'])
    out['checks_total']=len(out['checks'])
    out['failures']=[row for row in out['checks'] if not row['passed']]
    out['status']='FAIL' if out['failures'] else ('PASS' if complete else 'PARTIAL_PASS')
    out['limitations']=['Saved densities/projections are float32; saved powers were calculated before quantization.',
                        'Particles are not stored at each output. Force/momentum summaries can be checked '
                        'against saved histories, not rederived from unsaved final particle arrays.',
                        'Common initial positions and momenta are reconstructible; no per-run entrance '
                        'array hashes were originally stored. All available initial density snapshots '
                        'and momentum statistics are independently checked.',
                        'Timing records describe concurrent execution, not a serial performance benchmark.']
    REPORT.parent.mkdir(exist_ok=True)
    if REPORT.exists():
        previous=REPORT.read_bytes()
        digest=hashlib.sha256(previous).hexdigest()
        (REPORT.parent/f'clustering_independent_audit_previous_{digest[:12]}.json').write_bytes(previous)
    REPORT.write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    print(json.dumps({key:out[key] for key in ('status','completed_runs','checks_passed','checks_total','failures')},indent=2))


if __name__=='__main__':
    run()
