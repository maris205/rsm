#!/usr/bin/env python3
"""Targeted precision audit: rerun only N512 lattice with tighter tolerance.

Does not modify/import the earlier reference script or rerun other models.
"""
import os
for _name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[_name]='1'
from pathlib import Path
from time import perf_counter
import hashlib
import json
import platform
import numpy as np
import scipy
from scipy.integrate import solve_ivp

ROOT=Path(__file__).resolve().parents[1]
R=ROOT/'results'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rms_max(array):
    return float(np.max(np.sqrt(np.mean(array**2,axis=1))))


def relative_all_time_error(field,reference):
    return float(np.sqrt(np.sum((field-reference)**2)/np.sum(reference**2)))


def main():
    original_path=R/'reference_lattice_N512.npz'
    correction_path=R/'reference_corrected_N512.npz'
    baseline=np.load(original_path)
    correction=np.load(correction_path)
    N=512
    a=2*np.pi/N
    x=a*np.arange(N)
    t=baseline['t']
    initial=np.r_[.5*np.cos(x),np.zeros(N)]
    def rhs(_t,y):
        q,p=y[:N],y[N:]
        lap=(np.roll(q,1)+np.roll(q,-1)-2*q)/a**2
        return np.r_[p,lap-q-q**3]
    began=perf_counter()
    sol=solve_ivp(rhs,(0.,20.),initial,method='DOP853',
                  rtol=2e-13,atol=2e-15,t_eval=t)
    elapsed=perf_counter()-began
    if not sol.success:
        raise RuntimeError(sol.message)
    q,p=sol.y[:N].T.copy(),sol.y[N:].T.copy()
    energy=a*np.sum(.5*p*p+.5*((np.roll(q,-1,axis=1)-q)/a)**2+
                    .5*q*q+.25*q**4,axis=1)
    qhat=np.fft.fft(q,axis=1,norm='forward')
    vhat=np.fft.fft(p,axis=1,norm='forward')
    target=R/'reference_lattice_N512_tight.npz'
    np.savez_compressed(target,t=t,x=x,q=q,velocity=p,energy=energy,
                        qhat=qhat,velocity_hat=vhat,
                        k=np.fft.fftfreq(N,d=1/N).astype(int),spacing=np.array(a),
                        cos_modes=2*qhat[:,[1,3,5]].real,
                        cos_mode_numbers=np.array([1,3,5]))
    buf=np.zeros((len(t),N),complex)
    buf[:,correction['k']%N]=correction['qhat']
    qc=np.fft.ifft(buf,axis=1,norm='forward').real
    err_old=rms_max(baseline['q']-qc)
    err_tight=rms_max(q-qc)
    global_old=relative_all_time_error(baseline['q'],qc)
    global_tight=relative_all_time_error(q,qc)
    state=dict(max_field_rms_difference=rms_max(q-baseline['q']),
               max_velocity_rms_difference=rms_max(p-baseline['velocity']),
               all_time_relative_field_difference=relative_all_time_error(q,baseline['q']))
    effect=dict(original_max_rms_corrected_error=err_old,
                tight_max_rms_corrected_error=err_tight,
                fractional_change_in_max_rms_corrected_error=err_tight/err_old-1,
                original_all_time_relative_corrected_error=global_old,
                tight_all_time_relative_corrected_error=global_tight,
                fractional_change_in_all_time_relative_corrected_error=global_tight/global_old-1,
                max_field_rms_reference_change_over_tight_correction_error=state['max_field_rms_difference']/err_tight)
    drift=float(np.max(np.abs(energy/energy[0]-1)))
    checks=[dict(name='finite_states_and_energy',passed=bool(np.all(np.isfinite(q)) and np.all(np.isfinite(p)) and np.all(np.isfinite(energy)))),
            dict(name='tight_energy_drift',passed=bool(drift<1e-10),value=drift,limit=1e-10),
            dict(name='field_reference_change_under_1percent_of_correction_error',
                 passed=bool(state['max_field_rms_difference']<.01*err_tight),
                 value=state['max_field_rms_difference']/err_tight,limit=.01),
            dict(name='max_rms_correction_error_change_under_1percent',
                 passed=bool(abs(effect['fractional_change_in_max_rms_corrected_error'])<.01),
                 value=effect['fractional_change_in_max_rms_corrected_error'],limit=.01),
            dict(name='protocol_all_time_error_change_under_1percent',
                 passed=bool(abs(effect['fractional_change_in_all_time_relative_corrected_error'])<.01),
                 value=effect['fractional_change_in_all_time_relative_corrected_error'],limit=.01)]
    result=dict(scope='One targeted N512 lattice tolerance refinement, same fixed physical model and initial data; existing corrected PDE trajectory reused without rerun.',
                method='DOP853',N=N,rtol=2e-13,atol=2e-15,nfev=sol.nfev,
                elapsed_seconds=elapsed,max_relative_energy_drift=drift,
                state_differences=state,correction_error_budget=effect,
                input_hashes={original_path.name:sha(original_path),correction_path.name:sha(correction_path)},
                code_sha256=sha(Path(__file__)),artifact=dict(path=str(target.relative_to(ROOT)),sha256=sha(target)),
                software=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__),
                checks=checks,all_checks_passed=all(x['passed'] for x in checks))
    (R/'reference_precision_audit.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
    if not result['all_checks_passed']:
        raise SystemExit('Targeted precision audit failed.')

if __name__=='__main__':
    main()
