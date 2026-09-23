#!/usr/bin/env python3
"""Independent Galerkin and high-accuracy lattice reference trajectories.

No imports from the production lattice implementation. Normalization:
phi(x)=sum_k qhat[k] exp(i k x), Fourier forward FFT normalized by 1/M.
The +a^2 phi_xxxx/12 correction is only a cutoff, long-wave effective PDE.
"""
import os
for _name in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[_name] = '1'
from pathlib import Path
import csv
import hashlib
import json
import platform
from time import perf_counter
import numpy as np
import scipy
from scipy.integrate import solve_ivp

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / 'results'
L = 2*np.pi
A = .5
T = np.linspace(0., 20., 1001)
NS = (64, 128, 256, 512)
COS_MODES = np.array([1, 3, 5])


def save_npz(name, **arrays):
    path = RESULTS / name
    np.savez_compressed(path, **arrays)
    return dict(path=str(path.relative_to(ROOT)),
                sha256=hashlib.sha256(path.read_bytes()).hexdigest())


def project_cube(coeff, k, size):
    buf = np.zeros(size, dtype=complex)
    buf[k % size] = coeff
    field = np.fft.ifft(buf, norm='forward')
    return np.fft.fft(field.real**3, norm='forward')[k % size]


def solve_galerkin(K, rtol=1e-11, atol=1e-13, spacing=0., suffix=None):
    k = np.arange(-K, K+1)
    M = 4*K+1
    omega2 = 1+k*k-spacing**2*k**4/12
    if np.min(omega2) <= 0:
        raise ValueError('The retained correction band is not linearly stable.')
    qinit = np.zeros(2*K+1, dtype=complex)
    qinit[np.abs(k)==1] = A/2
    yinit = np.r_[qinit, np.zeros_like(qinit)]
    d = len(k)
    def rhs(_t, y):
        q, p = y[:d], y[d:]
        return np.r_[p, -omega2*q-project_cube(q, k, M)]
    began = perf_counter()
    sol = solve_ivp(rhs, (T[0], T[-1]), yinit, method='DOP853',
                    rtol=rtol, atol=atol, t_eval=T, dense_output=True)
    elapsed = perf_counter()-began
    if not sol.success:
        raise RuntimeError(sol.message)
    q, p = sol.y[:d].T.copy(), sol.y[d:].T.copy()
    buf = np.zeros((len(T), M), dtype=complex)
    buf[:, k % M] = q
    field = np.fft.ifft(buf, axis=1, norm='forward')
    energy = L*(.5*np.sum(np.abs(p)**2+omega2*np.abs(q)**2, axis=1)
                +.25*np.mean(field.real**4, axis=1))
    cos = np.stack([2*q[:, np.where(k==m)[0][0]].real for m in COS_MODES], axis=1)
    fname = suffix or f'reference_continuum_K{K}'
    artifact = save_npz(fname+'.npz', t=T, k=k, qhat=q, velocity_hat=p,
                        cos_modes=cos, cos_mode_numbers=COS_MODES,
                        energy=energy, omega_squared=omega2,
                        collocation_M=np.array(M), spacing=np.array(spacing))
    hermitian = max(np.max(np.abs(q-np.conj(q[:, ::-1]))),
                    np.max(np.abs(p-np.conj(p[:, ::-1]))))
    metadata = dict(name=fname, equation='continuum' if spacing==0 else 'cutoff_long_wave_correction',
                    K=K, M=M, spacing=spacing, max_abs_ka=float(spacing*K),
                    omega_squared_min=float(np.min(omega2)),
                    all_retained_modes_linearly_stable=bool(np.min(omega2)>0),
                    rtol=rtol, atol=atol, method='DOP853', nfev=sol.nfev,
                    elapsed_seconds=elapsed, hermitian_error=float(hermitian),
                    max_relative_energy_drift=float(np.max(np.abs(energy/energy[0]-1))),
                    all_finite=bool(np.all(np.isfinite(q)) and np.all(np.isfinite(p))),
                    initial_energy=float(energy[0]), artifact=artifact)
    if spacing:
        modal_energy=.5*(np.abs(p)**2+omega2*np.abs(q)**2)
        tail=np.abs(k)*spacing>.5
        metadata['max_quadratic_energy_fraction_above_abs_ka_0p5']=float(
            np.max(np.sum(modal_energy[:, tail], axis=1)/np.sum(modal_energy, axis=1)))
        metadata['validity_note']='Finite retained band only; positivity of omega^2 is not proof of short-wave accuracy. No untruncated +fourth-derivative PDE is solved.'
    # Analytic initial mode-3 coefficient: a3(t)=-t^2/64+O(t^4).
    accel = rhs(0., yinit)[d:]
    metadata['mode3_initial_t_squared_coefficient']=float(accel[np.where(k==3)[0][0]].real)
    return dict(metadata=metadata, k=k, q=q, p=p, cos=cos, energy=energy)


def solve_lattice(N, rtol=2e-12, atol=2e-14):
    a=L/N
    x=np.arange(N)*a
    initial=np.r_[A*np.cos(x), np.zeros(N)]
    def rhs(_t, y):
        q,p=y[:N],y[N:]
        lap=(np.roll(q,1)-2*q+np.roll(q,-1))/a**2
        return np.r_[p, lap-q-q**3]
    began=perf_counter()
    sol=solve_ivp(rhs, (T[0],T[-1]), initial, method='DOP853',
                  rtol=rtol, atol=atol, t_eval=T)
    elapsed=perf_counter()-began
    if not sol.success:
        raise RuntimeError(sol.message)
    q,p=sol.y[:N].T.copy(),sol.y[N:].T.copy()
    grad=(np.roll(q,-1,axis=1)-q)/a
    energy=a*np.sum(.5*p*p+.5*grad*grad+.5*q*q+.25*q**4,axis=1)
    allhat=np.fft.fft(q,axis=1,norm='forward')
    allvelhat=np.fft.fft(p,axis=1,norm='forward')
    cos=2*allhat[:,COS_MODES].real
    name=f'reference_lattice_N{N}'
    artifact=save_npz(name+'.npz',t=T,x=x,q=q,velocity=p,energy=energy,
                      qhat=allhat,velocity_hat=allvelhat,
                      k=np.fft.fftfreq(N,d=1/N).astype(int),
                      cos_modes=cos,cos_mode_numbers=COS_MODES,
                      spacing=np.array(a))
    expected=L*(A*A/4*(1+4*np.sin(a/2)**2/a**2)+3*A**4/32)
    metadata=dict(name=name,N=N,spacing=a,method='DOP853',rtol=rtol,atol=atol,
                  nfev=sol.nfev,elapsed_seconds=elapsed,
                  max_relative_energy_drift=float(np.max(np.abs(energy/energy[0]-1))),
                  initial_energy=float(energy[0]),initial_energy_analytic=float(expected),
                  initial_energy_relative_error=float(abs(energy[0]/expected-1)),
                  all_finite=bool(np.all(np.isfinite(q)) and np.all(np.isfinite(p))),
                  artifact=artifact)
    return dict(metadata=metadata,cos=cos,energy=energy)


def compare_spectra(left, right):
    kl,kr=left['k'],right['k']
    def pad(k, data):
        out=np.zeros_like(right['q'])
        out[:,k+int(kr[-1])]=data
        return out
    dq=pad(kl,left['q'])-right['q']
    dp=pad(kl,left['p'])-right['p']
    return dict(max_field_rms_difference=float(np.max(np.sqrt(np.sum(np.abs(dq)**2,axis=1)))),
                max_velocity_rms_difference=float(np.max(np.sqrt(np.sum(np.abs(dp)**2,axis=1)))),
                max_cos_mode_amplitude_difference=float(np.max(np.abs(left['cos']-right['cos']))))


def main():
    RESULTS.mkdir(parents=True,exist_ok=True)
    checks=[]
    def check(name, passed, **extra):
        checks.append(dict(name=name,passed=bool(passed),**extra))
    # Alias check against a larger quadrature for a nontrivial Hermitian state.
    rng=np.random.default_rng(2041)
    ktest=np.arange(-24,25)
    coeff=np.zeros(49,dtype=complex)
    coeff[24]=.11
    for k in range(1,25):
        coeff[24+k]=(rng.normal()+1j*rng.normal())/(1+k*k)
        coeff[24-k]=np.conj(coeff[24+k])
    projection_error=float(np.max(np.abs(project_cube(coeff,ktest,97)-project_cube(coeff,ktest,195))))
    check('cubic_projection_dealiasing',projection_error<1e-12,max_abs_difference=projection_error)
    base16=solve_galerkin(16)
    base24=solve_galerkin(24)
    tight=solve_galerkin(24,rtol=2e-13,atol=2e-15,suffix='reference_continuum_K24_tight')
    corrected=[solve_galerkin(24,rtol=2e-13,atol=2e-15,spacing=L/N,
                             suffix=f'reference_corrected_N{N}') for N in NS]
    lattice=[solve_lattice(N) for N in NS]
    galerkins=[base16,base24,tight]+corrected
    for run in galerkins:
        m=run['metadata']
        check(m['name']+'_finite_energy_symmetry', m['all_finite'] and
              m['max_relative_energy_drift']<1e-9 and m['hermitian_error']<1e-12,
              energy_drift=m['max_relative_energy_drift'],hermitian_error=m['hermitian_error'])
        check(m['name']+'_mode3_initial_acceleration',
              abs(m['mode3_initial_t_squared_coefficient']+1/64)<1e-13,
              value=m['mode3_initial_t_squared_coefficient'])
    for run in lattice:
        m=run['metadata']
        check(m['name']+'_finite_energy_initial',m['all_finite'] and
              m['max_relative_energy_drift']<1e-9 and m['initial_energy_relative_error']<1e-12,
              energy_drift=m['max_relative_energy_drift'],initial_energy_error=m['initial_energy_relative_error'])
    analytic=67*np.pi/256
    err=abs(tight['energy'][0]/analytic-1)
    check('continuum_initial_energy_67pi_over_256',err<1e-13,relative_error=float(err))
    spatial=compare_spectra(base16,base24)
    temporal=compare_spectra(base24,tight)
    for name,comp in [('Galerkin_K16_to_K24',spatial),('DOP853_tolerance_refinement',temporal)]:
        check(name,max(comp.values())<2e-8,**comp)
    allruns=galerkins+lattice
    result=dict(scope='Independent mathematical reference solutions; autonomous phi^4 example only; no observations, no variable constants.',
                parameters=dict(L=L,v=1.,omega0=1.,g=1.,A=A,T=20.,sample_dt=.02,samples=len(T)),
                Fourier_normalization='phi=sum qhat[k]*exp(i*k*x); fft(norm=forward); continuum k is sorted -K..K, lattice k is FFT order',
                corrected_equation='phi_tt=phi_xx-phi-phi^3+a^2 phi_xxxx/12, projected to |k|<=24; not an uncut PDE',
                recommended_continuum_reference='results/reference_continuum_K24_tight.npz',
                original_lattice_reference='results/reference_lattice_N128.npz (also N64, N256, N512)',
                software=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__),
                thread_environment={x:os.environ[x] for x in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS')},
                code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                runs=[r['metadata'] for r in allruns],
                convergence=dict(spatial_K16_to_K24=spatial,time_tolerance_K24=temporal),
                checks=checks,all_checks_passed=all(c['passed'] for c in checks))
    (RESULTS/'reference_summary.json').write_text(json.dumps(result,indent=2)+'\n')
    with (RESULTS/'reference_modes.csv').open('w',newline='') as stream:
        writer=csv.writer(stream)
        writer.writerow(['reference','t','cos1','cos3','cos5','energy'])
        for run in allruns:
            for i,t in enumerate(T):
                writer.writerow([run['metadata']['name'],t,*run['cos'][i],run['energy'][i]])
    print(json.dumps(dict(all_checks_passed=result['all_checks_passed'],checks=len(checks),
                         convergence=result['convergence'],
                         runs=[dict(name=r['metadata']['name'],nfev=r['metadata']['nfev'],
                                    seconds=r['metadata']['elapsed_seconds'],energy_drift=r['metadata']['max_relative_energy_drift']) for r in allruns]),indent=2))
    if not result['all_checks_passed']:
        raise SystemExit('A reference validation failed; inspect reference_summary.json.')

if __name__=='__main__':
    main()
