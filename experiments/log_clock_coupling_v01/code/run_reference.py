#!/usr/bin/env python3
"""Independent DOP853 references for an autonomous positive-energy R + phi^4 model.

The inverse-log mass response is an input assumption. R is a mechanical clock
coordinate, not the canonical exponential-potential chi model from v0.2.
No observation fitting. No imports from the main lattice implementation.
"""
import os
for _key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[_key]='1'
from pathlib import Path
from time import perf_counter
import csv
import hashlib
import json
import platform
import numpy as np
import scipy
from scipy.integrate import solve_ivp

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results'
L=2*np.pi
RSTAR=np.exp(-2.)
R0=1.
U0=.1
B=2.
A=.5
TIMES=np.linspace(0.,20.,1001)
CMODES=np.array([1,3,5])


def mass(R,b=B):
    chi=np.log(R/RSTAR)
    return 1+b*(chi**-2-.25),chi


def artifact(name,**arrays):
    path=OUT/(name+'.npz')
    np.savez_compressed(path,**arrays)
    return dict(path=str(path.relative_to(ROOT)),sha256=hashlib.sha256(path.read_bytes()).hexdigest())


def energy_checks(field,clock,work,mode):
    baseline=max(abs(float(field[0])),1.)
    out=dict(max_relative_field_work_residual=float(np.max(np.abs(field-field[0]-work))/baseline))
    if mode=='full':
        total=field+clock
        out.update(max_relative_total_energy_drift=float(np.max(np.abs(total-total[0]))/max(abs(total[0]),1.)),
                   max_relative_clock_work_residual=float(np.max(np.abs(clock-clock[0]+work))/max(abs(clock[0]),1.)))
    return out


def lattice(N,I=100.,mode='full',b=B,amplitude=A,name=None,rtol=2e-12,atol=2e-14):
    a=L/N
    x=a*np.arange(N)
    qi=amplitude*np.cos(x)
    y0=np.r_[qi,np.zeros(N)]
    if mode=='full':
        y0=np.r_[y0,R0,U0,0.]
    elif mode=='external':
        y0=np.r_[y0,0.]
    elif mode!='frozen':
        raise ValueError(mode)
    def rhs(t,y):
        q,p=y[:N],y[N:2*N]
        q2=a*np.dot(q,q)
        if mode=='full':
            R,u,_W=y[2*N:]
            m2,chi=mass(R,b)
        elif mode=='external':
            R,u=R0+U0*t,U0
            m2,chi=mass(R,b)
        else:
            m2=1.
        lap=(np.roll(q,1)-2*q+np.roll(q,-1))/a**2
        main=np.r_[p,lap-m2*q-q**3]
        if mode=='frozen':
            return main
        power=-b*u*q2/(R*chi**3)
        if mode=='external':
            return np.r_[main,power]
        accel=b*q2/(I*R*chi**3)
        return np.r_[main,u,accel,power]
    started=perf_counter()
    sol=solve_ivp(rhs,(0.,20.),y0,method='DOP853',rtol=rtol,atol=atol,t_eval=TIMES)
    seconds=perf_counter()-started
    if not sol.success:
        raise RuntimeError(sol.message)
    q,p=sol.y[:N].T,sol.y[N:2*N].T
    arrays=dict(t=TIMES,x=x,q=q,velocity=p,spacing=np.array(a))
    if mode=='full':
        R,u,W=sol.y[2*N:]
        m2,chi=mass(R,b)
        clock=.5*I*u*u
        arrays.update(R=R,Rdot=u,chi=chi,mass_squared=m2,work=W,clock_energy=clock)
    elif mode=='external':
        R=R0+U0*TIMES
        u=np.full_like(TIMES,U0)
        m2,chi=mass(R,b)
        W=sol.y[2*N]
        clock=None
        arrays.update(R=R,Rdot=u,chi=chi,mass_squared=m2,work=W)
    else:
        m2=np.ones_like(TIMES)
        W=np.zeros_like(TIMES)
        clock=None
        arrays.update(mass_squared=m2,work=W)
    grad=(np.roll(q,-1,axis=1)-q)/a
    field=a*np.sum(.5*p*p+.5*grad*grad+.5*m2[:,None]*q*q+.25*q**4,axis=1)
    allq=np.fft.fft(q,axis=1,norm='forward')
    allp=np.fft.fft(p,axis=1,norm='forward')
    cos=2*allq[:,CMODES].real
    arrays.update(qhat=allq,velocity_hat=allp,k=np.fft.fftfreq(N,d=1/N).astype(int),
                  cos_modes=cos,cos_mode_numbers=CMODES,field_energy=field,
                  energy=field+clock if mode=='full' else field)
    if mode=='full':
        arrays['total_energy']=field+clock
    if name is None:
        name=f'reference_full_N{N}_I{I:g}' if mode=='full' else f'reference_{mode}_N{N}'
    meta=dict(name=name,mode=mode,N=N,I=I if mode=='full' else None,b=b,amplitude=amplitude,
              method='DOP853',rtol=rtol,atol=atol,nfev=sol.nfev,elapsed_seconds=seconds,
              all_finite=bool(np.all(np.isfinite(sol.y))),mass_squared_min=float(np.min(m2)),
              initial_field_energy=float(field[0]),final_field_energy=float(field[-1]),
              final_work=float(W[-1]),**energy_checks(field,clock,W,mode),
              artifact=artifact(name,**arrays))
    if mode!='frozen':
        meta.update(R_min=float(np.min(R)),chi_min=float(np.min(chi)),
                    R_final=float(R[-1]),Rdot_final=float(u[-1]),
                    max_R_departure_from_free=float(np.max(np.abs(R-(R0+U0*TIMES)))))
    return dict(meta=meta,cos=cos,field=field,clock=clock,work=W)


def project_cube(q,k,M):
    buf=np.zeros(M,complex)
    buf[k%M]=q
    f=np.fft.ifft(buf,norm='forward').real
    return np.fft.fft(f**3,norm='forward')[k%M]


def galerkin(K,rtol=2e-12,atol=2e-14,name=None):
    I=100.
    k=np.arange(-K,K+1)
    d=len(k)
    M=4*K+1
    q0=np.zeros(d,complex)
    q0[np.abs(k)==1]=A/2
    y0=np.r_[q0,np.zeros(d,complex),complex(R0),complex(U0),0j]
    def rhs(_t,y):
        q,p=y[:d],y[d:2*d]
        R,u,_W=y[2*d:].real
        m2,chi=mass(R)
        Q2=L*np.sum(np.abs(q)**2)
        acc=-((k*k)+m2)*q-project_cube(q,k,M)
        racc=B*Q2/(I*R*chi**3)
        power=-B*u*Q2/(R*chi**3)
        return np.r_[p,acc,u,racc,power]
    started=perf_counter()
    sol=solve_ivp(rhs,(0.,20.),y0,method='DOP853',rtol=rtol,atol=atol,t_eval=TIMES)
    seconds=perf_counter()-started
    if not sol.success:
        raise RuntimeError(sol.message)
    q,p=sol.y[:d].T,sol.y[d:2*d].T
    R,u,W=sol.y[2*d:].real
    m2,chi=mass(R)
    buf=np.zeros((len(TIMES),M),complex)
    buf[:,k%M]=q
    fieldgrid=np.fft.ifft(buf,axis=1,norm='forward').real
    field=L*(.5*np.sum(np.abs(p)**2+(k*k+m2[:,None])*np.abs(q)**2,axis=1)
             +.25*np.mean(fieldgrid**4,axis=1))
    clock=.5*I*u*u
    cos=np.stack([2*q[:,np.where(k==m)[0][0]].real for m in CMODES],axis=1)
    if name is None:
        name=f'reference_continuum_I100_K{K}'
    meta=dict(name=name,mode='full_continuum',I=I,b=B,amplitude=A,K=K,M=M,
              method='DOP853',rtol=rtol,atol=atol,nfev=sol.nfev,elapsed_seconds=seconds,
              all_finite=bool(np.all(np.isfinite(sol.y))),
              mass_squared_min=float(np.min(m2)),R_min=float(np.min(R)),chi_min=float(np.min(chi)),
              R_final=float(R[-1]),Rdot_final=float(u[-1]),
              hermitian_error=float(max(np.max(np.abs(q-np.conj(q[:,::-1]))),np.max(np.abs(p-np.conj(p[:,::-1]))))),
              real_clock_imaginary_error=float(np.max(np.abs(sol.y[2*d:].imag))),
              initial_field_energy=float(field[0]),final_field_energy=float(field[-1]),
              final_work=float(W[-1]),**energy_checks(field,clock,W,'full'),
              artifact=artifact(name,t=TIMES,k=k,qhat=q,velocity_hat=p,R=R,Rdot=u,chi=chi,
                                mass_squared=m2,work=W,field_energy=field,clock_energy=clock,
                                total_energy=field+clock,energy=field+clock,cos_modes=cos,
                                cos_mode_numbers=CMODES,collocation_M=np.array(M)))
    return dict(meta=meta,k=k,q=q,p=p,R=R,Rdot=u,cos=cos,field=field,clock=clock,work=W)


def compare(a,b):
    off=int(b['k'][-1])
    aq=np.zeros_like(b['q']); ap=np.zeros_like(b['p'])
    aq[:,a['k']+off]=a['q']; ap[:,a['k']+off]=a['p']
    return dict(max_field_rms_difference=float(np.max(np.sqrt(np.sum(np.abs(aq-b['q'])**2,axis=1)))),
                max_velocity_rms_difference=float(np.max(np.sqrt(np.sum(np.abs(ap-b['p'])**2,axis=1)))),
                max_R_difference=float(np.max(np.abs(a['R']-b['R']))),
                max_Rdot_difference=float(np.max(np.abs(a['Rdot']-b['Rdot']))))


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    runs=[lattice(128,I=I) for I in (10.,100.,1000.)]
    runs.extend([lattice(N,I=100.) for N in (64,256)])
    runs.extend([lattice(128,mode=mode) for mode in ('external','frozen')])
    null_b=lattice(32,I=100.,b=0.,name='reference_null_b0_N32')
    null_q=lattice(16,I=100.,amplitude=0.,name='reference_null_q0_N16')
    runs.extend([null_b,null_q])
    c16=galerkin(16)
    c24=galerkin(24)
    ctight=galerkin(24,rtol=2e-13,atol=2e-15,name='reference_continuum_I100_K24_tight')
    runs.extend([c16,c24,ctight])
    checks=[]
    for run in runs:
        m=run['meta']
        item=dict(name=m['name']+'_physical_domain_and_work',
                  passed=bool(m['all_finite'] and m['mass_squared_min']>0 and
                              m.get('chi_min',1.)>0 and m['max_relative_field_work_residual']<1e-9),
                  max_field_work_residual=m['max_relative_field_work_residual'])
        checks.append(item)
        if 'max_relative_total_energy_drift' in m:
            checks.append(dict(name=m['name']+'_clock_total_energy',
                               passed=bool(m['max_relative_total_energy_drift']<1e-9 and m['max_relative_clock_work_residual']<1e-9),
                               total_energy_drift=m['max_relative_total_energy_drift'],
                               clock_work_residual=m['max_relative_clock_work_residual']))
    for case in (null_b,null_q):
        m=case['meta']
        checks.append(dict(name=m['name']+'_free_clock_analytic',
                           passed=bool(m['max_R_departure_from_free']<1e-12 and abs(m['final_work'])<1e-14),
                           max_R_error=m['max_R_departure_from_free'],work=m['final_work']))
    convergence=dict(spatial_K16_to_K24=compare(c16,c24),tolerance_K24=compare(c24,ctight))
    for label,comp in convergence.items():
        checks.append(dict(name=label,passed=bool(max(comp.values())<1e-9),**comp))
    field_analytic=67*np.pi/256
    checks.append(dict(name='continuum_initial_energy_analytic',
                       passed=bool(abs(ctight['field'][0]-field_analytic)<1e-13),
                       absolute_error=float(abs(ctight['field'][0]-field_analytic))))
    result=dict(scope='Autonomous mechanical R-clock + phi^4 toy model. Inverse-log mass response assumed, not derived; distinct from canonical exponential-potential chi completion. No observed data.',
                parameters=dict(L=L,v=1.,m0=1.,g=1.,A=A,R_initial=R0,Rdot_initial=U0,Rstar=RSTAR,
                                chi_initial=2.,Td=10.,b=B,t_end=20.,output_dt=.02),
                formulas=dict(mass_squared='1+b*(log(R/Rstar)**(-2)-0.25)',
                              R_acceleration='b*integral(q**2)/(I*R*chi**3)',
                              field_power='-b*Rdot*integral(q**2)/(R*chi**3)',
                              free_clock='R=1+t/Td; chi=2+log(1+t/Td)'),
                Fourier_normalization='phi=sum qhat[k]*exp(i*k*x); fft norm=forward. Continuum k sorted -K..K; lattice k in FFT order. Stored velocity is qdot, not canonical a*qdot.',
                energy_normalization='energy=total_energy for full coupling; energy=field_energy for external/frozen controls. work integrates field power independently.',
                reference_recommendation='results/reference_continuum_I100_K24_tight.npz',
                code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                software=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__),
                threads={x:os.environ[x] for x in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS')},
                runs=[r['meta'] for r in runs],convergence=convergence,
                checks=checks,all_checks_passed=all(x['passed'] for x in checks))
    (OUT/'reference_summary.json').write_text(json.dumps(result,indent=2)+'\n')
    with (OUT/'reference_modes.csv').open('w',newline='') as stream:
        writer=csv.writer(stream)
        writer.writerow(['reference','t','cos1','cos3','cos5','field_energy','clock_energy','work'])
        for r in runs:
            for i,t in enumerate(TIMES):
                ce='' if r['clock'] is None else r['clock'][i]
                writer.writerow([r['meta']['name'],t,*r['cos'][i],r['field'][i],ce,r['work'][i]])
    print(json.dumps(dict(all_checks_passed=result['all_checks_passed'],checks=len(checks),
                         convergence=convergence,runs=[r['meta'] for r in runs]),indent=2))
    if not result['all_checks_passed']:
        raise SystemExit('Reference validation failed; inspect reference_summary.json.')

if __name__=='__main__':
    main()
