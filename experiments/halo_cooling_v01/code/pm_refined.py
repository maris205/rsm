"""Refined periodic PM with a controlled compensated sphere.

Only the FFT storage differs from the archived CIC force operator. Positions
are in box units, tau=H0*t, and p=a^2 dx/dtau. Adaptive steps are not claimed
to define a strictly symplectic or time-reversible integration map.
"""
from __future__ import annotations

import hashlib
import importlib.util
from pathlib import Path

import numpy as np
from scipy import fft

ROOT = Path(__file__).resolve().parents[1]
LEGACY_PATH = ROOT.parent / 'cosmic_bridge_v01/code/pm.py'
LEGACY_SHA256 = '1b2607ac0b7ba35c5df9ceb7375910ba75b7bed65a1edec79c0583f067a5d982'
if hashlib.sha256(LEGACY_PATH.read_bytes()).hexdigest() != LEGACY_SHA256:
    raise RuntimeError('Archived PM source changed; refusing unrecorded algorithm change')
_spec = importlib.util.spec_from_file_location('archived_cosmic_bridge_pm', LEGACY_PATH)
legacy = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(legacy)
cic_deposit, cic_gather = legacy.cic_deposit, legacy.cic_gather
Background, advance_kdk = legacy.Background, legacy.advance_kdk


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def array_sha(value):
    return hashlib.sha256(np.ascontiguousarray(value).view(np.uint8)).hexdigest()


class PMGrid:
    """Same spectral kernel and CIC as legacy; real FFT, component-wise gather."""

    def __init__(self, nmesh, workers=4):
        self.nmesh, self.workers = int(nmesh), int(workers)
        if self.nmesh < 2 or self.workers < 1 or self.workers > 4:
            raise ValueError('Need mesh>=2 and 1<=FFT workers<=4')
        k = 2*np.pi*fft.fftfreq(self.nmesh, d=1/self.nmesh)
        kz = 2*np.pi*fft.rfftfreq(self.nmesh, d=1/self.nmesh)
        k2 = k[:,None,None]**2 + k[None,:,None]**2 + kz[None,None,:]**2
        self.inverse_k2 = np.zeros_like(k2)
        np.divide(1, k2, out=self.inverse_k2, where=k2>0)
        gx, gz = k.copy(), kz.copy()
        if self.nmesh % 2 == 0:
            gx[self.nmesh//2] = 0
            gz[-1] = 0
        self.gradient_k = (gx[:,None,None],gx[None,:,None],gz[None,None,:])

    def _components(self, delta):
        if delta.shape != (self.nmesh,)*3 or not np.all(np.isfinite(delta)):
            raise ValueError('Finite cubic density grid matching force mesh required')
        spectrum = fft.rfftn(delta, workers=self.workers)
        for wave in self.gradient_k:
            yield fft.irfftn(1j*wave*self.inverse_k2*spectrum,
                            s=(self.nmesh,)*3, workers=self.workers)

    def mesh_force(self, delta):
        """Validation API. irfft has real output; no measured imaginary residual."""
        force = np.empty(delta.shape+(3,), dtype=float)
        for axis, component in enumerate(self._components(delta)):
            force[...,axis] = component
        return force, 0.0

    def force(self, positions):
        mass = cic_deposit(positions, self.nmesh)
        mass_sum = float(np.sum(mass))
        mass *= self.nmesh**3
        mass -= 1
        force = np.empty_like(positions, dtype=float)
        for axis, component in enumerate(self._components(mass)):
            force[:,axis] = cic_gather(component,positions)
        if not np.all(np.isfinite(force)):
            raise FloatingPointError('Nonfinite refined PM force')
        mean = np.mean(force,axis=0)
        norm = np.linalg.norm(force,axis=1)
        absolute = float(np.linalg.norm(mean))
        typical = float(np.mean(norm))
        return force, {
            'raw_mass':mass_sum, 'mass_relative_error':abs(mass_sum-1),
            'density_contrast_mean':float(np.mean(mass)),
            'total_force_absolute':absolute,
            'total_force_normalized':absolute/max(typical,1e-300),
            'mean_force_magnitude':typical,'max_force_magnitude':float(np.max(norm)),
            'real_fft':True, 'imaginary_residual_status':'not observable from irfft; separately checked against full FFT',
        }


def minimum_vector(positions, center):
    return np.remainder(positions-np.asarray(center)+0.5,1)-0.5


def compensated_initial(nparticle, delta_i, reference_background,
                        shift_cells=(0,0,0), force_mesh=256):
    """Finite-ai nonlinear enclosed-density mapping; root protocol equation."""
    if int(nparticle)!=nparticle or nparticle<2 or not np.isfinite(delta_i) or delta_i<=0:
        raise ValueError('Positive finite delta_i and integer particle grid required')
    nparticle = int(nparticle)
    ai, ql, qc, box = .02,.25,.4,10.
    bg = reference_background.evaluate_a(ai)
    e_ref, f_ref = float(bg['E']), float(reference_background.f_reference)
    grid=(np.arange(nparticle)+.5)/nparticle
    lagrangian=np.stack(np.meshgrid(grid,grid,grid,indexing='ij'),axis=-1).reshape((-1,3))
    q=minimum_vector(lagrangian,np.full(3,.5))
    qr=np.linalg.norm(q,axis=1)
    s=np.clip((qr-ql)/(qc-ql),0,1)
    window=1-10*s**3+15*s**4-6*s**5
    window[qr>=qc]=0
    enclosed=delta_i*window
    scale=(1+enclosed)**(-1/3)
    yn=-f_ref*enclosed*scale/(3*(1+enclosed))
    shift=np.asarray(shift_cells,dtype=float)/force_mesh
    if shift.shape!=(3,) or not np.all(np.isfinite(shift)):
        raise ValueError('Finite three-component grid shift required')
    center=np.remainder(np.full(3,.5)+shift,1)
    positions=np.remainder(center+q*scale[:,None],1)
    momenta=ai**2*e_ref*q*yn[:,None]
    window_derivative=(-30*s**2+60*s**3-30*s**4)/(qc-ql)
    scale_derivative=-delta_i*window_derivative/(3*(1+enclosed)**(4/3))
    jacobian=scale**2*(scale+qr*scale_derivative)
    if np.min(jacobian)<=0 or not np.all(np.isfinite(momenta)):
        raise FloatingPointError('Initial mapping is noninvertible or nonfinite')
    core=(qr>=.25*ql)&(qr<=.75*ql)
    rho0=2.77536627e11*.315
    metadata={
        'nparticle':nparticle,'particles':nparticle**3,'box_mpc_h':box,
        'a_initial':ai,'delta_nl_core':float(delta_i),'scale_core':float((1+delta_i)**(-1/3)),
        'core_lagrangian_radius_box':ql,'compensation_radius_box':qc,
        'core_label_interval_q_over_qL':[.25,.75], 'core_label_count':int(core.sum()),
        'E_reference_initial':e_ref,'f_reference_initial':f_ref,
        'velocity_core_y_N':float(-f_ref*delta_i*(1+delta_i)**(-4/3)/3),
        'shift_box':shift.tolist(),'center_box':center.tolist(),
        'mass_excess_compensated':True,'mass_compensation_scope':'continuum enclosed map; particle discretization measured separately',
        'minimum_mapping_jacobian':float(np.min(jacobian)),
        'particle_mass_msun_h':float(rho0*box**3/nparticle**3),
        'continuum_core_mass_msun_h':float(4*np.pi/3*rho0*(box*ql)**3),
        'discrete_core_particle_count':int(np.sum(qr<=ql)),
        'discrete_core_mass_msun_h':float(np.sum(qr<=ql)*rho0*box**3/nparticle**3),
        'positions_sha256':array_sha(positions),'momenta_sha256':array_sha(momenta),
        'interpretation':'Prescribed compensated sphere; not observed or random-spectrum initial conditions',
    }
    return {'positions':positions,'momenta':momenta,'q':q,'q_radius':qr,
            'center':center,'core_mask':core,'metadata':metadata}


def core_diagnostics(positions,momenta,a,background,initial):
    mask=initial['core_mask']
    vectors=minimum_vector(positions[mask],initial['center'])
    radius=np.linalg.norm(vectors,axis=1)
    q=initial['q_radius'][mask]
    ratios=radius/q
    q16,median,q84=np.percentile(ratios,[16,50,84])
    if not np.isfinite(median) or median<=0:
        raise FloatingPointError('Nonpositive/nonfinite core scale')
    e=float(background.E(a))
    radial_p=np.sum(momenta[mask]*vectors,axis=1)/np.maximum(radius,1e-300)
    yn=radial_p/(a*a*e*q)
    transverse=momenta[mask]-radial_p[:,None]*vectors/np.maximum(radius[:,None],1e-300)
    eig=np.linalg.eigvalsh(vectors.T@vectors/len(vectors))
    axis_ratio=float(np.sqrt(max(eig[0],0)/max(eig[-1],1e-300)))
    proxy_radius=median*initial['metadata']['core_lagrangian_radius_box']
    all_radius=np.linalg.norm(minimum_vector(positions,initial['center']),axis=1)
    enclosed_count=int(np.sum(all_radius<=proxy_radius))
    enclosed_ratio=enclosed_count/(len(positions)*4*np.pi/3*proxy_radius**3)
    return {
        'a':float(a),'z':float(1/a-1),'y_median':float(median),
        'y_p16':float(q16),'y_p84':float(q84),'Delta_proxy':float(median**-3),
        'nonhomology_fraction':float((q84-q16)/(2*median)),
        'axis_ratio_min_max':axis_ratio,'axis_deviation':1-axis_ratio,
        'y_N_median':float(np.median(yn)),
        'physical_expansion_ratio_median':float(np.median(ratios+yn)),
        'transverse_momentum_rms':float(np.sqrt(np.mean(np.sum(transverse**2,axis=1)))),
        'radial_momentum_rms':float(np.sqrt(np.mean(radial_p**2))),
        'radius_proxy_comoving_box':float(proxy_radius),
        'actual_enclosed_particle_count_at_proxy_radius':enclosed_count,
        'actual_enclosed_mean_density_ratio':float(enclosed_ratio),
    }


def core_scale(positions,initial):
    mask=initial['core_mask']
    radius=np.linalg.norm(minimum_vector(positions[mask],initial['center']),axis=1)
    return float(np.median(radius/initial['q_radius'][mask]))


def time_factor(background,a1,a2):
    """Actual Delta(H0 t), independent of kick/drift integral powers."""
    left,right=np.log(a1),np.log(a2)
    nodes=(right+left)/2+(right-left)/2*background._gauss_nodes
    return float((right-left)/2*np.sum(background._gauss_weights/background.E(np.exp(nodes))))


def dynamical_time(a,y,omega_m0):
    return float(np.sqrt(2/omega_m0)*a**1.5*y**1.5)


def predict_positions(x,p,g,a1,a2,background):
    midpoint=np.sqrt(a1*a2)
    k1,_=background.factors(a1,midpoint)
    _,drift=background.factors(a1,a2)
    displacement=drift*(p+1.5*background.omega_m0*k1*g)
    return np.remainder(x+displacement,1),displacement


def select_step(x,p,g,a,background,initial,nmesh,a_max=.55,halfstep=False):
    factor=.5 if halfstep else 1.
    d_n=min(np.log(a_max/a),np.log(.55/.02)/256*factor)
    y=core_scale(x,initial)
    tdyn=dynamical_time(a,y,background.omega_m0)
    for reductions in range(48):
        end=min(float(a*np.exp(d_n)),a_max)
        _,displacement=predict_positions(x,p,g,a,end,background)
        cells=float(np.max(np.linalg.norm(displacement,axis=1))*nmesh)
        dt=time_factor(background,a,end)
        if cells<=.15*factor and dt/tdyn<=.03*factor:
            return end,{'predicted_max_drift_cells':cells,'predicted_dt_over_tdyn':dt/tdyn,
                        'step_halvings':reductions,'dt_H0':dt,'delta_ln_a':float(np.log(end/a))}
        d_n*=.5
    raise RuntimeError('Adaptive step failed to find positive bounded step')


def final_profiles(positions,momenta,a,background,initial,analysis_mesh=128):
    """Small audit arrays; density volume is never returned or archived."""
    vectors=minimum_vector(positions,initial['center'])
    radius=np.linalg.norm(vectors,axis=1)
    edges=np.linspace(0,.45,91)
    counts=np.histogram(radius,bins=edges)[0]
    volume=4*np.pi/3*np.diff(edges**3)
    enclosed=np.cumsum(counts)
    qr=initial['q_radius']
    qedges=np.linspace(0,.4,33)
    labels=np.searchsorted(qedges,qr,side='right')-1
    ratio=np.divide(radius,qr,out=np.zeros_like(radius),where=qr>0)
    yr=np.zeros((32,3))
    qcounts=np.zeros(32,dtype=int)
    for i in range(32):
        active=labels==i
        qcounts[i]=active.sum()
        yr[i]=np.percentile(ratio[active],[16,50,84]) if active.any() else np.zeros(3)
    # Same analysis mesh for all force/particle cases; projected CIC mass / mean.
    mass=cic_deposit(positions,analysis_mesh)
    projection=np.sum(mass,axis=2)*analysis_mesh**2
    return {'radius_edges_box':edges,'radial_particle_count':counts,
            'radial_density_over_mean':counts/(len(positions)*volume),
            'enclosed_particle_count':enclosed,
            'enclosed_density_over_mean':enclosed/(len(positions)*4*np.pi/3*edges[1:]**3),
            'q_edges_box':qedges,'q_shell_count':qcounts,'q_shell_y_percentiles':yr,
            'projection_density_over_mean':projection,
            'projection_axis_x_box':(np.arange(analysis_mesh)+.5)/analysis_mesh,
            'projection_axis_y_box':(np.arange(analysis_mesh)+.5)/analysis_mesh}
