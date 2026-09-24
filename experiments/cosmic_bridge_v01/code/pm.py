"""Periodic, equal-mass CIC particle-mesh gravity in a supplied FLRW background.

Units: the comoving box is one, tau=H0*t, p=a**2*dx/dtau.  The background
changes the drift/kick factors only.  There is no fifth force, varying
particle mass, double-well force, or claim of resolved halos/galaxies.
"""

from __future__ import annotations

import hashlib
from pathlib import Path

import numpy as np
from scipy import fft
from scipy.interpolate import PchipInterpolator


def file_sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


class Background:
    """Positive E(a), interpolated as log E versus log a without extrapolation."""

    def __init__(self, a, e, omega_m0, label="background", metadata=None):
        self.a = np.asarray(a, dtype=float)
        self.e = np.asarray(e, dtype=float)
        if (self.a.ndim != 1 or self.a.shape != self.e.shape or len(self.a) < 2
                or not np.all(np.isfinite(self.a)) or not np.all(np.isfinite(self.e))
                or np.any(self.a <= 0) or np.any(self.e <= 0)
                or np.any(np.diff(self.a) <= 0)):
            raise ValueError("Background requires ordered positive finite a and E arrays")
        self.omega_m0 = float(omega_m0)
        if not np.isfinite(self.omega_m0) or self.omega_m0 <= 0:
            raise ValueError("omega_m0 must be positive")
        self.label = str(label)
        self.metadata = {} if metadata is None else dict(metadata)
        self._interpolator = PchipInterpolator(np.log(self.a), np.log(self.e), extrapolate=False)
        self._gauss_nodes, self._gauss_weights = np.polynomial.legendre.leggauss(16)

    @classmethod
    def from_csv(cls, path, omega_m0, label="background"):
        table = np.genfromtxt(path, delimiter=",", names=True, dtype=float)
        fields = table.dtype.names or ()
        acol = "a" if "a" in fields else None
        ecol = "E" if "E" in fields else ("e" if "e" in fields else None)
        if acol is None or ecol is None:
            raise ValueError(f"Background CSV needs a,E columns, found {fields}")
        return cls(table[acol], table[ecol], omega_m0, label,
                   {"source_csv": str(path), "source_sha256": file_sha256(path)})

    @classmethod
    def eds(cls, a_min=0.02, a_max=1.0):
        a = np.geomspace(a_min, a_max, 65)
        return cls(a, a**(-1.5), 1.0, "Einstein-de Sitter",
                   {"analytic_E": "a^(-3/2)", "test_only": True})

    def E(self, a):
        values = np.asarray(a, dtype=float)
        tol = 32 * np.finfo(float).eps
        if np.any(values < self.a[0] * (1 - tol)) or np.any(values > self.a[-1] * (1 + tol)):
            raise ValueError("Requested a is outside background table")
        values = np.clip(values, self.a[0], self.a[-1])
        result = np.exp(self._interpolator(np.log(values)))
        return float(result) if result.ndim == 0 else result

    def factors(self, a1, a2):
        """Return integrals da/(a^2 E), da/(a^3 E), with signed intervals."""
        if a1 == a2:
            return 0.0, 0.0
        left, right = np.log(float(a1)), np.log(float(a2))
        nodes = (right + left) / 2 + (right - left) / 2 * self._gauss_nodes
        scale = np.exp(nodes)
        e = self.E(scale)
        weights = (right - left) / 2 * self._gauss_weights
        return (float(np.sum(weights / (scale * e))),
                float(np.sum(weights / (scale**2 * e))))


def _cic_weights(positions, nmesh):
    positions = np.asarray(positions, dtype=float)
    if positions.ndim != 2 or positions.shape[1] != 3 or not np.all(np.isfinite(positions)):
        raise ValueError("CIC positions must be a finite array with shape (particles,3)")
    scaled = np.remainder(positions, 1.0) * nmesh
    lower = np.floor(scaled).astype(np.int64)
    fraction = scaled - lower
    for dx in (0, 1):
        wx = fraction[:, 0] if dx else 1.0 - fraction[:, 0]
        ix = (lower[:, 0] + dx) % nmesh
        for dy in (0, 1):
            wy = fraction[:, 1] if dy else 1.0 - fraction[:, 1]
            iy = (lower[:, 1] + dy) % nmesh
            for dz in (0, 1):
                wz = fraction[:, 2] if dz else 1.0 - fraction[:, 2]
                iz = (lower[:, 2] + dz) % nmesh
                yield (ix * nmesh + iy) * nmesh + iz, wx * wy * wz


def cic_deposit(positions, nmesh):
    """Raw mass deposition: particle masses sum to one; never renormalized."""
    count = len(positions)
    if count == 0:
        raise ValueError("Need at least one particle")
    mass = np.zeros(nmesh**3, dtype=float)
    for index, weight in _cic_weights(positions, nmesh):
        mass += np.bincount(index, weights=weight / count, minlength=nmesh**3)
    return mass.reshape((nmesh,) * 3)


def cic_gather(grid, positions):
    """Same CIC shape as deposition; grid shape is (N,N,N) or (N,N,N,3)."""
    nmesh = grid.shape[0]
    components = grid.shape[3:]
    flat = grid.reshape((nmesh**3,) + components)
    result = np.zeros((len(positions),) + components, dtype=float)
    for index, weight in _cic_weights(positions, nmesh):
        if components:
            result += flat[index] * weight[:, None]
        else:
            result += flat[index] * weight
    return result


class PMGrid:
    def __init__(self, nmesh):
        self.nmesh = int(nmesh)
        k = 2 * np.pi * fft.fftfreq(nmesh, d=1.0 / nmesh)
        self.k = k
        self.k2 = (k[:, None, None]**2 + k[None, :, None]**2 + k[None, None, :]**2)
        self.inverse_k2 = np.zeros_like(self.k2)
        np.divide(1.0, self.k2, out=self.inverse_k2, where=self.k2 > 0)
        grad_k = k.copy()
        if nmesh % 2 == 0:
            grad_k[nmesh // 2] = 0.0
        self.gradient_k = (grad_k[:, None, None], grad_k[None, :, None], grad_k[None, None, :])

    def mesh_force(self, delta):
        spectrum = fft.fftn(delta, workers=1)
        force = np.empty(delta.shape + (3,), dtype=float)
        imaginary = 0.0
        for axis, wave in enumerate(self.gradient_k):
            transformed = fft.ifftn(1j * wave * self.inverse_k2 * spectrum, workers=1)
            imaginary = max(imaginary, float(np.max(np.abs(transformed.imag))))
            force[..., axis] = transformed.real
        return force, imaginary

    def force(self, positions):
        mass = cic_deposit(positions, self.nmesh)
        mass_sum = float(np.sum(mass))
        delta = mass * self.nmesh**3 - 1.0
        grid_force, imaginary = self.mesh_force(delta)
        particle_force = cic_gather(grid_force, positions)
        force_mean = np.mean(particle_force, axis=0)
        force_norm = float(np.mean(np.linalg.norm(particle_force, axis=1)))
        total_abs = float(np.linalg.norm(force_mean))
        diagnostics = {
            "raw_mass": mass_sum, "mass_relative_error": abs(mass_sum - 1.0),
            "density_contrast_mean": float(np.mean(delta)),
            "total_force_absolute": total_abs,
            "total_force_normalized": total_abs / max(force_norm, 1e-300),
            "mean_force_magnitude": force_norm,
            "force_ifft_max_imaginary": imaginary,
            "max_force_magnitude": float(np.max(np.linalg.norm(particle_force, axis=1))),
        }
        if not np.all(np.isfinite(particle_force)):
            raise FloatingPointError("Non-finite PM force; no clipping or repair applied")
        return particle_force, diagnostics


def advance_kdk(positions, momenta, a1, a2, background, mesh, initial_force=None):
    """One symmetric step; inputs are not mutated.  Endpoints may be reversed."""
    if initial_force is None:
        initial_force, _ = mesh.force(positions)
    midpoint = np.sqrt(a1 * a2)
    k1, _ = background.factors(a1, midpoint)
    k2, _ = background.factors(midpoint, a2)
    _, drift = background.factors(a1, a2)
    coefficient = 1.5 * background.omega_m0
    half_p = momenta + coefficient * k1 * initial_force
    step_displacement = drift * half_p
    following_x = np.remainder(positions + step_displacement, 1.0)
    final_force, force_diagnostics = mesh.force(following_x)
    following_p = half_p + coefficient * k2 * final_force
    diagnostics = dict(force_diagnostics)
    diagnostics["max_drift_cells"] = float(np.max(np.linalg.norm(step_displacement, axis=1)) * mesh.nmesh)
    if not np.all(np.isfinite(following_x)) or not np.all(np.isfinite(following_p)):
        raise FloatingPointError("Non-finite KDK state; no clipping or repair applied")
    return following_x, following_p, final_force, diagnostics


def generate_initial_modes(power_csv, seed=20260924, base_grid=32, box_mpc_h=50.0, cutoff_fraction=0.6):
    """One shared continuous Fourier field drawn from the supplied z=49 spectrum.

    Fourier coefficients obey E|delta_m|^2=P(2*pi*|m|/L)/L^3.  No realized
    variance rescaling is applied.  The returned coefficients are independent
    of the eventual particle/force mesh resolution.
    """
    table = np.genfromtxt(power_csv, delimiter=",", names=True, dtype=float)
    names = table.dtype.names or ()
    k_name = next((key for key in ("k_h_mpc", "k_h_Mpc") if key in names), None)
    p_name = next((key for key in ("pk_z49_mpc_over_h3", "P_lin_z49_Mpc_h3", "P_Mpc_h3") if key in names), None)
    if k_name is None or p_name is None:
        raise ValueError(f"Missing declared z49 spectrum columns in {names}")
    source_k, source_p = np.asarray(table[k_name]), np.asarray(table[p_name])
    if (source_k.ndim != 1 or source_k.shape != source_p.shape or len(source_k) < 2
            or not np.all(np.isfinite(source_k)) or not np.all(np.isfinite(source_p))
            or np.any(np.diff(source_k) <= 0) or np.any(source_k <= 0) or np.any(source_p <= 0)):
        raise ValueError("Power table must have finite increasing positive k and finite positive P")
    modes = fft.fftfreq(base_grid) * base_grid
    radius = np.sqrt(modes[:, None, None]**2 + modes[None, :, None]**2 + modes[None, None, :]**2)
    k = radius * 2 * np.pi / box_mpc_h
    kmax = cutoff_fraction * np.pi * base_grid / box_mpc_h
    retained = (k > 0) & (k <= kmax)
    if k[retained].min() < source_k.min() or k[retained].max() > source_k.max():
        raise ValueError("Power table does not cover retained modes; no extrapolation allowed")
    power = np.zeros_like(k)
    power[retained] = np.exp(np.interp(np.log(k[retained]), np.log(source_k), np.log(source_p)))
    white = np.random.default_rng(seed).standard_normal((base_grid,) * 3)
    coefficients = fft.fftn(white, workers=1) * np.sqrt(power / (box_mpc_h**3 * base_grid**3))
    coefficients[~retained] = 0.0
    if not np.all(np.isfinite(coefficients)):
        raise ValueError("Non-finite initial Fourier coefficient")
    metadata = {
        "seed": int(seed), "base_grid": int(base_grid), "box_mpc_h": float(box_mpc_h),
        "cutoff_fraction_base_nyquist": float(cutoff_fraction), "kmax_h_mpc": float(kmax),
        "retained_full_fft_modes": int(np.sum(retained)), "spectrum_redshift": 49,
        "source_csv": str(power_csv), "source_sha256": file_sha256(power_csv),
        "coefficient_sha256": hashlib.sha256(np.ascontiguousarray(coefficients).view(np.uint8)).hexdigest(),
        "normalization": "E|delta_m|^2 = P(k)/L^3; no realization rescaling",
        "realized_linear_density_variance": float(np.sum(np.abs(coefficients)**2)),
        "expected_linear_density_variance": float(np.sum(power) / box_mpc_h**3),
        "interpretation": "External standard matter transfer, not a fitted SMICA phase realization",
    }
    return {"coefficients": coefficients, "metadata": metadata}


def initial_conditions(modes, nmesh, a_initial, E_reference, f_reference):
    """Evaluate the shared bandlimited field at this mesh's particle centers."""
    coefficients = np.asarray(modes["coefficients"])
    if (coefficients.ndim != 3 or len(set(coefficients.shape)) != 1
            or not np.all(np.isfinite(coefficients))):
        raise ValueError("Initial coefficients must form a finite cubic three-dimensional array")
    if not np.all(np.isfinite([a_initial, E_reference, f_reference])) or a_initial <= 0 or E_reference <= 0:
        raise ValueError("Initial time and E must be positive and all velocity inputs finite")
    base = coefficients.shape[0]
    if nmesh < base:
        raise ValueError("Particle mesh must not be smaller than source mode grid")
    integers = np.rint(fft.fftfreq(base) * base).astype(int)
    ix, iy, iz = np.meshgrid(integers, integers, integers, indexing="ij")
    active = np.abs(coefficients) > 0
    mx, my, mz = ix[active], iy[active], iz[active]
    kvec = 2 * np.pi * np.stack([mx, my, mz], axis=-1)
    k2 = np.sum(kvec**2, axis=1)
    center_phase = np.exp(1j * np.pi * (mx + my + mz) / nmesh)
    index = (mx % nmesh, my % nmesh, mz % nmesh)
    displacement = np.empty((nmesh, nmesh, nmesh, 3))
    max_imaginary = 0.0
    for axis in range(3):
        embedded = np.zeros((nmesh,) * 3, dtype=complex)
        embedded[index] = 1j * kvec[:, axis] / k2 * coefficients[active] * center_phase
        field = fft.ifftn(embedded, workers=1) * nmesh**3
        if not np.all(np.isfinite(field)):
            raise ValueError("Non-finite initial displacement")
        imaginary = float(np.max(np.abs(field.imag)))
        hermitian_tolerance = 1e-14 + 1e-12 * float(np.max(np.abs(field.real)))
        if imaginary > hermitian_tolerance:
            raise ValueError(f"Initial field violates real/Hermitian condition: {imaginary} > {hermitian_tolerance}")
        displacement[..., axis] = field.real
        max_imaginary = max(max_imaginary, imaginary)
    grid = (np.arange(nmesh) + 0.5) / nmesh
    lagrangian = np.stack(np.meshgrid(grid, grid, grid, indexing="ij"), axis=-1).reshape((-1, 3))
    displacement = displacement.reshape((-1, 3))
    momentum_factor = float(a_initial)**2 * float(E_reference) * float(f_reference)
    metadata = dict(modes["metadata"])
    metadata.update({"nmesh": int(nmesh), "a_initial": float(a_initial),
                     "E_reference_initial": float(E_reference), "f_reference_initial": float(f_reference),
                     "common_momentum_factor": momentum_factor,
                     "initial_condition_order": "Zeldovich displacement and common reference growing velocity",
                     "max_ifft_imaginary_displacement": max_imaginary,
                     "displacement_rms_box_units": float(np.sqrt(np.mean(displacement**2))),
                     "max_displacement_cells": float(np.max(np.linalg.norm(displacement, axis=1)) * nmesh)})
    return {"positions": np.remainder(lagrangian + displacement, 1.0),
            "momenta": momentum_factor * displacement, "lagrangian": lagrangian,
            "displacement": displacement, "metadata": metadata}


def power_spectrum(delta, box_mpc_h=50.0, bin_edges=None):
    """Raw CIC-density auto power; full FFT mode count includes conjugate pairs."""
    n = delta.shape[0]
    modes = fft.fftfreq(n) * n
    radius = np.sqrt(modes[:, None, None]**2 + modes[None, :, None]**2 + modes[None, None, :]**2)
    k = radius * 2 * np.pi / box_mpc_h
    if bin_edges is None:
        bin_edges = np.arange(0.5, 17.5, 1.0) * 2 * np.pi / box_mpc_h
    edges = np.asarray(bin_edges, dtype=float)
    transformed = fft.fftn(delta, workers=1) / n**3
    mode_power = box_mpc_h**3 * np.abs(transformed)**2
    count, _ = np.histogram(k, bins=edges)
    psum, _ = np.histogram(k, bins=edges, weights=mode_power)
    ksum, _ = np.histogram(k, bins=edges, weights=k)
    p = np.divide(psum, count, out=np.full(len(count), np.nan), where=count > 0)
    center = np.divide(ksum, count, out=np.full(len(count), np.nan), where=count > 0)
    return {"k": center, "pk": p, "count": count, "edges": edges,
            "units": {"k": "h/Mpc", "pk": "(Mpc/h)^3"},
            "estimator": "raw CIC density; no window deconvolution or shot-noise subtraction",
            "mode_count_note": "Full Fourier count includes both members of conjugate pairs; not independent samples"}


def integrate(initial, background, nmesh, steps=256, a_final=1.0,
              a_outputs=(0.02, 0.1, 0.25, 0.5, 1.0), box_mpc_h=50.0,
              bin_edges=None, retain_particles=False):
    """Return states/diagnostics; no files, fits, normalization repairs, or plotting.

    Non-grid output times use a temporary step from the preceding main state.
    They do not replace or subdivide the main fixed-log-a trajectory.
    """
    x = np.array(initial["positions"], dtype=float, copy=True)
    p = np.array(initial["momenta"], dtype=float, copy=True)
    if x.shape != p.shape or x.ndim != 2 or x.shape[1] != 3 or not np.all(np.isfinite(x)) or not np.all(np.isfinite(p)):
        raise ValueError("Initial positions and momenta must be equally shaped finite (N,3) arrays")
    a_initial = float(initial["metadata"]["a_initial"])
    if not a_final > a_initial:
        raise ValueError("Final scale factor must exceed initial scale factor")
    boundaries = np.geomspace(a_initial, float(a_final), int(steps) + 1)
    requested = np.asarray(sorted(set(float(a) for a in a_outputs)))
    if np.any(requested < a_initial - 1e-14) or np.any(requested > a_final + 1e-14):
        raise ValueError("Output time outside integration interval")
    mesh = PMGrid(nmesh)
    force, start_diagnostics = mesh.force(x)
    initial_mean_p = np.mean(p, axis=0)
    history = []
    snapshots = []
    output_index = 0

    def snapshot(a, positions, momenta, kind):
        raw_mass = cic_deposit(positions, nmesh)
        delta = raw_mass * nmesh**3 - 1.0
        item = {"a": float(a), "z": float(1.0 / a - 1.0), "delta": delta,
                "projection": np.mean(1.0 + delta, axis=2),
                "spectrum": power_spectrum(delta, box_mpc_h, bin_edges),
                "output_step_kind": kind, "raw_mass": float(np.sum(raw_mass)),
                "density_variance": float(np.mean(delta**2)),
                "density_minimum": float(np.min(1.0 + delta)),
                "density_maximum": float(np.max(1.0 + delta))}
        if retain_particles:
            item["positions"] = positions.copy()
            item["momenta"] = momenta.copy()
        return item

    def record(a, diag, momenta):
        row = dict(diag)
        mean_p = np.mean(momenta, axis=0)
        row.update({"a": float(a), "mean_momentum": mean_p.tolist(),
                    "momentum_change_absolute": float(np.linalg.norm(mean_p - initial_mean_p)),
                    "momentum_rms": float(np.sqrt(np.mean(momenta**2))),
                    "max_drift_cells": float(diag.get("max_drift_cells", 0.0))})
        return row

    history.append(record(a_initial, start_diagnostics, p))
    while output_index < len(requested) and abs(requested[output_index] - a_initial) < 1e-13:
        snapshots.append(snapshot(requested[output_index], x, p, "initial"))
        output_index += 1
    for left, right in zip(boundaries[:-1], boundaries[1:]):
        while output_index < len(requested) and requested[output_index] < right - 1e-13:
            output_a = requested[output_index]
            out_x, out_p, _, _ = advance_kdk(x, p, left, output_a, background, mesh, force)
            snapshots.append(snapshot(output_a, out_x, out_p, "temporary_KDK_from_previous_main_state"))
            output_index += 1
        x, p, force, diag = advance_kdk(x, p, left, right, background, mesh, force)
        history.append(record(right, diag, p))
        while output_index < len(requested) and abs(requested[output_index] - right) < 1e-13:
            snapshots.append(snapshot(requested[output_index], x, p, "main_fixed_log_a_state"))
            output_index += 1
    if output_index != len(requested):
        raise RuntimeError("Not all requested output times were produced")
    diagnostics = {"max_mass_relative_error": max(r["mass_relative_error"] for r in history),
                   "max_total_force_normalized": max(r["total_force_normalized"] for r in history),
                   "max_total_force_absolute": max(r["total_force_absolute"] for r in history),
                   "max_momentum_change_absolute": max(r["momentum_change_absolute"] for r in history),
                   "max_drift_cells": max(r["max_drift_cells"] for r in history),
                   "max_force_ifft_imaginary": max(r["force_ifft_max_imaginary"] for r in history),
                   "finite": bool(np.all(np.isfinite(x)) and np.all(np.isfinite(p)))}
    metadata = {"background_label": background.label, "background": background.metadata,
                "omega_m0": background.omega_m0, "nmesh": int(nmesh), "particles": len(x),
                "steps": int(steps), "a_initial": a_initial, "a_final": float(a_final),
                "box_mpc_h": float(box_mpc_h), "mesh_length_mpc_h": float(box_mpc_h / nmesh),
                "integrator": "KDK, fixed logarithmic scale-factor boundaries, 16-point factor quadrature",
                "initial_conditions": dict(initial["metadata"]),
                "code_sha256": file_sha256(__file__),
                "gravity": "Periodic Newtonian PM, CIC deposit/gather, spectral Poisson and Nyquist-safe gradient",
                "interpretation": "Collisionless density clustering, not resolved halos or galaxies"}
    return {"snapshots": snapshots, "final_positions": x, "final_momenta": p,
            "history": history, "diagnostics": diagnostics, "metadata": metadata}
