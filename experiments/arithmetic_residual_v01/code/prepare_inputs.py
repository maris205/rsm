#!/usr/bin/env python3
"""Prepare fixed arithmetic inputs and conditional IAAFT controls.

This is mathematical source data, not astronomical measurements. Evaluation
blocks are held out from normalization, not blind forecasts: each surrogate is
conditioned on the complete block's values and Fourier amplitudes.
"""
from __future__ import annotations

import csv
import hashlib
import io
import json
import time
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import scipy
from scipy.interpolate import CubicSpline

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
RESULTS = ROOT / "results"
CACHE_URL = "https://raw.githubusercontent.com/maris205/riemann_logistic/e3419dda3d5515afd91a5dd6e8d8b8398f359dde/python/riemann_10k_true.npy"
CACHE_SHA = "6bb18ef724b76280ee91e56c9ad0669a8026c996bd969c710b3a708f28e57545"
OFFICIAL_URL = "https://www-users.cse.umn.edu/~odlyzko/zeta_tables/zeros1"
SEED_BASE = 20260922
N_SURROGATES = 99
BLOCK_SIZE = 128
N_BLOCKS = 8
MAX_ITERATIONS = 2000
QUALITY_LIMIT = 0.05
EPSILON = 0.02


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def acquire(url: str, path: Path, expected_sha: str | None = None) -> dict:
    if not path.exists():
        request = urllib.request.Request(url, headers={"User-Agent": "Research-reproducibility/1.0"})
        with urllib.request.urlopen(request, timeout=90) as response:
            body = response.read()
        temporary = path.with_suffix(path.suffix + ".part")
        temporary.write_bytes(body)
        temporary.replace(path)
    actual_sha = sha256(path)
    if expected_sha is not None and actual_sha != expected_sha:
        raise ValueError(f"SHA256 mismatch for {path.name}: {actual_sha}")
    return {"url": url, "path": str(path.relative_to(ROOT)), "sha256": actual_sha,
            "bytes": path.stat().st_size, "expected_sha256": expected_sha}


def project_spectrum(sequence: np.ndarray, target_fft: np.ndarray) -> np.ndarray:
    """Real Fourier amplitude projection, preserving the exact target mean."""
    spectrum = np.fft.rfft(sequence)
    projected = np.abs(target_fft) * np.exp(1j * np.angle(spectrum))
    projected[0] = target_fft[0]
    # For an even sample count, the Nyquist coefficient must be real.
    projected[-1] = np.copysign(abs(target_fft[-1]), spectrum[-1].real)
    return np.fft.irfft(projected, n=len(sequence))


def amplitude_error(sequence: np.ndarray, target_amplitude: np.ndarray) -> float:
    return float(np.linalg.norm(np.abs(np.fft.rfft(sequence))[1:] - target_amplitude[1:])
                 / np.linalg.norm(target_amplitude[1:]))


def iaaft(reference: np.ndarray, seed: int) -> tuple[np.ndarray, np.ndarray, dict]:
    """Keep the best exact-rank iterate; never optimize an output observable."""
    rng = np.random.default_rng(seed)
    ordered = np.sort(reference, kind="mergesort")
    target_fft = np.fft.rfft(reference)
    target_amplitude = np.abs(target_fft)
    current = rng.permutation(reference)
    best = current.copy()
    best_error = amplitude_error(best, target_amplitude)
    best_iteration = 0
    stop_reason = "max_iterations"
    previous_order = np.argsort(current, kind="mergesort")
    for iteration in range(1, MAX_ITERATIONS + 1):
        projected = project_spectrum(current, target_fft)
        order = np.argsort(projected, kind="mergesort")
        candidate = np.empty_like(current)
        candidate[order] = ordered
        error = amplitude_error(candidate, target_amplitude)
        if error < best_error:
            best = candidate.copy()
            best_error = error
            best_iteration = iteration
        if np.array_equal(order, previous_order):
            stop_reason = "unchanged_rank_permutation"
            break
        current = candidate
        previous_order = order
    companion = project_spectrum(best, target_fft)
    return best, companion, {"seed": seed, "iterations": iteration,
                            "best_iteration": best_iteration, "stop_reason": stop_reason}


def relative_error(values: np.ndarray, reference: np.ndarray) -> float:
    denominator = np.linalg.norm(reference)
    return float(np.linalg.norm(values - reference) / denominator) if denominator else 0.0


def acf(sequence: np.ndarray) -> np.ndarray:
    centered = sequence - np.mean(sequence)
    result = np.fft.irfft(np.abs(np.fft.rfft(centered)) ** 2, n=len(sequence)) / len(sequence)
    return result / result[0] if result[0] else result


def basic_quality(sequence: np.ndarray, reference: np.ndarray) -> dict:
    amplitude = np.abs(np.fft.rfft(sequence))[1:]
    ref_amplitude = np.abs(np.fft.rfft(reference))[1:]
    sorted_difference = np.sort(sequence) - np.sort(reference)
    return {
        "fft_amplitude_relative_l2_no_dc": relative_error(amplitude, ref_amplitude),
        "fft_power_relative_l2_no_dc": relative_error(amplitude**2, ref_amplitude**2),
        "sorted_max_abs_difference": float(np.max(np.abs(sorted_difference))),
        "sorted_relative_l2_centered_denominator": float(np.linalg.norm(sorted_difference)
                     / np.linalg.norm(reference - np.mean(reference))) if np.std(reference) else 0.0,
        "periodic_acf_relative_l2": relative_error(acf(sequence), acf(reference)),
        "mean": float(np.mean(sequence)), "variance_ddof0": float(np.var(sequence)),
        "mean_difference": float(np.mean(sequence) - np.mean(reference)),
        "variance_relative_difference": float(np.var(sequence) / np.var(reference) - 1.0)
                if np.var(reference) else 0.0,
    }


def dense_functions(spline: CubicSpline, dense_nodes: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    angle = np.pi * (dense_nodes - 1.0) / 2.0
    window = np.sin(angle)**4
    window_derivative = 2.0 * np.pi * np.sin(angle)**3 * np.cos(angle)
    values = spline(dense_nodes)
    return window * values, window_derivative * values + window * spline(dense_nodes, 1)


def dense_quality(values: np.ndarray, derivative: np.ndarray,
                  ref_values: np.ndarray, ref_derivative: np.ndarray) -> dict:
    result = {
        "dense_windowed_mean": float(np.mean(values)),
        "dense_windowed_rms": float(np.sqrt(np.mean(values**2))),
        "dense_derivative_mean": float(np.mean(derivative)),
        "dense_derivative_rms": float(np.sqrt(np.mean(derivative**2))),
        "dense_windowed_mean_difference": float(np.mean(values) - np.mean(ref_values)),
        "dense_windowed_rms_relative_difference": float(np.linalg.norm(values) / np.linalg.norm(ref_values) - 1.0)
            if np.linalg.norm(ref_values) else 0.0,
        "dense_derivative_rms_relative_difference": float(np.linalg.norm(derivative) / np.linalg.norm(ref_derivative) - 1.0)
            if np.linalg.norm(ref_derivative) else 0.0,
        "dense_windowed_max_abs": float(np.max(np.abs(values))),
        "dense_derivative_max_abs": float(np.max(np.abs(derivative))),
    }
    for label, vector, reference in (("windowed", values, ref_values), ("derivative", derivative, ref_derivative)):
        p = np.abs(np.fft.rfft(vector))[1:]**2
        ref_p = np.abs(np.fft.rfft(reference))[1:]**2
        result[f"dense_{label}_fft_power_relative_l2_no_dc"] = relative_error(p, ref_p)
        result[f"dense_{label}_periodic_acf_relative_l2"] = relative_error(acf(vector), acf(reference))
    return result


def main() -> None:
    started = time.perf_counter()
    (DATA / "raw").mkdir(parents=True, exist_ok=True)
    RESULTS.mkdir(parents=True, exist_ok=True)
    cache_path = DATA / "raw" / "riemann_10k_true.npy"
    official_path = DATA / "raw" / "odlyzko_zeros1.txt"
    with ThreadPoolExecutor(max_workers=2) as executor:
        futures = [executor.submit(acquire, CACHE_URL, cache_path, CACHE_SHA),
                   executor.submit(acquire, OFFICIAL_URL, official_path)]
        source_metadata = [future.result() for future in futures]
    zeros = np.load(cache_path, allow_pickle=False)
    if zeros.shape != (10000,) or zeros.dtype != np.float64 or not np.all(np.diff(zeros) > 0):
        raise ValueError("Expected 10000 strictly increasing float64 positive-zero ordinates")
    official = np.loadtxt(io.StringIO(official_path.read_text())).reshape(-1)
    if len(official) < 10000:
        raise ValueError("Official table has fewer than 10000 ordinates")
    official_max_difference = float(np.max(np.abs(zeros - official[:10000])))
    if official_max_difference > 3e-9:
        raise ValueError(f"Official cross-check failed: {official_max_difference}")
    local_path = DATA / "raw/mathematical_reference_zeros.csv"
    with local_path.open() as source:
        local_rows = list(csv.DictReader(source))
    local_indices = np.array([int(row["zero_index"]) for row in local_rows])
    local_zeros = np.array([float(row["ordinate_decimal"]) for row in local_rows])
    local_max_difference = float(np.max(np.abs(zeros[local_indices - 1] - local_zeros)))
    if local_max_difference > 3e-9:
        raise ValueError(f"Local mathematical-reference cross-check failed: {local_max_difference}")
    source_metadata.append({"path": str(local_path.relative_to(ROOT)), "sha256": sha256(local_path),
                            "original_workspace_path": "riemann_clock/data/processed/mathematical_reference_zeros.csv",
                            "rows": len(local_rows), "method": "mpmath.zetazero (existing project data)"})
    unfolded = zeros / (2 * np.pi) * np.log(zeros / (2 * np.pi)) - zeros / (2 * np.pi) + 7 / 8
    gaps = np.diff(unfolded)
    raw_residual = gaps - 1.0
    training_gap_indices = np.arange(1001, 2025)
    training = raw_residual[training_gap_indices - 1]
    mu, sigma = float(np.mean(training)), float(np.std(training, ddof=0))
    eta_all = np.tanh((raw_residual - mu) / (3 * sigma))
    block_gap_indices = np.array([np.arange(4097 + BLOCK_SIZE * block, 4097 + BLOCK_SIZE * (block + 1))
                                 for block in range(N_BLOCKS)])
    np.savez_compressed(DATA / "full_zero_gaps.npz", zeros=zeros,
                        zero_indices=np.arange(1, 10001), unfolded_zeros=unfolded,
                        gap_indices=np.arange(1, 10000), d=gaps, raw_r=raw_residual, eta=eta_all,
                        training_gap_indices=training_gap_indices,
                        training_d=gaps[training_gap_indices - 1], training_r=training,
                        training_eta=eta_all[training_gap_indices - 1], block_gap_indices=block_gap_indices,
                        block_d=gaps[block_gap_indices - 1], block_r=raw_residual[block_gap_indices - 1],
                        block_eta=eta_all[block_gap_indices - 1], normalization_mu=mu, normalization_sigma=sigma)
    nodes = np.linspace(1.0, 3.0, BLOCK_SIZE)
    dense_nodes = np.linspace(1.0, 3.0, 4097)
    sequences, coefficients, families, blocks, surrogate_numbers, identifiers, quality_rows = [], [], [], [], [], [], []
    max_dense_amplitude = 0.0
    for block in range(-1, N_BLOCKS):
        if block == -1:
            reference = np.zeros(BLOCK_SIZE)
            candidates = [("baseline", 0, reference, {})]
        else:
            reference = eta_all[block_gap_indices[block] - 1]
            primary, companions = [], []
            for j in range(1, N_SURROGATES + 1):
                exact_rank, exact_spectrum, details = iaaft(reference, SEED_BASE + 1000 * block + j)
                primary.append(("iaaft", j, exact_rank, details))
                companions.append(("exact_spectrum", j, exact_spectrum, details))
            candidates = [("real", 0, reference, {})] + primary + companions
        reference_spline = CubicSpline(nodes, reference, bc_type="natural")
        reference_dense, reference_derivative = dense_functions(reference_spline, dense_nodes)
        for family, number, sequence, algorithm_details in candidates:
            case_index = len(sequences)
            identifier = "baseline" if block < 0 else f"block{block + 1:02d}_{family}_{number:03d}"
            spline = CubicSpline(nodes, sequence, bc_type="natural")
            dense, derivative = dense_functions(spline, dense_nodes)
            quality = basic_quality(sequence, reference)
            passes = quality["fft_amplitude_relative_l2_no_dc"] <= QUALITY_LIMIT and quality["sorted_max_abs_difference"] <= 1e-14
            row = {"case_index": case_index, "case_id": identifier, "family": family, "block": block,
                   "surrogate": number, "seed": algorithm_details.get("seed", 0),
                   "iterations": algorithm_details.get("iterations", 0),
                   "best_iteration": algorithm_details.get("best_iteration", 0),
                   "stop_reason": algorithm_details.get("stop_reason", "not_applicable"),
                   "primary_node_quality_pass": bool(passes), **quality,
                   **dense_quality(dense, derivative, reference_dense, reference_derivative)}
            quality_rows.append(row)
            sequences.append(sequence)
            coefficients.append(spline.c)
            families.append(family)
            blocks.append(block)
            surrogate_numbers.append(number)
            identifiers.append(identifier)
            max_dense_amplitude = max(max_dense_amplitude, row["dense_windowed_max_abs"])
        print(f"prepared block {block + 1}/{N_BLOCKS}, cases={len(sequences)}", flush=True)
    quality_pass = np.array([row["primary_node_quality_pass"] for row in quality_rows])
    np.savez_compressed(DATA / "inputs.npz", eta=np.stack(sequences), spline_coefficients=np.stack(coefficients),
                        nodes=nodes, case_family=np.array(families), case_block=np.array(blocks, dtype=np.int64),
                        case_surrogate=np.array(surrogate_numbers, dtype=np.int64), case_id=np.array(identifiers),
                        primary_node_quality_pass=quality_pass, epsilon=np.array(EPSILON),
                        block_gap_indices=block_gap_indices,
                        normalization_mu=np.array(mu), normalization_sigma=np.array(sigma))
    with (RESULTS / "input_quality.csv").open("w", newline="") as output:
        writer = csv.DictWriter(output, fieldnames=list(quality_rows[0]))
        writer.writeheader()
        writer.writerows(quality_rows)
    primary_rows = [row for row in quality_rows if row["family"] == "iaaft"]
    companion_rows = [row for row in quality_rows if row["family"] == "exact_spectrum"]
    block_summaries = []
    for block in range(N_BLOCKS):
        rows = [row for row in primary_rows if row["block"] == block]
        block_summaries.append({"block": block, "gap_start_1based": int(block_gap_indices[block, 0]),
             "gap_end_1based": int(block_gap_indices[block, -1]),
             "primary_count": len(rows), "primary_quality_pass_count": sum(row["primary_node_quality_pass"] for row in rows),
             "amplitude_error_min": min(row["fft_amplitude_relative_l2_no_dc"] for row in rows),
             "amplitude_error_max": max(row["fft_amplitude_relative_l2_no_dc"] for row in rows)})
    checks = {
        "cache_expected_sha256": source_metadata[0]["sha256"] == CACHE_SHA,
        "zero_count_10000_float64_increasing": bool(zeros.shape == (10000,) and zeros.dtype == np.float64 and np.all(np.diff(zeros) > 0)),
        "official_10000_agreement_le_3e_9": official_max_difference <= 3e-9,
        "local_reference_agreement_le_3e_9": local_max_difference <= 3e-9,
        "local_reference_has_83_rows": len(local_rows) == 83,
        "normalization_training_disjoint_from_evaluation": bool(np.intersect1d(training_gap_indices, block_gap_indices).size == 0),
        "all_primary_sorted_max_difference_le_1e_14": max(row["sorted_max_abs_difference"] for row in primary_rows) <= 1e-14,
        "all_companion_amplitude_errors_le_1e_12": max(row["fft_amplitude_relative_l2_no_dc"] for row in companion_rows) <= 1e-12,
        "canonical_case_count_1593": len(sequences) == 1593,
        "all_cases_finite": bool(np.isfinite(np.stack(sequences)).all() and np.isfinite(np.stack(coefficients)).all()),
    }
    summary = {
        "experiment": "arithmetic_residual_v01", "script": "code/prepare_inputs.py",
        "data_kind": "numerical Riemann zero ordinates; mathematical reference, not physical observations",
        "case_count": len(sequences), "eta_shape": list(np.stack(sequences).shape),
        "spline_coefficients_shape": list(np.stack(coefficients).shape),
        "case_order": "baseline; for each block 0..7: real, iaaft 1..99, exact_spectrum 1..99",
        "family_counts": {name: families.count(name) for name in sorted(set(families))},
        "sources": source_metadata,
        "official_crosscheck": {"compared_count": 10000, "official_table_count": len(official),
                                 "max_abs_difference": official_max_difference, "threshold": 3e-9},
        "local_crosscheck": {"compared_count": len(local_rows), "max_abs_difference": local_max_difference},
        "normalization": {"training_gap_start_1based": 1001, "training_gap_end_1based": 2024,
                           "count": len(training), "mu_raw_r": mu, "sigma_raw_r_ddof0": sigma,
                           "eta_formula": "tanh((r-mu)/(3*sigma))", "unfolding": "gamma/(2*pi)*log(gamma/(2*pi))-gamma/(2*pi)+7/8"},
        "surrogates": {"per_block_per_family": N_SURROGATES, "seed_formula": "20260922 + 1000*block_0based + j_1based",
                       "max_iterations": MAX_ITERATIONS, "selection": "smallest non-DC relative FFT-amplitude L2 among exact-rank iterations, including initialization",
                       "early_stop": "unchanged rank permutation", "sort_kind": "mergesort",
                       "primary_quality_amplitude_threshold": QUALITY_LIMIT, "primary_quality_sorted_max_threshold": 1e-14,
                       "primary_quality_pass_count": sum(row["primary_node_quality_pass"] for row in primary_rows),
                       "primary_quality_total": len(primary_rows),
                       "failed_cases_retained": [row["case_id"] for row in primary_rows if not row["primary_node_quality_pass"]],
                       "companion_definition": "one Fourier amplitude projection of corresponding best exact-rank sequence; same Fourier phases, exact spectrum, approximate value distribution",
                       "interpretation": "conditional structure controls, not independent zero predictions; no output-aware surrogate selection"},
        "interpolation": {"type": "natural cubic spline", "nodes": 128, "range": [1.0, 3.0],
                          "linear_in_node_values": True, "coefficient_convention": "SciPy CubicSpline.c; descending local powers 3,2,1,0",
                          "epsilon": EPSILON, "window": "sin(pi*(R-1)/2)^4 inside [1,3], zero outside",
                          "dense_diagnostic_nodes": len(dense_nodes),
                          "max_dense_abs_windowed_eta_all_cases": max_dense_amplitude,
                          "max_dense_abs_delta_mass_squared_all_cases": EPSILON * max_dense_amplitude,
                          "diagnostics_only_no_selection": True},
        "blocks": block_summaries, "checks": checks,
        "versions": {"numpy": np.__version__, "scipy": scipy.__version__},
        "elapsed_seconds": time.perf_counter() - started,
    }
    summary["artifact_sha256"] = {str(path.relative_to(ROOT)): sha256(path) for path in
         [Path(__file__), DATA / "inputs.npz", DATA / "full_zero_gaps.npz", RESULTS / "input_quality.csv"]}
    (DATA / "source_metadata.json").write_text(json.dumps({"sources": source_metadata,
                 "warning": summary["data_kind"], "cross_checks": {"official": summary["official_crosscheck"],
                 "local": summary["local_crosscheck"]}}, indent=2) + "\n")
    (RESULTS / "input_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({"case_count": len(sequences), "primary_quality_pass_count": summary["surrogates"]["primary_quality_pass_count"],
                      "primary_quality_total": len(primary_rows), "official_max_abs_difference": official_max_difference,
                      "local_max_abs_difference": local_max_difference, "checks": checks,
                      "elapsed_seconds": summary["elapsed_seconds"]}, indent=2), flush=True)
    if not all(checks.values()):
        raise SystemExit("A required source or construction check failed; see input_summary.json")


if __name__ == "__main__":
    main()
