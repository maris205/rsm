#!/usr/bin/env python3
"""Small deterministic checks of the v0.1 equations; no observational fit.

Run from any directory.  Outputs a JSON record beside the mathematical review.
Requires NumPy.  All matrices are at most 14 by 14.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np


def step_matrix(k: np.ndarray, h: float) -> np.ndarray:
    eye = np.eye(k.shape[0])
    return np.block([[eye - h * h * k, h * eye], [-h * k, eye]])


def scalar_step(k: float, h: float = 1.0) -> np.ndarray:
    return step_matrix(np.array([[k]], dtype=float), h)


def energy(x: np.ndarray, k: np.ndarray) -> float:
    q, p = np.split(x, 2)
    return float((p @ p + q @ k @ q) / 2.0)


def modified_metric(k: np.ndarray, h: float) -> np.ndarray:
    return np.block([[k, -h * k / 2.0], [-h * k / 2.0, np.eye(len(k))]])


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).resolve().parents[1]
        / "reports"
        / "minimal_dynamics_v01_checks.json",
    )
    args = parser.parse_args()
    checks: list[dict] = []

    def residual_check(name: str, residual: float, tolerance: float) -> None:
        checks.append(
            {
                "name": name,
                "residual": float(residual),
                "tolerance": float(tolerance),
                "passed": bool(np.isfinite(residual) and residual <= tolerance),
            }
        )

    def condition_check(name: str, condition: bool, **details: object) -> None:
        checks.append({"name": name, "passed": bool(condition), **details})

    # A symmetric positive definite example distinct from the scalar examples.
    rng = np.random.default_rng(271828)
    a = rng.normal(size=(4, 4))
    k0 = a.T @ a + 0.7 * np.eye(4)
    h = 0.13
    m = step_matrix(k0, h)
    z = np.zeros_like(k0)
    eye = np.eye(4)
    j = np.block([[z, eye], [-eye, z]])
    residual_check("symmetric_k_symplectic", np.max(np.abs(m.T @ j @ m - j)), 2e-13)
    g = modified_metric(k0, h)
    residual_check("frozen_modified_energy", np.max(np.abs(m.T @ g @ m - g)), 2e-12)
    condition_check(
        "frozen_modified_metric_positive",
        np.linalg.eigvalsh(g).min() > 0 and h * h * np.linalg.eigvalsh(k0).max() < 4,
        minimum_metric_eigenvalue=float(np.linalg.eigvalsh(g).min()),
    )

    # Compare actual matrix eigenphases with the analytic arcsine expression.
    for k in (0.25, 1.0, 3.0, 3.9):
        eig = np.linalg.eigvals(scalar_step(k))
        theta = 2.0 * np.arcsin(np.sqrt(k) / 2.0)
        phase_error = abs(float(np.max(np.angle(eig))) - theta)
        residual_check(f"frozen_phase_k_{k:g}", phase_error, 3e-14)
        residual_check(f"frozen_unit_modulus_k_{k:g}", np.max(np.abs(np.abs(eig) - 1)), 3e-14)

    # Endpoints require separate treatment even though spectral radius is one.
    initial = np.array([0.0, 1.0])
    free20 = np.linalg.matrix_power(scalar_step(0.0), 20) @ initial
    residual_check("zero_mode_linear_drift", np.max(np.abs(free20 - [20.0, 1.0])), 0.0)
    edge20 = np.linalg.matrix_power(scalar_step(4.0), 20) @ initial
    residual_check("upper_endpoint_jordan_drift", np.max(np.abs(edge20 - [-20.0, -39.0])), 0.0)
    for k in (-1.0, 4.1):
        rho = float(np.max(np.abs(np.linalg.eigvals(scalar_step(k)))))
        condition_check(f"outside_stable_interval_k_{k:g}", rho > 1, spectral_radius=rho)

    # Integer-matrix counterexample to arbitrary time-varying stability.
    m1, m3 = scalar_step(1.0), scalar_step(3.0)
    product = m3 @ m1
    expected = np.array([[-1.0, -1.0], [-1.0, -2.0]])
    residual_check("alternating_product", np.max(np.abs(product - expected)), 0.0)
    rho = float(np.max(np.abs(np.linalg.eigvals(product))))
    residual_check("alternating_spectral_radius", abs(rho - (3 + np.sqrt(5)) / 2), 2e-14)
    j2 = np.array([[0.0, 1.0], [-1.0, 0.0]])
    for k, matrix in ((1.0, m1), (3.0, m3)):
        residual_check(f"integer_symplectic_k_{k:g}", np.max(np.abs(matrix.T @ j2 @ matrix - j2)), 0.0)
        metric = modified_metric(np.array([[k]]), 1.0)
        residual_check(f"integer_frozen_invariant_k_{k:g}", np.max(np.abs(matrix.T @ metric @ matrix - metric)), 0.0)

    # Exact bookkeeping at finite step size, with nonzero driven and numerical terms.
    x0 = rng.normal(size=8)
    delta_k = np.diag([0.15, -0.07, 0.08, -0.04])
    k1 = k0 + delta_k
    x1 = m @ x0
    q1 = x1[:4]
    frozen_increment = energy(x1, k0) - energy(x0, k0)
    switch_work = float(q1 @ delta_k @ q1 / 2)
    total_increment = energy(x1, k1) - energy(x0, k0)
    residual_check("discrete_energy_bookkeeping", abs(total_increment - frozen_increment - switch_work), 3e-13)
    condition_check("symplectic_step_need_not_preserve_H", abs(frozen_increment) > 1e-6, frozen_energy_increment=frozen_increment)

    # Independent central directional derivative of H along the continuous equations.
    q, p = np.split(x0, 2)
    tangent = np.concatenate([p, -k0 @ q])
    eps = 1e-6
    numerical_rate = (
        energy(x0 + eps * tangent, k0 + eps * delta_k)
        - energy(x0 - eps * tangent, k0 - eps * delta_k)
    ) / (2 * eps)
    analytic_rate = float(q @ delta_k @ q / 2)
    residual_check("continuous_external_work_directional_derivative", abs(numerical_rate - analytic_rate), 3e-8)

    # Path graph: no state Jacobian entry can jump more than one edge per step.
    nodes = 7
    adjacency = np.zeros((nodes, nodes))
    for index in range(nodes - 1):
        adjacency[index, index + 1] = adjacency[index + 1, index] = 1.0
    laplacian = np.diag(adjacency.sum(axis=1)) - adjacency
    graph_k = np.diag(np.linspace(0.7, 1.3, nodes)) + 0.3 * laplacian
    graph_step = step_matrix(graph_k, 0.2)
    distances = np.abs(np.arange(nodes)[:, None] - np.arange(nodes)[None, :])
    state_distances = np.tile(distances, (2, 2))
    for steps in (1, 2, 3):
        jacobian = np.linalg.matrix_power(graph_step, steps)
        forbidden = np.abs(jacobian[state_distances > steps])
        residual_check(f"local_dependence_at_{steps}_steps", float(forbidden.max()), 0.0)
    condition_check("local_test_has_nonzero_neighbour_transfer", abs(graph_step[0, 1]) > 0, transfer=float(graph_step[0, 1]))

    report = {
        "scope": "Deterministic small-matrix verification of the v0.1 model equations; not a data fit or a stability proof for the logarithmic schedule.",
        "seed": 271828,
        "numpy_version": np.__version__,
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "check_count": len(checks),
        "passed_count": sum(item["passed"] for item in checks),
        "all_passed": all(item["passed"] for item in checks),
        "checks": checks,
        "illustrative_values": {
            "alternating_spectral_radius": rho,
            "zero_mode_state_after_20_steps": free20.tolist(),
            "upper_endpoint_state_after_20_steps": edge20.tolist(),
            "frozen_energy_increment": frozen_increment,
            "parameter_switch_work": switch_work,
            "total_energy_increment": total_increment,
            "continuous_work_rate": analytic_rate,
        },
        "limits": [
            "Frozen-system stability is distinct from stability of a time-dependent schedule.",
            "Only fixed symmetric state-independent K is used in the symplectic checks.",
            "Local dependence is a graph statement, not a proof of relativistic symmetry.",
            "The atomic response is an assumed interface, not a microscopic derivation.",
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps({"all_passed": report["all_passed"], "passed": report["passed_count"], "checks": report["check_count"], "output": str(args.output)}, ensure_ascii=False))
    if not report["all_passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
