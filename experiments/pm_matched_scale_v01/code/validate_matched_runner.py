#!/usr/bin/env python3
"""Small, independent consistency checks for the matched-resolution runner.

No 64^3/128^3 production trajectory or production force benchmark is executed.
The new protocol must be frozen before this script may run. Archived experiment
directories are read only, including suppression of implicit bytecode writes.
"""

from __future__ import annotations

import sys

sys.dont_write_bytecode = True

import argparse
import ast
import hashlib
import importlib.util
import json
import platform
import tempfile
import time
import traceback
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import scipy


ROOT = Path(__file__).resolve().parents[1]
EXPERIMENTS = ROOT.parent
ARCHIVES = (EXPERIMENTS / "halo_cooling_v01", EXPERIMENTS / "cosmic_bridge_v01")
EXPECTED_SOURCE_HASHES = {
    "halo_cooling_v01/code/pm_refined.py": "f8e5bcbc32d726a0345282e228458894ac3b5e75baf0051f65cdbfb05a268602",
    "halo_cooling_v01/code/spherical_collapse.py": "4e0afba9b783575ec02ad60ec2f6acfaf87da8a0824f6a978b115e89ba522e47",
    "halo_cooling_v01/code/run_spherical_pm.py": "4c8f29acc2598935a30b027691ccf3f01ae34b6b80a27aac0f7cf5349fafa28b",
    "cosmic_bridge_v01/code/pm.py": "1b2607ac0b7ba35c5df9ceb7375910ba75b7bed65a1edec79c0583f067a5d982",
}
EXPECTED_MATRIX = (
    ("ref_matched64", "reference", 64, 64, True, (0, 0, 0)),
    ("ref_matched128", "reference", 128, 128, True, (0, 0, 0)),
)


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def serializable(value):
    if isinstance(value, dict):
        return {str(key): serializable(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [serializable(item) for item in value]
    if isinstance(value, np.ndarray):
        return serializable(value.tolist())
    if isinstance(value, np.generic):
        return serializable(value.item())
    if isinstance(value, float) and not np.isfinite(value):
        return str(value)
    return value


def archive_snapshot():
    """Hash every regular archived file; no result or __pycache__ exclusions."""
    snapshot = {}
    for archive in ARCHIVES:
        for path in sorted(archive.rglob("*")):
            if path.is_file():
                stat = path.stat()
                snapshot[str(path.relative_to(EXPERIMENTS))] = {
                    "bytes": stat.st_size, "mtime_ns": stat.st_mtime_ns, "sha256": sha(path)}
    return snapshot


def snapshot_difference(before, after):
    return {
        "created": sorted(after.keys() - before.keys()),
        "removed": sorted(before.keys() - after.keys()),
        "modified": sorted(key for key in before.keys() & after.keys() if before[key] != after[key]),
    }


def load_runner():
    path = ROOT / "code/run_matched_pm.py"
    spec = importlib.util.spec_from_file_location("matched_pm_validation_target", path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


class Checks:
    def __init__(self):
        self.rows = []

    def require(self, name, passed, **detail):
        self.rows.append({"name": name, "scope": "implementation_consistency",
                          "passed": bool(passed), **detail})

    def upper(self, name, value, threshold):
        self.require(name, np.isfinite(value) and value <= threshold,
                     value=float(value), inclusive_upper_threshold=float(threshold))


def relative_l2(actual, expected):
    expected_norm = np.linalg.norm(expected)
    return float(np.linalg.norm(actual - expected) / max(expected_norm, 1e-300))


def direct_compensated_initial(nparticle, delta_i, reference_background):
    """Independent closed-form mapping, with no call to archived IC helpers."""
    axis = (np.arange(nparticle) + 0.5) / nparticle
    original = np.stack(np.meshgrid(axis, axis, axis, indexing="ij"), axis=-1).reshape(-1, 3)
    q = original - 0.5
    radius = np.linalg.norm(q, axis=1)
    window = np.ones_like(radius)
    shell = (radius > 0.25) & (radius < 0.4)
    s = (radius[shell] - 0.25) / 0.15
    window[shell] = 1 - 10 * s**3 + 15 * s**4 - 6 * s**5
    window[radius >= 0.4] = 0
    scale = (1 + delta_i * window)**(-1 / 3)
    a_initial = 0.02
    expansion = float(reference_background.evaluate_a(a_initial)["E"])
    growth_rate = float(reference_background.f_reference)
    velocity_scale = -growth_rate * delta_i * window / (3 * (1 + delta_i * window)**(4 / 3))
    return {
        "positions": np.mod(0.5 + scale[:, None] * q, 1),
        "momenta": a_initial**2 * expansion * velocity_scale[:, None] * q,
        "q": q, "q_radius": radius,
        "core_mask": (radius >= 0.25 * 0.25) & (radius <= 0.75 * 0.25),
    }


def manual_eds_kdk(pm, x, p, a_left, a_right, nmesh):
    """Independent EdS factors plus old full-FFT force for one KDK step."""
    mesh = pm.legacy.PMGrid(nmesh)
    first_force, _ = mesh.force(x)
    middle = np.sqrt(a_left * a_right)
    first_kick = 2 * (np.sqrt(middle) - np.sqrt(a_left))
    second_kick = 2 * (np.sqrt(a_right) - np.sqrt(middle))
    drift = 2 * (a_left**(-0.5) - a_right**(-0.5))
    halfway = p + 1.5 * first_kick * first_force
    following_x = np.mod(x + drift * halfway, 1)
    second_force, _ = mesh.force(following_x)
    following_p = halfway + 1.5 * second_kick * second_force
    return following_x, following_p


def source_structure(path):
    tree = ast.parse(Path(path).read_text())
    return tree, {node.name: node for node in tree.body
                  if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))}


def report(result):
    lines = ["# 匹配粒子与力网格运行器：小规模一致性验证", "",
             f"状态：**{result['status']}，{result['checks_passed']}/{result['checks_total']} 项通过**。", "",
             "本验证只使用 8³ 粒子的短步测试及临时文件，不运行 64³／128³ 生产轨迹，也不重复巨大 DFT 或占用生产初力预检。实现一致性通过不代表连续物理或坍缩事件已经收敛。", "",
             f"冻结协议 SHA256：`{result['protocol_sha256']}`。", "",
             "| 检查 | 状态 | 数值或说明 |", "| --- | --- | --- |"]
    for row in result["checks"]:
        value = row.get("value", row.get("reason", "见 JSON 详情"))
        lines.append(f"| {row['name']} | {'PASS' if row['passed'] else 'FAIL'} | {value} |")
    lines += ["", "旧实验目录在导入运行器之前及测试结束后逐文件校验；禁止隐式生成旧目录字节码。分片测试只在临时目录进行，检查顺序、形状、float64、范围、逐片哈希及完整往返。旧实验的物理离散失败不因新运行器一致性通过而改变。", ""]
    if result["failures"]:
        lines += ["失败保持原样：", ""]
        lines += [f"- {row['name']}：{json.dumps(serializable(row), ensure_ascii=False)}" for row in result["failures"]]
    return "\n".join(lines) + "\n"


def run_checks(runner, checks, detail):
    pm = runner.pm
    checks.require("frozen_protocol_matches_runner", runner.FROZEN_PROTOCOL == sha(ROOT / "protocol.md"))
    checks.require("exact_two_case_matched_matrix", runner.MATRIX == EXPECTED_MATRIX,
                   actual_matrix=runner.MATRIX, expected_matrix=EXPECTED_MATRIX)
    checks.require("archived_PM_import_identity", Path(pm.__file__).resolve() == ARCHIVES[0] / "code/pm_refined.py")
    checks.require("archived_sphere_import_identity", Path(runner.sphere.__file__).resolve() == ARCHIVES[0] / "code/spherical_collapse.py")
    checks.require("old_canonical_KDK_identity", pm.advance_kdk is pm.legacy.advance_kdk)
    for relative_path, expected in runner.SOURCE_SHA256.items():
        checks.require(f"runner_dependency:{relative_path}", sha(EXPERIMENTS / relative_path) == expected)
    constants = {"EVENT_THRESHOLD": 0.02, "RADIUS_THRESHOLD": 0.03,
                 "PAIR_EVENT_THRESHOLD": 0.01, "SHAPE_THRESHOLD": 0.1,
                 "MAX_TOTAL_WALL_SECONDS": 14400, "MAX_STEPS": 10000,
                 "MAX_ARCHIVE_BYTES": 100_000_000}
    for name, value in constants.items():
        checks.require(f"frozen_constant:{name}", getattr(runner, name) == value,
                       actual=getattr(runner, name), expected=value)

    tree, functions = source_structure(ROOT / "code/run_matched_pm.py")
    check_calls = {}
    for node in ast.walk(functions["run_case"]):
        if (isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
                and node.func.id == "check" and node.args and isinstance(node.args[0], ast.Constant)):
            if len(node.args) > 2:
                threshold = node.args[2]
                check_calls[node.args[0].value] = (threshold.value if isinstance(threshold, ast.Constant)
                                                   else getattr(runner, threshold.id))
    expected_thresholds = {"mass": 1e-12, "total_force": 1e-11, "actual_drift_cells": 0.2,
                           "actual_dt_over_tdyn": 0.05, "nonhomology": 0.1, "axis_deviation": 0.1,
                           "PM_vs_ODE_event": 0.02, "PM_vs_ODE_radius": 0.03}
    checks.require("both_cases_apply_frozen_case_thresholds",
                   all(check_calls.get(name) == value for name, value in expected_thresholds.items()),
                   actual=check_calls, expected=expected_thresholds)
    checks.require("no_name_condition_around_ODE_thresholds",
                   not any(isinstance(node, ast.If)
                           and any(isinstance(child, ast.Name) and child.id == "name" for child in ast.walk(node.test))
                           and any(isinstance(child, ast.Constant) and child.value in ("PM_vs_ODE_event", "PM_vs_ODE_radius")
                                   for child in ast.walk(node)) for node in ast.walk(functions["run_case"])))
    checks.require("actual_tdyn_uses_minimum_of_both_endpoints",
                   any(isinstance(node, ast.Call) and isinstance(node.func, ast.Name) and node.func.id == "min"
                       and [item.id for item in node.args if isinstance(item, ast.Name)] == ["tdyn_left", "tdyn_right"]
                       for node in ast.walk(functions["run_case"])))

    reference = runner.sphere.load_backgrounds()["reference"]
    calibration = json.loads((ARCHIVES[0] / "results/sphere_summary.json").read_text())
    delta_i = float(calibration["delta_i"])
    checks.require("unchanged_calibrated_initial_amplitude", delta_i == 0.07295182597218097,
                   value=delta_i)
    initial = pm.compensated_initial(8, delta_i, reference, force_mesh=8)
    direct = direct_compensated_initial(8, delta_i, reference)
    for field in ("positions", "momenta", "q", "q_radius"):
        checks.upper(f"independent_initial_mapping:{field}",
                     float(np.max(np.abs(initial[field] - direct[field]))), 2e-14)
    checks.require("independent_fixed_core_labels", np.array_equal(initial["core_mask"], direct["core_mask"]))
    second_mesh_initial = pm.compensated_initial(8, delta_i, reference, force_mesh=16)
    for field in ("positions", "momenta", "q", "q_radius", "core_mask"):
        checks.require(f"initial_physics_independent_of_force_mesh:{field}",
                       np.array_equal(initial[field], second_mesh_initial[field]))
    detail["initial_conditions"] = {"nparticle": 8, "delta_i": delta_i,
                                     "reference_background_sha256": sha(reference.path),
                                     "calibration_sha256": sha(ARCHIVES[0] / "results/sphere_summary.json"),
                                     "production_particle_grids_constructed": False}

    eds = pm.Background.eds()
    mesh = pm.PMGrid(8, workers=1)
    x, p = initial["positions"].copy(), initial["momenta"].copy()
    reference_x, reference_p = x.copy(), p.copy()
    step_errors = []
    for left, right in zip((0.02, 0.0201, 0.0202), (0.0201, 0.0202, 0.0203)):
        x, p, _, _ = pm.advance_kdk(x, p, left, right, eds, mesh)
        reference_x, reference_p = manual_eds_kdk(pm, reference_x, reference_p, left, right, 8)
        position_error = float(np.max(np.abs(np.mod(x - reference_x + 0.5, 1) - 0.5)))
        momentum_error = relative_l2(p, reference_p)
        checks.upper(f"independent_EdS_KDK_positions:a={right}", position_error, 2e-13)
        checks.upper(f"independent_EdS_KDK_momenta:a={right}", momentum_error, 1e-10)
        step_errors.append({"a": right, "position_max_abs": position_error, "momentum_relative_l2": momentum_error})
    detail["independent_short_KDK"] = {"nparticle": 8, "nforce": 8, "accepted_steps": 3,
                                       "reference": "Closed-form EdS factors plus archived full-complex FFT force",
                                       "steps": step_errors}

    # Three artificial states exercise max Delta N, drift, and dynamical-time
    # bounds separately. No physical inference is made from these states.
    zero = np.zeros_like(initial["positions"])
    policies = [("ordinary", initial["positions"], initial["momenta"], zero),
                ("large_momentum", initial["positions"], initial["momenta"] + 0.2 * initial["q"], zero),
                ("short_dynamical_time", np.mod(0.5 + 2000**(-1/3) * initial["q"], 1), zero, zero)]
    policy_rows = []
    for label, x_test, p_test, force_test in policies:
        a_left = 0.02
        actual_a, diagnostics = pm.select_step(x_test, p_test, force_test, a_left, eds, initial, 8, 0.55, True)
        vectors = np.mod(x_test[initial["core_mask"]] - 0.5 + 0.5, 1) - 0.5
        y = float(np.median(np.linalg.norm(vectors, axis=1) / initial["q_radius"][initial["core_mask"]]))
        tdyn = np.sqrt(2) * a_left**1.5 * y**1.5
        delta_n = np.log(0.55 / 0.02) / 512
        halvings = 0
        while True:
            expected_a = a_left * np.exp(delta_n)
            middle = np.sqrt(a_left * expected_a)
            kick = 2 * (np.sqrt(middle) - np.sqrt(a_left))
            drift = 2 * (a_left**(-0.5) - expected_a**(-0.5))
            displacement = drift * (p_test + 1.5 * kick * force_test)
            cells = float(np.max(np.linalg.norm(displacement, axis=1)) * 8)
            delta_t = (2 / 3) * (expected_a**1.5 - a_left**1.5)
            if cells <= 0.075 and delta_t / tdyn <= 0.015:
                break
            delta_n /= 2
            halvings += 1
            if halvings > 48:
                raise RuntimeError("Independent small-state time-policy reference failed to terminate")
        checks.upper(f"halfstep_independent_endpoint:{label}", abs(actual_a / expected_a - 1), 1e-12)
        checks.require(f"halfstep_independent_halvings:{label}", diagnostics["step_halvings"] == halvings)
        checks.upper(f"halfstep_maximum_delta_ln_a:{label}", np.log(actual_a / a_left), np.log(0.55 / 0.02) / 512 * (1 + 1e-12))
        checks.upper(f"halfstep_drift_bound:{label}", cells, 0.075 * (1 + 1e-12))
        checks.upper(f"halfstep_dynamical_bound:{label}", delta_t / tdyn, 0.015 * (1 + 1e-12))
        policy_rows.append({"state": label, "a_right": actual_a, "step_halvings": halvings,
                            "drift_cells": cells, "dt_over_tdyn": delta_t / tdyn})
    checks.require("halfstep_test_exercises_reduction", any(row["step_halvings"] > 0 for row in policy_rows))
    detail["halfstep_policy"] = policy_rows

    # Zero initial force yields an exactly solvable EdS drift event, while the
    # final complete partial KDK may still compute a nonzero end force.
    a_left, a_right = 0.02, 0.0205
    scale_left, scale_right, scale_target = 199**(-1/3), 201**(-1/3), 200**(-1/3)
    full_drift = 2 * (a_left**(-0.5) - a_right**(-0.5))
    velocity_coefficient = (scale_right - scale_left) / full_drift
    x_event = np.mod(0.5 + scale_left * initial["q"], 1)
    p_event = velocity_coefficient * initial["q"]
    event_a, calls = runner.locate_event_a(x_event, p_event, zero, a_left, a_right, eds, initial)
    target_drift = (scale_target - scale_left) / velocity_coefficient
    expected_event_a = (a_left**(-0.5) - target_drift / 2)**(-2)
    checks.upper("partial_event_closed_form_scale_factor", abs(event_a / expected_event_a - 1), 1e-9)
    predicted, _ = pm.predict_positions(x_event, p_event, zero, a_left, event_a, eds)
    complete_x, _, _, _ = pm.advance_kdk(x_event, p_event, a_left, event_a, eds, mesh, zero)
    checks.upper("partial_event_prediction_matches_complete_KDK_position", float(np.max(np.abs(predicted - complete_x))), 2e-14)
    checks.upper("partial_event_reaches_density_proxy_200", abs(pm.core_scale(complete_x, initial)**(-3) / 200 - 1), 1e-8)
    checks.require("partial_event_reports_root_calls", isinstance(calls, int) and calls >= 2)
    try:
        runner.locate_event_a(initial["positions"], zero, zero, a_left, a_right, eds, initial)
        rejected = False
    except ValueError:
        rejected = True
    checks.require("unbracketed_event_is_not_fabricated", rejected)
    detail["event_helper"] = {"actual_a": event_a, "closed_form_a": expected_event_a,
                               "root_calls": calls, "scope": "Artificial near-threshold drift state, not a collapse simulation"}

    arrays = runner.corrected_final_profiles(initial["positions"], initial["momenta"], 0.02, eds, initial, 8)
    old_arrays = pm.final_profiles(initial["positions"], initial["momenta"], 0.02, eds, initial, 8)
    for axis in ("projection_axis_x_box", "projection_axis_y_box"):
        checks.require(f"projection_nodal_coordinate:{axis}", np.array_equal(arrays[axis], np.arange(8) / 8))
    unchanged = [key for key in old_arrays if not key.startswith("projection_axis_")]
    checks.require("coordinate_correction_preserves_all_density_and_profile_values",
                   all(np.array_equal(arrays[key], old_arrays[key]) for key in unchanged))
    checks.upper("projected_mass_normalization", abs(float(arrays["projection_density_over_mean"].mean()) - 1), 1e-12)

    # Scalar outcomes must retain missingness rather than silently assigning 0.
    for missing_name in ("ref_matched64", "ref_matched128"):
        runs = {name: {"status": "completed", "event200": {"a": 0.5}} for name in ("ref_matched64", "ref_matched128")}
        runs[missing_name]["event200"] = None
        comparison, flags = runner.pair_comparison(runs)
        checks.require(f"missing_event_is_null_and_failed:{missing_name}",
                       comparison["matched128_event_relative_to_matched64"] is None
                       and not flags[0]["passed"] and flags[0]["value"] is None)
    comparison, flags = runner.pair_comparison({})
    checks.require("missing_cases_explicitly_fail", comparison["matched128_event_relative_to_matched64"] is None
                   and all(not flag["passed"] for flag in flags))
    for event_a_test, should_pass in ((0.503, True), (0.506, False)):
        _, flags = runner.pair_comparison({"ref_matched64": {"status": "completed", "event200": {"a": 0.5}},
                                           "ref_matched128": {"status": "completed", "event200": {"a": event_a_test}}})
        checks.require(f"pair_event_threshold:a={event_a_test}", flags[0]["passed"] is should_pass)
    checks.require("missing_scalar_check_fails", not runner.check("required_event", None, 0.02)["passed"])

    with tempfile.TemporaryDirectory(prefix="matched_pm_validation_") as temporary:
        temporary = Path(temporary)
        positions = initial["positions"].astype(np.float64, copy=True)
        momenta = initial["momenta"].astype(np.float64, copy=True)
        manifest = runner.save_final_state(temporary, "tiny8", positions, momenta, initial["center"], 0.02)
        checks.require("separate_position_momentum_archives", set(manifest) == {"positions", "momenta"}
                       and manifest["positions"]["file"] != manifest["momenta"]["file"])
        for quantity, expected_array in (("positions", positions), ("momenta", momenta)):
            entry = manifest[quantity]
            path = temporary / entry["file"]
            checks.require(f"state_archive_hash:{quantity}", entry["sha256"] == sha(path))
            checks.require(f"state_archive_bytes:{quantity}", entry["bytes"] == path.stat().st_size
                           and entry["bytes"] < 100_000_000)
            checks.require(f"state_archive_shape_dtype_manifest:{quantity}", entry["shape"] == [512, 3]
                           and entry["dtype"] == "float64")
            with np.load(path, allow_pickle=False) as restored:
                checks.require(f"state_archive_exact_roundtrip:{quantity}",
                               restored[quantity].dtype == np.float64 and np.array_equal(restored[quantity], expected_array))
                checks.require(f"state_archive_event_metadata:{quantity}",
                               np.array_equal(restored["center_box"], initial["center"]) and float(restored["a"]) == 0.02)
        preserved_hashes = {path.name: sha(path) for path in temporary.iterdir()}
        try:
            runner.save_final_state(temporary, "tiny8", positions, momenta, initial["center"], 0.02)
            overwrite_rejected = False
        except RuntimeError:
            overwrite_rejected = True
        checks.require("state_archive_overwrite_rejected_without_changes", overwrite_rejected
                       and preserved_hashes == {path.name: sha(path) for path in temporary.iterdir()})
        try:
            runner.save_final_state(temporary, "wrong_dtype", positions.astype(np.float32), momenta, initial["center"], 0.02)
            dtype_rejected = False
        except ValueError:
            dtype_rejected = True
        checks.require("non_float64_final_state_rejected", dtype_rejected)
        invalid = positions.copy()
        invalid[0, 0] = np.nan
        try:
            runner.save_final_state(temporary, "nonfinite", invalid, momenta, initial["center"], 0.02)
            nonfinite_rejected = False
        except ValueError:
            nonfinite_rejected = True
        checks.require("nonfinite_final_state_rejected", nonfinite_rejected)
        runner.write_json(temporary / "missing.json", {"event200": None, "comparison": None})
        missing = json.loads((temporary / "missing.json").read_text())
        checks.require("missing_events_remain_JSON_null", missing == {"event200": None, "comparison": None})
        detail["terminal_archive"] = {"temporary_only": True, "particle_count": 512, "manifest": manifest,
                                      "storage": "One float64 positions file and one float64 momenta file"}

    vectors = pm.minimum_vector(initial["positions"], initial["center"])
    shape = runner.force_shape_diagnostics(initial["positions"], -vectors, initial)
    core_vector = vectors[initial["core_mask"]]
    direction = core_vector / np.linalg.norm(core_vector, axis=1)[:, None]
    expected_q4 = float(np.mean(np.sum(direction**4, axis=1)) - 3 / 5)
    checks.upper("angular_Q4_independent_definition", abs(shape["angular_cubic_Q4"] - expected_q4), 1e-14)
    checks.upper("radial_force_has_zero_transverse_component", shape["core_transverse_to_radial_force_L2_ratio"], 1e-13)
    zero_shape = runner.force_shape_diagnostics(initial["positions"], zero, initial)
    checks.require("zero_radial_force_ratio_is_missing_not_NaN_or_zero", zero_shape["core_transverse_to_radial_force_L2_ratio"] is None)
    detail["scope"] = {"production_cases_run": 0, "production_initial_force_calls": 0,
                        "largest_force_mesh": 8, "largest_particle_grid": 8,
                        "large_DFT_repeated": False, "scientific_acceptance_inferred": False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--protocol-sha256", required=True)
    args = parser.parse_args()
    if len(args.protocol_sha256) != 64 or any(char not in "0123456789abcdef" for char in args.protocol_sha256):
        parser.error("A lowercase SHA256 of the frozen protocol is required")
    protocol_hash = sha(ROOT / "protocol.md")
    if protocol_hash != args.protocol_sha256:
        raise RuntimeError("Protocol differs from the explicitly supplied frozen hash")
    started = time.perf_counter()
    checks, detail = Checks(), {}
    archive_before = archive_snapshot()
    runner = None
    try:
        for relative_path, expected in EXPECTED_SOURCE_HASHES.items():
            actual = sha(EXPERIMENTS / relative_path)
            checks.require(f"archived_source_hash:{relative_path}", actual == expected,
                           expected_sha256=expected, actual_sha256=actual)
        if any(not row["passed"] for row in checks.rows):
            raise RuntimeError("Archived source mismatch: cannot validate against changed physical code")
        runner = load_runner()
        run_checks(runner, checks, detail)
    except Exception:
        checks.require("execution_exception", False, traceback=traceback.format_exc())
    archive_after = archive_snapshot()
    differences = snapshot_difference(archive_before, archive_after)
    checks.require("archived_directories_unchanged", not any(differences.values()), **differences)
    detail["archive_audit"] = {"files_checked_before": len(archive_before),
                               "files_checked_after": len(archive_after),
                               "before": archive_before, "after": archive_after,
                               "differences": differences,
                               "bytecode_writes_disabled": sys.dont_write_bytecode}
    failures = [row for row in checks.rows if not row["passed"]]
    result_directory = ROOT / "results"
    report_directory = ROOT / "reports"
    result_directory.mkdir(parents=True, exist_ok=True)
    report_directory.mkdir(parents=True, exist_ok=True)
    attempts = sorted(result_directory.glob("matched_runner_validation_attempt_*.json"))
    attempt = max([int(path.stem.rsplit("_", 1)[1]) for path in attempts], default=0) + 1
    result = {"status": "FAIL" if failures else "PASS", "passed": not failures, "attempt": attempt,
              "created_utc": datetime.now(timezone.utc).isoformat(),
              "protocol_sha256": protocol_hash,
              "runner_sha256": sha(ROOT / "code/run_matched_pm.py"),
              "checks_passed": len(checks.rows) - len(failures), "checks_total": len(checks.rows),
              "checks": checks.rows, "failures": failures, "implementation_failures": failures,
              "detail": detail,
              "code_sha256": {"validate_matched_runner.py": sha(__file__),
                              "run_matched_pm.py": sha(ROOT / "code/run_matched_pm.py")},
              "previous_attempts": [{"path": path.name, "sha256": sha(path)} for path in attempts],
              "environment": {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__},
              "wall_seconds": time.perf_counter() - started,
              "scope": "Small runner/operator/archive consistency only; not production or physical convergence"}
    encoded = json.dumps(serializable(result), indent=2, ensure_ascii=False, allow_nan=False) + "\n"
    with (result_directory / f"matched_runner_validation_attempt_{attempt:03d}.json").open("x") as handle:
        handle.write(encoded)
    (result_directory / "matched_runner_validation.json").write_text(encoded)
    (report_directory / "matched_runner_validation_cn.md").write_text(report(result))
    print(json.dumps(serializable({"status": result["status"], "attempt": attempt,
                                   "checks_passed": result["checks_passed"], "checks_total": result["checks_total"],
                                   "wall_seconds": result["wall_seconds"], "failures": failures}),
                     ensure_ascii=False, indent=2), flush=True)
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
