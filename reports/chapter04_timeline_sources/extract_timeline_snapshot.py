#!/usr/bin/env python3
"""Extract a small, exact plotting snapshot from the fixed CMB archive.

This script selects saved arrays and checks provenance. It does not integrate
the lattice, optimize a sky, or assign physical ages to model steps.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess

import numpy as np


HERE = Path(__file__).resolve().parent
COMMIT = "f42db51748003fc83be064dade9bb52de5eb8319"
REPOSITORY = "https://github.com/maris205/riemann_cmb"
SOURCE_FILES = (
    "results/lattice_timeline.npz", "results/lattice_timeline.json",
    "reports/timeline_verification.json", "code/run_lattice_timeline.py",
    "code/verify_lattice_timeline.py", "code/lattice_transfer.py",
)
BRANCHES = ("p2_nonlinear", "p2_linear", "p0_nonlinear", "p2_nonlinear_undamped")


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def committed_digest(source, relative):
    """Stream the fixed Git blob; do not stage or change the source checkout."""
    command = ["git", "-C", str(source), "cat-file", "blob", f"{COMMIT}:{relative}"]
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    h = hashlib.sha256()
    size = 0
    while block := process.stdout.read(1024 * 1024):
        h.update(block)
        size += len(block)
    error = process.stderr.read().decode(errors="replace")
    if process.wait() != 0:
        raise RuntimeError(f"Cannot read fixed source blob {relative}: {error}")
    return h.hexdigest(), size


def require(condition, message):
    if not condition:
        raise ValueError(message)


def validate_arrays(data, summary):
    require(np.array_equal(data["steps"], np.arange(96)), "Unexpected step axis")
    require(np.array_equal(data["snapshot_steps"], [0, 45, 70, 95]), "Unexpected stored states")
    require(tuple(data["branch_names"]) == BRANCHES, "Unexpected branch order")
    require(data["slices_z48"].shape == (4, 96, 96), "Unexpected slice shape")
    require(data["histories"].shape == (4, 96, 5), "Unexpected history shape")
    require(np.isfinite(data["slices_z48"]).all(), "Nonfinite slice")
    require(np.isfinite(data["histories"]).all(), "Nonfinite history")
    require((data["histories"][:, :, 1] > 0).all(), "Nonpositive volume variance")
    for i, name in enumerate(BRANCHES):
        require(np.array_equal(data["histories"][0, :46], data["histories"][i, :46]),
                f"The prefix is not shared by {name}")
        require(data["histories"][i, -1, 1] == summary["branches"][name]["final_variance"],
                f"Final variance differs from source summary for {name}")
    require(np.isclose(summary["checks"]["late_phase_proxy_p0"],
                       summary["checks"]["late_phase_proxy_p2"], rtol=1e-14, atol=0),
            "Late constant-control propagation proxy differs")
    require(float(np.max(np.abs(data["slices_z48"]))) < 3.2,
            "The stated shared colour range would clip the shown slices")
    require(np.isnan(data["c2_p2"][0]) and np.isfinite(data["c2_p2"][1:]).all(),
            "Unexpected coefficient-array convention")


def extract(source):
    source = source.resolve()
    sources = {}
    for relative in SOURCE_FILES:
        path = source / relative
        file_hash = digest(path)
        git_hash, git_size = committed_digest(source, relative)
        require(file_hash == git_hash and path.stat().st_size == git_size,
                f"Working file does not match the fixed source commit: {relative}")
        sources[relative] = {"sha256": file_hash, "bytes": git_size,
                             "matches_fixed_commit": True,
                             "url": f"{REPOSITORY}/blob/{COMMIT}/{relative}"}

    summary = json.loads((source / "results/lattice_timeline.json").read_text())
    old_verification = json.loads((source / "reports/timeline_verification.json").read_text())
    require(old_verification["status"] == "PASS", "The archived source verification did not pass")
    for relative in ("results/lattice_timeline.npz", "results/lattice_timeline.json"):
        require(old_verification["inputs_sha256"][relative] == sources[relative]["sha256"],
                f"Archived verification refers to a different source: {relative}")

    with np.load(source / "results/lattice_timeline.npz", allow_pickle=False) as archive:
        volumes = archive["snapshots"]
        require(volumes.shape == (4, 96, 96, 96), "Unexpected original volume shape")
        data = {key: archive[key].copy() for key in
                ("steps", "snapshot_steps", "branch_names", "histories", "c2_p2")}
        data["slices_z48"] = volumes[:, :, :, 48].copy()
        volume_variances = np.var(volumes, axis=(1, 2, 3))
        require(np.allclose(volume_variances,
                            data["histories"][0, data["snapshot_steps"], 1],
                            rtol=1e-14, atol=1e-14),
                "Saved full-volume variances differ from history")
    validate_arrays(data, summary)

    metadata = {
        "source_repository": REPOSITORY, "source_commit": COMMIT,
        "source_npz_sha256": sources["results/lattice_timeline.npz"]["sha256"],
        "source_json_sha256": sources["results/lattice_timeline.json"]["sha256"],
        "selection": "snapshots[:, :, :, 48] copied exactly; no interpolation or normalization",
        "volume_shape": [96, 96, 96], "slice_axis": 2, "slice_index": 48,
        "snapshot_steps": [0, 45, 70, 95], "shared_colour_limits": [-3.2, 3.2],
        "units": {"q": "original dimensionless scalar-field units",
                  "variance": "population variance of the entire 96^3 field in original q^2 units",
                  "n": "dimensionless model update index; no cosmic-age calibration"},
        "moment_columns": summary["moment_columns"],
        "phase1": summary["phase1"], "phase2": summary["phase2"],
        "old_verification_status": "PASS", "old_verification_checks": [52, 52],
        "simulation_or_fit_executed_during_extraction": False,
        "not_included": ["complete 3D fields", "observed or fitted sky pixels",
                         "step-to-age mapping", "matter density or galaxy identification"],
    }
    data["metadata_json"] = np.array(json.dumps(metadata, ensure_ascii=False, sort_keys=True))
    np.savez_compressed(HERE / "timeline_snapshot.npz", **data)
    shutil.copyfile(source / "results/lattice_timeline.json", HERE / "lattice_timeline_source.json")
    shutil.copyfile(source / "reports/timeline_verification.json", HERE / "timeline_verification_source.json")

    output_names = ("timeline_snapshot.npz", "lattice_timeline_source.json",
                    "timeline_verification_source.json", "extract_timeline_snapshot.py", "README.md")
    manifest = {
        "format": "rsm-chapter04-timeline-plot-snapshot-v1", "source_repository": REPOSITORY,
        "source_commit": COMMIT, "source_inputs": sources,
        "extraction": metadata,
        "checks": {
            "source_bytes_match_fixed_commit": True,
            "source_verification_is_bound_to_input_hashes": True,
            "snapshot_selection_is_exact": True,
            "full_volume_variance_matches_saved_history": True,
            "four_branches_share_complete_prefix": True,
            "summary_final_variances_match_histories": True,
            "common_colour_scale_does_not_clip_saved_slices": True,
        },
        "local_files": {name: {"sha256": digest(HERE / name),
                                "bytes": (HERE / name).stat().st_size} for name in output_names},
        "scope": "Plotting snapshot and byte-identical old summaries; no new simulation or fit.",
    }
    (HERE / "manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n")
    verify_local()


def verify_local():
    manifest = json.loads((HERE / "manifest.json").read_text())
    require(manifest["source_commit"] == COMMIT, "Unexpected manifest source commit")
    for name, record in manifest["local_files"].items():
        path = HERE / name
        require(path.stat().st_size == record["bytes"] and digest(path) == record["sha256"],
                f"Local snapshot integrity failure: {name}")
    summary = json.loads((HERE / "lattice_timeline_source.json").read_text())
    verification = json.loads((HERE / "timeline_verification_source.json").read_text())
    require(digest(HERE / "lattice_timeline_source.json") ==
            manifest["source_inputs"]["results/lattice_timeline.json"]["sha256"],
            "Copied source summary differs")
    with np.load(HERE / "timeline_snapshot.npz", allow_pickle=False) as archive:
        data = {key: archive[key].copy() for key in archive.files}
    validate_arrays(data, summary)
    metadata = json.loads(str(data["metadata_json"]))
    require(metadata == manifest["extraction"], "Embedded provenance differs from manifest")
    require(metadata["source_npz_sha256"] ==
            verification["inputs_sha256"]["results/lattice_timeline.npz"],
            "Embedded source hash differs from old verification")
    print("PASS: 5 local file hashes, exact snapshot metadata, stored histories and source binding.")
    print("No simulation, fitting, network access or neighboring-repository access during --verify-local.")
    print(f"Plot snapshot: {(HERE / 'timeline_snapshot.npz').stat().st_size:,} bytes")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    choice = parser.add_mutually_exclusive_group(required=True)
    choice.add_argument("--source-root", type=Path, help="Local checkout containing the fixed CMB source commit")
    choice.add_argument("--verify-local", action="store_true", help="Check this snapshot without the CMB checkout")
    args = parser.parse_args()
    extract(args.source_root) if args.source_root else verify_local()
