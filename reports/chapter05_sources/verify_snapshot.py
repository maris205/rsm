"""Check the frozen Chapter 5 inputs locally, without running any fit."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    manifest = json.loads((ROOT / "manifest.json").read_text(encoding="utf-8"))
    entries = manifest["files"]
    if len(entries) != 11:
        raise ValueError("Expected exactly 11 frozen source files.")
    for entry in entries:
        path = (ROOT / entry["path"]).resolve()
        if not path.is_relative_to(ROOT):
            raise ValueError("Manifest path leaves the snapshot directory.")
        if path.stat().st_size != entry["bytes"]:
            raise ValueError(f"Byte count mismatch: {entry['path']}")
        if sha256(path) != entry["sha256"]:
            raise ValueError(f"SHA256 mismatch: {entry['path']}")

    alpha = json.loads(
        (ROOT / "alpha/raw/King2012_provenance.json").read_text(encoding="utf-8")
    )
    for entry in alpha["files"]:
        path = ROOT / "alpha/raw" / Path(entry["path"]).name
        if sha256(path) != entry["sha256"] or path.stat().st_size != entry["bytes"]:
            raise ValueError(f"King provenance mismatch: {path.name}")

    hubble = json.loads(
        (ROOT / "hubble/summary_provenance.json").read_text(encoding="utf-8")
    )
    for filename, field in [
        ("comparison.csv", "comparison_csv_sha256"),
        ("null_profile.json", "null_profile_sha256"),
    ]:
        if sha256(ROOT / "hubble" / filename) != hubble[field]:
            raise ValueError(f"Hubble provenance mismatch: {filename}")

    print(
        json.dumps(
            {
                "status": "passed",
                "source_files": len(entries),
                "source_bytes": sum(entry["bytes"] for entry in entries),
                "manifest": "reports/chapter05_sources/manifest.json",
                "fit_executed": False,
            },
            ensure_ascii=False,
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
