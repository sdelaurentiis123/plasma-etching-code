#!/usr/bin/env python3
"""Create a compact checksum-tree inventory of untracked local result artifacts.

The index makes a large ignored campaign auditable without committing every checkpoint and
replicate. It does not move, delete, or reinterpret any result.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
from hashlib import sha256
import json
from pathlib import Path
import subprocess


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"


def _tracked_result_paths():
    payload = subprocess.check_output(
        ["git", "ls-files", "-z", "results"], cwd=ROOT)
    return {
        ROOT / value.decode("utf-8")
        for value in payload.split(b"\0") if value
    }


def _file_sha256(path):
    digest = sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _tree_digest(records):
    digest = sha256()
    for record in records:
        digest.update(record["path"].encode("utf-8"))
        digest.update(b"\0")
        digest.update(str(record["bytes"]).encode("ascii"))
        digest.update(b"\0")
        digest.update(record["sha256"].encode("ascii"))
        digest.update(b"\n")
    return digest.hexdigest()


def build_index():
    tracked = _tracked_result_paths()
    grouped = {}
    for path in sorted(RESULTS.rglob("*")):
        if (
            not path.is_file()
            or path in tracked
            or ".plot-cache" in path.parts
            or path.name == ".DS_Store"
            or path == RESULTS / "README.md"
            or (RESULTS / "curated") in path.parents
        ):
            continue
        relative = path.relative_to(ROOT).as_posix()
        group = path.relative_to(RESULTS).parts[0]
        grouped.setdefault(group, []).append({
            "path": relative,
            "bytes": path.stat().st_size,
            "sha256": _file_sha256(path),
        })

    group_rows = []
    global_records = []
    for name in sorted(grouped):
        records = grouped[name]
        global_records.extend(records)
        group_rows.append({
            "name": name,
            "file_count": len(records),
            "bytes": sum(item["bytes"] for item in records),
            "tree_sha256": _tree_digest(records),
        })
    global_records.sort(key=lambda item: item["path"])
    revision = subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    return {
        "schema": "petch.local-result-artifact-index.v1",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "git_revision_before_consolidation": revision,
        "scope": "untracked non-curated files below results/, excluding plot caches",
        "file_count": len(global_records),
        "bytes": sum(item["bytes"] for item in global_records),
        "tree_sha256": _tree_digest(global_records),
        "groups": group_rows,
        "note": (
            "This compact index authenticates local raw evidence by top-level checksum tree. "
            "It is not an experimental-validation claim and contains no simulation arrays."
        ),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output", type=Path,
        default=ROOT / "LOCAL_RESULT_ARTIFACT_INDEX_2026-07-16.json")
    args = parser.parse_args()
    output = args.output if args.output.is_absolute() else ROOT / args.output
    payload = build_index()
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(output.suffix + ".tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    temporary.replace(output)
    print(json.dumps({
        "output": str(output),
        "file_count": payload["file_count"],
        "bytes": payload["bytes"],
        "tree_sha256": payload["tree_sha256"],
    }, indent=2))


if __name__ == "__main__":
    main()
