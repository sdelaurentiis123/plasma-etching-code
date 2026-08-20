#!/usr/bin/env python3
"""Create/check the label-free Bosch per-wafer process summary."""
from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from petch.bosch_process_data import (  # noqa: E402
    load_bosch_process_traces,
    process_ingestion_manifest,
    summarize_bosch_process_traces,
)


DEFAULT_DATA = ROOT / "data" / "experimental" / "zenodo_17122442"
DEFAULT_OUTPUT = DEFAULT_DATA / "process_wafer_summary.csv"
DEFAULT_MANIFEST = DEFAULT_DATA / "process_wafer_summary_manifest.json"


def _render(process_path: Path, dictionary_path: Path):
    traces = load_bosch_process_traces(process_path, dictionary_path)
    summaries = summarize_bosch_process_traces(traces)
    metric_names = sorted(summaries[0].metrics)
    if any(set(summary.metrics) != set(metric_names) for summary in summaries):
        raise RuntimeError("Bosch process summaries do not share one schema")
    header = [
        "experiment_key", "source_group", "process_date", "wafer_number",
        *metric_names,
    ]
    rows = [[
        summary.experiment_key,
        summary.source_group,
        summary.process_date,
        summary.wafer_number,
        *(summary.metrics[name] for name in metric_names),
    ] for summary in summaries]

    from io import StringIO
    stream = StringIO(newline="")
    writer = csv.writer(stream, lineterminator="\n")
    writer.writerow(header)
    writer.writerows(rows)
    csv_text = stream.getvalue()
    manifest = process_ingestion_manifest()
    manifest.update({
        "summary_schema": header,
        "summary_row_count": len(rows),
        "summary_sha256": sha256(csv_text.encode("utf-8")).hexdigest(),
        "source_feature_schema_sizes": sorted({len(trace.channels) for trace in traces}),
        "calculated_without_measurement_csv": True,
    })
    manifest_text = json.dumps(manifest, indent=2, sort_keys=True) + "\n"
    return csv_text, manifest_text


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--process", type=Path, default=DEFAULT_DATA / "Process_data.nc")
    parser.add_argument(
        "--dictionary", type=Path, default=DEFAULT_DATA / "Dictionary_process.nc")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args(argv)

    csv_text, manifest_text = _render(args.process, args.dictionary)
    if args.check:
        if (
            not args.output.exists()
            or args.output.read_text() != csv_text
            or not args.manifest.exists()
            or args.manifest.read_text() != manifest_text
        ):
            raise SystemExit("committed Bosch process extraction is stale")
        print("Bosch process extraction is current")
        return 0

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(csv_text)
    args.manifest.write_text(manifest_text)
    print(f"wrote {len(csv_text.splitlines()) - 1} rows to {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
