#!/usr/bin/env python3
"""Inspect generated Aero-FSI report schema manifests."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from aeroelast.solvers.aero import load_aero_report_schema, load_aero_report_table

_REPORT_TYPES = ("aero_report", "aero_spanwise_report", "aero_sectional_report")


def _render_summary(output_dir: Path) -> str:
    schema = load_aero_report_schema(output_dir)
    lines = [
        f"Manifest: {schema.manifest_path}",
        f"Backend: {schema.backend}",
        f"Blades: {schema.n_blades}",
        f"Sectional bins: {schema.sectional_bins}",
        "Reports:",
    ]
    for report in schema.reports.values():
        exists_label = "yes" if report.path.exists() else "no"
        lines.append(
            f"- {report.report_type}: {report.path.name}"
            f" [schema={report.schema_version}, exists={exists_label}]"
        )
        lines.append(f"  {report.granularity}")
    return "\n".join(lines)


def _render_report_summary(output_dir: Path, report_type: str) -> str:
    schema = load_aero_report_schema(output_dir)
    report = schema.get_report(report_type)
    lines = [
        f"Report: {report.report_type}",
        f"Path: {report.path}",
        f"Schema: {report.schema_version}",
        f"Granularity: {report.granularity}",
        f"Summary: {report.summary}",
    ]
    if report.path.exists():
        table = load_aero_report_table(output_dir, report_type)
        lines.append(f"Rows: {table.height}")
        lines.append(f"Columns: {len(table.columns)}")
    else:
        lines.append("Rows: unavailable (file not written yet)")
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="fem-shell-aero-report-info",
        description="Inspect the Aero-FSI report schema manifest in an output folder.",
    )
    parser.add_argument(
        "output_dir",
        type=Path,
        help="Output folder containing aero_report_schema.json",
    )
    parser.add_argument(
        "--json",
        action="store_true",
        help="Print the raw manifest JSON instead of a human-readable summary.",
    )
    parser.add_argument(
        "--report",
        choices=_REPORT_TYPES,
        help="Inspect a single report type instead of the whole manifest summary.",
    )
    parser.add_argument(
        "--fields",
        action="store_true",
        help="Print ordered field names for the selected report type.",
    )
    args = parser.parse_args(argv)

    if args.fields and args.report is None:
        print("Error: --fields requires --report.", file=sys.stderr)
        return 2

    try:
        if args.json:
            schema = load_aero_report_schema(args.output_dir)
            with schema.manifest_path.open(encoding="utf-8") as handle:
                payload = json.load(handle)
            if args.report is not None:
                payload = payload["reports"][args.report]
            print(json.dumps(payload, indent=2))
        elif args.report is not None and args.fields:
            schema = load_aero_report_schema(args.output_dir)
            print("\n".join(schema.get_report(args.report).fieldnames))
        elif args.report is not None:
            print(_render_report_summary(args.output_dir, args.report))
        else:
            print(_render_summary(args.output_dir))
    except FileNotFoundError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())