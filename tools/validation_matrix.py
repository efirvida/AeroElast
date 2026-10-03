#!/usr/bin/env python3
"""Query and maintain the AeroElast validation matrix store.

Usage:
    python tools/validation_matrix.py check [--store DIR] [--group ID] [--json]

The store lives in `docs/validation/`:

    references.yaml             canonical bibliography (rector)
    groups.yaml                 group registry + un_inventoried allowlist
    flags.yaml                  flag registry (judgement flags + derived flags)
    rows/<group-slug>.yaml      the rows of one group

`docs/validation-matrix.md` and `docs/references.md` are generated from the store;
they are the paper-facing views, never the input.

Model. A row is one test group: it owns an optional single `result` (our own
measured value) and one or more `comparisons`. A comparison pairs that result with
one independent reference, its tolerance and its margin, because a single test
routinely measures against several references at once -- CCX and an analytical bar,
the paper's N=16 cell and the Kirchhoff closed form, several papers in a modal
table. `asserted` says whether an assertion consumes the comparison; it is
orthogonal to `measured.status`, which says whether a usable residual was printed.

This module currently implements the store model and `check`. The read verbs
(`find`, `get`, `list`, `headline`), the generators (`references`, `render`) and the
migration pipeline (`extract`, `diff-against-md`) are added by the pilot work units
in `odd/tasks/validation-matrix-store.md` section 13.1.

Validation is hand-rolled on purpose: `jsonschema` is not importable in the pinned
environment, and CONTRIBUTING rule 6 makes this tool part of the suite's discipline
rather than an optional extra. The JSON Schemas under `schemas/` document the same
rules for editors and reviewers; `check` is the enforcement.

Exit codes: 0 clean, 1 findings, 2 usage or I/O error.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import yaml

REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_STORE = REPO_ROOT / "docs" / "validation"

EXIT_OK = 0
EXIT_FINDINGS = 1
EXIT_ERROR = 2

ID_RE = re.compile(r"^[a-z0-9]+(\.[a-z0-9_]+)+$")
SLUG_RE = re.compile(r"^[a-z][a-z0-9_]*$")
NODE_ID_RE = re.compile(r"^tests/[A-Za-z0-9_./-]+\.py::[A-Za-z0-9_]+(\[[^\]]*\])?$")
DATE_RE = re.compile(r"^\d{4}-\d{2}-\d{2}$")

REFERENCE_KINDS = {"paper", "code", "analytical", "self", "schema"}
TOLERANCE_KINDS = {"rtol", "atol", "rel_err", "sign", "exact", "subset"}
RATIO_TOLERANCE_KINDS = {"rtol", "rel_err"}
MEASURED_STATUSES = {"measured", "not_printed", "not_measured", "sign_only", "n/a"}
SUITE_TOLERANCE_RULE = 0.05

ROW_KEYS = {
    "id",
    "group",
    "title",
    "tests",
    "evidence",
    "validates",
    "result",
    "comparisons",
    "flags",
    "notes",
    "history",
}
COMPARISON_KEYS = {
    "label",
    "asserted",
    "reference",
    "tolerance",
    "measured",
    "expected",
    "notes",
}
RESULT_KEYS = {"text", "raw", "unit", "node"}
REFERENCE_KEYS = {"kind", "label", "citation", "also_asserts"}
TOLERANCE_KEYS = {"kind", "value", "source", "justified", "justification"}
MEASURED_KEYS = {"status", "run", "date", "raw", "margin_pct", "text"}
HISTORY_KEYS = {"rev", "note", "evidence"}
GROUP_KEYS = {
    "id",
    "title",
    "slug",
    "source_files",
    "common_tolerance",
    "headline",
    "provenance_note",
    "prose_after_table",
}
FLAG_KEYS = {"id", "derived", "label", "user_facing", "section", "legend"}


class StoreError(Exception):
    """The store could not be read at all."""


@dataclass
class Finding:
    level: str  # "error" or "warning"
    where: str
    message: str

    def as_text(self) -> str:
        return f"{self.level}: {self.where}: {self.message}"

    def as_dict(self) -> dict[str, str]:
        return {"level": self.level, "where": self.where, "message": self.message}


@dataclass
class Store:
    """The loaded store, plus the findings produced while loading it."""

    root: Path
    groups: dict[str, dict[str, Any]] = field(default_factory=dict)
    un_inventoried: dict[str, dict[str, Any]] = field(default_factory=dict)
    flag_entries: dict[str, dict[str, Any]] = field(default_factory=dict)
    derived_flags: set[str] = field(default_factory=set)
    rows: list[tuple[str, dict[str, Any]]] = field(default_factory=list)
    findings: list[Finding] = field(default_factory=list)

    def error(self, where: str, message: str) -> None:
        self.findings.append(Finding("error", where, message))

    def warn(self, where: str, message: str) -> None:
        self.findings.append(Finding("warning", where, message))


def load_yaml(path: Path) -> Any:
    try:
        with path.open("r", encoding="utf-8") as handle:
            return yaml.safe_load(handle)
    except FileNotFoundError as exc:
        raise StoreError(f"missing file: {path}") from exc
    except yaml.YAMLError as exc:
        raise StoreError(f"invalid YAML in {path}: {exc}") from exc


def _unknown_keys(data: dict[str, Any], allowed: set[str]) -> list[str]:
    return sorted(set(data) - allowed)


def load_flags(store: Store) -> None:
    path = store.root / "flags.yaml"
    data = load_yaml(path)
    if not isinstance(data, dict) or not isinstance(data.get("flags"), list):
        store.error("flags.yaml", "expected a mapping with a 'flags' list")
        return
    for index, entry in enumerate(data["flags"]):
        where = f"flags.yaml[{index}]"
        if not isinstance(entry, dict):
            store.error(where, "expected a mapping")
            continue
        unknown = _unknown_keys(entry, FLAG_KEYS)
        if unknown:
            store.error(where, f"unknown keys: {', '.join(unknown)}")
        flag_id = entry.get("id")
        if not isinstance(flag_id, str) or not SLUG_RE.match(flag_id):
            store.error(where, f"invalid flag id: {flag_id!r}")
            continue
        where = f"flags.yaml[{flag_id}]"
        if flag_id in store.flag_entries:
            store.error(where, "duplicate flag id")
            continue
        derived = entry.get("derived")
        if derived is not None:
            if not isinstance(derived, str) or not derived:
                store.error(where, "'derived' must be a non-empty expression string")
            if not isinstance(entry.get("user_facing"), bool):
                store.error(where, "a derived flag needs a boolean 'user_facing'")
            if entry.get("section") is not None:
                store.error(where, "a derived flag must not carry a 'section'")
            store.derived_flags.add(flag_id)
        else:
            if not isinstance(entry.get("section"), str):
                store.error(where, "a judgement flag needs a 'section'")
            if not isinstance(entry.get("label"), str):
                store.error(where, "a judgement flag needs a 'label'")
            if entry.get("user_facing") is not None:
                store.error(where, "'user_facing' belongs to derived flags only")
        store.flag_entries[flag_id] = entry


def load_groups(store: Store) -> None:
    path = store.root / "groups.yaml"
    data = load_yaml(path)
    if not isinstance(data, dict) or not isinstance(data.get("groups"), list):
        store.error("groups.yaml", "expected a mapping with a 'groups' list")
        return
    for index, entry in enumerate(data["groups"]):
        where = f"groups.yaml[{index}]"
        if not isinstance(entry, dict):
            store.error(where, "expected a mapping")
            continue
        unknown = _unknown_keys(entry, GROUP_KEYS)
        if unknown:
            store.error(where, f"unknown keys: {', '.join(unknown)}")
        group_id = entry.get("id")
        if not isinstance(group_id, str) or not group_id:
            store.error(where, f"invalid group id: {group_id!r}")
            continue
        where = f"groups.yaml[{group_id}]"
        if group_id in store.groups:
            store.error(where, "duplicate group id")
            continue
        slug = entry.get("slug")
        if not isinstance(slug, str) or not SLUG_RE.match(slug):
            store.error(where, f"invalid 'slug': {slug!r}")
        for source in entry.get("source_files") or []:
            if not isinstance(source, str) or not (REPO_ROOT / source).exists():
                store.error(where, f"source file does not exist: {source!r}")
        prose = entry.get("prose_after_table")
        if prose is not None and not (REPO_ROOT / str(prose).strip()).exists():
            store.warn(where, f"prose file not written yet: {str(prose).strip()}")
        store.groups[group_id] = entry

    for index, entry in enumerate(data.get("un_inventoried") or []):
        where = f"groups.yaml:un_inventoried[{index}]"
        if not isinstance(entry, dict) or not isinstance(entry.get("file"), str):
            store.error(where, "expected a mapping with a 'file'")
            continue
        if not (REPO_ROOT / entry["file"]).exists():
            store.error(where, f"file does not exist: {entry['file']}")
        if not isinstance(entry.get("reason"), str) or not entry["reason"]:
            store.error(where, "a documented drift entry needs a 'reason'")
        store.un_inventoried[entry["file"]] = entry


def _validate_reference(store: Store, where: str, ref: Any) -> None:
    if not isinstance(ref, dict):
        store.error(where, "'reference' must be a mapping")
        return
    unknown = _unknown_keys(ref, REFERENCE_KEYS)
    if unknown:
        store.error(where, f"reference: unknown keys: {', '.join(unknown)}")
    kind = ref.get("kind")
    if kind not in REFERENCE_KINDS:
        store.error(where, f"reference.kind must be one of {sorted(REFERENCE_KINDS)}")
    if not isinstance(ref.get("label"), str) or not ref["label"]:
        store.error(where, "reference.label must be a non-empty string")
    citation = ref.get("citation")
    if kind == "paper":
        if not isinstance(citation, str) or not SLUG_RE.match(citation):
            store.error(where, "reference.kind == 'paper' needs a 'citation' key")
    elif citation is not None:
        store.error(where, "reference.citation is only meaningful when kind == 'paper'")


def _validate_tolerance(store: Store, where: str, tol: Any) -> None:
    if not isinstance(tol, dict):
        store.error(where, "'tolerance' must be a mapping")
        return
    unknown = _unknown_keys(tol, TOLERANCE_KEYS)
    if unknown:
        store.error(where, f"tolerance: unknown keys: {', '.join(unknown)}")
    kind = tol.get("kind")
    if kind not in TOLERANCE_KINDS:
        store.error(where, f"tolerance.kind must be one of {sorted(TOLERANCE_KINDS)}")
    value = tol.get("value")
    if kind in {"sign", "exact", "subset"}:
        if value is not None:
            store.error(where, f"tolerance.value must be null for kind {kind!r}")
    elif not isinstance(value, (int, float)) or isinstance(value, bool) or value <= 0:
        store.error(where, f"tolerance.value must be a positive number for kind {kind!r}")
    if not isinstance(tol.get("source"), str) or not tol["source"]:
        store.error(where, "tolerance.source must state where the code sets it")
    if not isinstance(tol.get("justified"), bool):
        store.error(where, "tolerance.justified must be a boolean")
    if kind in RATIO_TOLERANCE_KINDS and isinstance(value, (int, float)):
        if value > SUITE_TOLERANCE_RULE:
            if not tol.get("justified"):
                store.error(
                    where,
                    f"tolerance {value} exceeds the {SUITE_TOLERANCE_RULE:.0%} suite rule "
                    "without justification",
                )
            justification = tol.get("justification")
            if not isinstance(justification, str) or not justification:
                store.error(where, "a tolerance above 5% needs a written 'justification'")


def _validate_measured(store: Store, where: str, measured: Any) -> None:
    if not isinstance(measured, dict):
        store.error(where, "'measured' must be a mapping")
        return
    unknown = _unknown_keys(measured, MEASURED_KEYS)
    if unknown:
        store.error(where, f"measured: unknown keys: {', '.join(unknown)}")
    status = measured.get("status")
    if status not in MEASURED_STATUSES:
        store.error(where, f"measured.status must be one of {sorted(MEASURED_STATUSES)}")
        return
    margin = measured.get("margin_pct")
    if status == "measured":
        if not isinstance(margin, (int, float)) or isinstance(margin, bool):
            store.error(where, "measured.status == 'measured' requires a numeric margin_pct")
    elif margin is not None:
        store.error(where, f"measured.margin_pct must be null when status == {status!r}")
    date = measured.get("date")
    if date is not None and not (isinstance(date, str) and DATE_RE.match(date)):
        store.error(where, f"measured.date must be YYYY-MM-DD, got {date!r}")
    raw = measured.get("raw")
    if raw is not None and (not isinstance(raw, (int, float)) or isinstance(raw, bool)):
        store.error(where, "measured.raw must be a number or null")


def _validate_result(store: Store, where: str, result: Any) -> None:
    if not isinstance(result, dict):
        store.error(where, "'result' must be a mapping")
        return
    unknown = _unknown_keys(result, RESULT_KEYS)
    if unknown:
        store.error(where, f"result: unknown keys: {', '.join(unknown)}")
    if not isinstance(result.get("text"), str) or not result["text"]:
        store.error(where, "result.text must be a non-empty string")
    raw = result.get("raw")
    if raw is not None and (not isinstance(raw, (int, float)) or isinstance(raw, bool)):
        store.error(where, "result.raw must be a number or null")


def _validate_comparison(store: Store, where: str, comparison: Any) -> bool:
    """Validate one comparison. Returns True when it is explicitly not asserted."""
    if not isinstance(comparison, dict):
        store.error(where, "expected a mapping")
        return False
    unknown = _unknown_keys(comparison, COMPARISON_KEYS)
    if unknown:
        store.error(where, f"unknown keys: {', '.join(unknown)}")
    if not isinstance(comparison.get("label"), str) or not comparison["label"]:
        store.error(
            where,
            "label must be a non-empty string naming which part of the test this "
            "comparison covers",
        )
    asserted = comparison.get("asserted")
    if not isinstance(asserted, bool):
        store.error(where, "'asserted' must be a boolean")
    _validate_reference(store, where, comparison.get("reference"))
    _validate_tolerance(store, where, comparison.get("tolerance"))
    _validate_measured(store, where, comparison.get("measured"))
    return isinstance(asserted, bool) and not asserted


def _validate_history(store: Store, where: str, history: Any) -> None:
    if history is None:
        return
    if not isinstance(history, list):
        store.error(where, "'history' must be a list")
        return
    for index, entry in enumerate(history):
        entry_where = f"{where}:history[{index}]"
        if not isinstance(entry, dict):
            store.error(entry_where, "expected a mapping")
            continue
        unknown = _unknown_keys(entry, HISTORY_KEYS)
        if unknown:
            store.error(entry_where, f"unknown keys: {', '.join(unknown)}")
        rev = entry.get("rev")
        if not isinstance(rev, str) or len(rev) < 4:
            store.error(entry_where, f"invalid 'rev': {rev!r}")
        if not isinstance(entry.get("note"), str) or not entry["note"]:
            store.error(entry_where, "a history entry needs a non-empty 'note'")


def validate_row(store: Store, where: str, row: Any) -> None:
    if not isinstance(row, dict):
        store.error(where, "expected a mapping")
        return
    unknown = _unknown_keys(row, ROW_KEYS)
    if unknown:
        store.error(where, f"unknown keys: {', '.join(unknown)}")

    row_id = row.get("id")
    if not isinstance(row_id, str) or not ID_RE.match(row_id):
        store.error(where, f"invalid id: {row_id!r}")
    else:
        where = f"{where}:{row_id}"

    group_id = row.get("group")
    group = store.groups.get(group_id) if isinstance(group_id, str) else None
    if group is None:
        store.error(where, f"unknown group: {group_id!r}")
    elif isinstance(row_id, str) and ID_RE.match(row_id):
        expected_prefix = group.get("slug")
        if expected_prefix and row_id.split(".")[0] != expected_prefix:
            store.error(
                where,
                f"id prefix {row_id.split('.')[0]!r} does not match the group slug "
                f"{expected_prefix!r}",
            )

    if not isinstance(row.get("title"), str) or not row["title"]:
        store.error(where, "title must be a non-empty string")
    if not isinstance(row.get("validates"), str) or not row["validates"]:
        store.error(where, "validates must be a non-empty string")

    evidence = row.get("evidence")
    if evidence not in (None, "out_of_band"):
        store.error(where, f"evidence must be null or 'out_of_band', got {evidence!r}")
    tests = row.get("tests")
    if not isinstance(tests, list):
        store.error(where, "'tests' must be a list of collected node ids")
    else:
        for node in tests:
            if not isinstance(node, str) or not NODE_ID_RE.match(node):
                store.error(where, f"not a collected node id: {node!r}")
            elif not (REPO_ROOT / node.split("::")[0]).exists():
                store.error(where, f"node id names a missing file: {node!r}")
        if not tests and evidence != "out_of_band":
            store.error(where, "an empty 'tests' list requires evidence: out_of_band")

    if row.get("result") is not None:
        _validate_result(store, where, row["result"])

    comparisons = row.get("comparisons")
    unasserted = False
    if not isinstance(comparisons, list) or not comparisons:
        store.error(
            where,
            "'comparisons' must be a non-empty list: one entry per independent reference",
        )
    else:
        for index, comparison in enumerate(comparisons):
            if _validate_comparison(store, f"{where}:comparisons[{index}]", comparison):
                unasserted = True

    _validate_history(store, where, row.get("history"))

    flags = row.get("flags")
    if not isinstance(flags, list):
        store.error(where, "'flags' must be a list")
        return
    for flag_id in flags:
        if flag_id in store.derived_flags:
            store.error(
                where,
                f"derived flag {flag_id!r} must not be stored: it is computed per "
                "comparison from tolerance and margin",
            )
        elif flag_id not in store.flag_entries:
            store.error(where, f"unknown flag id: {flag_id!r}")
    if unasserted and "info_only" not in flags:
        store.error(
            where,
            "a comparison with asserted: false requires the 'info_only' flag "
            "(check matrix invariant 8)",
        )


def load_rows(store: Store, only_group: str | None = None) -> None:
    rows_dir = store.root / "rows"
    if not rows_dir.exists():
        return
    for path in sorted(rows_dir.glob("*.yaml")):
        data = load_yaml(path)
        rel = path.relative_to(REPO_ROOT) if path.is_relative_to(REPO_ROOT) else path
        if isinstance(data, dict) and "rows" in data:
            file_group = data.get("group")
            unknown = _unknown_keys(data, {"group", "rows"})
            if unknown:
                store.error(str(rel), f"unknown keys: {', '.join(unknown)}")
            entries = data.get("rows")
            if not isinstance(entries, list):
                store.error(str(rel), "'rows' must be a list")
                continue
        elif isinstance(data, list):
            file_group = None
            entries = data
        else:
            store.error(str(rel), "expected a mapping with a 'rows' list")
            continue
        for index, row in enumerate(entries):
            where = f"{rel}[{index}]"
            if only_group is not None:
                if not isinstance(row, dict) or row.get("group") != only_group:
                    continue
            if file_group is not None and isinstance(row, dict) and row.get("group") != file_group:
                store.error(
                    where,
                    f"row group {row.get('group')!r} does not match the file group "
                    f"{file_group!r}",
                )
            validate_row(store, where, row)
            if isinstance(row, dict):
                store.rows.append((where, row))

    seen: dict[str, str] = {}
    for where, row in store.rows:
        row_id = row.get("id")
        if not isinstance(row_id, str):
            continue
        if row_id in seen:
            store.error(where, f"duplicate row id, first seen at {seen[row_id]}")
        else:
            seen[row_id] = where


def load_store(root: Path, only_group: str | None = None) -> Store:
    store = Store(root=root)
    if not root.exists():
        raise StoreError(f"store directory does not exist: {root}")
    load_flags(store)
    load_groups(store)
    if only_group is not None and only_group not in store.groups:
        raise StoreError(f"unknown group: {only_group}")
    load_rows(store, only_group=only_group)
    return store


def count_comparisons(store: Store) -> int:
    total = 0
    for _, row in store.rows:
        comparisons = row.get("comparisons")
        if isinstance(comparisons, list):
            total += len(comparisons)
    return total


def command_check(args: argparse.Namespace) -> int:
    store = load_store(args.store, only_group=args.group)
    errors = [finding for finding in store.findings if finding.level == "error"]
    warnings = [finding for finding in store.findings if finding.level == "warning"]
    if args.json:
        payload = {
            "store": str(store.root),
            "group": args.group,
            "rows": len(store.rows),
            "comparisons": count_comparisons(store),
            "errors": len(errors),
            "warnings": len(warnings),
            "findings": [finding.as_dict() for finding in store.findings],
        }
        print(json.dumps(payload, indent=2))
    else:
        for finding in store.findings:
            print(finding.as_text())
        scope = f" (group {args.group})" if args.group else ""
        print(
            f"check{scope}: {len(store.rows)} row(s), "
            f"{count_comparisons(store)} comparison(s), "
            f"{len(errors)} error(s), {len(warnings)} warning(s)"
        )
    return EXIT_FINDINGS if errors else EXIT_OK


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="validation_matrix",
        description="Query and maintain the validation matrix store.",
    )
    parser.add_argument(
        "--store",
        type=Path,
        default=DEFAULT_STORE,
        help="store directory (default: docs/validation)",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    check = subparsers.add_parser("check", help="validate the store")
    check.add_argument("--group", help="limit row validation to one group id")
    check.add_argument("--json", action="store_true", help="machine-readable output")
    check.set_defaults(func=command_check)

    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return args.func(args)
    except StoreError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return EXIT_ERROR


if __name__ == "__main__":
    raise SystemExit(main())
