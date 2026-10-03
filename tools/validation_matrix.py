#!/usr/bin/env python3
"""Query and maintain the AeroElast validation matrix store.

Usage:
    python tools/validation_matrix.py check  [--group ID] [--json]
    python tools/validation_matrix.py find   SUBSTRING [--group ID] [--file F] [--json]
    python tools/validation_matrix.py list   [--group ID] [--file F] [--id ID]
                                             [--flag FLAG] [--near] [--gt5]
                                             [--unit row|comparison]
                                             [--sort id|slack] [--fields A,B] [--json]
    python tools/validation_matrix.py get    ID [--comparisons] [--json]
    python tools/validation_matrix.py headline [--json]
    python tools/validation_matrix.py set    ID PATH=VALUE... [--unset PATH]...
                                             [--add-flag FLAG]... [--remove-flag FLAG]...

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

Derived values (`slack_pp`, `near`, `gt5`) are computed per comparison and are never
stored: `list --unit comparison` shows them, and `set` refuses them.

Still to come, per the pilot work units in `odd/tasks/validation-matrix-store.md`
section 13.1: the `references` subgroup with its generator, and the migration
pipeline (`extract`, `diff-against-md`, `render`).

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
from typing import Any, Iterable

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
PATH_SEGMENT_RE = re.compile(r"^([A-Za-z_][A-Za-z0-9_]*)(?:\[(\d+)\])?$")

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
# Bibliography entries (`docs/validation/references.yaml`). This is a different shape
# from the nested `reference` of a comparison: one entry per work, with the auditing
# fields that make the DOI and the citation sites checkable.
BIB_KEYS = {
    "key",
    "section",
    "kind",
    "authors",
    "title",
    "venue",
    "volume",
    "issue",
    "pages",
    "year",
    "doi",
    "doi_status",
    "held",
    "held_files",
    "verification",
    "verification_note",
    "cited_by_declared",
    "code_mentions",
    "bibtex_key",
    "notes",
    "orphan_ok",
}
BIB_KINDS = {"journal", "conference", "report", "manual", "software", "thesis", "book"}
DOI_STATUSES = {"verified", "printed_on_pdf", "to_verify", "not_applicable"}
VERIFICATIONS = {"verified_pdf", "repository_citation", "unverified"}
DOI_RE = re.compile(r"^10\.\d{4,9}/\S+$")
CITATION_SITE_RE = re.compile(r"^[A-Za-z0-9_./-]+(:\d+(-\d+)?)?$")
CODE_ROOTS = ("src", "tests", "crates", "tools")
CODE_SUFFIXES = {".py", ".rs"}
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

# The dotted paths `set` may write. Anything absent is refused, which is what keeps
# a typo from adding a field nothing reads. A list spec means "an element of this
# list", and must be addressed with an index.
LEAF = "leaf"
SET_SCHEMA: dict[str, Any] = {
    "group": LEAF,
    "title": LEAF,
    "tests": LEAF,
    "evidence": LEAF,
    "validates": LEAF,
    "notes": LEAF,
    "flags": LEAF,
    "result": {
        "text": LEAF,
        "raw": LEAF,
        "unit": LEAF,
        "node": LEAF,
    },
    "comparisons": [
        {
            "label": LEAF,
            "asserted": LEAF,
            "expected": LEAF,
            "notes": LEAF,
            "reference": {
                "kind": LEAF,
                "label": LEAF,
                "citation": LEAF,
                "also_asserts": LEAF,
            },
            "tolerance": {
                "kind": LEAF,
                "value": LEAF,
                "source": LEAF,
                "justified": LEAF,
                "justification": LEAF,
            },
            "measured": {
                "status": LEAF,
                "run": LEAF,
                "date": LEAF,
                "raw": LEAF,
                "margin_pct": LEAF,
                "text": LEAF,
            },
        }
    ],
    "history": [{"rev": LEAF, "note": LEAF, "evidence": LEAF}],
}
DERIVED_FIELD_NAMES = {"near", "gt5", "slack", "slack_pp"}

ROW_FIELDS = [
    "id",
    "group",
    "title",
    "tests",
    "comparisons",
    "worst_margin_pct",
    "worst_slack_pp",
    "flags",
]
COMPARISON_FIELDS = [
    "row_id",
    "group",
    "label",
    "kind",
    "tolerance_pct",
    "margin_pct",
    "slack_pp",
    "asserted",
    "flags",
]


class StoreError(Exception):
    """The store could not be read, or an argument cannot be honoured."""


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
class RowRef:
    """One row, with enough provenance to write it back."""

    id: str
    where: str
    path: Path | None
    index: int
    data: dict[str, Any]


@dataclass
class Store:
    """The loaded store, plus the findings produced while loading it."""

    root: Path
    groups: dict[str, dict[str, Any]] = field(default_factory=dict)
    un_inventoried: dict[str, dict[str, Any]] = field(default_factory=dict)
    flag_entries: dict[str, dict[str, Any]] = field(default_factory=dict)
    derived_flags: set[str] = field(default_factory=set)
    references: dict[str, dict[str, Any]] = field(default_factory=dict)
    rows: list[RowRef] = field(default_factory=list)
    findings: list[Finding] = field(default_factory=list)

    def error(self, where: str, message: str) -> None:
        self.findings.append(Finding("error", where, message))

    def warn(self, where: str, message: str) -> None:
        self.findings.append(Finding("warning", where, message))

    def by_id(self) -> dict[str, RowRef]:
        return {ref.id: ref for ref in self.rows}

    def errors(self) -> list[Finding]:
        return [finding for finding in self.findings if finding.level == "error"]

    def warnings(self) -> list[Finding]:
        return [finding for finding in self.findings if finding.level == "warning"]


def load_yaml(path: Path) -> Any:
    try:
        with path.open("r", encoding="utf-8") as handle:
            return yaml.safe_load(handle)
    except FileNotFoundError as exc:
        raise StoreError(f"missing file: {path}") from exc
    except yaml.YAMLError as exc:
        raise StoreError(f"invalid YAML in {path}: {exc}") from exc


def dump_yaml(path: Path, data: Any) -> None:
    text = yaml.safe_dump(
        data,
        allow_unicode=True,
        default_flow_style=False,
        sort_keys=False,
        width=1000,
    )
    path.write_text(text, encoding="utf-8")


def _unknown_keys(data: dict[str, Any], allowed: set[str]) -> list[str]:
    return sorted(set(data) - allowed)


def display_path(path: Path) -> str:
    """Repository-relative when possible; absolute otherwise (a fixture store)."""
    try:
        return str(path.relative_to(REPO_ROOT))
    except ValueError:
        return str(path)


def load_flags(store: Store) -> None:
    data = load_yaml(store.root / "flags.yaml")
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
    data = load_yaml(store.root / "groups.yaml")
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
        elif store.references and citation not in store.references:
            store.error(
                where,
                f"reference.citation {citation!r} does not resolve in references.yaml",
            )
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
        rel = display_path(path)
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
            if not isinstance(row, dict) or not isinstance(row.get("id"), str):
                continue
            store.rows.append(
                RowRef(id=row["id"], where=where, path=path, index=index, data=row)
            )

    seen: dict[str, str] = {}
    for ref in store.rows:
        if ref.id in seen:
            store.error(ref.where, f"duplicate row id, first seen at {seen[ref.id]}")
        else:
            seen[ref.id] = ref.where


def load_store(root: Path, only_group: str | None = None) -> Store:
    store = Store(root=root)
    if not root.exists():
        raise StoreError(f"store directory does not exist: {root}")
    load_flags(store)
    load_references(store)
    load_groups(store)
    if only_group is not None and only_group not in store.groups:
        raise StoreError(f"unknown group: {only_group}")
    load_rows(store, only_group=only_group)
    return store


# --------------------------------------------------------------------------- #
# Derived values
# --------------------------------------------------------------------------- #


def _number(value: Any) -> float | None:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        return None
    return float(value)


def comparison_views(ref: RowRef) -> list[dict[str, Any]]:
    """One derived view per comparison: tolerance, margin, slack, derived flags."""
    views: list[dict[str, Any]] = []
    comparisons = ref.data.get("comparisons")
    if not isinstance(comparisons, list):
        return views
    for index, comparison in enumerate(comparisons):
        if not isinstance(comparison, dict):
            continue
        tolerance = comparison.get("tolerance") or {}
        measured = comparison.get("measured") or {}
        kind = tolerance.get("kind")
        value = _number(tolerance.get("value"))
        margin = _number(measured.get("margin_pct"))
        tolerance_pct = 100.0 * value if kind in RATIO_TOLERANCE_KINDS and value else None
        slack_pp = None
        if tolerance_pct is not None and margin is not None:
            slack_pp = tolerance_pct - margin
        derived: list[str] = []
        if slack_pp is not None and slack_pp < 1.0:
            derived.append("near")
        if tolerance_pct is not None and tolerance_pct > 100.0 * SUITE_TOLERANCE_RULE:
            derived.append("gt5")
        views.append(
            {
                "index": index,
                "row_id": ref.id,
                "group": ref.data.get("group"),
                "label": comparison.get("label"),
                "kind": kind,
                "asserted": comparison.get("asserted"),
                "status": measured.get("status"),
                "tolerance_pct": tolerance_pct,
                "margin_pct": margin,
                "slack_pp": slack_pp,
                "flags": derived,
                "text": measured.get("text"),
            }
        )
    return views


def row_view(ref: RowRef) -> dict[str, Any]:
    """One derived view per row: counts, worst margin, tightest slack, flags."""
    views = comparison_views(ref)
    margins = [view["margin_pct"] for view in views if view["margin_pct"] is not None]
    slacks = [view["slack_pp"] for view in views if view["slack_pp"] is not None]
    derived = sorted({flag for view in views for flag in view["flags"]})
    stored = ref.data.get("flags") or []
    return {
        "id": ref.id,
        "group": ref.data.get("group"),
        "title": ref.data.get("title"),
        "tests": len(ref.data.get("tests") or []),
        "comparisons": len(views),
        "worst_margin_pct": max(margins) if margins else None,
        "worst_slack_pp": min(slacks) if slacks else None,
        "flags": sorted(set(stored) | set(derived)),
        "derived_flags": derived,
    }


def all_comparison_views(store: Store) -> list[dict[str, Any]]:
    return [view for ref in store.rows for view in comparison_views(ref)]


# --------------------------------------------------------------------------- #
# Formatting
# --------------------------------------------------------------------------- #


def fmt_number(value: Any) -> str:
    if value is None:
        return "-"
    return f"{value:g}"


def fmt_pct(value: Any, signed: bool = False) -> str:
    if value is None:
        return "-"
    return f"{value:+g}%" if signed else f"{value:g}%"


def _field_value(unit: str, view: dict[str, Any], name: str) -> Any:
    if name == "kind" and unit == "comparison":
        return view.get("kind")
    return view.get(name)


def render_views(
    unit: str,
    views: list[dict[str, Any]],
    fields: list[str] | None,
) -> list[str]:
    names = fields or (ROW_FIELDS if unit == "row" else COMPARISON_FIELDS)
    lines = []
    for view in views:
        cells = []
        for name in names:
            value = _field_value(unit, view, name)
            if value is None:
                cells.append("-")
            elif name.endswith("_pct"):
                cells.append(fmt_pct(value, signed=name == "slack_pp"))
            elif name.endswith("_pp"):
                cells.append(fmt_pct(value, signed=True))
            elif isinstance(value, list):
                cells.append(",".join(str(item) for item in value) or "-")
            elif isinstance(value, bool):
                cells.append("yes" if value else "no")
            else:
                cells.append(str(value))
        lines.append(" | ".join(cells))
    return lines


# --------------------------------------------------------------------------- #
# Commands
# --------------------------------------------------------------------------- #


def command_check(args: argparse.Namespace) -> int:
    store = load_store(args.store, only_group=args.group)
    errors = store.errors()
    warnings = store.warnings()
    if args.json:
        payload = {
            "store": str(store.root),
            "group": args.group,
            "rows": len(store.rows),
            "comparisons": len(all_comparison_views(store)),
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
            f"{len(all_comparison_views(store))} comparison(s), "
            f"{len(errors)} error(s), {len(warnings)} warning(s)"
        )
    return EXIT_FINDINGS if errors else EXIT_OK


def _matches_file(store: Store, ref: RowRef, needle: str) -> bool:
    group = store.groups.get(str(ref.data.get("group"))) or {}
    if needle in [str(item) for item in group.get("source_files") or []]:
        return True
    return any(needle in node.split("::")[0] for node in ref.data.get("tests") or [])


def select_rows(
    store: Store,
    *,
    group: str | None = None,
    file: str | None = None,
    row_id: str | None = None,
    flag: str | None = None,
    near: bool = False,
    gt5: bool = False,
) -> list[RowRef]:
    refs = list(store.rows)
    if group:
        refs = [ref for ref in refs if ref.data.get("group") == group]
    if file:
        refs = [ref for ref in refs if _matches_file(store, ref, file)]
    if row_id:
        refs = [ref for ref in refs if ref.id == row_id]
    if flag:
        refs = [
            ref
            for ref in refs
            if flag in (ref.data.get("flags") or [])
            or any(flag in view["flags"] for view in comparison_views(ref))
        ]
    if near:
        refs = [ref for ref in refs if "near" in row_view(ref)["derived_flags"]]
    if gt5:
        refs = [ref for ref in refs if "gt5" in row_view(ref)["derived_flags"]]
    return refs


def command_find(args: argparse.Namespace) -> int:
    store = load_store(args.store)
    needle = args.substring.lower()
    hits: list[RowRef] = []
    for ref in select_rows(store, group=args.group, file=args.file):
        haystack = " ".join(
            [
                ref.id,
                str(ref.data.get("group")),
                str(ref.data.get("title")),
                str(ref.data.get("validates")),
                " ".join(ref.data.get("tests") or []),
            ]
        ).lower()
        if needle in haystack:
            hits.append(ref)
    if args.limit:
        hits = hits[: args.limit]
    if args.json:
        print(json.dumps([row_view(ref) for ref in hits], indent=2))
    else:
        for ref in hits:
            view = row_view(ref)
            print(
                f"{view['id']} | §{view['group']} | {view['title']} | "
                f"{view['tests']} test(s), {view['comparisons']} comparison(s)"
            )
        print(f"find: {len(hits)} match(es)")
    return EXIT_OK


def command_list(args: argparse.Namespace) -> int:
    store = load_store(args.store)
    refs = select_rows(
        store,
        group=args.group,
        file=args.file,
        row_id=args.id,
        flag=args.flag,
        near=args.near,
        gt5=args.gt5,
    )
    unit = args.unit
    if unit == "comparison":
        views = [view for ref in refs for view in comparison_views(ref)]
        # The row-level filter is a coarse prefilter; narrowing to the comparisons that
        # actually carry the property is what makes `--near` answer "which comparison is
        # close to failing", not "which row contains one".
        if args.near:
            views = [view for view in views if "near" in view["flags"]]
        if args.gt5:
            views = [view for view in views if "gt5" in view["flags"]]
        if args.flag:
            views = [view for view in views if args.flag in view["flags"]]
    else:
        views = [row_view(ref) for ref in refs]
    if args.sort == "slack":
        key = "slack_pp" if unit == "comparison" else "worst_slack_pp"
        views.sort(key=lambda view: (view.get(key) is None, view.get(key)))
    elif args.sort == "id":
        key = "id" if unit == "row" else "row_id"
        views.sort(key=lambda view: str(view.get(key)))
    fields = [name.strip() for name in args.fields.split(",")] if args.fields else None
    if fields:
        known = set(ROW_FIELDS if unit == "row" else COMPARISON_FIELDS)
        unknown = [name for name in fields if name not in known]
        if unknown:
            raise StoreError(f"unknown field(s): {', '.join(unknown)}")
    if args.json:
        payload = [{name: view.get(name) for name in (fields or view)} for view in views]
        print(json.dumps(payload, indent=2))
    else:
        for line in render_views(unit, views, fields):
            print(line)
        print(f"list: {len(views)} {unit}(s)")
    return EXIT_OK


def command_get(args: argparse.Namespace) -> int:
    store = load_store(args.store)
    ref = store.by_id().get(args.id)
    if ref is None:
        raise StoreError(f"unknown row id: {args.id}")
    if args.comparisons:
        views = comparison_views(ref)
        if args.json:
            print(json.dumps(views, indent=2))
        else:
            for line in render_views("comparison", views, None):
                print(line)
        return EXIT_OK
    if args.json:
        print(json.dumps(ref.data, indent=2))
    else:
        print(yaml.safe_dump(ref.data, allow_unicode=True, sort_keys=False, width=1000).rstrip())
    return EXIT_OK


def command_headline(args: argparse.Namespace) -> int:
    store = load_store(args.store)
    rows: list[dict[str, Any]] = []
    for group_id, group in store.groups.items():
        group_refs = [ref for ref in store.rows if ref.data.get("group") == group_id]
        margins = [
            view["margin_pct"]
            for ref in group_refs
            for view in comparison_views(ref)
            if view["margin_pct"] is not None
        ]
        flagged = sum(
            1
            for ref in group_refs
            for view in comparison_views(ref)
            if "gt5" in view["flags"]
        )
        rows.append(
            {
                "group": group_id,
                "headline": " ".join(str(group.get("headline") or "").split()) or "-",
                "rows": len(group_refs),
                "comparisons": sum(len(comparison_views(ref)) for ref in group_refs),
                "worst_margin_pct": max(margins) if margins else None,
                "comparisons_above_5pct": flagged,
            }
        )
    if args.json:
        print(json.dumps(rows, indent=2))
    else:
        print("| § | headline | rows | cmp | worst margin | >5% |")
        print("| --- | --- | --- | --- | --- | --- |")
        for row in rows:
            print(
                f"| {row['group']} | {row['headline']} | {row['rows']} | "
                f"{row['comparisons']} | {fmt_pct(row['worst_margin_pct'])} | "
                f"{row['comparisons_above_5pct']} |"
            )
    return EXIT_OK


def resolve_path(path: str) -> list[tuple[str, int | None]]:
    """Walk a dotted path against SET_SCHEMA, refusing anything undeclared."""
    segments: list[tuple[str, int | None]] = []
    for raw in path.split("."):
        match = PATH_SEGMENT_RE.match(raw)
        if match is None:
            raise StoreError(f"invalid path segment: {raw!r}")
        index = int(match.group(2)) if match.group(2) is not None else None
        segments.append((match.group(1), index))
    node: Any = SET_SCHEMA
    for name, index in segments:
        if not isinstance(node, dict):
            raise StoreError(f"cannot descend past {name!r}: it is a scalar")
        if name not in node:
            raise StoreError(
                f"unknown field {name!r} in {path!r}: allowed here: "
                f"{', '.join(sorted(node))}"
            )
        child = node[name]
        if isinstance(child, list):
            if index is None:
                raise StoreError(f"{name!r} is a list: address an element, e.g. {name}[0]")
            node = child[0]
        else:
            if index is not None:
                raise StoreError(f"{name!r} is not a list")
            node = child
    return segments


def apply_path(row: dict[str, Any], path: str, value: Any, *, unset: bool = False) -> None:
    first = path.split(".", 1)[0].split("[", 1)[0]
    if first in DERIVED_FIELD_NAMES:
        raise StoreError(f"{first!r} is derived and must not be written")
    if first == "id":
        raise StoreError("'id' is immutable: rename with a history entry and a redirect")
    segments = resolve_path(path)
    cursor: Any = row
    for name, index in segments[:-1]:
        if not isinstance(cursor, dict):
            raise StoreError(f"cannot descend past {path!r}")
        if index is None:
            if name not in cursor:
                cursor[name] = {}
            cursor = cursor[name]
            continue
        container = cursor.get(name)
        if not isinstance(container, list):
            raise StoreError(f"no such list: {name}")
        if index >= len(container):
            raise StoreError(f"index out of range: {name}[{index}]")
        cursor = container[index]
    name, index = segments[-1]
    if not isinstance(cursor, dict):
        raise StoreError(f"cannot set {path!r}")
    if index is None:
        if unset:
            cursor.pop(name, None)
        else:
            cursor[name] = value
        return
    container = cursor.get(name)
    if not isinstance(container, list):
        raise StoreError(f"no such list: {name}")
    if index >= len(container):
        raise StoreError(f"index out of range: {name}[{index}]")
    if unset:
        raise StoreError("cannot unset a list element; unset a field inside it instead")
    container[index] = value


def parse_value(raw: str) -> Any:
    if len(raw) >= 2 and raw[0] == raw[-1] and raw[0] in {'"', "'"}:
        return raw[1:-1]
    lowered = raw.lower()
    if lowered in {"true", "false"}:
        return lowered == "true"
    if lowered in {"null", "none"}:
        return None
    if raw[:1] in {"[", "{"}:
        try:
            return json.loads(raw)
        except json.JSONDecodeError as exc:
            raise StoreError(f"invalid JSON value {raw!r}: {exc}") from exc
    try:
        return int(raw)
    except ValueError:
        pass
    try:
        return float(raw)
    except ValueError:
        return raw


def command_set(args: argparse.Namespace) -> int:
    store = load_store(args.store)
    ref = store.by_id().get(args.id)
    if ref is None:
        raise StoreError(f"unknown row id: {args.id}")
    if ref.path is None:
        raise StoreError(f"row {args.id} has no backing file to write")

    if store.errors():
        for finding in store.errors():
            print(finding.as_text(), file=sys.stderr)
        raise StoreError("refusing to write: the store already has errors")

    document = load_yaml(ref.path)
    rows = document.get("rows") if isinstance(document, dict) else document
    if not isinstance(rows, list) or ref.index >= len(rows):
        raise StoreError(f"{ref.path} is not a row file this tool can write back")
    row = rows[ref.index]
    if not isinstance(row, dict):
        raise StoreError(f"{ref.path}[{ref.index}] is not a mapping")

    for assignment in args.assign:
        if "=" not in assignment:
            raise StoreError(f"expected PATH=VALUE, got {assignment!r}")
        path, raw = assignment.split("=", 1)
        apply_path(row, path, parse_value(raw))
    for path in args.unset:
        apply_path(row, path, None, unset=True)
    for flag in args.add_flag:
        if flag in store.derived_flags:
            raise StoreError(f"{flag!r} is a derived flag and must not be stored")
        if flag not in store.flag_entries:
            raise StoreError(f"unknown flag id: {flag!r}")
        flags = row.setdefault("flags", [])
        if flag not in flags:
            flags.append(flag)
    for flag in args.remove_flag:
        flags = row.get("flags") or []
        if flag in flags:
            flags.remove(flag)

    probe = Store(root=store.root, groups=store.groups, un_inventoried=store.un_inventoried)
    probe.flag_entries = store.flag_entries
    probe.derived_flags = store.derived_flags
    validate_row(probe, str(ref.path), row)
    if probe.errors():
        for finding in probe.errors():
            print(finding.as_text(), file=sys.stderr)
        raise StoreError("refusing to write: the edited row is invalid")

    dump_yaml(ref.path, document)
    print(f"set {args.id}: wrote {display_path(ref.path)}")
    return EXIT_OK


# --------------------------------------------------------------------------- #
# Bibliography (`docs/validation/references.yaml`)
# --------------------------------------------------------------------------- #


def _validate_bib_entry(store: Store, where: str, entry: Any) -> None:
    if not isinstance(entry, dict):
        store.error(where, "expected a mapping")
        return
    unknown = _unknown_keys(entry, BIB_KEYS)
    if unknown:
        store.error(where, f"unknown keys: {', '.join(unknown)}")

    key = entry.get("key")
    if not isinstance(key, str) or not SLUG_RE.match(key):
        store.error(where, f"invalid key: {key!r}")
    else:
        where = f"references.yaml[{key}]"
        if key in store.references:
            store.error(where, "duplicate key")

    if not isinstance(entry.get("section"), str) or not entry["section"]:
        store.error(where, "section must be a non-empty string")
    if entry.get("kind") not in BIB_KINDS:
        store.error(where, f"kind must be one of {sorted(BIB_KINDS)}")
    authors = entry.get("authors")
    if not isinstance(authors, list) or not authors or not all(
        isinstance(author, str) and author for author in authors
    ):
        store.error(where, "authors must be a non-empty list of non-empty strings")
    if not isinstance(entry.get("title"), str) or not entry["title"]:
        store.error(where, "title must be a non-empty string")
    year = entry.get("year")
    if not isinstance(year, int) or isinstance(year, bool) or not 1900 <= year <= 2100:
        store.error(where, f"year must be an integer in 1900..2100, got {year!r}")

    doi = entry.get("doi")
    if doi is not None and (not isinstance(doi, str) or not DOI_RE.match(doi)):
        store.error(where, f"doi must be a bare DOI, got {doi!r}")
    status = entry.get("doi_status")
    if status not in DOI_STATUSES:
        store.error(where, f"doi_status must be one of {sorted(DOI_STATUSES)}")
        return
    held = entry.get("held")
    if not isinstance(held, bool):
        store.error(where, "held must be a boolean")
    held_files = entry.get("held_files")
    if held_files is not None and not isinstance(held_files, list):
        store.error(where, "held_files must be a list")

    if status in {"verified", "printed_on_pdf"}:
        if status == "verified" and not (isinstance(held, bool) and held):
            store.error(where, "doi_status: verified requires held: true")
        if held_files is not None and not held_files:
            store.error(where, f"doi_status: {status} requires a non-empty held_files")
        if entry.get("verification") not in {"verified_pdf", "repository_citation"}:
            store.error(where, f"doi_status: {status} requires a 'verification' method")
        if not isinstance(entry.get("verification_note"), str) or not entry["verification_note"]:
            store.error(where, f"doi_status: {status} requires a 'verification_note'")
    elif status == "to_verify":
        if doi is not None:
            store.error(where, "doi_status: to_verify requires doi: null")
    elif status == "not_applicable":
        if doi is not None:
            store.error(where, "doi_status: not_applicable requires doi: null")
        if not isinstance(entry.get("notes"), str) or not entry["notes"]:
            store.error(where, "doi_status: not_applicable requires 'notes' explaining why")

    for site in entry.get("cited_by_declared") or []:
        if not isinstance(site, str) or not CITATION_SITE_RE.match(site):
            store.error(where, f"not a file or file:line site: {site!r}")
            continue
        target = REPO_ROOT / site.split(":", 1)[0]
        if not target.exists():
            store.error(where, f"declared citation site does not exist: {site}")
    mentions = entry.get("code_mentions")
    if mentions is not None and (
        not isinstance(mentions, list)
        or not all(isinstance(item, str) and item for item in mentions)
    ):
        store.error(where, "code_mentions must be a list of non-empty strings")
    if entry.get("orphan_ok") is not None and not isinstance(entry.get("orphan_ok"), bool):
        store.error(where, "orphan_ok must be a boolean")


def load_references(store: Store) -> None:
    data = load_yaml(store.root / "references.yaml")
    if not isinstance(data, dict) or not isinstance(data.get("references"), list):
        store.error("references.yaml", "expected a mapping with a 'references' list")
        return
    for index, entry in enumerate(data["references"]):
        _validate_bib_entry(store, f"references.yaml[{index}]", entry)
        if isinstance(entry, dict) and isinstance(entry.get("key"), str):
            store.references[entry["key"]] = entry


def iter_code_files() -> Iterable[Path]:
    for root in CODE_ROOTS:
        base = REPO_ROOT / root
        if not base.exists():
            continue
        for path in sorted(base.rglob("*")):
            if path.is_file() and path.suffix in CODE_SUFFIXES:
                yield path


def find_code_mentions(store: Store) -> dict[str, list[str]]:
    """Sites that mention a work, by its key or by a declared `code_mentions` string.

    The key alone is useless as a needle: the code cites author-year text
    ("Ko2017 ratio-based mesh distortion"), never the store's key. `code_mentions`
    is therefore the auditable list of literal strings the code actually uses.
    """
    needles = {
        key: [key, *[str(item) for item in entry.get("code_mentions") or []]]
        for key, entry in store.references.items()
    }
    found: dict[str, list[str]] = {key: [] for key in store.references}
    for path in iter_code_files():
        try:
            lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
        except OSError:
            continue
        rel = display_path(path)
        for number, line in enumerate(lines, start=1):
            for key, phrases in needles.items():
                if any(phrase and phrase in line for phrase in phrases):
                    found[key].append(f"{rel}:{number}")
    return found


def rows_citing(store: Store, key: str) -> list[str]:
    """(row id, comparison label) pairs that cite a key."""
    users: list[str] = []
    for ref in store.rows:
        for view in comparison_views(ref):
            comparison = (ref.data.get("comparisons") or [])[view["index"]]
            citation = (comparison.get("reference") or {}).get("citation")
            if citation == key:
                users.append(f"{ref.id} [{view['label']}]")
    return users


def bibtex_entry(store: Store, key: str) -> str:
    entry = store.references[key]
    kind = entry.get("kind")
    template = {
        "journal": "article",
        "conference": "inproceedings",
        "report": "techreport",
        "book": "book",
    }.get(str(kind), "misc")
    fields: list[tuple[str, str]] = [
        ("author", " and ".join(entry.get("authors") or [])),
        ("title", str(entry.get("title"))),
    ]
    if entry.get("venue"):
        fields.append(("journal" if template == "article" else "booktitle", str(entry["venue"])))
    for name in ("volume", "issue", "pages", "year"):
        if entry.get(name) is not None:
            fields.append((name, str(entry[name])))
    if entry.get("doi"):
        fields.append(("doi", str(entry["doi"])))
    body = ".\n".join(f"  {name} = {{{value}}}" for name, value in fields)
    return f"@{template}{{{entry.get('bibtex_key') or key},\n{body}\n}}"


def reference_scan_findings(store: Store, found: dict[str, list[str]]) -> list[Finding]:
    findings: list[Finding] = []
    used = set()
    for ref in store.rows:
        for view in comparison_views(ref):
            comparison = (ref.data.get("comparisons") or [])[view["index"]]
            citation = (comparison.get("reference") or {}).get("citation")
            if isinstance(citation, str):
                used.add(citation)
    for key, entry in sorted(store.references.items()):
        where = f"references.yaml[{key}]"
        sites = [str(site) for site in entry.get("cited_by_declared") or []]
        for site in sites:
            target = REPO_ROOT / site.split(":", 1)[0]
            if ":" in site:
                line_no = int(site.rsplit(":", 1)[1].split("-", 1)[0])
                try:
                    lines = target.read_text(encoding="utf-8", errors="replace").splitlines()
                except OSError:
                    lines = []
                if line_no > len(lines):
                    findings.append(
                        Finding("error", where, f"declared site {site} is past end of file")
                    )
                    continue
            file_part, _, line_part = site.partition(":")
            if line_part:
                ok = site in found.get(key, [])
            else:
                ok = any(mention.split(":", 1)[0] == file_part for mention in found.get(key, []))
            if not ok:
                findings.append(
                    Finding(
                        "error",
                        where,
                        f"declared citation site {site} no longer mentions this work "
                        f"(key {key!r} or code_mentions)",
                    )
                )
        if not used and not sites and not entry.get("orphan_ok"):
            findings.append(
                Finding(
                    "warning",
                    where,
                    "orphan: no comparison cites it and no site is declared; set orphan_ok "
                    "with a reason if the bibliography legitimately holds it",
                )
            )
        for path in entry.get("held_files") or []:
            if not (REPO_ROOT / str(path)).exists():
                findings.append(
                    Finding("warning", where, f"held file absent (.sources is gitignored): {path}")
                )
    return findings


def command_references(args: argparse.Namespace) -> int:
    store = load_store(args.store)
    action = args.action
    keys = sorted(store.references)
    if action == "get":
        if not args.key:
            raise StoreError("references get needs a key")
        entry = store.references.get(args.key)
        if entry is None:
            raise StoreError(f"unknown reference key: {args.key}")
        if args.json:
            print(json.dumps(entry, indent=2))
        else:
            print(yaml.safe_dump(entry, allow_unicode=True, sort_keys=False, width=1000).rstrip())
        return EXIT_OK
    if action == "list" or action == "gaps":
        selected = []
        for key in keys:
            entry = store.references[key]
            if action == "gaps" and entry.get("doi_status") not in {"to_verify", "not_applicable"}:
                continue
            if args.doi_status and entry.get("doi_status") != args.doi_status:
                continue
            if args.kind and entry.get("kind") != args.kind:
                continue
            if args.held and not entry.get("held"):
                continue
            if args.not_held and entry.get("held"):
                continue
            selected.append(entry)
        if args.json:
            print(json.dumps(selected, indent=2))
        else:
            for entry in selected:
                print(
                    f"{entry['key']} | {entry.get('year')} | {entry.get('kind')} | "
                    f"{entry.get('doi') or entry.get('doi_status')} | {entry.get('title')}"
                )
            print(f"references {action}: {len(selected)} entr(ies)")
        return EXIT_OK
    if action == "bibtex":
        text = "\n\n".join(bibtex_entry(store, key) for key in keys) + "\n"
        if args.out:
            args.out.write_text(text, encoding="utf-8")
            print(f"wrote {display_path(args.out)}")
        else:
            print(text, end="")
        return EXIT_OK
    if action == "where-used":
        if not args.key:
            raise StoreError("references where-used needs a key")
        if args.key not in store.references:
            raise StoreError(f"unknown reference key: {args.key}")
        found = find_code_mentions(store)[args.key]
        payload = {"key": args.key, "rows": rows_citing(store, args.key), "code": found}
        if args.json:
            print(json.dumps(payload, indent=2))
        else:
            for user in payload["rows"]:
                print(f"row  | {user}")
            for site in payload["code"]:
                print(f"code | {site}")
        return EXIT_OK
    if action == "check":
        findings = store.findings + reference_scan_findings(store, find_code_mentions(store))
        errors = [finding for finding in findings if finding.level == "error"]
        if args.json:
            print(json.dumps([finding.as_dict() for finding in findings], indent=2))
        else:
            for finding in findings:
                print(finding.as_text())
            print(
                f"references check: {len(store.references)} entr(ies), "
                f"{len(errors)} error(s), {len(findings) - len(errors)} warning(s)"
            )
        return EXIT_FINDINGS if errors else EXIT_OK
    raise StoreError(f"unknown references action: {action}")


# --------------------------------------------------------------------------- #
# Entry point
# --------------------------------------------------------------------------- #


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

    find = subparsers.add_parser("find", help="search rows by id, title, group or node id")
    find.add_argument("substring")
    find.add_argument("--group")
    find.add_argument("--file")
    find.add_argument("--limit", type=int, default=0)
    find.add_argument("--json", action="store_true")
    find.set_defaults(func=command_find)

    listing = subparsers.add_parser("list", help="one line per row, or per comparison")
    listing.add_argument("--group")
    listing.add_argument("--file")
    listing.add_argument("--id")
    listing.add_argument("--flag")
    listing.add_argument("--near", action="store_true", help="only comparisons with slack < 1pp")
    listing.add_argument("--gt5", action="store_true", help="only comparisons above the 5% rule")
    listing.add_argument("--unit", choices=["row", "comparison"], default="row")
    listing.add_argument("--sort", choices=["id", "slack"], default="id")
    listing.add_argument("--fields", help="comma-separated field names")
    listing.add_argument("--json", action="store_true")
    listing.set_defaults(func=command_list)

    get = subparsers.add_parser("get", help="print one row verbatim")
    get.add_argument("id")
    get.add_argument("--comparisons", action="store_true", help="derived comparison table")
    get.add_argument("--json", action="store_true")
    get.set_defaults(func=command_get)

    headline = subparsers.add_parser("headline", help="derived section 2.1 summary table")
    headline.add_argument("--json", action="store_true")
    headline.set_defaults(func=command_headline)

    bibliography = subparsers.add_parser("references", help="query the bibliography store")
    bibliography.add_argument(
        "action", choices=["list", "gaps", "get", "check", "bibtex", "where-used"]
    )
    bibliography.add_argument("key", nargs="?", help="reference key for get and where-used")
    bibliography.add_argument("--doi-status", choices=sorted(DOI_STATUSES))
    bibliography.add_argument("--kind", choices=sorted(BIB_KINDS))
    bibliography.add_argument("--held", action="store_true", help="only entries whose PDF is held")
    bibliography.add_argument("--not-held", action="store_true")
    bibliography.add_argument("--out", type=Path, help="write bibtex to a file")
    bibliography.add_argument("--json", action="store_true")
    bibliography.set_defaults(func=command_references)

    setter = subparsers.add_parser("set", help="edit one row in place")
    setter.add_argument("id")
    setter.add_argument("assign", nargs="*", metavar="PATH=VALUE")
    setter.add_argument("--unset", action="append", default=[], metavar="PATH")
    setter.add_argument("--add-flag", action="append", default=[], metavar="FLAG")
    setter.add_argument("--remove-flag", action="append", default=[], metavar="FLAG")
    setter.set_defaults(func=command_set)

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
