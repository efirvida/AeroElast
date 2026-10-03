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

`docs/validation-matrix.md` and `docs/references.md` were the store's generated views, and
they are gone: the store is the only artifact. `references render` survives as an export
command, and `diff-against-md` as a migration-time cross-check, which is the one thing here
that still needs the old file.

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
import ast
import json
import re
import subprocess
import sys
from datetime import date
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

REFERENCE_KINDS = {"paper", "code", "analytical", "self"}
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
    "version",
    "doi",
    "doi_status",
    "held",
    "held_files",
    "verification",
    "verification_note",
    "cited_by_declared",
    "cited_by_stale",
    "code_mentions",
    "bibtex_key",
    "notes",
    "orphan_ok",
    "bibliographic_gaps",
}
BIB_GAPPABLE = {"authors", "title", "year", "venue", "volume", "pages"}
BIB_KINDS = {"journal", "conference", "report", "manual", "software", "thesis", "book"}
DOI_STATUSES = {"verified", "printed_on_pdf", "verified_externally", "to_verify", "not_applicable"}
VERIFICATIONS = {"verified_pdf", "repository_citation", "external_record", "unverified"}
DOI_RE = re.compile(r"^10\.\d{4,9}/\S+$")
CITATION_SITE_RE = re.compile(r"^[A-Za-z0-9_./-]+(:\d+(-\d+)?)?$")
CODE_ROOTS = ("src", "tests", "crates", "tools")
# Manifests count as a citation site: Cargo.toml names the crates that carry a work.
# Markdown does not, because the bibliography itself names every work and would match
# every entry trivially.
CODE_SUFFIXES = {".py", ".rs", ".toml"}
TOLERANCE_KEYS = {"kind", "value", "source", "justified", "justification"}
MEASURED_KEYS = {"status", "run", "date", "raw", "margin_pct", "text"}
HISTORY_KEYS = {"rev", "note", "evidence"}
GROUP_KEYS = {
    "id",
    "title",
    "slug",
    "citation",
    "source_files",
    "common_tolerance",
    "headline",
    "provenance_note",
    "prose_after_table",
}
# A group declares how its tests print their residuals, because the suite prints prose, not
# a format. `asserted` matches the line behind a real assertion, `unasserted` the line the
# test prints and never asserts.
# The patterns themselves live in `residual-patterns.json`, keyed by group id: a regex
# is machine data, and the YAML prose linter enforces a line budget a long pattern cannot
# satisfy. `asserted` matches the line behind a real assertion, `unasserted` the line the
# test prints and never asserts.
RESIDUAL_KEYS = {"asserted", "unasserted"}
RESIDUAL_FILE = "residual-patterns.json"
NODE_START_RE = re.compile(r"^(?P<node>tests/[^\s:]+\.py::\S+)(?:\s+(?P<rest>.*))?$")
NODE_STATUS_RE = re.compile(r"^(PASSED|FAILED|XFAIL|XPASS|SKIPPED|ERROR)\b")
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
    reference_header: str | None = None
    reference_sections: dict[str, dict[str, Any]] = field(default_factory=dict)
    residual_patterns: dict[str, dict[str, str]] = field(default_factory=dict)
    gaps: dict[str, dict[str, Any]] = field(default_factory=dict)
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


def is_iso_date(value: Any) -> bool:
    """A date as the store writes it, or as YAML hands it back.

    `yaml.safe_dump` quotes a date string so it round-trips as a string, but a date written
    by hand in a YAML file is unquoted and resolves to a `datetime.date`. Rejecting that
    would make the store's own files invalid depending on who wrote them.
    """
    if isinstance(value, date):
        return True
    return isinstance(value, str) and bool(DATE_RE.match(value))


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
    date_value = measured.get("date")
    if date_value is not None and not is_iso_date(date_value):
        store.error(where, f"measured.date must be YYYY-MM-DD, got {date_value!r}")
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
    load_residual_patterns(store)
    load_rows(store, only_group=only_group)
    load_adjudications(store)
    load_gaps(store)
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
    gaps = entry.get("bibliographic_gaps")
    if gaps is not None:
        if not isinstance(gaps, list) or not all(isinstance(gap, str) and gap for gap in gaps):
            store.error(where, "bibliographic_gaps must be a list of non-empty strings")
            gaps = None
        else:
            unknown = [gap for gap in gaps if gap not in BIB_GAPPABLE]
            if unknown:
                store.error(
                    where,
                    f"bibliographic_gaps may only name {sorted(BIB_GAPPABLE)}, got {unknown}",
                )
    declared_gaps = set(gaps or [])

    def absent(field: str, value: Any) -> bool:
        """A field the source does not state: allowed only when declared as a gap."""
        if value not in (None, "", [], {}):
            return False
        if field not in declared_gaps:
            store.error(
                where,
                f"{field!r} is missing and is not declared in bibliographic_gaps; "
                "an absence is recorded, never left implicit",
            )
            return False
        return True

    authors = entry.get("authors")
    if not isinstance(authors, list) or not all(isinstance(item, str) and item for item in authors):
        store.error(where, "authors must be a list of non-empty strings")
    elif not authors:
        absent("authors", "")

    title = entry.get("title")
    if title is not None and (not isinstance(title, str) or not title):
        store.error(where, "title must be a non-empty string or null")
    elif title is None:
        absent("title", None)

    year = entry.get("year")
    if year is not None and (
        not isinstance(year, int) or isinstance(year, bool) or not 1900 <= year <= 2100
    ):
        store.error(where, f"year must be an integer in 1900..2100 or null, got {year!r}")
    elif year is None:
        absent("year", None)
    version = entry.get("version")
    if version is not None and (not isinstance(version, str) or not version):
        store.error(where, "version must be a non-empty string or null")
    if entry.get("kind") == "software" and year is None and not version:
        store.error(
            where,
            "a software entry with no publication year must pin its identity with 'version'",
        )

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
        if status == "printed_on_pdf" and not (isinstance(held, bool) and held):
            store.error(
                where,
                "doi_status: printed_on_pdf claims a DOI on a held copy, so it requires "
                "held: true; use verified_externally when no copy is held",
            )
        if held_files is not None and not held_files:
            store.error(where, f"doi_status: {status} requires a non-empty held_files")
        if entry.get("verification") not in {"verified_pdf", "repository_citation"}:
            store.error(where, f"doi_status: {status} requires a 'verification' method")
        if not isinstance(entry.get("verification_note"), str) or not entry["verification_note"]:
            store.error(where, f"doi_status: {status} requires a 'verification_note'")
    elif status == "verified_externally":
        if doi is None:
            store.error(where, "doi_status: verified_externally requires the doi it verified")
        if entry.get("verification") != "external_record":
            store.error(
                where,
                "doi_status: verified_externally requires verification: external_record",
            )
        if not isinstance(entry.get("verification_note"), str) or not entry["verification_note"]:
            store.error(
                where,
                "doi_status: verified_externally requires a 'verification_note' naming "
                "the record it was read from",
            )
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
    stale = entry.get("cited_by_stale")
    if stale is not None:
        if not isinstance(stale, list) or not all(isinstance(site, str) and site for site in stale):
            store.error(where, "cited_by_stale must be a list of non-empty strings")
        else:
            for site in stale:
                if not CITATION_SITE_RE.match(site):
                    store.error(where, f"not a file or file:line site: {site!r}")
                elif not (REPO_ROOT / site.split(":", 1)[0]).exists():
                    store.error(where, f"stale citation site does not exist: {site}")
                elif site in (entry.get("cited_by_declared") or []):
                    store.error(
                        where,
                        f"{site} is both declared and stale; it is one or the other",
                    )
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
    header = data.get("header")
    if header is not None and (not isinstance(header, str) or not header.strip()):
        store.error("references.yaml", "'header' must be a non-empty string or absent")
    store.reference_header = header if isinstance(header, str) else None

    for index, section in enumerate(data.get("sections") or []):
        where = f"references.yaml:sections[{index}]"
        if not isinstance(section, dict):
            store.error(where, "expected a mapping")
            continue
        unknown = _unknown_keys(section, {"title", "intro"})
        if unknown:
            store.error(where, f"unknown keys: {', '.join(unknown)}")
        title = section.get("title")
        if not isinstance(title, str) or not title:
            store.error(where, "title must be a non-empty string")
            continue
        if title in store.reference_sections:
            store.error(f"references.yaml:sections[{title}]", "duplicate section title")
        intro = section.get("intro")
        if intro is not None and (not isinstance(intro, str) or not intro.strip()):
            store.error(f"references.yaml:sections[{title}]", "intro must be non-empty or absent")
        store.reference_sections[title] = section

    for index, entry in enumerate(data["references"]):
        _validate_bib_entry(store, f"references.yaml[{index}]", entry)
        if isinstance(entry, dict) and isinstance(entry.get("key"), str):
            store.references[entry["key"]] = entry

    declared = [
        str(entry.get("section")) for entry in store.references.values() if entry.get("section")
    ]
    for title in sorted(set(declared)):
        if title not in store.reference_sections:
            store.warn(
                "references.yaml",
                f"entries use section {title!r}, which `sections` does not declare: its "
                "intro cannot be rendered",
            )


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
        mentions_available = bool(found.get(key)) or bool(entry.get("code_mentions"))
        if sites and not mentions_available:
            # Without needles the scan can only look for the key string, which the code
            # never uses, so every declared site would read as stale. Report the inability
            # to check instead of manufacturing a verdict.
            findings.append(
                Finding(
                    "warning",
                    where,
                    "unverifiable: the entry declares citation sites but no "
                    "`code_mentions`, so nothing the scan can search for would appear in "
                    "the code. Declare the literal strings the code cites.",
                )
            )
        for site in sites:
            if not mentions_available:
                continue
            target = REPO_ROOT / site.split(":", 1)[0]
            file_part, _, line_part = site.partition(":")
            start, _, end = line_part.partition("-")
            if line_part and not start.isdigit():
                continue  # malformed; the loader already reported it
            try:
                lines = target.read_text(encoding="utf-8", errors="replace").splitlines()
            except OSError:
                lines = []
            mentioned = {
                mention.rsplit(":", 1)[-1]
                for mention in found.get(key, [])
                if mention.split(":", 1)[0] == file_part
            }
            if line_part:
                first = int(start)
                if first > len(lines):
                    findings.append(
                        Finding("error", where, f"declared site {site} is past end of file")
                    )
                    continue
                # A range is satisfied by a mention on ANY of its lines: reading only the
                # first line turned every range whose mention sits further down into a
                # false stale verdict.
                last = int(end) if end.isdigit() else first
                ok = any(str(number) in mentioned for number in range(first, last + 1))
            else:
                ok = bool(mentioned)
            if not ok:
                findings.append(
                    Finding(
                        "error",
                        where,
                        f"declared citation site {site} no longer mentions this work "
                        f"(key {key!r} or code_mentions)",
                    )
                )
        if not used and not sites and not entry.get("cited_by_stale") and not entry.get(
            "orphan_ok"
        ):
            findings.append(
                Finding(
                    "warning",
                    where,
                    "orphan: no comparison cites it and no site is declared; set orphan_ok "
                    "with a reason if the bibliography legitimately holds it",
                )
            )
        for site in entry.get("cited_by_stale") or []:
            if site in found.get(key, []):
                findings.append(
                    Finding(
                        "error",
                        where,
                        f"{site} is recorded as stale but the code does mention this work "
                        "there: move it to cited_by_declared",
                    )
                )
            else:
                findings.append(
                    Finding(
                        "warning",
                        where,
                        f"declared citation site {site} is known stale: the prose claims it "
                        "and the code no longer mentions the work there",
                    )
                )
        for path in entry.get("held_files") or []:
            if not (REPO_ROOT / str(path)).exists():
                findings.append(
                    Finding("warning", where, f"held file absent (.sources is gitignored): {path}")
                )
    return findings


def reference_heading(title: str) -> str:
    """Subsections of a numbered section (6.1, 6.2) render one level deeper."""
    return "###" if re.match(r"^\d+\.\d+\s", title) else "##"


def render_citation_line(entry: dict[str, Any]) -> str:
    authors = ", ".join(entry.get("authors") or [])
    title = entry.get("title")
    parts = [part for part in (authors, f'"{title}"' if title else None) if part]
    if entry.get("venue"):
        parts.append(f"*{entry['venue']}*")
    if entry.get("version"):
        parts.append(str(entry["version"]))
    volume = entry.get("volume")
    if volume:
        issue = entry.get("issue")
        parts.append(f"{volume}({issue})" if issue else str(volume))
    if entry.get("pages"):
        parts.append(str(entry["pages"]))
    line = ", ".join(parts)
    # A source that states no year prints none: never the string "None".
    tail = f", {entry['year']}." if entry.get("year") else "."
    if entry.get("doi"):
        tail += f" DOI: {entry['doi']}."
    elif entry.get("doi_status") == "to_verify":
        tail += " DOI: to verify."
    return line + tail


def render_reference_paragraphs(entry: dict[str, Any]) -> list[str]:
    paragraphs: list[str] = []
    note = entry.get("verification_note")
    held = (entry.get("held_files") or [None])[0]
    if entry.get("verification") == "verified_pdf" and held:
        detail = f" ({note})" if note else ""
        paragraphs.append(f"*Verified against `{held}`{detail}.*")
    elif entry.get("verification") == "external_record":
        # The note is already a full sentence describing what was read and where.
        paragraphs.append(
            f"*{note}*" if note else "*DOI verified against an external record; no copy held.*"
        )
    elif entry.get("verification") == "repository_citation":
        detail = f": {note}" if note else "."
        paragraphs.append(f"*Source: repository citation{detail}*")
    elif entry.get("verification") == "unverified":
        paragraphs.append("*Not verified against a held copy.*")
    if entry.get("cited_by_declared"):
        sites = ", ".join(f"`{site}`" for site in entry["cited_by_declared"])
        paragraphs.append(f"*Cited by the code: {sites}.*")
    if entry.get("cited_by_stale"):
        sites = ", ".join(f"`{site}`" for site in entry["cited_by_stale"])
        paragraphs.append(
            f"*Recorded stale: the prose cited {sites}, and the code no longer mentions "
            "this work there.*"
        )
    if entry.get("notes"):
        paragraphs.append(f"*{entry['notes']}*")
    return paragraphs


def render_references(store: Store) -> str:
    """Render the bibliography: an export, not a committed artifact."""
    lines = [
        "# Canonical bibliography",
        "",
        "Generated from `docs/validation/references.yaml` by",
        "`python tools/validation_matrix.py references render`. Do not edit by hand.",
        "",
    ]
    if store.reference_header:
        lines += [store.reference_header.strip(), ""]

    grouped: dict[str, list[dict[str, Any]]] = {}
    for entry in store.references.values():
        grouped.setdefault(str(entry.get("section")), []).append(entry)

    def emit_section(title: str, intro: str | None) -> None:
        lines.append(f"{reference_heading(title)} {title}")
        lines.append("")
        if intro:
            lines.extend([intro.strip(), ""])
        for entry in grouped.get(title, []):
            lines.append(f"- {render_citation_line(entry)}")
            for paragraph in render_reference_paragraphs(entry):
                lines.append(f"  {paragraph}")
            lines.append("")

    for title, section in store.reference_sections.items():
        emit_section(title, section.get("intro"))
    for title in grouped:
        if title not in store.reference_sections:
            emit_section(title, None)
    return "\n".join(lines).rstrip() + "\n"


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
    if action == "render":
        rendered = render_references(store)
        target = args.out or (REPO_ROOT / "docs" / "references.md")
        if args.check:
            current = target.read_text(encoding="utf-8") if target.exists() else ""
            if current != rendered:
                print(f"stale: {display_path(target)} differs from a fresh render", file=sys.stderr)
                return EXIT_FINDINGS
            print(f"render --check: {display_path(target)} is up to date")
            return EXIT_OK
        target.write_text(rendered, encoding="utf-8")
        print(f"wrote {display_path(target)}")
        return EXIT_OK
    raise StoreError(f"unknown references action: {action}")


# --------------------------------------------------------------------------- #
# Code extraction (pilot stage E1-E3)
# --------------------------------------------------------------------------- #


@dataclass
class NodeInfo:
    node: str
    file: str
    function: str
    params: str | None


@dataclass
class ToleranceSite:
    line: int
    kind: str
    value: float | None
    source: str
    asserts: bool


def collect_nodes(scope: str) -> list[NodeInfo]:
    """E1: the collected node set, verbatim from pytest.

    A zero-node result is an error, never a clean run: a module-level
    `pytest.importorskip` outside the pinned environment collects nothing, and a
    silent zero would read as "zero drift".
    """
    completed = subprocess.run(
        [sys.executable, "-m", "pytest", "-o", "addopts=", "--collect-only", "-q", scope],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
        check=False,
    )
    nodes: list[NodeInfo] = []
    for raw in completed.stdout.splitlines():
        line = raw.strip()
        if "::" not in line or line.startswith("="):
            continue
        if not NODE_ID_RE.match(line):
            continue
        file, _, rest = line.partition("::")
        function, _, params = rest.partition("[")
        nodes.append(
            NodeInfo(
                node=line,
                file=file,
                function=function,
                params=params.rstrip("]") if params else None,
            )
        )
    if not nodes:
        raise StoreError(
            f"collect-only returned no nodes for {scope!r} (exit {completed.returncode}); "
            "run it inside the pinned environment (conda activate aeroelast-dev)"
        )
    return nodes


def module_constants(tree: ast.Module) -> dict[str, Any]:
    consts: dict[str, Any] = {}
    for node in tree.body:
        if isinstance(node, ast.Assign) and len(node.targets) == 1:
            target = node.targets[0]
            if isinstance(target, ast.Name):
                try:
                    consts[target.id] = ast.literal_eval(node.value)
                except (ValueError, SyntaxError):
                    continue
    return consts


def resolve_literal(node: ast.expr, consts: dict[str, Any]) -> Any:
    """A literal, or a module constant holding one. Nothing else."""
    if isinstance(node, ast.Name) and node.id in consts:
        return consts[node.id]
    try:
        return ast.literal_eval(node)
    except (ValueError, SyntaxError):
        return None


def _is_number(text: str) -> bool:
    try:
        float(text)
    except ValueError:
        return False
    return True


def param_tokens(params: str | None, expected: int | None = None) -> list[str]:
    """Split a parametrisation id into parameter tokens.

    Ambiguous by construction: pytest joins parameters with `-`, an exponent carries
    its own minus (`1e-06`), and a value may too (`In-plane`). Knowing how many
    parameters the function declares resolves it: merge pieces until the count
    matches, merging either an exponent fragment or two adjacent non-numbers.
    """
    pieces = [piece for piece in (params or "").split("-") if piece]
    if expected is None:
        return pieces
    # Pass 1: rejoin an exponent fragment, because `1e-06` is one number and not two
    # tokens. Done first so that a value ending in `e` (`In-plane`) cannot claim the
    # exponent's minus as its own separator.
    rejoined: list[str] = []
    index = 0
    while index < len(pieces):
        joined = f"{pieces[index]}-{pieces[index + 1]}" if index + 1 < len(pieces) else ""
        if joined and not _is_number(pieces[index]) and _is_number(joined):
            rejoined.append(joined)
            index += 2
            continue
        rejoined.append(pieces[index])
        index += 1
    # Pass 2: rejoin the pieces of a value that itself contains the separator, in place
    # and only as long as the count still exceeds the number of declared parameters.
    tokens: list[str] = []
    index = 0
    while index < len(rejoined):
        piece = rejoined[index]
        needed = expected - len(tokens)
        if (
            index + 1 < len(rejoined)
            and len(rejoined) - index > needed
            and not _is_number(piece)
            and not _is_number(rejoined[index + 1])
        ):
            rejoined[index + 1] = f"{piece}-{rejoined[index + 1]}"
            index += 1
            continue
        tokens.append(piece)
        index += 1
    return tokens


def param_numbers(params: str | None, expected: int | None = None) -> set[float]:
    return {
        float(token) for token in param_tokens(params, expected) if _is_number(token)
    }


def _numbers_in(
    node: ast.AST,
    consts: dict[str, Any],
    functions: dict[str, ast.FunctionDef] | None,
    depth: int = 1,
) -> set[float]:
    numbers: set[float] = set()
    for sub in ast.walk(node):
        if isinstance(sub, ast.Name) and isinstance(consts.get(sub.id), list):
            for row in consts[sub.id]:
                for item in row if isinstance(row, tuple) else [row]:
                    value = eval_number(ast.Constant(item), consts)
                    if value is not None:
                        numbers.add(value)
        if depth and isinstance(sub, ast.Call) and isinstance(sub.func, ast.Name):
            # `_twisted_beam_params()` builds the case list, so the numbers live in
            # the helper, not in the test's own decorator.
            helper = (functions or {}).get(sub.func.id)
            if helper is not None:
                numbers |= _numbers_in(helper, consts, functions, depth - 1)
        value = eval_number(sub, consts)
        if value is not None:
            numbers.add(value)
    return numbers


def stated_numbers(
    func: ast.FunctionDef,
    consts: dict[str, Any],
    functions: dict[str, ast.FunctionDef] | None = None,
) -> set[float]:
    """Every number the function's parameters state, one helper call deep."""
    numbers: set[float] = set()
    for decorator in func.decorator_list:
        numbers |= _numbers_in(decorator, consts, functions)
    return numbers


def assertion_calls(func: ast.FunctionDef) -> list[ast.Call]:
    """Calls that assert, which is where a *validation* tolerance lives.

    `tol=` also appears on geometric helpers (`_find_node_by_xyz(..., tol=1e-4)`,
    `_twisted_beam_fixed(..., tol=1e-6)`), so scanning every keyword would import
    node-search tolerances into the matrix as if they were acceptance bounds.
    """
    asserted: set[int] = set()
    for sub in ast.walk(func):
        if isinstance(sub, ast.Assert):
            for inner in ast.walk(sub.test):
                asserted.add(id(inner))
    calls: list[ast.Call] = []
    for sub in ast.walk(func):
        if not isinstance(sub, ast.Call):
            continue
        if id(sub) in asserted:
            calls.append(sub)
            continue
        name = ast.unparse(sub.func).rsplit(".", 1)[-1].lower()
        if name.startswith(("assert", "verify", "check")):
            calls.append(sub)
    return calls


def parametrized_argnames(func: ast.FunctionDef) -> list[str]:
    """The parameter names in the order pytest writes them into the node id.

    Stacked decorators reach the id bottom-up (the decorator closest to the function
    first), which is the reverse of their source order, so the names are collected
    that way to line up positionally with the id's tokens.
    """
    names: list[str] = []
    for decorator in reversed(func.decorator_list):
        if isinstance(decorator, ast.Call) and "parametrize" in ast.unparse(decorator.func):
            if not decorator.args:
                continue
            raw = ast.unparse(decorator.args[0]).strip("'\"")
            names.extend(name.strip() for name in raw.split(",") if name.strip())
    return names


def param_value_by_name(
    func: ast.FunctionDef,
    consts: dict[str, Any],
    name: str,
    params: str | None,
) -> Any:
    """What a name holds for this node: a parametrise parameter or a loop variable.

    Case 3.5 asserts with `tol=tol` where `tol` is a parametrise parameter, so the
    value is only recoverable from the node's own id; the loop-over-a-constant case is
    kept as a fallback for parameters built inline.
    """
    argnames = parametrized_argnames(func)
    if name in argnames:
        tokens = param_tokens(params, len(argnames))
        index = argnames.index(name)
        if index < len(tokens):
            token = tokens[index]
            try:
                return float(token)
            except ValueError:
                return token
    wanted = param_numbers(params, len(argnames))
    for sub in ast.walk(func):
        if not isinstance(sub, ast.For):
            continue
        if isinstance(sub.target, ast.Tuple):
            targets = [item.id for item in sub.target.elts if isinstance(item, ast.Name)]
        elif isinstance(sub.target, ast.Name):
            targets = [sub.target.id]
        else:
            continue
        if name not in targets:
            continue
        rows = consts.get(sub.iter.id) if isinstance(sub.iter, ast.Name) else None
        if not isinstance(rows, list):
            continue
        position = targets.index(name)
        for row in rows:
            if not isinstance(row, tuple) or position >= len(row):
                continue
            numbers = set()
            for item in row:
                value = eval_number(ast.Constant(item), consts)
                if value is not None:
                    numbers.add(value)
            if wanted and wanted <= numbers:
                return row[position]
    return None


def tolerance_sites(
    func: ast.FunctionDef,
    consts: dict[str, Any],
    params: str | None,
) -> tuple[list[ToleranceSite], list[str]]:
    """The tolerances the function asserts on, plus the bounds it ignores.

    A zero bound (`atol=0.0`) asserts nothing, so it is reported as ignored rather
    than stored as if it bounded the comparison.
    """
    sites: list[ToleranceSite] = []
    ignored: list[str] = []
    calls = assertion_calls(func)
    for call in calls:
        for keyword in call.keywords:
            if keyword.arg not in {"tol", "rtol", "atol"}:
                continue
            value = resolve_literal(keyword.value, consts)
            if value is None and isinstance(keyword.value, ast.Name):
                value = param_value_by_name(func, consts, keyword.value.id, params)
            numeric = float(value) if isinstance(value, int | float) else None
            source = ast.unparse(call)[:120]
            if numeric == 0.0:
                ignored.append(f"{keyword.arg}=0 at line {call.lineno} ({source})")
                continue
            sites.append(
                ToleranceSite(
                    line=call.lineno,
                    # `tol` on this suite's helpers is the relative bound
                    # (`assert_relative_error` divides by the reference).
                    kind="rtol" if keyword.arg == "tol" else str(keyword.arg),
                    value=numeric,
                    source=source,
                    asserts=True,
                )
            )
    for sub in ast.walk(func):
        if isinstance(sub, ast.Assert) and isinstance(sub.test, ast.Compare):
            for op, comparator in zip(sub.test.ops, sub.test.comparators, strict=True):
                if not isinstance(op, ast.Lt | ast.LtE):
                    continue
                value = resolve_literal(comparator, consts)
                if not isinstance(value, int | float) or value == 0:
                    continue
                sites.append(
                    ToleranceSite(
                        line=sub.lineno,
                        kind="rel_err",
                        value=float(value),
                        source=ast.unparse(sub.test)[:120],
                        asserts=True,
                    )
                )
    unique: dict[tuple[int, str, float | None], ToleranceSite] = {}
    for site in sites:
        unique.setdefault((site.line, site.kind, site.value), site)
    return [unique[key] for key in sorted(unique, key=lambda item: item[0])], ignored


def parametrize_table(func: ast.FunctionDef, consts: dict[str, Any]) -> list[dict[str, Any]]:
    """Decorators, with literal case values where the code states them."""
    table: list[dict[str, Any]] = []
    for decorator in func.decorator_list:
        if not isinstance(decorator, ast.Call) or "parametrize" not in ast.unparse(decorator.func):
            continue
        argnames = ast.unparse(decorator.args[0]).strip("'\"") if decorator.args else ""
        cases: list[Any] = []
        if len(decorator.args) > 1:
            raw = decorator.args[1]
            literal = resolve_literal(raw, consts)
            if isinstance(literal, list):
                cases = [list(case) if isinstance(case, tuple) else case for case in literal]
        table.append(
            {
                "argnames": [name.strip() for name in argnames.split(",") if name.strip()],
                "cases": cases,
                "source": ast.unparse(decorator.args[1])[:80] if len(decorator.args) > 1 else "",
            }
        )
    return table


def eval_number(node: ast.AST, consts: dict[str, Any]) -> float | None:
    """The numeric value of a constant expression, or None.

    The suite states parameters as expressions (`1 / 100`, `1.0e2`), so
    `ast.literal_eval` alone cannot resolve them and a literal-only check would
    report every evaluated node-id value as a divergence.
    """
    if isinstance(node, ast.Constant):
        value = node.value
        if isinstance(value, int | float) and not isinstance(value, bool):
            return float(value)
        return None
    if isinstance(node, ast.UnaryOp) and isinstance(node.op, (ast.USub, ast.UAdd)):
        inner = eval_number(node.operand, consts)
        if inner is None:
            return None
        return -inner if isinstance(node.op, ast.USub) else inner
    if isinstance(node, ast.BinOp) and isinstance(
        node.op, (ast.Add, ast.Sub, ast.Mult, ast.Div, ast.Pow)
    ):
        left = eval_number(node.left, consts)
        right = eval_number(node.right, consts)
        if left is None or right is None:
            return None
        if isinstance(node.op, ast.Add):
            return left + right
        if isinstance(node.op, ast.Sub):
            return left - right
        if isinstance(node.op, ast.Mult):
            return left * right
        if isinstance(node.op, ast.Div):
            try:
                return left / right
            except ZeroDivisionError:
                return None
        try:
            return left**right
        except OverflowError:
            return None
    if isinstance(node, ast.Name) and isinstance(consts.get(node.id), int | float):
        return float(consts[node.id])
    return None


def cross_check_params(
    params: str | None,
    consts: dict[str, Any],
    func: ast.FunctionDef,
    functions: dict[str, ast.FunctionDef] | None = None,
) -> tuple[list[str], list[str], bool]:
    """E3: the numbers pytest put in the node id, checked against the code.

    The parametrisation id is a second witness for the values the code states: where
    the two disagree, one of them is wrong and the difference is reported, not
    silently reconciled. When the function states no readable numbers at all the
    check is reported as unavailable rather than as a divergence, because a
    divergence nobody can substantiate is a false finding.
    """
    if not params:
        return ([], [], True)
    numbers = stated_numbers(func, consts, functions)
    wanted = param_numbers(params, len(parametrized_argnames(func)))
    if not numbers:
        return (sorted(f"{value:g}" for value in wanted), [], False)
    agreed: list[str] = []
    diverged: list[str] = []
    for value in sorted(param_numbers(params, len(parametrized_argnames(func)))):
        match = any(abs(value - stated) <= 1e-9 * max(1.0, abs(stated)) for stated in numbers)
        (agreed if match else diverged).append(f"{value:g}")
    return (agreed, diverged, True)


def slugify(text: str, limit: int = 56) -> str:
    slug = re.sub(r"[^a-z0-9]+", "_", text.lower()).strip("_")
    return slug[:limit].strip("_") or "case"


def build_rows(
    group_id: str,
    slug: str,
    citation: str,
    scope: str,
    nodes: list[NodeInfo],
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """E1-E3 into rows: one row per collected node, none left unclaimed."""
    tree = ast.parse((REPO_ROOT / scope).read_text(encoding="utf-8"))
    consts = module_constants(tree)
    functions = {
        node.name: node for node in tree.body if isinstance(node, ast.FunctionDef)
    }
    report: dict[str, Any] = {"nodes": [], "unclaimed": [], "diverged": {}, "prints": {}}
    rows: list[dict[str, Any]] = []
    used: set[str] = set()
    for info in nodes:
        func = functions.get(info.function)
        if func is None:
            report["unclaimed"].append(info.node)
            continue
        docstring = ast.get_docstring(func) or ""
        validates = docstring.strip().splitlines()[0] if docstring.strip() else info.function
        sites, ignored = tolerance_sites(func, consts, info.params)
        if ignored:
            report.setdefault("ignored_bounds", {})[info.node] = ignored
        if not sites:
            report["unclaimed"].append(info.node)
            continue
        unresolved = [site.source for site in sites if site.value is None]
        if unresolved:
            # Never write a row whose tolerance cannot be read: an invalid row is
            # worse than a reported gap, because it looks like evidence.
            report.setdefault("unresolved", {})[info.node] = unresolved
            report["unclaimed"].append(info.node)
            continue
        row_id = f"{slug}.{slugify(info.function, 40)}.{slugify(info.params or 'single', 40)}"
        suffix = 2
        candidate = row_id
        while candidate in used:
            candidate = f"{row_id}_{suffix}"
            suffix += 1
        row_id = candidate
        used.add(row_id)

        agreed, diverged, checked = cross_check_params(info.params, consts, func, functions)
        if not checked:
            report.setdefault("cross_check_unavailable", []).append(info.node)
        elif diverged:
            report["diverged"][info.node] = diverged
        report.setdefault("parametrize", {})[info.function] = [
            table["argnames"] for table in parametrize_table(func, consts)
        ]
        prints = [
            ast.unparse(sub)[:100]
            for sub in ast.walk(func)
            if isinstance(sub, ast.Call)
            and isinstance(sub.func, ast.Name)
            and sub.func.id == "print"
            and re.search(
                r"expect|ref|err|margin|norm|ratio|delta|%", ast.unparse(sub), re.I
            )
        ]
        if prints:
            report["prints"][info.node] = prints

        comparisons = []
        for site in sites:
            kind = "paper" if re.search(r"paper|table|expected", site.source, re.I) else "analytical"
            reference: dict[str, Any] = {"kind": kind, "label": site.source}
            if kind == "paper":
                reference["citation"] = citation
            comparisons.append(
                {
                    "label": f"{site.kind} at line {site.line}",
                    "asserted": site.asserts,
                    "reference": reference,
                    "tolerance": {
                        "kind": site.kind,
                        "value": site.value,
                        "source": f"{scope}:{site.line}",
                        "justified": True,
                        "justification": (
                            "extracted from the code; whether the justification holds is "
                            "adjudicated in the pilot (T6)"
                        ),
                    },
                    "measured": {"status": "not_measured", "margin_pct": None},
                }
            )
        rows.append(
            {
                "id": row_id,
                "group": group_id,
                "title": info.node.split("::", 1)[1],
                "tests": [info.node],
                "validates": validates,
                "comparisons": comparisons,
                "flags": [],
                "notes": (
                    "Derived by `validation_matrix extract` from the code alone: the "
                    "human-facing description and the adjudicated reference land in the "
                    "pilot's T5/T6. No margin is measured yet."
                ),
            }
        )
        report["nodes"].append(info.node)
    return rows, report


def command_extract(args: argparse.Namespace) -> int:
    store = load_store(args.store)
    group = store.groups.get(args.group)
    if group is None:
        raise StoreError(f"unknown group: {args.group}")
    scope = args.scope or (group.get("source_files") or [None])[0]
    if not scope:
        raise StoreError(f"group {args.group} declares no source file; pass --scope")
    citation = args.citation or str(group.get("citation") or "")
    if not citation:
        raise StoreError(
            f"group {args.group} declares no 'citation'; pass --citation for paper references"
        )
    nodes = collect_nodes(str(scope))
    rows, report = build_rows(
        args.group,
        str(group.get("slug")),
        citation,
        str(scope),
        nodes,
    )
    payload = {
        "scope": str(scope),
        "collected": len(nodes),
        "rows": len(rows),
        "claimed": len(report["nodes"]),
        "unclaimed": report["unclaimed"],
        "diverged_params": report["diverged"],
        "unasserted_prints": report["prints"],
    }
    if args.json:
        print(json.dumps(payload, indent=2))
    else:
        print(f"scope      : {scope}")
        print(f"collected  : {len(nodes)}")
        print(f"rows       : {len(rows)}")
        print(f"claimed    : {len(report['nodes'])}")
        if report.get("ignored_bounds"):
            print(f"ignored bounds (assert nothing): {len(report['ignored_bounds'])} node(s)")
        if report.get("cross_check_unavailable"):
            print(
                "cross-check unavailable (no readable numbers in the code): "
                f"{len(report['cross_check_unavailable'])} node(s)"
            )
        for node, sources in (report.get("unresolved") or {}).items():
            print(f"UNRESOLVED tolerance (row not written): {node}")
            for source in sources:
                print(f"  - {source}")
        if report["unclaimed"]:
            print(f"UNCLAIMED  : {len(report['unclaimed'])}")
            for node in report["unclaimed"]:
                print(f"  - {node}")
        for node, tokens in report["diverged"].items():
            print(f"param id diverges from the code's literals: {node} -> {tokens}")
        for node, calls in report["prints"].items():
            print(f"printed and maybe never asserted (candidate info_only): {node}")
            for call in calls:
                print(f"  - {call}")
    if args.write:
        target = args.store / "rows" / f"{args.group}-{group.get('slug')}.yaml"
        target.parent.mkdir(parents=True, exist_ok=True)
        dump_yaml(target, {"group": args.group, "rows": rows})
        print(f"wrote {display_path(target)} ({len(rows)} rows)")
        probe = load_store(args.store)
        for finding in probe.errors():
            print(finding.as_text(), file=sys.stderr)
        if probe.errors():
            return EXIT_FINDINGS
    return EXIT_OK


# --------------------------------------------------------------------------- #
# Adjudications: the audit trail of a reconciled group
# --------------------------------------------------------------------------- #
#
# The log is data, not a document: `check` validates it, so an entry that references a row
# that no longer exists, or that states no verdict, is an error rather than a stale note.

ADJUDICATION_KEYS = {
    "id",
    "row",
    "field",
    "candidate_code",
    "candidate_md",
    "verdict",
    "resolution",
    "decided_by",
    "date",
    "status",
}
ADJUDICATION_STATUSES = {"open", "resolved", "normalized"}
ADJUDICATION_FILE_KEYS = {"group", "scope", "run", "date", "reconciliation", "adjudications"}
ADJUDICATION_ID_RE = re.compile(r"^ADJ-\d{4}$")


def load_adjudications(store: Store) -> None:
    directory = store.root / "adjudications"
    if not directory.exists():
        return
    known_rows = set(store.by_id())
    seen: set[str] = set()
    for path in sorted(directory.glob("*.yaml")):
        rel = display_path(path)
        data = load_yaml(path)
        if not isinstance(data, dict):
            store.error(rel, "expected a mapping")
            continue
        unknown = _unknown_keys(data, ADJUDICATION_FILE_KEYS)
        if unknown:
            store.error(rel, f"unknown keys: {', '.join(unknown)}")
        group_id = data.get("group")
        if not isinstance(group_id, str) or group_id not in store.groups:
            store.error(rel, f"unknown group: {group_id!r}")
        entries = data.get("adjudications")
        if entries is None:
            entries = []
        if not isinstance(entries, list):
            store.error(rel, "'adjudications' must be a list")
            continue
        for index, entry in enumerate(entries):
            where = f"{rel}[{index}]"
            if not isinstance(entry, dict):
                store.error(where, "expected a mapping")
                continue
            unknown = _unknown_keys(entry, ADJUDICATION_KEYS)
            if unknown:
                store.error(where, f"unknown keys: {', '.join(unknown)}")
            entry_id = entry.get("id")
            if not isinstance(entry_id, str) or not ADJUDICATION_ID_RE.match(entry_id):
                store.error(where, f"invalid id: {entry_id!r} (expected ADJ-NNNN)")
            elif entry_id in seen:
                store.error(where, f"duplicate adjudication id: {entry_id}")
            else:
                seen.add(entry_id)
                where = f"{rel}[{entry_id}]"
            row_id = entry.get("row")
            if row_id is not None and row_id not in known_rows:
                store.error(where, f"references a row that does not exist: {row_id!r}")
            for name in ("field", "verdict", "resolution", "decided_by"):
                if not isinstance(entry.get(name), str) or not entry[name]:
                    store.error(where, f"{name!r} must be a non-empty string")
            status = entry.get("status", "resolved")
            if status not in ADJUDICATION_STATUSES:
                store.error(
                    where, f"status must be one of {sorted(ADJUDICATION_STATUSES)}, got {status!r}"
                )
            date_value = entry.get("date")
            if date_value is not None and not is_iso_date(date_value):
                store.error(where, f"date must be YYYY-MM-DD, got {date_value!r}")


# --------------------------------------------------------------------------- #
# Gaps: what the suite does NOT validate
# --------------------------------------------------------------------------- #
#
# A gap is a claim about absence, so it cannot be derived from rows: a row that does not
# exist cannot be queried. The not-citable list is therefore data, and `check` validates it.

GAP_STATUSES = {"not_validated", "bounded", "open_defect"}
GAP_KEYS = {
    "id",
    "scope",
    "status",
    "reason",
    "evidence",
    "consequence",
    "citations_forbidden",
    "rows",
    "notes",
}
GAP_FILE_KEYS = {"version", "gaps"}


def load_residual_patterns(store: Store) -> None:
    """How each group's tests print their residuals, keyed by group id."""
    path = store.root / RESIDUAL_FILE
    if not path.exists():
        return
    rel = display_path(path)
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        store.error(rel, f"could not be read: {exc}")
        return
    patterns = data.get("patterns") if isinstance(data, dict) else None
    if not isinstance(patterns, dict):
        store.error(rel, "expected a mapping with a 'patterns' object")
        return
    for group_id, entry in patterns.items():
        where = f"{rel}[{group_id}]"
        if group_id not in store.groups:
            store.error(where, f"unknown group: {group_id!r}")
        if not isinstance(entry, dict):
            store.error(where, "expected a mapping")
            continue
        unknown = _unknown_keys(entry, RESIDUAL_KEYS)
        if unknown:
            store.error(where, f"unknown keys: {', '.join(unknown)}")
        if not isinstance(entry.get("asserted"), str):
            store.error(where, "an 'asserted' pattern is required")
        for name, pattern in entry.items():
            try:
                re.compile(str(pattern))
            except re.error as exc:
                store.error(where, f"{name} is not a valid regex: {exc}")
        store.residual_patterns[str(group_id)] = {key: str(value) for key, value in entry.items()}


def load_gaps(store: Store) -> None:
    path = store.root / "gaps.yaml"
    if not path.exists():
        return
    rel = display_path(path)
    data = load_yaml(path)
    if not isinstance(data, dict) or not isinstance(data.get("gaps"), list):
        store.error(rel, "expected a mapping with a 'gaps' list")
        return
    unknown = _unknown_keys(data, GAP_FILE_KEYS)
    if unknown:
        store.error(rel, f"unknown keys: {', '.join(unknown)}")
    known_rows = set(store.by_id())
    for index, entry in enumerate(data["gaps"]):
        where = f"{rel}[{index}]"
        if not isinstance(entry, dict):
            store.error(where, "expected a mapping")
            continue
        unknown = _unknown_keys(entry, GAP_KEYS)
        if unknown:
            store.error(where, f"unknown keys: {', '.join(unknown)}")
        gap_id = entry.get("id")
        if not isinstance(gap_id, str) or not SLUG_RE.match(gap_id):
            store.error(where, f"invalid id: {gap_id!r}")
        else:
            where = f"{rel}[{gap_id}]"
            if gap_id in store.gaps:
                store.error(where, "duplicate gap id")
            else:
                store.gaps[gap_id] = entry
        if entry.get("status") not in GAP_STATUSES:
            store.error(where, f"status must be one of {sorted(GAP_STATUSES)}")
        for name in ("scope", "reason", "evidence"):
            if not isinstance(entry.get(name), str) or not entry[name]:
                store.error(where, f"{name!r} must be a non-empty string")
        if not isinstance(entry.get("citations_forbidden"), bool):
            store.error(where, "citations_forbidden must be a boolean")
        for row_id in entry.get("rows") or []:
            if row_id not in known_rows:
                store.error(where, f"references a row that does not exist: {row_id!r}")


def command_gaps(args: argparse.Namespace) -> int:
    store = load_store(args.store)
    entries = list(store.gaps.values())
    if args.json:
        print(json.dumps(entries, indent=2))
    else:
        for entry in entries:
            forbidden = "NOT CITABLE" if entry.get("citations_forbidden") else "citable"
            print(f"{entry['id']:34s} {entry['status']:14s} {forbidden:11s} {entry['scope']}")
        print(f"gaps: {len(entries)}")
    return EXIT_OK


# --------------------------------------------------------------------------- #
# Regression: the store is the baseline
# --------------------------------------------------------------------------- #


def parse_prints(output: str) -> dict[str, list[str]]:
    """Map each collected node to the lines it printed.

    pytest prints a verbose node id with **no newline**, so the test's first print lands
    on that same line; the rest follow, and the status word comes on its own line after
    them. Parsing for `nodeid PASSED` on one line therefore matches nothing.
    """
    prints: dict[str, list[str]] = {}
    current: str | None = None
    for raw in output.splitlines():
        line = raw.rstrip()
        start = NODE_START_RE.match(line)
        if start:
            node = start.group("node")
            if not node:
                continue
            current = node
            prints.setdefault(node, [])
            rest = (start.group("rest") or "").strip()
            if rest:
                prints[node].append(rest)
            continue
        if NODE_STATUS_RE.match(line.strip()):
            current = None
            continue
        if current is not None and line.strip():
            prints[current].append(line.strip())
    return prints


def capture_prints(scope: str) -> dict[str, list[str]]:
    """Run the scope and read the prints out of its output."""
    completed = subprocess.run(
        [sys.executable, "-m", "pytest", "-o", "addopts=", "-s", "-v", scope],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
        check=False,
    )
    return parse_prints(completed.stdout)


def extract_residuals(
    lines: list[str], residual: dict[str, Any]
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    asserted: list[dict[str, str]] = []
    unasserted: list[dict[str, str]] = []
    for line in lines:
        match = re.search(str(residual["asserted"]), line)
        if match:
            asserted.append(match.groupdict())
            continue
        pattern = residual.get("unasserted")
        if pattern:
            match = re.search(str(pattern), line)
            if match:
                unasserted.append(match.groupdict())
    return asserted, unasserted


def _float(value: Any) -> float | None:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def compare_row(
    ref: RowRef,
    asserted: list[dict[str, str]],
    unasserted: list[dict[str, str]],
) -> list[dict[str, Any]]:
    """Pair each printed residual with the comparison it belongs to.

    The prints and the asserted comparisons are both in source order, so the Nth print
    is the Nth assertion. A count mismatch is reported as `unmapped` instead of guessed:
    attaching a margin to the wrong comparison would manufacture a baseline.
    """
    comparisons = [item for item in (ref.data.get("comparisons") or []) if isinstance(item, dict)]
    checked = [
        item for item in comparisons if isinstance(item.get("asserted"), bool) and item["asserted"]
    ]
    if len(checked) != len(asserted):
        return [
            {
                "row": ref.id,
                "verdict": "unmapped",
                "detail": (
                    f"{len(asserted)} printed residual(s) for {len(checked)} asserted "
                    "comparison(s)"
                ),
            }
        ]
    results: list[dict[str, Any]] = []
    for index, (comparison, printed) in enumerate(zip(checked, asserted, strict=True)):
        current = _float(printed.get("error"))
        baseline = _float((comparison.get("measured") or {}).get("margin_pct"))
        if baseline is None:
            verdict = "new_baseline"
        elif current is not None and abs(current - baseline) > 0.005:
            verdict = "drifted"
        else:
            verdict = "same"
        results.append(
            {
                "row": ref.id,
                "comparison": index,
                "label": comparison.get("label"),
                "verdict": verdict,
                "baseline": baseline,
                "current": current,
                "value": _float(printed.get("value")),
                "expected": printed.get("expected"),
                "text": f"{printed.get('value')} ({printed.get('error')}%)",
            }
        )
    for printed in unasserted:
        results.append(
            {
                "row": ref.id,
                "comparison": None,
                "label": "printed and never asserted",
                "verdict": "informational",
                "baseline": None,
                "current": _float(printed.get("error")),
                "value": _float(printed.get("value")),
                "expected": printed.get("expected"),
                "text": f"{printed.get('value')} ({printed.get('error')}%)",
            }
        )
    return results


def head_revision() -> str:
    completed = subprocess.run(
        ["git", "rev-parse", "--short", "HEAD"],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
        check=False,
    )
    return completed.stdout.strip() or "unknown"


def write_measurement(store: Store, ref: RowRef, index: int, result: dict[str, Any]) -> None:
    """Record one printed residual as the comparison's baseline."""
    path = ref.path
    if path is None:
        raise StoreError(f"row {ref.id} has no backing file to write")
    document = load_yaml(path)
    rows = document.get("rows") if isinstance(document, dict) else document
    if not isinstance(rows, list):
        raise StoreError(f"{ref.path} is not a row file this tool can write back")
    row = rows[ref.index]
    comparison = row["comparisons"][index]
    comparison["measured"] = {
        "status": "measured",
        "run": head_revision(),
        "date": date.today().isoformat(),
        "raw": result["value"],
        "margin_pct": result["current"],
        "text": result["text"],
    }
    if comparison.get("expected") is None and result.get("expected") is not None:
        comparison["expected"] = result["expected"]
    dump_yaml(path, document)


def command_regression(args: argparse.Namespace) -> int:
    store = load_store(args.store)
    group = store.groups.get(args.group)
    if group is None:
        raise StoreError(f"unknown group: {args.group}")
    residual = store.residual_patterns.get(str(args.group))
    if not residual:
        raise StoreError(
            f"group {args.group} declares no residual patterns, so its prints cannot be "
            f"read; declare them in docs/validation/{RESIDUAL_FILE}"
        )
    scope = args.scope or " ".join(str(item) for item in group.get("source_files") or [])
    prints = capture_prints(scope)
    if not prints:
        raise StoreError(f"no node output captured for {scope!r}; is the environment active?")

    results: list[dict[str, Any]] = []
    claimed: set[str] = set()
    for ref in store.rows:
        nodes = [node for node in (ref.data.get("tests") or []) if node in prints]
        if not nodes:
            continue
        claimed.update(nodes)
        lines = [line for node in nodes for line in prints[node]]
        asserted, unasserted = extract_residuals(lines, residual)
        results.extend(compare_row(ref, asserted, unasserted))
    for node in sorted(set(prints) - claimed):
        results.append({"row": None, "verdict": "unclaimed", "detail": node})

    if args.write:
        for result in results:
            index = result.get("comparison")
            if result["verdict"] not in {"same", "drifted", "new_baseline"}:
                continue
            if not isinstance(index, int):
                continue
            write_measurement(store, store.by_id()[result["row"]], index, result)
        print(f"wrote {len(results)} measurement(s)")

    drifted = [item for item in results if item["verdict"] in {"drifted", "unmapped", "unclaimed"}]
    if args.json:
        print(json.dumps(results, indent=2))
    else:
        for item in results:
            if item["verdict"] == "same":
                continue
            where = item.get("row") or item.get("detail")
            if item["verdict"] == "drifted":
                print(
                    f"DRIFTED  {where} [{item['label']}]: was {item['baseline']}%, "
                    f"now {item['current']}%"
                )
            elif item["verdict"] == "new_baseline":
                print(f"NEW      {where} [{item['label']}]: {item['current']}% ({item['text']})")
            elif item["verdict"] == "informational":
                print(f"INFO     {where}: {item['text']} (never asserted)")
            else:
                print(f"{item['verdict'].upper():8s} {where} {item.get('detail') or ''}")
        counts: dict[str, int] = {}
        for item in results:
            counts[item["verdict"]] = counts.get(item["verdict"], 0) + 1
        summary = ", ".join(f"{count} {name}" for name, count in sorted(counts.items()))
        print(f"regression: {summary}")
    return EXIT_FINDINGS if drifted else EXIT_OK


# --------------------------------------------------------------------------- #
# Cross-validating the Markdown view against the store
# --------------------------------------------------------------------------- #

MD_TABLE_MIN_CELLS = 5


def md_function(cell: str) -> str:
    """The test function a Markdown row names, or '...' for a continuation row."""
    return cell.strip().strip("`").split("[", 1)[0].strip()


def find_md_table(md_text: str, source_file: str) -> list[dict[str, Any]]:
    """The Markdown table of the section that documents this source file.

    The section is found by its own heading, so no path or anchor is hardcoded: the
    group names its source file and the heading names the same file.
    """
    lines = md_text.splitlines()
    start = next(
        (index for index, line in enumerate(lines) if line.startswith("## ") and source_file in line),
        None,
    )
    if start is None:
        raise StoreError(f"no Markdown section documents {source_file}")
    end = next(
        (index for index in range(start + 1, len(lines)) if lines[index].startswith("## ")),
        len(lines),
    )
    rows: list[dict[str, Any]] = []
    function = ""
    for line in lines[start:end]:
        if not line.startswith("|"):
            continue
        cells = [cell.strip() for cell in line.strip().strip("|").split("|")]
        if len(cells) < MD_TABLE_MIN_CELLS or cells[0] == "test" or set(cells[0]) <= {"-", " "}:
            continue
        named = md_function(cells[0])
        if named and not named.startswith("..."):
            function = named
        rows.append(
            {
                "test": cells[0],
                "function": function,
                "cases": md_case_count(cells[0]),
                "validates": cells[1],
                "reference": cells[2],
                "tolerance": cells[3],
                "margin": cells[4],
                "notes": cells[5] if len(cells) > 5 else "",
            }
        )
    return rows


def md_case_count(cell: str) -> int:
    """How many cases a Markdown row covers, in either phrasing the document uses.

    Most rows say `(3 cases)`; the circular-plate rows encode it inside the test cell as
    `[clamped, 3 t/L]`. This is the number that reconciles the two counts: a Markdown row
    covers a *comparison*, while a collected node covers a *test*, and a test asserting two
    tables is one node and two comparison slots.
    """
    for pattern in (r"\((\d+)\s+cases?\)", r"\b(\d+)\s+t/L\b"):
        match = re.search(pattern, cell)
        if match:
            return int(match.group(1))
    return 1


def md_tolerances(text: str) -> set[tuple[str, float]]:
    """The tolerance a Markdown cell states, in the store's own vocabulary.

    `tol=` on this suite's helpers is a relative bound, so it maps to `rtol` exactly as the
    extractor maps it; without that, every row would look like a conflict.
    """
    found: set[tuple[str, float]] = set()
    for match in re.finditer(r"(rtol|tol|atol|rel_err)\s*(?:<|=)\s*([0-9.]+)", text):
        kind = "rtol" if match.group(1) == "tol" else match.group(1)
        found.add((kind, float(match.group(2))))
    return found


def md_margins(text: str) -> set[float]:
    return {float(value) for value in re.findall(r"\(([0-9.]+)%\)", text)}


def store_function_expectations(refs: list[RowRef]) -> tuple[set[tuple[str, float]], set[float]]:
    tolerances: set[tuple[str, float]] = set()
    margins: set[float] = set()
    for ref in refs:
        for comparison in ref.data.get("comparisons") or []:
            if not isinstance(comparison, dict):
                continue
            tolerance = comparison.get("tolerance") or {}
            value = _number(tolerance.get("value"))
            if value is not None and tolerance.get("kind"):
                tolerances.add((str(tolerance["kind"]), value))
            margin = _number((comparison.get("measured") or {}).get("margin_pct"))
            if margin is not None:
                margins.add(round(margin, 2))
    return tolerances, margins


def command_diff_against_md(args: argparse.Namespace) -> int:
    store = load_store(args.store)
    group = store.groups.get(args.group)
    if group is None:
        raise StoreError(f"unknown group: {args.group}")
    source = str((group.get("source_files") or [""])[0])
    md_path = args.md or (REPO_ROOT / "docs" / "validation-matrix.md")
    if not md_path.exists():
        raise StoreError(f"{display_path(md_path)} is gone; there is no view to diff against")
    md_rows = find_md_table(md_path.read_text(encoding="utf-8"), source)

    md_by_function: dict[str, list[dict[str, Any]]] = {}
    for row in md_rows:
        md_by_function.setdefault(row["function"], []).append(row)
    store_by_function: dict[str, list[RowRef]] = {}
    for ref in store.rows:
        if ref.data.get("group") != args.group:
            continue
        store_by_function.setdefault(md_function(str(ref.data.get("title"))), []).append(ref)

    report: list[dict[str, Any]] = []
    for function in sorted(set(md_by_function) | set(store_by_function)):
        md_group = md_by_function.get(function, [])
        refs = store_by_function.get(function, [])
        if not refs:
            report.append({"function": function, "verdict": "only_in_md", "md_rows": len(md_group)})
            continue
        if not md_group:
            report.append(
                {"function": function, "verdict": "only_in_code", "nodes": len(refs)}
            )
            continue
        md_tol: set[tuple[str, float]] = set()
        for row in md_group:
            md_tol |= md_tolerances(row["tolerance"])
        md_margin: set[float] = set()
        for row in md_group:
            md_margin |= md_margins(row["margin"])
        store_tol, store_margin = store_function_expectations(refs)
        verdict = "match"
        detail: list[str] = []
        if md_tol != store_tol:
            verdict = "value_conflict"
            detail.append(
                f"tolerance: md {sorted(md_tol)} vs store {sorted(store_tol)}"
            )
        missing = {value for value in md_margin if value not in store_margin}
        if missing:
            verdict = "value_conflict"
            detail.append(
                f"margins in the Markdown and not in the store: {sorted(missing)} "
                f"(store: {sorted(store_margin)})"
            )
        report.append(
            {
                "function": function,
                "verdict": verdict,
                "md_rows": len(md_group),
                "nodes": len(refs),
                "detail": detail,
            }
        )

    conflicts = [item for item in report if item["verdict"] != "match"]
    if args.json:
        print(json.dumps(report, indent=2))
    else:
        print(f"{'function':52s} {'md':>3s} {'nodes':>5s}  verdict")
        for item in report:
            print(
                f"{item['function'][:52]:52s} {item.get('md_rows', 0):3d} "
                f"{item.get('nodes', 0):5d}  {item['verdict']}"
            )
            for line in item.get("detail") or []:
                print(f"      {line}")
        md_total = len(md_rows)
        md_slots = sum(row["cases"] for row in md_rows)
        nodes_total = sum(len(refs) for refs in store_by_function.values())
        comparisons = sum(
            len([item for item in (ref.data.get("comparisons") or []) if isinstance(item, dict)])
            for refs in store_by_function.values()
            for ref in refs
        )
        print(
            f"\ndiff: {md_total} Markdown row(s) covering {md_slots} case slot(s) against "
            f"{nodes_total} collected node(s); {len(conflicts)} function(s) need adjudication"
        )
        if md_slots != comparisons:
            print(
                f"NOTE: the Markdown's case slots ({md_slots}) do not equal the store's "
                f"comparison count ({comparisons}); a row and a comparison are the same unit, "
                "so this is a real discrepancy"
            )
    return EXIT_FINDINGS if conflicts else EXIT_OK


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
        "action",
        choices=["list", "gaps", "get", "check", "bibtex", "where-used", "render"],
    )
    bibliography.add_argument("key", nargs="?", help="reference key for get and where-used")
    bibliography.add_argument("--doi-status", choices=sorted(DOI_STATUSES))
    bibliography.add_argument("--kind", choices=sorted(BIB_KINDS))
    bibliography.add_argument("--held", action="store_true", help="only entries whose PDF is held")
    bibliography.add_argument("--not-held", action="store_true")
    bibliography.add_argument("--out", type=Path, help="write bibtex or the rendered file here")
    bibliography.add_argument(
        "--check", action="store_true", help="with render: fail instead of writing"
    )
    bibliography.add_argument("--json", action="store_true")
    bibliography.set_defaults(func=command_references)

    setter = subparsers.add_parser("set", help="edit one row in place")
    setter.add_argument("id")
    setter.add_argument("assign", nargs="*", metavar="PATH=VALUE")
    setter.add_argument("--unset", action="append", default=[], metavar="PATH")
    setter.add_argument("--add-flag", action="append", default=[], metavar="FLAG")
    setter.add_argument("--remove-flag", action="append", default=[], metavar="FLAG")
    setter.set_defaults(func=command_set)

    extract = subparsers.add_parser(
        "extract", help="derive rows from the collected nodes and the test source"
    )
    extract.add_argument("--group", default="3")
    extract.add_argument("--scope", help="test file to read (default: the group's source)")
    extract.add_argument("--citation", help="bibliography key for paper references")
    extract.add_argument("--write", action="store_true", help="write the rows file")
    extract.add_argument("--json", action="store_true")
    extract.set_defaults(func=command_extract)

    regression = subparsers.add_parser(
        "regression", help="re-run the scope and diff its printed residuals against the store"
    )
    regression.add_argument("--group", default="3")
    regression.add_argument("--scope", help="what to run (default: the group's source files)")
    regression.add_argument(
        "--write", action="store_true", help="record the printed residuals as the baseline"
    )
    regression.add_argument("--json", action="store_true")
    regression.set_defaults(func=command_regression)

    cross_check = subparsers.add_parser(
        "diff-against-md", help="cross-validate the Markdown view against the store"
    )
    cross_check.add_argument("--group", default="3")
    cross_check.add_argument(
        "--md", type=Path, help="the deleted Markdown view to diff against"
    )
    cross_check.add_argument("--json", action="store_true")
    cross_check.set_defaults(func=command_diff_against_md)

    gaps = subparsers.add_parser("gaps", help="what the suite does not validate")
    gaps.add_argument("--json", action="store_true")
    gaps.set_defaults(func=command_gaps)

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
