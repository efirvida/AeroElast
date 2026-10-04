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
    groups.yaml                 group registry (the directory is the classification)
    flags.yaml                  flag registry (judgement flags + derived flags)
    rows/<group-slug>.yaml      the rows of one group

`docs/validation-matrix.md` and `docs/references.md` were the store's generated views, and
they are gone: the store is the only artifact. `references render` survives as an export
command.

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
pipeline (`extract`, `render`).

Validation is hand-rolled on purpose: `jsonschema` is not importable in the pinned
environment, and CONTRIBUTING rule 6 makes this tool part of the suite's discipline
rather than an optional extra. The JSON Schemas under `docs/validation/schemas/` document the same
rules for editors and reviewers; `check` is the enforcement.

Exit codes: 0 clean, 1 findings, 2 usage or I/O error.
"""

from __future__ import annotations

import argparse
import ast
import contextlib
import hashlib
import io
from collections import Counter
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

# `ID_RE` and `SLUG_RE` have to agree, because a row id IS `<slug>.<test>.<params>`: whatever
# a group may call itself, the id grammar has to accept as the id's first segment. They did not
# agree, and a group called `orthotropic_shell_parity` produced three ids that no `set` could
# repair -- every write was refused for editing an invalid row.
ID_RE = re.compile(r"^[a-z][a-z0-9_]*(\.[a-z0-9_]+)+$")
SLUG_RE = re.compile(r"^[a-z][a-z0-9_]*$")
# A node id is `file::test[param]` for a module function and `file::Class::test[param]` for
# a method, so the class chain is part of the shape. Accepting only one segment after the file
# made every class-based test invisible to the extractor, which then reported "no nodes" for the
# file: 15 of the 42 validation files are written that way.
NODE_ID_RE = re.compile(
    r"^tests/[A-Za-z0-9_./-]+\.py::(?:[A-Za-z0-9_]+::)*[A-Za-z0-9_]+(?:\[[^\]]*\])?$"
)
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
    "cited_by_docs",
    "recovered_from",
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
MEASURED_KEYS = {"status", "run", "date", "raw", "margin_pct", "text", "digest"}
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
    # What the group's comparisons are against, when every one of them is against the same
    # thing. A mixed group declares none and each comparison is decided by hand: the extractor
    # reports a declared kind and never infers one, because a keyword guess is not evidence.
    "reference_kind",
    # Tests in this group's scope that are not validation rows, each with its reason. It lives
    # here, next to the scope it excludes from, and not in `gaps.yaml`: a gap is a claim about
    # evidence the suite does not have, while this is a classification of tests.
    "non_validation_tests",
    # Same-module helpers whose assertions ARE the comparison the test makes. Declared, because
    # structure cannot tell them from the machinery helpers beside them. See `followed_helpers`.
    "validation_helpers",
    # Assertions that are not comparisons against an independent reference: our own torque ruler, a
    # symmetry of our own model, a test's setup constant. Policy rules 1 and 6 leave them no place
    # in a row -- a row is evidence, and a reference must be independent -- so they are named here
    # by their `tolerance.source` and left out. Matched exactly, with stale reporting, like every
    # other declaration in this store: a line that moves makes the declaration visible, not silent.
    "non_reference_asserts",
}


def followed_helpers(
    func: ast.FunctionDef,
    functions: dict[str, ast.FunctionDef],
    declared: set[str],
) -> list[ast.FunctionDef]:
    """The declared same-module helpers the test reaches, and only those.

    A validation whose assertions live in a local helper is otherwise invisible: the tube torsion
    file has one, beside a sibling test whose own assert is machinery. Following every same-module
    helper is not the answer either -- group 3's tests reach `_assemble_global` and
    `_pinched_cylinder_load`, whose asserts check the load that was just applied, so the rows would
    fill with comparisons against nothing. Structure cannot tell the two apart, so the group declares
    which helpers hold reference comparisons, and this follows the declaration, transitively but only
    through declared names.
    """
    found: list[ast.FunctionDef] = []
    seen: set[str] = set()
    queue = [func]
    while queue:
        current = queue.pop(0)
        for node in ast.walk(current):
            if not isinstance(node, ast.Call) or not isinstance(node.func, ast.Name):
                continue
            name = node.func.id
            if name in seen or name not in declared:
                continue
            target = functions.get(name)
            if target is None or target is func:
                continue
            seen.add(name)
            found.append(target)
            queue.append(target)
    return found


def match_non_validation(group: dict[str, Any], node: str) -> tuple[str, str] | None:
    """The (pattern, reason) declaring `node` not a validation row, or None.

    A declaration matches three ways, and the third exists so that a declaration does not have to
    spell a node id that runs past the line budget with no way to wrap it:

      - the node id exactly;
      - a `::`-boundary prefix, so a whole class is declared once instead of twenty times;
      - a bare test name, matched after the last `::` of the node, so one test is named in twenty
        characters instead of a hundred.

    The last one can match more than intended only if two tests in the same group share a name and
    the declaration is not specific; the report prints the node each declaration matched, beside the
    reason, so a match that reached further than the author meant is visible rather than silent.
    """
    for entry in group.get("non_validation_tests") or []:
        if not isinstance(entry, dict):
            continue
        for pattern in entry.get("tests") or []:
            if not isinstance(pattern, str):
                continue
            # `test_x` and `test_x[Quad4]` are the same test, so the parametrised id is compared
            # by its name alone.
            bare = pattern.rsplit("::", 1)[-1].split("[", 1)[0]
            if (
                node == pattern
                or node.startswith(pattern.rstrip(":") + "::")
                or node.rsplit("::", 1)[-1].split("[", 1)[0] == bare
            ):
                return pattern, str(entry.get("reason") or "")
    return None


def classify_unclaimed(
    group: dict[str, Any], nodes: list[str]
) -> tuple[list[str], list[str], list[str]]:
    """Split unclaimed nodes into declared and undeclared, and name the stale declarations.

    A declaration that matches nothing is reported as well. It reads like coverage and is not,
    and a test that was renamed or deleted would otherwise leave a reason behind that nobody
    ever rechecks.
    """
    declared: list[str] = []
    undeclared: list[str] = []
    used: set[str] = set()
    for node in nodes:
        match = match_non_validation(group, node)
        if match is None:
            undeclared.append(node)
        else:
            declared.append(node)
            used.add(match[0])
    all_patterns = {
        pattern
        for entry in group.get("non_validation_tests") or []
        if isinstance(entry, dict)
        for pattern in entry.get("tests") or []
        if isinstance(pattern, str)
    }
    return declared, undeclared, sorted(all_patterns - used)
# A group declares how its tests print their residuals, because the suite prints prose, not
# a format. `asserted` matches the line behind a real assertion, `unasserted` the line the
# test prints and never asserts.
# The patterns themselves live in `residual-patterns.json`, keyed by group id: a regex
# is machine data, and the YAML prose linter enforces a line budget a long pattern cannot
# satisfy. `asserted` matches the line behind a real assertion, `unasserted` the line the
# test prints and never asserts.
RESIDUAL_KEYS = {"asserted", "unasserted"}
RESIDUAL_FILE = "residual-patterns.json"
# Digests of each group's source files, in JSON for the same reason: machine data, and the
# registries carry comments a YAML round-trip would drop.
SOURCES_FILE = "sources.json"
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
                "digest": LEAF,
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
    # Files under tests/validation/ that are not validation files, each with its reason.
    out_of_scope: list[dict[str, Any]] = field(default_factory=list)
    groups: dict[str, dict[str, Any]] = field(default_factory=dict)
    flag_entries: dict[str, dict[str, Any]] = field(default_factory=dict)
    derived_flags: set[str] = field(default_factory=set)
    references: dict[str, dict[str, Any]] = field(default_factory=dict)
    reference_header: str | None = None
    reference_sections: dict[str, dict[str, Any]] = field(default_factory=dict)
    residual_patterns: dict[str, dict[str, str]] = field(default_factory=dict)
    source_digests: dict[str, str] = field(default_factory=dict)
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
    # A validation file with no canonical comparison is not a row: rule 6 leaves a test that
    # compares nothing against an independent reference nowhere to sit. Declaring the file here is
    # what separates "out of scope" from "forgotten", and `status` counts only the undeclared ones
    # -- which the contract's own page promised before the tool could do it.
    out_of_scope = data.get("out_of_scope")
    if out_of_scope is not None and not isinstance(out_of_scope, list):
        store.error("groups.yaml", "expected 'out_of_scope' to be a list")
    for index, entry in enumerate(out_of_scope if isinstance(out_of_scope, list) else []):
        where = f"groups.yaml.out_of_scope[{index}]"
        if not isinstance(entry, dict):
            store.error(where, "expected a mapping")
            continue
        unknown = _unknown_keys(entry, {"files", "reason"})
        if unknown:
            store.error(where, f"unknown keys: {', '.join(unknown)}")
        files = entry.get("files")
        if (
            not isinstance(files, list)
            or not files
            or not all(isinstance(item, str) and item for item in files)
        ):
            store.error(where, "'files' must be a non-empty list of test paths")
            continue
        for item in files:
            if not (REPO_ROOT / item).exists():
                store.error(where, f"file does not exist: {item}")
        if not isinstance(entry.get("reason"), str) or not entry["reason"].strip():
            store.error(where, "'reason' must say why these are not validation files")
        store.out_of_scope.append(entry)

    for index, entry in enumerate(data["groups"]):
        where = f"groups.yaml[{index}]"
        if not isinstance(entry, dict):
            store.error(where, "expected a mapping")
            continue
        unknown = _unknown_keys(entry, GROUP_KEYS)
        if unknown:
            store.error(where, f"unknown keys: {', '.join(unknown)}")
        reference_kind = entry.get("reference_kind")
        if reference_kind is not None and reference_kind not in REFERENCE_KINDS:
            store.error(
                where,
                f"reference_kind must be one of {sorted(REFERENCE_KINDS)}, "
                f"got {reference_kind!r}",
            )
        assertions = entry.get("non_reference_asserts")
        if assertions is not None and not isinstance(assertions, list):
            store.error(where, "non_reference_asserts must be a list")
        for index, item in enumerate(assertions if isinstance(assertions, list) else []):
            at = f"{where}.non_reference_asserts[{index}]"
            if not isinstance(item, dict):
                store.error(at, "must be a mapping")
                continue
            unknown = _unknown_keys(item, {"at", "reason"})
            if unknown:
                store.error(at, f"unknown keys: {', '.join(unknown)}")
            source = item.get("at")
            if not isinstance(source, str) or not CITATION_SITE_RE.match(source):
                store.error(at, "'at' must be a <path>:<line> source, as the extractor writes it")
            if not isinstance(item.get("reason"), str) or not item["reason"].strip():
                store.error(at, "'reason' must say why this assertion is not a reference")
        helpers = entry.get("validation_helpers")
        if helpers is not None and (
            not isinstance(helpers, list)
            or not helpers
            or not all(isinstance(name, str) and name for name in helpers)
        ):
            store.error(where, "validation_helpers must be a non-empty list of function names")
        exclusions = entry.get("non_validation_tests")
        if exclusions is not None and not isinstance(exclusions, list):
            store.error(where, "non_validation_tests must be a list")
        for index, item in enumerate(exclusions if isinstance(exclusions, list) else []):
            at = f"{where}.non_validation_tests[{index}]"
            if not isinstance(item, dict):
                store.error(at, "must be a mapping")
                continue
            unknown = _unknown_keys(item, {"tests", "reason"})
            if unknown:
                store.error(at, f"unknown keys: {', '.join(unknown)}")
            patterns = item.get("tests")
            if (
                not isinstance(patterns, list)
                or not patterns
                or not all(isinstance(p, str) and p for p in patterns)
            ):
                store.error(at, "'tests' must be a non-empty list of node ids or prefixes")
            reason = item.get("reason")
            if not isinstance(reason, str) or not reason.strip():
                store.error(at, "'reason' must say why these tests are not validation rows")
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



def _validate_reference(store: Store, where: str, ref: Any) -> None:
    if not isinstance(ref, dict):
        store.error(where, "'reference' must be a mapping")
        return
    unknown = _unknown_keys(ref, REFERENCE_KEYS)
    if unknown:
        store.error(where, f"reference: unknown keys: {', '.join(unknown)}")
    kind = ref.get("kind")
    if kind is None:
        store.error(
            where,
            "reference.kind is not declared: the group declares what its comparisons are "
            "against, and the extractor never guesses one from the test's prose",
        )
    elif kind not in REFERENCE_KINDS:
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


# A ratio tolerance is a fraction. Nothing said so, so a percentage written where a fraction
# belongs - `5` for five percent - passed every check and then made `slack` and the `near` and
# `gt5` flags wrong by a factor of a hundred, silently, in the direction that looks like "far from
# the bound". The stored data is a fraction everywhere today, and now it has to stay one.
RATIO_TOLERANCE_MAX = 1.0


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
    value = tol.get("value") if isinstance(tol, dict) else None
    if kind in RATIO_TOLERANCE_KINDS and isinstance(value, int | float):
        if value > RATIO_TOLERANCE_MAX:
            store.error(
                where,
                f"tolerance.value is {value} for a {kind} tolerance, but a ratio tolerance is a "
                "fraction: five percent is 0.05, not 5. A percentage here would make slack and the "
                "near and gt5 flags wrong by a factor of a hundred",
            )
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
    # The class of the bug a reader found by eye: a reference side that is one number is stored
    # as a number, so that numeric queries on it work and it stays symmetric with measured.raw.
    # A string that parses as a number is exactly that mistake, and nothing else looks like it:
    # a cell list or a mode table does not parse.
    expected = comparison.get("expected")
    if isinstance(expected, str):
        try:
            float(expected)
        except ValueError:
            pass
        else:
            store.error(
                where,
                f"expected is the numeric string {expected!r}, but a reference side that is one "
                "number is stored as a number",
            )
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
            elif not node.startswith("tests/validation/"):
                store.error(
                    where,
                    f"a row may only claim a test under tests/validation/, and {node!r} is "
                    "not: a test that validates no physical quantity belongs in "
                    "tests/software/ and needs no row",
                )
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
    load_source_digests(store)
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


def row_errors(store: Store, path: Path, row: Any) -> set[str]:
    """What one row fails on by itself, so a write can refuse only to make it worse.

    `set` used to refuse any row that had an error, which made a row with two undeclared reference
    kinds unrepairable: fixing the first left the second undeclared, and every write was refused.
    The store-level `check` still reports whatever is left; this verb only declines to add to it.
    """
    scratch = Store(root=store.root, groups=store.groups)
    scratch.flag_entries = store.flag_entries
    scratch.derived_flags = store.derived_flags
    validate_row(scratch, str(path), row)
    return {finding.as_text() for finding in scratch.errors()}


def command_set(args: argparse.Namespace) -> int:
    store = load_store(args.store)
    ref = store.by_id().get(args.id)
    if ref is None:
        raise StoreError(f"unknown row id: {args.id}")
    if ref.path is None:
        raise StoreError(f"row {args.id} has no backing file to write")

    if store.errors():
        # A store that is broken is not one to pile edits onto, but a null `reference.kind` is
        # exactly the error `set` exists to repair: the extractor writes a row whose kind its group
        # did not declare, and declaring it is this verb's job. So the question is not whether the
        # store has errors but whether they are inside the reach of this verb. An error in a row is
        # repairable here; one in `groups.yaml` or `references.yaml` is not, and writing anyway
        # would be editing on top of a registry that does not hold.
        rows_root = display_path(args.store / "rows") + "/"
        foreign = [f for f in store.errors() if not f.where.startswith(rows_root)]
        if foreign:
            for finding in foreign:
                print(finding.as_text(), file=sys.stderr)
            raise StoreError(
                "refusing to write: the store has errors outside the row files, which this "
                "verb cannot repair"
            )

    document = load_yaml(ref.path)
    rows = document.get("rows") if isinstance(document, dict) else document
    if not isinstance(rows, list) or ref.index >= len(rows):
        raise StoreError(f"{ref.path} is not a row file this tool can write back")
    row = rows[ref.index]
    # What the row already failed on, taken before the edit is touched. The returned set is text,
    # so later mutations of `row` cannot reach it.
    before_errors = row_errors(store, ref.path, row)
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

    # The guard asks whether this write makes the row worse, not whether the row is already
    # imperfect. A row with two undeclared reference kinds cannot be repaired one kind at a time
    # otherwise: fixing the first leaves the second undeclared and every write is refused, which is
    # the same trap that stopped `set` from repairing a null kind at all.
    introduced = row_errors(store, ref.path, row) - before_errors
    if introduced:
        for finding in sorted(introduced):
            print(finding, file=sys.stderr)
        raise StoreError(
            "refusing to write: the edit would leave an error the row did not have"
        )

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
    # Two relations that rot the same way as a code site, and are checked the same way:
    # a document that cites the work, and the site the bibliographic data was recovered
    # from. Kept out of `notes`, where nothing validates them and a deleted document leaves
    # a sentence pointing nowhere.
    for name in ("cited_by_docs", "recovered_from"):
        for site in entry.get(name) or []:
            if not isinstance(site, str) or not CITATION_SITE_RE.match(site):
                store.error(where, f"{name}: not a file or file:line site: {site!r}")
                continue
            if not (REPO_ROOT / site.split(":", 1)[0]).exists():
                store.error(where, f"{name}: site does not exist: {site}")

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
    if entry.get("cited_by_docs"):
        sites = ", ".join(f"`{site}`" for site in entry["cited_by_docs"])
        paragraphs.append(f"*Cited by the documentation: {sites}.*")
    if entry.get("recovered_from"):
        sites = ", ".join(f"`{site}`" for site in entry["recovered_from"])
        paragraphs.append(f"*Recovered from: {sites}.*")
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
    # The class chain pytest prints between the file and the test, for a method; None for a
    # module-level function.
    owner: str | None = None

    @property
    def qualname(self) -> str:
        """The test as pytest names it, with its class chain: `Class::test`."""
        return f"{self.owner}::{self.function}" if self.owner else self.function


@dataclass
class ToleranceSite:
    line: int
    kind: str
    value: float | None
    source: str
    asserts: bool
    # The reference side of the comparison: `case.expected_normalized`, or the literal bound of
    # a bare assert. It names the reference, where the whole call only names the assertion.
    reference_expr: str | None = None
    # A canonical call (tests/support/assertions.py) states both of these, and then the store takes
    # them from the call instead of from the group: the comparison says what it is against, which is
    # the whole point of the contract, and a group-level default is only a fallback for the files
    # written before it.
    reference_name: str | None = None
    reference_kind: str | None = None


def parse_node_id(line: str) -> NodeInfo | None:
    """One collected node id as its parts, or None when the line is not a node id.

    The test is the last `::` segment and everything before it is the class chain, because that
    is the order pytest prints. Reading only the first segment after the file dropped every
    class-based test on the floor, and the extractor blamed the environment for it.
    """
    line = line.strip()
    if "::" not in line or not NODE_ID_RE.match(line):
        return None
    file, _, rest = line.partition("::")
    path, _, params = rest.partition("[")
    *owners, function = path.split("::")
    return NodeInfo(
        node=line,
        file=file,
        function=function,
        params=params.rstrip("]") if params else None,
        owner="::".join(owners) or None,
    )


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
        info = parse_node_id(line)
        if info is not None:
            nodes.append(info)
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
    extra: list[ast.FunctionDef] | None = None,
) -> tuple[list[ToleranceSite], list[str]]:
    """The tolerances the function asserts on, plus the bounds it ignores.

    Only a canonical call is a site: one that names its reference and its kind and
    carries the bound in `tol`, `rtol` or `atol`. A bare relational assert
    (`assert err < 0.05`) says a value came in below a bound and never says what it
    was measured against, so it is not a comparison the store can attribute; such a
    comparison is converted to the canonical call, not read. A zero bound
    (`atol=0.0`, or `x > 0.0`) still asserts nothing, so it is reported as ignored
    rather than stored as if it bounded anything.
    """
    sites: list[ToleranceSite] = []
    ignored: list[str] = []
    # The test first, then the declared helpers it reaches, in that order: appending after the
    # test's own sites keeps every existing row byte-identical when a group declares no helper.
    scopes = [func, *(extra or [])]
    calls = [call for scope in scopes for call in assertion_calls(scope)]
    for call in calls:
        stated_name: str | None = None
        stated_kind: str | None = None
        for keyword in call.keywords:
            if keyword.arg == "reference_name" and isinstance(keyword.value, ast.Constant):
                stated_name = str(keyword.value.value)
            elif keyword.arg == "kind" and isinstance(keyword.value, ast.Constant):
                stated_kind = str(keyword.value.value)
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
                    reference_expr=(
                        ast.unparse(call.args[1])[:80] if len(call.args) > 1 else None
                    ),
                    reference_name=stated_name,
                    reference_kind=stated_kind,
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


def name_without_prefix(name: str) -> str:
    """The test's own name, minus pytest's prefix, verbatim.

    Verbatim on purpose. A name is not a description, and composing one would mean the extractor
    inventing a summary instead of recording what the code states. Replacing the underscores with
    spaces read well for `test_axial_extension_matches_ccx` and destroyed the numbering for
    `test_3_1_square_plate_tables_2_to_5`, which became "3 1 square plate tables 2 to 5".
    """
    return name[len("test_") :] if name.startswith("test_") else name


def slugify(text: str, limit: int = 56) -> str:
    slug = re.sub(r"[^a-z0-9]+", "_", text.lower()).strip("_")
    return slug[:limit].strip("_") or "case"


def build_rows(
    group_id: str,
    slug: str,
    citation: str,
    scope: str,
    nodes: list[NodeInfo],
    reference_kind: str | None = None,
    skip: dict[str, str] | None = None,
    validation_helpers: set[str] | None = None,
    drop: set[str] | None = None,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """E1-E3 into rows: one row per collected node, none left unclaimed."""
    tree = ast.parse((REPO_ROOT / scope).read_text(encoding="utf-8"))
    consts = module_constants(tree)
    functions = {
        node.name: node for node in tree.body if isinstance(node, ast.FunctionDef)
    }
    # A method lives inside its class, not in the module body, so key it by the qualified name
    # pytest prints. Without this the node resolves to nothing and is reported as unclaimed.
    for parent in tree.body:
        if isinstance(parent, ast.ClassDef):
            for node in parent.body:
                if isinstance(node, ast.FunctionDef):
                    functions.setdefault(f"{parent.name}::{node.name}", node)
    report: dict[str, Any] = {
        "nodes": [],
        "unclaimed": [],
        # Why a node is not a row. The three reasons exist in the code already, and lumping
        # them together hid the difference between "this test compares nothing" and "the
        # tolerance could not be read": the first is a property of the test, the second is a
        # defect of the extractor, and the exclusion rule has to be written against the first.
        "unclaimed_reasons": {},
        "diverged": {},
        "prints": {},
    }
    rows: list[dict[str, Any]] = []
    reached: set[str] = set()
    matched_drops: set[str] = set()

    def unclaim(node: str, reason: str) -> None:
        report["unclaimed"].append(node)
        report["unclaimed_reasons"][node] = reason
    used: set[str] = set()
    for info in nodes:
        # A declared test is not a row, and it may not be claimed either. The declaration has to
        # work in both directions: it used to be consulted only for nodes the extractor could not
        # claim, so a declaration could not suppress a row whose assertions are machinery. The
        # tube torsion file is the case that forced it -- four torque-ruler and symmetry asserts
        # sit in rows of their own, beside a row holding the Bredt limit, and only the last is a
        # comparison against an independent reference.
        if skip and info.node in skip:
            report.setdefault("declared_not_rows", {})[info.node] = skip[info.node]
            continue
        func = functions.get(info.qualname)
        if func is None:
            unclaim(info.node, "not_in_ast")
            continue
        docstring = ast.get_docstring(func) or ""
        validates = (
            docstring.strip().splitlines()[0]
            if docstring.strip()
            else name_without_prefix(info.function)
        )
        # The test plus the declared helpers it reaches: a validation whose assertions live in a
        # local helper is otherwise invisible, and structure cannot say which helpers those are.
        extra = followed_helpers(func, functions, validation_helpers or set())
        reached.update(scope.name for scope in extra)
        sites, ignored = tolerance_sites(func, consts, info.params, extra)
        if ignored:
            report.setdefault("ignored_bounds", {})[info.node] = ignored
        if not sites:
            unclaim(info.node, "no_comparison")
            continue
        # Left out before anything else looks at them, so a dropped assertion cannot be counted,
        # paired with a print, or keep a row alive on its own.
        dropped = [site for site in sites if f"{scope}:{site.line}" in (drop or set())]
        if dropped:
            matched_drops.update(f"{scope}:{site.line}" for site in dropped)
            report.setdefault("dropped_asserts", {})[info.node] = [
                f"{scope}:{site.line}" for site in dropped
            ]
            sites = [site for site in sites if site not in dropped]
        if not sites:
            # Every assertion this test makes is machinery, so there is no comparison to record. A
            # row with no comparisons is not a smaller row; it is not a row.
            unclaim(info.node, "no_reference_comparison" if dropped else "no_comparison")
            continue
        resolved = [site for site in sites if site.value is not None]
        if not resolved:
            # Never write a comparison whose tolerance cannot be read: an invalid comparison is
            # worse than a reported gap, because it looks like evidence.
            report.setdefault("unresolved", {})[info.node] = [s.source for s in sites]
            unclaim(info.node, "tolerance_unresolved")
            continue
        # The unit of refusal is the comparison, not the test. An unreadable site is left out of
        # the row and reported, because it is often not a reference comparison at all: both
        # Scordelis-Lo smoothed tests were dropped whole by a symmetry assert carrying a computed
        # `atol`, while the comparison they exist for -- `rel < 0.05` against Lee & Lee Table 6 --
        # reads perfectly well.
        unreadable = [site.source for site in sites if site.value is None]
        if unreadable:
            report.setdefault("unresolved", {})[info.node] = unreadable
        # The owner is part of the id: two classes can hold same-named methods.
        row_id = (
            f"{slug}.{slugify(info.qualname, 40)}."
            f"{slugify(info.params or 'single', 40)}"
        )
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
            for scope in [func, *extra]
            for sub in ast.walk(scope)
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
        for site in resolved:
            # The group declares it. Guessing from the prose -- "paper" when the source text
            # happens to say paper, table or expected, "analytical" otherwise -- mislabelled the
            # hook test that compares against the paper's Table 14, because the assert reads
            # `rel_err < 0.03` and names none of those words, and it would have stamped
            # "analytical" on every symmetry and consistency check in the suite.
            kind = site.reference_kind or reference_kind
            reference: dict[str, Any] = {
                "kind": kind,
                "label": site.reference_name or site.reference_expr or site.source,
            }
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
    # A declared helper that no test in the scope reaches is a declaration that reads like coverage
    # and is not: the helper was renamed, or the test that used it went away.
    report["declared_helpers_used"] = sorted(reached)
    report["declared_helpers_stale"] = sorted((validation_helpers or set()) - reached)
    report["dropped_declared_stale"] = sorted((drop or set()) - matched_drops)
    return rows, report


def preserve_measurements(target: Path, rows: list[dict[str, Any]]) -> int:
    """Carry each row's measured evidence across a re-extraction.

    Re-extracting is how a row gains a better description or a new comparison, and it must not
    quietly discard the baselines: the digests *are* the regression record, and losing them
    would turn every comparison into a `new_baseline` on the next run, which reads like nothing
    was ever measured.

    Rows are matched by their `tests` set, which is stable as long as the test does not move, and
    comparisons are matched by their `tolerance.source` -- `<file>:<line>` -- and not by position.
    Position was the wrong key and dangerously so: dropping or reordering one comparison would have
    slid every later measurement one slot over, attaching a margin to a reference it was never
    measured against. A comparison whose source is not found keeps no measurement, and the count of
    those is reported rather than passed over, because a moved line means someone should decide
    whether the number still applies.
    """
    if not target.exists():
        return 0
    document = load_yaml(target)
    existing = document.get("rows") if isinstance(document, dict) else None
    if not isinstance(existing, list):
        return 0
    by_tests = {
        tuple(row.get("tests") or []): row for row in existing if isinstance(row, dict)
    }
    carried = 0
    unmatched: list[str] = []
    for row in rows:
        previous = by_tests.get(tuple(row.get("tests") or []))
        if previous is None:
            continue
        # Everything that is decided or measured rather than derived from the code: a human
        # flag, the history, the recorded expectation. A re-extraction improves the derived
        # fields, and must not drop these.
        if previous.get("history"):
            row["history"] = previous["history"]
        if previous.get("flags"):
            row["flags"] = previous["flags"]
        old_by_source: dict[str, dict[str, Any]] = {}
        for old in previous.get("comparisons") or []:
            if not isinstance(old, dict):
                continue
            source = str((old.get("tolerance") or {}).get("source") or "")
            if source:
                old_by_source[source] = old
        for comparison in row.get("comparisons") or []:
            if not isinstance(comparison, dict):
                continue
            source = str((comparison.get("tolerance") or {}).get("source") or "")
            old = old_by_source.get(source)
            if old is None:
                if comparison.get("tolerance"):
                    unmatched.append(source or "<no source>")
                continue
            for name in ("measured", "expected"):
                value = old.get(name)
                if value is not None:
                    comparison[name] = value
                    carried += 1
    if unmatched:
        print(
            f"no stored measurement matches {len(unmatched)} comparison(s) by source; a "
            "re-measure is due: "
            + ", ".join(sorted(unmatched)[:5]),
            file=sys.stderr,
        )
    return carried


def command_extract(args: argparse.Namespace) -> int:
    store = load_store(args.store)
    group = store.groups.get(args.group)
    if group is None:
        raise StoreError(f"unknown group: {args.group}")
    scope = args.scope or (group.get("source_files") or [None])[0]
    if not scope:
        raise StoreError(f"group {args.group} declares no source file; pass --scope")
    citation = args.citation or str(group.get("citation") or "")
    # A citation names a bibliography entry, and the store resolves one only for a paper
    # reference: `_validate_reference` refuses `citation` on any other kind. Requiring it of every
    # group would mean inventing a paper for CCX or OpenFAST -- the false reference this store
    # exists to prevent -- and it blocked every cross-code group, whose rows name their reference
    # through `reference.label` and carry the version in the group's `provenance_note`.
    if group.get("reference_kind") == "paper" and not citation:
        raise StoreError(
            f"group {args.group} is a paper group and declares no 'citation'; a paper reference "
            "has to resolve in references.yaml"
        )
    nodes = collect_nodes(str(scope))
    skip = {
        info.node: match[1]
        for info in nodes
        if (match := match_non_validation(group, info.node)) is not None
    }
    rows, report = build_rows(
        args.group,
        str(group.get("slug")),
        citation,
        str(scope),
        nodes,
        reference_kind=group.get("reference_kind"),
        skip=skip,
        validation_helpers=set(group.get("validation_helpers") or []),
        drop={str(item.get("at")) for item in (group.get("non_reference_asserts") or []) if isinstance(item, dict)},
    )
    suppressed = report.get("declared_not_rows") or {}
    # The declared nodes that were suppressed count as used patterns just like the unclaimed ones,
    # or a declaration that did its job would be reported as stale.
    declared_nodes, undeclared, stale = classify_unclaimed(
        group, list(report["unclaimed"]) + list(suppressed)
    )
    payload = {
        "scope": str(scope),
        "collected": len(nodes),
        "rows": len(rows),
        "claimed": len(report["nodes"]),
        "unclaimed": report["unclaimed"],
        "unclaimed_reasons": report["unclaimed_reasons"],
        "declared_non_validation": declared_nodes,
            "declared_not_rows": suppressed,
        "declared_helpers_used": report["declared_helpers_used"],
        "declared_helpers_stale": report["declared_helpers_stale"],
        "dropped_asserts": report.get("dropped_asserts") or {},
        "dropped_declared_stale": report["dropped_declared_stale"],
        "undeclared": undeclared,
        "stale_declarations": stale,
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
            # The row is written without it: say so, because a message that claims otherwise
            # sends the reader looking for a missing row that is not missing.
            print(f"UNREADABLE comparison (left out of the row): {node}")
            for source in sources:
                print(f"  - {source}")
        for node, reason in suppressed.items():
            print(f"declared not a row, no row written: {node}")
            print(f"  because: {reason}")
        if declared_nodes:
            print(f"not validation rows (declared): {len(declared_nodes)}")
            for node in declared_nodes:
                print(f"  - {node}")
        if undeclared:
            # Not an error any more. A test that makes no canonical comparison against a named
            # reference is a physics property or a machinery check, which is a legitimate thing for
            # it to be, and the contract says so by shape instead of asking a maintainer to declare
            # it. It stays in the report because a comparison written the old way -- a bare assert
            # against a literal -- lands here, and that is what a reviewer has to see.
            counts = Counter(
                report["unclaimed_reasons"].get(node, "unknown") for node in undeclared
            )
            breakdown = ", ".join(f"{r}: {n}" for r, n in counts.most_common())
            print(
                f"not a comparison against an independent reference: {len(undeclared)} "
                f"({breakdown})"
            )
            for node in undeclared:
                print(f"  - {node}")
        for pattern in stale:
            print(
                "DECLARED BUT MATCHING NOTHING: "
                f"{pattern} -- the test was renamed or deleted, or the prefix is wrong"
            )
        for source in report["dropped_declared_stale"]:
            print(
                f"DECLARED NON-REFERENCE ASSERTION NOT FOUND: {source} -- the line moved, or the "
                "assertion changed"
            )
        for name in report["declared_helpers_stale"]:
            print(
                f"DECLARED HELPER REACHED BY NO TEST: {name} -- renamed, or the tests that "
                "used it are gone"
            )
        for node, tokens in report["diverged"].items():
            print(f"param id diverges from the code's literals: {node} -> {tokens}")
        for node, calls in report["prints"].items():
            print(f"printed and maybe never asserted (candidate info_only): {node}")
            for call in calls:
                print(f"  - {call}")
    if args.write:
        target = args.store / "rows" / f"{args.group}-{group.get('slug')}.yaml"
        target.parent.mkdir(parents=True, exist_ok=True)
        carried = preserve_measurements(target, rows)
        dump_yaml(target, {"group": args.group, "rows": rows})
        print(
            f"wrote {display_path(target)} ({len(rows)} rows, {carried} measured "
            "comparison(s) carried over)"
        )
        probe = load_store(args.store)
        for finding in probe.errors():
            print(finding.as_text(), file=sys.stderr)
        if probe.errors():
            return EXIT_FINDINGS
    # A captured node that is neither claimed nor declared is the drift this store exists to
    # prevent, and a declaration that matches nothing is coverage that is not there. Neither is
    # reported by `check`: both need the collected set, and collecting means running pytest.
    if stale or report["declared_helpers_stale"] or report["dropped_declared_stale"]:
        return EXIT_FINDINGS
    return EXIT_OK


def command_coherence(args: argparse.Namespace) -> int:
    """Re-derive every group and report the row files that no longer match the code.

    A row file cites the line each comparison sits on. Editing a test moves those lines,
    and nothing else sees it: `check` cannot read the code, and comparing how many rows a
    group extracts compares a count rather than the bytes of the file, which is how a
    stale citation hid until `git status` showed it. This writes what the code says and
    compares in memory -- the same thing as `extract --write` followed by `git diff`,
    without needing git and with the groups named. A stale file is left refreshed, because
    the diff it prints is the fix.
    """
    store = load_store(args.store)
    groups = [args.group] if args.group else sorted(store.groups, key=int)
    stale: list[str] = []
    for group_id in groups:
        group = store.groups.get(group_id)
        if group is None:
            raise StoreError(f"unknown group: {group_id}")
        scope = (group.get("source_files") or [None])[0]
        target = args.store / "rows" / f"{group_id}-{group.get('slug')}.yaml"
        if not scope:
            # Group 10's rows are written by hand; there is nothing to derive them from.
            print(f"{target.name}: no source file to derive from, left alone")
            continue
        before = target.read_text(encoding="utf-8") if target.exists() else None
        quiet = argparse.Namespace(
            store=args.store, group=group_id, scope=scope, citation=None, json=False, write=True
        )
        with contextlib.redirect_stdout(io.StringIO()):
            code = command_extract(quiet)
        if code != EXIT_OK:
            print(f"{target.name}: the extraction itself did not finish cleanly")
            return code
        after = target.read_text(encoding="utf-8") if target.exists() else None
        if before != after:
            stale.append(target.name)
            print(f"{target.name}: the rows on disk were not what the code says, and are refreshed now")
    if stale:
        print(
            f"{len(stale)} of {len(groups)} row file(s) were stale. Read the change, and commit it "
            "if it is right."
        )
        return EXIT_FINDINGS
    print(f"coherent: {len(groups)} group(s) re-derive to the rows on disk")
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
        scope = data.get("scope")
        if isinstance(scope, str) and "/" in scope and not (REPO_ROOT / scope).exists():
            store.error(rel, f"scope points at a path that does not exist: {scope!r}")
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


def command_status(args: argparse.Namespace) -> int:
    """What the store covers, and what it does not yet.

    `check` cannot notice a validation file with no rows at all: enumerating collected nodes
    needs a pytest run and the built environment, which is why the row-to-node reconciliation is
    a separate step. The filesystem answers the coarser and more useful question on its own —
    which validation files have a group — so the state of the migration is a report instead of a
    guess.
    """
    store = load_store(args.store)
    grouped = {
        str(source)
        for group in store.groups.values()
        for source in group.get("source_files") or []
    }
    root = REPO_ROOT / "tests" / "validation"
    # Declared out of scope is not the same as forgotten, and only the second is a to-do.
    declared_out = {
        str(item) for entry in store.out_of_scope for item in entry.get("files") or []
    }
    ungrouped = sorted(
        str(path.relative_to(REPO_ROOT))
        for path in root.rglob("test_*.py")
        if str(path.relative_to(REPO_ROOT)) not in grouped
        and str(path.relative_to(REPO_ROOT)) not in declared_out
    )
    rows = []
    # By the number the group carries, not by where it sits in groups.yaml: the file's order is an
    # accident of editing, and someone looking for group 10 should find it between 9 and 11.
    for group_id, group in sorted(store.groups.items(), key=lambda item: int(item[0])):
        refs = [ref for ref in store.rows if ref.data.get("group") == group_id]
        views = [view for ref in refs for view in comparison_views(ref)]
        measured = [view for view in views if view["margin_pct"] is not None]
        rows.append(
            {
                "group": group_id,
                "files": len(group.get("source_files") or []),
                "rows": len(refs),
                "comparisons": len(views),
                "measured": len(measured),
                "near": sum(1 for view in views if "near" in view["flags"]),
                "gt5": sum(1 for view in views if "gt5" in view["flags"]),
            }
        )
    payload = {
        "groups": rows,
        "ungrouped_validation_files": ungrouped,
        "out_of_scope_files": sorted(declared_out),
        "gaps": len(store.gaps),
    }
    if args.json:
        print(json.dumps(payload, indent=2))
    else:
        print(f"{'group':8s} {'src':>5s} {'rows':>5s} {'cmp':>5s} {'measured':>8s} {'near':>4s} {'gt5':>4s}")
        for row in rows:
            # A group with no source file is not a broken group: its rows are driven by hand, and
            # no re-derivation can refresh them. Showing "hand" keeps a zero from reading as a gap.
            source = str(row["files"]) if row["files"] else "hand"
            print(
                f"{row['group']:8s} {source:>5s} {row['rows']:5d} {row['comparisons']:5d} "
                f"{row['measured']:8d} {row['near']:4d} {row['gt5']:4d}"
            )
        hand = [row["group"] for row in rows if not row["files"]]
        if hand:
            print(
                f"\ngroups with no source file: {', '.join(hand)}. Their rows are driven by hand, "
                "so src reads hand\n  and re-deriving cannot refresh them."
            )
        if not any(row["measured"] for row in rows if row["files"]):
            print(
                "\nthe only measured comparisons are the hand-driven ones of group 10. The rest of "
                "the store gets\n  its margins when regression --write records them, which is a "
                "run of the suite and not of this command."
            )
        print(f"\ngaps declared: {len(store.gaps)}")
        print(f"validation files declared out of scope: {len(declared_out)}")
        print(f"validation files neither grouped nor declared: {len(ungrouped)}")
        for path in ungrouped[:8]:
            print(f"  - {path}")
        if len(ungrouped) > 8:
            print(f"  ... and {len(ungrouped) - 8} more")
    return EXIT_OK


OUTCOME_RE = re.compile(
    r"^(?P<node>\S+::\S+)\s+(?P<status>PASSED|FAILED|ERROR|SKIPPED|XFAIL|XPASS)\b"
)
FAILURE_HEADER_RE = re.compile(r"^_{3,}\s*(?P<name>.+?)\s*_{3,}$")
# The suite's canonical assertions (tests/support/assertions.py) print where they landed and what
# bound they missed, so a failure says how far past the bound it went without anyone opening a file.
CANONICAL_MISS_RE = re.compile(
    r"(?P<measured>[\d.]+)% against .*?, above the (?P<bound>[\d.]+)% bound"
)


def triage_report(stdout: str) -> dict[str, Any]:
    """What a pytest run says, reduced to what a caller has to act on.

    Pure, so it is tested on captured output. The numbers it reports for a failure come from the
    canonical assertion's own message (tests/support/assertions.py), and that message format is part
    of the contract: a change to it should fail a test here rather than quietly stop reporting how
    far past a bound a failure went.
    """
    counts: Counter[str] = Counter()
    failed: list[dict[str, Any]] = []
    skipped: list[dict[str, str]] = []
    # pytest's own last summary line, which is authoritative. Per-test parsing is a convenience for
    # naming the failures, not a way to count them: `-v` wraps a long node id onto a second line and
    # a wrapped line matches no outcome, so counting with it under-reported a 114-test run as 62.
    summary_line = ""
    for raw in stdout.splitlines():
        stripped = raw.strip().strip("=").strip()
        if re.match(r"^\d+ (passed|failed|error|skipped|xfailed|xpassed)", stripped):
            summary_line = stripped
    for raw in stdout.splitlines():
        # pytest writes `<STATUS> [<n>] <file>:<line>: <reason>`: the path may carry its own colon
        # and a line range, so the path is matched lazily and the line suffix is optional, which
        # leaves the last colon as the separator before the reason.
        short = re.match(
            r"^(?P<status>SKIPPED|XFAIL|XPASS)\s+\[\d+\]\s+"
            r"(?P<where>[^\s:]+(?::\d+(?:-\d+)?)?):\s*(?P<reason>.*)$",
            raw.strip(),
        )
        if short is not None:
            skipped.append(
                {
                    "status": short.group("status"),
                    "where": short.group("where"),
                    "reason": short.group("reason").strip(),
                }
            )
        match = OUTCOME_RE.match(raw.strip())
        if match is not None:
            status = match.group("status")
            counts[status.lower()] += 1
            if status in {"FAILED", "ERROR"}:
                failed.append({"node": match.group("node"), "status": status})
    # The assertion message lives under the FAILURES header, and pytest wraps it across as many
    # `E ` lines as it needs -- a long canonical message is exactly the case that wraps, so taking
    # only the first line loses the "above the bound" half and the report goes quiet about the
    # distance. All of them are joined.
    messages: dict[str, list[str]] = {}
    current: str | None = None
    for raw in stdout.splitlines():
        header = FAILURE_HEADER_RE.match(raw.strip())
        if header is not None:
            current = header.group("name").strip()
            continue
        if current and raw.strip().startswith("E "):
            messages.setdefault(current, []).append(raw.strip()[2:].strip())
    for item in failed:
        name = item["node"].split("::")[-1]
        parts = messages.get(name) or messages.get(item["node"]) or []
        message = " ".join(parts).strip()
        item["message"] = message
        miss = CANONICAL_MISS_RE.search(message)
        if miss is not None:
            measured = float(miss.group("measured"))
            bound = float(miss.group("bound"))
            item["measured_pct"] = measured
            item["bound_pct"] = bound
            item["over_pp"] = round(measured - bound, 4)
            item["times_bound"] = round(measured / bound, 3) if bound else None
    return {
        "summary": summary_line,
        "counts": dict(counts),
        "failed": failed,
        "skipped": skipped,
    }


def command_triage(args: argparse.Namespace) -> int:
    """Run a scope and report what failed and by how much, in one place.

    `regression` also runs a scope, but it asks whether the stored numbers moved rather than whether
    the tests passed: it ignores pytest's exit status on purpose, because a failed test still prints
    its residual, which is all it needs. This verb is the other half, for the caller who has to act
    on a red suite. It names the test, what it compared and against what reference, how far past the
    bound it went, and which store row the comparison belongs to -- so nobody has to open the file to
    find out, which is the whole cost of a red run.
    """
    scope = args.scope or (
        " ".join(
            str(item)
            for item in (load_store(args.store).groups.get(str(args.group), {}) or {}).get(
                "source_files"
            )
            or []
        )
        if args.group
        else "tests"
    )
    if not scope:
        raise StoreError(f"group {args.group} declares no source_files; pass --scope")
    # `-ra` prints the short summary for everything that did not pass, which is where a skipped
    # module says why it skipped -- and in this suite a missing external tool makes a row skip rather
    # than fail, so the reason is the answer, not a footnote.
    # `--first` stops at the first failure. A scope can be minutes long -- the blade mesh study
    # solves CCX on every mesh in it -- and the first failure is usually the one to act on, so the
    # caller who wants a verdict rather than a full census should not have to wait for the rest.
    flags = ["-o", "addopts=", "-s", "-v", "-ra", *(["-x"] if args.first else [])]
    completed = subprocess.run(
        [sys.executable, "-m", "pytest", *flags, *scope.split()],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
        check=False,
    )
    report = triage_report(completed.stdout)
    failed = report["failed"]
    skipped = report["skipped"]
    for ref in load_store(args.store).rows:
        for item in failed:
            if item["node"] in (ref.data.get("tests") or []):
                item["row"] = ref.data.get("id")
    if args.json:
        print(json.dumps({"scope": scope, **report}, indent=2))
    else:
        # pytest's line first: it counts what actually ran, and the per-test parse is only used to
        # name failures. A scope that printed nothing still says its exit status rather than nothing.
        summary = report["summary"] or f"no report (pytest exit {completed.returncode})"
        print(f"{scope}: {summary}")
        for item in skipped:
            print(f"  {item['status']} {item['where']} -- {item['reason'] or 'no reason given'}")
        if not summary and not skipped and not failed:
            for line in completed.stdout.splitlines()[-4:]:
                print(f"  {line}")
        for item in failed:
            print(f"\n{item['status']}  {item['node']}")
            if item.get("row"):
                print(f"  store row : {item['row']}")
            if item.get("measured_pct") is not None:
                print(
                    f"  measured  : {item['measured_pct']}% against a "
                    f"{item['bound_pct']}% bound -> over by {item['over_pp']}pp "
                    f"({item['times_bound']}x the bound)"
                )
            if item["message"]:
                print(f"  assertion : {item['message'][:160]}")
    return EXIT_FINDINGS if failed else EXIT_OK


def command_kinds(args: argparse.Namespace) -> int:
    """Every comparison whose reference kind is undeclared, with the code that decides it.

    Policy rule 7 reads the kind from the comparison, never from the file, so this verb prints the
    comparison and the line of code behind it. The operator decides; the decision is written with
    `set`. It exists because a mixed file needs one decision per comparison -- and every parity
    module examined so far is mixed -- which makes reading them out of the tests by hand the whole
    cost of extending the store. This is the review sheet for that cost.
    """
    # No error guard here on purpose. This verb is read-only, and the state it exists to review is
    # an error state: an undeclared reference kind IS a `check` error. Refusing to run whenever the
    # store has errors would make the review sheet unavailable exactly when it is needed -- the same
    # trap that made `set` unable to repair a null kind.
    store = load_store(args.store)
    pending: list[dict[str, Any]] = []
    for ref in store.rows:
        if args.group and str(ref.data.get("group")) != args.group:
            continue
        for comparison in ref.data.get("comparisons") or []:
            if not isinstance(comparison, dict):
                continue
            reference = comparison.get("reference")
            if not isinstance(reference, dict) or reference.get("kind") is not None:
                continue
            tolerance = comparison.get("tolerance") or {}
            source = str(tolerance.get("source") or "")
            path, _, line = source.rpartition(":")
            code = ""
            if path and line.isdigit():
                target = REPO_ROOT / path
                if target.exists():
                    body = target.read_text(encoding="utf-8").splitlines()
                    number = int(line)
                    if 0 < number <= len(body):
                        code = body[number - 1].strip()
            pending.append(
                {
                    "row": ref.data.get("id"),
                    "comparison": comparison.get("label"),
                    "reference_label": reference.get("label"),
                    "source": source,
                    "code": code,
                }
            )
    if args.json:
        print(json.dumps(pending, indent=2))
    elif not pending:
        print("every comparison declares its reference kind")
    for item in pending:
        print(f"{item['row']}")
        print(f"  comparison : {item['comparison']}")
        print(f"  reference  : {item['reference_label']}")
        print(f"  code       : {str(item['code'])[:92]}")
        print(f"  at         : {item['source']}")
    return EXIT_OK


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


def evidence_digest(printed: dict[str, Any]) -> str:
    """A digest of one comparison's printed evidence.

    The digest is the *detector*: any change in the printed digits changes it, so the check is
    exact and needs no tolerance. A tolerance in the detector would swallow exactly the small
    movements a regression check exists to catch; the tolerance belongs in the test's
    assertion, and the numeric delta beside the digest is what *measures* the change once the
    digest has said which comparison moved.

    It hashes the **printed strings**, not the parsed numbers, and it is computed once and
    stored by the same code that compares it. Hashing parsed values from one side and formatted
    strings from the other made every comparison report as changed with a delta of zero — the
    detector catching its own inconsistency, which is the whole point of having one.
    """
    payload = json.dumps([printed.get(key) for key in ("value", "error", "expected")])
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()[:16]



def sources_digest(group: dict[str, Any]) -> str:
    """A digest of a group's source files: the test's own contract.

    Kept beside the rows it separates the two kinds of change a regression report has to tell
    apart: the numbers moved while the test stood still, or the test itself changed and the
    numbers are expected to move.
    """
    digest = hashlib.sha256()
    for source in sorted(str(item) for item in group.get("source_files") or []):
        digest.update(source.encode("utf-8"))
        path = REPO_ROOT / source
        if path.exists():
            digest.update(path.read_bytes())
    return digest.hexdigest()[:16]


def load_source_digests(store: Store) -> None:
    """The stored source digests, keyed by group id."""
    path = store.root / SOURCES_FILE
    if not path.exists():
        return
    rel = display_path(path)
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        store.error(rel, f"could not be read: {exc}")
        return
    digests = data.get("digests") if isinstance(data, dict) else None
    if not isinstance(digests, dict):
        store.error(rel, "expected a mapping with a 'digests' object")
        return
    for group_id, digest in digests.items():
        if group_id not in store.groups:
            store.error(f"{rel}[{group_id}]", f"unknown group: {group_id!r}")
        if not isinstance(digest, str) or not digest:
            store.error(f"{rel}[{group_id}]", f"invalid digest: {digest!r}")
        store.source_digests[str(group_id)] = str(digest)


def write_source_digest(store: Store, group_id: str, digest: str) -> None:
    path = store.root / SOURCES_FILE
    data: dict[str, Any] = {"version": 1, "digests": {}}
    if path.exists():
        try:
            loaded = json.loads(path.read_text(encoding="utf-8"))
            if isinstance(loaded, dict) and isinstance(loaded.get("digests"), dict):
                data = loaded
        except json.JSONDecodeError:
            pass
    data["digests"][group_id] = digest
    path.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def scalar(text: Any) -> Any:
    """A printed value as a number when it is one, and as text when it is not.

    Every residual pattern captures with `\\S+`, so what it produces is always a string. Storing
    a numeric expectation as text made `expected` asymmetric with `measured.raw`, which is a
    number, and made any numeric query on it compare strings. A reference side that is not one
    number -- a mode table, a cell list -- cannot become one, so it stays text.
    """
    if text is None:
        return None
    try:
        return float(text)
    except (TypeError, ValueError):
        return text


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
        measured = comparison.get("measured") or {}
        baseline = _float(measured.get("margin_pct"))
        digest = evidence_digest(printed)
        stored = measured.get("digest")
        if baseline is None or not stored:
            verdict = "new_baseline"
        elif stored == digest:
            verdict = "same"
        else:
            verdict = "changed"
        delta = None
        if current is not None and baseline is not None:
            delta = round(current - baseline, 6)
        results.append(
            {
                "row": ref.id,
                "comparison": index,
                "label": comparison.get("label"),
                "verdict": verdict,
                "baseline": baseline,
                "current": current,
                "delta": delta,
                "digest": digest,
                "stored_digest": stored,
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
        "digest": result["digest"],  # computed once, by the code that compares it
    }
    if comparison.get("expected") is None and result.get("expected") is not None:
        comparison["expected"] = scalar(result["expected"])
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

    sources_now = sources_digest(group)
    sources_stored = store.source_digests.get(str(args.group))
    sources_changed = bool(sources_stored) and sources_stored != sources_now

    if args.write:
        for result in results:
            index = result.get("comparison")
            if result["verdict"] not in {"same", "changed", "new_baseline"}:
                continue
            if not isinstance(index, int):
                continue
            write_measurement(store, store.by_id()[result["row"]], index, result)
        write_source_digest(store, str(args.group), sources_now)
        print(f"wrote {len(results)} measurement(s) and the source digest")

    failing = [
        item for item in results if item["verdict"] in {"changed", "unmapped", "unclaimed"}
    ]
    if sources_changed:
        failing.append(
            {
                "verdict": "test_changed",
                "row": None,
                "detail": ", ".join(str(item) for item in group.get("source_files") or []),
            }
        )
    if args.json:
        print(json.dumps(results, indent=2))
    else:
        for item in results:
            if item["verdict"] == "same":
                continue
            where = item.get("row") or item.get("detail")
            if item["verdict"] == "changed":
                delta = item.get("delta")
                moved = f"{delta:+g}" if isinstance(delta, (int, float)) else "n/a"
                print(
                    f"CHANGED  {where} [{item['label']}]: {item['baseline']}% -> "
                    f"{item['current']}% (delta {moved})"
                )
            elif item["verdict"] == "test_changed":
                print(
                    f"TEST MOVED  {where}: the source digest changed, so its numbers are "
                    "expected to move too"
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
        print(
            "source digest: "
            + ("moved, the test itself changed" if sources_changed else "unchanged")
            if sources_stored
            else "source digest: none stored yet; --write records it"
        )
    return EXIT_FINDINGS if failing else EXIT_OK


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

    coherence = subparsers.add_parser(
        "coherence",
        help="re-derive every group and report the row files that no longer match the code",
    )
    coherence.add_argument("--group", help="one group instead of all of them")
    coherence.set_defaults(func=command_coherence)

    gaps = subparsers.add_parser("gaps", help="what the suite does not validate")
    gaps.add_argument("--json", action="store_true")
    gaps.set_defaults(func=command_gaps)

    triage = subparsers.add_parser(
        "triage",
        help="run a scope and say which tests failed and by how much",
    )
    triage.add_argument("--group", help="run this group's source_files")
    triage.add_argument("--scope", help="what to run (default: tests)")
    triage.add_argument(
        "--first", action="store_true", help="stop at the first failure instead of running the rest"
    )
    triage.add_argument("--json", action="store_true")
    triage.set_defaults(func=command_triage)

    kinds = subparsers.add_parser(
        "kinds",
        help="undeclared reference kinds, with the code that decides each one",
    )
    kinds.add_argument("--group", help="only rows of this group")
    kinds.add_argument("--json", action="store_true")
    kinds.set_defaults(func=command_kinds)

    status = subparsers.add_parser(
        "status", help="what the store covers, and which validation files it does not"
    )
    status.add_argument("--json", action="store_true")
    status.set_defaults(func=command_status)

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
