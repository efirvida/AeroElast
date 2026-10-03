"""Tests for tools/validation_matrix.py.

These tests are deliberately not part of the scientific suite: `pyproject.toml`
sets `testpaths = ["tests"]`, so a bare `pytest` never collects them. Run them
explicitly:

    python -m pytest tools/tests

Keeping them out matters: the validation matrix reconciles its row count against
`pytest --collect-only`, and a tool test inside `tests/` would have to appear in
the matrix as a validation row.
"""

from __future__ import annotations

import copy
import importlib.util
import ast
import json
import subprocess
import sys
from pathlib import Path
from typing import Any

import pytest
import yaml

REPO_ROOT = Path(__file__).resolve().parents[2]
TOOL = REPO_ROOT / "tools" / "validation_matrix.py"
REAL_STORE = REPO_ROOT / "docs" / "validation"


def _load_tool_module() -> Any:
    """Import the tool from its path: tools/ is not a package."""
    spec = importlib.util.spec_from_file_location("validation_matrix_tool", TOOL)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    # @dataclass looks its owner up in sys.modules by cls.__module__.
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


GROUP_SOURCE = "tests/validation/benchmarks/test_ko2017_performance.py"
NODE_ID = (
    "tests/validation/benchmarks/test_ko2017_performance.py::"
    "test_3_1_square_plate_tables_2_to_5"
    "[expected_table2_30-expected_table4_50-False-0.01-100.0]"
)

GROUPS = {
    "version": 1,
    "groups": [
        {
            "id": "3",
            "title": "`tests/validation/benchmarks/test_ko2017_performance.py`",
            "slug": "ko2017",
            "source_files": [GROUP_SOURCE],
        }
    ],
}

FLAGS = {
    "version": 1,
    "flags": [
        {"id": "near", "derived": "slack_pp < 1.0", "user_facing": True},
        {"id": "gt5", "derived": "value > 0.05", "user_facing": True},
        {"id": "info_only", "section": "9.5", "label": "Printed and never asserted"},
        {"id": "tautology", "section": "9.1", "label": "Tautological reference"},
    ],
}

REFERENCES = {
    "version": 1,
    "header": "This file is the canonical bibliography for AeroElast.",
    "sections": [
        {"title": "About the recovered PDFs", "intro": "Verification convention used below."},
        {"title": "1. Shell element formulations"},
        {"title": "5. Verification benchmarks, reference models and validation literature"},
    ],
    "references": [
        {
            "key": "ko2017_perf",
            "section": "1. Shell element formulations",
            "kind": "journal",
            "authors": ["Ko, Y.", "Lee, Y.", "Lee, P.-S.", "Bathe, K.-J."],
            "title": "Performance of the MITC3+ and MITC4+ shell elements",
            "venue": "Computers and Structures",
            "volume": "193",
            "pages": "187-206",
            "year": 2017,
            "doi": "10.1016/j.compstruc.2017.08.003",
            "doi_status": "verified",
            "held": True,
            "held_files": [GROUP_SOURCE],
            "verification": "verified_pdf",
            "verification_note": "read off the held copy",
            "cited_by_declared": [],
            "code_mentions": ["Ko2017"],
        },
        {
            "key": "escalera2023",
            "section": "5. Verification benchmarks",
            "kind": "conference",
            "authors": ["Escalera Mendoza, A."],
            "title": "An open-source NuMAD model for the IEA 15 MW blade",
            "year": 2023,
            "doi": None,
            "doi_status": "to_verify",
            "held": False,
            "cited_by_declared": [],
        },
    ],
}


def make_row(**overrides: Any) -> dict[str, Any]:
    """A valid row: one own result compared against two references."""
    row: dict[str, Any] = {
        "id": "ko2017.square_plate.reg_clamped",
        "group": "3",
        "title": "test_3_1_square_plate_tables_2_to_5[reg clamped]",
        "tests": [NODE_ID],
        "validates": "clamped square plate centre deflection, N=16",
        "result": {"text": "0.9996 normalized", "raw": 0.9996, "unit": "w / w_Kirchhoff"},
        "flags": [],
        "comparisons": [
            {
                "label": "paper cell, N=16",
                "asserted": True,
                "reference": {
                    "kind": "paper",
                    "label": "Table 2, MITC4/MITC4+ N=16: 0.9984",
                    "citation": "ko2017_perf",
                },
                "tolerance": {
                    "kind": "rtol",
                    "value": 0.05,
                    "source": "assert_relative_error",
                    "justified": True,
                },
                "expected": 0.9984,
                "measured": {
                    "status": "measured",
                    "run": "34e2328",
                    "date": "2026-10-02",
                    "raw": 0.9996,
                    "margin_pct": 0.12,
                    "text": "0.9996 (0.12%)",
                },
            },
            {
                "label": "Kirchhoff closed form used for normalization",
                "asserted": True,
                "reference": {"kind": "analytical", "label": "alpha p L^4 / D, alpha = 1.267e-3"},
                "tolerance": {"kind": "rtol", "value": 0.05, "source": "same assert", "justified": True},
                "measured": {"status": "measured", "margin_pct": 0.0, "text": "1.0000 (0.00%)"},
            },
        ],
    }
    for key, value in overrides.items():
        row[key] = copy.deepcopy(value)
    return row


def write_store(
    tmp_path: Path,
    rows: list[dict[str, Any]],
    groups: dict[str, Any] | None = None,
    flags: dict[str, Any] | None = None,
    references: dict[str, Any] | None = None,
) -> Path:
    store = tmp_path / "validation"
    (store / "rows").mkdir(parents=True, exist_ok=True)
    (store / "groups.yaml").write_text(yaml.safe_dump(groups or GROUPS), encoding="utf-8")
    (store / "flags.yaml").write_text(yaml.safe_dump(flags or FLAGS), encoding="utf-8")
    (store / "references.yaml").write_text(
        yaml.safe_dump(references or REFERENCES), encoding="utf-8"
    )
    if rows:
        (store / "rows" / "3-ko2017.yaml").write_text(
            yaml.safe_dump({"group": "3", "rows": rows}), encoding="utf-8"
        )
    return store


def run(store: Path, *args: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, str(TOOL), "--store", str(store), *args],
        capture_output=True,
        text=True,
        check=False,
    )


def run_check(store: Path, *args: str) -> subprocess.CompletedProcess[str]:
    return run(store, "check", *args)


def test_three_valid_rows_pass(tmp_path: Path) -> None:
    rows = [
        make_row(),
        make_row(id="ko2017.square_plate.dist_clamped", flags=["tautology"]),
        make_row(id="ko2017.circular_plate.clamped"),
    ]
    store = write_store(tmp_path, rows)
    completed = run_check(store)
    assert completed.returncode == 0, completed.stdout + completed.stderr
    assert "3 row(s), 6 comparison(s)" in completed.stdout


def test_real_store_is_clean() -> None:
    """The shipped store must satisfy the same rules as the fixture."""
    completed = run_check(REAL_STORE)
    assert completed.returncode == 0, completed.stdout + completed.stderr


def test_multi_reference_row_keeps_two_tolerances(tmp_path: Path) -> None:
    """One result against CCX and against an analytical bar, as section 4.3 does."""
    row = make_row(
        comparisons=[
            {
                "label": "the CCX leg",
                "asserted": True,
                "reference": {"kind": "code", "label": "CCX 2.23, S8R"},
                "tolerance": {
                    "kind": "rtol",
                    "value": 0.015,
                    "source": "CCX_MEMBRANE_TOL",
                    "justified": True,
                },
                "measured": {"status": "measured", "margin_pct": 0.65, "text": "41.94 um (0.65%)"},
            },
            {
                "label": "the analytical bar",
                "asserted": True,
                "reference": {"kind": "analytical", "label": "delta = a22 (P/B) L"},
                "tolerance": {
                    "kind": "rtol",
                    "value": 0.02,
                    "source": "CLT_ANALYTICAL_TOL",
                    "justified": True,
                },
                "measured": {"status": "measured", "margin_pct": 1.25, "text": "42.19 um (1.25%)"},
            },
        ]
    )
    store = write_store(tmp_path, [row])
    completed = run_check(store, "--json")
    assert completed.returncode == 0, completed.stdout + completed.stderr
    payload = yaml.safe_load(completed.stdout)
    assert payload["rows"] == 1
    assert payload["comparisons"] == 2


@pytest.mark.parametrize(
    ("mutation", "expected"),
    [
        ("derived_flag_stored", "derived flag 'near' must not be stored"),
        ("unknown_flag", "unknown flag id: 'bogus'"),
        ("measured_without_margin", "requires a numeric margin_pct"),
        ("margin_without_measured", "must be null when status"),
        ("tolerance_above_rule", "without justification"),
        ("unasserted_without_info_only", "requires the 'info_only' flag"),
        ("no_comparisons", "'comparisons' must be a non-empty list"),
        ("empty_tests", "requires evidence: out_of_band"),
        ("unknown_group", "unknown group"),
        ("bad_id_slug", "does not match the group slug"),
        ("missing_label", "label must be a non-empty string"),
        ("unknown_row_key", "unknown keys"),
        ("bad_date", "must be YYYY-MM-DD"),
        ("expected_as_string", "expected is the numeric string"),
        ("no_reference_kind", "reference.kind is not declared"),
        ("citation_on_non_paper", "only meaningful when kind == 'paper'"),
    ],
)
def test_check_reports_one_error(tmp_path: Path, mutation: str, expected: str) -> None:
    row = make_row()
    if mutation == "derived_flag_stored":
        row["flags"] = ["near"]
    elif mutation == "unknown_flag":
        row["flags"] = ["bogus"]
    elif mutation == "measured_without_margin":
        del row["comparisons"][0]["measured"]["margin_pct"]
    elif mutation == "margin_without_measured":
        row["comparisons"][0]["measured"]["status"] = "not_printed"
    elif mutation == "tolerance_above_rule":
        row["comparisons"][0]["tolerance"] = {
            "kind": "rtol",
            "value": 0.15,
            "source": "15% bound",
            "justified": False,
        }
    elif mutation == "unasserted_without_info_only":
        row["comparisons"][1]["asserted"] = False
        row["comparisons"][1]["measured"] = {"status": "n/a"}
    elif mutation == "no_comparisons":
        row["comparisons"] = []
    elif mutation == "empty_tests":
        row["tests"] = []
    elif mutation == "unknown_group":
        row["group"] = "99"
    elif mutation == "bad_id_slug":
        row["id"] = "nope.row.here"
    elif mutation == "missing_label":
        del row["comparisons"][0]["label"]
    elif mutation == "unknown_row_key":
        row["surprise"] = True
    elif mutation == "bad_date":
        row["comparisons"][0]["measured"]["date"] = "02/10/2026"
    elif mutation == "expected_as_string":
        row["comparisons"][0]["expected"] = "0.9984"
    elif mutation == "no_reference_kind":
        row["comparisons"][0]["reference"]["kind"] = None
    elif mutation == "citation_on_non_paper":
        row["comparisons"][1]["reference"]["citation"] = "ko2017_perf"

    store = write_store(tmp_path, [row])
    completed = run_check(store)
    assert completed.returncode == 1, completed.stdout + completed.stderr
    assert expected in completed.stdout


def test_unasserted_comparison_with_info_only_passes(tmp_path: Path) -> None:
    """The section 9.5 shape: printed, never asserted, and flagged as such."""
    row = make_row(flags=["info_only"])
    row["comparisons"][0]["asserted"] = False
    row["comparisons"][0]["measured"] = {
        "status": "n/a",
        "margin_pct": None,
        "text": "printed errors 523% to 6.2e8%",
    }
    store = write_store(tmp_path, [row])
    completed = run_check(store)
    assert completed.returncode == 0, completed.stdout + completed.stderr


def test_row_group_must_match_file_group(tmp_path: Path) -> None:
    second = make_row(id="ko2017.square_plate.other")
    second["group"] = "3"
    store = write_store(tmp_path, [make_row()])
    (store / "rows" / "3-ko2017.yaml").write_text(
        yaml.safe_dump({"group": "4", "rows": [second]}), encoding="utf-8"
    )
    completed = run_check(store)
    assert completed.returncode == 1
    assert "does not match the file group" in completed.stdout


def test_duplicate_row_id_is_reported(tmp_path: Path) -> None:
    store = write_store(tmp_path, [make_row(), make_row()])
    completed = run_check(store)
    assert completed.returncode == 1
    assert "duplicate row id" in completed.stdout


def test_group_scope_limits_rows(tmp_path: Path) -> None:
    store = write_store(tmp_path, [make_row()])
    completed = run_check(store, "--group", "3")
    assert completed.returncode == 0
    assert "1 row(s)" in completed.stdout
    completed = run_check(store, "--group", "4")
    assert completed.returncode == 2
    assert "unknown group: 4" in completed.stderr


def test_missing_store_is_a_usage_error(tmp_path: Path) -> None:
    completed = run_check(tmp_path / "nowhere")
    assert completed.returncode == 2
    assert "does not exist" in completed.stderr


# --------------------------------------------------------------------------- #
# T2: read verbs and `set`
# --------------------------------------------------------------------------- #


def _mixed_store(tmp_path: Path) -> Path:
    """One row whose three comparisons are ordinary, near, and above 5%."""
    row = make_row()
    row["comparisons"] = [
        {
            "label": "ordinary",
            "asserted": True,
            "reference": {"kind": "code", "label": "CCX 2.23, S4"},
            "tolerance": {"kind": "rtol", "value": 0.05, "source": "tol=0.05", "justified": True},
            "measured": {"status": "measured", "margin_pct": 0.12, "text": "0.12%"},
        },
        {
            "label": "tight",
            "asserted": True,
            "reference": {"kind": "analytical", "label": "F L^3/(3 E I)"},
            "tolerance": {"kind": "rtol", "value": 0.01, "source": "tol=0.01", "justified": True},
            "measured": {"status": "measured", "margin_pct": 0.5, "text": "0.50%"},
        },
        {
            "label": "wide",
            "asserted": True,
            "reference": {"kind": "paper", "label": "Escalera modes", "citation": "escalera2023"},
            "tolerance": {
                "kind": "rtol",
                "value": 0.15,
                "source": "15% bound",
                "justified": True,
                "justification": "the article's own scatter is 11.8%",
            },
            "measured": {"status": "measured", "margin_pct": 3.0, "text": "3.00%"},
        },
    ]
    return write_store(tmp_path, [row])


def test_list_row_unit_shows_derived_worst_values(tmp_path: Path) -> None:
    store = _mixed_store(tmp_path)
    completed = run(store, "list", "--unit", "row")
    assert completed.returncode == 0, completed.stdout + completed.stderr
    line = completed.stdout.splitlines()[0]
    assert "ko2017.square_plate.reg_clamped" in line
    assert "3" in line  # three comparisons
    assert "3%" in line  # worst margin is the wide comparison
    assert "+0.5%" in line  # tightest slack: 1% tolerance - 0.5% margin
    assert "gt5" in line
    assert "near" in line


def test_list_comparison_unit_flags_near_and_gt5(tmp_path: Path) -> None:
    store = _mixed_store(tmp_path)
    near = run(store, "list", "--unit", "comparison", "--near")
    assert near.returncode == 0, near.stdout + near.stderr
    assert "tight" in near.stdout
    assert "ordinary" not in near.stdout
    assert "wide" not in near.stdout

    gt5 = run(store, "list", "--unit", "comparison", "--gt5")
    assert gt5.returncode == 0, gt5.stdout + gt5.stderr
    assert "wide" in gt5.stdout
    assert "tight" not in gt5.stdout


def test_list_sorts_by_slack_and_selects_fields(tmp_path: Path) -> None:
    store = _mixed_store(tmp_path)
    completed = run(
        store, "list", "--unit", "comparison", "--sort", "slack", "--fields", "label,slack_pp"
    )
    assert completed.returncode == 0, completed.stdout + completed.stderr
    lines = [line for line in completed.stdout.splitlines() if "|" in line]
    assert "tight" in lines[0]
    assert "ordinary" in lines[1]
    assert "wide" in lines[2]

    unknown = run(store, "list", "--fields", "nope")
    assert unknown.returncode == 2
    assert "unknown field(s): nope" in unknown.stderr


def test_find_matches_title_and_narrows_by_group(tmp_path: Path) -> None:
    store = _mixed_store(tmp_path)
    completed = run(store, "find", "clamped")
    assert completed.returncode == 0
    assert "1 match(es)" in completed.stdout
    assert "ko2017.square_plate.reg_clamped" in completed.stdout

    assert "0 match(es)" in run(store, "find", "clamped", "--group", "9").stdout
    assert "1 match(es)" in run(store, "find", "clamped", "--group", "3").stdout
    empty = run(store, "find", "zzz")
    assert "0 match(es)" in empty.stdout


def test_get_round_trips_the_stored_record(tmp_path: Path) -> None:
    store = _mixed_store(tmp_path)
    as_json = run(store, "get", "ko2017.square_plate.reg_clamped", "--json")
    assert as_json.returncode == 0
    stored = yaml.safe_load(as_json.stdout)
    assert stored["comparisons"][2]["label"] == "wide"

    as_yaml = run(store, "get", "ko2017.square_plate.reg_clamped")
    assert yaml.safe_load(as_yaml.stdout) == stored

    derived = run(store, "get", "ko2017.square_plate.reg_clamped", "--comparisons")
    assert "tight" in derived.stdout
    assert "gt5" in derived.stdout

    missing = run(store, "get", "ko2017.nope.here")
    assert missing.returncode == 2
    assert "unknown row id" in missing.stderr


def test_headline_reports_worst_margin_and_above_rule_count(tmp_path: Path) -> None:
    store = _mixed_store(tmp_path)
    groups = {
        "version": 1,
        "groups": [
            {
                "id": "3",
                "title": "`tests/validation/benchmarks/test_ko2017_performance.py`",
                "slug": "ko2017",
                "source_files": [GROUP_SOURCE],
                "headline": "31 benchmark cases; worst 2.73%",
            }
        ],
    }
    store = write_store(tmp_path, [make_row()], groups=groups)
    completed = run(store, "headline")
    assert completed.returncode == 0
    assert "| § | headline | rows | cmp | worst margin | >5% |" in completed.stdout
    assert "31 benchmark cases; worst 2.73%" in completed.stdout
    assert "| 3 |" in completed.stdout

    as_json = yaml.safe_load(run(store, "headline", "--json").stdout)
    assert as_json[0]["headline"] == "31 benchmark cases; worst 2.73%"
    assert as_json[0]["comparisons"] == 2


def test_set_edits_one_field_and_leaves_the_store_valid(tmp_path: Path) -> None:
    store = _mixed_store(tmp_path)
    completed = run(
        store,
        "set",
        "ko2017.square_plate.reg_clamped",
        "comparisons[2].measured.margin_pct=2.5",
        "comparisons[2].measured.text=2.50%",
        "measured.note=ignored",
    )
    assert completed.returncode == 2
    assert "unknown field 'measured'" in completed.stderr
    # The refused write must not have touched the file.
    assert yaml.safe_load(run(store, "get", "ko2017.square_plate.reg_clamped", "--json").stdout)[
        "comparisons"
    ][2]["measured"]["margin_pct"] == 3.0

    completed = run(
        store,
        "set",
        "ko2017.square_plate.reg_clamped",
        "comparisons[2].measured.margin_pct=2.5",
        "comparisons[2].measured.text=2.50%",
    )
    assert completed.returncode == 0, completed.stdout + completed.stderr
    after = run(store, "get", "ko2017.square_plate.reg_clamped", "--json")
    assert yaml.safe_load(after.stdout)["comparisons"][2]["measured"]["margin_pct"] == 2.5
    assert run_check(store).returncode == 0


def test_set_nested_mapping_and_json_list(tmp_path: Path) -> None:
    store = _mixed_store(tmp_path)
    completed = run(
        store,
        "set",
        "ko2017.square_plate.reg_clamped",
        'comparisons[0].reference.label="CCX 2.23, S8"',
        "comparisons[0].measured.status=measured",
        "comparisons[0].measured.margin_pct=1",
        f'tests=["{NODE_ID}"]',
        "comparisons[1].tolerance.justification=justified in the module comment",
    )
    assert completed.returncode == 0, completed.stdout + completed.stderr
    after = yaml.safe_load(run(store, "get", "ko2017.square_plate.reg_clamped", "--json").stdout)
    assert after["comparisons"][0]["reference"]["label"] == "CCX 2.23, S8"
    assert after["tests"] == [NODE_ID]


@pytest.mark.parametrize(
    ("assignment", "expected"),
    [
        ("near=true", "is derived and must not be written"),
        ("gt5=true", "is derived and must not be written"),
        ("id=ko2017.other.case", "'id' is immutable"),
        ("comparisons.bogus=1", "is a list: address an element"),
        ("comparisons[9].label=x", "index out of range"),
        ("comparisons[0].reference.bogus=x", "unknown field 'bogus'"),
        ("title", "expected PATH=VALUE"),
        ("comparisons[0].measured.status=nonsense", "must be one of"),
    ],
)
def test_set_refuses_bad_paths_and_values(tmp_path: Path, assignment: str, expected: str) -> None:
    store = _mixed_store(tmp_path)
    before = run(store, "get", "ko2017.square_plate.reg_clamped", "--json").stdout
    completed = run(store, "set", "ko2017.square_plate.reg_clamped", assignment)
    assert completed.returncode == 2, completed.stdout + completed.stderr
    assert expected in completed.stderr
    assert run(store, "get", "ko2017.square_plate.reg_clamped", "--json").stdout == before


def test_set_refuses_to_write_an_invalid_row(tmp_path: Path) -> None:
    store = _mixed_store(tmp_path)
    before = run(store, "get", "ko2017.square_plate.reg_clamped", "--json").stdout
    # margin_pct null while status is still 'measured' breaks a check invariant.
    completed = run(
        store, "set", "ko2017.square_plate.reg_clamped", "comparisons[0].measured.margin_pct=null"
    )
    assert completed.returncode == 2
    assert "refusing to write: the edited row is invalid" in completed.stderr
    assert run(store, "get", "ko2017.square_plate.reg_clamped", "--json").stdout == before


def test_set_flag_add_and_remove(tmp_path: Path) -> None:
    store = _mixed_store(tmp_path)
    completed = run(store, "set", "ko2017.square_plate.reg_clamped", "--add-flag", "tautology")
    assert completed.returncode == 0, completed.stdout + completed.stderr
    after = yaml.safe_load(run(store, "get", "ko2017.square_plate.reg_clamped", "--json").stdout)
    assert after["flags"] == ["tautology"]

    completed = run(store, "set", "ko2017.square_plate.reg_clamped", "--remove-flag", "tautology")
    assert completed.returncode == 0
    after = yaml.safe_load(run(store, "get", "ko2017.square_plate.reg_clamped", "--json").stdout)
    assert after["flags"] == []

    refused = run(store, "set", "ko2017.square_plate.reg_clamped", "--add-flag", "near")
    assert refused.returncode == 2
    assert "derived flag" in refused.stderr
    refused = run(store, "set", "ko2017.square_plate.reg_clamped", "--add-flag", "bogus")
    assert refused.returncode == 2
    assert "unknown flag id" in refused.stderr


def test_set_unset_removes_an_optional_field(tmp_path: Path) -> None:
    store = _mixed_store(tmp_path)
    completed = run(
        store, "set", "ko2017.square_plate.reg_clamped", "--unset", "comparisons[2].notes"
    )
    assert completed.returncode == 0, completed.stdout + completed.stderr
    after = yaml.safe_load(run(store, "get", "ko2017.square_plate.reg_clamped", "--json").stdout)
    assert "notes" not in after["comparisons"][2]


def test_real_store_answers_read_verbs() -> None:
    """The shipped store must be queryable, even with no rows yet."""
    assert run(REAL_STORE, "list").returncode == 0
    headline = run(REAL_STORE, "headline")
    assert headline.returncode == 0
    assert "| 3 |" in headline.stdout
    assert run(REAL_STORE, "find", "ko2017").returncode == 0


def test_set_schema_matches_the_validation_key_sets() -> None:
    """The write schema must not drift from the read-side key sets."""
    module = _load_tool_module()
    schema = module.SET_SCHEMA
    assert set(schema) == module.ROW_KEYS - {"id"}
    assert set(schema["comparisons"][0]) == module.COMPARISON_KEYS
    assert set(schema["result"]) == module.RESULT_KEYS
    assert set(schema["comparisons"][0]["reference"]) == module.REFERENCE_KEYS
    assert set(schema["comparisons"][0]["tolerance"]) == module.TOLERANCE_KEYS
    assert set(schema["comparisons"][0]["measured"]) == module.MEASURED_KEYS
    assert set(schema["history"][0]) == module.HISTORY_KEYS


# --------------------------------------------------------------------------- #
# T3a: the bibliography store
# --------------------------------------------------------------------------- #


def test_references_render_is_deterministic_and_checkable(tmp_path: Path) -> None:
    store = _mixed_store(tmp_path)
    target = tmp_path / "references.md"
    first = run(store, "references", "render", "--out", str(target))
    assert first.returncode == 0, first.stdout + first.stderr
    rendered = target.read_text(encoding="utf-8")

    assert rendered.startswith("# Canonical bibliography\n")
    assert "Do not edit by hand" in rendered
    assert "This file is the canonical bibliography for AeroElast." in rendered
    assert "## About the recovered PDFs" in rendered
    assert "Verification convention used below." in rendered
    assert "### " not in rendered
    assert "DOI: 10.1016/j.compstruc.2017.08.003." in rendered
    assert "DOI: to verify." in rendered
    assert "*Verified against `" in rendered

    unchanged = run(store, "references", "render", "--out", str(target), "--check")
    assert unchanged.returncode == 0, unchanged.stdout + unchanged.stderr

    target.write_text(rendered + "\ntampered\n", encoding="utf-8")
    stale = run(store, "references", "render", "--out", str(target), "--check")
    assert stale.returncode == 1
    assert "differs from a fresh render" in stale.stderr


def test_references_render_declares_undeclared_sections(tmp_path: Path) -> None:
    references = copy.deepcopy(REFERENCES)
    references["references"][0]["section"] = "9. A section with no intro"
    store = write_store(tmp_path, [], references=references)
    warned = run(store, "references", "check")
    assert "which `sections` does not declare" in warned.stdout
    rendered = run(store, "references", "render", "--out", str(tmp_path / "out.md"))
    assert rendered.returncode == 0
    assert "## 9. A section with no intro" in (tmp_path / "out.md").read_text(encoding="utf-8")


def test_references_render_emits_numeral_subsection_headings(tmp_path: Path) -> None:
    references = copy.deepcopy(REFERENCES)
    references["sections"].append({"title": "6.1 BEM theory"})
    references["references"][0]["section"] = "6.1 BEM theory"
    store = write_store(tmp_path, [], references=references)
    target = tmp_path / "references.md"
    assert run(store, "references", "render", "--out", str(target)).returncode == 0
    assert "### 6.1 BEM theory" in target.read_text(encoding="utf-8")


def test_param_tokens_survive_the_ids_pytest_actually_builds() -> None:
    """Ids join parameters with `-`, exponents carry one, and `In-plane` ends in `e`."""
    module = _load_tool_module()
    thick = "0.02667-In-plane-1.0-0.005424-0.001754-0.9971-0.01"
    assert module.param_tokens(thick, 7) == [
        "0.02667",
        "In-plane",
        "1.0",
        "0.005424",
        "0.001754",
        "0.9971",
        "0.01",
    ]
    assert module.param_numbers(thick, 7) == {0.02667, 1.0, 0.005424, 0.001754, 0.9971, 0.01}

    thin = "0.0002667-Out-of-plane-1e-06-0.005256-0.001294-0.9982-0.01"
    assert module.param_tokens(thin, 7) == [
        "0.0002667",
        "Out-of-plane",
        "1e-06",
        "0.005256",
        "0.001294",
        "0.9982",
        "0.01",
    ]
    assert module.param_numbers(thin, 7) == {0.0002667, 1e-06, 0.005256, 0.001294, 0.9982, 0.01}

    # The separator is not a sign, and a parameter *name* is not a number.
    table = "expected_table2_30-expected_table4_50-False-0.01-100.0"
    assert module.param_numbers(table, 5) == {0.01, 100.0}
    assert module.param_numbers("0.004-2.0-expected_mitc40-True", 4) == {0.004, 2.0}


def test_assertion_calls_exclude_geometric_tolerances() -> None:
    """`tol=` on a node-search helper is not an acceptance bound."""
    module = _load_tool_module()
    tree = ast.parse(
        (REPO_ROOT / "tests/validation/benchmarks/test_ko2017_performance.py").read_text(
            encoding="utf-8"
        )
    )
    consts = module.module_constants(tree)
    func = next(
        node
        for node in tree.body
        if isinstance(node, ast.FunctionDef)
        and node.name == "test_3_1_square_plate_tables_2_to_5"
    )
    texts = [ast.unparse(call) for call in module.assertion_calls(func)]
    assert any("isclose" in text for text in texts)
    assert not any("_find_node_by_xyz" in text for text in texts)
    assert not any("_clamped_edge_fixed_dofs" in text for text in texts)

    sites, ignored = module.tolerance_sites(func, consts, "0.01-100.0-False-expected_table2_30")
    assert {site.kind for site in sites} == {"rtol"}
    assert all(site.value == 0.05 for site in sites)
    # `atol=0.0` bounds nothing, so it is reported rather than stored.
    assert len(ignored) == 2
    assert all("atol=0" in entry for entry in ignored)


def test_extract_claims_every_collected_node() -> None:
    """T4 acceptance: the section 3 scope extracts to one row per collected node."""
    completed = run(REAL_STORE, "extract", "--json")
    assert completed.returncode == 0, completed.stdout + completed.stderr
    payload = json.loads(completed.stdout)
    assert payload["collected"] == 31
    assert payload["rows"] == 31
    assert payload["claimed"] == 31
    assert payload["unclaimed"] == []
    assert payload["diverged_params"] == {}


def test_reference_gaps_make_an_absence_explicit(tmp_path: Path) -> None:
    """A missing year is allowed only when declared, and software pins a version instead."""
    references = copy.deepcopy(REFERENCES)
    software = {
        "key": "openfoam",
        "section": "4. Software and vendored code",
        "kind": "software",
        "authors": [],
        "title": "OpenFOAM",
        "year": None,
        "version": "v2406",
        "doi": None,
        "doi_status": "not_applicable",
        "held": False,
        "bibliographic_gaps": ["authors", "year"],
        "notes": "the repository pins v2406; the source states no publication year",
    }
    references["references"].append(copy.deepcopy(software))
    store = write_store(tmp_path, [], references=references)
    completed = run_check(store)
    assert completed.returncode == 0, completed.stdout + completed.stderr

    # Same entry, without pinning a version: a software entry needs one identity.
    unpinned = copy.deepcopy(software)
    unpinned["version"] = None
    references2 = copy.deepcopy(REFERENCES)
    references2["references"].append(unpinned)
    store = write_store(tmp_path / "unpinned", [], references=references2)
    completed = run_check(store)
    assert completed.returncode == 1
    assert "must pin its identity with 'version'" in completed.stdout

    # Same entry, with the year absent but not declared: the absence must be stated.
    undeclared = copy.deepcopy(software)
    undeclared["bibliographic_gaps"] = ["authors"]
    references3 = copy.deepcopy(REFERENCES)
    references3["references"].append(undeclared)
    store = write_store(tmp_path / "undeclared", [], references=references3)
    completed = run_check(store)
    assert completed.returncode == 1
    assert "'year' is missing and is not declared" in completed.stdout


def test_cited_by_stale_is_visible_and_never_silently_true(tmp_path: Path) -> None:
    """A stale prose claim is recorded as a warning, and a live one as an error."""
    references = copy.deepcopy(REFERENCES)
    references["references"][0]["cited_by_stale"] = ["src/aeroelast/core/mesh/generators.py:1"]
    store = write_store(tmp_path, [], references=references)
    stale = run(store, "references", "check")
    assert stale.returncode == 0, stale.stdout + stale.stderr
    assert "is known stale" in stale.stdout

    # generators.py:65 does mention "Ko2017", so recording it as stale is wrong.
    references["references"][0]["cited_by_stale"] = ["src/aeroelast/core/mesh/generators.py:65"]
    store = write_store(tmp_path / "alive", [], references=references)
    alive = run(store, "references", "check")
    assert alive.returncode == 1
    assert "move it to cited_by_declared" in alive.stdout


def test_references_get_shows_the_doi(tmp_path: Path) -> None:
    store = _mixed_store(tmp_path)
    completed = run(store, "references", "get", "ko2017_perf")
    assert completed.returncode == 0, completed.stdout + completed.stderr
    entry = yaml.safe_load(completed.stdout)
    assert entry["doi"] == "10.1016/j.compstruc.2017.08.003"
    assert entry["doi_status"] == "verified"

    missing = run(store, "references", "get", "nope")
    assert missing.returncode == 2
    assert "unknown reference key" in missing.stderr


def test_references_list_gaps_and_filters(tmp_path: Path) -> None:
    store = _mixed_store(tmp_path)
    every = run(store, "references", "list")
    assert every.returncode == 0
    assert "2 entr(ies)" in every.stdout

    gaps = run(store, "references", "gaps")
    assert "escalera2023" in gaps.stdout
    assert "ko2017_perf" not in gaps.stdout

    unheld = run(store, "references", "list", "--not-held")
    assert "escalera2023" in unheld.stdout
    assert "ko2017_perf" not in unheld.stdout


def test_real_store_serves_the_section_3_key() -> None:
    """T3 acceptance: the shipped store answers with the printed DOI."""
    completed = run(REAL_STORE, "references", "get", "ko2017_perf")
    assert completed.returncode == 0, completed.stdout + completed.stderr
    entry = yaml.safe_load(completed.stdout)
    assert entry["doi"] == "10.1016/j.compstruc.2017.08.003"
    assert entry["doi_status"] == "verified"
    assert entry["held"] is True


def test_references_check_reports_a_stale_declared_site(tmp_path: Path) -> None:
    """A declared site that no longer mentions the work is a finding; a live one is not."""
    references = copy.deepcopy(REFERENCES)
    references["references"][0]["cited_by_declared"] = [
        "src/aeroelast/core/mesh/generators.py:65",
    ]
    store = write_store(tmp_path, [], references=references)
    clean = run(store, "references", "check")
    assert clean.returncode == 0, clean.stdout + clean.stderr
    assert "no longer mentions" not in clean.stdout

    references["references"][0]["cited_by_declared"] = ["src/aeroelast/core/mesh/generators.py:1"]
    store = write_store(tmp_path / "stale", [], references=references)
    stale = run(store, "references", "check")
    assert stale.returncode == 1
    assert "no longer mentions this work" in stale.stdout
    assert "generators.py:1" in stale.stdout


def test_references_check_flags_past_end_of_file(tmp_path: Path) -> None:
    references = copy.deepcopy(REFERENCES)
    references["references"][0]["cited_by_declared"] = [
        "src/aeroelast/core/mesh/generators.py:999999"
    ]
    store = write_store(tmp_path, [], references=references)
    completed = run(store, "references", "check")
    assert completed.returncode == 1
    assert "past end of file" in completed.stdout


def test_references_where_used_lists_rows_and_code(tmp_path: Path) -> None:
    store = _mixed_store(tmp_path)
    citing = run(store, "references", "where-used", "escalera2023")
    assert citing.returncode == 0
    assert "ko2017.square_plate.reg_clamped [wide]" in citing.stdout

    by_code = run(store, "references", "where-used", "ko2017_perf")
    assert by_code.returncode == 0
    assert "code | src/aeroelast/core/mesh/generators.py:65" in by_code.stdout


def test_references_bibtex_exports_the_doi(tmp_path: Path) -> None:
    store = _mixed_store(tmp_path)
    completed = run(store, "references", "bibtex")
    assert completed.returncode == 0
    assert "@article{ko2017_perf," in completed.stdout
    assert "doi = {10.1016/j.compstruc.2017.08.003}" in completed.stdout
    assert "@inproceedings{escalera2023," in completed.stdout


def test_rows_citing_an_unknown_key_are_reported(tmp_path: Path) -> None:
    references = copy.deepcopy(REFERENCES)
    references["references"] = [references["references"][1]]  # drop ko2017_perf
    store = write_store(tmp_path, [make_row()], references=references)
    completed = run_check(store)
    assert completed.returncode == 1
    assert "does not resolve in references.yaml" in completed.stdout


@pytest.mark.parametrize(
    ("mutation", "expected"),
    [
        ("verified_without_held", "verified requires held: true"),
        ("verified_without_note", "requires a 'verification_note'"),
        ("to_verify_with_doi", "to_verify requires doi: null"),
        ("not_applicable_without_notes", "not_applicable requires 'notes'"),
        ("bad_doi_format", "doi must be a bare DOI"),
        ("bad_doi_status", "doi_status must be one of"),
        ("declared_site_missing", "declared citation site does not exist"),
        ("docs_site_missing", "cited_by_docs: site does not exist"),
        ("docs_site_malformed", "cited_by_docs: not a file or file:line site"),
        ("recovered_site_missing", "recovered_from: site does not exist"),
        ("bad_year", "year must be an integer"),
        ("no_authors", "is not declared in bibliographic_gaps"),
        ("duplicate_key", "duplicate key"),
        ("unknown_bib_key", "unknown keys"),
    ],
)
def test_reference_entry_validation(tmp_path: Path, mutation: str, expected: str) -> None:
    references = copy.deepcopy(REFERENCES)
    entry = references["references"][0]
    if mutation == "verified_without_held":
        entry["held"] = False
    elif mutation == "verified_without_note":
        entry.pop("verification_note")
    elif mutation == "to_verify_with_doi":
        entry["doi_status"] = "to_verify"
    elif mutation == "not_applicable_without_notes":
        entry["doi_status"] = "not_applicable"
        entry["doi"] = None
    elif mutation == "bad_doi_format":
        entry["doi"] = "https://doi.org/10.1016/j.compstruc.2017.08.003"
    elif mutation == "bad_doi_status":
        entry["doi_status"] = "maybe"
    elif mutation == "declared_site_missing":
        entry["cited_by_declared"] = ["src/aeroelast/nope.py:1"]
    elif mutation == "docs_site_missing":
        entry["cited_by_docs"] = ["docs/nope.md:1"]
    elif mutation == "docs_site_malformed":
        entry["cited_by_docs"] = ["a site with spaces:1"]
    elif mutation == "recovered_site_missing":
        entry["recovered_from"] = ["docs/gone.md"]
    elif mutation == "bad_year":
        entry["year"] = "2017"
    elif mutation == "no_authors":
        entry["authors"] = []
    elif mutation == "duplicate_key":
        references["references"].append(copy.deepcopy(entry))
    elif mutation == "unknown_bib_key":
        entry["surprise"] = True
    store = write_store(tmp_path, [], references=references)
    completed = run_check(store)
    assert completed.returncode == 1, completed.stdout + completed.stderr
    assert expected in completed.stdout


# --------------------------------------------------------------------------- #
# Regression: the store is the baseline
# --------------------------------------------------------------------------- #

SYNTHETIC_PYTEST = """\
collecting ... collected 1 item

tests/validation/benchmarks/test_ko2017_performance.py::test_3_6_hook_table_14_minimal_fix[0.9782] DEBUG: Running hook_MITC4
  Nodes: 441, Elements: 384
  Norm vs Kirchhoff: 0.9814 (expected: 0.9782, error: 0.33%)
Hook MITC4: normalized = 0.9814 (expected 0.9782)
PASSED

============================== 1 passed in 5.32s ===============================
"""
HOOK_NODE = (
    "tests/validation/benchmarks/test_ko2017_performance.py::"
    "test_3_6_hook_table_14_minimal_fix[0.9782]"
)


def test_parse_prints_reads_the_shape_pytest_actually_emits() -> None:
    """pytest puts the node id and the test's first print on the same line."""
    module = _load_tool_module()
    prints = module.parse_prints(SYNTHETIC_PYTEST)
    assert list(prints) == [HOOK_NODE]
    assert "Norm vs Kirchhoff: 0.9814 (expected: 0.9782, error: 0.33%)" in prints[HOOK_NODE]
    # The status word and the summary belong to the run, not to the node's output.
    assert "PASSED" not in prints[HOOK_NODE]
    assert not any("1 passed" in line for line in prints[HOOK_NODE])


def test_extract_residuals_uses_the_declared_group_patterns() -> None:
    module = _load_tool_module()
    patterns = json.loads(
        (REAL_STORE / module.RESIDUAL_FILE).read_text(encoding="utf-8")
    )
    residual = patterns["patterns"]["3"]
    asserted, unasserted = module.extract_residuals(
        [
            "Norm vs Kirchhoff: 0.9814 (expected: 0.9782, error: 0.33%)",
            "[x] Norm vs Paper 3D: 6249.0122 (expected: 0.0, error: 622932.13%)",
            "Nodes: 441, Elements: 384",
        ],
        residual,
    )
    assert asserted == [{"value": "0.9814", "expected": "0.9782", "error": "0.33"}]
    assert unasserted == [{"value": "6249.0122", "expected": "0.0", "error": "622932.13"}]


def _row_ref(comparisons: list[dict[str, Any]]) -> Any:
    module = _load_tool_module()
    return module.RowRef(
        id="ko2017.toy.case",
        where="toy",
        path=None,
        index=0,
        data={"id": "ko2017.toy.case", "tests": [], "comparisons": comparisons},
    )


def _asserted(baseline: float | None) -> dict[str, Any]:
    measured: dict[str, Any] = {"status": "not_measured", "margin_pct": None}
    if baseline is not None:
        measured = {"status": "measured", "margin_pct": baseline}
    return {
        "label": "rtol",
        "asserted": True,
        "reference": {"kind": "analytical", "label": "x"},
        "tolerance": {"kind": "rtol", "value": 0.05, "source": "s", "justified": True},
        "measured": measured,
    }


def _printed(error: str, value: str = "1.0000", expected: str = "1.0000") -> dict[str, str]:
    return {"value": value, "expected": expected, "error": error}


def test_a_group_slug_always_produces_a_valid_row_id() -> None:
    """The two grammars have to agree: a row id starts with its group's slug.

    They did not agree. `SLUG_RE` accepts an underscore and `ID_RE` did not accept one in the id's
    first segment, so a group called `orthotropic_shell_parity` produced three row ids that no
    `set` could repair -- every write was refused for editing an invalid row.
    """
    module = _load_tool_module()
    for slug in ("ko2017", "outofband", "orthotropic_shell_parity", "a1_b2"):
        assert module.SLUG_RE.match(slug), slug
        assert module.ID_RE.match(f"{slug}.test_something.single"), slug


def test_an_unclaimed_node_needs_a_declaration_to_stay_silent() -> None:
    """The closure rule: a grouped file cannot gain a test that is neither row nor declaration."""
    module = _load_tool_module()
    group = {
        "id": "4",
        "non_validation_tests": [
            {
                "tests": [
                    "tests/validation/element/test_quad.py::TestStiffnessMatrix",
                    "tests/validation/element/test_quad.py::test_api_shape",
                ],
                "reason": "structural properties and API shape, not validation",
            }
        ],
    }
    nodes = [
        "tests/validation/element/test_quad.py::TestStiffnessMatrix::test_symmetry[Quad4]",
        "tests/validation/element/test_quad.py::test_api_shape",
        "tests/validation/element/test_quad.py::test_something_new",
    ]
    declared, undeclared, stale = module.classify_unclaimed(group, nodes)
    assert len(declared) == 2  # the class prefix and the exact id
    assert undeclared == ["tests/validation/element/test_quad.py::test_something_new"]
    assert stale == []

    # A declaration that matches nothing reads like coverage and is not.
    _, _, stale = module.classify_unclaimed(group, nodes[:2])
    assert stale == []  # both patterns were used
    _, _, stale = module.classify_unclaimed(group, [])
    assert len(stale) == 2

    # A prefix matches at a `::` boundary only: `test_api` must not shadow `test_api_shape2`.
    assert module.match_non_validation(group, "tests/x.py::TestStiffnessMatrix2::t") is None


def test_a_class_based_node_id_is_understood() -> None:
    """pytest prints `file::Class::test[param]` for a method.

    Reading only the first segment after the file silently dropped every class-based test and
    made the extractor report "no nodes" for the file; 15 of the 42 validation files are written
    that way, so the store could not see a third of the suite.
    """
    module = _load_tool_module()
    info = module.parse_node_id(
        "tests/validation/element/test_quad_elements.py::TestStiffnessMatrix::test_symmetry[Quad4]"
    )
    assert info is not None
    assert info.file == "tests/validation/element/test_quad_elements.py"
    assert info.owner == "TestStiffnessMatrix"
    assert info.function == "test_symmetry"
    assert info.qualname == "TestStiffnessMatrix::test_symmetry"
    assert info.params == "Quad4"
    assert info.node.endswith("TestStiffnessMatrix::test_symmetry[Quad4]")


def test_a_module_level_node_id_is_understood() -> None:
    module = _load_tool_module()
    info = module.parse_node_id("tests/validation/parity/test_a.py::test_case[3]")
    assert info is not None
    assert info.owner is None
    assert info.function == "test_case"
    assert info.qualname == "test_case"
    assert info.params == "3"
    assert module.parse_node_id("not-a-node-id") is None
    # No test name to take. Which `tests/` subtree a row may claim is `check`'s rule, not the
    # shape of a node id.
    assert module.parse_node_id("tests/other/test_a.py::") is None


def test_a_printed_expectation_is_stored_as_a_number() -> None:
    r"""Every residual pattern captures with `\S+`, so the print always yields a string."""
    module = _load_tool_module()
    assert module.scalar("0.9978") == 0.9978
    assert isinstance(module.scalar("0.9978"), float)
    # A reference side that is not one number stays text, rather than becoming a wrong float.
    assert module.scalar("1.1017/1.0323/1.0075") == "1.1017/1.0323/1.0075"
    assert module.scalar(None) is None


def test_re_extraction_keeps_what_was_measured_and_flagged(tmp_path: Path) -> None:
    """A better label must not cost the baselines, the digests or a human flag."""
    module = _load_tool_module()
    target = tmp_path / "rows.yaml"
    previous_row = {
        "id": "ko2017.toy.case",
        "group": "3",
        "title": "old",
        "tests": ["tests/a.py::b"],
        "validates": "old",
        "flags": ["out_of_band"],
        "history": [{"rev": "abc1234", "note": "earlier"}],
        "comparisons": [
            {
                "label": "rtol",
                "asserted": True,
                "reference": {"kind": "analytical", "label": "old"},
                "tolerance": {
                    "kind": "rtol",
                    "value": 0.05,
                    "source": "s",
                    "justified": True,
                },
                "measured": {"status": "measured", "margin_pct": 0.1, "digest": "abc"},
                "expected": 0.9978,
            }
        ],
    }
    target.write_text(
        yaml.safe_dump({"group": "3", "rows": [previous_row]}), encoding="utf-8"
    )
    new_rows = [
        {
            "id": "ko2017.toy.case",
            "group": "3",
            "title": "new",
            "tests": ["tests/a.py::b"],
            "validates": "better",
            "flags": [],
            "comparisons": [
                {
                    "label": "rtol",
                    "asserted": True,
                    "reference": {"kind": "analytical", "label": "better"},
                    "tolerance": {
                        "kind": "rtol",
                        "value": 0.05,
                        "source": "s",
                        "justified": True,
                    },
                    "measured": {"status": "not_measured", "margin_pct": None},
                }
            ],
        }
    ]
    carried = module.preserve_measurements(target, new_rows)
    assert carried == 2  # the measurement and the expectation
    assert new_rows[0]["flags"] == ["out_of_band"]
    assert new_rows[0]["history"][0]["rev"] == "abc1234"
    assert new_rows[0]["comparisons"][0]["measured"]["digest"] == "abc"
    assert new_rows[0]["comparisons"][0]["expected"] == 0.9978
    # Derived fields are refreshed, which is the whole point of re-extracting.
    assert new_rows[0]["validates"] == "better"
    assert new_rows[0]["comparisons"][0]["reference"]["label"] == "better"


def test_a_reordered_comparison_keeps_its_own_measurement(tmp_path: Path) -> None:
    """Comparisons are matched by their source, not by their position.

    Position was the key, and it was dangerous: dropping or reordering one comparison slid every
    later measurement one slot over, attaching a margin to a reference it was never measured
    against, silently. Here the two comparisons swap places between the stored row and the
    re-extracted one, so each measurement has to follow its own `tolerance.source`.
    """
    module = _load_tool_module()

    def comparison(source: str, margin: float | None) -> dict:
        return {
            "label": f"rel_err at {source}",
            "asserted": True,
            "reference": {"kind": "paper", "label": "x", "citation": "k"},
            "tolerance": {
                "kind": "rel_err",
                "value": 0.05,
                "source": source,
                "justified": True,
                "justification": "j",
            },
            "measured": {"status": "measured", "margin_pct": margin, "digest": None},
        }

    stored = comparison("tests/a.py:10", 1.0)
    stored["measured"] = {"status": "measured", "margin_pct": 1.0, "digest": "ten"}
    other = comparison("tests/a.py:20", 2.0)
    other["measured"] = {"status": "measured", "margin_pct": 2.0, "digest": "twenty"}
    target = tmp_path / "rows.yaml"
    target.write_text(
        yaml.safe_dump(
            {
                "rows": [
                    {
                        "id": "g.t.s",
                        "group": "g",
                        "title": "t",
                        "tests": ["tests/a.py::b"],
                        "validates": "v",
                        "comparisons": [stored, other],
                    }
                ]
            }
        ),
        encoding="utf-8",
    )
    swapped = comparison("tests/a.py:20", None)
    swapped["measured"] = {"status": "not_measured", "margin_pct": None}
    moved = comparison("tests/a.py:10", None)
    moved["measured"] = {"status": "not_measured", "margin_pct": None}
    new_rows = [
        {
            "id": "g.t.s",
            "group": "g",
            "title": "t",
            "tests": ["tests/a.py::b"],
            "validates": "v",
            "comparisons": [swapped, moved],
        }
    ]
    module.preserve_measurements(target, new_rows)
    first, second = new_rows[0]["comparisons"]
    assert first["measured"]["digest"] == "twenty"  # :20 keeps its own
    assert second["measured"]["digest"] == "ten"  # :10 keeps its own


def test_a_movement_too_small_for_a_tolerance_is_still_detected() -> None:
    """The digest is exact, so the detector needs no tolerance and hides nothing.

    A tolerance in the detector would swallow precisely the small movements a regression check
    exists to catch; the tolerance belongs in the test's own assertion.
    """
    module = _load_tool_module()
    ref = _row_ref([_asserted(0.1000)])
    baseline = module.compare_row(ref, [_printed("0.1000")], [])[0]
    comparison = ref.data["comparisons"][0]
    comparison["measured"]["digest"] = module.evidence_digest(_printed("0.1000"))

    moved = module.compare_row(ref, [_printed("0.1004")], [])[0]
    assert moved["verdict"] == "changed", "a 0.004 point movement must not go unnoticed"
    assert abs(moved["delta"]) < 0.005  # and it is smaller than the tolerance I used before
    assert baseline["digest"] != moved["digest"]


def test_compare_row_classifies_every_outcome() -> None:
    module = _load_tool_module()
    ref = _row_ref([_asserted(0.10), _asserted(None)])
    # A stored digest is what makes a comparison comparable at all: without one the verdict is
    # new_baseline, which is what the first --write records.
    ref.data["comparisons"][0]["measured"]["digest"] = module.evidence_digest(_printed("0.10"))

    same = module.compare_row(ref, [_printed("0.10"), _printed("0.20")], [])
    assert [item["verdict"] for item in same] == ["same", "new_baseline"]

    drifted = module.compare_row(ref, [_printed("0.90"), _printed("0.20")], [])
    assert drifted[0]["verdict"] == "changed"
    assert drifted[0]["baseline"] == 0.10
    assert drifted[0]["current"] == 0.90
    assert drifted[0]["delta"] == 0.8

    # A count mismatch is reported, never guessed: attaching a margin to the wrong
    # comparison would manufacture a baseline.
    unmapped = module.compare_row(ref, [_printed("0.10")], [])
    assert len(unmapped) == 1
    assert unmapped[0]["verdict"] == "unmapped"

    informational = module.compare_row(ref, [_printed("0.10"), _printed("0.20")], [_printed("622932.13")])
    assert informational[-1]["verdict"] == "informational"


def write_residual_patterns(store: Path, patterns: dict[str, Any]) -> None:
    (store / "residual-patterns.json").write_text(
        json.dumps({"version": 1, "patterns": patterns}), encoding="utf-8"
    )


def test_residual_patterns_must_be_valid_regexes(tmp_path: Path) -> None:
    store = write_store(tmp_path, [make_row()])
    write_residual_patterns(store, {"3": {"asserted": "([unclosed"}})
    completed = run_check(store)
    assert completed.returncode == 1
    assert "not a valid regex" in completed.stdout


def test_a_pattern_for_an_unknown_group_is_reported(tmp_path: Path) -> None:
    store = write_store(tmp_path, [make_row()])
    write_residual_patterns(store, {"99": {"asserted": "x"}})
    completed = run_check(store)
    assert completed.returncode == 1
    assert "unknown group" in completed.stdout


def test_regression_needs_a_declared_pattern(tmp_path: Path) -> None:
    store = write_store(tmp_path, [make_row()])
    completed = run(store, "regression", "--group", "3")
    assert completed.returncode == 2
    assert "declares no residual patterns" in completed.stderr


# --------------------------------------------------------------------------- #
# T5: cross-validating the Markdown view
# --------------------------------------------------------------------------- #

def write_adjudications(store: Path, entries: list[dict[str, Any]], group: str = "3") -> None:
    directory = store / "adjudications"
    directory.mkdir(exist_ok=True)
    (directory / "3-ko2017.yaml").write_text(
        yaml.safe_dump({"group": group, "adjudications": entries}), encoding="utf-8"
    )


def adjudication(**overrides: Any) -> dict[str, Any]:
    entry: dict[str, Any] = {
        "id": "ADJ-0001",
        "row": "ko2017.square_plate.reg_clamped",
        "field": "tolerance.kind",
        "verdict": "the same bound under two names",
        "resolution": "the extractor maps tol to rtol",
        "decided_by": "maintainer",
        "date": "2026-10-03",
    }
    entry.update(overrides)
    return entry


def test_an_adjudication_record_is_validated_not_trusted(tmp_path: Path) -> None:
    """The audit trail is data: an entry naming a dead row is an error, not a note."""
    store = write_store(tmp_path, [make_row()])
    write_adjudications(store, [adjudication()])
    assert run_check(store).returncode == 0

    write_adjudications(store, [adjudication(row="ko2017.gone.away")])
    completed = run_check(store)
    assert completed.returncode == 1
    assert "references a row that does not exist" in completed.stdout


def test_the_real_adjudication_record_passes() -> None:
    """The shipped audit trail must satisfy its own validator."""
    assert (REAL_STORE / "adjudications" / "3-ko2017.yaml").exists()
    assert run_check(REAL_STORE).returncode == 0


@pytest.mark.parametrize(
    ("entry", "expected"),
    [
        (adjudication(id="ADJ-1"), "invalid id"),
        (adjudication(status="maybe"), "status must be one of"),
        (adjudication(verdict=""), "'verdict' must be a non-empty string"),
        (adjudication(resolution=""), "'resolution' must be a non-empty string"),
        (adjudication(decided_by=""), "'decided_by' must be a non-empty string"),
        (adjudication(date="03/10/2026"), "date must be YYYY-MM-DD"),
        (adjudication(surprise=True), "unknown keys"),
    ],
)
def test_a_malformed_adjudication_is_reported(
    tmp_path: Path, entry: dict[str, Any], expected: str
) -> None:
    store = write_store(tmp_path, [make_row()])
    write_adjudications(store, [entry])
    completed = run_check(store)
    assert completed.returncode == 1, completed.stdout + completed.stderr
    assert expected in completed.stdout


def test_duplicate_adjudication_ids_are_reported(tmp_path: Path) -> None:
    store = write_store(tmp_path, [make_row()])
    write_adjudications(store, [adjudication(), adjudication()])
    completed = run_check(store)
    assert completed.returncode == 1
    assert "duplicate adjudication id" in completed.stdout


# --------------------------------------------------------------------------- #
# Gaps: what the suite does not validate
# --------------------------------------------------------------------------- #


def write_gaps(store: Path, entries: list[dict[str, Any]]) -> None:
    (store / "gaps.yaml").write_text(
        yaml.safe_dump({"version": 1, "gaps": entries}), encoding="utf-8"
    )


def gap(**overrides: Any) -> dict[str, Any]:
    entry: dict[str, Any] = {
        "id": "composite_stress_recovery",
        "scope": "composite outer-fibre stress recovery",
        "status": "not_validated",
        "reason": "CCX ignores OUTPUT=3D for composite sections",
        "evidence": "issue #3",
        "consequence": "only the ABD matrices may be cited",
        "citations_forbidden": True,
        "rows": [],
    }
    entry.update(overrides)
    return entry


def test_the_shipped_gap_list_is_present_and_valid() -> None:
    """The not-citable list must survive the Markdown it came from.

    The count is read out of the store rather than pinned here: a number written into this test
    churns on every gap the maintainer adds, and the pinned one went stale the moment a defect was
    recorded. What is worth asserting is that every entry the store holds reaches the output, and
    that the list is not empty.
    """
    completed = run(REAL_STORE, "gaps")
    assert completed.returncode == 0, completed.stdout + completed.stderr
    entries = yaml.safe_load((REAL_STORE / "gaps.yaml").read_text(encoding="utf-8"))["gaps"]
    assert entries, "the not-citable list is empty"
    assert f"gaps: {len(entries)}" in completed.stdout
    for entry in entries:
        assert entry["id"] in completed.stdout
    # The marker tracks the entrys own field, so the invariant is the relation between them and
    # not a number that someone has to remember to bump.
    not_citable = [entry for entry in entries if entry.get("citations_forbidden")]
    assert not_citable, "no entry is marked not citable"
    assert completed.stdout.count("NOT CITABLE") == len(not_citable)


@pytest.mark.parametrize(
    ("entry", "expected"),
    [
        (gap(id="hyphen-not-allowed"), "invalid id"),
        (gap(status="maybe"), "status must be one of"),
        (gap(scope=""), "'scope' must be a non-empty string"),
        (gap(reason=""), "'reason' must be a non-empty string"),
        (gap(evidence=""), "'evidence' must be a non-empty string"),
        (gap(citations_forbidden="yes"), "citations_forbidden must be a boolean"),
        (gap(rows=["ko2017.gone.away"]), "references a row that does not exist"),
        (gap(surprise=True), "unknown keys"),
    ],
)
def test_a_malformed_gap_is_reported(
    tmp_path: Path, entry: dict[str, Any], expected: str
) -> None:
    store = write_store(tmp_path, [make_row()])
    write_gaps(store, [entry])
    completed = run_check(store)
    assert completed.returncode == 1, completed.stdout + completed.stderr
    assert expected in completed.stdout


def test_a_document_that_cites_the_work_is_declared_and_validated(tmp_path: Path) -> None:
    """Deleting or moving the document must break check, not leave a note that rots."""
    references = copy.deepcopy(REFERENCES)
    entry = references["references"][0]
    entry["cited_by_docs"] = ["docs/validation-policy.md"]
    entry["recovered_from"] = ["README.md", "docs/reading-sources.md"]
    store = write_store(tmp_path, [], references=references)
    assert run_check(store).returncode == 0

    entry["cited_by_docs"] = ["docs/deleted.md"]
    store = write_store(tmp_path / "moved", [], references=references)
    completed = run_check(store)
    assert completed.returncode == 1
    assert "cited_by_docs: site does not exist: docs/deleted.md" in completed.stdout


def test_a_gap_may_point_at_the_rows_that_evidence_it(tmp_path: Path) -> None:
    store = write_store(tmp_path, [make_row()])
    write_gaps(store, [gap(rows=["ko2017.square_plate.reg_clamped"])])
    assert run_check(store).returncode == 0


def test_status_reports_coverage_and_pending_migration(tmp_path: Path) -> None:
    store = write_store(tmp_path, [make_row()])
    completed = run(store, "status", "--json")
    assert completed.returncode == 0, completed.stdout + completed.stderr
    payload = json.loads(completed.stdout)
    assert payload["groups"][0]["group"] == "3"
    assert payload["groups"][0]["rows"] == 1
    assert isinstance(payload["ungrouped_validation_files"], list)


def test_the_real_store_status_names_the_pending_migration() -> None:
    """`check` cannot see a validation file with no rows; this is what makes it visible."""
    completed = run(REAL_STORE, "status", "--json")
    assert completed.returncode == 0, completed.stdout + completed.stderr
    payload = json.loads(completed.stdout)
    groups = {row["group"]: row for row in payload["groups"]}
    assert groups["3"]["rows"] == 31
    assert groups["10"]["rows"] == 1
    assert len(payload["ungrouped_validation_files"]) > 0


def test_duplicate_gap_ids_are_reported(tmp_path: Path) -> None:
    store = write_store(tmp_path, [make_row()])
    write_gaps(store, [gap(), gap()])
    completed = run_check(store)
    assert completed.returncode == 1
    assert "duplicate gap id" in completed.stdout
