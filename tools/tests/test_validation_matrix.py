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
import shutil
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
                "tolerance": {
                    "kind": "rtol",
                    "value": 0.05,
                    "source": "same assert",
                    "justified": True,
                },
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
    assert (
        yaml.safe_load(run(store, "get", "ko2017.square_plate.reg_clamped", "--json").stdout)[
            "comparisons"
        ][2]["measured"]["margin_pct"]
        == 3.0
    )

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
    assert "the edit would leave an error the row did not have" in completed.stderr
    # The introduced error is named, so the operator does not have to hunt for what broke.
    assert "numeric margin_pct" in completed.stderr
    assert run(store, "get", "ko2017.square_plate.reg_clamped", "--json").stdout == before


def test_set_may_repair_a_row_that_is_already_imperfect(tmp_path: Path) -> None:
    """A write is refused for making a row worse, not for the row being imperfect.

    A row with two undeclared reference kinds cannot be repaired one kind at a time otherwise: fixing
    the first leaves the second undeclared, and refusing every write made `set` useless exactly where
    it is needed. The store-level `check` still reports what is left.
    """
    store = _mixed_store(tmp_path)
    broken = run(
        store, "set", "ko2017.square_plate.reg_clamped", "comparisons[0].measured.margin_pct=null"
    )
    assert broken.returncode == 2, broken.stdout + broken.stderr
    # The row is imperfect: `status: measured` with a null margin. One edit that leaves it invalid
    # for a reason the row did not have is refused, and the repair that resolves the reason it does
    # have is allowed -- both fields in one call, because either alone would introduce a new error.
    assert (
        run(
            store,
            "set",
            "ko2017.square_plate.reg_clamped",
            "comparisons[0].measured.status=not_measured",
        ).returncode
        == 2
    )
    fixed = run(
        store,
        "set",
        "ko2017.square_plate.reg_clamped",
        "comparisons[0].measured.status=not_measured",
        "--unset",
        "comparisons[0].measured.margin_pct",
    )
    assert fixed.returncode == 0, fixed.stdout + fixed.stderr


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
        if isinstance(node, ast.FunctionDef) and node.name == "test_3_1_square_plate_tables_2_to_5"
    )
    texts = [ast.unparse(call) for call in module.assertion_calls(func)]
    # The file's comparisons go through the suite's own helper, which is what the criterion reads.
    assert any("assert_residual_below" in text for text in texts)
    assert not any("_find_node_by_xyz" in text for text in texts)
    assert not any("_clamped_edge_fixed_dofs" in text for text in texts)

    sites, ignored = module.tolerance_sites(func, consts, "0.01-100.0-False-expected_table2_30")
    assert {site.kind for site in sites} == {"rtol"}
    assert all(site.value == 0.05 for site in sites)
    assert not ignored


def test_a_zero_bound_bounds_nothing() -> None:
    """An absolute bound of zero is reported as ignored, not stored as if it bounded anything.

    It used to be covered through a file that happened to pass atol=0.0 to a numpy helper.
    That file was converted, so the behaviour is pinned here on a source of its own rather
    than on the shape of a file that has already changed once.
    """
    module = _load_tool_module()
    tree = ast.parse("def test_x():\n    assert_allclose(value, reference, rtol=1e-6, atol=0.0)\n")
    func = next(node for node in tree.body if isinstance(node, ast.FunctionDef))
    sites, ignored = module.tolerance_sites(func, module.module_constants(tree), None)
    assert [site.kind for site in sites] == ["rtol"]
    assert len(ignored) == 1
    assert "atol=0" in ignored[0]


def test_the_site_criterion_reads_canonical_calls_only() -> None:
    """A bare relational assert is not a site: it carries a bound and no reference.

    A relational assert says a value came in below a bound; it does not say what
    the value was measured against, so the store has nothing to attribute. The
    comparison becomes visible when it is written through the helper that names
    its reference and its kind, and not before.
    """
    module = _load_tool_module()
    source = "\n".join(
        [
            "def test_x():",
            "    err = compute()",
            "    assert err < 0.05",
            "    assert_relative_error(",
            "        err, ref, tol=0.05,",
            '        reference_name="r", kind="analytical", what="w",',
            "    )",
        ]
    )
    tree = ast.parse(source)
    func = next(node for node in tree.body if isinstance(node, ast.FunctionDef))
    sites, ignored = module.tolerance_sites(func, module.module_constants(tree), None)
    assert [site.kind for site in sites] == ["rtol"]
    assert not ignored


def test_extract_claims_every_collected_node() -> None:
    """T4 acceptance: the section 3 scope extracts to one row per collected node."""
    completed = run(REAL_STORE, "extract", "--group", "3", "--json")
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
    patterns = json.loads((REAL_STORE / module.RESIDUAL_FILE).read_text(encoding="utf-8"))
    residual = patterns["patterns"]["3"]
    asserted, unasserted = module.extract_residuals(
        [
            "Norm vs Kirchhoff: 0.9814 (expected: 0.9782, error: 0.33%)",
            "[x] Norm vs Paper 3D: 6249.0122 (expected: 0.0, error: 622932.13%)",
            "Nodes: 441, Elements: 384",
        ],
        residual,
    )
    assert asserted == [
        {
            "value": "0.9814",
            "expected": "0.9782",
            "error": "0.33",
            "line": "Norm vs Kirchhoff: 0.9814 (expected: 0.9782, error: 0.33%)",
        }
    ]
    assert unasserted == [
        {
            "value": "6249.0122",
            "expected": "0.0",
            "error": "622932.13",
            "line": "[x] Norm vs Paper 3D: 6249.0122 (expected: 0.0, error: 622932.13%)",
        }
    ]


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

    # A bare test name is enough, so a declaration does not have to spell a node id too long to
    # wrap. It matches the name after the last `::` and nothing else.
    named = {"non_validation_tests": [{"tests": ["test_api_shape"], "reason": "r"}]}
    assert module.match_non_validation(named, "tests/x.py::test_api_shape") is not None
    assert module.match_non_validation(named, "tests/x.py::TestA::test_api_shape[1]") is not None
    assert module.match_non_validation(named, "tests/x.py::test_api_shape2") is None


def test_regression_does_not_count_a_declared_non_validation_test_as_unclaimed() -> None:
    """`regression` must split its leftovers the way `extract` does.

    The command built its leftover list by comparing the collected nodes against the store's rows
    and never consulted `non_validation_tests`, so a group whose declarations are correct was
    reported as `unclaimed` - and `unclaimed` is in the failing set, so `regression` exited
    non-zero on exactly the groups that had done the declaring. Measured on the shipped store:
    `--group 30` reported its 29 non-validation tests and `--group 34` its one.
    """
    module = _load_tool_module()
    group = {
        "id": "30",
        "non_validation_tests": [
            {"tests": ["test_omega_provider"], "reason": "a property, not a comparison"}
        ],
    }
    nodes = [
        "tests/x.py::test_omega_provider",
        "tests/x.py::test_something_undeclared",
    ]
    verdicts = module.regression_leftover_verdicts(group, nodes)
    by_node = {item["detail"]: item["verdict"] for item in verdicts}
    assert by_node["tests/x.py::test_omega_provider"] == "declared_non_validation"
    assert by_node["tests/x.py::test_something_undeclared"] == "unclaimed"

    # The declared test is the group saying so on purpose, so it is not a regression finding.
    assert "declared_non_validation" not in module.REGRESSION_FAILING_VERDICTS
    assert "stale_declaration" not in module.REGRESSION_FAILING_VERDICTS

    # A declaration that matched nothing is surfaced, not left to read like coverage.
    stale = module.regression_leftover_verdicts(group, [])
    assert [item["verdict"] for item in stale] == ["stale_declaration"]
    assert "test_omega_provider" in stale[0]["detail"]


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
    target.write_text(yaml.safe_dump({"group": "3", "rows": [previous_row]}), encoding="utf-8")
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


PYTEST_FAILURE = """\
============================= test session starts =============================
collected 2 items

tests/validation/parity/test_x.py::test_isotropic PASSED                [ 50%]
tests/validation/parity/test_x.py::test_laminate FAILED                 [100%]

=================================== FAILURES ===================================
_____________________________ test_laminate __________________________________

    def test_laminate():
>       assert_residual_below(
E       AssertionError: laminate theta_z rate on the Bredt rate: 3.1200% against Bredt T L / GJ
E       is 3.1200%, above the 2.0000% bound (kind analytical)

tests/validation/parity/test_x.py:12: AssertionError
=========================== short test summary info ============================
FAILED tests/validation/parity/test_x.py::test_laminate - AssertionError: lami
SKIPPED [1] tests/validation/parity/test_y.py:36: PETSc not available
========================= 1 failed, 1 passed, 1 skipped in 3.10s ===============
"""


def test_triage_says_which_test_failed_and_by_how_much() -> None:
    """The report a caller acts on: the test, its reference, and the distance past the bound.

    These numbers come from the canonical assertion's message, and that message is part of the
    contract. Pinning the parse here means a change to the message fails this test instead of quietly
    turning the report into "something failed, somewhere".
    """
    module = _load_tool_module()
    report = module.triage_report(PYTEST_FAILURE)
    assert report["counts"] == {"passed": 1, "failed": 1}
    assert len(report["failed"]) == 1
    item = report["failed"][0]
    assert item["node"] == "tests/validation/parity/test_x.py::test_laminate"
    assert item["status"] == "FAILED"
    assert item["measured_pct"] == 3.12
    assert item["bound_pct"] == 2.0
    assert item["over_pp"] == 1.12
    assert item["times_bound"] == 1.56
    assert "Bredt" in item["message"]
    # A skipped module says why: in this suite a missing external tool skips a row rather than
    # failing it, so the reason is the answer.
    assert report["skipped"] == [
        {
            "status": "SKIPPED",
            "where": "tests/validation/parity/test_y.py:36",
            "reason": "PETSc not available",
        }
    ]


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

    informational = module.compare_row(
        ref, [_printed("0.10"), _printed("0.20")], [_printed("622932.13")]
    )
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


def test_capture_prints_runs_every_scope_argument(monkeypatch: Any) -> None:
    """A multi-file group is one pytest call with several paths, not one path with a space.

    `command_regression` joins the group's `source_files` with spaces, and `capture_prints`
    passed that string as a *single* argv element, so pytest looked for a file literally
    named `a.py b.py` and the group failed with "no node output captured". `triage` already
    splats its scope (`*scope.split()`); this is the same contract.
    """
    module = _load_tool_module()
    calls: list[list[str]] = []

    def fake_run(argv: list[str], **kwargs: Any) -> subprocess.CompletedProcess[str]:
        calls.append(list(argv))
        return subprocess.CompletedProcess(argv, 0, stdout="", stderr="")

    monkeypatch.setattr(module.subprocess, "run", fake_run)
    module.capture_prints("tests/a.py tests/b.py")

    assert calls[0][-2:] == ["tests/a.py", "tests/b.py"], calls[0]


def test_group_is_required_and_never_assumed(tmp_path: Path) -> None:
    """`--group` defaulted to 3, so `--scope <file>` ran an unrelated group's pattern.

    The message names the groups that declare the scope, because that is the one token the
    caller is missing; assuming a group is how `source digest: moved` got printed against a
    group the run had never read.
    """
    store = write_store(tmp_path, [make_row()])

    completed = run(store, "regression", "--scope", GROUP_SOURCE)
    assert completed.returncode == 2, completed.stdout + completed.stderr
    assert "--group is required" in completed.stderr
    assert GROUP_SOURCE in completed.stderr
    assert "declared by group 3" in completed.stderr

    extracted = run(store, "extract", "--scope", GROUP_SOURCE)
    assert extracted.returncode == 2, extracted.stdout + extracted.stderr
    assert "--group is required" in extracted.stderr


def test_group_is_required_even_when_no_scope_is_given(tmp_path: Path) -> None:
    """The default was the whole problem: a bare `regression` ran group 3 unasked."""
    store = write_store(tmp_path, [make_row()])
    completed = run(store, "regression")
    assert completed.returncode == 2, completed.stdout + completed.stderr
    assert "--group is required" in completed.stderr
    assert "has no default" in completed.stderr


#: Two real validation files, because a fixture cannot be one: `NODE_ID_RE` admits only node ids
#: under `tests/`, so a temporary file collects to nothing and the tool blames the environment.
#: These are the group 36 pair -- the multi-file group this defect actually broke -- and
#: collecting them is the cheapest of the candidates.
TWO_FILE_SOURCES = (
    "tests/validation/parity/test_composite_ply_stress_parity.py",
    "tests/validation/parity/test_composite_stress_ccx_parity.py",
)


def two_file_groups() -> dict[str, Any]:
    """One group over two files.

    `reference_kind` keeps the derived rows valid without the `non_validation_tests` and `drop`
    declarations the real group carries: a test the real group declares is not a row here.
    """
    return {
        "version": 1,
        "groups": [
            {
                "id": "3",
                "title": "two files",
                "slug": "ko2017",
                "reference_kind": "analytical",
                "source_files": list(TWO_FILE_SOURCES),
            }
        ],
    }


def test_extract_derives_every_declared_source_file(tmp_path: Path) -> None:
    """A group's scope is every file it declares, not only the first.

    `command_extract` read `source_files[0]`, and `--write` dumped just those rows, so
    re-deriving a multi-file group deleted the rows of its other files -- which is how group 36
    lost its two CalculiX-parity rows. `coherence` reads the same path, which is why it reported
    `the extraction itself did not finish cleanly` for a group whose second file declares tests.
    """
    store = write_store(tmp_path, [], groups=two_file_groups())

    completed = run(store, "extract", "--group", "3", "--json")
    assert completed.returncode == 0, completed.stdout + completed.stderr
    payload = json.loads(completed.stdout)
    assert payload["scope"].split() == list(TWO_FILE_SOURCES)
    assert payload["collected"] > 0
    assert payload["rows"] == payload["claimed"]
    assert isinstance(payload["dynamic_multiplicity"], dict)


def test_a_loop_body_assertion_becomes_one_comparison_per_iteration(tmp_path: Path) -> None:
    """The derived row for `for station in (...)`: four comparisons, in execution order.

    The measured file (`test_composite_ply_stress_parity.py`) prints one residual per iteration
    from two call sites, so the extraction reported `unmapped` and the store's gate refused to
    attach a baseline to either: 4 printed residuals against 2 comparisons.
    """
    store = write_store(tmp_path, [], groups=two_file_groups())

    completed = run(store, "extract", "--group", "3", "--write")
    assert completed.returncode == 0, completed.stdout + completed.stderr
    rows = yaml.safe_load((store / "rows" / "3-ko2017.yaml").read_text(encoding="utf-8"))["rows"]
    loop_row = next(row for row in rows if "bending_top_and_bottom_match_clt_at" in row["id"])

    labels = [comparison["label"] for comparison in loop_row["comparisons"]]
    assert len(labels) == 4
    # Execution order, and not one site expanded at a time: the same iteration's two sites sit
    # together (`[line A, line B, line A, line B]`, not `[A, A, B, B]`).
    numbers = [int(label.split("line ")[1].split(" ")[0]) for label in labels]
    assert numbers[0] == numbers[2] and numbers[1] == numbers[3] and numbers[0] != numbers[1]
    assert labels[0].endswith("(station=centre)") and labels[1].endswith("(station=centre)")
    assert labels[2].endswith("(station=three_quarter)")
    assert labels[3].endswith("(station=three_quarter)")
    assert len({comparison["tolerance"]["source"] for comparison in loop_row["comparisons"]}) == 2


def test_extract_write_keeps_the_other_source_file_rows(tmp_path: Path) -> None:
    """The write is the whole row file: a file left out of the derivation loses its rows."""
    store = write_store(tmp_path, [], groups=two_file_groups())

    completed = run(store, "extract", "--group", "3", "--write")
    assert completed.returncode == 0, completed.stdout + completed.stderr
    rows = yaml.safe_load((store / "rows" / "3-ko2017.yaml").read_text(encoding="utf-8"))["rows"]
    derived = {str(row["tests"][0]).split("::", 1)[0] for row in rows}
    assert derived == set(TWO_FILE_SOURCES)


# --------------------------------------------------------------------------- #
# The comparison is one printed residual, not one call site
# --------------------------------------------------------------------------- #

_TWO_SITES_IN_ONE_LOOP = """\
for station in ("centre", "three_quarter"):
    assert_relative_error(
        1.0, 2.0, tol=0.05, kind="analytical", reference_name="ref", what=f"{station} TOP"
    )
    assert_relative_error(
        3.0, 4.0, tol=0.05, kind="analytical", reference_name="ref", what=f"{station} BOTTOM"
    )
"""

_NESTED_LOOPS = """\
for station in ("centre", "three_quarter"):
    for comp in (0, 1):
        assert_relative_error(
            1.0, 2.0, tol=0.05, kind="analytical", reference_name="ref", what="k"
        )
"""

_COMPUTED_LOOP = """\
for step in steps():
    assert_relative_error(
        1.0, 2.0, tol=0.05, kind="analytical", reference_name="ref", what="x"
    )
"""


def _sites_of(body: str, module: Any = None) -> list[Any]:
    """The tolerance sites of a test function whose body is `body`."""
    module = module or _load_tool_module()
    source = 'def test_case() -> None:\n    """One case."""\n' + "".join(
        f"    {line}\n" if line.strip() else "\n" for line in body.splitlines()
    )
    tree = ast.parse(source)
    func = next(node for node in tree.body if isinstance(node, ast.FunctionDef))
    sites, _ = module.tolerance_sites(func, module.module_constants(tree), None)
    return list(sites)


def test_a_literal_loop_yields_one_comparison_per_execution() -> None:
    """The prints come out in execution order, so the comparisons must be emitted in it too.

    Expanding site by site (`for site: for iteration`) would give line 1, line 1, line 2, line 2,
    and `compare_row` pairs the Nth print with the Nth comparison -- so the second print would be
    attached to the wrong site's baseline. The walk has to re-enter the loop body once per
    iteration instead.
    """
    sites = _sites_of(_TWO_SITES_IN_ONE_LOOP)

    assert len(sites) == 4
    lines = [site.line for site in sites]
    assert lines[0] == lines[2] and lines[1] == lines[3] and lines[0] != lines[1]
    assert [site.binding["station"] for site in sites] == [
        "centre",
        "centre",
        "three_quarter",
        "three_quarter",
    ]


def test_nested_literal_loops_multiply_in_nesting_order() -> None:
    """One call in two loops is one comparison per pair, innermost varying fastest."""
    sites = _sites_of(_NESTED_LOOPS)

    assert len(sites) == 4
    assert [(site.binding["station"], site.binding["comp"]) for site in sites] == [
        ("centre", "0"),
        ("centre", "1"),
        ("three_quarter", "0"),
        ("three_quarter", "1"),
    ]


def test_a_computed_loop_is_flagged_not_guessed() -> None:
    """A sequence the code computes has no knowable multiplicity: one site, and it says so."""
    sites = _sites_of(_COMPUTED_LOOP)

    assert len(sites) == 1
    assert sites[0].dynamic is True
    assert sites[0].binding is None


def test_a_single_execution_site_carries_no_binding() -> None:
    """Every existing row derives byte-identically: no loop, no suffix, no change."""
    sites = _sites_of(
        "assert_relative_error(\n"
        '    1.0, 2.0, tol=0.05, kind="analytical", reference_name="ref", what="w"\n'
        ")\n"
    )

    assert len(sites) == 1
    assert sites[0].binding is None
    assert sites[0].dynamic is False


def test_a_geometric_tolerance_is_not_a_comparison() -> None:
    """`atol=` on a node search bounds nothing about a result: only assertion calls are read.

    The walk collects calls in execution order, so it needs the same filter `assertion_calls`
    applies -- without it `isclose(..., atol=1e-12)` landed in the store as a comparison with no
    reference, which is a `check` error.
    """
    sites = _sites_of(
        "tip = [n for n in nodes if isclose(float(n.z), length, atol=1e-12)]\n"
        "assert_relative_error(\n"
        '    1.0, 2.0, tol=0.05, kind="analytical", reference_name="ref", what="w"\n'
        ")\n"
    )

    assert len(sites) == 1


def test_preserve_measurements_carries_duplicate_sources_in_order(tmp_path: Path) -> None:
    """Two executions of one call share `tolerance.source`; each keeps its own measurement.

    Keying the on-disk comparisons by source collapsed them into a dict, so the last one
    overwrote the rest and every execution came back carrying the same baseline.
    """
    module = _load_tool_module()
    target = tmp_path / "rows" / "3-ko2017.yaml"
    target.parent.mkdir(parents=True)
    stored = make_row()
    fresh = make_row()
    for row in (stored, fresh):
        row["comparisons"][1]["tolerance"]["source"] = row["comparisons"][0]["tolerance"]["source"]
    stored["comparisons"][0]["measured"]["margin_pct"] = 0.10
    stored["comparisons"][1]["measured"]["margin_pct"] = 0.20
    module.dump_yaml(target, {"group": "3", "rows": [stored]})

    for comparison in fresh["comparisons"]:
        comparison["measured"] = {"status": "not_measured", "margin_pct": None}
    carried = module.preserve_measurements(target, [fresh])

    # One carried block per comparison, in order: the property under test is which one each
    # comparison gets, not how many fields came with it.
    assert carried >= 2
    assert [item["measured"]["margin_pct"] for item in fresh["comparisons"]] == [0.10, 0.20]


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
    # The marker tracks each entry's own field, so the invariant is the relation between them
    # and not a number that someone has to remember to bump.
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
def test_a_malformed_gap_is_reported(tmp_path: Path, entry: dict[str, Any], expected: str) -> None:
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


def test_measurement_text_never_says_none() -> None:
    """The stored text is built from what the pattern captured, not from what it did not.

    A canonical print carries an error and a bound and no separate value, and rendering the
    missing part wrote "None (0.0000%)" into the store as the readable half of a measurement.
    """
    module = _load_tool_module()

    assert module.measurement_text({"error": "0.0023", "bound": "1.0000"}) == (
        "0.0023% bound 1.0000%"
    )
    assert module.measurement_text({"value": "0.9996", "error": "0.08"}) == "0.9996 0.08%"
    assert "None" not in module.measurement_text({"error": "0.0023"})
    assert module.measurement_text({}) == "the pattern captured no number"
    # The suffix comes from the comparison: a percent for a relative residual, the declared
    # unit for an absolute one, and nothing when the comparison does not say.
    assert module.residual_suffix({"tolerance": {"kind": "rtol", "value": 0.05}}) == "%"
    assert module.residual_suffix({"tolerance": {"kind": "atol", "unit": "%"}}) == "%"
    assert module.residual_suffix({"tolerance": {"kind": "atol", "unit": "m"}}) == "m"
    assert module.residual_suffix({"tolerance": {"kind": "atol", "unit": None}}) == ""
    assert module.measurement_text({"error": "1e-10"}, "") == "1e-10"


def test_the_stored_text_is_the_line_the_test_printed() -> None:
    """The reader keeps the printed line, so a text never guesses a unit.

    Rebuilding the text from the parsed groups lost what the print wrote -- a percent sign on
    group 3s unasserted diagnostics, which the reader dropped while the print carried it. A
    raw line cannot drift from itself, and the digest is computed from the three named groups
    rather than from this key, so keeping it changes no stored digest.
    """
    module = _load_tool_module()
    lines = ["[x] Norm vs Paper 3D: 6.2475 (expected: 0.9984, error: 525.75%)"]
    residual = {
        "asserted": "(?P<nothing>zzz)",
        "unasserted": r"Norm\svs\sPaper\s3D:\s(?P<value>\S+)\s\(expected:\s(?P<expected>\S+),\serror:\s(?P<error>\S+)%\)",
    }
    _, unasserted = module.extract_residuals(lines, residual)

    assert unasserted[0]["line"] == lines[0]
    assert module.evidence_digest(unasserted[0]) == module.evidence_digest(
        {"value": "6.2475", "expected": "0.9984", "error": "525.75"}
    )


def test_a_measurement_that_is_not_numeric_is_refused(tmp_path: Path) -> None:
    """A margin that is not a number is not a measurement, and the writer says so.

    Recording it would store a null and leave check to reject the file afterwards, which is
    how a run that could not read a print ends up looking like a run that measured something.
    The group residual pattern not reading the print is the finding, and it should surface as
    one instead of as a broken row file.
    """
    module = _load_tool_module()
    store = module.load_store(write_store(tmp_path, [make_row()]))
    ref = next(item for item in store.rows if item.data.get("group") == "3")
    result = {
        "verdict": "new_baseline",
        "row": ref.id,
        "comparison": 0,
        "current": None,
        "value": None,
        "text": "None (0.0824%%)",
        "digest": "d",
    }

    with pytest.raises(module.StoreError, match="not a number"):
        module.write_measurement(store, ref, 0, result)


def test_a_comparison_with_no_counterpart_says_that(tmp_path: Path, capsys) -> None:
    """A moved line leaves a fresh comparison with no counterpart, and the message says so.

    The message used to claim that no stored measurement matched, which is a different
    claim: the condition is a source with no counterpart on disk, so it fired whether or not
    a measurement existed, and a row that had never been measured still told the reader a
    re-measure was due.
    """
    module = _load_tool_module()
    target = tmp_path / "rows" / "1-x.yaml"
    target.parent.mkdir(parents=True)
    module.dump_yaml(
        target,
        {
            "group": "1",
            "rows": [
                {
                    "id": "x",
                    "comparisons": [
                        {
                            "label": "rtol at line 10",
                            "tolerance": {"kind": "rtol", "value": 0.05, "source": "f.py:10"},
                            "measured": {"status": "measured", "margin_pct": 3.0},
                        }
                    ],
                }
            ],
        },
    )
    fresh = [
        {
            "id": "x",
            "comparisons": [
                {
                    "label": "rtol at line 13",
                    "tolerance": {"kind": "rtol", "value": 0.05, "source": "f.py:13"},
                }
            ],
        }
    ]

    carried = module.preserve_measurements(target, fresh)

    err = capsys.readouterr().err
    assert carried == 0
    assert "no counterpart" in err
    assert "re-measure is due" not in err


def test_coherence_names_a_row_file_whose_citations_moved(tmp_path: Path) -> None:
    """A row file cites the line each comparison sits on, and editing a test moves it.

    Nothing else sees that. `check` cannot read the code, and comparing how many rows a
    group extracts compares a count rather than the bytes of the file, which is the way a
    stale citation hid until git status showed it. This re-derives the group and names it.
    """
    store = tmp_path / "store"
    shutil.copytree(REAL_STORE, store)
    row = next((store / "rows").glob("4-*.yaml"))
    data = yaml.safe_load(row.read_text(encoding="utf-8"))
    data["rows"][0]["comparisons"][0]["tolerance"]["source"] = "tests/nowhere.py:1"
    row.write_text(yaml.safe_dump(data), encoding="utf-8")

    completed = run(store, "coherence", "--group", "4")

    assert completed.returncode == 1, completed.stdout + completed.stderr
    assert row.name in completed.stdout, completed.stdout
    assert "coherent" not in completed.stdout


def test_the_real_store_is_coherent(real_store_copy: Path) -> None:
    """A cheap spot check: one group re-derives to the rows on disk.

    It runs on a copy: `coherence` refreshes a stale file as it finds it, and a test must not
    write the shipped store (#32).
    """
    completed = run(real_store_copy, "coherence", "--group", "4")

    assert completed.returncode == 0, completed.stdout + completed.stderr
    assert "coherent" in completed.stdout


@pytest.mark.slow
def test_every_group_re_derives_to_the_rows_on_disk(real_store_copy: Path) -> None:
    """The guard: no row file in the store cites a line the code does not have.

    This is the whole sweep, and it is slow because it collects every validation file. It
    belongs with the tool's own tests rather than with the physics suite: a row file goes
    stale when a test is edited, so the person it protects is the one working on the store.
    Deselect it with -m "not slow".

    The sweep runs on a copy, because `coherence` refreshes every stale file it finds and a
    test must not write the shipped store (#32). Refreshing the real store stays the
    maintainer's explicit act: `python tools/validation_matrix.py coherence`.
    """
    completed = run(real_store_copy, "coherence")

    assert completed.returncode == 0, (
        "stale row files; refresh them with `python tools/validation_matrix.py coherence`\n"
        + completed.stdout
        + completed.stderr
    )
    assert "were stale" not in completed.stdout
    assert "re-derive to the rows on disk" in completed.stdout


def test_the_real_store_status_reports_no_pending_migration() -> None:
    """`check` cannot see a validation file with no rows; this is what makes it visible.

    The sweep is complete: every file under the validation tree is grouped, declared
    out of scope, or was deleted as software. The list is empty now and has to stay
    that way, so the test asserts the empty list rather than its non-emptiness, and
    its name changed with the fact.
    """
    completed = run(REAL_STORE, "status", "--json")
    assert completed.returncode == 0, completed.stdout + completed.stderr
    payload = json.loads(completed.stdout)
    groups = {row["group"]: row for row in payload["groups"]}
    assert groups["3"]["rows"] == 31
    assert groups["10"]["rows"] == 1
    assert payload["ungrouped_validation_files"] == []


def test_duplicate_gap_ids_are_reported(tmp_path: Path) -> None:
    store = write_store(tmp_path, [make_row()])
    write_gaps(store, [gap(), gap()])
    completed = run_check(store)
    assert completed.returncode == 1
    assert "duplicate gap id" in completed.stdout


# --------------------------------------------------------------------------- #
# The contract layer: levels.yaml and contract.yaml
# --------------------------------------------------------------------------- #

ROW_ID = "ko2017.square_plate.reg_clamped"


def _git(*args: str) -> str:
    completed = subprocess.run(
        ["git", *args], cwd=REPO_ROOT, capture_output=True, text=True, check=True
    )
    return completed.stdout.strip()


def write_contract(
    store: Path,
    verdicts: list[dict[str, Any]] | None = None,
    levels: list[dict[str, Any]] | None = None,
) -> None:
    level = {"id": "L1", "title": "elements", "scope": "the fixture group", "groups": ["3"]}
    (store / "levels.yaml").write_text(
        yaml.safe_dump({"version": 1, "ceiling": 0.05, "levels": levels or [level]}),
        encoding="utf-8",
    )
    (store / "contract.yaml").write_text(
        yaml.safe_dump({"version": 1, "verdicts": verdicts or []}), encoding="utf-8"
    )


def contract_report(store: Path) -> dict[str, Any]:
    completed = run(store, "contract", "--json")
    assert completed.returncode == 0, completed.stderr
    return json.loads(completed.stdout)


def trusted(**overrides: Any) -> dict[str, Any]:
    entry = {
        "target": ROW_ID,
        "class": "physics",
        "verdict": "trusted",
        "audited_rev": _git("rev-parse", "HEAD"),
    }
    entry.update(overrides)
    return entry


def test_contract_counts_every_row_unaudited_without_a_verdict(tmp_path: Path) -> None:
    store = write_store(tmp_path, [make_row()])
    write_contract(store)
    level = contract_report(store)["levels"][0]
    assert level["counts"] == {"unaudited": 1}
    assert level["closed"] is False


def test_contract_closes_a_level_whose_rows_are_all_trusted(tmp_path: Path) -> None:
    store = write_store(tmp_path, [make_row()])
    write_contract(store, [trusted()])
    level = contract_report(store)["levels"][0]
    assert level["counts"] == {"trusted": 1}
    assert level["closed"] is True


def test_contract_marks_a_verdict_stale_when_a_bound_file_changed(tmp_path: Path) -> None:
    # AGENTS.md was added after its creating commit's parent, so a verdict audited at that
    # parent and bound to it has at least one commit against it.
    added = _git("log", "--format=%H", "--diff-filter=A", "-n1", "--", "AGENTS.md")
    store = write_store(tmp_path, [make_row()])
    write_contract(
        store, [trusted(audited_rev=_git("rev-parse", f"{added}^"), files=["AGENTS.md"])]
    )
    entry = contract_report(store)["levels"][0]["entries"][0]
    assert entry["state"] == "stale"
    assert entry["stale_by"]


def test_contract_marks_an_unresolvable_revision_stale(tmp_path: Path) -> None:
    store = write_store(tmp_path, [make_row()])
    write_contract(store, [trusted(audited_rev="0000000")])
    assert contract_report(store)["levels"][0]["entries"][0]["state"] == "stale"


def _widened_row(value: float) -> dict[str, Any]:
    row = make_row()
    row["comparisons"][0]["tolerance"]["value"] = value
    # The store itself demands prose for a band over 5 %; extraction fills it in for every
    # real row, which is exactly why that prose is not a verdict.
    row["comparisons"][0]["tolerance"]["justification"] = "extracted from the code"
    return row


def test_contract_flags_a_physics_tolerance_over_the_ceiling(tmp_path: Path) -> None:
    store = write_store(tmp_path, [_widened_row(0.30)])
    write_contract(store, [trusted()])
    level = contract_report(store)["levels"][0]
    assert level["over_ceiling"] == 1
    assert level["closed"] is False


def test_contract_exception_covers_only_up_to_its_bound(tmp_path: Path) -> None:
    exception = {
        "bound": 0.10,
        "source": "two independent beam references disagree by 10 %",
        "justification": "the shell cannot be held tighter than its references agree",
        "fixed_before_measuring": True,
    }
    for value, covered in ((0.08, True), (0.30, False)):
        store = write_store(tmp_path / str(value), [_widened_row(value)])
        write_contract(store, [trusted(verdict="exception", exception=exception)])
        level = contract_report(store)["levels"][0]
        assert (level["over_ceiling"] == 0) is covered, value


def test_contract_ignores_the_ceiling_for_an_identity(tmp_path: Path) -> None:
    store = write_store(tmp_path, [_widened_row(0.30)])
    write_contract(store, [trusted(**{"class": "identity"})])
    assert contract_report(store)["levels"][0]["over_ceiling"] == 0


def test_contract_keeps_an_absent_anchor_open(tmp_path: Path) -> None:
    store = write_store(tmp_path, [make_row()])
    level = {
        "id": "L1",
        "title": "elements",
        "scope": "the fixture group",
        "groups": ["3"],
        "absent": [{"evidence": "tests/test_anchor.py", "reason": "not imported"}],
    }
    write_contract(store, [trusted()], levels=[level])
    level_report = contract_report(store)["levels"][0]
    assert level_report["counts"] == {"trusted": 1, "absent": 1}
    assert level_report["closed"] is False


@pytest.mark.parametrize(
    ("verdicts", "levels", "expected"),
    [
        (
            [
                {
                    "target": ROW_ID,
                    "class": "physics",
                    "verdict": "exception",
                    "audited_rev": "abcdef0",
                    "exception": {
                        "bound": 0.1,
                        "source": "s",
                        "justification": "j",
                        "fixed_before_measuring": False,
                    },
                }
            ],
            None,
            "fixed before measuring",
        ),
        ([{"target": ROW_ID, "class": "physics", "verdict": "trusted"}], None, "audited_rev"),
        (
            [{"target": "no.such.row", "class": "physics", "verdict": "undecided"}],
            None,
            "neither a row id",
        ),
        ([], [{"id": "L1", "title": "t", "scope": "s", "groups": ["99"]}], "unknown group"),
    ],
)
def test_check_validates_the_contract_layer(
    tmp_path: Path, verdicts: list[dict[str, Any]], levels: Any, expected: str
) -> None:
    store = write_store(tmp_path, [make_row()])
    write_contract(store, verdicts, levels=levels)
    completed = run_check(store)
    assert completed.returncode == 1
    assert expected in completed.stdout


def test_contract_command_writes_nothing(tmp_path: Path) -> None:
    store = write_store(tmp_path, [make_row()])
    write_contract(store, [trusted()])
    before = {path: path.read_bytes() for path in store.rglob("*") if path.is_file()}
    run(store, "contract", "--rows")
    after = {path: path.read_bytes() for path in store.rglob("*") if path.is_file()}
    assert before == after


def test_contract_verdict_expires_when_its_levels_code_changes(tmp_path: Path) -> None:
    # The fixture's test file did not change after AGENTS.md's parent commit, AGENTS.md did:
    # standing in for the element code, it must expire the verdict only through depends_on.
    added = _git("log", "--format=%H", "--diff-filter=A", "-n1", "--", "AGENTS.md")
    rev = _git("rev-parse", f"{added}^")
    assert not _git("log", "--format=%h", f"{rev}..HEAD", "--", NODE_ID.split("::")[0])
    level = {"id": "L1", "title": "elements", "scope": "the fixture group", "groups": ["3"]}
    for depends_on, expected in (([], "trusted"), (["AGENTS.md"], "stale")):
        store = write_store(tmp_path / expected, [make_row()])
        write_contract(
            store, [trusted(audited_rev=rev)], levels=[{**level, "depends_on": depends_on}]
        )
        assert contract_report(store)["levels"][0]["entries"][0]["state"] == expected


def test_contract_a_level_is_blocked_by_an_open_level_it_requires(tmp_path: Path) -> None:
    levels = [
        {"id": "L1", "title": "elements", "scope": "the fixture group", "groups": ["3"]},
        {
            "id": "L2",
            "title": "solver",
            "scope": "one anchor",
            "absent": [{"evidence": "tests/test_anchor.py::test_it", "reason": "outside"}],
        },
    ]
    anchor = trusted(
        target="tests/test_anchor.py::test_it",
        arbiters=[{"kind": "analytical", "label": "closed form"}],
    )
    for row_verdict, l2_closed in (("trusted", True), ("undecided", False)):
        store = write_store(tmp_path / row_verdict, [make_row()])
        write_contract(store, [trusted(verdict=row_verdict), anchor], levels=levels)
        l2 = contract_report(store)["levels"][1]
        assert l2["own_closed"] is True
        assert l2["closed"] is l2_closed
        assert l2["blocked_by"] == ([] if l2_closed else ["L1"])


def _code_only_row() -> dict[str, Any]:
    row = make_row()
    row["comparisons"] = [row["comparisons"][0]]
    row["comparisons"][0]["reference"] = {"kind": "code", "label": "CalculiX 2.20 S8R"}
    return row


def test_contract_a_single_numerical_arbiter_is_not_validation(tmp_path: Path) -> None:
    store = write_store(tmp_path, [_code_only_row()])
    write_contract(store, [trusted()])
    level = contract_report(store)["levels"][0]
    assert level["counts"] == {"unarbitrated": 1}
    assert level["closed"] is False


def test_contract_two_independent_numerical_arbiters_validate(tmp_path: Path) -> None:
    arbiters = [
        {"kind": "produced_numerical", "label": "CalculiX 2.20 S8R"},
        {"kind": "published_numerical", "label": "BeamDyn deck, IEA 15 MW"},
    ]
    store = write_store(tmp_path, [_code_only_row()])
    write_contract(store, [trusted(arbiters=arbiters)])
    assert contract_report(store)["levels"][0]["closed"] is True


def test_contract_experimental_data_arbitrates_alone(tmp_path: Path) -> None:
    store = write_store(tmp_path, [_code_only_row()])
    write_contract(store, [trusted(arbiters=[{"kind": "experimental", "label": "test rig"}])])
    assert contract_report(store)["levels"][0]["closed"] is True


@pytest.mark.parametrize(("status", "closed"), [("bounded", True), ("not_validated", False)])
def test_contract_only_an_unbounded_gap_blocks_its_level(
    tmp_path: Path, status: str, closed: bool
) -> None:
    store = write_store(tmp_path, [make_row()])
    entry = gap()
    entry["status"] = status
    write_gaps(store, [entry])
    level = {"id": "L1", "title": "elements", "scope": "s", "groups": ["3"], "gaps": [entry["id"]]}
    write_contract(store, [trusted()], levels=[level])
    assert contract_report(store)["levels"][0]["closed"] is closed
