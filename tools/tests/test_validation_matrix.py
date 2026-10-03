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
import subprocess
import sys
from pathlib import Path
from typing import Any

import pytest
import yaml

REPO_ROOT = Path(__file__).resolve().parents[2]
TOOL = REPO_ROOT / "tools" / "validation_matrix.py"
REAL_STORE = REPO_ROOT / "docs" / "validation"

GROUP_SOURCE = "tests/test_ko2017_performance.py"
NODE_ID = (
    "tests/test_ko2017_performance.py::test_3_1_square_plate_tables_2_to_5"
    "[expected_table2_30-expected_table4_50-False-0.01-100.0]"
)

GROUPS = {
    "version": 1,
    "groups": [
        {
            "id": "3",
            "title": "`tests/test_ko2017_performance.py`",
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
                "expected": "0.9984",
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
) -> Path:
    store = tmp_path / "validation"
    (store / "rows").mkdir(parents=True)
    (store / "groups.yaml").write_text(yaml.safe_dump(groups or GROUPS), encoding="utf-8")
    (store / "flags.yaml").write_text(yaml.safe_dump(flags or FLAGS), encoding="utf-8")
    if rows:
        (store / "rows" / "3-ko2017.yaml").write_text(
            yaml.safe_dump({"group": "3", "rows": rows}), encoding="utf-8"
        )
    return store


def run_check(store: Path, *args: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, str(TOOL), "--store", str(store), "check", *args],
        capture_output=True,
        text=True,
        check=False,
    )


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
