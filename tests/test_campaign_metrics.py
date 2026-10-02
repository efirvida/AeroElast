"""Tests for the campaign before/after metrics tool (``tools/campaign_metrics.py``).

The tool is the only thing that turns the G4 re-run into a readable before/after,
so its name resolution and delta convention are pinned here.  The module is loaded
from its path because ``tools/`` is not an importable package.
"""

from __future__ import annotations

import csv
import importlib.util
from pathlib import Path

import pytest

MODULE_PATH = Path(__file__).resolve().parents[1] / "tools" / "campaign_metrics.py"
_spec = importlib.util.spec_from_file_location("campaign_metrics", MODULE_PATH)
campaign_metrics = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(campaign_metrics)


def _row(campaign: str, case: str, **metrics) -> dict:
    return {"campaign": campaign, "case": case, **metrics}


def _write_csv(path: Path, rows: list[dict]) -> None:
    fields = ["campaign", "case", *campaign_metrics.METRICS]
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def test_alias_maps_a_renamed_campaign():
    old = [_row("c_mitc3fix", "yaw_0")]
    new = [_row("c", "yaw_0")]
    pairs, unmatched = campaign_metrics.resolve_pairs(new, old, {"c": "c_mitc3fix"})
    assert [key for key, _, _ in pairs] == ["c_mitc3fix/yaw_0"]
    assert unmatched == []


def test_unrelated_campaign_names_do_not_match_by_case_alone():
    old = [_row("conv_results_official", "h_fine")]
    new = [_row("conv_rerun", "h_fine")]
    # "conv_rerun" does not prefix-match "conv_results_official"; the exact
    # ``case`` is not enough on its own, so this must stay unmatched.
    pairs, unmatched = campaign_metrics.resolve_pairs(new, old, {})
    assert pairs == []
    assert unmatched == ["conv_rerun/h_fine"]


def test_prefix_suffix_campaign_matches_when_unique():
    old = [_row("campaign_official", "h_fine")]
    new = [_row("campaign", "h_fine")]
    pairs, unmatched = campaign_metrics.resolve_pairs(new, old, {})
    assert [key for key, _, _ in pairs] == ["campaign_official/h_fine"]
    assert unmatched == []


def test_two_recorded_states_are_ambiguous_and_stay_unmatched():
    old = [_row("c_mitc3fix", "yaw_0"), _row("c_twistfix", "yaw_0")]
    new = [_row("c", "yaw_0")]
    pairs, unmatched = campaign_metrics.resolve_pairs(new, old, {})
    assert pairs == []
    assert unmatched == ["c/yaw_0"]
    # ... and the explicit alias resolves the ambiguity.
    pairs, unmatched = campaign_metrics.resolve_pairs(new, old, {"c": "c_twistfix"})
    assert [key for key, _, _ in pairs] == ["c_twistfix/yaw_0"]
    assert unmatched == []


def test_alias_pointing_at_a_missing_campaign_leaves_the_row_unmatched():
    old = [_row("c_mitc3fix", "yaw_0")]
    new = [_row("c", "yaw_0")]
    pairs, unmatched = campaign_metrics.resolve_pairs(new, old, {"c": "c_does_not_exist"})
    assert pairs == []
    assert unmatched == ["c/yaw_0"]


def test_parse_aliases_rejects_a_malformed_value():
    with pytest.raises(SystemExit):
        campaign_metrics.parse_aliases(["c_mitc3fix"])


def test_compare_prints_delta_relative_to_before(tmp_path, capsys):
    before = tmp_path / "before.csv"
    after = tmp_path / "after.csv"
    _write_csv(before, [_row("c_mitc3fix", "yaw_0", flap_mean=8.0)])
    _write_csv(after, [_row("c", "yaw_0", flap_mean=7.2)])

    rc = campaign_metrics.compare(after, before, {"c": "c_mitc3fix"})
    out = capsys.readouterr().out

    assert rc == 0
    assert "-10.00%" in out  # (7.2 - 8.0) / 8.0
    assert "matched 1 of 1 cases" in out


def test_compare_reports_unmatched_rows(tmp_path, capsys):
    before = tmp_path / "before.csv"
    after = tmp_path / "after.csv"
    _write_csv(before, [_row("some_other_campaign", "yaw_0", flap_mean=8.0)])
    _write_csv(after, [_row("c", "yaw_0", flap_mean=7.2)])

    rc = campaign_metrics.compare(after, before, {})
    out = capsys.readouterr().out

    assert rc == 0
    assert "matched 0 of 1 cases" in out
    assert "--alias NEW=OLD" in out
