"""blade_mass_report on the official windIO v2 blade (T18, S8/S9): reference is the yaml's own declared mass."""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import blade_mass_report as bmr  # noqa: E402
import yaml_layup_independent as ev  # noqa: E402


def test_blade_mass_is_within_2pct_of_declared_inertia_matrix():
    """Measured 2026-11-23: layup 67740.1 kg vs declared 67015.1 kg (+1.08 %); 2 % is a sanity bound, not an anchor."""
    layup, declared = ev.blade_mass(ev.load())
    assert declared == pytest.approx(67015.1, rel=1e-3)
    assert layup / declared - 1 == pytest.approx(0.0108, abs=2e-3)


def test_report_lines_and_v1_rejected(tmp_path):
    lines = bmr.report(ev.load(), stations=2)
    assert any(line.startswith("BLADE_MASS") for line in lines)
    v1 = tmp_path / "v1.yaml"
    v1.write_text("windIO_version: '1.0'\n")
    with pytest.raises(ValueError, match="windIO 2.x required"):
        ev.load(v1)
