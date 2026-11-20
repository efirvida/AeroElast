"""Pin the windIO v2 port of the independent layup evaluator (T18, S9).

The square skin+web+anchored-band section is the only exact expectation (perimeter 8 m,
web height 2 m, 140 kg/m, TE-style width residual 0); a windIO 1.x file must be rejected.
"""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import yaml_layup_independent as ev  # noqa: E402


def test_self_check_square_section_passes(capsys):
    ev.self_check()
    assert "SELF-CHECK PASS" in capsys.readouterr().out


def test_windio_1x_file_is_rejected(tmp_path):
    v1 = tmp_path / "v1.yaml"
    v1.write_text("windIO_version: '1.0'\ncomponents: {}\n")
    with pytest.raises(ValueError, match="windIO 2.x required"):
        ev.load(v1)


def test_official_v2_blade_loads():
    data = ev.load()
    assert str(data["windIO_version"]).startswith("2")
