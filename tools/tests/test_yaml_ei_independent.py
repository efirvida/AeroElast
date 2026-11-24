"""yaml_ei_independent on the official windIO v2 blade (T18): positive principal EI, edge >> flap at mid-span."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import yaml_ei_independent as ei  # noqa: E402
import yaml_layup_independent as ev  # noqa: E402


def test_principal_ei_v2_midspan():
    flap, edge, ea = ei.principal_ei(ev.load(), 0.525)
    assert 0 < flap < edge and ea > 0
