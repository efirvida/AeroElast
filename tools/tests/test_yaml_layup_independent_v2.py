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


def test_no_inverted_interval_with_positive_thickness():
    import numpy as np

    st = ev.load()["components"]["blade"]["structure"]
    tab = ev.anchor_table(st)
    for lay in st["layers"]:
        if "web" in lay:
            continue
        for eta in np.linspace(0, 1, 201):
            s, e = (ev.arc(tab, lay[k], eta) for k in ("start_nd_arc", "end_nd_arc"))
            if e < s:
                assert ev.interp(lay["thickness"], eta) == 0.0, (lay["name"], eta)


def test_inverted_interval_with_thickness_raises():
    const = lambda v: dict(grid=[0.0, 1.0], values=[v, v])  # noqa: E731
    tab = {"a": dict(start_nd_arc=const(0.6), end_nd_arc=const(0.4))}
    layer = dict(name="x", thickness=const(0.01),
                 start_nd_arc=dict(anchor=dict(name="a", handle="start_nd_arc")),
                 end_nd_arc=dict(anchor=dict(name="a", handle="end_nd_arc")))
    with pytest.raises(ValueError, match="Invalid nd_arc interval"):
        ev.layer_interval(tab, layer, 0.5)


def test_te_declared_width_is_reported_not_used():
    """L-S61 finding: the v2 TE_reinforcement nd_arc is ~6.8 % wider than its declared 0.75 m at eta 0.30."""
    r = ev.section(ev.load(), 0.30)
    te = next(lay for lay in r["layers"] if lay["name"] == "TE_reinforcement_SS")
    assert te["declared_width_m"] == pytest.approx(0.75)
    assert te["width_residual_m"] / 0.75 == pytest.approx(0.0676, abs=2e-3)
