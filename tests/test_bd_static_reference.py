"""Deterministic BeamDyn static input/output checks; never launches BeamDyn."""
import importlib.util
import sys
from pathlib import Path

import pytest

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'tools'))
SPEC = importlib.util.spec_from_file_location(
    "bd_static_reference", Path(__file__).resolve().parents[1] / "tools/bd_static_reference.py"
)
bd = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(bd)


@pytest.mark.parametrize("case,index,value", [("LC1", 0, -9.81), ("LC2", 1, 1e6), ("LC4", 2, 5e3)])
def test_loads(case, index, value):
    loads = bd.case_loads(case, 0.01)
    assert loads[index] == pytest.approx((value * 0.01, 0, 0))
    assert all(v == (0, 0, 0) for i, v in enumerate(loads) if i != index)


@pytest.mark.parametrize("scale", [0, -1, float("nan"), float("inf")])
def test_invalid_scale(scale):
    with pytest.raises(ValueError):
        bd.case_loads("LC2", scale)


def test_parser_reordered_channels_and_fortran_exponent():
    text = "header\nTime TipTDyr TipTDzr TipTDxr\n(s) (m) (m) (m)\n0 0 0 0\n1 -2D-3 3e-4 7.5e-2\n"
    assert bd.parse_output(text) == {"TipTDxr": 0.075, "TipTDyr": -0.002, "TipTDzr": 0.0003}


@pytest.mark.parametrize("text", ["", "Time TipTDxr\n(s) (m)\n1 2", "Time TipTDxr TipTDyr TipTDzr\n(s) (m) (m) (m)\n1 NaN 0 0"])
def test_parser_rejects_incomplete_or_nonfinite(text):
    with pytest.raises(ValueError):
        bd.parse_output(text)


@pytest.mark.parametrize("units,row", [
    ("(s) (N) (m) (m)", "1 1 0 0"),
    ("(s) (m) (m) (m)", "0 0 0 0"),
])
def test_parser_rejects_wrong_units_or_unsolved_initial_state(units, row):
    with pytest.raises(ValueError):
        bd.parse_output(f"Time TipTDxr TipTDyr TipTDzr\n{units}\n{row}\n")


def test_build_preserves_geometry_and_matrices(tmp_path, monkeypatch):
    monkeypatch.setattr(bd, "INPUTS", tmp_path)
    name = bd.build("LC4", 15, 0.001)
    beam = (tmp_path / name / 'beam.dat').read_text()
    source = bd.GEOMETRY.read_text()
    assert beam.split('(deg)\n')[1].split('---')[0] == source.split('(deg)\n')[1].split('---')[0]
    assert '1   member_total' in beam
    assert '15   order_elem' in beam
    assert 'PITCH ACTUATOR' not in beam
    prop = (tmp_path / name / 'blade.dat').read_text()
    assert prop.count('MODAL DAMPING') == 1
    assert prop.split('DISTRIBUTED PROPERTIES')[1] == bd.PROPERTIES.read_text().split('DISTRIBUTED PROPERTIES')[1]
    with pytest.raises(FileExistsError):
        bd.build("LC4", 15, 0.001)


def test_driver_static_clamped_identity():
    text = bd.driver_text("LC2", 0.01, Path("beam.dat"))
    assert "False DynamicSolve" in text
    assert "1 0 0\n0 1 0\n0 0 1" in text
    assert "10000 TipFx" in text
    assert "0 omega_x\n0 omega_y\n0 omega_z" in text
