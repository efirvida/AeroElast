"""The windIO 2.x rigid blade variant differs from the official upstream yaml only by its transform.

``tests/IEA15MW/reference/upstream/IEA-15-240-RWT_rigid.yaml`` replaces the retired windIO 1.x
``tests/IEA-15-240-RWT_rigid.yaml``: every material stiffness ``E`` and ``G`` is multiplied by 1e4,
``rho``, ``nu`` and every other key are unchanged (gmsh-official-mesher T16 data, S8/S9).
"""

from __future__ import annotations

from pathlib import Path

import pytest
import yaml

TESTS = Path(__file__).resolve().parents[1]
UPSTREAM = TESTS / "IEA15MW" / "reference" / "upstream"
BASE = UPSTREAM / "IEA-15-240-RWT.yaml"
RIGID = UPSTREAM / "IEA-15-240-RWT_rigid.yaml"
CASE = TESTS / "IEA15MW" / "frontiersin_2025_yaw_anim" / "solid_corotational_rigid_yaw_0.yaml"
FACTOR = 1.0e4


def _load(path: Path) -> dict:
    return yaml.safe_load(path.read_text())


def _scaled(value):
    return [v * FACTOR for v in value] if isinstance(value, list) else value * FACTOR


def test_rigid_variant_is_the_upstream_yaml_with_E_and_G_scaled_by_1e4():
    base, rigid = _load(BASE), _load(RIGID)
    expected = {**base, "materials": []}
    for mat in base["materials"]:
        expected["materials"].append({k: (_scaled(v) if k in ("E", "G") else v) for k, v in mat.items()})
    assert rigid.keys() == base.keys()
    for key in base:
        if key != "materials":
            assert rigid[key] == base[key], key
    assert len(rigid["materials"]) == len(base["materials"]) == 11
    for got, want in zip(rigid["materials"], expected["materials"], strict=True):
        assert got.keys() == want.keys()
        for k, v in want.items():
            if k in ("E", "G"):
                assert got[k] == pytest.approx(v, rel=1e-12), (want["name"], k)
            else:
                assert got[k] == v, (want["name"], k)


def test_rigid_variant_names_its_source_and_transform_in_the_header():
    header = "\n".join(RIGID.read_text().splitlines()[:6])
    assert "IEA15MW/reference/upstream/IEA-15-240-RWT.yaml" in header
    assert "1e4" in header


def test_rigid_variant_is_read_by_the_windio2_reader():
    from aeroelast.models.blade.windio import is_windio2, normalize_windio

    data = _load(RIGID)
    assert is_windio2(data)
    assert normalize_windio(data)["materials"][0]["E"] == pytest.approx(3.44e9 * FACTOR)


def test_rigid_yaw0_case_reads_the_variant_through_the_gmsh_mesher():
    gen = _load(CASE)["mesh"]["generator"]
    assert gen["type"] == "GmshBladeMesh"
    assert (CASE.parent / gen["params"]["yaml_file"]).resolve() == RIGID.resolve()
    assert not {"n_samples", "airfoil_spacing", "span_grading", "airfoil_dir"} & gen["params"].keys()
