"""Blade/rotor meshing without shear webs (CFD outer-mold-line surface)."""

import os

import numpy as np
import pytest

from aeroelast.core.mesh.generators import BladeMesh, RotorMesh

_HERE = os.path.dirname(os.path.abspath(__file__))
from tests.support.paths import DATA_DIR  # noqa: E402
_YAML = str(DATA_DIR / "IEA15MW" / "reference" / "upstream" / "IEA-15-240-RWT.yaml")  # windIO v2 (T6c)

# T6c: BladeMesh meshes with gmsh; chord_divisions 10 keeps the with-webs reference cheap (AGENTS rule 9).
_CHORD_DIVISIONS = 10


def _is_web_name(name: str) -> bool:
    lowered = name.lower()
    return "web" in lowered or lowered.endswith("_sw")


@pytest.fixture(scope="module")
def with_webs():
    return BladeMesh(_YAML, chord_divisions=_CHORD_DIVISIONS, include_webs=True).generate(verbose=False)


@pytest.fixture(scope="module")
def no_webs():
    return BladeMesh(_YAML, chord_divisions=_CHORD_DIVISIONS, include_webs=False).generate(verbose=False)


def test_no_webs_has_no_web_sets(no_webs):
    for name in list(no_webs.element_sets) + list(no_webs.node_sets):
        assert not _is_web_name(name), f"web set leaked into the surface mesh: {name}"
    assert "allWebEls" not in no_webs.element_sets
    assert "allShearWebNods" not in no_webs.node_sets


def test_no_webs_keeps_the_outer_shell(no_webs):
    assert "allOuterShellEls" in no_webs.element_sets
    assert len(no_webs.element_sets["allOuterShellEls"].elements) > 0
    assert len(no_webs.elements) > 0
    assert len(no_webs.nodes) > 0


def test_outer_shell_is_unchanged(with_webs, no_webs):
    assert len(no_webs.element_sets["allOuterShellEls"].elements) == len(
        with_webs.element_sets["allOuterShellEls"].elements
    )


def test_web_elements_are_exactly_the_difference(with_webs, no_webs):
    web_elements = len(with_webs.element_sets["allWebEls"].elements)
    assert web_elements > 0
    assert with_webs.elements_count - no_webs.elements_count == web_elements


def test_no_webs_is_smaller(with_webs, no_webs):
    assert no_webs.elements_count < with_webs.elements_count
    assert no_webs.node_count < with_webs.node_count


def test_no_webs_has_no_orphan_nodes_and_the_skin_is_the_webbed_one(with_webs, no_webs):
    used = {id(n) for e in no_webs.elements for n in e.nodes}
    assert len(used) == no_webs.node_count
    skin = with_webs.element_sets["allOuterShellEls"].elements
    skin_nodes = {id(n) for e in skin for n in e.nodes}
    assert len(skin_nodes) == no_webs.node_count
    ref = np.array(sorted(map(tuple, np.round([n.coords for e in skin for n in e.nodes], 9))))
    got = np.array(sorted(map(tuple, np.round([n.coords for e in no_webs.elements for n in e.nodes], 9))))
    np.testing.assert_array_equal(got, ref)
    assert {"RootNodes", "TipNodes", "allOuterShellNods"} <= set(no_webs.node_sets)


def test_no_webs_sections_all_have_elements(no_webs):
    sections = [s["elementSet"] for s in no_webs.numad_mesh_data["sections"]]
    assert sections and all(name in no_webs.element_sets for name in sections)


def test_rotor_forwards_include_webs():
    rotor = RotorMesh(_YAML, n_blades=2, hub_radius=3.97, chord_divisions=_CHORD_DIVISIONS, include_webs=False)
    mesh = rotor.generate(verbose=False)
    assert not any(_is_web_name(n) for n in mesh.element_sets)
    assert mesh.elements_count > 0
