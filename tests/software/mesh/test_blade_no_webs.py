"""Blade/rotor meshing without shear webs (CFD outer-mold-line surface)."""

import os

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


# T7 (odd/tasks/gmsh-official-mesher.md, S7e): the surface-only gmsh blade is not built yet; BladeMesh(include_webs=False) raises
# NotImplementedError.  The tests below that need it are the T7 target contract: strict xfail, so finishing T7 flips them to XPASS.
_T7 = pytest.mark.xfail(raises=NotImplementedError, strict=True, reason="T7: include_webs=False on gmsh")


@pytest.fixture(scope="module")
def no_webs():
    return BladeMesh(_YAML, chord_divisions=_CHORD_DIVISIONS, include_webs=False).generate(verbose=False)


@_T7
def test_no_webs_has_no_web_sets(no_webs):
    for name in list(no_webs.element_sets) + list(no_webs.node_sets):
        assert not _is_web_name(name), f"web set leaked into the surface mesh: {name}"
    assert "allWebEls" not in no_webs.element_sets
    assert "allShearWebNods" not in no_webs.node_sets


@_T7
def test_no_webs_keeps_the_outer_shell(no_webs):
    assert "allOuterShellEls" in no_webs.element_sets
    assert len(no_webs.element_sets["allOuterShellEls"].elements) > 0
    assert len(no_webs.elements) > 0
    assert len(no_webs.nodes) > 0


@_T7
def test_outer_shell_is_unchanged(with_webs, no_webs):
    assert len(no_webs.element_sets["allOuterShellEls"].elements) == len(
        with_webs.element_sets["allOuterShellEls"].elements
    )


@_T7
def test_web_elements_are_exactly_the_difference(with_webs, no_webs):
    web_elements = len(with_webs.element_sets["allWebEls"].elements)
    assert web_elements > 0
    assert with_webs.elements_count - no_webs.elements_count == web_elements


@_T7
def test_no_webs_is_smaller(with_webs, no_webs):
    assert no_webs.elements_count < with_webs.elements_count
    assert no_webs.node_count < with_webs.node_count


def test_blade_include_webs_false_surfaces_the_t7_error():
    with pytest.raises(NotImplementedError, match="T7"):
        BladeMesh(_YAML, chord_divisions=_CHORD_DIVISIONS, include_webs=False)


def test_rotor_forwards_include_webs():
    # a forwarded include_webs=False must reach BladeMesh and surface its T7 error (it is not silently dropped)
    rotor = RotorMesh(
        _YAML,
        n_blades=2,
        hub_radius=3.97,
        chord_divisions=_CHORD_DIVISIONS,
        include_webs=False,
    )
    with pytest.raises(NotImplementedError, match="T7"):
        rotor.generate(verbose=False)
