"""Blade/rotor meshing without shear webs (CFD outer-mold-line surface)."""

import os

import pytest

from aeroelast.core.mesh.generators import BladeMesh, RotorMesh

_HERE = os.path.dirname(os.path.abspath(__file__))
_YAML = os.path.join(_HERE, "IEA-15-240-RWT.yaml")

# Coarse enough to keep the test fast, fine enough that the shear web elements
# are a meaningful fraction of the total.
_ELEMENT_SIZE = 1.0


def _is_web_name(name: str) -> bool:
    lowered = name.lower()
    return "web" in lowered or lowered.endswith("_sw")


@pytest.fixture(scope="module")
def with_webs():
    return BladeMesh(_YAML, element_size=_ELEMENT_SIZE, include_webs=True).generate(verbose=False)


@pytest.fixture(scope="module")
def no_webs():
    return BladeMesh(_YAML, element_size=_ELEMENT_SIZE, include_webs=False).generate(verbose=False)


def test_no_webs_has_no_web_sets(no_webs):
    for name in list(no_webs.element_sets) + list(no_webs.node_sets):
        assert not _is_web_name(name), f"web set leaked into the surface mesh: {name}"
    assert "allShearWebEls" not in no_webs.element_sets
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
    web_elements = len(with_webs.element_sets["allShearWebEls"].elements)
    assert web_elements > 0
    assert with_webs.elements_count - no_webs.elements_count == web_elements


def test_no_webs_is_smaller(with_webs, no_webs):
    assert no_webs.elements_count < with_webs.elements_count
    assert no_webs.node_count < with_webs.node_count


def test_rotor_forwards_include_webs():
    rotor = RotorMesh(
        _YAML,
        n_blades=2,
        hub_radius=3.97,
        element_size=_ELEMENT_SIZE * 2,
        include_webs=False,
    )
    mesh = rotor.generate(verbose=False)
    assert "allShearWebEls" not in mesh.element_sets
    assert "allOuterShellEls_blade_1" in mesh.element_sets
    assert "allOuterShellEls_blade_2" in mesh.element_sets
