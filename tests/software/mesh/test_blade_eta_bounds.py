"""``BladeMesh.generate(eta_bounds=...)`` meshes only a spanwise sub-range of the blade."""

import pytest

from aeroelast.core.mesh.generators import BladeMesh
from tests.support.paths import V2_BLADE_YAML

_YAML = str(V2_BLADE_YAML)
_CHORD_DIVISIONS = 10  # reference chord-scaled sizing, cheap (AGENTS rule 9: nothing >= 0.3 m)
_BOUNDS = (0.4, 0.6)


@pytest.fixture(scope="module")
def full():
    return BladeMesh(_YAML, chord_divisions=_CHORD_DIVISIONS, include_webs=False).generate()


@pytest.fixture(scope="module")
def patch():
    return BladeMesh(_YAML, chord_divisions=_CHORD_DIVISIONS, include_webs=False).generate(eta_bounds=_BOUNDS)


def test_sub_range_spans_only_those_stations(full, patch):
    z_full = full.coords_array[:, 2]
    length = z_full.max() - z_full.min()
    z = patch.coords_array[:, 2]
    lo, hi = z_full.min() + _BOUNDS[0] * length, z_full.min() + _BOUNDS[1] * length
    # the straight-span coordinate; prebend/sweep do not move z by more than the section offset
    assert z.min() == pytest.approx(lo, abs=0.5)
    assert z.max() == pytest.approx(hi, abs=0.5)
    assert 0 < patch.node_count < full.node_count
    assert patch.elements_count < full.elements_count


def test_sub_range_keeps_the_contract_sets(patch):
    assert {"RootNodes", "TipNodes", "allOuterShellNods"} <= set(patch.node_sets)
    assert len(patch.element_sets["allOuterShellEls"].elements) == patch.elements_count


@pytest.mark.parametrize("bounds", [(0.6, 0.4), (-0.1, 0.5), (0.2, 1.2), (0.5, 0.5)])
def test_invalid_eta_bounds_raise(bounds):
    with pytest.raises(ValueError, match="eta_bounds"):
        BladeMesh(_YAML, chord_divisions=_CHORD_DIVISIONS).generate(eta_bounds=bounds)
