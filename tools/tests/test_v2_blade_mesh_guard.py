"""Pin the reference-mesh rule (AGENTS.md rule 9) in tools/v2_blade_mesh.py.

element_size >= 0.3 m is not converged for the IEA-15 blade and must be rejected;
without a size the helper builds the chord-scaled reference mesh N20 AR 1.
"""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import v2_blade_mesh  # noqa: E402


@pytest.mark.parametrize("es", [0.3, 0.5, 1.0])
def test_element_size_at_or_above_0p3_is_rejected(es):
    with pytest.raises(ValueError, match="0.3"):
        v2_blade_mesh.blade_generator(element_size=es)


def test_reference_mesh_constants():
    assert (v2_blade_mesh.REFERENCE_CHORD_DIVISIONS, v2_blade_mesh.REFERENCE_SPAN_ASPECT) == (20, 1.0)
