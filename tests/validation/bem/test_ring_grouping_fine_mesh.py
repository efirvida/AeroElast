"""Physical-ring grouping must follow the mesh, not a tolerance tuned on the 1.0 m mesh.

The IEA-15 MW blade is prebent: one physical ring's nodes do not lie in one plane.  At
``element_size`` 0.25 the ring at z = 107.44 m has 2 nodes at 107.4369 and 26 at 107.449
(skew gap 0.0121 m), while the span-extent gap tolerance (``RING_GAP_FRACTION``, 1e-4 of
117 m = 0.0117 m) was tuned on the 1.0 m mesh, whose skew is 0.0072 m.  Single linkage on
the span coordinate therefore split those rings (26 + 2 nodes) and the elastic twist,
fitted on the 26-node arc, differed from the whole ring's by 0.03-0.04 deg
(``test_elastic_twist_matches_the_nodal_section_rotation``, strips 44-47).  The same mesh
is graded at the tip, ring spacing 0.0746-0.109 m, which is below the 0.119 m (5 % of the
strip width) the twist ring binder used before: it fused two rings (strips 46/49: 0.1808 deg).

Measured on that mesh (SLURM 11622584): every element edge inside a ring has
``|dz| / length <= 0.088``; every edge between rings has ``|dz| / length >= 0.559``.  The
element topology separates rings at any resolution; these tests drive that on synthetic
meshes (no solve) and on the strip-ring selection of the twist.
"""

from __future__ import annotations

import numpy as np
import pytest

from aeroelast.core.mesh.entities import ElementType, MeshElement, Node
from aeroelast.core.mesh.model import MeshModel
from aeroelast.solvers.bem.force_projection import RING_GAP_FRACTION, physical_ring_labels
from aeroelast.solvers.bem.section_contour import rings_in_band

N_PER_RING = 12
SKEW = 0.0121  # span gap inside a ring (measured 0.0121-0.0126 m at element_size 0.25)
SPAN = np.array([0.0, 0.0, 1.0])


def _stack(spacing: float, n_rings: int = 6):
    """Rings of 12 nodes in the xy plane; the last two nodes of each sit ``SKEW`` higher.

    Quads join consecutive rings.  A root ring at span 0 and a tip ring at 117 m make the
    mesh span extent that of the real blade, so the span-extent tolerance is 0.0117 m.
    """
    mesh = MeshModel()
    angles = np.linspace(0.0, 2.0 * np.pi, N_PER_RING, endpoint=False)
    z0 = np.concatenate([[0.0], 100.0 + spacing * np.arange(n_rings), [117.0]])
    rings = []
    for z in z0:
        ring = []
        for j, a in enumerate(angles):
            node = Node(
                [0.4 * np.cos(a), 0.1 * np.sin(a), z + (SKEW if j >= N_PER_RING - 2 else 0.0)]
            )
            mesh.add_node(node)
            ring.append(node)
        rings.append(ring)
    for lower, upper in zip(rings[:-1], rings[1:], strict=True):
        for j in range(N_PER_RING):
            k = (j + 1) % N_PER_RING
            mesh.add_element(
                MeshElement([lower[j], lower[k], upper[k], upper[j]], ElementType.quad)
            )
    expected = np.repeat(np.arange(len(z0)), N_PER_RING)
    return mesh, expected


def _same_partition(labels: np.ndarray, expected: np.ndarray) -> bool:
    return len(np.unique(labels)) == len(np.unique(expected)) and all(
        len(set(labels[expected == r])) == 1 for r in np.unique(expected)
    )


@pytest.mark.parametrize("spacing", [0.30, 0.25, 0.109, 0.0746])
def test_rings_are_the_physical_stations_whatever_the_spacing(spacing):
    """One label per station, the skewed pair kept with its ring, at every spacing.

    The 0.0121 m skew gap exceeds the 0.0117 m span-extent tolerance, so span gap alone
    splits every ring (precondition asserted); 0.0746-0.109 m spacings are below 5 % of
    a 2.388 m strip (the former twist binder) yet are separate rings.
    """
    mesh, expected = _stack(spacing)
    span = mesh.coords_array @ SPAN
    gap = RING_GAP_FRACTION * float(np.ptp(span))
    assert len(rings_in_band(span, gap)) > len(np.unique(expected))  # the defect's premise

    labels = physical_ring_labels(mesh, span, gap)

    assert _same_partition(labels, expected)


def test_rings_without_elements_fall_back_to_the_span_gap():
    """Nodes with no element edge still group by span proximity (previous behaviour)."""
    mesh = MeshModel()
    for z in (0.0, 0.001, 0.002, 1.0, 1.001, 117.0):
        mesh.add_node(Node([z * 0.1, 0.0, z]))
    span = mesh.coords_array @ SPAN
    labels = physical_ring_labels(mesh, span, RING_GAP_FRACTION * float(np.ptp(span)))
    assert labels[0] == labels[1] == labels[2]
    assert labels[3] == labels[4]
    assert len({labels[0], labels[3], labels[5]}) == 3
