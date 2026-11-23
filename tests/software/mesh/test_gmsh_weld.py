"""Surfaces that share an attachment line share nodes (odd/tasks/gmsh-official-mesher.md T19, S10).

web0 and web1 of the IEA-15 blade share an attachment line at the root (z 9.525-9.578 m); each web was meshed on its own and
the mesh carried 16 node pairs at distance 0 with distinct ids.  ``weld_coincident_nodes`` merges them for any pair of surfaces.
"""

import numpy as np
import pytest

pytest.importorskip("gmsh", reason="gmsh not available")

from aeroelast.models.blade.numad.mesh_gen import gmsh_skin as gs  # noqa: E402
from tests.support import v2_blade  # noqa: E402


def _two_quads_sharing_a_line():
    # quad A: 0 1 2 3, quad B: 4 5 6 7 with 4 == 0-position... B sits on top of A's edge 1-2 (positions of 4,5 = 1,2)
    nodes = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [1, 0, 0], [1, 1, 0], [2, 1, 0], [2, 0, 0], [5, 5, 5]], float)
    elements = np.array([[0, 1, 2, 3], [4, 7, 6, 5], [1, 7, 8, -1]])
    return nodes, elements


def test_coincident_ids_are_merged_and_connectivity_remapped():
    nodes, elements = _two_quads_sharing_a_line()
    attach = {("w", "start"): (0, 4, 5, 8)}
    n2, e2, a2, removed = gs.weld_coincident_nodes(nodes, elements, attach)
    assert removed == 2 and len(n2) == 7
    assert not gs.conformity_report(n2, e2)["coincident_distinct_id_pairs"]
    assert gs.conformity_report(n2, e2)["unused_nodes"] == 0
    assert np.array_equal(n2[e2[0]], nodes[elements[0]]) and np.array_equal(n2[e2[1]], nodes[elements[1]])
    assert e2[2, 3] == -1 and np.array_equal(n2[e2[2, :3]], nodes[elements[2, :3]])  # padding kept, positions unchanged
    assert tuple(n2[list(a2[("w", "start")])].ravel()) == tuple(nodes[[0, 4, 5, 8]].ravel())
    assert a2[("w", "start")][1] == e2[0][1]  # the attachment id is the surviving id


def test_no_coincident_pair_leaves_the_mesh_untouched():
    nodes, elements = _two_quads_sharing_a_line()
    n2, e2, a2, removed = gs.weld_coincident_nodes(nodes[[0, 1, 2, 3]], elements[:1], {"k": (0, 1)})
    assert removed == 0 and np.array_equal(e2, elements[:1]) and a2 == {"k": (0, 1)} and len(n2) == 4


@pytest.mark.slow
def test_iea15_gmsh_mesh_has_no_two_node_ids_at_one_position():
    from aeroelast.models.blade.numad.objects.blade import Blade

    blade = Blade()
    blade.read_yaml(str(v2_blade.DEFAULT_YAML))
    raw, info = gs.generate_mesh_gmsh(blade.windio2_topology, 0.25, with_materials=True, return_info=True)
    xyz, el = np.asarray(raw.nodes, float), np.asarray(raw.elements)
    assert not gs.conformity_report(xyz, el)["coincident_distinct_id_pairs"]
    assert info["welded_nodes"] == 16  # web0/web1 root line, z 9.525-9.578 m
    assert len(xyz) == 24996 and len(el) == 26015
    assert len(np.unique(el[el >= 0])) == len(xyz)  # every node used
