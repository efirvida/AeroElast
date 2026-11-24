"""``sectional.box_node_mask`` on a REAL windIO v2 blade mesh (l3-outboard-section-stiffness T14, S2, D4; L20 risk).

L20 gave ``box_node_mask`` a v2 branch (SEC<n> laminates, no SPAR/_SW names) but pinned it only on synthetic strips.
This test pins it on the official ``tests/IEA15MW/reference/upstream/IEA-15-240-RWT.yaml`` meshed with the reference mesher (gmsh chord-scaled
N20 AR 1, skin reference ``mid``; ``AGENTS.md`` rule 9), via ``tests/support/v2_blade.build_v2_blade`` (disk-cached under
``$SCRATCH/tmp/v2_blade_cache``; a cold build is slow, a warm one takes ~3 s).

At three named stations (z = 30, 60, 90 m, slab |z - z0| < 0.5 m) the expected box is derived INDEPENDENTLY of the
mask from the mesh's layups (``numad_mesh_data``): spar cap = an element whose layup carries CarbonUD, shear web =
``allWebEls``, TE reinforcement = a glass_uni laminate without CarbonUD, LE/TE panels = every other skin laminate.
Asserted through the public function only: all cap and web nodes are selected; a node touched only by LE/TE-panel or
TE-reinforcement elements is not; and the box fraction of the slab is that of a two-cap-plus-webs box, not the full ring.
"""

import numpy as np
import pytest

pytest.importorskip("gmsh")
from tests.support.v2_blade import build_v2_blade  # noqa: E402

from aeroelast.postprocess.sectional import BOX_MIN_NODES, box_node_mask  # noqa: E402

STATIONS_Z = (30.0, 60.0, 90.0)
HALF_SLAB = 0.5


@pytest.fixture(scope="module")
def blade():
    mesh, _ = build_v2_blade(
        0.25, skin_reference="mid", mesher="gmsh", chord_divisions=20
    )  # element size unused (N20)
    return mesh, box_node_mask(mesh)


def _node_indices(mesh, elements):
    idx = mesh.node_id_to_index
    return {idx[n] for e in elements for n in e.node_ids}


def _classes(mesh):
    """Node-index sets (cap, web, te_reinforcement, panel) from the layups, independent of ``box_node_mask``."""
    cap, te, panel = set(), set(), set()
    for sec in mesh.numad_mesh_data["sections"]:
        nodes = _node_indices(mesh, mesh.element_sets[sec["elementSet"]].elements)
        mats = {ply[0] for ply in sec["layup"]}
        if "allWebEls" in mesh.element_sets and "CarbonUD" in mats:
            cap |= nodes
        elif "glass_uni" in mats:
            te |= nodes
        else:
            panel |= (
                nodes  # skin panels AND web laminates (glass_biax); web nodes are removed below
            )
    web = _node_indices(mesh, mesh.element_sets["allWebEls"].elements)
    return cap, web, te, panel


@pytest.mark.parametrize("z0", STATIONS_Z)
def test_box_mask_selects_caps_and_webs_not_panels_on_real_v2_mesh(blade, z0):
    mesh, mask = blade
    assert not any("SPAR" in n or n.endswith("_SW") for n in mesh.element_sets), "not a v2 mesh"
    cap, web, te, panel = _classes(mesh)
    z = mesh.coords_array[:, 2]
    slab = set(np.nonzero(np.abs(z - z0) < HALF_SLAB)[0].tolist())
    sel = {i for i in slab if mask[i]}

    box = (cap | web) & slab
    assert len(cap & slab) > 0 and len(web & slab) > 0 and len(sel) >= BOX_MIN_NODES
    assert box <= sel, "a cap or web node was left out of the box"

    only_outside = slab - cap - web  # touched by no cap and no web element
    assert only_outside, "station has no LE/TE panel nodes"
    assert not (only_outside & sel), (
        "a node of only LE/TE-panel / TE-reinforcement elements was selected"
    )
    te_only = (te - cap - web) & slab
    assert te_only and not (te_only & sel)

    # The box is a strict subset of the ring: caps + webs are 30-60 % of the slab nodes, not all of them.
    frac = len(sel) / len(slab)
    assert 0.30 < frac < 0.60, f"box fraction {frac:.2f} at z={z0}"


def test_box_mask_has_no_nodes_outside_box_stations_overall(blade):
    mesh, mask = blade
    cap, web, _, _ = _classes(mesh)
    assert set(np.nonzero(mask)[0].tolist()) == cap | web
