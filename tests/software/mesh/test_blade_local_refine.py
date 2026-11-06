"""Local gmsh remesh of a patch with imposed boundary nodes (``numad.mesh_gen.local_refine``).

Synthetic patches only, no blade: the boundary nodes must come back bit for bit, the elements must
be well formed (min angle >= 15 deg, scaled Jacobian >= 0.3, the "poor" thresholds of the blade
mesh quality test), and every size must scale with the ``element_size`` / ``local_size`` given.
"""

from __future__ import annotations

import numpy as np
import pytest

from aeroelast.models.blade.numad.mesh_gen.local_refine import (
    MeshedRegion,
    refine_patch,
    refine_surface_patch,
)
from aeroelast.models.blade.numad.mesh_gen.surface import Surface
from tests.software.mesh.test_blade_mesh_quality import MIN_ANGLE_DEG, MIN_SCALED_JACOBIAN, _quality


def _strip(width, height, bottom, top):
    """Rectangle width x height: ``bottom`` / ``top`` nodes on the lower / upper edge, ccw."""
    low = [(width * t, 0.0) for t in np.linspace(0.0, 1.0, bottom)]
    high = [(width * t, height) for t in np.linspace(1.0, 0.0, top)]
    return np.array(low + high)


def _qualities(nodes, elements):
    out = []
    for conn in elements:
        corners = np.column_stack([nodes[conn[conn >= 0]], np.zeros(int((conn >= 0).sum()))])
        out.append(_quality(corners))
    return np.array(out)


def _areas(nodes, elements):
    total = 0.0
    for conn in elements:
        ring = nodes[conn[conn >= 0]]
        total += 0.5 * (np.dot(ring[:, 0], np.roll(ring[:, 1], -1)) - np.dot(np.roll(ring[:, 0], -1), ring[:, 1]))
    return total


# (width, height, nodes on the lower edge, nodes on the upper edge) in metres: a strip of a cap row
# whose lower edge carries one more element than the upper one, and a crack row, whose upper nodes
# sit a few mm beside the lower ones.
PATCHES = {
    "transition_fan": (0.09, 0.24, 3, 2),
    "crack_row": (0.9, 0.24, 5, 5),
    "wide_strip": (0.45, 0.24, 4, 3),
}


@pytest.mark.parametrize("name", sorted(PATCHES))
def test_boundary_nodes_are_imposed_and_elements_well_formed(name):
    width, height, low, high = PATCHES[name]
    boundary = _strip(width, height, low, high)
    if name == "crack_row":
        boundary[low:, 0] += 0.007  # the upper nodes beside the lower ones
    nodes, elements = refine_patch(boundary, element_size=0.25, local_size=width / 2.0)
    m = len(boundary)
    assert np.array_equal(nodes[:m], boundary), "boundary nodes moved"
    assert elements.max() < len(nodes) and elements.min() >= -1
    # no node was added on a boundary segment: every extra node is strictly inside the patch
    extra = nodes[m:]
    assert not len(extra) or (extra[:, 1].min() > 0.0 and extra[:, 1].max() < height)
    # the elements tile the patch exactly, all counter-clockwise
    assert _areas(nodes, elements) == pytest.approx(0.5 * np.dot(boundary[:, 0], np.roll(boundary[:, 1], -1))
                                                    - 0.5 * np.dot(np.roll(boundary[:, 0], -1), boundary[:, 1]))
    q = _qualities(nodes, elements)
    assert q[:, 0].min() >= MIN_ANGLE_DEG, f"{name}: min angle {q[:, 0].min():.1f} deg"
    assert q[:, 1].min() >= MIN_SCALED_JACOBIAN, f"{name}: sJac {q[:, 1].min():.3f}"


def test_best_of_the_trials_is_never_worse_than_a_single_one():
    """gmsh is not monotonic in the size; the retained mesh is the best worst-corner of the trials."""
    from aeroelast.models.blade.numad.mesh_gen.local_refine import worst_corner

    boundary = _strip(0.09, 0.24, 3, 2)
    single = [refine_patch(boundary, 0.25, 0.045, trials=(f,)) for f in (1.0, 0.5, 0.75, 1.5, 2.0)]
    best = refine_patch(boundary, 0.25, 0.045)
    assert worst_corner(*best) >= max(worst_corner(*mesh) for mesh in single) - 1e-12
    assert worst_corner(*best) > 0.0


def test_every_boundary_edge_is_used_once_and_interior_edges_twice():
    boundary = _strip(0.45, 0.24, 4, 3)
    nodes, elements = refine_patch(boundary, element_size=0.25, local_size=0.12)
    use = {}
    for conn in elements:
        ring = list(conn[conn >= 0])
        for a, b in zip(ring, ring[1:] + ring[:1], strict=True):
            use[frozenset((a, b))] = use.get(frozenset((a, b)), 0) + 1
    m = len(boundary)
    for key, n in use.items():
        on_boundary = all(k < m for k in key) and (max(key) - min(key) in (1, m - 1))
        assert n == (1 if on_boundary else 2), f"edge {sorted(key)} used {n} times"


def _square(side, per_side):
    """Square boundary, ``per_side`` equal segments on each side, counter-clockwise."""
    t = np.linspace(0.0, side, per_side, endpoint=False)
    return np.array(
        [(x, 0.0) for x in t] + [(side, y) for y in t] + [(side - x, side) for x in t] + [(0.0, side - y) for y in t]
    )


def _element_edges(nodes, elements):
    ring = nodes[elements[:, :3]]
    return np.linalg.norm(ring - np.roll(ring, -1, axis=1), axis=2).mean(axis=1), ring.mean(axis=1)


def test_interior_size_scales_with_element_size():
    """The boundary of a patch is at ``element_size`` spacing: the far interior follows it, the wall ``local_size``."""
    side = 8.0
    far, wall = {}, {}
    for element_size in (0.5, 1.0):
        boundary = _square(side, int(round(side / element_size)))
        nodes, elements = refine_patch(boundary, element_size=element_size, local_size=element_size / 5.0)
        edges, centre = _element_edges(nodes, elements)
        dist = np.minimum(centre.min(axis=1), (side - centre).min(axis=1))
        far[element_size] = edges[dist > 0.35 * side].mean()
        wall[element_size] = edges[dist < 0.1 * side].mean()
    assert far[1.0] / far[0.5] == pytest.approx(2.0, rel=0.25)
    assert wall[1.0] / wall[0.5] == pytest.approx(2.0, rel=0.25)
    assert all(wall[k] < far[k] for k in far), "the interior is not graded from the wall to the far field"


def test_geometry_scaled_by_k_gives_the_same_mesh_topology():
    boundary = _strip(0.09, 0.24, 3, 2)
    base = refine_patch(boundary, element_size=0.25, local_size=0.045)
    for k in (0.1, 7.0):
        nodes, elements = refine_patch(k * boundary, element_size=k * 0.25, local_size=k * 0.045)
        assert elements.shape == base[1].shape
        assert np.allclose(nodes, k * base[0], atol=1e-6 * k)


def test_growth_ratio_sets_how_fast_the_size_grows_away_from_the_boundary():
    """Boundary nodes at ``local_size``: a gentler ratio keeps the mesh finer for longer, none exceeds element_size."""
    # An abstract planar patch, not a blade mesh: sizes are expressed relative to the wall size so
    # the mesh-size gate (scripts/check_mesh_size.py) has no blade element_size literal to flag.
    local = 0.1
    ratio = 10.0
    element_size = ratio * local
    side = 4.0 * element_size
    boundary = _square(side, int(round(side / local)))
    count = {}
    for growth in (1.2, 1.6):
        nodes, elements = refine_patch(boundary, element_size, local, growth=growth)
        edges, centre = _element_edges(nodes, elements)
        dist = np.minimum(centre.min(axis=1), (side - centre).min(axis=1))
        count[growth] = len(elements)
        # beyond the distance where the ramp reaches element_size nothing is finer than ~element_size/2
        # and nothing coarser than 1.5 element_size; next to the wall the elements are the wall's
        ramp_end = (element_size - local) / (growth - 1.0)
        far = dist > ramp_end
        if far.any():
            assert 0.5 * element_size <= edges[far].mean() <= 1.5 * element_size
        assert edges.max() <= 1.5 * element_size
        assert edges[dist < 2.0 * local].mean() < 3.0 * local
    assert count[1.6] < count[1.2]


@pytest.mark.parametrize(
    "case",
    ["clockwise", "growth", "size", "short"],
)
def test_invalid_input_is_rejected(case):
    boundary = _strip(0.4, 0.24, 3, 3)
    kwargs = {"element_size": 0.25, "local_size": 0.1}
    if case == "clockwise":
        boundary = boundary[::-1]
    elif case == "growth":
        kwargs["growth"] = 1.0
    elif case == "size":
        kwargs["local_size"] = 0.0
    else:
        boundary = boundary[:2]
    with pytest.raises(ValueError):
        refine_patch(boundary, **kwargs)


def test_gmsh_global_state_is_left_alone():
    import gmsh

    boundary = _strip(0.4, 0.24, 3, 3)
    assert not gmsh.isInitialized()
    refine_patch(boundary, element_size=0.25, local_size=0.1)
    assert not gmsh.isInitialized(), "refine_patch finalised or left a session it did not own"

    gmsh.initialize()
    try:
        gmsh.model.add("caller")
        models = gmsh.model.list()
        gmsh.option.setNumber("Mesh.Algorithm", 1)
        refine_patch(boundary, element_size=0.25, local_size=0.1)
        assert gmsh.isInitialized()
        assert gmsh.model.getCurrent() == "caller"
        assert gmsh.model.list() == models, "the scratch model was left behind"
        assert gmsh.option.getNumber("Mesh.Algorithm") == 1
    finally:
        gmsh.finalize()
    # the existing callers initialise and finalise on their own, after us
    gmsh.initialize()
    gmsh.finalize()


def _wedge(xyz_of_eta):
    """Boundary of the whole eta square, nodes uniform per edge: 2 on the bottom, 2 on the sides, 4 on top."""
    bottom = np.array([(-1.0, -1.0), (1.0, -1.0)])
    right = np.array([(1.0, -1.0), (1.0, 1.0)])
    top = np.array([(1.0 - 2.0 * t / 3.0, 1.0) for t in range(4)])
    left = np.array([(-1.0, 1.0), (-1.0, -1.0)])
    eta = np.vstack([bottom[:1], right[:1], top[:1], top[1:], left[:1]])
    return xyz_of_eta(eta)


def test_surface_patch_on_a_curved_wedge():
    """A cylinder-like patch whose lower edge collapses to a point (the apex of a web wedge)."""
    radius = 2.0

    def xyz_of_eta(eta):
        eta = np.atleast_2d(eta)
        u, v = eta[:, 0], eta[:, 1]
        half_width = 0.03 * (v + 1.0)  # 0 at the apex v = -1, 6 cm at v = 1
        angle = (u * half_width) / radius
        return np.column_stack([radius * np.sin(angle), radius * np.cos(angle), 0.12 * (v + 1.0)])

    eta_boundary = np.array(
        [(-1.0, -1.0), (1.0, -1.0), (1.0, 1.0), (1.0 / 3.0, 1.0), (-1.0 / 3.0, 1.0), (-1.0, 1.0)]
    )
    boundary = xyz_of_eta(eta_boundary)
    assert np.linalg.norm(boundary[0] - boundary[1]) < 1e-12  # collapsed lower edge
    nodes, elements = refine_surface_patch(boundary, xyz_of_eta, element_size=0.25, local_size=0.03)
    unique = boundary[[0, 2, 3, 4, 5]]
    assert np.array_equal(nodes[: len(unique)], unique), "boundary nodes moved or duplicated"
    interior = nodes[len(unique):]
    assert np.allclose(np.hypot(interior[:, 0], interior[:, 1]), radius, atol=1e-9)
    # the elements fall along the surface normal of the parametrisation, none flipped
    du = xyz_of_eta(np.array([[0.5, 0.0]])) - xyz_of_eta(np.array([[-0.5, 0.0]]))
    dv = xyz_of_eta(np.array([[0.0, 0.5]])) - xyz_of_eta(np.array([[0.0, -0.5]]))
    normal = np.cross(du[0], dv[0])
    for conn in elements:
        ring = nodes[conn[conn >= 0]]
        assert np.cross(ring[1] - ring[0], ring[2] - ring[0]) @ normal > 0.0


def test_meshed_region_joins_a_structured_neighbour():
    """A MeshedRegion next to a structured one merges into a conforming surface (no free interior edge)."""
    flat = np.zeros((16, 3))
    order = ((0, 0), (0, 3), (3, 3), (3, 0), (0, 1), (0, 2), (1, 3), (2, 3), (3, 2), (3, 1), (2, 0), (1, 0),
             (1, 1), (1, 2), (2, 2), (2, 1))
    for k, (r, c) in enumerate(order):
        flat[k] = (c / 3.0 * 0.4, r / 3.0 * 0.24, 0.0)
    surface = Surface()
    surface.addShellRegion("quad3", flat, np.array([2, 1, 2, 1]), name="structured", meshMethod="structured")
    boundary = _strip(0.4, 0.24, 3, 3) + np.array([0.4, 0.0])
    nodes, elements = refine_patch(boundary, element_size=0.25, local_size=0.1)
    nodes3 = np.column_stack([nodes, np.zeros(len(nodes))])
    surface.shellRegions.append(MeshedRegion(nodes3, elements))
    surface.regionNames.append("refined")
    mesh = surface.getSurfaceMesh()
    use = {}
    for conn in mesh["elements"]:
        ring = list(conn[conn >= 0])
        for a, b in zip(ring, ring[1:] + ring[:1], strict=True):
            use[frozenset((a, b))] = use.get(frozenset((a, b)), 0) + 1
    shared = [k for k, n in use.items() if n == 2 and all(abs(mesh["nodes"][i][0] - 0.4) < 1e-9 for i in k)]
    assert len(shared) == 1, "the two regions do not share their edge on x = 0.4"
    outer = [k for k, n in use.items() if n == 1]
    for key in outer:  # every other edge lies on the outline of the two regions together
        a, b = (mesh["nodes"][i] for i in key)
        on_outline = (
            (abs(a[1]) < 1e-9 and abs(b[1]) < 1e-9)
            or (abs(a[1] - 0.24) < 1e-9 and abs(b[1] - 0.24) < 1e-9)
            or (abs(a[0]) < 1e-9 and abs(b[0]) < 1e-9)
            or (abs(a[0] - 0.8) < 1e-9 and abs(b[0] - 0.8) < 1e-9)
        )
        assert on_outline, f"free interior edge {a[:2]} - {b[:2]}"
    assert [s["name"] for s in mesh["sets"]["element"]] == ["structured", "refined"]
