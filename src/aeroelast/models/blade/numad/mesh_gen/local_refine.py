"""Local remesh of a shell patch whose boundary nodes are given and must not move.

The blade mesher builds every region as a structured grid in its natural plane ``eta`` (nodes
uniform on each edge, triangles where the two chordwise edges carry different counts). Where the
counts of a narrow piece disagree, or where the neighbouring row cuts the cap at another set of
lines, that grid has no way to follow the nodes the neighbours dictate on a shared line, and the
coordinate merge then leaves a crack (or the transition collapses into slivers).

``refine_patch`` meshes such a patch in the plane with gmsh: the boundary nodes are imposed
exactly (one element per boundary segment, no node is added or moved on the boundary) and the
interior size grows from ``local_size`` at the boundary to ``element_size`` with a bounded growth
ratio. ``refine_surface_patch`` does the same for a patch on a curved region: it unfolds the
boundary onto the mean plane, meshes there and maps the interior nodes back with the region's own
parametrisation. Neither function holds a length of its own: every size comes from the caller.
"""

from __future__ import annotations

from contextlib import contextmanager

import numpy as np

_GMSH_FRONTAL_DELAUNAY = 6  # Mesh.Algorithm: Frontal-Delaunay, triangles (paired into quads on request)
_GMSH_TRIANGLE = 2  # element type codes of gmsh
_GMSH_QUAD = 3


@contextmanager
def _gmsh_model():
    """A scratch gmsh model that leaves no global state behind.

    gmsh is initialised here only when nobody has done it; an already running session keeps its
    current model and its options, and is not finalised.
    """
    import gmsh

    owner = not gmsh.isInitialized()
    if owner:
        gmsh.initialize(interruptible=False)
    names = (
        "General.Terminal",
        "Mesh.Algorithm",
        "Mesh.RecombineAll",
        "Mesh.MeshSizeFromPoints",
        "Mesh.MeshSizeExtendFromBoundary",
        "Mesh.MeshSizeFromCurvature",
    )
    saved = {name: gmsh.option.getNumber(name) for name in names}
    previous = None if owner else gmsh.model.getCurrent()
    gmsh.option.setNumber("General.Terminal", 0)
    name = "local_refine"
    while name in gmsh.model.list():
        name += "_"
    gmsh.model.add(name)
    try:
        yield gmsh
    finally:
        gmsh.model.remove()
        if owner:
            gmsh.finalize()
        else:
            for name, value in saved.items():
                gmsh.option.setNumber(name, value)
            if previous:
                gmsh.model.setCurrent(previous)


def _signed_area(loop):
    x, y = loop[:, 0], loop[:, 1]
    return 0.5 * float(np.dot(x, np.roll(y, -1)) - np.dot(np.roll(x, -1), y))


def _mesh_once(boundary, element_size, local_size, growth, recombine):
    """One gmsh run: ``(nodes (n, 2), elements (e, 4))``, boundary nodes first."""
    m = len(boundary)
    with _gmsh_model() as gmsh:
        geo = gmsh.model.geo
        points = [geo.addPoint(x, y, 0.0, element_size) for x, y in boundary]
        lines = [geo.addLine(points[k], points[(k + 1) % m]) for k in range(m)]
        surface = geo.addPlaneSurface([geo.addCurveLoop(lines)])
        geo.synchronize()
        for line in lines:
            gmsh.model.mesh.setTransfiniteCurve(line, 2)  # one element: the boundary is imposed

        field = gmsh.model.mesh.field
        distance = field.add("Distance")
        field.setNumbers(distance, "CurvesList", lines)
        field.setNumber(distance, "Sampling", 2 * m)
        ramp = field.add("Threshold")
        field.setNumber(ramp, "InField", distance)
        field.setNumber(ramp, "SizeMin", local_size)
        field.setNumber(ramp, "SizeMax", element_size)
        field.setNumber(ramp, "DistMin", 0.0)
        field.setNumber(ramp, "DistMax", (element_size - local_size) / (growth - 1.0))
        field.setAsBackgroundMesh(ramp)
        gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
        gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
        gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
        gmsh.option.setNumber("Mesh.Algorithm", _GMSH_FRONTAL_DELAUNAY)
        gmsh.option.setNumber("Mesh.RecombineAll", 1 if recombine else 0)
        gmsh.model.mesh.generate(2)

        tags, coords, _ = gmsh.model.mesh.getNodes()
        position = {int(t): c[:2] for t, c in zip(tags, coords.reshape(-1, 3), strict=True)}
        boundary_tags = [int(gmsh.model.mesh.getNodes(0, p)[0][0]) for p in points]
        interior_tags = sorted(set(position) - set(boundary_tags))
        order = {t: k for k, t in enumerate(boundary_tags + interior_tags)}
        nodes = np.array([position[t] for t in boundary_tags + interior_tags])
        nodes[:m] = boundary  # exactly the imposed nodes, not their gmsh copy

        elements = []
        types, _, node_tags = gmsh.model.mesh.getElements(2, surface)
        for kind, flat in zip(types, node_tags, strict=True):
            size = 3 if kind == _GMSH_TRIANGLE else 4
            for conn in np.asarray(flat, dtype=int).reshape(-1, size):
                elements.append([order[int(t)] for t in conn] + [-1] * (4 - size))
    elements = np.array(elements, dtype=int).reshape(-1, 4)
    for e, conn in enumerate(elements):
        ring = nodes[conn[conn >= 0]]
        if _signed_area(ring) < 0.0:
            elements[e, : len(ring)] = conn[: len(ring)][::-1]
    return nodes, elements


def worst_corner(nodes, elements):
    """Smallest sine of an element corner angle (the scaled Jacobian of a planar element)."""
    worst = 1.0
    for conn in elements:
        ring = nodes[conn[conn >= 0]]
        ahead = np.roll(ring, -1, axis=0) - ring
        behind = np.roll(ring, 1, axis=0) - ring
        cross = ahead[:, 0] * behind[:, 1] - ahead[:, 1] * behind[:, 0]
        worst = min(worst, float(np.min(np.abs(cross) / (np.linalg.norm(ahead, axis=1) * np.linalg.norm(behind, axis=1)))))
    return worst


def refine_patch(boundary, element_size, local_size, growth=1.3, recombine=False, trials=(1.0, 0.5, 0.75, 1.5, 2.0)):
    """Mesh the interior of a planar patch whose boundary nodes are imposed.

    Parameters
    ----------
    boundary : (m, 2) array
        The boundary nodes in counter-clockwise order, without repeating the first one. Each pair
        of consecutive nodes is one element edge: no node is added or moved on the boundary.
    element_size : float
        Size of the regular mesh around the patch; the interior never grows beyond it.
    local_size : float
        Interior size next to the boundary (for a narrow strip, a fraction of its width).
    growth : float
        Bounded growth ratio ``dh/dd > 0`` of the interior size with the distance ``d`` from the
        boundary: ``h(d) = local_size + growth_rate * d`` up to ``element_size``, where
        ``growth_rate = growth - 1``.
    recombine : bool
        Pair triangles into quadrilaterals where gmsh can (mixed meshes are allowed).
    trials : sequence of float
        Factors on ``local_size`` to try. The Delaunay insertion of gmsh is not monotonic in the
        size (a finer size can leave a worse triangle next to a long boundary edge), so the mesh
        with the best worst corner (``worst_corner``) of the trials is returned.

    Returns
    -------
    nodes : (n, 2) array
        The ``m`` boundary nodes first, unchanged and in the given order, then the interior ones.
    elements : (e, 4) int array
        Counter-clockwise connectivity, triangles padded with -1.
    """
    boundary = np.asarray(boundary, dtype=float)
    if boundary.ndim != 2 or boundary.shape[1] != 2 or len(boundary) < 3:
        raise ValueError("boundary must be an (m >= 3, 2) array")
    if not (element_size > 0.0 and local_size > 0.0):
        raise ValueError("element_size and local_size must be positive")
    if not growth > 1.0:
        raise ValueError("growth must exceed 1 (it is the growth ratio of the interior size)")
    if _signed_area(boundary) <= 0.0:
        raise ValueError("boundary must be counter-clockwise")
    best = None
    for factor in trials:
        size = min(local_size * factor, element_size)
        nodes, elements = _mesh_once(boundary, element_size, size, growth, recombine)
        quality = worst_corner(nodes, elements)
        if best is None or quality > best[0]:
            best = (quality, nodes, elements)
    return best[1], best[2]


def unfold(points, u_direction, v_direction):
    """Orthonormal frame ``(origin, basis (2, 3))`` of the plane that best follows a patch.

    The first axis follows ``u_direction`` and the second is ``v_direction`` made orthogonal to
    it, so the plane has the orientation of the patch parametrisation ``(u, v)``.
    """
    points = np.asarray(points, dtype=float)
    first = np.asarray(u_direction, dtype=float)
    first = first / np.linalg.norm(first)
    second = np.asarray(v_direction, dtype=float)
    second = second - np.dot(second, first) * first
    second = second / np.linalg.norm(second)
    return points.mean(axis=0), np.vstack((first, second))


def locate(targets, plane_of_eta, grid=9, tolerance=1.0e-10, iterations=30):
    """Parameters ``eta`` in [-1, 1]^2 whose image ``plane_of_eta(eta)`` is each target point.

    Newton iteration with a finite-difference Jacobian, started from the nearest of a ``grid`` x
    ``grid`` sampling of the patch; ``tolerance`` is relative to the patch extent.
    """
    targets = np.asarray(targets, dtype=float)
    axis = np.linspace(-1.0, 1.0, grid)
    samples = np.array([(u, v) for v in axis for u in axis])
    image = plane_of_eta(samples)
    extent = float(np.ptp(image, axis=0).max())
    eta = np.empty_like(targets)
    step = 1.0 / (grid * grid)
    for k, target in enumerate(targets):
        guess = samples[np.argmin(np.linalg.norm(image - target, axis=1))].copy()
        for _ in range(iterations):
            here = plane_of_eta(guess[None, :])[0]
            residual = target - here
            if np.linalg.norm(residual) <= tolerance * extent:
                break
            jac = np.empty((2, 2))
            for axis_id in range(2):
                bumped = guess.copy()
                bumped[axis_id] += step
                jac[:, axis_id] = (plane_of_eta(bumped[None, :])[0] - here) / step
            guess = np.clip(guess + np.linalg.solve(jac, residual), -1.0, 1.0)
        else:
            raise ValueError("locate: no convergence for a patch interior node")
        eta[k] = guess
    return eta


def refine_surface_patch(
    boundary_xyz,
    xyz_of_eta,
    element_size,
    local_size,
    growth=1.3,
    recombine=False,
    coincide=1.0e-6,
    record=None,
    replay=None,
):
    """``refine_patch`` for a patch on a curved region ``xyz_of_eta`` (eta in [-1, 1]^2 -> 3D).

    ``boundary_xyz`` are the boundary nodes in 3D, counter-clockwise in ``eta``; consecutive
    nodes closer than ``coincide`` times the patch extent are one node (a collapsed edge, the apex
    of a wedge). They are kept exactly. The interior nodes are meshed in the unfolded plane and
    placed on the region at the ``eta`` of their plane position.

    ``record`` (a list) receives ``(keep, eta, elements, count)``: the indices of ``boundary_xyz``
    that were kept, the ``eta`` of the interior nodes, the connectivity and the boundary node
    count. Passing that tuple back as ``replay`` rebuilds the same topology on another boundary of
    the same node count (another parametrisation of the same patch) without running gmsh: only the
    node positions follow it. A boundary with another node count is meshed afresh.

    Returns ``(nodes (n, 3), elements (e, 4))``: the unique boundary nodes first, then interior.
    """
    boundary_xyz = np.asarray(boundary_xyz, dtype=float)
    extent = float(np.ptp(boundary_xyz, axis=0).max())
    keep = [0]
    for k in range(1, len(boundary_xyz)):
        if np.linalg.norm(boundary_xyz[k] - boundary_xyz[keep[-1]]) > coincide * extent:
            keep.append(k)
    if len(keep) > 1 and np.linalg.norm(boundary_xyz[keep[-1]] - boundary_xyz[0]) <= coincide * extent:
        keep.pop()
    if replay is not None and len(boundary_xyz) == replay[3]:
        keep, eta, elements, _ = replay
        loop = boundary_xyz[keep]
        return np.vstack((loop, xyz_of_eta(np.asarray(eta, dtype=float).reshape(-1, 2)))), np.array(elements)
    loop = boundary_xyz[keep]

    half = np.array([[-0.5, 0.0], [0.5, 0.0], [0.0, -0.5], [0.0, 0.5]])
    end = xyz_of_eta(half)
    origin, basis = unfold(loop, end[1] - end[0], end[3] - end[2])

    def plane_of_eta(eta):
        return (xyz_of_eta(eta) - origin) @ basis.T

    nodes, elements = refine_patch((loop - origin) @ basis.T, element_size, local_size, growth, recombine)
    nb = len(loop)
    out = np.empty((len(nodes), 3))
    out[:nb] = loop
    eta = np.empty((0, 2))
    if len(nodes) > nb:
        eta = locate(nodes[nb:], plane_of_eta)
        out[nb:] = xyz_of_eta(eta)
    if record is not None:
        record.append((np.array(keep), eta, elements.copy(), len(boundary_xyz)))
    return out, elements


class MeshedRegion:
    """A region whose mesh is already made; ``Surface.getSurfaceMesh`` takes it like a ShellRegion."""

    def __init__(self, nodes, elements):
        self.nodes = np.asarray(nodes, dtype=float)
        self.elements = np.asarray(elements, dtype=int)

    def createShellMesh(self):
        return {"nodes": self.nodes.copy(), "elements": self.elements.copy()}
