"""Isolated T11 skin/web mesher and opt-in material mapping. No BladeMesh dispatch.

Mesh in (native eta, source nd_arc), then evaluate every node on the actual
completed loft. Internal triangles are split into three four-distinct-node
quads; triangles/degenerate quads are never returned to the shell adapter.
"""
from collections import defaultdict, deque
from dataclasses import dataclass, field, replace
from types import SimpleNamespace

from aeroelast.models.blade.numad.objects.windio2_topology import BoundaryPoint
from aeroelast.models.blade.numad.objects.windio2_precision import (
    Pchip, extended, extended_array, linspace, rational,
)
from functools import lru_cache

import numpy as np


@dataclass(frozen=True)
class SkinMesh:
    nodes: np.ndarray  # metres, NuMAD axes
    elements: np.ndarray  # zero-based four-node shells
    parameters: np.ndarray  # native (eta, nd_arc); seam represented by nd_arc=0
    element_parameters: np.ndarray  # per-corner chart coordinates, preserving nd=1 at seam
    cells: tuple  # source laminate sampled at each element's parametric centroid
    eta_bounds: tuple
    attachments: dict = field(default_factory=dict)


def generate_skin_mesh(topology, element_size=.25, eta_bounds=(0., 1.)):
    return _generate_chart(topology, element_size, eta_bounds, 'skin')


def _generate_chart(topology, element_size, eta_bounds, surface):
    """Return boundary-conforming skin ONLY, including empty-laminate cells.

    ``element_size`` is a resolution target, not a maximum-edge guarantee on a
    warped loft. Native knots/crossings/dropoffs are mandatory stations. Source
    arc widths and reference length set extra subdivisions. All horizontal
    interfaces are reconciled before quadrangulation, including partial-domain
    boundaries. No geometric welding tolerance or material averaging is used.

    Connectivity follows positive (eta, nd_arc) orientation: increasing span
    crossed with source suction->LE->pressure contour gives the skin outward
    normal. T12 must preserve this physical direction. Layer angles in cells
    remain SOURCE angles; the opt-in raw material adapter stores separate angles.
    """
    if not np.isfinite(element_size) or not .15 <= element_size <= .25:
        raise ValueError('element_size must be in [0.15, 0.25] m')
    lo, hi = eta_bounds
    if not np.isfinite([lo, hi]).all() or not 0 <= lo < hi <= 1:
        raise ValueError('eta_bounds must satisfy 0 <= lo < hi <= 1')
    topology = topology.with_extended_precision()
    lo, hi = rational(lo), rational(hi)

    @lru_cache(None)
    def section(eta):
        raw = topology.section(eta)
        if surface == 'skin':
            return raw
        cells = tuple(c for c in raw.web_cells if c.surface == surface)
        if not cells:
            return SimpleNamespace(points=(), intervals=(), skin_cells=())
        positions = topology._positions(eta, surface)
        positions.update({('perimeter', surface, 'start'): rational(0),
                          ('perimeter', surface, 'end'): rational(1)})
        ends = dict(raw.webs)[surface]
        start, end = np.array(ends[0].xyz), np.array(ends[1].xyz)
        width = np.linalg.norm(end-start)
        points = tuple(BoundaryPoint(q, extended(q)*width,
                                     tuple(start*(1-extended(q))+end*extended(q)),
                                     tuple(k for k, v in positions.items() if v == q))
                       for q in sorted(set(positions.values())))
        return SimpleNamespace(points=points, intervals=tuple(zip(points[:-1], points[1:], strict=True)),
                               skin_cells=cells)

    # The definition's span is the authoritative curvilinear reference length.
    length = topology._definition.span[-1]
    stations = sorted({lo, hi, *(e for e in topology.event_stations if lo < e < hi)})
    refined = [stations[0]]
    for a, b in zip(stations[:-1], stations[1:], strict=True):
        n = max(1, int(np.ceil((b-a)*length/element_size)))
        refined.extend(linspace(a, b, n+1)[1:])
    rings = {e: {p.nd_arc for p in section(e).points} for e in refined}
    strips = []

    def precision_error(detail):
        raise ValueError('Unsupported extended chart subdivision precision: ' + detail
                         + '; distinct source intervals were not merged; '
                         'this is a representation limit, not invalid source geometry')

    for a, b in zip(refined[:-1], refined[1:], strict=True):
        midpoint = a+(b-a)/2
        if not a < midpoint < b:
            precision_error(f'no interior span sample in [{a!r}, {b!r}]')
        mid = section(midpoint)
        if not mid.points:
            continue
        ends = [{identity: p.nd_arc for p in section(e).points for identity in p.identities}
                for e in (a, b)]
        perimeter = max(abs(section(e).points[-1].arc-section(e).points[0].arc)
                        for e in (a, (a+b)/2, b))
        for left, right in mid.intervals:
            # Inside one event strip identities cannot change order. Their
            # endpoint values include analytically unified crossing identities.
            bounds = [[mapping[p.identities[0]] for p in (left, right)] for mapping in ends]
            n = max(1, int(np.ceil(extended(max(y-x for x, y in bounds))*perimeter/element_size)))
            qa, qb = [linspace(x, y, n+1) for x, y in bounds]
            rings[a].update(qa)
            rings[b].update(qb)
            strips.extend((a, b, qa[j], qa[j+1], qb[j], qb[j+1]) for j in range(n))

    # Reconcile hanging cuts by transport across each noncollapsed strip.
    # Exact normalized positions keep both horizontal edges compatible, even
    # when neighboring stations requested different contour subdivisions.
    # Transport reaches a fixed point before adding any element: repairing one
    # side after quadrangulation would leave a crack in its neighbor.
    incident = defaultdict(list)
    for index, (a, b, *_) in enumerate(strips):
        incident[a].append(index)
        incident[b].append(index)
    pending = deque(range(len(strips)))
    queued = set(pending)
    while pending:
        index = pending.popleft()
        queued.remove(index)
        a, b, al, ar, bl, br = strips[index]
        if al == ar or bl == br:
            # An exact crossing has a genuinely collapsed edge. Its polygon
            # keeps the existing triangle-to-quad construction below.
            continue
        fractions = {(q-al)/(ar-al) for q in rings[a] if al <= q <= ar}
        fractions.update((q-bl)/(br-bl) for q in rings[b] if bl <= q <= br)
        for eta, left, right in ((a, al, ar), (b, bl, br)):
            extra = {left+t*(right-left) for t in fractions}-rings[eta]
            if extra:
                rings[eta].update(extra)
                for neighbor in incident[eta]:
                    if neighbor not in queued:
                        pending.append(neighbor)
                        queued.add(neighbor)

    uv, elements, ids = [], [], {}

    def node(p):
        p = tuple(p)
        if p not in ids:
            ids[p] = len(uv)
            uv.append(p)
        return ids[p]

    def average(*points):
        return tuple(np.mean(points, axis=0))

    def append_quad(quad):
        if len(set(quad)) != 4:
            precision_error(f'subdivision collapsed distinct corners in strip '
                            f'{(a, b, al, ar, bl, br)!r}')
        elements.append([node(v) for v in quad])

    for a, b, al, ar, bl, br in strips:
        if al < ar and bl < br:
            fractions = sorted((q-al)/(ar-al) for q in rings[a] if al <= q <= ar)
            for left, right in zip(fractions[:-1], fractions[1:], strict=True):
                # Pair equal strip coordinates, not diagonal contour samples.
                # Subdivide the tensor quad to retain the same midpoint nodes
                # used by collapsed-edge neighbors and web attachments.
                polygon = [(a, al+left*(ar-al)), (b, bl+left*(br-bl)),
                           (b, bl+right*(br-bl)), (a, al+right*(ar-al))]
                center = average(*polygon)
                for j, corner in enumerate(polygon):
                    prev, nex = polygon[(j-1) % 4], polygon[(j+1) % 4]
                    append_quad([corner, average(corner, nex), center, average(prev, corner)])
            continue
        # CCW parameter polygon. Include ALL cuts on shared horizontal edges;
        # adjacent strips may have different subdivision counts or identities.
        polygon = [(a, q) for q in sorted(rings[a]) if al <= q <= ar]
        polygon += [(b, q) for q in sorted(rings[b], reverse=True) if bl <= q <= br]
        # Our desired order is eta then nd (reverse the above boundary order).
        polygon.reverse()
        center = average(*polygon)
        for p, q in zip(polygon, polygon[1:]+polygon[:1], strict=True):
            # A collapsed source interval contributes no horizontal edge, but
            # the fan's side triangles remain valid at the crossing event.
            if p == q:
                continue
            triangle = (center, p, q)
            centroid = average(*triangle)
            for j, corner in enumerate(triangle):
                prev, nex = triangle[(j-1) % 3], triangle[(j+1) % 3]
                quad = [corner, average(corner, nex), centroid, average(prev, corner)]
                append_quad(quad)

    # Evaluate on the true loft, not on linear chords between mesh vertices.
    parameters = extended_array(uv)
    if len({tuple(p) for p in parameters}) != len(uv):
        precision_error('distinct rational chart vertices share a longdouble image')
    nodes = np.empty((len(uv), 3), dtype=np.longdouble)
    groups = defaultdict(list)
    for i, (eta, _) in enumerate(uv):
        groups[eta].append(i)
    for eta, indices in groups.items():
        if surface == 'skin':
            nodes[indices] = topology.skin_points(eta, [uv[i][1] for i in indices])
        else:
            points = section(eta).points
            start, end = np.array(points[0].xyz), np.array(points[-1].xyz)
            q = parameters[indices, 1, None]
            nodes[indices] = start*(1-q)+end*q

    # Identify only the exact topological TE seam. No spatial tolerance welds.
    remap = np.arange(len(nodes))
    for (eta, q), i in ids.items():
        if surface == 'skin' and q == 1.:
            remap[i] = ids[(eta, 0.)]
    connectivity = remap[np.asarray(elements, dtype=int)]
    cells = []
    for quad in elements:
        eta, q = average(*(uv[i] for i in quad))
        cell = next(c for c in section(eta).skin_cells if c.nd_bounds[0] <= q <= c.nd_bounds[1])
        cells.append(cell)
    used = np.unique(connectivity)
    compact = np.empty(len(nodes), dtype=int)
    compact[used] = np.arange(len(used))
    connectivity = compact[connectivity]
    nodes, parameters = nodes[used], parameters[used]
    xyz = nodes[connectivity]
    cross1 = np.cross(xyz[:, 1]-xyz[:, 0], xyz[:, 3]-xyz[:, 0])
    cross2 = np.cross(xyz[:, 2]-xyz[:, 1], xyz[:, 0]-xyz[:, 1])
    if (not np.isfinite(nodes).all() or np.any(np.einsum('ij,ij->i', cross1, cross2) <= 0)
            or any(len(set(q)) != 4 for q in connectivity)):
        raise ValueError('Nonpositive skin quadrilateral on source geometry')
    attachments = defaultdict(dict)
    if surface == 'skin':
        for a, b in zip(refined[:-1], refined[1:], strict=True):
            for name, _ in section((a+b)/2).webs:
                ends = [dict(section(e).webs)[name] for e in (a, b)]
                for side, handle in enumerate(('start_nd_arc', 'end_nd_arc')):
                    p, q = (a, ends[0][side].nd_arc), (b, ends[1][side].nd_arc)
                    for eta, nd in (p, average(p, q), q):
                        attachments[(name, handle)][eta] = int(compact[remap[ids[(eta, nd)]]])
    chains = {key: tuple(mapping[e] for e in sorted(mapping)) for key, mapping in attachments.items()}
    return SkinMesh(nodes, connectivity, parameters, extended_array(uv)[elements],
                    tuple(cells), (extended(lo), extended(hi)), chains)


@dataclass(frozen=True)
class RawMesh:
    nodes: np.ndarray
    elements: np.ndarray
    cells: tuple
    surfaces: tuple
    element_parameters: np.ndarray
    bases: np.ndarray  # rows: geometric shell x, y, normal; NOT source material axes
    attachments: dict
    material_stacks: tuple | None = None

    @property
    def material_angles(self):
        """Element angles in mapped bottom-to-top order, or None if not requested."""
        if self.material_stacks is None:
            return None
        return tuple(tuple(layer.angle_deg for layer in stack.layers)
                     for stack in self.material_stacks)


def generate_mesh(topology, element_size=.25, eta_bounds=(0., 1.), *, with_materials=False):
    """Isolated skin+web raw mesh; not the BladeMesh/MeshModel dispatch.

    Exact source attachment identities share skin nodes and whole edges. The
    geometric basis is suitable for explicit physical-vector angle conversion.
    ``with_materials=True`` opts into the restricted S12 SONATA-derived policy;
    unsupported web stacks raise ValueError. Default geometry-only behavior and
    source cells are unchanged. OML/reference nodes are never relocated.
    """
    from .windio2_basis import shell_bases
    topology = topology.with_extended_precision()
    skin = generate_skin_mesh(topology, element_size, eta_bounds)
    nodes, elements = list(skin.nodes), list(skin.elements)
    cells, surfaces = list(skin.cells), ['skin']*len(skin.cells)
    parameters = list(skin.element_parameters)
    for name in dict(topology.section(sum(eta_bounds)/2).webs) | {
            key[0]: None for key in skin.attachments}:
        if (name, 'start_nd_arc') not in skin.attachments:
            continue
        web = _generate_chart(topology, element_size, eta_bounds, name)
        mapping = []
        shared = [{skin.parameters[i, 0]: i for i in skin.attachments[(name, handle)]}
                  for handle in ('start_nd_arc', 'end_nd_arc')]
        for xyz, (eta, q) in zip(web.nodes, web.parameters, strict=True):
            if q == 0. or q == 1.:
                mapping.append(shared[int(q)][eta])
            else:
                mapping.append(len(nodes))
                nodes.append(xyz)
        elements.extend(np.asarray(mapping)[web.elements])
        parameters.extend(web.element_parameters)
        cells.extend(web.cells)
        surfaces.extend([name]*len(web.cells))
    nodes, elements = np.asarray(nodes), np.asarray(elements, dtype=int)
    mesh = RawMesh(nodes, elements, tuple(cells), tuple(surfaces), np.asarray(parameters),
                   shell_bases(nodes, elements), skin.attachments)
    return map_mesh_materials(topology, mesh) if with_materials else mesh


def map_mesh_materials(topology, mesh):
    """Return raw mesh with physical stack positions and element-frame angles.

    Reuse Geometry's native section-frame construction, not a parallel rotation
    recipe. Canonical normals follow positive chart orientation independently
    of the caller's element bases; web normals are oriented toward local LE.
    This is not the public BladeMesh adapter and performs no constitutive solve.
    """
    from aeroelast.models.blade.numad.objects.geometry import Geometry
    from .windio2_basis import map_stack, shell_bases

    eta = mesh.element_parameters[:, :, 0].mean(axis=1)
    definition = topology._definition
    frames = Geometry()
    frames.ispan = eta * definition.span[-1]
    twist = definition.windio2_outer_shape['twist']
    frames._windio2_eta = eta
    frames._v2_pchip = Pchip
    frames.idegreestwist = Pchip(twist['grid'], twist['values'])(eta)
    frames._build_windio2_frames(definition)
    canonical = shell_bases(mesh.nodes, mesh.elements)[:, 2]
    # Signed chart Jacobian, centred to avoid cancellation in tiny event strips.
    uv = mesh.element_parameters
    du = uv[:, 1]-uv[:, 0]+uv[:, 2]-uv[:, 3]
    dv = uv[:, 3]-uv[:, 0]+uv[:, 2]-uv[:, 1]
    determinant = du[:, 0]*dv[:, 1]-du[:, 1]*dv[:, 0]
    if np.any(determinant == 0):
        raise ValueError('Degenerate material mapping chart')
    canonical *= np.sign(determinant)[:, None]
    stacks = []
    for cell, frame, normal, basis in zip(mesh.cells, frames._windio2_frames,
                                          canonical, mesh.bases, strict=True):
        if cell.surface != 'skin':
            toward_le = -frame[:, 1]
            alignment = normal@toward_le
            if abs(alignment) <= 1e-12:
                raise ValueError('Web normal has no unambiguous LE side')
            normal = normal * np.sign(alignment)
        stacks.append(map_stack(cell, frame[:, 2], normal, basis))
    return replace(mesh, material_stacks=tuple(stacks))
