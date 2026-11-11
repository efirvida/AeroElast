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


COLLAPSE_EDGE_M = 1e-9  # physical edge length treated as a geometric collapse
MIN_RING_DIVISIONS = 16  # about 45 degrees of turning per element on average


def _centroid(parameters, elements):
    """Mean chart coordinate over each element's true corners (triangles repeat one)."""
    live = (elements >= 0)[:, :, None]
    return (parameters*live).sum(axis=1)/live.sum(axis=1)


def _has_material(cells):
    return any(cell.thickness_m > 0 for cell in cells)


@dataclass(frozen=True)
class SkinMesh:
    nodes: np.ndarray  # metres, NuMAD axes
    elements: np.ndarray  # zero-based four-node shells
    parameters: np.ndarray  # native (eta, nd_arc); seam represented by nd_arc=0
    element_parameters: np.ndarray  # per-corner chart coordinates, preserving nd=1 at seam
    cells: tuple  # source laminate sampled at each element's parametric centroid
    eta_bounds: tuple
    attachments: dict = field(default_factory=dict)
    layout: tuple = None  # (stations, {(e0, e1): sub-stations}) reused by the web charts


def generate_skin_mesh(topology, element_size=.25, eta_bounds=(0., 1.), min_strip_m=0.):
    return _generate_chart(topology, element_size, eta_bounds, 'skin', min_strip_m=min_strip_m)


def _generate_chart(topology, element_size, eta_bounds, surface, layout=None, min_strip_m=0.):
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
    if layout is not None:
        stations = list(layout[0])  # the skin's stations, so web rows land on skin nodes
    substations = {}

    def precision_error(detail):
        raise ValueError('Unsupported extended chart subdivision precision: ' + detail
                         + '; distinct source intervals were not merged; '
                         'this is a representation limit, not invalid source geometry')

    LE = ('le', 'landmark', 'arc0')

    @lru_cache(None)
    def idmap(eta):
        result = {identity: p.nd_arc for p in section(eta).points for identity in p.identities}
        if surface == 'skin':
            # The loft's leading-edge landmark (arc 0) is a mandatory ring cut:
            # without it one element wraps the highest-curvature region.
            points = section(eta).points
            for p1, p2 in zip(points[:-1], points[1:], strict=True):
                if p1.arc*p2.arc < 0:
                    result[LE] = p1.nd_arc+(p2.nd_arc-p1.nd_arc)*rational((0-p1.arc)/(p2.arc-p1.arc))
        return result

    def perimeter(*etas):
        return max(abs(section(e).points[-1].arc-section(e).points[0].arc) for e in etas)

    def count(width, reference):
        return max(1, int(np.ceil(extended(width)*reference/element_size)))

    # Rows: the ordered nd positions of every station. Inside one event strip no
    # identity changes order, so every interval keeps ONE subdivision count for
    # the whole strip and neighbouring rows conform without any transport. An
    # event station is a row of its own (its section carries the identities of
    # both sides); the strips next to it join it with transition elements.
    if surface == 'skin':
        # A moving boundary that sweeps over the leading-edge landmark is an
        # event of its own: add a station where it crosses (secant on the loft).
        for _ in range(4):
            added = []
            for e0, e1 in zip(stations[:-1], stations[1:], strict=True):
                m0, m1 = idmap(e0), idmap(e1)
                if LE not in m0 or LE not in m1:
                    continue
                for k in m0.keys() & m1.keys() - {LE}:
                    d0, d1 = m0[k]-m0[LE], m1[k]-m1[LE]
                    if d0*d1 < 0:
                        t = rational(extended(d0)/extended(d0-d1))
                        guess = e0+(e1-e0)*t
                        for _ in range(3):
                            m = idmap(guess)
                            if LE not in m:
                                break
                            d = m[k]-m[LE]
                            slope = extended((d1-d0)/(e1-e0))
                            guess = guess-rational(extended(d)/slope)
                            guess = min(max(guess, e0), e1)
                        # The landmark is computed, not a source event: a root
                        # within a quarter element of a station joins it, since a
                        # strip thinner than that cannot hold a usable element.
                        margin = rational(element_size/4/length)
                        if e0+margin < guess < e1-margin:
                            added.append(guess)
            if not added:
                break
            stations = sorted({*stations, *added})
    if min_strip_m > 0 and layout is None:
        # Resolution limit: two stations closer than min_strip_m of span are one
        # station. The strip between them is thinner than any element can be and
        # the sign of its Jacobian is numerical noise. The first one is kept.
        kept = [stations[0]]
        for e in stations[1:-1]:
            if (e-kept[-1])*length >= min_strip_m and (stations[-1]-e)*length >= min_strip_m:
                kept.append(e)
        stations = kept+[stations[-1]]
    rows, strips, event_strips, finest = {}, [], [], {}
    for e in stations:
        points = section(e).points
        rows[e] = {p.nd_arc for p in points}
    for e0, e1 in zip(stations[:-1], stations[1:], strict=True):
        midpoint = e0+(e1-e0)/2
        if not e0 < midpoint < e1:
            precision_error(f'no interior span sample in [{e0!r}, {e1!r}]')
        mid = section(midpoint)
        if not mid.points:
            continue
        if surface != 'skin' and not _has_material(mid.skin_cells):
            continue  # a web whose laminate is zero thick here does not exist
        reference = perimeter(e0, midpoint, e1)
        for use_le in (True, False):
            intervals = []
            for left, right in mid.intervals:
                keys = (left.identities[0], right.identities[0])
                width = max(idmap(e)[keys[1]]-idmap(e)[keys[0]] for e in (e0, e1))
                n = count(width, reference)
                if (use_le and LE in idmap(midpoint) and idmap(midpoint)[keys[0]] < idmap(midpoint)[LE] < idmap(midpoint)[keys[1]]
                        and all(LE in idmap(e)
                                and idmap(e)[LE]-idmap(e)[keys[0]] >= rational(element_size/4/perimeter(e))
                                and idmap(e)[keys[1]]-idmap(e)[LE] >= rational(element_size/4/perimeter(e))
                                for e in (e0, e1))):
                    # Same interval over the whole strip: split it at the landmark.
                    n = max(n, 2)
                    w1 = max(idmap(e)[LE]-idmap(e)[keys[0]] for e in (e0, e1))
                    n1 = min(n-1, max(1, int(round(n*extended(w1)/extended(width)))))
                    intervals.append(((keys[0], LE), n1))
                    intervals.append(((LE, keys[1]), n-n1))
                else:
                    intervals.append((keys, n))
            if surface == 'skin':
                # A closed profile needs a minimum angular resolution: below
                # MIN_RING_DIVISIONS a coarse size on a small chord wraps one
                # element around the leading edge or collapses the ring.
                widths = [max(idmap(e)[k[1]]-idmap(e)[k[0]] for e in (e0, e1)) for k, _ in intervals]
                counts = [n for _, n in intervals]
                while sum(counts) < MIN_RING_DIVISIONS:
                    counts[max(range(len(counts)), key=lambda i: widths[i]/counts[i])] += 1
                intervals = [(k, n) for (k, _), n in zip(intervals, counts, strict=True)]
            # A boundary that slides along the contour must not shear a row of
            # elements beyond ~45 degrees: the span step also bounds the slide.
            slide = max(abs(idmap(e1)[k]-idmap(e0)[k]) for pair, _ in intervals for k in pair)
            n_span = max(1, int(np.ceil((e1-e0)*length/element_size)),
                         int(np.ceil(extended(slide)*reference/element_size)))
            sub = layout[1][(e0, e1)] if layout is not None else linspace(e0, e1, n_span+1)
            # The landmark split holds only if the order is kept at every row.
            if all(idmap(s).get(k[0], 0) <= idmap(s).get(k[1], 1) and k[0] in idmap(s) and k[1] in idmap(s)
                   for s in sub for k, _ in intervals):
                break
        substations[(e0, e1)] = sub
        for s in sub[1:-1]:
            row = {idmap(s)[intervals[0][0][0]]}
            for (first, last), n in intervals:
                row.update(linspace(idmap(s)[first], idmap(s)[last], n+1))
            rows[s] = row
        for e in (e0, e1):
            for (first, last), n in intervals:
                finest[(e, first, last)] = max(n, finest.get((e, first, last), 0))
        event_strips.append((midpoint, sub))
        for a, b in zip(sub[:-1], sub[1:], strict=True):
            for first, last in (k for k, _ in intervals):
                strips.append((a, b, idmap(a)[first], idmap(a)[last], idmap(b)[first], idmap(b)[last]))

    # An event station takes the finer count of its two sides per interval, so
    # a neighbouring strip with the same count meets it without hanging nodes.
    for (e, first, last), n in finest.items():
        rows[e].update(linspace(idmap(e)[first], idmap(e)[last], n+1))

    uv, elements, ids = [], [], {}

    def node(p):
        p = tuple(p)
        if p not in ids:
            ids[p] = len(uv)
            uv.append(p)
        return ids[p]

    def average(*points):
        return tuple(np.mean(points, axis=0))

    def span_row(eta, lo_q, hi_q):
        return [(eta, q) for q in sorted(rows[eta]) if lo_q <= q <= hi_q]

    def fractions(row, lo_q, hi_q):
        width = extended(hi_q-lo_q)
        return [0. if width == 0 else float(extended(q-lo_q)/width) for _, q in row]

    def append(corners):
        if len(set(corners)) != len(corners):
            precision_error(f'subdivision collapsed distinct corners {corners!r}')
        elements.append([node(v) for v in corners])

    for a, b, al, ar, bl, br in strips:
        if al == ar and bl == br:
            continue
        lower, upper = span_row(a, al, ar), span_row(b, bl, br)
        if not lower or not upper:
            precision_error(f'strip {(a, b, al, ar, bl, br)!r} lost its boundary nodes '
                            f'(lower {len(lower)}, upper {len(upper)})')
        tl, tu = fractions(lower, al, ar), fractions(upper, bl, br)
        i = j = 0
        # Merge the two rows: a quad where the next cuts agree to within a
        # fraction of the local spacing, otherwise one transition triangle.
        while i < len(lower)-1 or j < len(upper)-1:
            if i == len(lower)-1:
                append([lower[i], upper[j], upper[j+1]]); j += 1
            elif j == len(upper)-1:
                append([lower[i], upper[j], lower[i+1]]); i += 1
            else:
                step = min(tl[i+1]-tl[i], tu[j+1]-tu[j])
                if abs(tl[i+1]-tu[j+1]) <= .3*step:
                    append([lower[i], upper[j], upper[j+1], lower[i+1]]); i += 1; j += 1
                elif tl[i+1] < tu[j+1]:
                    append([lower[i], upper[j], lower[i+1]]); i += 1
                else:
                    append([lower[i], upper[j], upper[j+1]]); j += 1

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
    # An edge of (numerically) zero physical length is a geometric collapse of
    # the source, e.g. a zero-width TE leg: its two chart nodes are one point.
    # Merge exactly those nodes (edge < COLLAPSE_EDGE_M, far below any band) and
    # let the element lose a corner; an element that loses its area disappears.
    parent = list(range(len(nodes)))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for e in elements:
        for u, v in zip(e, e[1:]+e[:1], strict=True):
            if uv[u][0] != uv[v][0]:
                continue  # only along a section: a thin span strip is a band, not a collapse
            a_, b_ = u, v
            u, v = remap[u], remap[v]
            if (u != v and np.linalg.norm(nodes[u]-nodes[v]) < COLLAPSE_EDGE_M
                    and abs(float(uv[a_][1]-uv[b_][1])) > 1e-6):  # a real chart gap, zero length on the loft
                parent[find(max(u, v))] = find(min(u, v))
    remap = np.array([find(int(remap[i])) for i in range(len(nodes))])
    kept, kept_ids, cells = [], [], []
    for e in elements:
        mapped = [remap[i] for i in e]
        keep = [k for k in range(len(e)) if mapped[k] not in mapped[:k]]
        if len(keep) < 3:
            continue
        kept.append([mapped[k] for k in keep])
        kept_ids.append([e[k] for k in keep])
        eta, q = average(*(uv[e[k]] for k in keep))
        cells.append(next(c for c in section(eta).skin_cells if c.nd_bounds[0] <= q <= c.nd_bounds[1]))
    padded = np.array([r+[-1]*(4-len(r)) for r in kept_ids], dtype=int)
    live = padded >= 0
    connectivity = np.array([r+[-1]*(4-len(r)) for r in kept], dtype=int)
    used = np.unique(connectivity[live])
    compact = np.full(len(nodes), -1, dtype=int)
    compact[used] = np.arange(len(used))
    connectivity = np.where(live, compact[np.where(live, connectivity, 0)], -1)
    nodes, parameters = nodes[used], parameters[used]
    if not np.isfinite(nodes).all():
        raise ValueError('Nonpositive skin quadrilateral on source geometry')
    for row in connectivity:
        ids4 = row[row >= 0]
        xyz = nodes[ids4]
        if len(set(ids4)) != len(ids4):
            raise ValueError('Nonpositive skin quadrilateral on source geometry')
        if len(ids4) == 4:
            bad = np.einsum('i,i', np.cross(xyz[1]-xyz[0], xyz[3]-xyz[0]),
                            np.cross(xyz[2]-xyz[1], xyz[0]-xyz[1])) <= 0
        else:
            bad = not np.linalg.norm(np.cross(xyz[1]-xyz[0], xyz[2]-xyz[0])) > 0
        if bad:
            raise ValueError('Nonpositive skin quadrilateral on source geometry')
    attachments = defaultdict(dict)
    if surface == 'skin':
        for midpoint, sub in event_strips:
            for name, _ in section(midpoint).webs:
                if not _has_material(c for c in topology.section(midpoint).web_cells
                                     if c.surface == name):
                    continue
                for handle in ('start_nd_arc', 'end_nd_arc'):
                    for eta in sub:
                        nd = idmap(eta)[('web', name, handle)]
                        attachments[(name, handle)][eta] = int(compact[remap[ids[(eta, nd)]]])
    chains = {key: tuple(mapping[e] for e in sorted(mapping)) for key, mapping in attachments.items()}
    corner_uv = extended_array(uv)
    corner_uv = np.stack([corner_uv[np.where(r >= 0, r, r[0])] for r in padded])
    return SkinMesh(nodes, connectivity, parameters, corner_uv,
                    tuple(cells), (extended(lo), extended(hi)), chains,
                    (tuple(stations), substations))


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


def generate_mesh(topology, element_size=.25, eta_bounds=(0., 1.), *, with_materials=False,
                  min_strip_m=0., collapse_ratio=None, collapse_max_shift=0.05, collapse_optimize=True):
    """Isolated skin+web raw mesh; not the BladeMesh/MeshModel dispatch.

    Exact source attachment identities share skin nodes and whole edges. The
    geometric basis is suitable for explicit physical-vector angle conversion.
    ``with_materials=True`` opts into the restricted S12 SONATA-derived policy;
    unsupported web stacks raise ValueError. Default geometry-only behavior and
    source cells are unchanged. OML/reference nodes are never relocated, except by the opt-in
    ``collapse_ratio`` (e.g. 0.15): strips thinner than ``collapse_ratio * element_size`` are cut
    out and end in triangles (``windio2_collapse``); that moves the merged lines by up to a few
    percent of the element size and drops the thin bands' elements, which the caller must accept. The
    cumulative shift of a merged line is capped at ``collapse_max_shift`` and, with ``collapse_optimize``, merged
    apexes slide along their own line (up to 0.6 es) to make the closing triangles near-equilateral
    (``collapse_opt``). Recommended: ``collapse_ratio=0.25``, ``collapse_max_shift=0.05``.
    """
    from .windio2_basis import shell_bases
    topology = topology.with_extended_precision()
    skin = generate_skin_mesh(topology, element_size, eta_bounds, min_strip_m)
    nodes, elements = list(skin.nodes), list(skin.elements)
    cells, surfaces = list(skin.cells), ['skin']*len(skin.cells)
    parameters = list(skin.element_parameters)
    for name in dict(topology.section(sum(eta_bounds)/2).webs) | {
            key[0]: None for key in skin.attachments}:
        if (name, 'start_nd_arc') not in skin.attachments:
            continue
        web = _generate_chart(topology, element_size, eta_bounds, name, skin.layout)
        mapping = []
        shared = [{skin.parameters[i, 0]: i for i in skin.attachments[(name, handle)]}
                  for handle in ('start_nd_arc', 'end_nd_arc')]
        for xyz, (eta, q) in zip(web.nodes, web.parameters, strict=True):
            if q == 0. or q == 1.:
                mapping.append(shared[int(q)][eta])
            else:
                mapping.append(len(nodes))
                nodes.append(xyz)
        mapped = np.asarray(mapping)[np.where(web.elements >= 0, web.elements, 0)]
        elements.extend(np.where(web.elements >= 0, mapped, -1))
        parameters.extend(web.element_parameters)
        cells.extend(web.cells)
        surfaces.extend([name]*len(web.cells))
    nodes, elements = np.asarray(nodes), np.asarray(elements, dtype=int)
    mesh = RawMesh(nodes, elements, tuple(cells), tuple(surfaces), np.asarray(parameters),
                   shell_bases(nodes, elements), skin.attachments)
    if collapse_ratio:
        from .collapse_opt import collapse_and_optimize_raw_mesh
        mesh, _ = collapse_and_optimize_raw_mesh(mesh, element_size, collapse_ratio,
                                                 max_shift=collapse_max_shift, optimize=collapse_optimize)
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

    eta = _centroid(mesh.element_parameters, mesh.elements)[:, 0]
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
    tri = mesh.elements[:, 3] < 0
    du = np.where(tri[:, None], uv[:, 1]-uv[:, 0], uv[:, 1]-uv[:, 0]+uv[:, 2]-uv[:, 3])
    dv = np.where(tri[:, None], uv[:, 2]-uv[:, 0], uv[:, 3]-uv[:, 0]+uv[:, 2]-uv[:, 1])
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
