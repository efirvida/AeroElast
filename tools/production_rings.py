"""Node rings of a gmsh blade mesh, grouped exactly as production does (gmsh-official-mesher T15).

The gmsh mesh has no constant-z node rings (24 952 distinct z of 25 012 nodes at 6 dp), so grouping by rounded z or by a z
band does not see a section.  Production (``ForceProjector.__init__`` -> ``_strip_ring_groups``,
``BEMFSIParticipant._find_strip_rings``) labels the nodes with ``physical_ring_labels`` (span gap
``RING_GAP_FRACTION * span_extent`` plus in-ring element edges, ``RING_EDGE_SLOPE``) and weights a ring's nodes by
``section_contour.ring_tributary_lengths``.  This module only calls those functions; it holds no estimator of its own.
"""

from __future__ import annotations

import numpy as np

SPAN_Z = np.array([0.0, 0.0, 1.0])


def ring_labels(mesh, span_dir=SPAN_Z) -> np.ndarray:
    """Physical ring label of every mesh node (production gap tolerance)."""
    from aeroelast.solvers.bem.force_projection import RING_GAP_FRACTION, physical_ring_labels

    span = np.asarray(mesh.coords_array, float) @ np.asarray(span_dir, float)
    return physical_ring_labels(mesh, span, RING_GAP_FRACTION * float(span.max() - span.min()))


def mesh_rings(mesh, span_dir=SPAN_Z, min_nodes: int = 1) -> list[np.ndarray]:
    """Physical rings (global node indices) in increasing mean span; rings below ``min_nodes`` nodes are dropped."""
    labels = ring_labels(mesh, span_dir)
    span = np.asarray(mesh.coords_array, float) @ np.asarray(span_dir, float)
    groups = [np.flatnonzero(labels == lab) for lab in np.unique(labels)]
    groups = [g for g in groups if len(g) >= min_nodes]
    groups.sort(key=lambda g: float(span[g].mean()))
    return groups


def ring_profile(mesh, z_min, span_dir=SPAN_Z, wall_edges=None):
    """``(ring mean z, node indices, wall-length weights)`` of every production ring whose mean z is >= ``z_min``.

    Weights are ``ring_tributary_lengths`` over all element corner edges of the mesh (skin and webs), or over
    ``wall_edges`` when given.
    """
    from aeroelast.solvers.bem.force_projection import mesh_corner_edges
    from aeroelast.solvers.bem.section_contour import ring_tributary_lengths

    coords = np.asarray(mesh.coords_array, float)
    edges = mesh_corner_edges(mesh) if wall_edges is None else wall_edges
    span = coords @ np.asarray(span_dir, float)
    out = []
    for ring in mesh_rings(mesh, span_dir):
        zk = float(span[ring].mean())
        if zk < z_min - 1e-6:
            continue
        w = ring_tributary_lengths(ring, edges, coords)
        keep = w > 0
        out.append((zk, ring[keep], w[keep]))
    return out


_CACHE: dict[int, tuple[object, list[np.ndarray]]] = {}


def cached_rings(mesh) -> list[np.ndarray]:
    """``mesh_rings(mesh)`` memoised per mesh object (the grouping walks every element)."""
    hit = _CACHE.get(id(mesh))
    if hit is None or hit[0] is not mesh:
        hit = (mesh, mesh_rings(mesh))
        _CACHE.clear()
        _CACHE[id(mesh)] = hit
    return hit[1]


def tip_ring(mesh) -> np.ndarray:
    """Node indices of the tip ring: the ``TipNodes`` set when the mesh has one, else the outermost production ring."""
    try:
        node_set = mesh.get_node_set("TipNodes")
    except (KeyError, ValueError, AttributeError):
        node_set = None
    if node_set is not None:
        idx = mesh.node_id_to_index
        return np.array(sorted(idx[nid] for nid in node_set.nodes), dtype=int)
    return cached_rings(mesh)[-1]
