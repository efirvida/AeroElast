"""The production ring grouping of a gmsh blade mesh, for tests (gmsh-official-mesher T14).

The gmsh mesh has no constant-z node rings (25 012 nodes, 24 952 distinct z at 6 dp), so a test that groups nodes by exact z
or by a z strip does not see a section.  Production groups them through ``physical_ring_labels`` (span gap + in-ring element
edges, ``RING_EDGE_SLOPE``) with the mesh-scaled gap ``RING_GAP_FRACTION * span_extent`` (``ForceProjector.__init__``,
which feeds ``_strip_ring_groups`` and ``BEMFSIParticipant._find_strip_rings``), and weights a ring's nodes by
``section_contour.ring_tributary_lengths``.  This module only calls those functions; it holds no estimator of its own.
"""

from __future__ import annotations

import numpy as np

from aeroelast.solvers.bem.force_projection import RING_GAP_FRACTION, mesh_corner_edges, physical_ring_labels
from aeroelast.solvers.bem.section_contour import ring_tributary_lengths

SPAN_Z = np.array([0.0, 0.0, 1.0])


def ring_labels(mesh, span_dir=SPAN_Z) -> np.ndarray:
    """Physical ring label of every mesh node (``physical_ring_labels`` with the production gap tolerance)."""
    span = np.asarray(mesh.coords_array, float) @ np.asarray(span_dir, float)
    return physical_ring_labels(mesh, span, RING_GAP_FRACTION * float(span.max() - span.min()))


def mesh_rings(mesh, span_dir=SPAN_Z, min_nodes: int = 1) -> list[np.ndarray]:
    """The physical rings (global node indices), in increasing mean span; rings with fewer than ``min_nodes`` are dropped."""
    labels = ring_labels(mesh, span_dir)
    span = np.asarray(mesh.coords_array, float) @ np.asarray(span_dir, float)
    groups = [np.flatnonzero(labels == lab) for lab in np.unique(labels)]
    groups = [g for g in groups if len(g) >= min_nodes]
    groups.sort(key=lambda g: float(span[g].mean()))
    return groups


def skin_wall_edges(mesh) -> np.ndarray:
    """Corner edges of the skin: every element edge except the shear webs' (``*_SW`` sets), as the torsion guard builds them."""
    edges = mesh_corner_edges(mesh)
    web_elements = [e for name, members in mesh.element_sets.items() if name.endswith("_SW") for e in members.elements]
    web = mesh_corner_edges(mesh, web_elements)
    if web.size == 0:
        return edges
    n = len(mesh.nodes)
    return edges[~np.isin(edges[:, 0] * n + edges[:, 1], web[:, 0] * n + web[:, 1])]


def ring_ids_and_weights(mesh, span_dir=SPAN_Z, edges=None):
    """``(ring_id, weight)`` per mesh node: ring index in ``mesh_rings`` order (-1: none), skin tributary length."""
    coords = np.asarray(mesh.coords_array, float)
    edges = skin_wall_edges(mesh) if edges is None else edges
    ring_id = np.full(len(coords), -1)
    weights = np.zeros(len(coords))
    for k, ring in enumerate(mesh_rings(mesh, span_dir)):
        ring_id[ring] = k
        weights[ring] = ring_tributary_lengths(ring, edges, coords)
    return ring_id, weights


def aerodyn_blade_aero():
    """The AeroDyn-deck ``BladeAero`` of the IEA-15 reference case (the strip grid production builds its projector on)."""
    from tests.support.openfast_bem import build_blade_aero_from_aerodyn
    from tests.support.paths import DATA_DIR

    return build_blade_aero_from_aerodyn(
        DATA_DIR / "reference" / "iea15mw_openfast" / "case" / "IEA-15-240-RWT_AeroDyn15.dat"
    )


def strip_station_rings(mesh, blade_aero=None, span_dir=SPAN_Z) -> list[np.ndarray]:
    """The ring production assigns to every BEM strip, as global node indices (``BEMFSIParticipant._find_strip_rings``).

    Built on the production object itself: ``ForceProjector`` groups every strip's nodes into physical rings
    (``_strip_ring_groups``) and the participant takes, per strip, the group holding the strip node nearest the strip
    centre.  Production never uses a ring that no strip centre selects, so a guard of "the section the solver reads"
    iterates these rings, not every ``physical_ring_labels`` group (a gmsh mesh carries partial groups of 1-25 nodes at
    z 19.3, 27.4-27.8 m and a few more; they are not sections).
    """
    from aeroelast.solvers.bem.force_projection import ForceProjector

    aero = aerodyn_blade_aero() if blade_aero is None else blade_aero
    projector = ForceProjector(mesh, aero, span_direction=np.asarray(span_dir, float))
    mesh_span = np.asarray(mesh.coords_array, float) @ np.asarray(span_dir, float)
    rings: list[np.ndarray] = []
    for strip, groups in zip(projector._strips, projector._strip_ring_groups, strict=True):
        idx = strip.node_indices
        if len(idx) == 0:
            continue
        anchor = int(np.argmin(np.abs(mesh_span[idx] - strip.r_center)))
        group = next(g for g in groups if anchor in g)
        rings.append(idx[np.sort(group)])
    return rings
