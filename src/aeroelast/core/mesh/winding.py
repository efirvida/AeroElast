"""Element-winding canonicalisation for shell meshes.

The composite assembler rotates every ply by a per-element angle offset
(:func:`element_angle_offset` in the Rust mesh crate): a *signed* angle
measured about the element normal, and the normal comes from the node
winding with no canonicalisation.  Meshers that emit mixed windings (blade
skins vs shear webs, triangles vs quads) therefore rotate the plies of some
elements the wrong way -- invisible for isotropic materials and for meshes
whose first edge is aligned with the span (offset ±90°, equivalent modulo
180), but real for oblique elements with off-axis plies.

Strategy (robust, no mesh-topology traversal):

1. Score each element with the physical "outward" test of a closed section:
   ``n · (element_centroid − section_centroid)``, projected perpendicular to
   the span axis.  Elements with ``|score| > 0.5`` are trusted (typically the
   skins) and are oriented so the normal points outward -- the blade's layup
   stacking starts at the outer mould line.
2. Every untrusted element (shear webs, TE band, tip cap, elements whose
   normal is nearly parallel to the span) inherits the orientation of the
   nearest trusted element of its own cross-section.

Open surfaces (a flat plate) have no intrinsic outward direction: no element
is trusted and the mesh is left untouched -- such meshes must define their
normal convention by construction.
"""

from __future__ import annotations

import numpy as np

__all__ = ["canonicalize_windings"]


def canonicalize_windings(mesh, span_axis: int = 2, decimals: int = 6) -> int:
    """Re-order nodes so element normals follow one consistent orientation.

    Returns the number of elements whose node order was reversed.
    """
    coords = np.asarray(mesh.coords_array, dtype=float)
    nid2idx = mesh.node_id_to_index
    span = np.zeros(3)
    span[span_axis] = 1.0

    # ── section centroids (group nodes by their rounded span coordinate) ────
    s = coords[:, span_axis]
    uniq, inv = np.unique(np.round(s, decimals), return_inverse=True)
    counts = np.bincount(inv, minlength=len(uniq)).astype(float)
    cent = np.zeros((len(uniq), 3))
    for d in range(3):
        cent[:, d] = np.bincount(inv, weights=coords[:, d], minlength=len(uniq)) / counts

    elements = list(mesh.elements)
    n_el = len(elements)

    # ── per-element: normal, outward score, section id, in-plane centroid ───
    score = np.zeros(n_el)
    section = np.full(n_el, -1, dtype=int)
    c_inplane = np.zeros((n_el, 3))
    for i, el in enumerate(elements):
        idx = [nid2idx[nd.id] for nd in el.nodes]
        p = coords[idx]
        if len(p) < 3:
            continue
        c = p.mean(axis=0)
        sec = int(inv[idx[0]])
        section[i] = sec
        c_inplane[i] = c
        n = np.cross(p[1] - p[0], p[2] - p[0])
        nn = float(np.linalg.norm(n))
        if nn < 1e-14:
            continue
        n /= nn
        out = c - cent[sec]
        out -= out[span_axis] * span
        no = float(np.linalg.norm(out))
        if no < 1e-9:
            continue
        score[i] = float(np.dot(n, out / no))

    # ── trusted orientations, then inherit from the nearest trusted element ──
    orient = np.zeros(n_el, dtype=bool)
    trusted = np.nonzero(np.abs(score) > 0.5)[0]
    for i in trusted:
        orient[i] = score[i] > 0.0

    flipped = 0
    for sec in np.unique(section[section >= 0]):
        sec_idx = np.nonzero(section == sec)[0]
        trusted_sec = sec_idx[np.abs(score[sec_idx]) > 0.5]
        if trusted_sec.size == 0:
            continue  # open surface / degenerate section: leave untouched
        pts = c_inplane[trusted_sec]
        for i in sec_idx:
            if abs(score[i]) > 0.5:
                continue
            d = np.linalg.norm(pts - c_inplane[i], axis=1)
            orient[i] = orient[trusted_sec[int(np.argmin(d))]]
    for i in np.nonzero(section >= 0)[0]:
        if not orient[i]:
            elements[i].nodes = list(reversed(elements[i].nodes))
            flipped += 1
    return flipped
