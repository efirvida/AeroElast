"""Which end of a blade ring is the leading edge? Measured from the skin only.

`_section_ends` calls the blunter end the leading edge, using the in-plane spread
of the nodes inside the outer quarter of the chord. On this mesh that slab also
contains the WEB nodes (which sit across the thickness where the webs meet the
skin) and the trailing edge is expanded, so the spread is not an airfoil
thickness and the call can invert. An inverted call puts the aerodynamic centre
on the wrong side of the moment reference, which flips the sign of the lever-arm
transfer - measured at +7.607e5 N.m nose-up against the polars' -2.766e5 N.m.

This diagnostic identifies the leading edge **independently**: the skin outline's
thickness is maximal at 25-35% of the chord from the leading edge for any real
aerofoil, so the end nearer the maximum-thickness station is the leading edge.
Skin edges are the ones a single cell owns; web edges are owned by two.
"""

import numpy as np

import tests.validation.blade.test_blade_rated_twist as t
from aeroelast.models.blade.model import Blade
from aeroelast.solvers.bem.force_projection import ForceProjector, ring_section
from tests.support.openfast_bem import build_blade_aero_from_aerodyn

SPAN = np.array([0.0, 0.0, 1.0])

blade = Blade(str(t.YAML), element_size=1.0)
blade.generate_mesh()
mesh = blade.mesh
if mesh is None:
    raise SystemExit("no mesh")
props = blade.get_element_properties()
aero = build_blade_aero_from_aerodyn(t.AD_PRIMARY)
proj = ForceProjector(mesh, aero, span_direction=SPAN, element_properties=props)

print(
    f"{'k':>3} {'c[m]':>7} {'max thk/c':>10} {'at x/c':>7} {'LE side':>9} {'_section_ends LE':>17} {'agree':>6}"
)
for k in range(0, len(proj._strips), max(1, len(proj._strips) // 8)):
    strip = proj._strips[k]
    ring = next((g for g in proj._strip_ring_groups[k] if len(g) >= 3), None)
    if ring is None:
        continue
    pts = strip.centroid + strip.offsets[ring]
    ch = proj._strip_chord_dirs[k]
    in_plane = np.cross(SPAN, ch)
    in_plane /= np.linalg.norm(in_plane)
    x = pts @ ch
    q = pts @ in_plane
    lo, hi = float(x.min()), float(x.max())
    chord = hi - lo

    cells, adjacency, edge_length, edge_S = ring_section(
        mesh, strip.node_indices[ring], SPAN, props
    )
    skin_nodes = set()
    for key, owners in adjacency.items():
        if len(owners) == 1:
            skin_nodes.update(key)
    if len(skin_nodes) < 4:
        continue

    # thickness profile of the SKIN outline only, on 10 chordwise bins
    nbins = 10
    prof = np.zeros(nbins)
    edges = (hi - lo) / nbins
    for b in range(nbins):
        sel = [
            i
            for i in range(len(pts))
            if i in skin_nodes and lo + b * edges <= x[i] < lo + (b + 1) * edges
        ]
        if len(sel) >= 2:
            prof[b] = q[sel].max() - q[sel].min()
    bmax = int(np.argmax(prof))
    x_max = lo + (bmax + 0.5) * edges
    frac_hi = (x_max - lo) / chord
    # the leading edge is the end nearer the max-thickness station
    le_side = "lo-end" if frac_hi < 0.5 else "hi-end"

    le_i, te_i = ForceProjector._section_ends(pts, ch, SPAN)
    picked = "hi-end" if le_i == int(np.argmax(x)) else "lo-end"
    print(
        f"{k:3d} {chord:7.2f} {prof[bmax] / chord:10.4f} {frac_hi:7.3f} "
        f"{le_side:>9} {picked:>17} {str(le_side == picked):>6}"
    )
