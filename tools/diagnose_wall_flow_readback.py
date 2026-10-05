"""Read the applied nodal field back, wall by wall, for one blade ring.

If the multi-cell flow were converted to nodal forces correctly, then for every
wall edge the force contributed by one cell must be q_i * ell along the wall
tangent, so the edge's net force is (sum of the owning cells' signed flows) * ell
along the tangent. This measures that, plus the section's total moment and net
force, on a real blade ring.
"""

import numpy as np

import tests.validation.blade.test_blade_rated_twist as t
from aeroelast.models.blade.model import Blade
from aeroelast.solvers.bem.force_projection import (
    ForceProjector,
    multi_cell_shear_flow,
    ring_section,
)
from tests.support.openfast_bem import build_blade_aero_from_aerodyn

SPAN = np.array([0.0, 0.0, 1.0])
T_REF = 1.0e6

blade = Blade(str(t.YAML), element_size=1.0)
blade.generate_mesh()
mesh = blade.mesh
if mesh is None:
    raise SystemExit("Blade.generate_mesh() produced no mesh")
props = blade.get_element_properties()
coords = mesh.coords_array
aero = build_blade_aero_from_aerodyn(t.AD_PRIMARY)
proj = ForceProjector(mesh, aero, span_direction=SPAN, element_properties=props)

k = len(proj._strips) // 2
strip = proj._strips[k]
groups = [g for g in proj._strip_ring_groups[k] if len(g) > 0]
sections = proj._strip_ring_sections[k]
means = [float(((strip.centroid + strip.offsets[g]) @ SPAN).mean()) for g in groups]
gi = int(np.argmin([abs(m - strip.r_center) for m in means]))
global_nodes = strip.node_indices[groups[gi]]
pts = coords[global_nodes]

cells, adjacency, edge_length, edge_S = ring_section(mesh, global_nodes, SPAN, props)
flows, theta = multi_cell_shear_flow(cells, adjacency, edge_length, edge_S, T_REF)

centroid = pts.mean(axis=0)
d = pts - centroid
n = len(pts)
f = np.zeros((n, 3))
edge_force: dict[tuple[int, int], np.ndarray] = {}
for i, cell in enumerate(cells):
    for a, b in cell.edges:
        pa, pb = pts[a], pts[b]
        e = pb - pa
        ell = float(np.linalg.norm(e))
        if ell == 0.0:
            continue
        f[a] += 0.5 * flows[i] * e
        f[b] += 0.5 * flows[i] * e
        key = (min(a, b), max(a, b))
        edge_force[key] = edge_force.get(key, np.zeros(3)) + flows[i] * e

M = float(np.cross(d, f).sum(axis=0) @ SPAN)
F = f.sum(axis=0)
print(f"ring: {n} nodes, {len(cells)} cells, {len(edge_force)} wall edges")
print(f"solver flows q_i [N/m]: {np.round(flows, 3).tolist()}")
print(f"total moment about the ring centroid / T_REF = {M / T_REF:+.12f}")
print(f"net force |sum f| / T_REF = {np.linalg.norm(F) / T_REF:.3e}")

print("\nper-wall read-back, grouped by how many cells own the wall:")
owners_hist = {1: [], 2: []}
for key, vec in edge_force.items():
    owners = adjacency.get(key, ())
    ell = edge_length[key]
    flow_back = float(np.linalg.norm(vec)) / ell
    q_expected = sum(flows[i] for i in owners)
    # a shared wall is traversed in opposite directions, so the expected flow is
    # q_i - q_j; with the unsigned read-back we compare magnitudes and the sign below
    owners_hist.setdefault(len(owners), []).append((key, flow_back, abs(q_expected), owners))

for nown, rows in sorted(owners_hist.items()):
    if not rows:
        continue
    err = [abs(rb - abs(qe)) / max(abs(qe), 1e-30) for _, rb, qe, _ in rows]
    S = [edge_S.get(kk, float("nan")) for kk, _, _, _ in rows]
    print(
        f"  walls owned by {nown} cell(s): n={len(rows)}  "
        f"|read-back flow| vs |sum q_i|: median {np.median(err):.3%}  max {max(err):.3%}  "
        f"| S range {min(S):.3e}..{max(S):.3e}"
    )
    for kk, rb, qe, owners in rows[:3]:
        print(f"      {kk}: read-back {rb:12.2f}  |sum q_i| {qe:12.2f}  owners {owners}")

# does the direction of the circulating field match the sign of the moment?
print(
    f"\nsign check: total moment {M:+.4e} vs T_REF {T_REF:+.4e} -> "
    f"{'SAME' if np.sign(M) == np.sign(T_REF) else 'OPPOSITE'}"
)
