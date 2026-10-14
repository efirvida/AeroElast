"""Diagnose which tip-deflection measure matches the D-Tube beam reference.

Solves the same linear problem as tests/validation/benchmarks/test_dtube_camarena_bending.py and reports the
tip uy as: arc-node mean (the test's measure), web-node mean, all-node mean, and the
section-centroid-equivalent (area-weighted over the tip ring). Read-only on the test.
"""
import importlib.util
import sys

import numpy as np

path = sys.argv[1]
spec = importlib.util.spec_from_file_location("dt", path)
dt = importlib.util.module_from_spec(spec)
spec.loader.exec_module(dt)

from aeroelast.core.bc import DirichletCondition
from aeroelast.core.material import IsotropicMaterial
from aeroelast.core.properties import ShellProperty
from aeroelast.elements import ElementFamily
from aeroelast.solvers.elasticity.static_linear import StaticLinearSolver

mesh, arc_nodes, arc_thetas = dt._build_dtube_mesh()
prop = ShellProperty(material=IsotropicMaterial(name="Steel", E=dt.E, nu=dt.NU, rho=7850.0), thickness=dt.T)
cfg = {"solver": {}, "elements": {"element_family": ElementFamily.SHELL,
       "properties": {"tube": prop}, "span_direction": (0.0, 0.0, 1.0)}}
s = StaticLinearSolver(mesh, cfg)
dpn = s.domain.dofs_per_node
cl = sorted(mesh.node_id_to_index[n.id] * dpn + d for n in mesh.get_node_set("clamped").nodes.values() for d in range(dpn))
s.add_dirichlet_conditions([DirichletCondition(cl, 0.0)])
fy = dt._apply_traction(mesh, s, arc_nodes, arc_thetas)
u = np.asarray(s.solve(), dtype=float)
uy = lambda n: u[mesh.node_id_to_index[n.id] * dpn + 1]
tip = [n for n in mesh.nodes if np.isclose(n.coords[2], dt.L)]
arc = [n for n in arc_nodes[dt.NZ] if n is not None]
web = [n for n in tip if np.isclose(n.coords[1], 0.0)]
print(f"ref beam                {dt.DELTA_REF:.4f}")
print(f"resultant / (w L)       {fy/(dt.W_RES*dt.L):.6f}")
for name, ns in (("arc nodes (test)", arc), ("web nodes", web), ("all tip nodes", tip)):
    v = np.mean([uy(n) for n in ns]); print(f"{name:22s}  {v:.4f}  rel {abs(v-dt.DELTA_REF)/dt.DELTA_REF:.2%}")
# segment-length-weighted ring mean (wall-area weighting, trapezoid on the polygon)
ring = [(n, n.coords) for n in tip]
x = np.array([c[0] for _, c in ring]); y = np.array([c[1] for _, c in ring]); v = np.array([uy(n) for n, _ in ring])
# arc polygon + web chords, node tributary length
def trib(nodes):
    return None
wts = np.zeros(len(ring)); idx = {n.id: i for i, (n, _) in enumerate(ring)}
chains = [[n for n in arc_nodes[dt.NZ]]]
xs = sorted(web, key=lambda n: n.coords[0]); chains.append(xs)
for ch in chains:
    for a, b in zip(ch[:-1], ch[1:]):
        l = np.hypot(*(np.array(a.coords[:2]) - np.array(b.coords[:2])))
        wts[idx[a.id]] += l / 2; wts[idx[b.id]] += l / 2
vw = np.sum(wts * v) / np.sum(wts)
print(f"wall-length weighted    {vw:.4f}  rel {abs(vw-dt.DELTA_REF)/dt.DELTA_REF:.2%}")
# mid-span
mid = [n for n in mesh.nodes if np.isclose(n.coords[2], dt.L / 2)]
print(f"midspan all nodes       {np.mean([uy(n) for n in mid]):.4f}")
