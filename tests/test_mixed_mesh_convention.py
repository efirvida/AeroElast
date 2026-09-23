"""Mixed-mesh regression: MITC3 and MITC4 must share the same physical
rotation convention.

MITC3's covariant shear B-matrix used an inverted director
(``V3 ≈ e3 − θy·e1 + θx·e2``) while MITC4 used the physical one
(``V3 = e3 + θy·e1 − θx·e2``).  Each element was self-consistent, so pure
meshes passed every benchmark; the inconsistency only shows when triangles
and quads are mixed in the same mesh (the same θy means opposite things to
each family).

This test is the regression guard: a cantilever plate meshed with quads,
with triangles, and with mixed rows must give the same tip deflection.
Pre-fix the mixed patterns were off by 27-100%; post-fix they agree within
0.2%.
"""

from __future__ import annotations

import numpy as np
import pytest

from aeroelast.core.bc import DirichletCondition, NodalLoad
from aeroelast.core.material import IsotropicMaterial
from aeroelast.core.mesh.entities import ElementType, MeshElement, Node
from aeroelast.core.mesh.model import ElementSet, MeshModel
from aeroelast.core.properties import ShellProperty
from aeroelast.elements import ElementFamily
from aeroelast.solvers.elasticity.static_linear import StaticLinearSolver

L, B, H = 1.0, 0.1, 0.001
E, NU = 2.1e11, 0.3
P_TIP = 10.0
NX, NY = 20, 4
TOL = 0.01  # 1 % vs the all-quad solution


def _plate(pattern: str) -> MeshModel:
    """Cantilever plate; `pattern` selects which rows are triangulated."""
    mesh = MeshModel()
    xs = np.linspace(0.0, L, NX + 1)
    ys = np.linspace(0.0, B, NY + 1)
    grid: dict[tuple[float, float], Node] = {}
    for y in ys:
        for x in xs:
            n = Node([float(x), float(y), 0.0], geometric_node=False)
            mesh.add_node(n)
            grid[(round(x, 12), round(y, 12))] = n

    def triangulated(row: int) -> bool:
        if pattern == "quads":
            return False
        if pattern == "tris":
            return True
        if pattern == "one_row":
            return row == NY // 2
        if pattern == "alternating":
            return row % 2 == 1
        raise ValueError(pattern)

    for j in range(NY):
        for i in range(NX):
            n00 = grid[(round(xs[i], 12), round(ys[j], 12))]
            n10 = grid[(round(xs[i + 1], 12), round(ys[j], 12))]
            n11 = grid[(round(xs[i + 1], 12), round(ys[j + 1], 12))]
            n01 = grid[(round(xs[i], 12), round(ys[j + 1], 12))]
            if triangulated(j):
                mesh.add_element(MeshElement(nodes=[n00, n10, n11],
                                             element_type=ElementType.triangle))
                mesh.add_element(MeshElement(nodes=[n00, n11, n01],
                                             element_type=ElementType.triangle))
            else:
                mesh.add_element(MeshElement(nodes=[n00, n10, n11, n01],
                                             element_type=ElementType.quad))
    mesh.add_element_set(ElementSet("plate", set(mesh.elements)))
    return mesh


def _tip_deflection(mesh: MeshModel) -> float:
    prop = ShellProperty(
        material=IsotropicMaterial(name="steel", E=E, nu=NU, rho=0.0),
        thickness=H,
    )
    cfg = {
        "solver": {},
        "elements": {
            "element_family": ElementFamily.SHELL,
            "span_direction": (1.0, 0.0, 0.0),
            "properties": {"plate": prop},
        },
    }
    solver = StaticLinearSolver(mesh, cfg)
    dpn = solver.domain.dofs_per_node
    nid2i = mesh.node_id_to_index
    root = [n for n in mesh.nodes if abs(n.coords[0]) < 1e-12]
    solver.add_dirichlet_conditions([DirichletCondition(
        sorted(nid2i[n.id] * dpn + d for n in root for d in range(dpn)), 0.0)])
    tip = [n for n in mesh.nodes if abs(n.coords[0] - L) < 1e-12]
    dofs, vals = [], []
    for n in tip:
        dofs.append(nid2i[n.id] * dpn + 2)
        vals.append(-P_TIP / len(tip))
    solver.add_nodal_loads([NodalLoad(dofs, vals)])
    u = np.asarray(solver.solve(), dtype=np.float64).reshape(-1, dpn)
    return float(np.mean([u[nid2i[n.id], 2] for n in tip]))


@pytest.mark.parametrize("pattern", ["tris", "one_row", "alternating"])
def test_mixed_mesh_matches_all_quad(pattern):
    """Triangulated and mixed-row meshes must match the all-quad solution."""
    ref = _tip_deflection(_plate("quads"))
    val = _tip_deflection(_plate(pattern))
    rel = abs(val - ref) / abs(ref)
    assert rel < TOL, (
        f"{pattern}: w_tip={val:.6e} vs all-quads {ref:.6e} (rel {rel:.2%} > "
        f"{TOL:.0%}) — MITC3/MITC4 rotation-convention mismatch?"
    )
