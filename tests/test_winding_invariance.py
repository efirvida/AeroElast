"""Element-winding canonicalisation (ply-orientation invariance).

The per-element ply-angle offset is a SIGNED angle about the element normal,
and the normal comes from the node winding.  Without canonicalisation the
material orientation of an oblique element (triangle, skewed quad) flips when
its node order is reversed, corrupting the extension-shear / bend-twist
coupling of off-axis laminates.

The fix lives in :func:`aeroelast.core.mesh.winding.canonicalize_windings`,
called by the blade mesh generator.  Contract pinned here on a CLOSED section
(where "outward" is physical, unlike a flat plate):

* flipping one element's winding changes the response (documents the bug);
* after canonicalisation the flipped and un-flipped meshes agree;
* canonicalisation is idempotent;
* the blade generator emits canonical windings.
"""

from __future__ import annotations

import numpy as np
import pytest

from aeroelast.core.bc import DirichletCondition, NodalLoad
from aeroelast.core.mesh.entities import ElementType, MeshElement, Node
from aeroelast.core.mesh.model import ElementSet, MeshModel
from aeroelast.core.mesh.winding import canonicalize_windings
from aeroelast.elements import ElementFamily
from aeroelast.solvers.elasticity.static_linear import StaticLinearSolver

E1, E2, NU12, G12, RHO = 4.46e10, 1.70e10, 0.262, 3.27e9, 1900.0
BX, BY, BL, H = 0.4, 0.2, 0.5, 0.002
NX, NY = 4, 4


def _tube(flip_one: bool):
    """Closed triangular tube; all elements are triangles (oblique edges).

    An unbalanced +30 laminate makes the response sensitive to the ply sign,
    and the closed section gives a physical outward normal for the
    canonicalisation seeds.
    """
    from _aeroelast import Laminate, OrthotropicMaterial, Ply

    mesh = MeshModel()
    xs = np.linspace(-BX / 2, BX / 2, NX + 1)
    ys = np.linspace(-BY / 2, BY / 2, NY + 1)
    zs = np.linspace(0.0, BL, NY + 1)
    grid: dict[tuple[float, float, float], Node] = {}

    def node(x, y, z):
        key = (round(x, 12), round(y, 12), round(z, 12))
        if key not in grid:
            grid[key] = Node([float(x), float(y), float(z)], geometric_node=False)
            mesh.add_node(grid[key])
        return grid[key]

    def wall(pts, n_u):
        for j in range(NY):
            for i in range(n_u):
                n00 = node(*pts(i, j), zs[j])
                n10 = node(*pts(i + 1, j), zs[j])
                n11 = node(*pts(i + 1, j + 1), zs[j + 1])
                n01 = node(*pts(i, j + 1), zs[j + 1])
                mesh.add_element(MeshElement(nodes=[n00, n10, n11],
                                             element_type=ElementType.triangle))
                mesh.add_element(MeshElement(nodes=[n00, n11, n01],
                                             element_type=ElementType.triangle))

    wall(lambda i, j: (xs[i], np.full(1, -BY / 2)[0]), NX)
    wall(lambda i, j: (xs[i], np.full(1, +BY / 2)[0]), NX)
    wall(lambda i, j: (np.full(1, -BX / 2)[0], ys[i]), NY)
    wall(lambda i, j: (np.full(1, +BX / 2)[0], ys[i]), NY)

    if flip_one:
        # reverse the LAST element's winding (an oblique diagonal element)
        last = list(mesh.elements)[-1]
        last.nodes = list(reversed(last.nodes))

    mesh.add_element_set(ElementSet("tube", set(mesh.elements)))
    mat = OrthotropicMaterial(E1, E2, E2, G12, G12, G12, NU12, NU12, NU12, RHO)
    plies = [Ply(mat, H / 4, +30.0) for _ in range(4)]
    return mesh, {"tube": Laminate(plies)}


def _solve(mesh, props) -> np.ndarray:
    cfg = {
        "solver": {},
        "elements": {
            "element_family": ElementFamily.SHELL,
            "span_direction": (0.0, 0.0, 1.0),
            "properties": props,
        },
    }
    solver = StaticLinearSolver(mesh, cfg)
    dpn = solver.domain.dofs_per_node
    nid2i = mesh.node_id_to_index
    root = [n for n in mesh.nodes if abs(n.coords[2]) < 1e-12]
    solver.add_dirichlet_conditions([DirichletCondition(
        sorted(nid2i[n.id] * dpn + d for n in root for d in range(dpn)), 0.0)])
    tip = [n for n in mesh.nodes if abs(n.coords[2] - BL) < 1e-12]
    dofs, vals = [], []
    for n in tip:  # axial tip load: +30 laminate couples it to shear/twist
        dofs.append(nid2i[n.id] * dpn + 2)
        vals.append(1.0e2 / len(tip))
    solver.add_nodal_loads([NodalLoad(dofs, vals)])
    return np.asarray(solver.solve(), dtype=np.float64)


def test_winding_flip_changes_response_before_canonicalisation():
    """Documents the bug: without canonicalisation the winding matters."""
    u_ref = _solve(*_tube(flip_one=False))
    u_flip = _solve(*_tube(flip_one=True))
    rel = np.linalg.norm(u_ref - u_flip) / np.linalg.norm(u_ref)
    assert rel > 1e-3, (
        f"expected the flipped winding to change the response (>0.1%); got "
        f"{rel:.2e} — the ply orientation no longer depends on the winding?"
    )


def test_canonicalisation_restores_winding_invariance():
    mesh_ref, props_ref = _tube(flip_one=False)
    mesh_flip, props_flip = _tube(flip_one=True)
    canonicalize_windings(mesh_ref)
    canonicalize_windings(mesh_flip)

    u_ref = _solve(mesh_ref, props_ref)
    u_flip = _solve(mesh_flip, props_flip)
    rel = np.linalg.norm(u_ref - u_flip) / np.linalg.norm(u_ref)
    assert rel < 1e-10, (
        f"responses differ by {rel:.2e} after canonicalisation"
    )


def test_canonicalisation_is_idempotent():
    mesh, _ = _tube(flip_one=True)
    first = canonicalize_windings(mesh)
    second = canonicalize_windings(mesh)
    assert first >= 1
    assert second == 0, "second pass should find nothing to flip"


@pytest.mark.slow
def test_blade_generator_emits_canonical_windings():
    """The blade mesh generator must leave no element disagreeing with the
    propagated orientation (i.e. no inward-pointing skin element)."""
    from aeroelast.core.mesh.generators import BladeMesh

    gen = BladeMesh(
        yaml_file="tests/IEA-15-240-RWT.yaml",
        airfoil_dir="tests/airfoils",
        element_size=4.0,
        n_samples=120,
        airfoil_spacing="constant",
        span_grading="chord",
    )
    mesh = gen.generate(renumber="rcm", verbose=False)
    assert canonicalize_windings(mesh) == 0
