"""Production shell fibre mapping against independent plane-stress strain energy.

Unlike direct MITC constructors, StaticLinearSolver prepares laminates from a
requested span direction. A zero-degree ply must follow its in-plane projection
for both triangles and quads. The reference is CLT Q11, not a material transform
from the implementation. The 5% bound is the existing material-suite physical
limit; the uniform membrane patch should in fact be nearly exact.
"""

import numpy as np
import pytest

from aeroelast.core.mesh.entities import ElementSet, ElementType, MeshElement, Node
from aeroelast.core.mesh.model import MeshModel
from aeroelast.elements import ElementFamily
from aeroelast.models.blade.model import build_rust_properties
from aeroelast.solvers.elasticity.static_linear import StaticLinearSolver
from tests.support.assertions import assert_relative_error

E1, E2, NU12 = 120e9, 10e9, 0.3
H, EPS = 0.002, 1e-4
TOL = 0.05


def _patch_energy(kind, *, a=(1., 0., 0.), b=(0., 0., 1.),
                  span=(0., 0., 1.), fibre=(0., 0., 1.),
                  order=None, isotropic=False):
    """Assemble through the production solver; impose affine fibre strain."""
    a, b = np.asarray(a), np.asarray(b)
    coords = np.array([np.zeros(3), a, b] if kind == "tri"
                      else [np.zeros(3), a, a + b, b])
    if order is not None:
        coords = coords[order]
    mesh = MeshModel()
    nodes = [Node(point.tolist()) for point in coords]
    for node in nodes:
        mesh.add_node(node)
    element_type = ElementType.triangle if kind == "tri" else ElementType.quad
    mesh.add_element(MeshElement(nodes, element_type))
    mesh.add_element_set(ElementSet("patch", mesh.elements))
    elastic = ({"E": E1, "nu": NU12} if isotropic else
               {"E": [E1, E2, E2], "G": [5e9, 3e9, 5e9],
                "nu": [NU12, NU12, NU12]})
    props = build_rust_properties({
        "materials": [{"name": "ud", "density": 1500., "elastic": elastic}],
        "sections": [{"elementSet": "patch", "layup": [["ud", H, 0.]]}],
    })
    elements = {"element_family": ElementFamily.SHELL, "properties": props}
    if span is not None:
        elements["span_direction"] = span
    solver = StaticLinearSolver(mesh, {"solver": {}, "elements": elements})
    rows, cols, values = solver.domain._rust.assemble_k()
    fibre = np.asarray(fibre, dtype=float)
    fibre /= np.linalg.norm(fibre)
    u = np.zeros((len(nodes), 6))
    u[:, :3] = EPS * np.outer(coords @ fibre, fibre)
    u = u.ravel()
    energy = float(0.5 * np.sum(values * u[rows] * u[cols]))
    area = np.linalg.norm(np.cross(a, b)) * (0.5 if kind == "tri" else 1.)
    return energy, area


@pytest.mark.parametrize("kind", ["tri", "quad"])
def test_production_zero_ply_follows_span(kind):
    energy, area = _patch_energy(kind)
    expected = 0.5 * area * H * E1 / (1 - NU12**2 * E2 / E1) * EPS**2
    assert_relative_error(
        energy, expected, tol=TOL,
        reference_name="CLT UD plane-stress Q11 uniform-strain energy",
        kind="analytical", what=f"{kind} production +Z fibre energy [J]",
    )


@pytest.mark.parametrize("geometry", [
    pytest.param({"order": [0, 2, 1]}, id="reversed-normal"),
    pytest.param({"order": [1, 2, 0]}, id="diagonal-local-edge"),
    pytest.param({"a": (0., 1., 0.)}, id="global-X-normal-fallback-Z"),
    pytest.param({"a": (2**-0.5, 0., 2**-0.5), "b": (0., 1., 0.),
                  "span": (0., 1., 2**0.5),
                  "fibre": (0.5, 2**-0.5, 0.5)}, id="tilted-projected-span"),
])
def test_triangle_span_uses_material_reference_not_local_edge(geometry):
    energy, area = _patch_energy("tri", **geometry)
    expected = 0.5 * area * H * E1 / (1 - NU12**2 * E2 / E1) * EPS**2
    assert_relative_error(
        energy, expected, tol=TOL,
        reference_name="CLT UD plane-stress Q11 uniform-strain energy",
        kind="analytical", what="triangle projected-span fibre energy [J]",
    )


@pytest.mark.parametrize("kind", ["tri", "quad"])
@pytest.mark.parametrize("span", [None, (0., 1., 0.)], ids=["default", "normal-span"])
def test_existing_reference_without_in_plane_span_is_preserved(kind, span):
    # With no in-plane span, both constructors on this patch use global X.
    energy, area = _patch_energy(kind, span=span, fibre=(1., 0., 0.))
    expected = 0.5 * area * H * E1 / (1 - NU12**2 * E2 / E1) * EPS**2
    assert_relative_error(
        energy, expected, tol=TOL,
        reference_name="CLT UD plane-stress Q11 default-reference energy",
        kind="analytical", what=f"{kind} unchanged default fibre energy [J]",
    )


@pytest.mark.parametrize("kind", ["tri", "quad"])
def test_isotropic_production_span_behavior_is_preserved(kind):
    energy, area = _patch_energy(kind, isotropic=True)
    expected = 0.5 * area * H * E1 / (1 - NU12**2) * EPS**2
    assert_relative_error(
        energy, expected, tol=TOL,
        reference_name="isotropic plane-stress uniform-strain energy",
        kind="analytical", what=f"{kind} unchanged isotropic energy [J]",
    )
