"""Centrifugal shell resultant: mass/area, never volumetric density.

Run with scripts/aeroenv.sh python -m pytest -q <this file>.
The 2 x 4 m panel has centroid (3, 0, 0); rotation about y gives
Fx = rho*t*8*omega**2*3. Composite ply densities integrate to rho*t.
"""

import numpy as np
import pytest

from aeroelast.core.assembler import MeshAssembler
from aeroelast.core.laminate import create_laminate_from_angles
from aeroelast.core.material import IsotropicMaterial, Material
from aeroelast.core.mesh.entities import ElementSet, ElementType, MeshElement, Node
from aeroelast.core.mesh.model import MeshModel
from aeroelast.core.properties import CompositeShellProperty, ShellProperty
from aeroelast.elements import ElementFamily


@pytest.mark.parametrize("path", ["legacy", "property", "dict", "composite"])
@pytest.mark.parametrize("thickness", [0.05, 0.2])
def test_centrifugal_resultant_uses_areal_density(path, thickness):
    mesh = MeshModel()
    nodes = [Node(list(xyz), geometric_node=False)
             for xyz in [(2., -2., 0.), (4., -2., 0.), (4., 2., 0.), (2., 2., 0.)]]
    for node in nodes:
        mesh.add_node(node)
    mesh.add_element(MeshElement(nodes=nodes, element_type=ElementType.quad))
    mesh.add_element_set(ElementSet("panel", set(mesh.elements)))
    rho, omega = 2700., 3.
    material = IsotropicMaterial(name="metal", E=70e9, nu=.3, rho=rho)
    config = {"element_family": ElementFamily.SHELL}
    if path == "legacy":
        config.update(material=material, thickness=thickness)
    else:
        if path == "composite":
            orthotropic = Material(name="ply", E=(70e9, 10e9, 10e9),
                                   G=(5e9, 5e9, 5e9), nu=(.3, .3, .3), rho=rho)
            prop = CompositeShellProperty(create_laminate_from_angles(
                orthotropic, thickness / 2, [0., 90.]))
        elif path == "dict":
            prop = {"type": "isotropic", "e": 70e9, "nu": .3,
                    "rho": rho, "thickness": thickness}
        else:
            prop = ShellProperty(material, thickness)
        config["properties"] = {"panel": prop}
    assembler = MeshAssembler(mesh, {"elements": config})
    load = assembler.assemble_centrifugal_load(omega, np.array([0., 1., 0.]), np.zeros(3))
    expected = rho * thickness * 8 * omega**2 * 3
    np.testing.assert_allclose(load.reshape(-1, 6).sum(axis=0),
                               [expected, 0., 0., 0., 0., 0.], rtol=1e-12, atol=1e-10)
    # Existing zero-speed behavior and reference-axis translation must survive.
    np.testing.assert_array_equal(
        assembler.assemble_centrifugal_load(0., np.array([0., 1., 0.]), np.zeros(3)), 0.)
    centered = assembler.assemble_centrifugal_load(
        omega, np.array([0., 1., 0.]), np.array([3., 0., 0.]))
    np.testing.assert_allclose(centered.reshape(-1, 6).sum(axis=0), 0., atol=1e-10)
