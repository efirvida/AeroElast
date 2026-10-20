"""Geometry-only tests of the optional blade skin reference surface."""

from collections import defaultdict
from pathlib import Path

import numpy as np
import pytest

from aeroelast.core.mesh.generators import BladeMesh
from aeroelast.models.blade.model import build_rust_properties

BLADE_YAML = Path(__file__).parent / "IEA-15-240-RWT.yaml"


def generator(**kwargs):
    return BladeMesh(yaml_file=str(BLADE_YAML), element_size=0.25, **kwargs)


@pytest.fixture(scope="module")
def oml():
    gen = generator()
    mesh = gen.generate(renumber="rcm", verbose=False)
    return mesh, build_rust_properties(gen.numad_mesh_data)


def test_default_is_identical_to_explicit_oml(oml):
    explicit = generator(skin_reference="oml").generate(renumber="rcm", verbose=False)
    np.testing.assert_array_equal(explicit.coords_array, oml[0].coords_array)
    assert [e.node_ids for e in explicit.elements] == [e.node_ids for e in oml[0].elements]


def test_mid_insets_skin_only_without_inverting_elements(oml):
    original, props = oml
    mid = generator(skin_reference="mid").generate(renumber="rcm", verbose=False)
    x0, x1 = original.coords_array, mid.coords_array
    idx = original.node_id_to_index
    outer = {idx[n] for n in original.get_node_set("allOuterShellNods").nodes}
    te_nodes, te_elements = set(), set()
    thickness = {}
    for name in original.element_sets_names:
        elements = original.get_element_set(name).elements
        if name.endswith("_TE_FLAT"):
            for element in elements:
                te_elements.add(element.id)
                te_nodes.update(idx[n] for n in element.node_ids)
        if name in props and hasattr(props[name], "total_thickness"):
            for element in elements:
                thickness[element.id] = props[name].total_thickness

    adjacent = defaultdict(list)
    skin_normals = defaultdict(lambda: np.zeros(3))
    for element in original.elements:
        ids = [idx[n] for n in element.node_ids]

        def normal(coords):
            p = coords[ids]
            return (np.cross(p[2] - p[0], p[3] - p[1]) if len(ids) == 4
                    else np.cross(p[1] - p[0], p[2] - p[0]))

        if all(i in outer for i in ids) and element.id in thickness:
            unit_normal = normal(x0) / np.linalg.norm(normal(x0))
            for i in ids:
                adjacent[i].append(thickness[element.id])
                skin_normals[i] += unit_normal
        if element.id not in te_elements:
            assert np.dot(normal(x0), normal(x1)) > 0

    moved = sorted(set(adjacent) - te_nodes)
    unchanged = sorted(set(range(len(x0))) - set(moved))
    assert moved
    assert set(range(len(x0))) - outer  # The fixture includes interior/web nodes.
    np.testing.assert_array_equal(x1[unchanged], x0[unchanged])
    np.testing.assert_array_equal(x1[:, 2], x0[:, 2])
    shift = np.linalg.norm(x1[moved] - x0[moved], axis=1)
    expected = np.array([0.5 * np.mean(adjacent[i]) for i in moved])
    np.testing.assert_allclose(shift, expected, rtol=0, atol=1e-9)
    assert abs(shift.mean() - expected.mean()) < 1e-9
    assert shift.max() < max(thickness.values()) / 2 + 1e-9
    sections = defaultdict(list)
    for i in outer:
        sections[round(x0[i, 2], 6)].append(i)
    # A concave profile's inward normal need not point toward its centroid.
    # Determine the mesh winding globally, then check each local normal strictly.
    centers = {z: x0[ids, :2].mean(axis=0) for z, ids in sections.items()}
    votes = sum(
        np.sign(np.dot(n[:2], centers[round(x0[i, 2], 6)] - x0[i, :2]))
        for i, n in skin_normals.items()
    )
    inward_sign = np.sign(votes)
    assert inward_sign != 0
    for i in moved:
        assert np.dot(x1[i, :2] - x0[i, :2], inward_sign * skin_normals[i][:2]) > 0


def test_invalid_skin_reference():
    with pytest.raises(ValueError, match="skin_reference"):
        generator(skin_reference="invalid")
