"""Analytic material proofs for the explicit S12 adapter, not a windIO-wide rule."""
from copy import deepcopy
from dataclasses import replace
from types import MappingProxyType, SimpleNamespace

import numpy as np
import pytest

from aeroelast.core.laminate import Laminate, Ply
from aeroelast.core.material import OrthotropicMaterial
from aeroelast.models.blade.numad.objects.windio2_laminate import LayerSample, LaminateCell
from tests.software.mesh.test_windio2_web_mesh import short as short
from tests.software.mesh.test_windio2_topology import official as official


ORTHO = MappingProxyType({'name': 'ud', 'E': (135e9, 9e9, 9e9),
                         'G': (5e9, 3e9, 5e9), 'nu': (.27, .3, .27), 'rho': 1600.})
CORE = MappingProxyType({'name': 'foam', 'E': 80e6, 'G': 30e6, 'nu': .3, 'rho': 90.})


def sample(name, thickness, angle, material=ORTHO):
    return LayerSample(name, material['name'], material, thickness, angle, None)


def cell(layers, surface='skin'):
    return LaminateCell(surface, (0., 1.), ((0., 0., 0.), (1., 0., 0.)), tuple(layers))


def mapping(*args):
    from aeroelast.models.blade.numad.mesh_gen.windio2_basis import map_stack
    return map_stack(*args)


def consumer(stack, centered=False):
    """Existing CLT consumer; explicitly translate its centered reference if needed."""
    plies = []
    for layer in stack.layers:
        m = layer.source.material
        if np.ndim(m['E']) == 0:
            e, nu = m['E'], m['nu']
            material = OrthotropicMaterial(m['name'], (e, e, e),
                                           (e/(2*(1+nu)),)*3, (nu,)*3, m['rho'])
        else:
            material = OrthotropicMaterial(m['name'], m['E'], m['G'], m['nu'], m['rho'])
        plies.append(Ply(material, layer.source.thickness_m, layer.angle_deg))
    laminate = Laminate(plies)
    offset = 0. if centered else stack.midplane_offset_m
    a, b, d = laminate.A, laminate.B, laminate.D
    return a, b + offset*a, d + 2*offset*b + offset**2*a


def energy_oracle(layers, basis, strain, curvature):
    """Independent tensor energy and exact two-point thickness integration.

    No production angle transform, Q/Qbar, ABD, or z-position helper is used.
    Fiber vectors and z bounds in `layers` are independently specified facts.
    """
    value = 0.
    for source, fiber, lo, hi in layers:
        m = source.material
        e1, e2 = (m['E'], m['E']) if np.ndim(m['E']) == 0 else m['E'][:2]
        nu = m['nu'] if np.ndim(m['nu']) == 0 else m['nu'][0]
        g = e1/(2*(1+nu)) if np.ndim(m['E']) == 0 else m['G'][0]
        compliance = np.array([[1/e1, -nu/e1, 0.], [-nu/e1, 1/e2, 0.], [0., 0., 1/g]])
        f = np.asarray(fiber)
        transverse = np.cross(basis[2], f)
        for xi in (-1/np.sqrt(3), 1/np.sqrt(3)):
            z = (lo+hi)/2 + xi*(hi-lo)/2
            engineering = strain + z*curvature
            tensor = np.array([[engineering[0], engineering[2]/2],
                               [engineering[2]/2, engineering[1]]])
            tensor3 = basis[:2].T @ tensor @ basis[:2]
            local = np.array([f@tensor3@f, transverse@tensor3@transverse,
                              2*f@tensor3@transverse])
            value += .5 * local @ np.linalg.solve(compliance, local) * (hi-lo)/2
    return value


def abd_energy(abd, strain, curvature):
    a, b, d = abd
    return .5*(strain@a@strain + 2*strain@b@curvature + curvature@d@curvature)


@pytest.mark.parametrize('normal,positive', [([0., 1., 0.], [1., 0., 0.]),
                                            ([0., -1., 0.], [-1., 0., 0.])])
def test_signed_skin_directions_and_orthotropic_response(normal, positive):
    # NuMAD: +x toward LE, +y suction, +z outboard. These analytic plane
    # directions are fixed independently of the element's arbitrary x-axis.
    n, transverse = np.array(normal), np.array(positive)
    beam = np.array([0., 0., 1.])
    basis = np.array([beam, transverse, n])
    responses = {}
    for angle in (0., 31., -31.):
        source = sample('outer', .002, angle)
        stack = mapping(cell([source]), beam, n, basis)
        layer, = stack.layers
        radians = np.deg2rad(angle)
        physical = np.cos(radians)*beam + np.sin(radians)*transverse
        np.testing.assert_allclose(layer.fiber_vector, physical, atol=2e-14)
        assert layer.angle_deg == pytest.approx(angle, abs=2e-13)
        assert layer.source is source and layer.source.material is ORTHO
        np.testing.assert_allclose(layer.z_bounds_m, [-.002, 0.], atol=1e-16)
        assert stack.reference == 'outer_mold_line'
        assert stack.midplane_offset_m == pytest.approx(-.001)
        abd = consumer(stack, centered=True)
        for strain in (np.array([.001, .0002, -.0003]), np.array([0., 0., .001])):
            analytic = energy_oracle([(source, physical, -.001, .001)], basis,
                                     strain, np.zeros(3))
            assert abd_energy(abd, strain, np.zeros(3)) == pytest.approx(analytic, rel=2e-13)
        responses[angle] = abd[0]
    assert responses[31.][0, 2] > 0
    np.testing.assert_allclose(responses[31.][:2, 2], -responses[-31.][:2, 2], rtol=1e-14)
    assert responses[0.][0, 2] == pytest.approx(0., abs=1e-9)


def test_asymmetric_reversal_same_center_and_oml_offset():
    beam, n = np.array([1., 0., 0.]), np.array([0., 0., 1.])
    basis = np.eye(3)
    layers = [sample('outer', .001, 0.), sample('inner', .003, 37.)]
    forward = mapping(cell(layers), beam, n, basis)
    reverse = mapping(cell(layers[::-1]), beam, n, basis)
    assert [l.source.name for l in forward.layers] == ['inner', 'outer']
    np.testing.assert_allclose([l.z_bounds_m for l in forward.layers],
                               [[-.004, -.001], [-.001, 0.]], atol=1e-16)
    f, r = consumer(forward, True), consumer(reverse, True)
    assert np.linalg.norm(f[1]) > 1e4  # asymmetric, not an all-zero B proof
    for index, sign in enumerate((1., -1., 1.)):
        np.testing.assert_allclose(r[index], sign*f[index], rtol=3e-14, atol=1e-8)
    strain, curvature = np.array([.001, -.0003, .0004]), np.array([.12, -.08, .15])
    physical = [(layers[1], [np.cos(np.deg2rad(37)), np.sin(np.deg2rad(37)), 0.], -.004, -.001),
                (layers[0], [1., 0., 0.], -.001, 0.)]
    expected = energy_oracle(physical, basis, strain, curvature)
    assert abd_energy(consumer(forward), strain, curvature) == pytest.approx(expected, rel=3e-13)
    assert not np.isclose(abd_energy(f, strain, curvature), expected, rtol=.01)


def test_projection_and_opposite_bases_preserve_physical_stack_and_energy():
    # Analytic sloped shell with oblique beam: project, do not use chart x.
    n = np.array([0., -.6, .8])
    zero = np.array([0., .8, .6])
    beam = np.array([0., .8, 1.6])  # zero + 1.25*n, up to a scalar
    transverse = np.array([-1., 0., 0.])
    phi = np.deg2rad(23.)
    x = np.cos(phi)*zero + np.sin(phi)*transverse
    y = -np.sin(phi)*zero + np.cos(phi)*transverse
    basis = np.array([x, y, n])
    opposite = np.array([x, -y, -n])
    sources = [sample('outside', .001, 17.), sample('inside', .002, -29.)]
    direct = mapping(cell(sources), beam, n, basis)
    flipped = mapping(cell(sources), beam, n, opposite)
    np.testing.assert_allclose(basis@basis.T, np.eye(3), atol=1e-14)
    assert np.linalg.det(opposite) == pytest.approx(1.)
    assert [l.source.name for l in flipped.layers] == ['outside', 'inside']
    for stack, frame in ((direct, basis), (flipped, opposite)):
        for l in stack.layers:
            angle = np.deg2rad(l.source.fiber_angle_deg)
            physical = np.cos(angle)*zero + np.sin(angle)*transverse
            np.testing.assert_allclose(l.fiber_vector, physical, atol=2e-14)
            recovered = np.cos(np.deg2rad(l.angle_deg))*frame[0] + np.sin(np.deg2rad(l.angle_deg))*frame[1]
            np.testing.assert_allclose(recovered, physical, atol=2e-14)
    for a, b in zip(direct.layers, reversed(flipped.layers), strict=True):
        np.testing.assert_allclose(np.outer(a.z_bounds_m, basis[2]),
                                   np.outer(-np.array(b.z_bounds_m)[::-1], basis[2]), atol=1e-16)
    strain, curvature = np.array([.001, -.0002, .0003]), np.array([.1, -.2, .3])
    # y and normal both reverse: engineering shear strain reverses; normal
    # reversal also reverses curvature's normal components, not its shear.
    other_strain = strain*np.array([1., 1., -1.])
    other_curvature = curvature*np.array([-1., -1., 1.])
    expected = energy_oracle([(sources[1], np.cos(np.deg2rad(-29))*zero + np.sin(np.deg2rad(-29))*transverse, -.003, -.001),
                              (sources[0], np.cos(np.deg2rad(17))*zero + np.sin(np.deg2rad(17))*transverse, -.001, 0.)],
                             basis, strain, curvature)
    assert abd_energy(consumer(direct), strain, curvature) == pytest.approx(expected, rel=3e-13)
    assert abd_energy(consumer(flipped), other_strain, other_curvature) == pytest.approx(expected, rel=3e-13)


def test_web_triple_faces_signed_directions_and_reference():
    beam, le = np.array([0., 0., 1.]), np.array([1., 0., 0.])
    basis = np.array([beam, [0., -1., 0.], le])
    sources = [sample('LEskin', .001, 31.), sample('core', .01, 0., CORE),
               sample('TEskin', .003, 31.)]
    stack = mapping(cell(sources, 'web'), beam, le, basis)
    assert stack.reference == 'web_core_center'
    assert [l.source.name for l in stack.layers] == ['TEskin', 'core', 'LEskin']
    np.testing.assert_allclose([l.z_bounds_m for l in stack.layers],
                               [[-.008, -.005], [-.005, .005], [.005, .006]], atol=1e-16)
    assert stack.midplane_offset_m == pytest.approx(-.001)
    assert [l.angle_deg for l in stack.layers] == pytest.approx([-31., 0., 31.])
    assert stack.layers[0].fiber_vector[1] > 0  # TE face: toward suction
    assert stack.layers[2].fiber_vector[1] < 0  # LE face: toward pressure
    strain, curvature = np.array([.001, .0002, -.0003]), np.array([.1, -.2, .3])
    angle = np.deg2rad(31.)
    physical = [(sources[2], [0., np.sin(angle), np.cos(angle)], -.008, -.005),
                (sources[1], beam, -.005, .005),
                (sources[0], [0., -np.sin(angle), np.cos(angle)], .005, .006)]
    assert abd_energy(consumer(stack), strain, curvature) == pytest.approx(
        energy_oracle(physical, basis, strain, curvature), rel=3e-13)
    flipped = mapping(cell(sources, 'web'), beam, le, np.array([basis[0], -basis[1], -basis[2]]))
    assert [l.source.name for l in flipped.layers] == ['LEskin', 'core', 'TEskin']
    assert abd_energy(consumer(flipped), strain*[1, 1, -1], curvature*[-1, -1, 1]) == pytest.approx(
        abd_energy(consumer(stack), strain, curvature), rel=3e-13)


@pytest.mark.parametrize('layers', [[], [sample('a', .001, 0.)],
    [sample('a', .001, 0.), sample('b', .001, 0.)],
    [sample('a', .001, 0.), sample('b', .001, 0.), sample('c', .001, 0.)],
    [sample('a', .001, 0.), sample('c', .01, 0., CORE), sample('b', .001, 0.), sample('extra', .001, 0.)]])
def test_unsupported_web_stacks_rejected_without_source_mutation(layers):
    original = cell(layers, 'web')
    with pytest.raises(ValueError, match='requires one orthotropic LE skin, isotropic core, orthotropic TE skin'):
        mapping(original, [1., 0., 0.], [0., 0., 1.], np.eye(3))
    assert original.layers == tuple(layers)


def test_degenerate_projection_rejected_and_empty_skin_preserved():
    with pytest.raises(ValueError, match='beam-axis projection'):
        mapping(cell([sample('a', .001, 0.)]), [0., 0., 1.], [0., 0., 1.], np.eye(3))
    assert mapping(cell([]), [1., 0., 0.], [0., 0., 1.], np.eye(3)).layers == ()


def test_raw_native_oblique_beam_projection_not_chart_tangent(short):
    from aeroelast.models.blade.numad.mesh_gen.windio2_mesh import RawMesh, map_mesh_materials
    definition = deepcopy(short[0]._definition)
    # Known rigid rotation, independent of the production section-frame recipe.
    a, b = np.deg2rad([19., 27.])
    rotation = np.array([[np.cos(b), -np.sin(b), 0.], [np.sin(b), np.cos(b), 0.], [0., 0., 1.]]) @ np.array([
        [1., 0., 0.], [0., np.cos(a), -np.sin(a)], [0., np.sin(a), np.cos(a)]])
    beam = rotation[:, 2]
    normal = rotation @ np.array([0., 1., .25])
    normal /= np.linalg.norm(normal)
    x = rotation[:, 0]
    y = np.cross(normal, x)
    basis = np.array([x, y, normal])
    nodes = np.array([-x-y, x-y, x+y, -x+y])*.05
    uv = np.array([[[.2, .3], [.3, .3], [.3, .4], [.2, .4]]])
    # The native source axis has the exact analytic beam above, expressed in
    # windIO axes, with nonzero translation and twist. Neither changes its zero.
    windio_slope = np.array([beam[1], -beam[0], beam[2]])
    definition.windio2_reference_axis = {axis: {'grid': [0., 1.], 'values': [i, i+v]}
                                        for i, (axis, v) in enumerate(zip('xyz', windio_slope, strict=True))}
    definition.windio2_outer_shape['twist'] = {'grid': [0., 1.], 'values': [23., 23.]}
    source = sample('outer', .002, 31.)
    raw = RawMesh(nodes, np.array([[0, 1, 2, 3]]), (cell([source]),), ('skin',), uv, basis[None], {})
    mapped = map_mesh_materials(SimpleNamespace(_definition=definition), raw)
    layer, = mapped.material_stacks[0].layers
    zero = beam - (beam@normal)*normal
    zero /= np.linalg.norm(zero)
    angle = np.deg2rad(31.)
    expected = np.cos(angle)*zero + np.sin(angle)*np.cross(normal, zero)
    np.testing.assert_allclose(layer.fiber_vector, expected, atol=3e-14)
    assert abs(expected@x) > .1  # nonzero coupling in this arbitrary chart basis
    assert abs(zero@x) < 1e-14  # chart x is NOT the beam-axis zero
    np.testing.assert_array_equal(mapped.nodes, nodes)
    np.testing.assert_allclose(layer.z_bounds_m, [-.002, 0.], atol=1e-16)


def test_official_raw_triple_mapping_keeps_reference_nodes(official):
    from aeroelast.models.blade.numad.mesh_gen.windio2_mesh import generate_mesh
    blade, _, _ = official
    eta = .08163265306122448
    geometry = generate_mesh(blade.windio2_topology, eta_bounds=(eta, eta+.0001))
    mapped = generate_mesh(blade.windio2_topology, eta_bounds=(eta, eta+.0001), with_materials=True)
    np.testing.assert_array_equal(mapped.nodes, geometry.nodes)
    np.testing.assert_array_equal(mapped.elements, geometry.elements)
    assert geometry.material_angles is None
    assert len(mapped.material_stacks) == len(mapped.elements)
    for source_cell, stack, basis in zip(mapped.cells, mapped.material_stacks, mapped.bases, strict=True):
        assert {l.source.name for l in stack.layers} == {l.name for l in source_cell.layers}
        if source_cell.surface != 'skin':
            assert stack.reference == 'web_core_center'
            assert len(stack.layers) == 3
            # Source local LE direction has +x at zero twist; use the stored
            # physical normal only to establish the signed element z-axis.
            for layer in stack.layers:
                lo, hi = layer.z_bounds_m
                face_position = (lo+hi)/2 * np.dot(stack.physical_normal, basis[2])
                if layer.source.name.endswith('skinLE'):
                    assert face_position > 0
                elif layer.source.name.endswith('skinTE'):
                    assert face_position < 0
                assert layer.source.fiber_angle_deg == 0.


def test_raw_mapping_opt_in_preserves_geometry_and_source(short):
    from aeroelast.models.blade.numad.mesh_gen.windio2_mesh import generate_mesh, map_mesh_materials
    topology, _, path, before = short
    geometry = generate_mesh(topology)
    assert geometry.material_angles is None
    with pytest.raises(ValueError, match='requires one orthotropic LE skin, isotropic core, orthotropic TE skin'):
        map_mesh_materials(topology, geometry)
    assert geometry.material_angles is None and path.read_bytes() == before
    # This fixture intentionally has unsupported webs; use its unchanged source
    # geometry with skin cells only to verify the production raw entry point.
    select = np.array(geometry.surfaces) == 'skin'
    skin = replace(geometry, elements=geometry.elements[select],
                   cells=tuple(c for c, yes in zip(geometry.cells, select, strict=True) if yes),
                   surfaces=tuple('skin' for _ in range(select.sum())),
                   element_parameters=geometry.element_parameters[select], bases=geometry.bases[select])
    mapped = map_mesh_materials(topology, skin)
    np.testing.assert_array_equal(mapped.nodes, skin.nodes)
    np.testing.assert_array_equal(mapped.elements, skin.elements)
    assert mapped.cells is skin.cells
    assert len(mapped.material_angles) == len(skin.elements)
    assert all(s.reference == 'outer_mold_line' for s in mapped.material_stacks)
    for stack, source_cell in zip(mapped.material_stacks, mapped.cells, strict=True):
        assert [l.source for l in stack.layers] == list(reversed(source_cell.layers))
        assert stack.midplane_offset_m == pytest.approx(-source_cell.thickness_m/2)
    # Changing only the right-handed element basis must not change physical
    # face positions or fibers in the production raw mapper either.
    flipped_bases = skin.bases * np.array([[1.], [-1.], [-1.]])
    flipped = map_mesh_materials(topology, replace(skin, bases=flipped_bases))
    for direct, other in zip(mapped.material_stacks, flipped.material_stacks, strict=True):
        assert [l.source for l in other.layers] == [l.source for l in reversed(direct.layers)]
        for a, b in zip(direct.layers, reversed(other.layers), strict=True):
            np.testing.assert_allclose(a.fiber_vector, b.fiber_vector, atol=2e-14)
            np.testing.assert_allclose(a.z_bounds_m, -np.array(b.z_bounds_m)[::-1], atol=1e-16)
    assert path.read_bytes() == before
