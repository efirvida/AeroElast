"""Analytic material proofs for the explicit S12 adapter, not a windIO-wide rule."""
from types import MappingProxyType

import numpy as np
import pytest

from aeroelast.core.laminate import Laminate, Ply
from aeroelast.core.material import OrthotropicMaterial
from aeroelast.models.blade.numad.objects.windio2_laminate import LayerSample, LaminateCell
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
