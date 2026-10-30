"""Characterize the public coordinate boundary, not an extended-loft solution.

Independent affine surface: F(eta, q) = (q, 0, eta). Its Jacobian is
constant and positive between the two exact parsed-binary64 crossing roots.
No production interpolation, topology predicates or meshing oracle is reused.
"""
from decimal import Decimal, localcontext
from fractions import Fraction

import numpy as np
import pytest

from aeroelast.core.mesh.entities import Node
from aeroelast.models.blade.numad.mesh_gen.windio2_basis import shell_bases
from tests.software.mesh.test_windio2_topology import official  # noqa: F401


def test_le_landmark_brackets_interior_minimum_not_te_extrema():
    from aeroelast.models.blade.numad.objects.windio2_precision import minimum_location

    # A single interior minimum at .4, with near-TE maxima at .02 and .98.
    # Endpoint slope signs are the reverse of the minimum's local bracket.
    def derivative(s):
        return -(s-np.longdouble('.02'))*(s-np.longdouble('.4'))*(s-np.longdouble('.98'))

    location = minimum_location(derivative)
    assert abs(location-np.longdouble('.4')) <= 2*np.finfo(np.longdouble).eps
    with pytest.raises(ValueError, match='leading-edge landmark is not bracketed'):
        minimum_location(lambda s: np.ones_like(s))


@pytest.mark.parametrize('translation', [0, 1, 100])
def test_exact_finite_band_collapses_at_public_node_boundary(translation):
    # Crossings eta=.3 and .7*eta=.21, not decimal-exact reconstruction.
    a = Fraction(.3)
    b = Fraction(.21)/Fraction(.7)
    lo, hi = sorted((a, b))
    assert lo < hi and float(lo) == float(hi)
    # Connectivity is positive in (eta,q), hence toward +y for this surface.
    exact = [(Fraction(0), Fraction(0), lo+translation),
             (Fraction(0), Fraction(0), hi+translation),
             (Fraction(1), Fraction(0), hi+translation),
             (Fraction(1), Fraction(0), lo+translation)]
    assert (hi-lo)/4 > 0  # exact physical bilinear Jacobian at ALL (r,s)
    # Decimal is only an independent conversion oracle, not a selected runtime
    # precision policy. NumPy longdouble is demonstrably sufficient here.
    with localcontext() as context:
        context.prec = 80
        extended = np.array([[np.longdouble(str(Decimal(v.numerator)/Decimal(v.denominator)))
                              for v in point] for point in exact])
    assert np.unique(extended, axis=0).shape == (4, 3)
    local = extended-extended[0]
    for r in (-1, 0, 1):
        for s in (-1, 0, 1):
            dr = np.array([-(1-s), 1-s, 1+s, -(1+s)], dtype=np.longdouble)/4
            ds = np.array([-(1-r), -(1+r), 1+r, 1-r], dtype=np.longdouble)/4
            assert np.cross(dr@local, ds@local)[1] > 0
    # The actual public entity constructor, not a simulation of its cast.
    public = np.array([Node(point).coords for point in extended])
    assert public.dtype == np.dtype('float64')
    assert np.unique(public, axis=0).shape == (2, 3)
    np.testing.assert_array_equal(public[0], public[1])
    np.testing.assert_array_equal(public[2], public[3])
    np.testing.assert_array_equal(np.cross(public[1]-public[0], public[3]-public[0]), 0)
    with pytest.raises(ValueError, match='Nonpositive shell quadrilateral on source geometry'):
        shell_bases(public, np.array([[0, 1, 2, 3]]))


@pytest.mark.parametrize('values', [[0., 1.], [0., 1., 1.5, 4.], [0., 1., -1., 2.]])
def test_extended_pchip_value_derivative_and_loft_chain_rule(values):
    from scipy.interpolate import PchipInterpolator
    from aeroelast.models.blade.numad.objects.windio2_precision import Pchip
    grid = np.linspace(0., 1., len(values))
    query = np.linspace(0., 1., 23)
    reference = PchipInterpolator(grid, values)
    curve = Pchip(grid, values)
    np.testing.assert_allclose(curve(query), reference(query), atol=2e-14, rtol=0)
    np.testing.assert_allclose(curve.derivative()(query), reference.derivative()(query), atol=2e-14, rtol=0)
    # Positive affine ordinate transformations commute with PCHIP. This gives
    # an analytic directional derivative without differencing any small step.
    tangent = 2*np.array(values)+3
    np.testing.assert_allclose(curve.directional(query, tangent),
                               2*reference(query)+3, atol=2e-14, rtol=0)


def test_node_precision_is_explicit_and_default_is_unchanged():
    value = np.longdouble('0.30000000000000001')
    ordinary = Node([value, 2])
    assert ordinary.coords.dtype == np.dtype(float)
    np.testing.assert_array_equal(ordinary.coords, [float(value), 2., 0.])
    extended = Node([value, 2], coordinate_dtype=np.longdouble)
    assert extended.coords.dtype == np.dtype(np.longdouble)
    assert extended.x == value and extended.x != ordinary.x
    assert extended.y == 2 and extended.z == 0
    assert extended.id == ordinary.id+1
    with pytest.raises(ValueError, match='coordinate_dtype must be float64 or longdouble'):
        Node([1, 2, 3], coordinate_dtype=np.int64)
    assert Node([0, 0, 0]).id == extended.id+1


@pytest.mark.parametrize('case', ['arc_precision', 'span_precision', 'same_float_roots'])
def test_extended_public_mesh_keeps_finite_bands(tmp_path, official, case):
    import copy
    import hashlib
    import yaml
    from collections import Counter
    from aeroelast.core.mesh.windio2 import geometry_mesh
    from aeroelast.models.blade.numad.objects.windio2_topology import WindIO2Topology
    from tests.software.mesh.test_windio2_mesh import synthetic_chart

    seed, source, path = synthetic_chart(tmp_path, official, case)
    definition = copy.deepcopy(seed._definition)
    # Constant section, straight axis: an extruded source contour, not a
    # chart-only oracle. Keep all material/source records and explicit cuts.
    definition.windio2_outer_shape['rthick']['values'] = [
        .5]*len(definition.windio2_outer_shape['rthick']['grid'])
    if case == 'same_float_roots':
        source['components']['blade']['structure']['anchors'] = [
            {'name': name, 'start_nd_arc': {'grid': [0., 1.], 'values': values}}
            for name, values in [('diagonal', [0., 1.]), ('level', [.3, .3]),
                                 ('scaled', [0., .7]), ('scaled_level', [.21, .21])]]
    path.write_text(yaml.safe_dump(source))
    before = path.read_bytes()
    topology = WindIO2Topology(source, hashlib.sha256(before).hexdigest(), definition)
    mesh = geometry_mesh(topology, .25, (0., 1.))
    raw = mesh.windio2_metadata.raw_mesh
    assert mesh.coords_array.dtype == np.dtype(np.longdouble)
    assert np.isfinite(mesh.coords_array).all()
    np.testing.assert_array_equal(mesh.coords_array, raw.nodes)
    xyz = raw.nodes[raw.elements]
    local = xyz-xyz[:, :1]
    normals = raw.bases[:, 2]
    for r in (-1, 0, 1):
        for s in (-1, 0, 1):
            dr = np.array([-(1-s), 1-s, 1+s, -(1+s)], dtype=np.longdouble)/4
            ds = np.array([-(1-r), -(1+r), 1+r, 1-r], dtype=np.longdouble)/4
            jac = np.sum(np.cross(dr@local, ds@local)*normals, axis=1)
            assert np.all(jac > 0)
    edges = Counter(tuple(sorted((a, b))) for quad in raw.elements
                    for a, b in zip(quad, np.roll(quad, -1), strict=True))
    assert set(edges.values()) <= {1, 2}
    for (a, b), count in edges.items():
        if count == 1:
            assert raw.nodes[a, 2] == raw.nodes[b, 2]
            assert raw.nodes[a, 2] in (0, .5)
    uv = raw.element_parameters-raw.element_parameters[:, :1]
    area = np.sum(uv[:, :, 0]*np.roll(uv[:, :, 1], -1, axis=1)
                  -uv[:, :, 1]*np.roll(uv[:, :, 0], -1, axis=1), axis=1)/2
    assert np.all(area > 0)
    np.testing.assert_allclose(area.sum(), 1, atol=32*np.finfo(np.longdouble).eps, rtol=0)
    probe = [any(layer.name == 'probe' for layer in cell.layers) for cell in raw.cells]
    width = Fraction(.3+2**-48)-Fraction(.3) if case == 'arc_precision' else Fraction(.7)-Fraction(.3)
    np.testing.assert_allclose(area[probe].sum(), np.longdouble(width.numerator)/width.denominator,
                               atol=32*np.finfo(np.longdouble).eps, rtol=0)
    assert all(cell.layers[0].name == source['components']['blade']['structure']['layers'][0]['name']
               for cell in raw.cells)
    if case == 'same_float_roots':
        roots = [Fraction(.3), Fraction(.21)/Fraction(.7)]
        extended = topology.with_extended_precision()
        for root, names in zip(roots, [('diagonal', 'level'), ('scaled', 'scaled_level')], strict=True):
            point = next(p for p in extended.section(root).points
                         if ('anchor', names[0], 'start_nd_arc') in p.identities)
            assert ('anchor', names[1], 'start_nd_arc') in point.identities
        assert all(np.any(raw.element_parameters[:, :, 0] == np.longdouble(str(root.numerator))/np.longdouble(str(root.denominator)))
                   for root in roots)
    assert path.read_bytes() == before
    assert mesh.windio2_metadata.source_sha256 == hashlib.sha256(before).hexdigest()
    assert mesh.windio2_metadata.structural_ready is False
