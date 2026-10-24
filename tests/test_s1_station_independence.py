"""S-1 extraction does not depend on the number of stations the mesher inserts (blade-input-audit T7d).

``BladeMesh.max_bay_length = 0.25`` turns the ~50 yaml-grid stations into ~490, and the tangent-frame fit
window (twice the median spacing) shrank from 4.8 m to 0.48 m, the spline smoothing bounds act on a sum over
more points, so the S-1 numbers moved with the mesher setting, not with the structure.
``ei_from_static_response(stations=...)`` fits on the mesh stations nearest to a fixed reference grid.

Synthetic cantilever with tapered EI(z): tip force F, curvature F (L - z) / EI(z), rigid section rotations.
"""

import numpy as np

from aeroelast.postprocess.sectional import SectionalExtractor, select_stations

L = 10.0
F = 1.0e3
HALF_W = (0.3, 0.2)  # ring half extents in x, y


def ei_true(z):
    return 5.0e6 * (1.0 - 0.6 * z / L)


def _extractor(zs):
    """Extractor shell with rings of 4 nodes per station and the exact rigid-section displacement."""
    zs = np.asarray(zs, dtype=float)
    ring = np.array([[1, 1], [-1, 1], [-1, -1], [1, -1]]) * np.array(HALF_W)
    coords = np.vstack([np.column_stack([ring, np.full(4, z)]) for z in zs])
    nodes = [np.arange(4 * i, 4 * i + 4) for i in range(zs.size)]
    zf = np.linspace(0.0, L, 20001)
    kap = F * (9.5 - zf) / ei_true(zf)  # f'' (the tip force acts at the last station, z = 9.5)
    slope = np.concatenate([[0.0], np.cumsum(0.5 * (kap[1:] + kap[:-1]) * np.diff(zf))])  # f'
    defl = np.concatenate([[0.0], np.cumsum(0.5 * (slope[1:] + slope[:-1]) * np.diff(zf))])  # f
    sol = {}
    for name in ("flap", "edge"):
        u = np.zeros((coords.shape[0], 6))
        for i, z in enumerate(zs):
            s, f = np.interp(z, zf, slope), np.interp(z, zf, defl)
            d = coords[nodes[i]] - [0.0, 0.0, z]
            if name == "flap":  # tip force in +Y: u_y = f, theta_x = -f'
                u[nodes[i], 1] = f
                u[nodes[i], 2] = -s * d[:, 1]
            else:  # tip force in +X: u_x = f, theta_y = +f'
                u[nodes[i], 0] = f
                u[nodes[i], 2] = -s * d[:, 0]
        sol[name] = u.ravel()
    ext = SectionalExtractor.__new__(SectionalExtractor)
    ext.coords, ext.dpn = coords, 6
    return ext, nodes, sol


def _ei(zs, stations=None):
    ext, nodes, sol = _extractor(zs)
    return ext._ei_tangent_frame(zs, nodes, sol, F, stations)


COARSE = np.linspace(0.5, 9.5, 10)
DENSE = np.linspace(0.5, 9.5, 91)  # 9x the stations, contains COARSE


def test_dense_mesh_with_reference_stations_equals_coarse_mesh():
    zc, fc, ec = _ei(COARSE)
    zd, fd, ed = _ei(DENSE, stations=COARSE)
    np.testing.assert_allclose(zd, zc, rtol=0, atol=1e-12)
    np.testing.assert_allclose(fd, fc, rtol=1e-9, atol=0)
    np.testing.assert_allclose(ed, ec, rtol=1e-9, atol=0)


def test_without_reference_stations_the_density_matters():
    """The defect: all stations fitted gives a different value at the coarse stations."""
    _, fc, _ = _ei(COARSE)
    zd, fd, _ = _ei(DENSE)
    at = np.isin(np.round(zd, 9), np.round(COARSE, 9))
    assert np.max(np.abs(fd[at][1:-1] / fc[1:-1] - 1)) > 1e-4


def test_coarse_set_alone_reproduces_current_method():
    """Identity on a mesh whose stations are the reference grid (the 9f72c6b extraction)."""
    base = _ei(COARSE)
    sel = _ei(COARSE, stations=COARSE)
    for a, b in zip(base, sel, strict=True):
        np.testing.assert_allclose(a, b, rtol=0, atol=1e-12)
    # sanity of the synthetic field: the fit tracks the analytic EI (a 2 m window is coarse, hence 6 %)
    np.testing.assert_allclose(base[1][1:-1], ei_true(COARSE)[1:-1], rtol=6e-2)


def test_select_stations_nearest_unique_sorted():
    zs = np.array([0.0, 0.1, 0.2, 1.0, 1.1, 2.0])
    np.testing.assert_array_equal(select_stations(zs, [0.12, 1.04, 1.9, 2.5]), [1, 3, 5])
    np.testing.assert_array_equal(select_stations(zs, [0.1, 0.1]), [1])
    np.testing.assert_array_equal(select_stations(zs, None), np.arange(zs.size))
    np.testing.assert_array_equal(select_stations(zs, zs), np.arange(zs.size))
