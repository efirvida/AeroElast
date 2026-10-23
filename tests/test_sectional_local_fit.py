"""The tangent-frame local fit survives non-uniform stations (blade-input-audit T7a).

With ``BladeMesh.max_bay_length`` the stations are no longer uniform, and a window of twice the median
spacing left one station alone: ``polyfit`` fell to degree 0 and S-1 died with an IndexError.
"""

import numpy as np

from aeroelast.postprocess.sectional import _local_quadratic


def test_isolated_station_gets_three_neighbours():
    zs = np.array([0.0, 0.1, 0.2, 0.3, 0.4, 2.0, 2.1, 2.2])
    y = 3.0 * zs**2 + 2.0 * zs + 1.0
    half = 2.0 * float(np.median(np.diff(zs)))
    np.testing.assert_allclose(_local_quadratic(zs, y, half, True), 6.0 * zs + 2.0, atol=1e-9)
    np.testing.assert_allclose(_local_quadratic(zs, y, half, False), y, atol=1e-9)


def test_uniform_stations_keep_the_median_window():
    zs = np.linspace(0.0, 10.0, 21)
    y = np.sin(zs)
    half = 2.0 * float(np.median(np.diff(zs)))
    expected = []
    for z0 in zs:
        m = np.abs(zs - z0) <= half + 1e-9
        expected.append(np.polyfit(zs[m] - z0, y[m], 2)[-2])
    np.testing.assert_allclose(_local_quadratic(zs, y, half, True), expected, rtol=0, atol=1e-12)
