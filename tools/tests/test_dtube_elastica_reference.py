import pytest

from tests.validation.benchmarks import test_dtube_camarena_bending as dtube
from tools.dtube_elastica_reference import brazier_secant_factor, tip_deflection


def test_brazier_factor_tends_to_one_for_small_curvature():
    assert brazier_secant_factor(0.0, dtube.R, dtube.T, dtube.NU) == 1.0
    assert brazier_secant_factor(1e-7, dtube.R, dtube.T, dtube.NU) == pytest.approx(1.0, abs=1e-6)


def test_brazier_factor_softens_monotonically():
    k = [1e-3, 3e-3, 6e-3]
    f = [brazier_secant_factor(x, dtube.R, dtube.T, dtube.NU) for x in k]
    assert 1.0 > f[0] > f[1] > f[2] > 0.0


def test_small_load_recovers_the_linear_cantilever():
    scale = 0.01
    assert tip_deflection(scale, "rigid") == pytest.approx(scale * dtube.DELTA_REF, rel=1e-3)
    assert tip_deflection(scale, "brazier") == pytest.approx(scale * dtube.DELTA_REF, rel=2e-3)
