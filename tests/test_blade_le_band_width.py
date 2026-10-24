"""The leading-edge reinforcement band keeps the width the yaml gives it.

Same class of defect as ``test_blade_te_band_width.py`` (blade input audit, finding G5,
``odd/tasks/blade-input-audit.md``). NuMAD capped keypoint ``a`` at 10 % of each side's arc, whatever
room there was before the spar cap, so where the side arc is under 4 m the 0.4 m per side of
``LE_reinforcement`` (width 0.8 m centred on the LE) was clipped and the strip got the a-b panel
laminate. The band may only yield to keep ``a`` 5 % of the side arc clear of the spar-cap edge ``b``.
"""

from pathlib import Path

import numpy as np
import pytest

from aeroelast.models.blade.numad.objects.blade import Blade

YAML = Path(__file__).parent / "IEA-15-240-RWT.yaml"
HP_B, HP_A, LE, LP_A, LP_B = 4, 5, 6, 7, 8  # key_arcs rows: te, e, d, c, b, a, le, a, b, c, d, e, te


@pytest.fixture(scope="module")
def blade():
    return Blade(str(YAML))


def _widths(blade):
    kp = blade.keypoints
    return kp.key_arcs[LE] - kp.key_arcs[HP_A], kp.key_arcs[LP_A] - kp.key_arcs[LE]


def test_le_band_matches_yaml_where_the_panel_has_room(blade):
    hp, lp = _widths(blade)
    kp, g = blade.keypoints, blade.geometry
    n = g.arclength.shape[0]
    ns, nf = 1, n - 2
    leband = 1e-3 * np.asarray(blade.definition.leband, dtype=float)
    # stations where the band stays 5 % of the side arc clear of b and above the 1 % floor
    room_hp = (leband - (kp.key_arcs[LE] - kp.key_arcs[HP_B])) < -0.05 * -g.arclength[ns]
    room_lp = (leband - (kp.key_arcs[LP_B] - kp.key_arcs[LE])) < -0.05 * g.arclength[nf]
    floor_hp = leband > 0.01 * -g.arclength[ns]
    floor_lp = leband > 0.01 * g.arclength[nf]
    sel_hp, sel_lp = room_hp & floor_hp, room_lp & floor_lp
    assert sel_hp.sum() > 30 and sel_lp.sum() > 30
    np.testing.assert_allclose(hp[sel_hp], leband[sel_hp], rtol=0, atol=1e-9)
    np.testing.assert_allclose(lp[sel_lp], leband[sel_lp], rtol=0, atol=1e-9)


def test_le_band_reaches_the_yaml_width_outboard(blade):
    """The clipped range (side arc < 4 m, r/R ~0.59-0.80) now carries the full 0.4 m per side."""
    hp, lp = _widths(blade)
    leband = 1e-3 * np.asarray(blade.definition.leband, dtype=float)
    span = np.asarray(blade.definition.ispan, dtype=float)
    rr = (span - span[0]) / (span[-1] - span[0])
    sel = (rr > 0.62) & (rr < 0.78) & (leband > 0.35)
    assert sel.sum() > 3
    np.testing.assert_allclose(hp[sel], leband[sel], rtol=0, atol=1e-9)
    np.testing.assert_allclose(lp[sel], leband[sel], rtol=0, atol=1e-9)


def test_band_stays_forward_of_the_spar_cap(blade):
    kp = blade.keypoints
    assert np.all(kp.key_arcs[HP_A] > kp.key_arcs[HP_B])
    assert np.all(kp.key_arcs[LP_A] < kp.key_arcs[LP_B])
    assert np.all(kp.key_arcs[HP_A] < kp.key_arcs[LE])
    assert np.all(kp.key_arcs[LP_A] > kp.key_arcs[LE])
