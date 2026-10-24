"""The trailing-edge reinforcement band keeps the width the yaml gives it.

Candidate A of the blade input audit (``odd/tasks/blade-input-audit.md``, T6). NuMAD capped keypoint
``d`` at 15 % of each side's arc (``keypoints.py``), whatever room the TE panel had. On the official
IEA 15 MW blade that clipped the 0.75 m TE band to 0.59 m at r/R 0.60 and dropped 23 % of the TE
glass there. The band may only yield to keep ``d`` aft of the spar-cap keypoint ``c``.

T15 item 4: the yaml band width is measured from the trailing edge, which for a flat-back section is
the corner (``arclength[ns]`` / ``arclength[nf]``) because the audit's ``nd_arc`` excludes the flat
face. ``d`` used to be measured from the face midpoint (``arclength[0]`` / ``[-1]``), so every band
was shortened by the half-face height (the official blade has a flat back up to r/R ~0.4).
"""

from pathlib import Path

import numpy as np
import pytest

from aeroelast.models.blade.numad.objects.blade import Blade

YAML = Path(__file__).parent / "IEA-15-240-RWT.yaml"
HP_D, HP_C, LP_C, LP_D = 2, 3, 9, 10  # key_arcs rows: te, e, d, c, b, a, le, a, b, c, d, e, te


@pytest.fixture(scope="module")
def blade():
    return Blade(str(YAML))


def _band_widths(blade):
    """Band width from the TE corner (arclength[ns] / arclength[nf]) per side, per station."""
    g, kp = blade.geometry, blade.keypoints
    ns, nf = 1, g.arclength.shape[0] - 2
    hp = kp.key_arcs[HP_D] - g.arclength[ns]
    lp = g.arclength[nf] - kp.key_arcs[LP_D]
    return hp, lp


def test_te_band_matches_yaml_where_the_panel_has_room(blade):
    hp, lp = _band_widths(blade)
    kp, g = blade.keypoints, blade.geometry
    n = g.arclength.shape[0]
    ns, nf = 1, n - 2
    teband = 1e-3 * np.asarray(blade.definition.teband, dtype=float)
    # stations where the band stays 5 % of the side arc clear of c and above the 2 %/4 % floors
    room_hp = (kp.key_arcs[HP_C] - (g.arclength[ns] + teband)) > 0.05 * -g.arclength[ns]
    room_lp = ((g.arclength[nf] - teband) - kp.key_arcs[LP_C]) > 0.05 * g.arclength[nf]
    floor_hp = teband > 0.02 * -g.arclength[ns]
    floor_lp = teband > 0.04 * g.arclength[nf]
    sel_hp, sel_lp = room_hp & floor_hp, room_lp & floor_lp
    assert sel_hp.sum() > 30 and sel_lp.sum() > 30
    np.testing.assert_allclose(hp[sel_hp], teband[sel_hp], rtol=0, atol=1e-9)
    np.testing.assert_allclose(lp[sel_lp], teband[sel_lp], rtol=0, atol=1e-9)


def test_band_stays_aft_of_the_spar_cap(blade):
    kp = blade.keypoints
    assert np.all(kp.key_arcs[HP_D] < kp.key_arcs[HP_C])
    assert np.all(kp.key_arcs[LP_D] > kp.key_arcs[LP_C])
