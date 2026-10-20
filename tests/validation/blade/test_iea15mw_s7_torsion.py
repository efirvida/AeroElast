"""S-7 — global torsional response: shell tip torque vs BeamDyn GJ tables.

The shell's global torsional response (tip torque via force couples, section
rotation from the displacement field) is compared against the BeamDyn 6×6 tables:
the station-by-station twist ratio in the 0.3-0.9 span band measures **0.8670** at
an element size of 0.25, trending to about **0.861** at 0.125 — the shell twists
about 13-14 % LESS than the beam, i.e. its global GJ runs about 15 % ABOVE the deck,
on the same side as the sectional S-1 comparison.

That direction has been corrected twice.  Until 2026-10-10 the tip couple
under-delivered the applied moment by a mesh-dependent factor κ (0.8978 at 0.25);
fixing it moved the ratio to 1.0466.  Until 2026-10-11 the beam reference integrated
1/GJ only up to the last BeamDyn station below each section and dropped the partial
segment, under-stating the beam twist by 16 % on average over the band; integrating to
``z`` moves the ratio to 0.8670.  Neither correction touches the shell model, and the
second one dissolves the sectional-versus-global sign disagreement.  The size of the
excess is not settled: about 2-3 % of it depends on how BeamDyn interpolates GJ between
key points (unverified), and most of the rest tracks the skin's reference surface — the
nodes sit on the outer mould line with the laminate centred on them, and moving them to
the laminate mid-surface closes the ratio to about 0.97.  See
``odd/tasks/s7-torsion-ratio.md`` sections E7 and E8.

The bound is the L3 ceiling, ±5 % around exact agreement, so the test is red while the
ratio reads about 0.91 at es 0.25.  The former ±30 % window had no quantified source and
could detect neither the 10.2 % load error nor the 16 % reference error, which is why the
couple construction and the reference integral each carry their own exact test.
"""

from __future__ import annotations

import sys

import numpy as np
import pytest

pytest.importorskip("petsc4py", reason="PETSc not available")

from tests.support.assertions import assert_relative_error  # noqa: E402
from tests.support.paths import REPO_ROOT  # noqa: E402

sys.path.insert(0, str(REPO_ROOT / "tools"))

from run_s7_torsion import (  # noqa: E402
    beam_twist_profile_deg,
    build_mesh_and_properties,
    run_torsion_case,
    section_twists_deg,
    tip_couple,
)

MOMENT_NM = 1000.0
# Span resolution.  The ratios below are the *measured* shell/beam twist ratio with the tip couple
# delivering exactly MOMENT_NM (see tip_couple) and the beam twist integrated to each section's z
# (see beam_twist_profile_deg).  es=2.0/1.0/0.5/0.25 were measured on 2026-10-11 with
# tools/diagnose_s7_beam_reference.py; es=0.125 is the 2026-10-10 load-corrected 1.039011 times
# the reference correction measured at 0.25 (0.8670/1.0466), so it is derived and good to about
# +-0.002, not measured.
#
#   es=2.000   1460 nodes  ratio=0.9251
#   es=1.000   3040 nodes  ratio=0.9227
#   es=0.500   9271 nodes  ratio=0.8874
#   es=0.250  32325 nodes  ratio=0.8670   <- the fixture
#   es=0.125 120352 nodes  ratio~0.861    <- derived (230 s, 12 GB peak RSS to measure)
#
# 0.250 sits about 0.7% above 0.125, which is what justifies it as a regression guard.  This is a
# slow test (about 1 min of PETSc per run), which is why it is marked slow.
#
# History: with the staircase reference these read 1.118319 / 1.116392 / 1.069529 / 1.046588 /
# 1.039011 (2026-10-10), and before the couple fix 1.0237 / 1.0219 / 0.9790 / 0.939603 / 0.933641.
# A table committed between 2026-09-30 and 2026-10-04 read 1.386 .. 1.273 and does not reproduce
# at HEAD at any mesh size.  See `odd/tasks/s7-torsion-ratio.md`.
ELEMENT_SIZE = 0.25

#: L3 ceiling (maintainer rule 2026-10-16): a relative bound of 0.05 on the ratio, i.e. the
#: band [0.95, 1.05] around exact agreement with the beam. It was 0.30 ([0.7, 1.3]), a band
#: with no quantified source; tightened in B3-2 batch 1 (odd/tasks/validation-l3-b3-audit.md),
#: red allowed.
_RATIO_TOL = 0.05


@pytest.mark.slow
def test_s7_global_torsion_ratio_within_band():
    mesh, properties = build_mesh_and_properties(ELEMENT_SIZE)
    u = run_torsion_case(mesh, properties, MOMENT_NM)
    u_arr = np.asarray(u, dtype=np.float64)

    zs, theta_s = section_twists_deg(mesh, u_arr)
    theta_beam = beam_twist_profile_deg(MOMENT_NM, zs)

    mid = (zs >= 0.3 * 117.0) & (zs <= 0.9 * 117.0)
    ratios = theta_s[mid] / theta_beam[mid]
    ratio_mean = float(np.nanmean(ratios))

    # The shell twists less than the beam: at the correct load and against the integral-to-z
    # reference the global ratio is 0.8670, i.e. the shell's global GJ is about 15 % above the
    # deck -- the same side as S-1's sectional comparison (+10.8 % signed at HEAD).  Most of the
    # excess tracks the skin's OML reference surface (E8); this test reports the measurement and
    # does not adjudicate the model difference.
    #
    # The bound is the L3 ceiling, [0.95, 1.05] around exact agreement, which is what
    # assert_relative_error expresses symmetrically (it was [0.7, 1.3] before B3-2 batch 1).
    #
    # The reference kind is `code`: the values compared against come from the
    # BeamDyn blade deck (a code's tables); the closed form is only how those
    # tables are integrated into a twist profile, so the comparison is not
    # `analytical`.
    assert_relative_error(
        ratio_mean,
        1.0,
        tol=_RATIO_TOL,
        kind="code",
        reference_name=(
            "BeamDyn blade deck of the official IEA 15 MW reference turbine "
            "(IEA-15-240-RWT_BeamDyn_blade.dat, K[5,5] per station), integrated "
            "as a Timoshenko torsion beam"
        ),
        what="shell/beam tip twist ratio over the 0.3-0.9 span band",
    )


def test_s7_beam_reference_integrates_to_the_query_point():
    """The beam twist at ``z`` is ``M * int_0^z dz'/GJ`` up to ``z`` itself, not to the station below.

    The reference used to stop at the last BeamDyn station at or below each query point and drop
    the partial segment.  With stations about 5.85 m apart and GJ falling fast, that staircase
    under-stated the beam twist by 16 % on average (37 % at worst) over the 0.3-0.9 band and
    carried the S-7 ratio from 0.867 to 1.0466, on the wrong side of 1.  A reference that
    integrates to ``z`` is continuous and grows at ``M/GJ(z)`` between stations; both are checked
    against the deck's own stations, linear in ``1/GJ`` (the trapezoid the reference uses).
    See ``odd/tasks/s7-torsion-ratio.md`` section E8.
    """
    from tools.beam_reference import load_beamdyn_blade

    bd = load_beamdyn_blade(
        REPO_ROOT / "tests" / "IEA15MW" / "reference" / "IEA-15-240-RWT_BeamDyn_blade.dat"
    )
    z = np.asarray(bd.z, dtype=np.float64)
    inv = 1.0 / np.asarray(bd.gj, dtype=np.float64)
    k = int(np.searchsorted(z, 0.5 * 117.0))  # first station above mid-span
    z_lo, z_hi = z[k - 1], z[k]
    z_mid = 0.5 * (z_lo + z_hi)
    theta = beam_twist_profile_deg(MOMENT_NM, np.array([z_lo, z_mid, z_hi]))
    inv_mid = np.interp(z_mid, z, inv)
    expected_increment = np.rad2deg(MOMENT_NM * 0.5 * (inv[k - 1] + inv_mid) * (z_mid - z_lo))
    assert theta[1] - theta[0] == pytest.approx(expected_increment, rel=1e-12)
    assert theta[0] < theta[1] < theta[2]


@pytest.mark.slow
def test_s7_applied_moment_matches_request():
    """The force-couple construction must transmit exactly the requested torque.

    Not a window.  A couple construction either delivers ``M`` or it does not, and the tolerance
    that used to live here (``< 0.3``) was wide enough to hide a systematic 10.2 % under-delivery
    whose factor carried the global twist ratio of the case above across 1.0 and inverted the
    direction this module reports -- the bound was five times the effect it was meant to guard.
    See ``odd/tasks/s7-torsion-ratio.md``.
    """
    mesh, _ = build_mesh_and_properties(ELEMENT_SIZE)
    tip_nodes, f_y = tip_couple(mesh, MOMENT_NM)
    x_off = mesh.coords_array[tip_nodes, 0] - mesh.coords_array[tip_nodes, 0].mean()
    mz_applied = float((x_off * f_y).sum())
    assert abs(mz_applied - MOMENT_NM) / MOMENT_NM < 1e-12
    # a couple carries no net force
    assert abs(float(f_y.sum())) < 1e-12 * float(np.abs(f_y).sum())
