"""S-7 — global torsional response: shell tip torque vs BeamDyn GJ tables.

The shell's global torsional response (tip torque via force couples, section
rotation from the displacement field) is compared against the BeamDyn 6×6 tables:
the station-by-station twist ratio in the 0.3-0.9 span band measures **1.0466** at
``ELEMENT_SIZE = 0.25``, converging to **1.0390** at ``0.125`` — the shell's global
GJ runs about **3.8 % BELOW** the beam's at convergence.

That direction is a correction (2026-10-10).  The ratios recorded here used to read
0.9396 / 0.9336 -- the shell *stiffer* -- because ``run_torsion_case`` normalised its
tip couple by ``2*sum(x_off**2 for x_off > 0)``, which equals the full chordwise second
moment only for an x-symmetric ring and under-delivered the applied moment by a
mesh-dependent factor κ (0.8978 at 0.25, 0.9154 on the coarser three).  The ratio was
therefore a *stiffness* ratio times κ, and κ is what carried it across 1.0.  Corrected
by dividing each recorded ratio by its own κ -- exact, because the response is linear
in the applied moment -- the trend is monotone to 1.039 and never crosses 1.0.

The sectional S-1 comparison still reports the shell's GJ *above* the BeamDyn tables, so
sectional and global now disagree in sign.  That tension is real and open; this module
reports the measurement and does not settle it.  See ``odd/tasks/s7-torsion-ratio.md``.

The bound is a ±30 % window around exact agreement, which is the documented S-1 model
difference with room for the coarse-mesh envelope, and not a fit to what this run
happens to produce.  Note it is five times the effect it names: the bound cannot detect
the 10.2 % load error that used to invert the direction, which is why the couple
construction carries its own exact test.
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
# delivering exactly MOMENT_NM (see tip_couple).  es=2.0/1.0/0.5/0.25 were re-measured on
# 2026-10-10 with that construction; es=0.125 is the HEAD-measured 0.933641 divided by its
# geometry-only kappa = 0.898586, which is exact because the response is linear in the applied
# moment.  The analytic kappa correction reproduces the four re-measured rows to 5-6 digits.
#
#   es=2.000   1460 nodes  ratio=1.118319
#   es=1.000   3040 nodes  ratio=1.116392
#   es=0.500   9271 nodes  ratio=1.069529
#   es=0.250  32325 nodes  ratio=1.046588   <- the fixture
#   es=0.125 120352 nodes  ratio=1.039011   <- converged (230 s, 12 GB peak RSS)
#
# 0.250 sits 0.73% above 0.125, which is what justifies it as a regression guard.  This is a
# slow test (about 1 min of PETSc per run), which is why it is marked slow.
#
# History: before 2026-10-10 these read 1.0237 / 1.0219 / 0.9790 / 0.939603 / 0.933641 -- the
# ratio crossing 1.0 between 0.500 and 0.250, which is what made "the shell is stiffer" look
# settled.  A table committed here between 2026-09-30 and 2026-10-04 read 1.386 / 1.406 / 1.303 /
# 1.286 / 1.273; its *direction* was right, and its magnitude still does not reproduce at HEAD at
# any mesh size, with or without the element-winding canonicalisation that `22d3ccc` removed.
# See `odd/tasks/s7-torsion-ratio.md`.
ELEMENT_SIZE = 0.25

#: The band [0.7, 1.3] is symmetric around exact agreement with the beam, so it is
#: the same statement as a relative bound of 0.30 on the ratio.
_RATIO_TOL = 0.30


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

    # S-1's sectional shell GJ sits +17-24% ABOVE the BeamDyn tables, so a shell that twisted
    # less would sit below 1 here.  It does not: at the correct load the global ratio is 1.0466,
    # i.e. the shell twists 4.7% MORE -- globally softer while sectionally stiffer.  The two
    # comparisons are not the same quantity and they now disagree in sign; this test reports the
    # global measurement, it does not adjudicate the tension (see the module docstring).
    #
    # The bound is the documented S-1 model difference plus the coarse-mesh envelope, and not a
    # fit to what this run happens to produce: [0.7, 1.3] is a +-30% window around exact
    # agreement, which is what assert_relative_error expresses symmetrically.
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
