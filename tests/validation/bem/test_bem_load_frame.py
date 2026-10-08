"""Which frame the BEM loads are in, ours and the published reference's (issue #15).

Two facts that item 3 of the roadmap (#18) turned on, pinned so they cannot drift:

1. **Our production BEM emits the rotor-plane pair.** ``BEMSolver`` wraps ccblade's
   ``distributedAeroLoads``, whose own source computes
   ``cn = cl cos(phi) + cd sin(phi)`` with ``phi`` the inflow angle measured from the
   rotor plane (``phi = alpha + theta``). The identity that says so, in the output's own
   numbers, is ``atan2(Tp, Np) + atan2(Cd, Cl) - theta == alpha``, and it holds to machine
   precision. ``test_force_projection_load_frame.py``'s docstring calls the same pair
   "normal and tangential to the section chord"; that is the *section* frame, and this
   file is what distinguishes them.

2. **Zhou et al. 2025's Fig. 11 pair is in the section frame.** Their published pair,
   read with the official AeroDyn polar at their published angle of attack (Fig. 10),
   reproduces that same angle of attack over the blade core,
   ``alpha == atan2(Tp, Np) + atan2(Cd, Cl)``; read in the rotor plane it is off by up to
   ``9 deg``. The two figures are therefore mutually consistent, and only one frame says
   so. This is what makes the ``+31%`` tangential comparison in the report a frame-mixed
   one rather than a physical difference.

The file is an investigation of frames, not a measurement of our model against an
independent physical reference, so it is declared `out_of_scope` in
`docs/validation/groups.yaml`. The measured numbers, the closure and the citation ban
live in `docs/validation_closures.md` and `docs/validation/gaps.yaml`.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

pytest.importorskip("ccblade", reason="ccblade not installed (pip install -e '.[bem]')")

from aeroelast.solvers.bem.engine import BEMSolver  # noqa: E402

from tests.support.openfast_bem import build_blade_aero_from_aerodyn  # noqa: E402
from tests.support.paths import DATA_DIR  # noqa: E402

AD_PRIMARY = DATA_DIR / "reference" / "iea15mw_openfast" / "case" / "IEA-15-240-RWT_AeroDyn15.dat"
ZHOU_LOADS = DATA_DIR.parent / "docs" / "validation_data" / "zhou_2025_fig11_loads.csv"
ZHOU_AOA = DATA_DIR.parent / "docs" / "validation_data" / "zhou_2025_fig10_aoa.csv"

V_RATED = 10.59
RPM_ZHOU = 7.55
PITCH_RATED = 0.0

#: The blade core, where the digitized Fig. 10 angle of attack is well resolved and the
#: polar is not being extrapolated. The digitized tip is separately excluded.
CORE_LO = 0.26
CORE_HI = 0.80

#: Bounds from the reference's own stated digitization uncertainty (Fig. 10, +/-0.5 deg)
#: with room for polar interpolation, not fitted to the measurement. The measured core
#: residuals are 0.37 deg (section) and 6.83 deg (rotor plane).
SECTION_BOUND_DEG = 1.5
ROTOR_FLOOR_DEG = 3.0

FRAME_IDENTITY_BOUND_DEG = 1e-9


@pytest.fixture(scope="module")
def blade_aero():
    if not AD_PRIMARY.exists():
        pytest.skip(f"AeroDyn reference not present: {AD_PRIMARY}")
    return build_blade_aero_from_aerodyn(AD_PRIMARY)


def _polar_at(aero, index: int, alpha_deg: float) -> tuple[float, float, float]:
    polar = aero.stations[index].airfoil.polars[0]
    alpha = np.degrees(polar.alpha)
    return (
        float(np.interp(alpha_deg, alpha, polar.cl)),
        float(np.interp(alpha_deg, alpha, polar.cd)),
        float(np.interp(alpha_deg, alpha, polar.cm)),
    )


def _zhou_frame_residuals(aero) -> tuple[np.ndarray, np.ndarray]:
    """Per-station (section-reading, rotor-plane-reading) alpha residual [deg]."""
    loads = pd.read_csv(ZHOU_LOADS)
    aoa = pd.read_csv(ZHOU_AOA)
    r_R = loads["r_R"].to_numpy(dtype=float)
    np_zhou = loads["normal_force_kN_m"].to_numpy(dtype=float) * 1.0e3
    tp_zhou = loads["tangential_force_kN_m"].to_numpy(dtype=float) * 1.0e3
    alpha_zhou = np.interp(
        r_R, aoa["r_R"].to_numpy(dtype=float), aoa["angle_of_attack_deg"].to_numpy(dtype=float)
    )
    r_hub = np.asarray(aero.r, dtype=float)
    twist_deg = np.degrees([st.twist for st in aero.stations])
    rotor_radius = float(aero.rotor_radius)

    section: list[float] = []
    rotor: list[float] = []
    for i, rr in enumerate(r_R):
        if not (CORE_LO <= rr <= CORE_HI):
            continue
        r_abs = rr * rotor_radius
        if r_abs <= r_hub[0] or r_abs >= r_hub[-1]:
            continue
        j = int(np.argmin(np.abs(r_hub - r_abs)))
        cl, cd, _ = _polar_at(aero, j, alpha_zhou[i])
        gamma = np.degrees(np.arctan2(cd, cl))
        psi = np.degrees(np.arctan2(tp_zhou[i], np_zhou[i]))
        section.append(psi + gamma - alpha_zhou[i])
        rotor.append(psi + gamma - twist_deg[j] - alpha_zhou[i])
    return np.asarray(section), np.asarray(rotor)


def test_the_production_bem_emits_the_rotor_plane_pair(blade_aero):
    """``atan2(Tp, Np) + atan2(Cd, Cl) - theta == alpha`` on the production output."""
    bem = BEMSolver(blade_aero, rho=1.225, mu=1.81206e-5, hub_height=150.0, shear_exp=0.2).compute(
        V_RATED, RPM_ZHOU, PITCH_RATED
    )
    twist_deg = np.degrees([st.twist for st in blade_aero.stations])
    psi = np.degrees(np.arctan2(bem.Tp, bem.Np))
    gamma = np.degrees(np.arctan2(bem.cd, bem.cl))
    residual = float(np.max(np.abs(psi + gamma - twist_deg - bem.alpha)))
    assert residual < FRAME_IDENTITY_BOUND_DEG, (
        f"production BEM is no longer the rotor-plane pair: max residual {residual:.3e} deg"
    )


def test_zhou_fig11_is_a_section_frame_pair(blade_aero):
    """Their published pair fits their own Fig. 10 angle of attack only in that frame."""
    section, rotor = _zhou_frame_residuals(blade_aero)
    assert section.size >= 8, "the digitized pair lost its blade-core stations"
    section_max = float(np.max(np.abs(section)))
    rotor_max = float(np.max(np.abs(rotor)))
    assert section_max < SECTION_BOUND_DEG, (
        f"their pair no longer reads as a section-frame pair: max {section_max:.2f} deg"
    )
    assert rotor_max > ROTOR_FLOOR_DEG, (
        f"the rotor-plane reading is now as good as the section reading: max {rotor_max:.2f} deg"
    )
