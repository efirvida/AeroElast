"""Arbitrate the twist-feedback sense against the frame production actually runs.

`BEMFSIParticipant` derives its mesh-rotation -> BEM-twist factor from
``normal_direction x (span_direction x tangential_direction)``.  That conflates two
different things the configured ``tangential_direction`` carries: the **sense of the
tangential force** ``Tp`` (line ~1020, `F = F_n*n + F_t*t`) and a **proxy for the
chord axis** used by the twist sense.

The static validation fixtures never pass ``tangential_direction``, so they run the
sign that makes the derivation's premise true (``span x tangential = normal``).  Every
production case YAML (62 of them) sets ``tangential_direction = [-1, 0, 0]``, which
makes ``span x tangential = -normal`` and therefore *flips* the derived factor.

This prints, on the real blade mesh and the real AeroDyn deck, for both frames:

  * the derived ``_twist_mesh_to_bem``;
  * ``_strip_chord_dirs`` and the deck-matched ``LE . c_hat``;
  * the elastic twist a rigid +2 deg section rotation about +span produces, i.e.
    which way CCBlade's ``theta`` moves;
  * the thrust CCBlade returns for that rotation, so de-loading vs re-loading is
    measured, not argued.

The physical statement being arbitrated: the deck's leading edge sits at +x, the
load frame puts downwind (thrust) at +y, so a **positive** rotation about +span moves
the leading edge downwind - the section pitches **nose-down**, alpha falls, and
CCBlade's ``alpha = phi - theta`` requires ``theta`` to **rise**.  The factor must
therefore be ``+1`` in the production frame, and anything else turns the corrected
nose-down rotation into the re-loading the code comment warns about.

Usage (login node):
    module load glu gcc/14.2.0_sequana
    python tools/diagnose_twist_sense_frame.py [--element-size 1.0]
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from tests.support.openfast_bem import build_blade_aero_from_aerodyn  # noqa: E402
from tests.support.paths import DATA_DIR  # noqa: E402

YAML = DATA_DIR / "IEA-15-240-RWT.yaml"
AD_PRIMARY = DATA_DIR / "reference" / "iea15mw_openfast" / "case" / "IEA-15-240-RWT_AeroDyn15.dat"

V_RATED, RPM_RATED, PITCH_RATED = 10.59, 7.56, 0.0
BASE_BEM = {
    "wind_speed": V_RATED,
    "omega": RPM_RATED * 2.0 * np.pi / 60.0,
    "pitch": PITCH_RATED,
    "azimuth": 0.0,
    "air_density": 1.225,
    "dynamic_viscosity": 1.81206e-5,
    "hub_height": 150.0,
    "shear_exp": 0.0,
}

#: The frame every production case YAML carries (62 files under tests/).
PRODUCTION_FRAME = {
    "span_direction": [0.0, 0.0, 1.0],
    "normal_direction": [0.0, 1.0, 0.0],
    "tangential_direction": [-1.0, 0.0, 0.0],
}
#: The frame the static fixtures run by default.
DEFAULT_FRAME: dict[str, list[float]] = {}

ANGLE_DEG = 2.0


def rigid_section_rotation(coords: np.ndarray, span: np.ndarray, angle_rad: float) -> np.ndarray:
    """Displacement field of a rigid rotation about ``+span`` through the origin.

    A pure rotation about the span axis is span-invariant and, for any ring, its
    in-plane part is exactly a section rotation - so ``_section_rotation`` must
    return ``+angle_rad`` and nothing else.
    """
    axis = np.asarray(span, dtype=float)
    axis = axis / np.linalg.norm(axis)
    return angle_rad * np.cross(np.broadcast_to(axis, coords.shape), coords)


def describe(label: str, bem_config: dict) -> dict:
    from aeroelast.models.blade.model import Blade
    from aeroelast.solvers.bem.fsi_participant import BEMFSIParticipant

    model = Blade(str(YAML), element_size=ARGS.element_size)
    model.generate_mesh()
    mesh = model.mesh
    assert mesh is not None, "blade mesh generation failed"
    blade_aero = build_blade_aero_from_aerodyn(AD_PRIMARY)
    participant = BEMFSIParticipant(mesh, blade_aero, {**BASE_BEM, **bem_config})

    coords = np.array([[nd.x, nd.y, nd.z] for nd in mesh.nodes], dtype=float)
    span = participant._span_dir
    normal = participant._normal_dir

    displacement = rigid_section_rotation(coords, span, np.deg2rad(ANGLE_DEG))
    r_def, twist_def = participant._compute_deformed_geometry(displacement)
    delta = twist_def - participant._ref_twist

    # deck side: the canonical WindIO aerofoil tables put x = 0 at the leading edge,
    # so a ring's leading-edge direction is the one toward small deck x.
    proj = participant._projector
    k_mid = len(proj._strips) // 2
    chord_hat = np.asarray(proj._strip_chord_dirs[k_mid], dtype=float)
    ring = proj._station_ring_points(k_mid, proj._strips[k_mid], span)
    x_proj = (ring - ring.mean(axis=0)) @ chord_hat
    le_dir = -chord_hat if x_proj.max() > -x_proj.min() else chord_hat

    # physical indicator: does +span rotation move the leading edge downwind?
    nose_down_indicator = float(le_dir @ np.cross(normal, span))

    n_valid = sum(
        1
        for k in range(len(proj._strips))
        if participant._section_rotation(displacement, participant._strip_ring_indices[k])
        is not None
    )

    solver, _ = participant._rebuild_bem_solver(r_def, twist_def)
    result = solver.compute(V_RATED, RPM_RATED, PITCH_RATED)
    ref_solver, _ = participant._rebuild_bem_solver(participant._ref_r, participant._ref_twist)
    ref_result = ref_solver.compute(V_RATED, RPM_RATED, PITCH_RATED)

    print(f"\n=== {label} ===")
    print(f"  span_direction       = {span}")
    print(f"  normal_direction     = {normal}")
    print(f"  tangential_direction = {participant._tangential_dir}")
    print(
        "  normal . (span x tangential) = "
        f"{normal @ np.cross(span, participant._tangential_dir):+.6f}"
    )
    print(f"  _twist_mesh_to_bem   = {participant._twist_mesh_to_bem:+.0f}")
    print(f"  chord_hat (mid ring) = {chord_hat}")
    print(f"  leading-edge dir     = {le_dir}")
    print(
        f"  LE . (normal x span) = {nose_down_indicator:+.6f}  (+1 => +span rotation is nose-down)"
    )
    print(f"  strips with a usable ring: {n_valid}/{len(proj._strips)}")
    print(
        f"  rigid +{ANGLE_DEG} deg section rotation -> twist_def - ref_twist "
        f"spans [{delta.min():+.4f}, {delta.max():+.4f}] deg"
    )
    print(
        f"  thrust: reference {ref_result.thrust * 1e-6:.4f} MN  |  deformed {result.thrust * 1e-6:.4f} MN"
    )
    print(f"  delta thrust: {100.0 * (result.thrust / ref_result.thrust - 1.0):+.3f} %")
    return {
        "factor": participant._twist_mesh_to_bem,
        "delta_thrust": result.thrust / ref_result.thrust - 1.0,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--element-size", type=float, default=1.0)
    ARGS = parser.parse_args()

    print(f"blade mesh element_size = {ARGS.element_size} m")
    default = describe("DEFAULT frame (what the static fixtures run)", DEFAULT_FRAME)
    production = describe("PRODUCTION frame (tangential_direction = -x)", PRODUCTION_FRAME)

    print("\n=== verdict ===")
    print(
        f"  default frame   factor {default['factor']:+.0f}  delta thrust {100 * default['delta_thrust']:+.3f} %"
    )
    print(
        f"  production frame factor {production['factor']:+.0f}  "
        f"delta thrust {100 * production['delta_thrust']:+.3f} %"
    )
