"""Does the production realisation of a section moment deliver Bredt's torsion on a closed tube? (P8)

The blade's elastic twist is the section's rigid rotation (``1df95fe``, verdict §22.8). Two
measurements pointed at *how the load is realised on the nodes* rather than at the element: the
production load path's ``distortion/|omega| = 9.36`` (§20.5) and the statically-equivalent end load
on a closed thin-walled tube that bought 125x the section distortion at the same torque
(``d4fec33``). The production ``ForceProjector._distribute`` realises a pure torsional moment as the
constrained **minimum-norm** nodal field, which on the tube is exactly the circle-tangential
``f_j = omega x d_j`` about the section centroid. It is not the Saint-Venant wall shear flow
``q = T/(2A)`` a closed thin-walled section carries, and it over-delivers the Bredt twist rate by
**31.6933x** in this test's interior window - a localised, long-decaying end distortion.

That minimum-norm field is what the fix replaced. ``realise_section_load`` in
``aeroelast.solvers.bem.force_projection`` groups a strip's nodes into physical rings, orders each
ring and applies Bredt's ``q = T/(2A)`` per edge, then feeds any residual through the unchanged
minimum-norm solve. The test is run on the one case whose right answer is exact and independent:
Bredt torsion of the closed thin-walled tube already validated in
``tests/test_thin_walled_tube_torsion.py``. The fixtures, the ``theta_fit`` / ``theta_z`` estimators
and the hand-written reference ``T/GJ`` are reused from that module (``import
test_thin_walled_tube_torsion as tube``); only the realisation changes. Both cases share **one mesh,
one torque and one solver**:

* **case A** - the validated hand-written shear flow ``q = T/(2A)`` at the tip ring and ``-q`` at
  the root ring (``tube._ring_shear_flow``, the independent reference). It stays the independent
  guard: it is not the production code.
* **case B** - the **same** net torque through the production entry point
  ``realise_section_load(strip, F=0, M=T*z_hat, span_dir=z_hat)`` at the tip and root rings. It is
  now the production guard, and it is the code the fix targets: the minimum-norm solve is no longer
  what realises the span moment.

Both cases are compared against the independent Bredt ``T/GJ`` rate inside the suite's 5% rule. The
invariances are asserted without tolerance: the realised nodal system is self-equilibrated (net
force vanishes), its torque read back with the independent ruler ``sum(x Fy - y Fx)`` equals the
requested torque to 1e-9 relative, and its sign follows the requested torque.
"""

from __future__ import annotations

import numpy as np
import pytest

pytest.importorskip("_aeroelast", reason="Rust backend not available")

from aeroelast.solvers.bem.force_projection import _Strip, realise_section_load  # noqa: E402

import tests.validation.parity.test_thin_walled_tube_torsion as tube  # noqa: E402
from tests.support.assertions import assert_relative_error  # noqa: E402

# The suite's 5% rule (CONTRIBUTING), written as a module-level literal because the store's
# extractor resolves a bound declared in this module and not an alias imported from another one:
# `TOL = tube.TOL` left the comparison below unreadable. Same value as
# tests/validation/parity/test_thin_walled_tube_torsion.py::TOL.
TOL = 0.05
WINDOW = (0.4 * tube.L, 0.9 * tube.L)  # the validated interior window of case A


def _ring_strip(coords: np.ndarray, ring: list[int]) -> _Strip:
    """The production ``_Strip`` container for one mesh ring (planar section at one z)."""
    idx = np.asarray(ring, dtype=int)
    pts = coords[idx]
    centroid = pts.mean(axis=0)
    return _Strip(
        node_indices=idx,
        r_center=float(centroid[2]),
        dr=tube.DZ,
        centroid=centroid,
        offsets=pts - centroid,
    )


def _production_tip_moment(coords: np.ndarray, ring: list[int], torque: float) -> np.ndarray:
    """Pure span-axis couple ``torque`` about the ring centroid, realised by production.

    ``F_strip = 0`` and ``M_strip = torque * z_hat``; the returned nodal force system is the
    production :func:`realise_section_load` output, embedded in the assembler's 6-DOF layout.
    """
    strip = _ring_strip(coords, ring)
    nodes = realise_section_load(
        strip,
        np.zeros(3),
        np.asarray([0.0, 0.0, torque]),
        np.asarray([0.0, 0.0, 1.0]),
    )
    f = np.zeros(6 * len(coords))
    for j, nd in enumerate(ring):
        f[6 * nd : 6 * nd + 3] = nodes[j]
    return f


def _measure(coords: np.ndarray, u: np.ndarray, rings: list[list[int]], ref_rate: float) -> dict:
    """Rigid-rotation rate, tip rigid rotation and section distortion for one solution."""
    z = np.asarray([coords[r[0], 2] for r in rings])
    theta_fit = np.asarray([tube._theta_fit(coords, u, r) for r in rings])
    theta_z = np.asarray([tube._theta_z(u, r) for r in rings])
    station = 3 * len(rings) // 4  # a ring at ~0.75 L, away from the self-equilibrium null
    alpha, beta, dev = tube._parallelogram(coords, u, rings[station])
    shear = 0.5 * (alpha + beta)  # the section parallelogram amplitude (in-plane distortion)
    return {
        "z": z,
        "slope_fit": tube._lsq_slope(z, theta_fit, *WINDOW),
        "slope_z": tube._lsq_slope(z, theta_z, *WINDOW),
        "tip_fit": float(theta_fit[-1]),
        "distortion": abs(shear),
        "deviation": dev,
        "distortion_over_rotation": abs(shear) / abs(theta_fit[station]),
    }


def test_production_shear_flow_moment_realisation_vs_bredt():
    """Case A (independent hand-written shear flow) and case B (production realisation) on one tube.

    The decider is the interior rigid-rotation rate against the exact Bredt rate ``T/GJ``. Both
    cases must land inside the 5% rule; a miss is the finding and is never accommodated. Case A is
    the independent guard, case B the production guard.
    """
    coords, conn, rings, n_ring = tube._tube_mesh()
    _, K = tube._assemble(coords, conn, 4, tube._iso_prop())
    modes = tube._rigid_modes(coords)
    gj, ref_rate = tube._bredt_isotropic()

    # ---- case A: the validated hand-written self-equilibrated shear flow (±T at the end rings) ----
    f_a_tip = tube._ring_shear_flow(coords, rings[-1], tube.TORQUE)
    f_a_root = tube._ring_shear_flow(coords, rings[0], -tube.TORQUE)
    assert abs(tube._realised_torque(coords, f_a_tip) - tube.TORQUE) / tube.TORQUE < 1e-9
    assert abs(tube._realised_torque(coords, f_a_root) + tube.TORQUE) / tube.TORQUE < 1e-9
    u_a = tube._solve_rigid_removed(K, f_a_tip + f_a_root, modes)

    # ---- case B: the same net torque through the production shear-flow realisation under test ----
    f_b_tip = _production_tip_moment(coords, rings[-1], +tube.TORQUE)
    f_b_root = _production_tip_moment(coords, rings[0], -tube.TORQUE)
    tip_nodes = f_b_tip.reshape(-1, 6)[np.asarray(rings[-1]), :3]
    root_nodes = f_b_root.reshape(-1, 6)[np.asarray(rings[0]), :3]
    # Tolerance-free invariances: the realised system is self-equilibrated, its torque is what the
    # independent ruler reads, and its sign follows the requested torque. None of these is a
    # comparison against a reference, so none of them is a store row.
    assert np.abs(tip_nodes.sum(axis=0)).max() < 1e-9, "tip net force must vanish"
    assert np.abs(root_nodes.sum(axis=0)).max() < 1e-9, "root net force must vanish"
    realised_tip = tube._realised_torque(coords, f_b_tip)
    realised_root = tube._realised_torque(coords, f_b_root)
    assert abs(realised_tip - tube.TORQUE) / tube.TORQUE < 1e-9, "realised tip torque"
    assert abs(realised_root + tube.TORQUE) / tube.TORQUE < 1e-9, "realised root torque"
    assert realised_tip > 0.0, "a positive M_z about +z must realise a positive torque"
    assert realised_root < 0.0, "a negative M_z about +z must realise a negative torque"
    u_b = tube._solve_rigid_removed(K, f_b_tip + f_b_root, modes)

    m_a = _measure(coords, u_a, rings, ref_rate)
    m_b = _measure(coords, u_b, rings, ref_rate)

    print(f"\nBredt reference: GJ = {gj:.6e} N.m^2, theta' = T/GJ = {ref_rate:.6e} rad/m")
    for label, m in (("A / ring shear flow      ", m_a), ("B / production shear flow", m_b)):
        print(
            f"[{label}] rate ratio = {m['slope_fit'] / ref_rate:.4f} "
            f"(theta_z {m['slope_z'] / ref_rate:.4f}), tip/Bredt(L/2) = "
            f"{m['tip_fit'] / (ref_rate * tube.L / 2):.4f}, "
            f"distortion/|rotation| = {m['distortion_over_rotation']:.3f}, "
            f"deviation-from-rigid = {m['deviation']:.3f}"
        )

    # Case A is the independent guard: the hand-written shear flow, not the production code, against
    # Bredt's T/GJ rate derived from the tube's geometry and isotropic GJ. A miss outside 5% is the
    # finding and is never accommodated.
    assert_relative_error(
        m_a["slope_fit"],
        ref_rate,
        tol=TOL,
        kind="analytical",
        reference_name="Bredt's T/GJ twist rate for the closed tube",
        what="shear-flow moment realisation interior twist rate",
    )
    # Case B is the production guard: the same independent Bredt reference, now through
    # `realise_section_load`. The fix must bring it inside the same 5% rule.
    assert_relative_error(
        m_b["slope_fit"],
        ref_rate,
        tol=TOL,
        kind="analytical",
        reference_name="Bredt's T/GJ twist rate for the closed tube",
        what="production shear-flow moment realisation interior twist rate",
    )
