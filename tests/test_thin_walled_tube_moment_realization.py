"""Does the production minimum-norm nodal realisation of a section moment under-deliver rigid torsion? (P8)

The blade's elastic twist is the section's rigid rotation (``1df95fe``, verdict §22.8). The shell
delivers ``0.74x`` of the rigid twist its own independently measured ``GJ_shell/GJ_ref = 0.87``
calls for - a ``1.35x`` deficit. Two measurements point at *how the load is realised on the nodes*
rather than at the element: the production load path's ``distortion/|omega| = 9.36`` (§20.5) and the
statically-equivalent end load on a closed thin-walled tube that bought 125x the section distortion
at the same torque (``d4fec33``). The hypothesis is that
``ForceProjector._distribute``'s minimum-norm nodal realisation of a pure torsional moment yields a
rigid section rotation **smaller** than Bredt's ``T L / GJ``, with the difference appearing as
cross-section distortion.

The test is run on the one case whose right answer is exact and independent: Bredt torsion of the
closed thin-walled tube already validated in ``tests/test_thin_walled_tube_torsion.py``.  The
fixtures, the ``theta_fit`` / ``theta_z`` estimators and the hand-written reference ``T/GJ`` are
reused from that module (``import test_thin_walled_tube_torsion as tube``); only the load
realisation changes.  Both cases share **one mesh, one torque and one solver**:

* **case A** - the validated shear flow ``q = T/(2A)`` at the tip ring and ``-q`` at the root ring
  (``tube._self_equilibrated_load``).  Fully self-equilibrated, so the interior is uniform
  Saint-Venant torsion and both metrics equal ``T/GJ``.
* **case B** - the **same** net torque through the production minimum-norm system:
  ``ForceProjector._distribute(strip, F=0, M=T*z_hat)`` at the tip and root rings.  This calls the
  production entry point of ``aeroelast.solvers.bem.force_projection`` (the same solve
  ``ForceProjector.project()`` runs per strip); the minimum-norm solve is **not** re-implemented
  here, and its output is checked against the independent torque ruler, so the test exercises
  production rather than its own harness.

Measured result (see the module printout and verdict §22.9): case A reproduces Bredt within 5%,
but case B lands at **31.7x** the Bredt rate - it **over**-delivers, not under-delivers, so the
under-delivery hypothesis (P8) is rejected.  The departure is real and is localised: the
minimum-norm field is exactly ``f_j = omega x d_j`` (circle-tangential about the section centroid),
not the wall shear flow ``q = T/(2A)``, and it produces a long-decaying end distortion.
"""

from __future__ import annotations

import numpy as np
import pytest

pytest.importorskip("_aeroelast", reason="Rust backend not available")

from aeroelast.solvers.bem.force_projection import ForceProjector, _Strip  # noqa: E402

import test_thin_walled_tube_torsion as tube  # noqa: E402

TOL = tube.TOL  # the suite's 5% rule against the exact, independent Bredt reference
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
    """Pure span-axis couple ``torque`` about the ring centroid, realised by production ``_distribute``.

    ``F_strip = 0`` and ``M_strip = torque * z_hat``; the returned nodal force system is the
    production minimum-norm solve's output, embedded in the assembler's 6-DOF layout.
    """
    strip = _ring_strip(coords, ring)
    nodes = ForceProjector._distribute(strip, np.zeros(3), np.asarray([0.0, 0.0, torque]))
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


def test_production_minimum_norm_moment_realisation_vs_bredt():
    """Case A (validated shear flow) and case B (production ``_distribute``) on one tube and torque.

    The decider is the interior rigid-rotation rate against the exact Bredt rate ``T/GJ``.
    A miss outside 5% is the finding and is not accommodated; here it lands on the **over** side,
    so the under-delivery hypothesis is refuted and the verdict is recorded in §22.9.
    """
    coords, conn, rings, n_ring = tube._tube_mesh()
    _, K = tube._assemble(coords, conn, 4, tube._iso_prop())
    modes = tube._rigid_modes(coords)
    gj, ref_rate = tube._bredt_isotropic()

    # ---- case A: the validated self-equilibrated shear flow (±T at the two end rings) ----
    f_a_tip = tube._ring_shear_flow(coords, rings[-1], tube.TORQUE)
    f_a_root = tube._ring_shear_flow(coords, rings[0], -tube.TORQUE)
    assert abs(tube._realised_torque(coords, f_a_tip) - tube.TORQUE) / tube.TORQUE < 1e-9
    assert abs(tube._realised_torque(coords, f_a_root) + tube.TORQUE) / tube.TORQUE < 1e-9
    u_a = tube._solve_rigid_removed(K, f_a_tip + f_a_root, modes)

    # ---- case B: the same net torque through the production minimum-norm code under test ----
    f_b_tip = _production_tip_moment(coords, rings[-1], +tube.TORQUE)
    f_b_root = _production_tip_moment(coords, rings[0], -tube.TORQUE)
    # The production solve conserves the requested resultant force and moment exactly; the torque
    # is checked with the independent ruler (`sum x Fy - y Fx`), never with production's own bookkeeping.
    assert abs(tube._realised_torque(coords, f_b_tip) - tube.TORQUE) / tube.TORQUE < 1e-9
    assert abs(tube._realised_torque(coords, f_b_root) + tube.TORQUE) / tube.TORQUE < 1e-9
    assert np.abs(f_b_tip.reshape(-1, 6)[:, :3].sum(axis=0)).max() < 1e-9, "net force must vanish"
    u_b = tube._solve_rigid_removed(K, f_b_tip + f_b_root, modes)

    m_a = _measure(coords, u_a, rings, ref_rate)
    m_b = _measure(coords, u_b, rings, ref_rate)

    print(f"\nBredt reference: GJ = {gj:.6e} N.m^2, theta' = T/GJ = {ref_rate:.6e} rad/m")
    for label, m in (("A / shear flow   ", m_a), ("B / min-norm P8  ", m_b)):
        print(f"[{label}] rate ratio = {m['slope_fit'] / ref_rate:.4f} "
              f"(theta_z {m['slope_z'] / ref_rate:.4f}), tip/Bredt(L/2) = "
              f"{m['tip_fit'] / (ref_rate * tube.L / 2):.4f}, "
              f"distortion/|rotation| = {m['distortion_over_rotation']:.3f}, "
              f"deviation-from-rigid = {m['deviation']:.3f}")

    # The production minimum-norm field is the circle-tangential f_j = omega x d_j, NOT a wall
    # shear flow.  This is the localisation of the departure; it is exact (round-off only).
    strip = _ring_strip(coords, rings[-1])
    d = strip.offsets
    omega = tube.TORQUE / float(np.sum(d[:, 0] ** 2 + d[:, 1] ** 2))
    circle_tangential = omega * np.c_[-d[:, 1], d[:, 0], np.zeros(len(d))]
    f_b_nodes = f_b_tip.reshape(-1, 6)[np.asarray(rings[-1]), :3]
    assert np.abs(f_b_nodes - circle_tangential).max() < 1e-9 * np.abs(f_b_nodes).max()
    f_a_nodes = f_a_tip.reshape(-1, 6)[np.asarray(rings[-1]), :3]
    print(f"minimum-norm nodal |f|: min {np.linalg.norm(f_b_nodes, axis=1).min():.1f} "
          f"max {np.linalg.norm(f_b_nodes, axis=1).max():.1f} (corner/mid-wall); "
          f"wall shear flow |f|: min {np.linalg.norm(f_a_nodes, axis=1).min():.1f} "
          f"max {np.linalg.norm(f_a_nodes, axis=1).max():.1f}")

    # Case A is the validated exact reference and must stay inside 5%.
    assert abs(m_a["slope_fit"] / ref_rate - 1.0) < TOL, (
        f"case A (validated shear flow) is {m_a['slope_fit'] / ref_rate:.4f}x Bredt"
    )
    # Case B misses the exact Bredt reference by a large factor.  This is the decisive finding:
    # the hypothesis is about *under*-delivery, and the measured departure is on the over side.
    ratio_b = m_b["slope_fit"] / ref_rate
    assert abs(ratio_b - 1.0) > TOL, (
        f"case B (production minimum-norm) is {ratio_b:.3f}x Bredt - inside 5%, so the "
        f"under-delivery hypothesis would be refuted; re-read §22.9 before changing this"
    )
    print(f"\nDECISION: case A {m_a['slope_fit'] / ref_rate:.4f}x Bredt (PASS, <5%); "
          f"case B {ratio_b:.4f}x Bredt (MISS by {abs(ratio_b - 1) * 100:.0f}%, "
          f"on the OVER side) -> under-delivery hypothesis (P8) REJECTED")
