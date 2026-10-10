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

from aeroelast.solvers.bem.force_projection import (  # noqa: E402
    ForceProjector,
    _ring_tributary_weights,
    _Strip,
    realise_section_load,
)

import tests.validation.parity.test_thin_walled_tube_torsion as tube  # noqa: E402
from tests.support.assertions import assert_relative_error  # noqa: E402

# The suite's 5% rule (CONTRIBUTING), written as a module-level literal because the store's
# extractor resolves a bound declared in this module and not an alias imported from another one:
# `TOL = tube.TOL` left the comparison below unreadable. Same value as
# tests/validation/parity/test_thin_walled_tube_torsion.py::TOL.
TOL = 0.05
WINDOW = (0.4 * tube.L, 0.9 * tube.L)  # the validated interior window of case A

# The legacy minimum-norm realisation's recorded over-delivery of Bredt on this fixture
# (`ForceProjector._distribute` realising a pure span couple as the circle-tangential field
# `f_j = omega x d_j`), promoted from the module docstring to an executed regression pin by T3 of
# issue #30: `31.69331` with the rigid modes removed and `10.40654` with a clamped root.  These are
# **recorded measurements of a defect**, not validated bounds - Bredt is the reference the field
# *misses* - so the pin is deliberately kept out of the store's reference machinery and the defect
# stays recorded in `docs/validation/gaps.yaml` (`moment_realization_over_delivers`).
MIN_NORM_OVER_DELIVERY_SELF_EQUILIBRATED = 31.69331
MIN_NORM_OVER_DELIVERY_CLAMPED = 10.40654
PIN_TOL = 1e-3

# T5 of issue #30: production does not apply the shear flow at the end rings.  `project()` hands
# `realise_section_load` a whole **strip** and `_realise_multi_cell_section_load` splits its torsion
# over the strip's physical rings with `_ring_tributary_weights`, each ring realising its own share
# with its own wall flow - so the internal torque ramps and a self-equilibrated couple is applied at
# every ring.  That mode is not covered by the end-couple bound above, and the two metrics disagree
# at a fixed mesh: the estimator-free energy balance `W / integral M_cum^2/(2 GJ)` sits a few percent
# over Bredt and **converges with refinement** (4.356 % at `n_z = 30`, 2.484 % at `n_z = 60`), while
# the local `theta'(z) * GJ / M_cum(z)` blows up near the tip where `M_cum -> 0`.  Neither number is
# therefore registered as a validated bound; they are bare regression pins, and the parked
# minimum-norm field's non-converging `+338 %` in the same mode is pinned next to them.
DISTRIBUTED_ENERGY_EXCESS = {30: 0.04356, 60: 0.02484}
MIN_NORM_DISTRIBUTED_ENERGY_EXCESS = {30: 3.38184, 60: 3.34399}


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


def _min_norm_ring_couple(coords, ring, torque):
    """The legacy minimum-norm realisation of a pure span couple about one ring's centroid.

    What ``ForceProjector._distribute`` returns for ``M = torque * z_hat`` on that ring, embedded in
    the mesh's 6-DOF layout. This is the parked production default, i.e. the field the shear-flow
    fix replaced.
    """
    nodes = ForceProjector._distribute(
        _ring_strip(coords, ring), np.zeros(3), np.asarray([0.0, 0.0, torque])
    )
    f = np.zeros(6 * len(coords))
    for j, node in enumerate(ring):
        f[6 * node : 6 * node + 3] = nodes[j]
    return f


def test_minimum_norm_realisation_over_delivers_bredt_by_the_recorded_factor():
    """Do the two moment realisations stay apart by the factor they were recorded at? (T3, #30)

    The companion assertion to ``test_production_shear_flow_moment_realisation_vs_bredt``: the same
    tube, the same torque and the same solver, but the span couple realised through the legacy
    ``ForceProjector._distribute`` instead of the wall shear flow. That field is the
    circle-tangential ``f_j = omega x d_j`` about the section centroid, so on a closed section it
    over-delivers Bredt by **31.69331x** with the rigid modes removed and **10.40654x** with a
    clamped root - a localised, long-decaying end distortion, not a torsion rate.

    This is a **regression pin**, not a validated bound: it asserts that the recorded defect has not
    moved, so a change to ``_distribute`` (or to the fixture) that silently shifted it makes noise.
    It deliberately does **not** go through ``assert_relative_error`` - the store's comparisons need
    an independent reference the measurement must *meet*, and Bredt is the one this field misses -
    so the defect stays recorded in ``docs/validation/gaps.yaml``
    (``moment_realization_over_delivers``) rather than dressed as validation.
    """
    coords, conn, rings, _n_ring = tube._tube_mesh()
    _, K = tube._assemble(coords, conn, 4, tube._iso_prop())
    modes = tube._rigid_modes(coords)
    _gj, ref_rate = tube._bredt_isotropic()

    # ---- the tolerance-free invariances the field must still satisfy: it is a self-equilibrated
    #      couple whose torque is what the independent ruler reads, with the requested sign.
    f_tip = _min_norm_ring_couple(coords, rings[-1], +tube.TORQUE)
    f_root = _min_norm_ring_couple(coords, rings[0], -tube.TORQUE)
    tip_nodes = f_tip.reshape(-1, 6)[np.asarray(rings[-1]), :3]
    root_nodes = f_root.reshape(-1, 6)[np.asarray(rings[0]), :3]
    assert np.abs(tip_nodes.sum(axis=0)).max() < 1e-9, "tip net force must vanish"
    assert np.abs(root_nodes.sum(axis=0)).max() < 1e-9, "root net force must vanish"
    assert abs(tube._realised_torque(coords, f_tip) - tube.TORQUE) / tube.TORQUE < 1e-9
    assert abs(tube._realised_torque(coords, f_root) + tube.TORQUE) / tube.TORQUE < 1e-9

    # ---- self-equilibrated (+T tip / -T root), rigid modes removed: the recorded configuration
    u_eq = tube._solve_rigid_removed(K, f_tip + f_root, modes)
    m_eq = _measure(coords, u_eq, rings, ref_rate)
    ratio_eq = m_eq["slope_fit"] / ref_rate

    # ---- clamped root, +T at the tip: the same fixture as case A/B's other end condition
    clamped = tube._clamped_dofs(rings)
    u_cl = tube._solve_clamped(K, f_tip, clamped)
    m_cl = _measure(coords, u_cl, rings, ref_rate)
    ratio_cl = m_cl["slope_fit"] / ref_rate

    print(
        f"\nminimum-norm over-delivery of Bredt: self-equilibrated {ratio_eq:.5f}x "
        f"(pinned {MIN_NORM_OVER_DELIVERY_SELF_EQUILIBRATED}), clamped {ratio_cl:.5f}x "
        f"(pinned {MIN_NORM_OVER_DELIVERY_CLAMPED}); "
        f"distortion/|rotation| = {m_eq['distortion_over_rotation']:.3f}"
    )
    assert ratio_eq == pytest.approx(MIN_NORM_OVER_DELIVERY_SELF_EQUILIBRATED, rel=PIN_TOL), (
        "the minimum-norm over-delivery moved: it is a recorded defect (gaps.yaml "
        "moment_realization_over_delivers), so a shift is a finding to report, never a bound "
        "to widen"
    )
    assert ratio_cl == pytest.approx(MIN_NORM_OVER_DELIVERY_CLAMPED, rel=PIN_TOL), (
        "the clamped minimum-norm over-delivery moved: same recorded defect as above"
    )


def _ring_couple(coords, ring, torque):
    """The validated single-ring production route, embedded in the mesh's 6-DOF layout."""
    idx = np.asarray(ring, dtype=int)
    nodes = _production_tip_moment(coords, ring, torque).reshape(-1, 6)
    f = np.zeros(6 * len(coords))
    for node in idx:
        f[6 * node : 6 * node + 3] = nodes[node, :3]
    return f


def _distributed_ring_by_ring(coords, rings, z):
    """Production's structure: one tributary share of the total torsion per ring, its own wall flow."""
    weights = _ring_tributary_weights(z)
    f = np.zeros(6 * len(coords))
    for ring, weight in zip(rings, weights, strict=True):
        f += _ring_couple(coords, ring, tube.TORQUE * float(weight))
    return f


def _distributed_one_solve(coords, dr):
    """The parked default's branch: one minimum-norm solve for the whole moment over all nodes."""
    centroid = coords.mean(axis=0)
    strip = _Strip(
        node_indices=np.arange(len(coords), dtype=np.intp),
        r_center=float(centroid[2]),
        dr=dr,
        centroid=centroid,
        offsets=coords - centroid,
    )
    nodes = ForceProjector._distribute(strip, np.zeros(3), np.array([0.0, 0.0, tube.TORQUE]))
    f = np.zeros(6 * len(coords))
    f[0::6], f[1::6], f[2::6] = nodes[:, 0], nodes[:, 1], nodes[:, 2]
    return f


def _cumulative_moment(coords, f, rings) -> np.ndarray:
    """Outboard cumulative applied moment per ring, read back from the applied field itself."""
    force = f.reshape(-1, 6)[:, :3]
    zz = coords[:, 2]
    profile = []
    for ring in rings:
        centre = coords[ring].mean(axis=0)
        outboard = zz > float(coords[ring][:, 2].mean())
        profile.append(float(np.cross(coords[outboard] - centre, force[outboard]).sum(axis=0)[2]))
    return np.asarray(profile)


def _energy_excess(coords, f, u, rings, z, gj) -> float:
    """``W / integral M_cum^2/(2 GJ) - 1``: estimator-free, so it cannot repeat a local bias."""
    work = 0.5 * float(f @ u)
    bredt = float(np.trapezoid(_cumulative_moment(coords, f, rings) ** 2, z)) / (2.0 * gj)
    return work / bredt - 1.0


def test_production_distributed_application_is_bredt_to_a_converging_residual():
    """T5 of #30: is the application mode production actually runs Bredt-consistent? (bare pins)

    The bound in ``test_production_shear_flow_moment_realisation_vs_bredt`` is measured with the
    shear flow on the **end rings**, so the internal torque is constant and the interior window sits
    far from the application. Production applies a tributary share at **every** ring of the strip,
    which is a different mode: the torque ramps and the per-ring self-equilibrated couples overlap
    (a Saint-Venant disturbance decays over ~one section dimension, here ~1 m, against `dz = 0.2`).

    This test measures that mode with the estimator-free energy balance, and pins two facts:

    * the wall flow's excess **shrinks with refinement** (4.356 % -> 2.484 % when `dz` halves): its
      residual is the discretisation of applying a discrete couple per ring, not a defect of the
      realisation. It is deliberately **not** turned into a validated bound - at a fixed mesh the
      local ``theta'(z) * GJ / M_cum(z)`` metric disagrees with the energy one and blows up near the
      tip where ``M_cum -> 0``, so a bound here would be metric-shopping;
    * the parked minimum-norm field's excess in the same mode is **+338 % and does not move**
      (3.38184 -> 3.34399), i.e. its over-delivery is not a discretisation effect at all.

    Read with the bound above: the validated end-couple number and the pin here answer different
    questions, and the production mode sits between them.
    """
    measured = {}
    for n_z in sorted(DISTRIBUTED_ENERGY_EXCESS):
        coords, conn, rings, _n_ring = tube._tube_mesh(n_z=n_z)
        _, K = tube._assemble(coords, conn, 4, tube._iso_prop())
        gj, _ref_rate = tube._bredt_isotropic()
        clamped = tube._clamped_dofs(rings)
        z = np.asarray([coords[ring][:, 2].mean() for ring in rings])

        wall = _distributed_ring_by_ring(coords, rings, z)
        assert abs(tube._realised_torque(coords, wall) - tube.TORQUE) / tube.TORQUE < 1e-9
        wall_excess = _energy_excess(
            coords, wall, tube._solve_clamped(K, wall, clamped), rings, z, gj
        )
        plain = _distributed_one_solve(coords, tube.L / n_z)
        plain_excess = _energy_excess(
            coords, plain, tube._solve_clamped(K, plain, clamped), rings, z, gj
        )
        measured[n_z] = (wall_excess, plain_excess)
        print(
            f"[n_z={n_z}, dz={tube.L / n_z:.3f} m] distributed wall flow energy excess "
            f"{wall_excess:+.5f}, minimum-norm {plain_excess:+.5f}"
        )

    for n_z, (wall_excess, plain_excess) in measured.items():
        assert wall_excess == pytest.approx(DISTRIBUTED_ENERGY_EXCESS[n_z], rel=0.02), (
            "the distributed wall-flow energy excess moved: it is a pinned discretisation residual "
            "of the production application mode (T5 of #30), not a bound to widen"
        )
        assert plain_excess == pytest.approx(MIN_NORM_DISTRIBUTED_ENERGY_EXCESS[n_z], rel=0.02), (
            "the distributed minimum-norm energy excess moved: this is the parked default's "
            "recorded defect (gaps.yaml moment_realization_over_delivers)"
        )
    coarse, fine = sorted(measured)
    assert measured[fine][0] < measured[coarse][0], (
        "the distributed wall-flow excess must shrink with refinement - that is what makes it a "
        "discretisation residual of a discrete per-ring application rather than a defect"
    )
    assert measured[fine][1] > 0.95 * measured[coarse][1], (
        "the minimum-norm excess is not expected to shrink: it is the realisation's own defect"
    )
