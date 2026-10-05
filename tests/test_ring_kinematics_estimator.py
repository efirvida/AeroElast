"""Instrument validation for the suite's section-rotation estimators - GitHub issue #11.

This file measures the *measuring instrument*.  It does **not** change any estimator and it
does **not** test the physics: every case is built from an analytic displacement field whose
answer is known exactly, independent of any production projector, mesh generator or solver.

Why it exists
-------------
The blade's reported tip section rotation is dominated by distortion (the module
``tests/validation/blade/test_blade_rated_twist.py`` prints ``distortion / |omega| = 42``,
where the same estimator previously reported ``8.6``).  Before any physics conclusion is drawn
from that number we must know whether the estimator recovers a rigid rotation when there is one
and refuses to invent one when there is not.

Estimator inventory (``path:line`` and calling convention)
---------------------------------------------------------
``tests/validation/blade/test_blade_rated_twist.py:459``  ``_ring_kinematics(coords, u, ring)``
    ``u`` is the global 6-DOF vector ``[ux, uy, uz, rx, ry, rz]`` per node.  It uses **only the
    in-plane translations** (DOF 0 and 1).  Full 6-parameter affine least-squares fit with an
    intercept; ``omega = 0.5 * (a21 - a12)``.  The nodal rotations are never read.  Because the
    fit carries an intercept, ``omega`` is the translation-invariant antisymmetric part of the
    displacement gradient and is independent of the ring's location.

``tests/validation/blade/test_blade_rated_twist.py:881``  ``_measured_tip_axes(coords, ring)``
    Geometry only: the SVD principal axes of the ring outline.  It takes **no displacement** and
    returns axes, not a rotation, so it is not a section-rotation estimator.

``tests/validation/parity/test_thin_walled_tube_torsion.py:344``  ``_theta_fit(coords, u, ring)``
    ``u`` is the global 6-DOF vector; uses only DOF 0 and 1.  Centroid estimate
    ``theta = sum(x*uy - y*ux) / sum(x**2 + y**2)`` with the mean removed from ``u`` but **not**
    from the coordinates, so it is exact only for a ring centred at the coordinate origin (the
    tube axis).  Nodal rotations are unused.

``tests/validation/parity/test_thin_walled_tube_torsion.py:355``  ``_theta_z(u, ring)``
    Reads the nodal rotation DOF 5 (``rz``) directly and averages it.  It never looks at the
    translations; it reports whatever rotation was prescribed or solved for.

``tests/validation/parity/test_thin_walled_tube_torsion.py:360``  ``_parallelogram(coords, u, ring)``
    ``u`` is the global 6-DOF vector; uses DOF 0 and 1.  Fits ``ux = alpha*y`` and ``uy = beta*x``
    after removing the mean of ``u``; returns ``(alpha, beta, deviation)``.  It is the suite's
    distortion metric; its internal rotation is the same centroid estimate as ``_theta_fit``.

``tests/validation/blade/test_blade_section_torsion_stiffness.py:145``
    ``_ring_section_rotation(coords, u, ring)``.  Global 6-DOF ``u``; affine fit on DOF 0 and 1;
    identical estimator to ``_ring_kinematics`` (imported lazily - the module pulls in
    ``openfast_toolbox``).

``tests/validation/blade/test_blade_twist_anchor_beam.py:85`` and
``tests/validation/rotor/test_bem_fsi_deformed_geometry.py:149``
    ``_ring_section_rotation(coords, disp, ring)``.  ``disp`` is an ``(N, 3)`` displacement
    array; affine fit on columns 0 and 1.

``tests/validation/blade/test_laminate_bend_twist.py:182/189/193``
    ``twist_lsq`` / ``twist_edge`` / ``twist_rotation``.  Plate conventions: a strip edge, not an
    in-plane ring.  ``twist_lsq`` and ``twist_edge`` use DOF 2 (``w``) against the ``y``
    coordinate; ``twist_rotation`` averages DOF 3 (``rx``).  They cannot see in-plane shear at
    all.

``src/aeroelast/postprocess/sectional.py``
    No field-to-rotation estimator.  ``_section_axes`` measures axes from geometry (like
    ``_measured_tip_axes``) and ``_projected_rots`` *prescribes* a rotation; nothing there
    recovers a section rotation from a displacement field.

Every ``print`` below is one case line: ``name | prescribed | recovered | error | distortion``.
The assertions are properties of our own instrumentation (plain ``assert``), not store
comparisons.  A failing case is the finding: the estimator is not adjusted and no bound is
weakened to make it pass.
"""

from __future__ import annotations

import importlib

import numpy as np
import pytest

from tests.validation.blade.test_blade_rated_twist import (
    _measured_tip_axes,
    _ring_kinematics,
)
from tests.validation.blade.test_blade_twist_anchor_beam import (
    _ring_section_rotation as _ring_section_rotation_anchor,
)
from tests.validation.blade.test_laminate_bend_twist import (
    twist_edge,
    twist_lsq,
    twist_rotation,
)
from tests.validation.parity.test_thin_walled_tube_torsion import (
    B,
    H,
    _midline_ring,
    _parallelogram,
    _theta_fit,
    _theta_z,
)
from tests.validation.rotor.test_bem_fsi_deformed_geometry import (
    _ring_section_rotation as _ring_section_rotation_bem,
)

#: Round-off scale for a unit-geometry, few-dozen-node analytic field.  Used absolutely for the
#: zero-rotation cases (never divided by a value that is itself zero) and as the *relative*
#: bound for the exact rigid-rotation cases.
ROUNDOFF = 1e-12

RAD = np.deg2rad


# --------------------------------------------------------------------------------------
# Estimator adapters.  Each returns a dict with a "rotation" (radians) and, when the
# estimator exposes one, its distortion metrics.
# --------------------------------------------------------------------------------------


def _adapt_ring_kinematics(coords, u, ring):
    d = _ring_kinematics(coords, u, np.asarray(ring))
    return {
        "rotation": d["omega"],
        "shear": d["shear"],
        "dilatation": d["dilatation"],
        "distortion": d["distortion"],
        "rigid_rotation": d["rigid_rotation"],
    }


def _adapt_theta_fit(coords, u, ring):
    return {"rotation": _theta_fit(coords, u, list(ring))}


def _adapt_theta_z(coords, u, ring):
    return {"rotation": _theta_z(u, list(ring))}


def _adapt_parallelogram(coords, u, ring):
    alpha, beta, deviation = _parallelogram(coords, u, list(ring))
    return {
        "rotation": 0.5 * (beta - alpha),
        "alpha": alpha,
        "beta": beta,
        "distortion": deviation,
    }


def _adapt_affine_u(coords, u, ring):
    mod = importlib.import_module("tests.validation.blade.test_blade_section_torsion_stiffness")
    return {"rotation": mod._ring_section_rotation(coords, u, np.asarray(ring))}


def _adapt_affine_disp(fn):
    def _call(coords, u, ring):
        disp = u.reshape(-1, 6)[:, :3]
        return {"rotation": fn(coords, disp, np.asarray(ring))}

    return _call


class _Estimator:
    def __init__(self, name, path_line, fn, affine=False):
        self.name = name
        self.path_line = path_line
        self.fn = fn
        self.affine = affine


RING_ESTIMATORS = [
    _Estimator(
        "_ring_kinematics",
        "tests/validation/blade/test_blade_rated_twist.py:459",
        _adapt_ring_kinematics,
    ),
    _Estimator(
        "_theta_fit",
        "tests/validation/parity/test_thin_walled_tube_torsion.py:344",
        _adapt_theta_fit,
    ),
    _Estimator(
        "_theta_z",
        "tests/validation/parity/test_thin_walled_tube_torsion.py:355",
        _adapt_theta_z,
    ),
    _Estimator(
        "_parallelogram",
        "tests/validation/parity/test_thin_walled_tube_torsion.py:360",
        _adapt_parallelogram,
    ),
    _Estimator(
        "_ring_section_rotation (u)",
        "tests/validation/blade/test_blade_section_torsion_stiffness.py:145",
        _adapt_affine_u,
        affine=True,
    ),
    _Estimator(
        "_ring_section_rotation (disp) [anchor]",
        "tests/validation/blade/test_blade_twist_anchor_beam.py:85",
        _adapt_affine_disp(_ring_section_rotation_anchor),
        affine=True,
    ),
    _Estimator(
        "_ring_section_rotation (disp) [bem_fsi]",
        "tests/validation/rotor/test_bem_fsi_deformed_geometry.py:149",
        _adapt_affine_disp(_ring_section_rotation_bem),
        affine=True,
    ),
]


# --------------------------------------------------------------------------------------
# Analytic ring / field builders.  No mesh generation, no solver: a few dozen nodes.
# --------------------------------------------------------------------------------------


def _circle_ring(n=24, radius=1.0, centre=(0.0, 0.0), z=0.0, warp=0.0):
    """A circular ring in the xy plane.  ``warp`` breaks planarity with a 1e-3-style offset."""
    phi = 2.0 * np.pi * np.arange(n) / n
    x = centre[0] + radius * np.cos(phi)
    y = centre[1] + radius * np.sin(phi)
    zz = np.full(n, z) + (warp * np.cos(phi) if warp else 0.0)
    return np.column_stack([x, y, zz])


def _tube_rect_ring():
    """The tube test's real mid-line rectangle (B x H), at z = 0."""
    return np.asarray([[x, y, 0.0] for x, y in _midline_ring()], dtype=float)


def _inplane_u(coords, ring, omega, g=0.0, centre=(0.0, 0.0), tx=0.0, ty=0.0):
    """The analytic field ``ux = (g - omega)(y - cy) + tx``, ``uy = (g + omega)(x - cx) + ty``.

    Its antisymmetric part is exactly ``omega`` (rigid rotation) and its symmetric part is the
    pure shear ``g`` (a parallelogram that stretches one diagonal and compresses the other).  The
    nodal rotation DOF 5 is set to ``omega`` so the DOF-reading estimators see the same rigid
    rotation.
    """
    u = np.zeros(6 * len(coords))
    idx = np.asarray(ring)
    dx = coords[idx, 0] - centre[0]
    dy = coords[idx, 1] - centre[1]
    u[6 * idx] = (g - omega) * dy + tx
    u[6 * idx + 1] = (g + omega) * dx + ty
    u[6 * idx + 5] = omega
    return u


def _report(label, est, prescribed, out):
    recovered = out["rotation"]
    err = abs(recovered - prescribed) / abs(prescribed) if prescribed else abs(recovered)
    distortion = out.get("distortion")
    dist_str = f"{distortion:.6e}" if distortion is not None else "n/a"
    print(
        f"  {label:34s} | {est.name:40s} | prescribed={prescribed:+.9e} "
        f"recovered={recovered:+.9e} error={err:.3e} distortion={dist_str}"
    )
    return err


# --------------------------------------------------------------------------------------
# 1. Rigid rotation only
# --------------------------------------------------------------------------------------


@pytest.mark.parametrize("angle_deg", [0.5, -0.5, 5.0, -5.0])
def test_rigid_rotation_only(angle_deg):
    """A ring rigidly rotated about its own centre: every estimator must return that angle."""
    omega = RAD(angle_deg)
    coords = _circle_ring(centre=(0.0, 0.0))
    ring = np.arange(len(coords))
    u = _inplane_u(coords, ring, omega, g=0.0, centre=(0.0, 0.0), tx=1e-3, ty=-2e-3)
    print(f"\nRIGID ROTATION ONLY: prescribed omega = {angle_deg:+.4f} deg = {omega:+.9e} rad")
    for est in RING_ESTIMATORS:
        out = est.fn(coords, u, ring)
        err = _report("rigid", est, omega, out)
        assert err < ROUNDOFF, (
            f"{est.name} ({est.path_line}) recovered {out['rotation']!r} for a rigid rotation "
            f"of {omega!r} rad; relative error {err:.3e} exceeds {ROUNDOFF:.0e}"
        )
        if "distortion" in out:
            assert abs(out["distortion"]) < ROUNDOFF, (
                f"{est.name} reported distortion {out['distortion']!r} for a pure rigid rotation"
            )
        if "shear" in out:
            assert abs(out["shear"]) < ROUNDOFF
            assert abs(out["dilatation"]) < ROUNDOFF


# --------------------------------------------------------------------------------------
# 2. Pure distortion only, zero net rotation
# --------------------------------------------------------------------------------------


@pytest.mark.parametrize("g", [1e-3, -1e-3])
def test_pure_distortion_only_zero_rotation(g):
    """A symmetric parallelogram shear with zero net rotation: no estimator may invent a rotation."""
    coords = _circle_ring(centre=(0.0, 0.0))
    ring = np.arange(len(coords))
    u = _inplane_u(coords, ring, omega=0.0, g=g, centre=(0.0, 0.0))
    print(f"\nPURE DISTORTION ONLY (circle): prescribed omega = 0, shear g = {g:+.6e}")
    for est in RING_ESTIMATORS:
        out = est.fn(coords, u, ring)
        print(
            f"  {'pure-distortion':34s} | {est.name:40s} | prescribed={0.0:+.9e} "
            f"recovered={out['rotation']:+.9e} abs_error={abs(out['rotation']):.3e} "
            f"distortion={out.get('distortion', float('nan')):.6e}"
        )
        assert abs(out["rotation"]) < ROUNDOFF, (
            f"{est.name} ({est.path_line}) invented a rotation of {out['rotation']!r} rad from a "
            f"pure symmetric distortion with zero net rotation"
        )
        if est.name == "_ring_kinematics":
            assert abs(out["shear"] - g) / abs(g) < ROUNDOFF, (
                f"_ring_kinematics shear {out['shear']!r} does not recover the prescribed {g!r}"
            )
            expected_distortion = np.sqrt(2.0) * abs(g)
            assert abs(out["distortion"] - expected_distortion) / expected_distortion < ROUNDOFF
        if est.name == "_parallelogram":
            assert abs(out["alpha"] - g) / abs(g) < ROUNDOFF
            assert abs(out["beta"] - g) / abs(g) < ROUNDOFF


@pytest.mark.xfail(
    strict=True,
    reason=(
        "measured instrument limit: on the tube's own B x H rectangle a pure shear g = 1e-3 "
        "reads as omega = 4.7059e-4 rad where the true rotation is zero. Strict, so fixing the "
        "estimator turns this into an XPASS and the marker must go."
    ),
)
def test_pure_distortion_on_the_real_tube_section():
    """The tube's own ``B x H`` rectangle: ``_theta_fit`` must not read shear as a rotation.

    This is the case that matters for the tube inspector, because the tube is *not* circular.
    """
    coords = _tube_rect_ring()
    ring = np.arange(len(coords))
    g = 1e-3
    u = _inplane_u(coords, ring, omega=0.0, g=g, centre=(0.0, 0.0))
    print(
        f"\nPURE DISTORTION ONLY (tube rectangle B={B} H={H}): prescribed omega = 0, "
        f"shear g = {g:+.6e}"
    )
    failures = []
    for est in RING_ESTIMATORS:
        if est.name == "_theta_z":
            continue  # reads rz only; covered by the circle case
        out = est.fn(coords, u, ring)
        print(
            f"  {'pure-distortion (rect)':34s} | {est.name:40s} | prescribed={0.0:+.9e} "
            f"recovered={out['rotation']:+.9e} abs_error={abs(out['rotation']):.3e} "
            f"distortion={out.get('distortion', float('nan')):.6e}"
        )
        if abs(out["rotation"]) >= ROUNDOFF:
            failures.append((est.name, est.path_line, float(out["rotation"])))
    assert not failures, (
        "these estimators invented a rotation from a pure symmetric distortion (zero net "
        "rotation) on the tube's own rectangular section: "
        + "; ".join(f"{n} ({p}) -> {r!r} rad" for n, p, r in failures)
    )


# --------------------------------------------------------------------------------------
# 3. Rigid rotation + distortion
# --------------------------------------------------------------------------------------


@pytest.mark.parametrize("angle_deg", [5.0, -5.0])
def test_rigid_plus_distortion(angle_deg):
    """Adding distortion must not move the recovered rigid-rotation part."""
    omega = RAD(angle_deg)
    g = 1e-3
    coords = _circle_ring(centre=(0.0, 0.0))
    ring = np.arange(len(coords))
    u = _inplane_u(coords, ring, omega=omega, g=g, centre=(0.0, 0.0))
    print(f"\nRIGID + DISTORTION: prescribed omega = {omega:+.9e} rad, shear g = {g:+.6e}")
    for est in RING_ESTIMATORS:
        out = est.fn(coords, u, ring)
        err = _report("rigid+distortion", est, omega, out)
        assert err < ROUNDOFF, (
            f"{est.name} ({est.path_line}) recovered {out['rotation']!r} for a rigid rotation of "
            f"{omega!r} rad combined with distortion; relative error {err:.3e}"
        )


# --------------------------------------------------------------------------------------
# 4. Twisted ring: a rigid rotation whose angle varies linearly along the span
# --------------------------------------------------------------------------------------


def test_twisted_ring_local_angles():
    """Several rings, each rigidly rotated by its own angle; each must recover its own angle."""
    thetas = RAD(np.array([0.0, 1.0, -2.5, 4.0]))
    rings, coords_parts = [], []
    off = 0
    for k in range(len(thetas)):
        part = _circle_ring(centre=(0.0, 0.0), z=float(k) * 0.5)
        rings.append(np.arange(off, off + len(part)))
        coords_parts.append(part)
        off += len(part)
    coords = np.vstack(coords_parts)
    u = np.zeros(6 * len(coords))
    for ring, theta in zip(rings, thetas, strict=True):
        dx = coords[ring, 0]
        dy = coords[ring, 1]
        u[6 * ring] = -theta * dy
        u[6 * ring + 1] = theta * dx
        u[6 * ring + 5] = theta
    print("\nTWISTED RING: prescribed omega varies linearly along the span")
    for k, (ring, theta) in enumerate(zip(rings, thetas, strict=True)):
        for est in RING_ESTIMATORS:
            out = est.fn(coords, u, ring)
            err = _report(f"twisted ring {k}", est, theta, out)
            assert err < ROUNDOFF, (
                f"{est.name} ({est.path_line}) recovered {out['rotation']!r} at twisted ring {k} "
                f"whose own angle is {theta!r} rad; relative error {err:.3e}"
            )


# --------------------------------------------------------------------------------------
# 5. Non-planar (prebent) ring: a small out-of-plane offset
# --------------------------------------------------------------------------------------


@pytest.mark.parametrize("angle_deg", [0.5, 5.0])
def test_prebent_nonplanar_ring(angle_deg):
    """A 1e-3 m out-of-plane warp must leave the in-plane rotation estimate unchanged."""
    omega = RAD(angle_deg)
    coords = _circle_ring(centre=(0.0, 0.0), warp=1e-3)
    ring = np.arange(len(coords))
    u = _inplane_u(coords, ring, omega=omega, g=0.0, centre=(0.0, 0.0))
    print(f"\nPREBENT / NON-PLANAR RING (warp = 1e-3 m): prescribed omega = {omega:+.9e} rad")
    for est in RING_ESTIMATORS:
        out = est.fn(coords, u, ring)
        err = _report("prebent", est, omega, out)
        assert err < ROUNDOFF, (
            f"{est.name} ({est.path_line}) recovered {out['rotation']!r} for a rigid rotation "
            f"{omega!r} with a 1e-3 m out-of-plane warp; relative error {err:.3e}"
        )


# --------------------------------------------------------------------------------------
# 6. Off-axis ring: origin dependence of the centroid estimators
# --------------------------------------------------------------------------------------


@pytest.mark.xfail(
    strict=True,
    reason=(
        "measured instrument limit: the centroid-based tube estimators are origin-dependent. "
        "A ring centred at (2.0, 1.5) and rotated rigidly by 3.491e-2 rad reads 4.815e-3 "
        "(_theta_fit, 86% error) and 5.113e-3 (_parallelogram, 85% error), while the affine "
        "estimators are exact. Strict, so a fix turns this into an XPASS."
    ),
)
def test_off_axis_ring_origin_dependence():
    """A ring centred away from the coordinate origin, rigidly rotated about its own centre.

    This is the blade's situation (its rings are not centred on the global origin).  The affine
    estimators and the DOF reader are origin-independent; the centroid estimators are not.
    """
    omega = RAD(2.0)
    centre = (2.0, 1.5)
    coords = _circle_ring(centre=centre, radius=1.0)
    ring = np.arange(len(coords))
    u = _inplane_u(coords, ring, omega=omega, g=0.0, centre=centre)
    print(f"\nOFF-AXIS RING (centre = {centre}): prescribed omega = {omega:+.9e} rad")
    failures = []
    for est in RING_ESTIMATORS:
        out = est.fn(coords, u, ring)
        err = _report("off-axis rigid", est, omega, out)
        if err >= ROUNDOFF:
            failures.append((est.name, est.path_line, float(out["rotation"]), err))
    assert not failures, (
        "these estimators are not invariant to the ring's location and mis-recover an off-axis "
        "rigid rotation about the ring's own centre: "
        + "; ".join(f"{n} ({p}) -> {r!r} rad (rel err {e:.3e})" for n, p, r, e in failures)
    )


# --------------------------------------------------------------------------------------
# 7. Plate-strip twist estimators
# --------------------------------------------------------------------------------------


def _strip(n_y=9, half=0.5, xs=(0.0, 1.0, 2.0, 3.0), warp=0.0):
    coords, stations = [], []
    for x in xs:
        ids = []
        for y in np.linspace(-half, half, n_y):
            z = warp * np.cos(np.pi * y / (2.0 * half)) if warp else 0.0
            ids.append(len(coords))
            coords.append([x, y, z])
        stations.append(np.asarray(ids))
    return np.asarray(coords, float), stations


def _strip_u(coords, stations, thetas, warp_amp=0.0):
    u = np.zeros(6 * len(coords))
    for st, theta in zip(stations, thetas, strict=True):
        y = coords[st, 1]
        u[6 * st + 2] = theta * y + warp_amp * y * y
        u[6 * st + 3] = theta
    return u


TWIST_ESTIMATORS = [
    (
        "twist_lsq",
        "tests/validation/blade/test_laminate_bend_twist.py:182",
        lambda coords, u, st: twist_lsq(coords, u, list(st)),
    ),
    (
        "twist_edge",
        "tests/validation/blade/test_laminate_bend_twist.py:189",
        lambda coords, u, st: twist_edge(coords, u, list(st)),
    ),
    (
        "twist_rotation",
        "tests/validation/blade/test_laminate_bend_twist.py:193",
        lambda coords, u, st: twist_rotation(u, list(st)),
    ),
]


def test_twist_estimators_rigid_and_distortion():
    """The plate strip: rigid twist is recovered; a zero-slope warping is not read as a twist.

    These estimators read only ``w`` (DOF 2) or ``rx`` (DOF 3), so the in-plane parallelogram
    mode is invisible to them.  The strongest zero-rotation case they accept is a symmetric
    (even) out-of-plane warping ``w = a * y**2``, whose linear slope is exactly zero.
    """
    thetas = RAD(np.array([0.5, -0.5, 5.0, -5.0]))
    coords, stations = _strip()
    print("\nPLATE STRIP - RIGID TWIST (each station's own angle)")
    u = _strip_u(coords, stations, thetas)
    for st, theta in zip(stations, thetas, strict=True):
        for name, path, fn in TWIST_ESTIMATORS:
            rec = fn(coords, u, st)
            print(
                f"  {'strip rigid twist':34s} | {name:40s} ({path}) | prescribed={theta:+.9e} "
                f"recovered={rec:+.9e} error={abs(rec - theta) / abs(theta):.3e}"
            )
            assert abs(rec - theta) / abs(theta) < ROUNDOFF

    print("\nPLATE STRIP - ZERO-NET-ROTATION WARPING (even mode, prescribed twist = 0)")
    a = 1e-3
    u = _strip_u(coords, stations, [0.0] * len(stations), warp_amp=a)
    for st in stations:
        for name, path, fn in TWIST_ESTIMATORS:
            rec = fn(coords, u, st)
            print(
                f"  {'strip even warping':34s} | {name:40s} | prescribed={0.0:+.9e} "
                f"recovered={rec:+.9e} abs_error={abs(rec):.3e}"
            )
            assert abs(rec) < ROUNDOFF, (
                f"{name} ({path}) invented a twist {rec!r} rad from a symmetric (even) "
                f"out-of-plane warping with zero net rotation"
            )


# --------------------------------------------------------------------------------------
# 8. Geometry-only helper: _measured_tip_axes
# --------------------------------------------------------------------------------------


def test_measured_tip_axes_recovers_known_ellipse_axes():
    """``_measured_tip_axes`` is geometry only; it must recover a known major-axis direction."""
    phi = np.linspace(0.0, 2.0 * np.pi, 32, endpoint=False)
    a, b = 2.0, 0.5
    for alpha_deg in (0.0, 30.0, -60.0):
        alpha = RAD(alpha_deg)
        x = a * np.cos(phi) * np.cos(alpha) - b * np.sin(phi) * np.sin(alpha)
        y = a * np.cos(phi) * np.sin(alpha) + b * np.sin(phi) * np.cos(alpha)
        coords = np.column_stack([x, y, np.zeros_like(x)])
        ring = np.arange(len(coords))
        chord, flap = _measured_tip_axes(coords, ring)
        major = np.array([np.cos(alpha), np.sin(alpha), 0.0])
        alignment = abs(float(np.dot(chord, major)))
        print(
            f"\nMEASURED TIP AXES: major axis at {alpha_deg:+.1f} deg, |chord . major| = "
            f"{alignment:.12f}, chord . flap = {float(np.dot(chord, flap)):.3e}"
        )
        assert abs(alignment - 1.0) < ROUNDOFF
        assert abs(float(np.dot(chord, flap))) < ROUNDOFF


def test_no_section_rotation_estimator_in_sectional_module():
    """Document that ``sectional.py`` has no field-to-rotation estimator (inventory guard)."""
    import aeroelast.postprocess.sectional as sectional

    source_has_recovery = any(
        name in vars(sectional)
        for name in ("_section_rotation", "section_rotation", "omega_from_field")
    )
    assert not source_has_recovery, (
        "sectional.py gained a field-to-rotation estimator; add it to the instrument inventory "
        "and extend this test file"
    )
    assert hasattr(sectional.SectionalExtractor, "_projected_rots"), (
        "sectional.py's _projected_rots (a rotation *prescriber*, not an estimator) moved; "
        "update the inventory docstring"
    )
