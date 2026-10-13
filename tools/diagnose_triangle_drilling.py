"""Does the composite triangle's drilling penalty (#38) change a physical response? (Line B.)

Gates were committed beforehand (42c645c, odd/tasks/shell-twist-mixed-mesh.md "Line B
pre-registration"):

* B0  a flat composite triangle under a rigid rotation about its normal stores no energy
      (ratio to a uniform stretch <= 1e-10, at thickness h and 2h).
* B1  self-equilibrated torsion of the closed [45,-45]s tube of
      tests/validation/parity/test_thin_walled_tube_torsion.py against its hand-written Bredt
      rate, on quad / triangle / mixed meshes, with the triangle drilling stiffness scaled
      through ``e_equiv`` by s = 0.1, 1, 10, 100.
      Gate 1: triangle and mixed rate / Bredt within 5 % at s = 1.
      Gate 2: triangle ratio changes <= 1 % between s = 0.1 and s = 10 (>5 % = #38 matters).

Run: ``scripts/aeroenv.sh env PYTHONPATH=$PWD/src python tools/diagnose_triangle_drilling.py``.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from _aeroelast import PyMeshAssembler
from scipy.sparse import coo_matrix

_REPO_ROOT = Path(__file__).resolve().parents[1]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

import tests.validation.parity.test_thin_walled_tube_torsion as tw  # noqa: E402
from aeroelast.core.laminate import create_laminate_from_angles  # noqa: E402
from aeroelast.core.material import OrthotropicMaterial  # noqa: E402

QUAD, TRI = 44, 33
B0_GATE = 1e-10
GATE1 = 0.05
GATE2_QUIET, GATE2_LOUD = 0.01, 0.05
SCALES = (0.1, 1.0, 10.0, 100.0)


def _prop(ply_t: float = tw.PLY_T, e_scale: float = 1.0) -> dict:
    mat = OrthotropicMaterial(
        "CFRP", E=(tw.E1, tw.E2, tw.E2), G=(tw.G12, tw.G12, tw.G12), nu=(tw.NU12, tw.NU12, tw.NU12), rho=tw.RHO_PLY
    )
    lam = create_laminate_from_angles(mat, ply_t, [45.0, -45.0, -45.0, 45.0])
    h = lam.total_thickness
    return {
        "type": "composite",
        "cm": lam.A.ravel().tolist(),
        "b_coupling": lam.B.ravel().tolist(),
        "cb": lam.D.ravel().tolist(),
        "cs": lam.Cs.ravel().tolist(),
        "cs_uncorrected": lam.shear_stiffness_uncorrected().ravel().tolist(),
        "thickness": h,
        "e_equiv": float(np.trace(lam.A) / (3.0 * h)) * e_scale,
        "mass_per_area": sum(p.material.rho * p.thickness for p in lam.plies),
        "rotational_inertia": sum(
            p.material.rho * (p.z_top**3 - p.z_bottom**3) / 3.0 for p in lam.plies
        ),
    }


def _k(coords, conn, codes, prop):
    asm = PyMeshAssembler(np.asarray(coords, float), conn, codes, [prop] * len(conn))
    r, c, v = asm.assemble_k()
    n = asm.dofs_count
    return coo_matrix((v, (r, c)), shape=(n, n)).tocsr()


def b0() -> bool:
    print("B0  rigid rotation about the normal, flat composite triangle")
    xyz = np.array([[0.0, 0.0, 0.0], [0.1, 0.0, 0.0], [0.0, 0.1, 0.0]])
    eps = 1e-3
    ok = True
    for label, ply_t in (("h", tw.PLY_T), ("2h", 2 * tw.PLY_T)):
        K = _k(xyz, [[0, 1, 2]], [TRI], _prop(ply_t))
        u_rig, u_str = np.zeros(18), np.zeros(18)
        for i, (x, y, _) in enumerate(xyz):
            u_rig[6 * i : 6 * i + 2] = (-eps * y, eps * x)
            u_rig[6 * i + 5] = eps
            u_str[6 * i] = eps * x
        ratio = float(u_rig @ K @ u_rig) / float(u_str @ K @ u_str)
        print(f"    thickness {label:>2}: u_rig K u_rig / u_str K u_str = {ratio:.3e}")
        ok &= abs(ratio) <= B0_GATE
    print(f"    GATE B0 (<= {B0_GATE:g}): {'PASS' if ok else 'FAIL'}")
    return ok


def _tube(pattern: str):
    """Tube connectivity and element codes; ``pattern`` in quad | tri | mixed."""
    coords, quads, rings, n_ring = tw._tube_mesh()
    conn, codes = [], []
    for e, q in enumerate(quads):
        k = e // n_ring  # z-row
        if pattern == "quad" or (pattern == "mixed" and k % 2 == 1):
            conn.append(q)
            codes.append(QUAD)
        else:
            conn += [[q[0], q[1], q[2]], [q[0], q[2], q[3]]]
            codes += [TRI, TRI]
    return coords, conn, codes, rings


def _ratio(pattern: str, s: float, ref: float) -> float:
    coords, conn, codes, rings = _tube(pattern)
    # e_equiv is read only by the composite triangle's k_drill; the quad does not take it.
    K = _k(coords, conn, codes, _prop(e_scale=s))
    f = tw._ring_shear_flow(coords, rings[-1], tw.TORQUE) + tw._ring_shear_flow(
        coords, rings[0], -tw.TORQUE
    )
    modes = tw._rigid_modes(coords)
    u = tw._solve_rigid_removed(K, f, modes)
    z = np.asarray([coords[r[0], 2] for r in rings])
    th = np.asarray([tw._theta_fit(coords, u, r) for r in rings])
    return tw._lsq_slope(z, th, 0.4 * tw.L, 0.9 * tw.L) / ref


def control_scaling_reaches_k_drill() -> bool:
    """e_equiv is the only knob B1 turns; check it really scales the triangle's theta_z stiffness."""
    xyz = np.array([[0.0, 0.0, 0.0], [0.1, 0.0, 0.0], [0.0, 0.1, 0.0]])
    d1 = _k(xyz, [[0, 1, 2]], [TRI], _prop(e_scale=1.0)).diagonal()[5::6]
    d10 = _k(xyz, [[0, 1, 2]], [TRI], _prop(e_scale=10.0)).diagonal()[5::6]
    ratio = d10 / d1
    ok = bool(np.allclose(ratio, 10.0, rtol=1e-6))
    print(f"CONTROL  triangle K[theta_z,theta_z] at s=10 / s=1 = {np.array2string(ratio, precision=6)}  (must be 10): {'OK' if ok else 'BROKEN'}")
    return ok


def b1() -> bool:
    _, _, _, ref = tw._bredt_laminate()
    print(f"\nB1  tube torsion, rate / Bredt (T/GJ = {ref:.6e} rad/m)")
    print(f"    {'mesh':>6} " + " ".join(f"s={s:<7g}" for s in SCALES))
    table = {}
    for pattern in ("quad", "tri", "mixed"):
        table[pattern] = [_ratio(pattern, s, ref) for s in SCALES]
        print(f"    {pattern:>6} " + " ".join(f"{r:9.6f}" for r in table[pattern]))
    i1 = SCALES.index(1.0)
    g1 = all(abs(table[p][i1] - 1.0) <= GATE1 for p in ("tri", "mixed"))
    tri = table["tri"]
    d2 = abs(tri[SCALES.index(10.0)] - tri[SCALES.index(0.1)]) / abs(tri[SCALES.index(0.1)])
    print(f"    GATE 1 (tri and mixed at s=1 within {GATE1:.0%} of Bredt): {'PASS' if g1 else 'FAIL'}")
    verdict = "invisible" if d2 <= GATE2_QUIET else ("MATTERS" if d2 > GATE2_LOUD else "in between")
    print(f"    GATE 2 triangle ratio change s=0.1 -> 10: {d2:.3%} (<= {GATE2_QUIET:.0%} quiet, > {GATE2_LOUD:.0%} loud): {verdict}")
    return g1 and d2 <= GATE2_QUIET


if __name__ == "__main__":
    a = b0()
    if not control_scaling_reaches_k_drill():
        raise SystemExit("the e_equiv knob does not scale k_drill; B1 would be void")
    b = b1()
    print(f"\nsummary: B0 {'PASS' if a else 'FAIL'}; B1 {'PASS' if b else 'see gates above'}")
