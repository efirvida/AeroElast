"""Are quads (MITC4) and triangles (MITC3) consistent for composites? Issue #34, element level.

Run from a checkout with ``PYTHONPATH=$PWD/src`` (and not ``python -I``)::

    scripts/aeroenv.sh python tools/diagnose_mixed_mesh_consistency.py --element
    scripts/aeroenv.sh python tools/diagnose_mixed_mesh_consistency.py --blade

``--element`` (E1): a cantilever composite strip along Y, meshed all-quad, all-triangle and with
alternating rows. A quad takes its first edge (the span) as laminate axis 1; a triangle takes global
X projected on its plane (``mitc3.rs::rotate_constitutive_to_local``), i.e. the *transverse* axis.
Torsion cannot show this (``Qbar66`` is invariant under a 90 degree swap), so the metrics are the
tip axial, bending and bend-twist responses. The CLT prediction of the triangle/quad ratio, for a
laminate rotated by 90 degrees, is printed BEFORE anything is solved.

``--blade`` (E2): triangles of the production blade mesh by span station, their area share and the
angle between their effective axis 1 and the span-projected axis the quads follow.

The blade is meshed at 0.25 m. 0.5 m is NOT converged for the blade (torsion and twist move
between 0.5 and 0.25) and must not be used for any result here; it only appears in the
already-done mesh-convergence studies.

Read-only diagnostic: no element code is touched and nothing here is a pass/fail test.
"""

import sys
from pathlib import Path

import numpy as np

_REPO_ROOT = Path(__file__).resolve().parents[1]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

# ---- strip definition -------------------------------------------------------
L, W = 1.0, 0.1  # length along Y, width along X [m]
NXC, NY = 4, 20  # elements across the width, along the span
E1, E2, G12, NU12 = 120e9, 10e9, 5e9, 0.3  # same CFRP as the tube test
PLY_T = 0.125e-3
P = 10.0  # tip load [N]
TOL = 0.05  # the suite's 5 % rule
LAMINATES = {"UD [0]4": [0.0] * 4, "[45,0]s": [45.0, 0.0, 0.0, 45.0]}


def _laminate(angles):
    from aeroelast.core.laminate import create_laminate_from_angles
    from aeroelast.core.material import OrthotropicMaterial

    mat = OrthotropicMaterial(
        "CFRP", E=(E1, E2, E2), G=(G12, G12, G12), nu=(NU12, NU12, NU12), rho=1500.0
    )
    return create_laminate_from_angles(mat, PLY_T, angles)


def _comp_prop(lam):
    h = lam.total_thickness
    return {
        "type": "composite",
        "cm": lam.A.ravel().tolist(),
        "b_coupling": lam.B.ravel().tolist(),
        "cb": lam.D.ravel().tolist(),
        "cs": lam.Cs.ravel().tolist(),
        "cs_uncorrected": lam.shear_stiffness_uncorrected().ravel().tolist(),
        "thickness": h,
        "e_equiv": float(np.trace(lam.A) / (3.0 * h)),
        "mass_per_area": 0.0,
        "rotational_inertia": 0.0,
    }


def _iso_prop():
    return {"type": "isotropic", "e": 70e9, "nu": 0.3, "rho": 2700.0,
            "thickness": 0.5e-3, "shear_correction": 5.0 / 6.0}


def _grid():
    nid = {}
    coords = []
    for i in range(NXC + 1):
        for j in range(NY + 1):
            nid[(i, j)] = len(coords)
            coords.append([W * i / NXC - W / 2.0, L * j / NY, 0.0])
    return np.asarray(coords), nid


def _mesh(pattern, quad_code, tri_code):
    """Connectivity and element codes. First edge of every element runs along +Y (the span)."""
    coords, nid = _grid()
    conn, codes = [], []
    for i in range(NXC):
        for j in range(NY):
            a, b = nid[(i + 1, j)], nid[(i + 1, j + 1)]
            c, d = nid[(i, j + 1)], nid[(i, j)]
            tri = pattern == "tri" or (pattern == "rows" and j % 2 == 0)
            if tri:
                conn += [[a, b, c], [a, c, d]]
                codes += [tri_code, tri_code]
            else:
                conn.append([a, b, c, d])
                codes.append(quad_code)
    return coords, conn, codes, nid


def _solve(pattern, prop, quad_code, tri_code, load_dir):
    from _aeroelast import PyMeshAssembler
    from scipy.sparse import coo_matrix
    from scipy.sparse.linalg import spsolve

    coords, conn, codes, nid = _mesh(pattern, quad_code, tri_code)
    asm = PyMeshAssembler(coords, conn, codes, [prop] * len(conn))
    r, c, v = asm.assemble_k()
    n = asm.dofs_count
    K = coo_matrix((v, (r, c)), shape=(n, n)).tocsr()
    tip = [nid[(i, NY)] for i in range(NXC + 1)]
    f = np.zeros(n)
    for t in tip:
        f[6 * t + load_dir] = P / len(tip)
    clamped = {6 * nid[(i, 0)] + d for i in range(NXC + 1) for d in range(6)}
    free = np.array(sorted(set(range(n)) - clamped))
    u = np.zeros(n)
    u[free] = spsolve(K[np.ix_(free, free)].tocsc(), f[free])
    ut = np.array([u[6 * t:6 * t + 3] for t in tip])
    x = coords[tip, 0] - coords[tip, 0].mean()
    twist = np.rad2deg(-np.sum(x * ut[:, 2]) / np.sum(x * x))  # rotation about +Y
    return ut[:, 1].mean(), ut[:, 2].mean(), twist


def _clt(lam):
    """Closed forms in the laminate frame (axis 1 = strip span; free edges)."""
    a = np.linalg.inv(lam.A)
    d = np.linalg.inv(lam.D)
    axial = a[0, 0] * (P / W) * L
    bend = d[0, 0] * (P / W) * L**3 / 3.0
    twist_rate = d[2, 0] * (P / W) * L**2 / 2.0  # same convention as tests/test_laminate_bend_twist
    return axial, bend, twist_rate


def element_study():
    print("E1: composite strip, quad vs triangle laminate axis (issue #34)\n")
    print("--- PRE-REGISTERED (CLT only, nothing solved yet) ---")
    pred = {}
    for name, ang in LAMINATES.items():
        q = _clt(_laminate(ang))
        t = _clt(_laminate([a + 90.0 for a in ang]))  # triangle axis 1 = transverse
        pred[name] = (q, t)
        print(f"{name:>9}: closed form quad  axial {q[0]:.4e} m  bend {q[1]:.4e} m  "
              f"twist-coupling {q[2]:.3e}")
        print(f"{'':>9}  predicted tri/quad   axial {t[0] / q[0]:.3f}  bend {t[1] / q[1]:.3f}  "
              f"twist {t[2] / q[2]:.3f}   rows(axial/bend) {(1 + t[0] / q[0]) / 2:.3f} / "
              f"{(1 + t[1] / q[1]) / 2:.3f}")
    print("isotropic control: predicted tri/quad 1.000 (tolerance 1 %)\n")
    print("--- MEASURED ---")
    print(f"{'case':>9} {'mesh':>5} {'uy axial':>11} {'uz bend':>11} {'twist[deg]':>11}"
          f"  {'tri/quad: axial':>15} {'bend':>7} {'twist':>7}")
    cases = [(nm, _comp_prop(_laminate(a)), 44, 33) for nm, a in LAMINATES.items()]
    cases.append(("isotropic", _iso_prop(), 4, 3))
    for name, prop, qc, tc in cases:
        base = None
        for pattern in ("quad", "tri", "rows"):
            ax = _solve(pattern, prop, qc, tc, 1)[0]
            _, bz, _ = _solve(pattern, prop, qc, tc, 2)
            tw = _solve(pattern, prop, qc, tc, 2)[2]
            if pattern == "quad":
                base = (ax, bz, tw)
            ratio = (ax / base[0], bz / base[1], tw / base[2] if base[2] else float("nan"))
            print(f"{name:>9} {pattern:>5} {ax:11.4e} {bz:11.4e} {tw:11.4e}  "
                  f"{ratio[0]:15.3f} {ratio[1]:7.3f} {ratio[2]:7.3f}")
        if name in pred:
            q = pred[name][0]
            print(f"{'':>9} closed form quad: axial {q[0]:.4e} (measured/closed "
                  f"{base[0] / q[0]:.3f}), bend {q[1]:.4e} (measured/closed {base[1] / q[1]:.3f})")


def blade_study():
    from aeroelast.models.blade.model import Blade
    import tests.validation.blade.test_blade_rated_twist as t

    model = Blade(str(t.YAML), element_size=0.25)
    model.generate_mesh()
    mesh = model.mesh
    xyz = np.array([[n.x, n.y, n.z] for n in mesh.nodes])
    idx = mesh.node_id_to_index
    R = xyz[:, 2].max()
    rows = []
    n_quad = n_tri = 0
    for el in mesh.elements:
        ids = [idx[i] for i in el.node_ids]
        p = xyz[ids]
        if len(ids) == 3:
            n_tri += 1
        elif len(ids) == 4:
            n_quad += 1
        else:
            continue
        nrm = np.cross(p[1] - p[0], p[2] - p[0]) if len(ids) == 3 else np.cross(
            p[2] - p[0], p[3] - p[1])
        area = 0.5 * np.linalg.norm(nrm)
        nrm = nrm / np.linalg.norm(nrm)
        xg = np.array([1.0, 0.0, 0.0])
        xp = xg - xg.dot(nrm) * nrm
        if np.linalg.norm(xp) > 1e-10:
            ref = xp / np.linalg.norm(xp)
        else:
            zp = np.array([0.0, 0.0, 1.0]) - nrm[2] * nrm
            ref = zp / np.linalg.norm(zp)
        sp = np.array([0.0, 0.0, 1.0]) - nrm[2] * nrm
        sn = np.linalg.norm(sp)
        err = np.degrees(np.arccos(min(1.0, abs(ref @ sp / sn)))) if sn > 1e-10 else np.nan
        rows.append((len(ids), p[:, 2].mean() / R, area, err))
    a = np.array(rows)
    tri = a[a[:, 0] == 3]
    print(f"E2: blade mesh {len(mesh.nodes)} nodes, {n_quad} quads, {n_tri} triangles "
          f"({100 * n_tri / (n_tri + n_quad):.2f} % of elements, "
          f"{100 * tri[:, 2].sum() / a[:, 2].sum():.2f} % of area)")
    bins = np.linspace(0, 1, 11)
    print("r/R bin      tri  tri area[m2]  tri-area share  area-wtd |axis err| [deg]  err>45deg share")
    for lo, hi in zip(bins[:-1], bins[1:]):
        m = (a[:, 1] >= lo) & (a[:, 1] < hi + (1e-9 if hi == 1 else 0))
        tm = m & (a[:, 0] == 3)
        if not tm.any():
            continue
        ar = a[tm, 2]
        e = np.nan_to_num(a[tm, 3])
        print(f"{lo:.1f}-{hi:.1f}  {tm.sum():7d}  {ar.sum():11.3f}  "
              f"{100 * ar.sum() / a[m, 2].sum():13.2f} %  {np.average(e, weights=ar):22.1f}"
              f"  {100 * ar[e > 45].sum() / ar.sum():14.1f} %")
    mid = (tri[:, 1] >= 0.3) & (tri[:, 1] <= 0.9)
    print(f"\ntriangle area inside 0.3-0.9 R: {tri[mid, 2].sum():.3f} m2 "
          f"({100 * tri[mid, 2].sum() / a[(a[:, 1] >= .3) & (a[:, 1] <= .9), 2].sum():.2f} % of the "
          f"area there)")
    sets = getattr(mesh, "element_sets", None)
    if isinstance(sets, dict):
        tri_ids = {el.id for el in mesh.elements if len(el.node_ids) == 3}
        by_set = {}
        for name, es in sets.items():
            members = getattr(es, "elements", es)
            ids = set(members.keys()) if isinstance(members, dict) else {
                getattr(m, "id", m) for m in members}
            k = len(ids & tri_ids)
            if k:
                by_set[name] = (k, len(ids))
        print(f"\nelement sets holding triangles: {len(by_set)} of {len(sets)}")
        for name, (k, tot) in sorted(by_set.items(), key=lambda kv: -kv[1][0] / kv[1][1])[:20]:
            print(f"  {name:<28} {k:3d} tri / {tot:4d} elements")
        comp = {}
        for name, (k, tot) in by_set.items():
            comp[name.split("_", 2)[-1] if name[:2].isdigit() else name] = comp.get(
                name.split("_", 2)[-1] if name[:2].isdigit() else name, 0) + k
        print("triangles by component suffix:", dict(sorted(comp.items(), key=lambda kv: -kv[1])))


# ---- E3: triangles moved onto the quad contract (axis 1 = projected span) ---------------
E3_THRESHOLD_DEG = 0.42 * 3.6  # the band width used in T7d: Zhou to BeamDyn, relative to Zhou


def psi_deg(p, span):
    """Signed angle [deg] about the element normal from projected global X to the projected span.

    A triangle's laminate axis 1 is projected global X (global Z when X is normal). Shifting every
    ply angle by this value points the fibres where the quad contract (axis 1 = projected span)
    would.  e3 follows the node order, ``(p1 - p0) x (p2 - p0)``.
    """
    e3 = np.cross(p[1] - p[0], p[2] - p[0])
    e3 = e3 / np.linalg.norm(e3)
    xg = np.array([1.0, 0.0, 0.0])
    ref = xg - xg.dot(e3) * e3
    if np.linalg.norm(ref) <= 1e-10:
        ref = np.array([0.0, 0.0, 1.0]) - e3[2] * e3
    ref = ref / np.linalg.norm(ref)
    sp = np.asarray(span, float) - np.dot(span, e3) * e3
    sp = sp / np.linalg.norm(sp)
    return float(np.degrees(np.arctan2(e3 @ np.cross(ref, sp), ref @ sp)))


def calibrate(alpha_deg=30.0):
    """Check the sign and the method of the shift on a strip turned by alpha about Z.

    At 90 deg the sign is invisible (a ply is symmetric under 180 deg), so the strip is turned: the
    quad strip (axis 1 = span) is the reference, the triangle strip must match it only with +psi.
    """
    from _aeroelast import PyMeshAssembler
    from scipy.sparse import coo_matrix
    from scipy.sparse.linalg import spsolve

    a = np.radians(alpha_deg)
    rot = np.array([[np.cos(a), -np.sin(a), 0.0], [np.sin(a), np.cos(a), 0.0], [0.0, 0.0, 1.0]])
    span = rot @ np.array([0.0, 1.0, 0.0])
    ang = LAMINATES["[45,0]s"]
    print(f"calibration: strip turned {alpha_deg:.0f} deg about Z, laminate [45,0]s")
    out = {}
    for label, pattern, sign in (("quad", "quad", 0.0), ("tri unshifted", "tri", 0.0),
                                 ("tri shifted +psi", "tri", 1.0), ("tri shifted -psi", "tri", -1.0)):
        coords, conn, codes, nid = _mesh(pattern, 44, 33)
        coords = coords @ rot.T
        psi = 0.0
        if pattern == "tri":
            psi = psi_deg(coords[conn[0]], span)
        prop = _comp_prop(_laminate([x + sign * psi for x in ang]))
        asm = PyMeshAssembler(coords, conn, codes, [prop] * len(conn))
        r, c, v = asm.assemble_k()
        n = asm.dofs_count
        K = coo_matrix((v, (r, c)), shape=(n, n)).tocsr()
        tip = [nid[(i, NY)] for i in range(NXC + 1)]
        clamped = {6 * nid[(i, 0)] + d for i in range(NXC + 1) for d in range(6)}
        free = np.array(sorted(set(range(n)) - clamped))
        res = []
        for vec in (span, np.array([0.0, 0.0, 1.0])):
            f = np.zeros(n)
            for t in tip:
                f[6 * t:6 * t + 3] = P / len(tip) * vec
            u = np.zeros(n)
            u[free] = spsolve(K[np.ix_(free, free)].tocsc(), f[free])
            res.append(np.mean([u[6 * t:6 * t + 3] @ vec for t in tip]))
        out[label] = res
        print(f"  {label:>17}  psi = {psi:+8.2f} deg   axial {res[0]:+.4e}   bend {res[1]:+.4e}")
    q = out["quad"]
    for label in ("tri unshifted", "tri shifted +psi", "tri shifted -psi"):
        t = out[label]
        print(f"  {label:>17}  / quad:  axial {t[0] / q[0]:.3f}   bend {t[1] / q[1]:.3f}")


if __name__ == "__main__":
    if "--calibrate" in sys.argv:
        calibrate()
    elif "--element" in sys.argv:
        element_study()
    elif "--blade" in sys.argv:
        blade_study()
    else:
        print(__doc__)
