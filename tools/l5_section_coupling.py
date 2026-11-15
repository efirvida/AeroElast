"""L5 T9: which model is right about the flap-edge bending cross-coupling? (odd/tasks/l5-one-way.md L8, L9)

Three views of the 2x2 bending compliance per station, all in ONE common frame:

    mesh x = edge direction, mesh y = flap direction, NOT rotated by the section twist.
    Moment/curvature components are about mesh x ("flap": M_f = Mx, kappa_f = kappa_x = d theta_x/dz) and about
    mesh y ("edge": M_e = My, kappa_e = kappa_y); C = d(kappa_f, kappa_e)/d(M_f, M_e), pure couple (F = 0, Mz = 0).
    C_fe = C[0, 1] is the flap curvature per edge moment, C_ef = C[1, 0] the reverse.  Principal angle psi = angle of
    the softest-moment axis (the one carrying the largest compliance, the flap-like axis) from mesh +x toward +y,
    range (-90, 90]; a section whose chord lies along mesh x has psi = 0, a twisted one psi ~ its chord angle.

(a) model level, linear parked solves under the SAME loads: a pure tip couple (Mx, then My) on the shell (moment-equivalent
    z-force couple over the tip-ring nodes, zero net force) and on BeamDyn (Mx/My point load at eta = 1, mapped to the BD
    frame).  A tip couple gives a UNIFORM internal moment, so the curvature is the compliance directly (no moment balance
    of distributed rows).  Shell: u_z = c + e dz + (theta_x + kappa_x dz) y' - (theta_y + kappa_y dz) x' fitted by least
    squares to the nodes of a slab |z - z0| <= W (rotation about x and y of the section and its z-rate in one fit).
    BeamDyn: Wiener-Milenkovic RD{x,y}r -> rotation vector -> mesh frame (v_mesh = v_bd @ MESH_TO_BD, twist_p1.py:37),
    slope of theta over the nodes within W of z0.
(b) deck: the BD blade 6x6 compliance S = K^-1 (rotational block = curvature per pure couple, all couplings condensed),
    rotated from the section frame to the blade frame by R = tilt * Rz(-kp_twist) (see DECK FRAME below) and to the mesh
    frame; per arc length -> per z by 1/t_z.
(c) independent reference: thin-wall Steiner integral of the windIO layup (tools/yaml_ei_independent.strips, mesh-free)
    with its cross term kept, expressed in the mesh frame by an ICP-like alignment of the yaml contour to the mesh ring,
    plus the mesh-cut thin-wall tensor (tools/l3_w3_ei_edge.py, mesh frame natively) as a second reading.

DECK FRAME (BeamDyn source, OpenFAST main; fetched 2026-11-15):
  * https://openfast.readthedocs.io/en/main/source/user/beamdyn/input_files.html section 4.7.3.4 (Blade Input File): "these
    matrices are defined in a local coordinate system along the blade axis with Z_l directing toward the unit tangent
    vector of the blade reference axis"; section 4.7.3.3.2: the reference axis "locates the origin and orientation of each
    a local coordinate system where the cross-sectional 6x6 stiffness and mass matrices are defined".
  * theory.html 4.7.5.1.4: the local blade system has Z_l tangent to the (deflected) blade axis.
  * modules/beamdyn/src/BeamDyn_IO.f90:898  ``kp_coordinate(i,4) = -TmpReAry(4)``  (the file's initial_twist is NEGATED on
    reading: IEC twist is about -Z_l).
  * BeamDyn_Subs.f90:1030-1096 (BD_ComputeIniNodalCrv): Rr = Rtwist(+phi about the tangent e3) * [e1 e2 e3], i.e. a
    right-hand rotation of the minimal-tilt frame by the (already negated) angle.
  * BeamDyn.f90:2692-2710 (Calc_Stif_betaC, called from BD_StifAtDeformedQP) Stif = RR0 Stif0 RR0^T, RR0 = R(deformed) R0:
    the 6x6 in the deck is the SECTION-frame matrix and BeamDyn rotates it at run time by the twisted frame.  The deck K is
    therefore NOT pre-rotated by the structural twist.
  ``tools/twist_p1.py:section_frames`` takes the twist angle as given; passing +kp_twist (as ``Deck`` in twist_p1 is
  built in ``run_case``) is the opposite sense of BeamDyn.  View (b) is evaluated with both senses and the BD run
  (a) arbitrates (the cross-check below).

Usage (each view bounded, results under --out):
    scripts/aeroenv.sh python tools/l5_section_coupling.py --view bd deck ref shell --out logs/check/l5_oneway/t9
    scripts/aeroenv.sh python tools/l5_section_coupling.py --report --out logs/check/l5_oneway/t9
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "tools"))

import twist_p1 as tp  # noqa: E402

R_M = 117.0  # blade length (arc z), L7
FRACTIONS = (0.5, 0.6, 0.7, 0.8, 0.9, 0.95)
HALF_WINDOWS_SHELL = (2.5, 3.5, 5.0)  # m, slab half widths; the table uses WINDOW
HALF_WINDOWS_BD = (3.5, 5.0, 7.0)  # BeamDyn outputs ~51 nodes (2.34 m apart): >= 3 of them inside the window
WINDOW = 3.5
FIT_ORDER = 2  # polynomial degree in dz of T, theta_x, theta_y in the shell slab fit
STABLE_TOL = 0.05  # a shell number is reported as resolved when C_ff and C_ee agree to this over the windows
M0 = 1.0e4  # N.m, the unit tip couple
HUB_RADIUS = 3.97
TIP_TOL = 0.02  # m, tip-ring node selection


# --- pure helpers (unit-tested) ----------------------------------------------------------------------------------


def wm_to_rotvec(c):
    """Rotation vector (rad) of Wiener-Milenkovic parameters ``c = 4 tan(phi/4) n``; rows or a single vector."""
    c = np.asarray(c, float)
    norm = np.linalg.norm(c, axis=-1, keepdims=True)
    safe = np.where(norm == 0.0, 1.0, norm)
    return np.where(norm == 0.0, 0.0, 4.0 * np.arctan(norm / 4.0) * c / safe)


def bd_to_mesh(v):
    """BeamDyn blade-frame vector(s) -> mesh frame: v_mesh = v_bd @ MESH_TO_BD (twist_p1.py:37)."""
    return np.asarray(v, float) @ tp.MESH_TO_BD


def mesh_to_bd(v):
    return np.asarray(v, float) @ tp.MESH_TO_BD.T


def tip_couple_forces(xy, Mx, My):
    """z-forces on tip-ring nodes with zero net force and moments (about x, about y) = (Mx, My).

    Moment of a z-force: Mx = sum y f, My = -sum x f.  f = alpha (y - ybar) + beta (x - xbar).
    """
    xy = np.asarray(xy, float)
    x, y = xy[:, 0] - xy[:, 0].mean(), xy[:, 1] - xy[:, 1].mean()
    sxx, syy, sxy = (x * x).sum(), (y * y).sum(), (x * y).sum()
    alpha, beta = np.linalg.solve([[syy, sxy], [-sxy, -sxx]], [Mx, My])
    return alpha * y + beta * x


def fit_rotation_curvature(X, uz, z0, half, order=2):
    """Section rotation (theta_x, theta_y) at z0 and its z-rate (kappa_x, kappa_y) from the u_z of the slab |z - z0| <= half.

    u_z = T(dz) + theta_x(dz) y' - theta_y(dz) x',  T, theta_x, theta_y polynomials of degree ``order`` in dz = z - z0
    (a rigid rotation about x gives u_z = theta_x y, about y gives u_z = -theta_y x), (x', y') relative to the slab centroid.
    The rate is the dz coefficient of theta; with order >= 2 it is unbiased by the quadratic variation of theta and by the
    spanwise variation of the translation (which the tilt of the reference axis puts into u_z).
    """
    X, uz = np.asarray(X, float), np.asarray(uz, float)
    m = np.abs(X[:, 2] - z0) <= half
    if m.sum() < 6 * (order + 1):
        raise ValueError(f"slab at z={z0} +- {half} has only {int(m.sum())} nodes")
    xp, yp = X[m, 0] - X[m, 0].mean(), X[m, 1] - X[m, 1].mean()
    dz = X[m, 2] - z0
    pw = [dz**k for k in range(order + 1)]
    A = np.column_stack([*pw, *[p * yp for p in pw], *[-p * xp for p in pw]])
    coef, *_ = np.linalg.lstsq(A, uz[m], rcond=None)
    n = order + 1
    rms = float(np.sqrt(np.mean((A @ coef - uz[m]) ** 2)))
    rot = A[:, n:] @ coef[n:]  # the rotation part of the fit (what a rigid-section model explains)
    return {"theta": np.array([coef[n], coef[2 * n]]), "kappa": np.array([coef[n + 1], coef[2 * n + 1]]),
            "n": int(m.sum()), "rms": rms, "distortion": float(rms / max(np.std(rot), 1e-300))}


def slope_fit(z, y, z0, half):
    """Least-squares slope d y/d z of the samples with |z - z0| <= half (>= 3 of them); y may be (n, k)."""
    z, y = np.asarray(z, float), np.asarray(y, float)
    m = np.abs(z - z0) <= half
    if m.sum() < 3:
        raise ValueError(f"only {int(m.sum())} samples within {half} m of z={z0}")
    A = np.column_stack([np.ones(m.sum()), z[m] - z0])
    return np.linalg.lstsq(A, y[m], rcond=None)[0][1]


def principal_angle_deg(C):
    """Angle (deg, (-90, 90]) from mesh +x toward +y of the principal axis of the LARGEST compliance of sym(C)."""
    C = np.asarray(C, float)
    w, v = np.linalg.eigh(0.5 * (C + C.T))
    a = np.degrees(np.arctan2(v[1, np.argmax(w)], v[0, np.argmax(w)]))
    a = (a + 90.0) % 180.0 - 90.0
    return float(90.0 if a == -90.0 else a)


def compliance_from_responses(kappa_Mx, kappa_My, moment):
    """C (2x2) from the curvature pairs (kappa_x, kappa_y) under unit-direction couples Mx and My of size ``moment``."""
    return np.column_stack([np.asarray(kappa_Mx, float), np.asarray(kappa_My, float)]) / moment


def deck_bending_compliance(K, R, tz=1.0):
    """2x2 bending compliance in the mesh frame from a section-frame 6x6 ``K`` and the frame ``R`` (columns = section axes).

    S = K^-1; the rotational block S[3:, 3:] is the curvature per pure couple (force and torque free), R S R^T puts it in the
    blade frame, MESH_TO_BD^T into the mesh frame; /tz converts per arc length to per z.
    """
    S = np.linalg.inv(np.asarray(K, float))[3:, 3:]
    Cb = R @ S @ R.T
    P = tp.MESH_TO_BD
    return (P.T @ Cb @ P)[:2, :2] / tz


def stiffness_block_to_compliance(block):
    return np.linalg.inv(np.asarray(block, float))


def thin_wall_block(strips, A=None):
    """[[Ixx, -Ixy], [-Ixy, Iyy]] (curvature -> moment) about the EA-weighted centroid of (EA, x, y) point masses.

    ``A`` (2x2 orthogonal) maps the strip coordinates into the target frame first.
    """
    s = np.asarray(strips, float)
    xy = s[:, 1:3] if A is None else s[:, 1:3] @ np.asarray(A, float).T
    w = s[:, 0]
    d = xy - (w[:, None] * xy).sum(0) / w.sum()
    ixx, iyy, ixy = (w * d[:, 1] ** 2).sum(), (w * d[:, 0] ** 2).sum(), (w * d[:, 0] * d[:, 1]).sum()
    return np.array([[ixx, -ixy], [-ixy, iyy]])


def farthest_pair(xy):
    xy = np.asarray(xy, float)
    d = np.linalg.norm(xy[:, None] - xy[None], axis=2)
    i, j = np.unravel_index(np.argmax(d), d.shape)
    return xy[i], xy[j]


def align_candidates(contour, ring):
    """Orthogonal maps (A, t) taking yaml-contour points p to A p + t onto the mesh ring, with their mean-NN distance.

    Four candidates from the chord endpoints (farthest pair): two endpoint assignments x reflection.  The ring chord
    sets rotation and translation; the reflection and the LE/TE sense are decided by the geometry residual alone.
    """
    from scipy.spatial import cKDTree

    a, b = farthest_pair(contour[:: max(1, len(contour) // 500)])
    c, d = farthest_pair(ring)
    e1 = (b - a) / np.linalg.norm(b - a)
    e = np.array([e1, [-e1[1], e1[0]]])  # rows
    tree_pts = contour
    out = []
    for swap in (False, True):
        p, q = (c, d) if not swap else (d, c)
        f1 = (q - p) / np.linalg.norm(q - p)
        for refl in (1.0, -1.0):
            f = np.array([f1, refl * np.array([-f1[1], f1[0]])])
            A = f.T @ e
            t = p - A @ a
            moved = tree_pts @ A.T + t
            dist = cKDTree(moved).query(ring)[0].mean()
            out.append({"A": A, "t": t, "residual": float(dist), "reflect": bool(refl < 0), "swap": swap})
    return sorted(out, key=lambda r: r["residual"])


def chord_angle_deg(ring):
    """Line angle (deg, (-90, 90]) of the chord (farthest pair of the ring) from mesh +x toward +y."""
    a, b = farthest_pair(ring)
    ang = np.degrees(np.arctan2(*(b - a)[::-1]))
    ang = (ang + 90.0) % 180.0 - 90.0
    return float(90.0 if ang == -90.0 else ang)


def table_row(C):
    C = np.asarray(C, float)
    return {"C": C.tolist(), "C_ff": float(C[0, 0]), "C_ee": float(C[1, 1]), "C_fe": float(C[0, 1]), "C_ef": float(C[1, 0]),
            "ratio_fe_ff": float(C[0, 1] / C[0, 0]), "ratio_ef_ee": float(C[1, 0] / C[1, 1]),
            "asymmetry": float((C[0, 1] - C[1, 0]) / (0.5 * abs(C[0, 1] + C[1, 0])) if (C[0, 1] + C[1, 0]) else 0.0),
            "psi_deg": principal_angle_deg(C)}


def stations(z_max=R_M):
    return [f * z_max for f in FRACTIONS]


# --- view (a) BeamDyn ---------------------------------------------------------------------------------------------


def run_bd_couples(work, moment=M0):
    """BeamDyn static pure tip couples Mx, My (mesh frame) at eta = 1; returns z and the mesh-frame rotation vectors."""
    import diagnose_radii_beamdyn as dmr
    import l5_oneway_bd_parity as P

    tp.HUB_RADIUS = HUB_RADIUS
    ctx = {"work": Path(work), "deck_dir": P.DECK_DIR, "beamdyn": P.BEAMDYN, "dmr": dmr, "hashes": {}}
    out = {}
    for name, M in (("Mx", [1.0, 0.0, 0.0]), ("My", [0.0, 1.0, 0.0])):
        rows = np.array([[1.0, 0.0, 0.0, 0.0, *mesh_to_bd(np.array(M) * moment)]])
        deck, names, data = tp.run_case(ctx, f"couple_{name}", dynamic=False, rows=rows)
        _, rd = tp.node_fields(names, tp.pick_row(data, False), len(deck.eta))
        out[name] = bd_to_mesh(wm_to_rotvec(rd))
        out["z"] = deck.qp[:, 2]
    return out, ctx["hashes"]


def view_bd(work, out_dir):
    res, hashes = run_bd_couples(Path(work) / "bd")
    res_lin, _ = run_bd_couples(Path(work) / "bd_lin", moment=0.1 * M0)  # linearity: same compliance at a tenth of the load
    rows = {}
    for half in HALF_WINDOWS_BD:
        per = []
        for z0 in stations():
            kx = [slope_fit(res["z"], res[n][:, :2], z0, half) for n in ("Mx", "My")]
            C = compliance_from_responses(kx[0], kx[1], M0)
            kl = [slope_fit(res_lin["z"], res_lin[n][:, :2], z0, half) for n in ("Mx", "My")]
            Cl = compliance_from_responses(kl[0], kl[1], 0.1 * M0)
            per.append({"z": z0, **table_row(C), "linearity_max_rel": float(np.max(np.abs(Cl / C - 1.0)))})
        rows[f"{half:g}"] = per
    doc = {"view": "beamdyn_couples", "moment_Nm": M0, "window_table": WINDOW, "windows": rows, "generated_sha256": hashes,
           "revision": tp.revision(REPO)}
    _write(out_dir / "bd.json", doc)
    return doc


# --- view (b) deck -----------------------------------------------------------------------------------------------


def deck_at(deck, z0):
    """Linear interpolation of K and the (re-orthonormalised) frame R and the tangent z-component at span z0."""
    z = deck.qp[:, 2]
    K = np.array([[np.interp(z0, z, deck.K[:, i, j]) for j in range(6)] for i in range(6)])
    R = np.array([[np.interp(z0, z, deck.R[:, i, j]) for j in range(3)] for i in range(3)])
    u, _, vt = np.linalg.svd(R)
    return K, u @ vt, float(np.interp(z0, z, deck.tangent[:, 2]))


def view_deck(out_dir):
    import l5_oneway_bd_parity as P

    primary = (P.DECK_DIR / "IEA-15-240-RWT_BeamDyn.dat").read_text()
    kpt = tp.read_keypoints(primary)
    blade = P.DECK_DIR / "IEA-15-240-RWT_BeamDyn_blade.dat"
    res = {}
    for label, sign in (("bd_sense_minus", -1.0), ("harness_sense_plus", 1.0)):
        deck = tp.Deck(kpt[:, :3], blade, sign * kpt[:, 3])
        per = []
        for z0 in stations():
            K, R, tz = deck_at(deck, z0)
            C = deck_bending_compliance(K, R, tz)
            Cs = deck_bending_compliance(K, np.eye(3), tz)  # section frame (no rotation), for the record
            per.append({"z": z0, "kp_twist_deg": float(np.interp(z0, kpt[:, 2], kpt[:, 3])), "tz": tz,
                        **table_row(C), "section_frame_ratio_fe_ff": float(Cs[0, 1] / Cs[0, 0]),
                        "K44": float(K[3, 3]), "K55": float(K[4, 4]), "K45": float(K[3, 4]), "K34": float(K[2, 3]),
                        "K35": float(K[2, 4]), "K46": float(K[3, 5]), "K56": float(K[4, 5])})
        res[label] = per
    doc = {"view": "deck_6x6", "twist_sign_in_R": {"bd_sense_minus": "Rz(-kp_twist), BeamDyn_IO.f90:898",
                                                   "harness_sense_plus": "Rz(+kp_twist), tools/twist_p1.py run_case"},
           "results": res, "blade_file_sha256": tp.sha256(blade)}
    _write(out_dir / "deck.json", doc)
    return doc


# --- view (c) reference ---------------------------------------------------------------------------------------------


def yaml_compliance(data, eta, X, z0):
    """(c): 2x2 compliance of the windIO thin-wall section at ``eta`` in the mesh frame (contour aligned to the mesh ring at z0)."""
    from yaml_ei_independent import strips
    from yaml_layup_independent import contour

    ring = X[np.abs(X[:, 2] - z0) <= 0.25, :2]
    _, cxy, _ = contour(data, eta)
    cand = align_candidates(cxy, ring)
    sp = strips(data, eta)
    K2 = thin_wall_block(sp, cand[0]["A"])
    return stiffness_block_to_compliance(K2), K2, cand, sp, ring


def view_ref(out_dir):
    import l5_oneway_bd_parity as P
    from yaml_layup_independent import interp, load

    mesh, _, _ = P.build_shell(20)
    X = np.array([n.coords for n in mesh.nodes], float)
    data = load()
    twist_grid = data["components"]["blade"]["outer_shape_bem"]["twist"]
    rows = []
    for z0 in stations():
        eta = z0 / R_M
        C, K2, cand, sp, ring = yaml_compliance(data, eta, X, z0)
        best = cand[0]
        K2_yaml = thin_wall_block(sp)  # in the yaml frame, for the record
        rows.append({"z": z0, "eta": eta, **table_row(C), "EI_flap_Nm2": float(K2[0, 0]), "EI_edge_Nm2": float(K2[1, 1]),
                     "EI_cross_Nm2": float(-K2[0, 1]), "chord_angle_mesh_deg": chord_angle_deg(ring),
                     "windio_twist_deg": float(np.degrees(interp(twist_grid, eta))),
                     "alignment": {"best_residual_m": best["residual"], "second_residual_m": cand[1]["residual"],
                                   "reflect": best["reflect"], "swap": best["swap"],
                                   "rotation_deg": float(np.degrees(np.arctan2(best["A"][1, 0], best["A"][0, 0]))),
                                   "n_ring_nodes": int(len(ring))},
                     "principal_angle_yaml_frame_deg": principal_angle_deg(stiffness_block_to_compliance(K2_yaml)),
                     "n_strips": int(len(sp))})
    cut_path = out_dir / "l3_w3_mesh_cut.json"
    subprocess.run([sys.executable, "-B", str(REPO / "tools" / "l3_w3_ei_edge.py"), str(cut_path), "20"], check=True,
                   capture_output=True, text=True, cwd=REPO, timeout=500)
    cut = json.loads(cut_path.read_text())
    mesh_cut = []
    for z0 in stations():
        eta = round(z0 / R_M, 2)
        r = next(r for r in cut["rows"] if abs(r["eta"] - eta) < 1e-9)
        for ref in ("mid", "oml"):
            m = r[ref]
            ang = np.radians(m["ang"])  # direction of the max-EI eigenvector in the (kappa_x, kappa_y) plane
            Rm = np.array([[np.cos(ang), -np.sin(ang)], [np.sin(ang), np.cos(ang)]])
            K2 = Rm @ np.diag([m["EI_edge"], m["EI_flap"]]) @ Rm.T  # Rm columns: [max-axis, min-axis]
            mesh_cut.append({"z": z0, "surface": ref, **table_row(stiffness_block_to_compliance(K2)),
                             "EI_flap_Nm2": m["EI_flap"], "EI_edge_Nm2": m["EI_edge"], "ang_max_axis_deg": m["ang"]})
    doc = {"view": "reference", "yaml_thin_wall": rows, "mesh_cut_thin_wall": mesh_cut,
           "tools": ["tools/yaml_ei_independent.py:strips", "tools/yaml_layup_independent.py", "tools/l3_w3_ei_edge.py"]}
    _write(out_dir / "ref.json", doc)
    return doc


# --- view (a) shell -----------------------------------------------------------------------------------------------


def view_shell(out_dir):
    import l5_oneway_bd_parity as P
    from run_s5_oneway import solve_static

    t0 = time.time()
    mesh, _, asm = P.build_shell(20)
    dpn = asm.dofs_per_node
    X = np.array([n.coords for n in mesh.nodes], float)
    K, free = P.rotating_keff(asm, mesh, 0.0, HUB_RADIUS)  # omega = 0: K only
    (Path(out_dir) / "work").mkdir(parents=True, exist_ok=True)
    tip = np.flatnonzero(X[:, 2] >= X[:, 2].max() - TIP_TOL)
    print(f"shell ready ({mesh.node_count} nodes, {len(tip)} tip nodes, {time.time() - t0:.0f} s)", flush=True)
    uz = {}
    for name, (mx, my) in (("Mx", (M0, 0.0)), ("My", (0.0, M0))):
        f = np.zeros(asm.dofs_count)
        f[tip * dpn + 2] = tip_couple_forces(X[tip, :2], mx, my)
        assert abs(f[2::dpn].sum()) < 1e-6 * M0 and abs((f[2::dpn] * X[:, 1]).sum() - mx) < 1e-6 * M0
        u = solve_static(asm, K, free, f)
        uz[name] = u[2::dpn]
        np.save(Path(out_dir) / "work" / f"shell_uz_{name}.npy", uz[name])
        print(f"solved {name} ({time.time() - t0:.0f} s)", flush=True)
    rows = {}
    for half in HALF_WINDOWS_SHELL:
        per = []
        for z0 in stations():
            fits = {n: fit_rotation_curvature(X, uz[n], z0, half, FIT_ORDER) for n in uz}
            C = compliance_from_responses(fits["Mx"]["kappa"], fits["My"]["kappa"], M0)
            per.append({"z": z0, **table_row(C), "n_slab_nodes": fits["Mx"]["n"],
                        "fit_rms_m": {n: fits[n]["rms"] for n in fits},
                        "plane_fit_distortion": {n: fits[n]["distortion"] for n in fits},
                        "theta_Mx_rad": fits["Mx"]["theta"].tolist(), "theta_My_rad": fits["My"]["theta"].tolist()})
        rows[f"{half:g}"] = per
    doc = {"view": "shell_couples", "moment_Nm": M0, "n_nodes": mesh.node_count, "n_tip_nodes": int(len(tip)),
           "tip_tol_m": TIP_TOL, "window_table": WINDOW, "windows": rows, "revision": tp.revision(REPO),
           "elapsed_s": time.time() - t0}
    _write(out_dir / "shell.json", doc)
    return doc


# --- harness validation: uniform-section BeamDyn blade ---------------------------------------------------------------


class _UniformDeck:
    """Stand-in for ``diagnose_radii_beamdyn`` in ``twist_p1.run_case``: the real deck written, then made UNIFORM.

    Every station carries the 6x6 ``K``; the key points are a straight blade (x = y = 0) with the constant initial twist
    ``tw0`` (deg, file sense).  One order-10 element then represents the solution exactly (constant curvature), so the BD
    rotation slope is the discretisation-free compliance of K in the frame BeamDyn builds from ``tw0``.
    """

    def __init__(self, real, K, tw0):
        self.real, self.K, self.tw0 = real, np.asarray(K, float), tw0

    def write_deck(self, deck, work):
        self.real.write_deck(deck, work)
        prim = (Path(work) / "bd_primary.dat").read_text().splitlines()
        i0 = next(i for i, ln in enumerate(prim) if "kp_xr" in ln) + 2
        kp = []
        for i in range(i0, i0 + 50):
            z = float(prim[i].split()[2])
            prim[i] = f"0.0 0.0 {z:.8e} {self.tw0:.10f}"
            kp.append([0.0, 0.0, z])
        (Path(work) / "bd_primary.dat").write_text("\n".join(prim) + "\n")
        blade = Path(work) / "IEA-15-240-RWT_BeamDyn_blade.dat"
        st, _, M = tp.parse_blade(blade)
        lines = blade.read_text().splitlines()
        k0 = next(i for i, ln in enumerate(lines) if "DISTRIBUTED PROPERTIES" in ln) + 1
        fmt = lambda A: ["\t " + " \t ".join(f"{v:.16e}" for v in r) for r in A]  # noqa: E731
        out = lines[:k0]
        for e, m in zip(st, M, strict=True):
            out += [f"\t {e:.6f}", *fmt(self.K), "", *fmt(m), ""]
        blade.write_text("\n".join(out) + "\n")
        return np.array(kp)


def view_bdcheck(work, out_dir):
    """BD on a uniform blade built from the deck's K at each station vs view (b) with both twist senses."""
    import diagnose_radii_beamdyn as dmr
    import l5_oneway_bd_parity as P

    tp.HUB_RADIUS = HUB_RADIUS
    kpt = tp.read_keypoints((P.DECK_DIR / "IEA-15-240-RWT_BeamDyn.dat").read_text())
    deck = tp.Deck(kpt[:, :3], P.DECK_DIR / "IEA-15-240-RWT_BeamDyn_blade.dat", -kpt[:, 3])
    z_axis = np.array([0.0, 0.0, 1.0])
    rows = []
    for z0 in stations():
        K, _, _ = deck_at(deck, z0)
        tw0 = float(np.interp(z0, kpt[:, 2], kpt[:, 3]))
        ctx = {"work": Path(work) / "uniform", "deck_dir": P.DECK_DIR, "beamdyn": P.BEAMDYN,
               "dmr": _UniformDeck(dmr, K, tw0), "hashes": {}}
        resp = {}
        for name, M in (("Mx", [1.0, 0.0, 0.0]), ("My", [0.0, 1.0, 0.0])):
            r = np.array([[1.0, 0.0, 0.0, 0.0, *mesh_to_bd(np.array(M) * M0)]])
            d, names, data = tp.run_case(ctx, f"u{z0:.0f}_{name}", dynamic=False, rows=r)
            _, rd = tp.node_fields(names, tp.pick_row(data, False), len(d.eta))
            resp[name] = slope_fit(d.qp[:, 2], bd_to_mesh(wm_to_rotvec(rd))[:, :2], 0.5 * d.qp[-1, 2], 1e9)
        C = compliance_from_responses(resp["Mx"], resp["My"], M0)
        entry = {"z": z0, "kp_twist_deg": tw0, "bd_uniform": table_row(C)}
        for label, sign in (("bd_sense_minus", -1.0), ("harness_sense_plus", 1.0)):
            Cb = deck_bending_compliance(K, tp.rotation_matrix(z_axis, np.deg2rad(sign * tw0)))
            entry[label] = {**table_row(Cb), "max_rel_diff_vs_bd": float(np.max(np.abs(Cb / C - 1.0)))}
        rows.append(entry)
    doc = {"view": "bd_uniform_blade_check", "rows": rows,
           "note": "BD on a straight uniform blade (K of the station, constant twist): exact in one element"}
    _write(out_dir / "bdcheck.json", doc)
    return doc


# --- 1D prediction of the L8 deflections from the deck compliance ------------------------------------------------------


def tip_deflections(z, C, fx, fy):
    """Tip (u_x, u_y) of a cantilever with per-z compliance ``C`` (n, 2, 2) in the mesh frame under distributed forces.

    f_x, f_y are N/m along mesh x, y at the span positions ``z``; Mx = -int_z^L (s - z) f_y ds, My = +int_z^L (s - z) f_x ds
    (moment of a transverse force about the axes through the section, F_z = 0); kappa = C (Mx, My); small rotations:
    u_x'' = kappa_y, u_y'' = -kappa_x, clamped at the root.
    """
    z = np.asarray(z, float)

    def tail(g):  # int_z^L g ds on the grid
        seg = 0.5 * (g[1:] + g[:-1]) * np.diff(z)
        return np.concatenate([np.cumsum(seg[::-1])[::-1], [0.0]])

    s0, s1 = tail(np.asarray(fy)), tail(z * np.asarray(fy))
    t0, t1 = tail(np.asarray(fx)), tail(z * np.asarray(fx))
    Mx, My = -(s1 - z * s0), (t1 - z * t0)
    kap = np.einsum("nij,nj->ni", C, np.column_stack([Mx, My]))
    w = z[-1] - z
    u = np.array([np.trapezoid(w * kap[:, 1], z), -np.trapezoid(w * kap[:, 0], z)])
    return u


def view_predict(out_dir):
    """Edge/flap tip deflection under the L8 flap-only and edge-only AeroDyn loads from the deck compliance, both senses."""
    import l5_oneway_bd_parity as P

    kpt = tp.read_keypoints((P.DECK_DIR / "IEA-15-240-RWT_BeamDyn.dat").read_text())
    blade = P.DECK_DIR / "IEA-15-240-RWT_BeamDyn_blade.dat"
    r_m, fn, ft, twist = np.loadtxt(P.AD_LOADS, delimiter=",", skiprows=1).T
    th = np.deg2rad(twist)
    # force per length in the mesh frame: normal n = (-sin th, cos th), tangential t = (cos th, sin th) (run_s5_oneway.apply_ad_loads)
    fx_all, fy_all = -fn * np.sin(th) + ft * np.cos(th), fn * np.cos(th) + ft * np.sin(th)
    res = {}
    for label, sign in (("bd_sense_minus", -1.0), ("harness_sense_plus", 1.0)):
        deck = tp.Deck(kpt[:, :3], blade, sign * kpt[:, 3])
        eta, arc, K, R = deck.dense(16)
        z = np.interp(eta, deck.arc_kp / deck.length, deck.kp[:, 2])
        tz = np.interp(eta, deck.arc_kp / deck.length, np.gradient(deck.kp[:, 2], deck.arc_kp))
        C = np.array([deck_bending_compliance(K[i], R[i], tz[i]) for i in range(len(eta))])
        keep = np.concatenate([[True], np.diff(z) > 0])
        z, C = z[keep], C[keep]
        fxi, fyi = np.interp(z, r_m, fx_all), np.interp(z, r_m, fy_all)
        zero = np.zeros_like(z)
        res[label] = {"flap_only_load(u_x edge, u_y flap)": tip_deflections(z, C, zero, fyi).tolist(),
                      "edge_only_load(u_x edge, u_y flap)": tip_deflections(z, C, fxi, zero).tolist()}
    # (c): thin-wall yaml compliance on eta 0.25..0.97 (mesh frame, aligned to the mesh ring at each station); the deck with the
    # MESH twist sense (Rz(+tw)) outside that range, where the load arm is largest but the twist cross terms are weakest
    import l5_oneway_bd_parity as P2
    from yaml_layup_independent import load

    mesh, _, _ = P2.build_shell(20)
    X = np.array([n.coords for n in mesh.nodes], float)
    data = load()
    deck = tp.Deck(kpt[:, :3], blade, kpt[:, 3])
    eta, arc, K, R = deck.dense(16)
    z = np.interp(eta, deck.arc_kp / deck.length, deck.kp[:, 2])
    tz = np.interp(eta, deck.arc_kp / deck.length, np.gradient(deck.kp[:, 2], deck.arc_kp))
    C = np.array([deck_bending_compliance(K[i], R[i], tz[i]) for i in range(len(eta))])
    ge = np.linspace(0.25, 0.97, 37)
    Cy = np.array([yaml_compliance(data, e, X, e * R_M)[0] for e in ge])
    inside = (z >= 0.25 * R_M) & (z <= 0.97 * R_M)
    for i in np.flatnonzero(inside):
        C[i] = [[np.interp(z[i], ge * R_M, Cy[:, a, b]) for b in range(2)] for a in range(2)]
    keep = np.concatenate([[True], np.diff(z) > 0])
    zz, Cc = z[keep], C[keep]
    fxi, fyi, zero = np.interp(zz, r_m, fx_all), np.interp(zz, r_m, fy_all), np.zeros(keep.sum())
    res["reference_c_yaml_0.25-0.97R_deck_plus_outside"] = {
        "flap_only_load(u_x edge, u_y flap)": tip_deflections(zz, Cc, zero, fyi).tolist(),
        "edge_only_load(u_x edge, u_y flap)": tip_deflections(zz, Cc, fxi, zero).tolist()}
    doc = {"view": "deck_1d_prediction", "results": res,
           "l8_reference_tip_m": {"flap_only_load_edge_cross": {"shell": 0.048, "beamdyn": 0.591},
                                  "edge_only_load_flap_cross": {"shell": 0.026, "beamdyn": 0.081},
                                  "source": "odd/tasks/l5-one-way.md L8"}}
    _write(out_dir / "predict.json", doc)
    return doc


# --- report ------------------------------------------------------------------------------------------------------


def _write(path, doc):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(doc, indent=1) + "\n")


def _load(out_dir, name):
    p = Path(out_dir) / f"{name}.json"
    return json.loads(p.read_text()) if p.exists() else None


def rel(a, b):
    return float(a / b - 1.0)


def verdict(model, ref, tol_c=0.05, tol_ang=1.0):
    """Per-station and overall verdict of ``model`` against ``ref`` (rows with C_ff, C_ee, psi_deg)."""
    rows = []
    for m, r in zip(model, ref, strict=True):
        dff, dee, dpsi = rel(m["C_ff"], r["C_ff"]), rel(m["C_ee"], r["C_ee"]), m["psi_deg"] - r["psi_deg"]
        rows.append({"z": m["z"], "dC_ff": dff, "dC_ee": dee, "dpsi_deg": dpsi,
                     "ok": abs(dff) <= tol_c and abs(dee) <= tol_c and abs(dpsi) <= tol_ang})
    return {"rows": rows, "all_ok": all(r["ok"] for r in rows), "tol_C": tol_c, "tol_psi_deg": tol_ang}


def shell_stability(shell):
    """Per station: spread of C_ff and C_ee over the shell windows (relative to the table window) and the plane-fit distortion."""
    out = []
    key = f"{WINDOW:g}"
    for i, base in enumerate(shell["windows"][key]):
        spread = {}
        for name in ("C_ff", "C_ee"):
            vals = [w[i][name] for w in shell["windows"].values()]
            spread[name] = float((max(vals) - min(vals)) / abs(base[name]))
        dist = max(base["plane_fit_distortion"].values())
        out.append({"z": base["z"], "spread": spread, "distortion": dist,
                    "resolved": bool(max(spread.values()) <= STABLE_TOL)})
    return out


def report(out_dir):
    bd, deck, ref, shell = (_load(out_dir, n) for n in ("bd", "deck", "ref", "shell"))
    chk, pred = _load(out_dir, "bdcheck"), _load(out_dir, "predict")
    key = f"{WINDOW:g}"
    lines = ["# L5 T9 flap-edge cross-coupling arbiter", "",
             "Frame: mesh x = edge, mesh y = flap, not rotated by twist; C = d(kappa_x, kappa_y)/d(Mx, My) [1/(N m^2)];",
             "psi = angle of the largest-compliance (flap-like) principal axis from mesh +x toward +y.",
             f"Window (a): +-{WINDOW:g} m (shell slab fit, order {FIT_ORDER}; BeamDyn linear slope). Tip couple {M0:g} N m.", ""]
    doc = {"window_m": WINDOW}
    views = {}
    stab = shell_stability(shell) if shell else None
    if shell:
        views["shell (a)"] = shell["windows"][key]
    if bd:
        views["BeamDyn (a)"] = bd["windows"][key]
    if deck:
        views["deck (b) R=Rz(-tw) [BeamDyn sense]"] = deck["results"]["bd_sense_minus"]
        views["deck (b) R=Rz(+tw) [mesh sense]"] = deck["results"]["harness_sense_plus"]
    if ref:
        views["reference (c) yaml thin-wall"] = ref["yaml_thin_wall"]
        views["reference (c2) mesh-cut mid"] = [r for r in ref["mesh_cut_thin_wall"] if r["surface"] == "mid"]
        views["reference (c2) mesh-cut OML"] = [r for r in ref["mesh_cut_thin_wall"] if r["surface"] == "oml"]
    lines += ["| z [m] | view | C_ff | C_ee | C_fe | C_ef | C_fe/C_ff | psi [deg] |", "|---|---|---|---|---|---|---|---|"]
    for i, z0 in enumerate(stations()):
        for name, rows in views.items():
            r = rows[i]
            flag = " (UNRESOLVED: window spread / plane-fit distortion)" if name == "shell (a)" and not stab[i]["resolved"] else ""
            lines.append(f"| {z0:.2f} | {name}{flag} | {r['C_ff']:.4e} | {r['C_ee']:.4e} | {r['C_fe']:.3e} | {r['C_ef']:.3e} |"
                         f" {r['ratio_fe_ff']:+.4f} | {r['psi_deg']:+.2f} |")
    doc["views"] = views
    if stab:
        lines += ["", "## Shell (a) fit quality", "", "| z [m] | spread C_ff over windows | spread C_ee | plane-fit distortion (rms/rotation signal) | resolved |",
                  "|---|---|---|---|---|"]
        for s in stab:
            lines.append(f"| {s['z']:.2f} | {s['spread']['C_ff']:.3f} | {s['spread']['C_ee']:.3f} | {s['distortion']:.3f} | {s['resolved']} |")
        doc["shell_stability"] = stab
    if ref:
        lines += ["", "## Chord/twist angle per station", "",
                  ("| z [m] | mesh chord angle [deg] | windIO twist [deg] | deck kp twist [deg] | yaml-frame principal angle [deg] | "
                   "alignment: reflect / 2nd-best residual ratio |"), "|---|---|---|---|---|---|"]
        for i, r in enumerate(ref["yaml_thin_wall"]):
            kt = deck["results"]["bd_sense_minus"][i]["kp_twist_deg"] if deck else float("nan")
            lines.append(f"| {r['z']:.2f} | {r['chord_angle_mesh_deg']:+.2f} | {r['windio_twist_deg']:+.2f} | {kt:+.2f} |"
                         f" {r['principal_angle_yaml_frame_deg']:+.2f} | {r['alignment']['reflect']} / "
                         f"{r['alignment']['second_residual_m'] / max(r['alignment']['best_residual_m'], 1e-12):.1f}x |")
    verdicts = {}
    if bd and deck:
        for label in ("bd_sense_minus", "harness_sense_plus"):
            verdicts[f"BD(a) vs deck(b) {label}"] = verdict(bd["windows"][key], deck["results"][label])
    if ref:
        target = ref["yaml_thin_wall"]
        for name, rows in (("shell (a)", shell and shell["windows"][key]), ("BeamDyn (a)", bd and bd["windows"][key]),
                           ("deck (b) BeamDyn sense", deck and deck["results"]["bd_sense_minus"]),
                           ("deck (b) mesh sense", deck and deck["results"]["harness_sense_plus"])):
            if rows:
                v = verdict(rows, target)
                if name == "shell (a)":
                    v["rows"] = [r for r, s in zip(v["rows"], stab, strict=True) if s["resolved"]]
                    v["all_ok"] = all(r["ok"] for r in v["rows"])
                    v["note"] = "resolved stations only"
                verdicts[f"{name} vs (c)"] = v
    if shell and deck:
        for label, name in (("harness_sense_plus", "mesh sense"), ("bd_sense_minus", "BeamDyn sense")):
            v = verdict(shell["windows"][key], deck["results"][label])
            v["rows"] = [r for r, s in zip(v["rows"], stab, strict=True) if s["resolved"]]
            v["all_ok"], v["note"] = all(r["ok"] for r in v["rows"]), "resolved stations only"
            verdicts[f"shell (a) vs deck (b) {name}"] = v
    doc["verdicts"] = verdicts
    if ref:
        spread = {}
        c1, cm, co = (np.array([r["C_ff"] for r in x]) for x in (ref["yaml_thin_wall"], views["reference (c2) mesh-cut mid"],
                                                                  views["reference (c2) mesh-cut OML"]))
        e1, em, eo = (np.array([r["C_ee"] for r in x]) for x in (ref["yaml_thin_wall"], views["reference (c2) mesh-cut mid"],
                                                                  views["reference (c2) mesh-cut OML"]))
        spread = {"C_ff mid/yaml-1": (cm / c1 - 1).tolist(), "C_ff oml/yaml-1": (co / c1 - 1).tolist(),
                  "C_ee mid/yaml-1": (em / e1 - 1).tolist(), "C_ee oml/yaml-1": (eo / e1 - 1).tolist()}
        doc["reference_internal_spread"] = spread
    lines += ["", "## Verdicts (5 % on C_ff and C_ee, 1 deg on psi; per-station detail in summary.json)", "",
              "| comparison | stations | all ok | worst dC_ff | worst dC_ee | worst dpsi [deg] |", "|---|---|---|---|---|---|"]
    for name, v in verdicts.items():
        rr = v["rows"]
        lines.append(f"| {name} | {len(rr)} | {v['all_ok']} | {max(abs(r['dC_ff']) for r in rr):+.3f} | "
                     f"{max(abs(r['dC_ee']) for r in rr):+.3f} | {max(abs(r['dpsi_deg']) for r in rr):.2f} |")
    if ref:
        lines += ["", "(c) is a thin-wall estimate with its own spread: mesh-cut/yaml - 1 on C_ff (mid, OML) = "
                  + ", ".join(f"{a:+.2f}/{b:+.2f}" for a, b in zip(spread["C_ff mid/yaml-1"], spread["C_ff oml/yaml-1"], strict=True))
                  + "; on C_ee = " + ", ".join(f"{a:+.2f}/{b:+.2f}" for a, b in
                                              zip(spread["C_ee mid/yaml-1"], spread["C_ee oml/yaml-1"], strict=True))]
    if chk:
        lines += ["", "## Harness check: BD on a uniform blade (K of the station, constant twist) vs view (b)", "",
                  "| z [m] | kp twist [deg] | BD C_fe/C_ff | psi_BD | max rel diff (b) Rz(-tw) | max rel diff (b) Rz(+tw) |", "|---|---|---|---|---|---|"]
        for r in chk["rows"]:
            lines.append(f"| {r['z']:.2f} | {r['kp_twist_deg']:+.2f} | {r['bd_uniform']['ratio_fe_ff']:+.4f} | "
                         f"{r['bd_uniform']['psi_deg']:+.2f} | {r['bd_sense_minus']['max_rel_diff_vs_bd']:.1e} | "
                         f"{r['harness_sense_plus']['max_rel_diff_vs_bd']:.1e} |")
        doc["bd_uniform_check"] = chk["rows"]
    if pred:
        lines += ["", "## 1D prediction of the L8 tip deflections [m] (parked, AeroDyn rows, deck/yaml compliance)", "",
                  "| compliance | flap-only load: u_edge / u_flap | edge-only load: u_edge / u_flap |", "|---|---|---|"]
        for k, v in pred["results"].items():
            a, b = v["flap_only_load(u_x edge, u_y flap)"], v["edge_only_load(u_x edge, u_y flap)"]
            lines.append(f"| {k} | {a[0]:+.3f} / {a[1]:+.3f} | {b[0]:+.3f} / {b[1]:+.3f} |")
        lines.append("| L8 shell (3-D, magnitudes) | 0.048 / - | - / 0.026 |")
        lines.append("| L8 BeamDyn (magnitudes) | 0.591 / - | - / 0.081 |")
        doc["prediction"] = pred
    _write(Path(out_dir) / "summary.json", doc)
    (Path(out_dir) / "summary.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--view", nargs="*", default=[], choices=("bd", "deck", "ref", "shell", "bdcheck", "predict"))
    ap.add_argument("--report", action="store_true")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--work", type=Path, default=None, help="BeamDyn case directory (default: <out>/work)")
    args = ap.parse_args()
    work = args.work or args.out / "work"
    for v in args.view:
        t0 = time.time()
        {"bd": lambda: view_bd(work, args.out), "deck": lambda: view_deck(args.out), "ref": lambda: view_ref(args.out),
         "shell": lambda: view_shell(args.out), "bdcheck": lambda: view_bdcheck(work, args.out),
         "predict": lambda: view_predict(args.out)}[v]()
        print(f"view {v} done ({time.time() - t0:.0f} s)", flush=True)
    if args.report:
        report(args.out)


if __name__ == "__main__":
    main()
