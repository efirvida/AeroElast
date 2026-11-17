"""L5 T9 arbiter: who is right about the flap-edge coupling, principal axes and torsion, our shell or the deck 6x6?
(odd/tasks/l5-one-way.md T9, L9-L10, L28; solid segment tooling reused from tools/l3_solid_section_ccx.py, L3 T4)

Three models of the SAME cut, one post-processing, one frame, one 3x3 section compliance per station:

  deck   BeamDyn deck 6x6 at the station: C = (K^-1)[3:, 3:], rotated by Rz(-initial_twist) (BeamDyn_IO.f90:898, the
         section frame is rotated at run time, BeamDyn.f90:2692-2710) into the blade frame and by MESH_TO_BD
         (tools/twist_p1.py:37) into the mesh frame; per unit z by 1/t_z (tangent z component).  Force-free.
  shell  segment cut from the production gmsh chord-scaled N20 AR 1 mid-surface mesh (tests/support/v2_blade.py,
         production MITC4+ elements via MeshAssembler): the elements whose centroid lies within +-L/2 of the station
         (L = n_chord * local chord).  End A clamped (all 6 dofs of the nodes within 0.05 L of the end); end B a rigid
         body (all 6 dofs of the nodes within 0.05 L tied to a reference node at the end centroid), the unit couple
         applied at the reference node = zero net force, a uniform internal moment, warping restrained at both ends.
  solid  CalculiX C3D20R layered solid of the yaml section (tools/l3_solid_section_ccx.py: layers grown inward from
         the OML, webs, orthotropic local frames), prismatic, same clamped / rigid-body ends and reference-node couple.
         The solid section is the yaml section at the station; skin_points() are already in the mesh frame (LE +x,
         suction +y, chord at -twist, checked against the mesh ring at z 60/80/95: psi within 0.2 deg), so no rotation
         is applied (``frame_deg`` exists for a check).  No TE filler unless --te-filler (the shell mesh has none, the
         deck's SONATA section does: L23).

Frame, units, signs (all three): MESH frame, x = leading-edge direction, y = suction / downwind (flap), z = span (m).
  loads    M_flap = moment about +x, M_edge = moment about +y, T = torque about +z, right-handed, N m.
  strains  kappa_flap = d theta_x / dz, kappa_edge = d theta_y / dz, twist rate = d theta_z / dz, rad/m, rotations
           u = theta x r (u_z = theta_x y - theta_y x, u_x = -theta_z y, u_y = theta_z x).
  C        3x3 compliance, rad/m per N m: C[i, j] = strain_i per unit load_j, order (flap, edge, torque), force and the
           other two couples zero.  ratio_fe_ff = C[0,1]/C[0,0] (flap curvature per edge moment over the flap compliance,
           the table_row convention of tools/l5_section_coupling.py), ratio_ef_ee = C[1,0]/C[1,1].
  psi      principal bending angle, deg in (-90, 90]: direction (from mesh +x toward +y, i.e. toward the edge axis for
           +psi) of the largest-compliance (flap-like) axis of sym(C[:2,:2]); a section with the chord along x has 0.
  GJ       1/C[2,2] (pure torque).  EI_flap = 1/C[0,0], EI_edge = 1/C[1,1] (diagonal compliance inverses, the pure-couple
           response); EI_min / EI_max = eigenvalues of inv(sym C[:2,:2]) (frame invariant).
  Section kinematics (shell and solid, the same function): least-squares fit over the nodes of the window |z - z_c| <= half
  of u = rigid section motion with translation cubic and rotation quadratic polynomials in dz; the rates at dz = 0 are the
  strains.  Window ``third`` = L/6 (the middle third of the segment), ``half`` = L/4.

Not evidence by itself: the solid and shell segments are SHORT (a few chords), the shell one is tapered and twisted (the
cut keeps the production geometry) while the solid and the deck are per-station; read the three-way table, not one row.

  scripts/aeroenv.sh env PYTHONPATH=$PWD OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 python tools/l5_t9_arbiter.py deck  --out DIR
  ...                                                                                     shell --z 60 80 95 --out DIR
  ...                                                                                     solid --z 95 --level A --out DIR  [--cmd build|run|analyze]
  ...                                                                                     table --out DIR [--solid-dirs D ...]
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "tools"))

import twist_p1 as tp  # noqa: E402

YAML = REPO / "tests/IEA15MW/reference/upstream/IEA-15-240-RWT.yaml"
DECK_DIR = Path("/scratch/leahk/eduardo.donestevez/ofruns/OpenFAST/IEA-15-240-RWT")
CCX = os.environ.get("CCX_BIN", "/scratch/leahk/eduardo.donestevez/venv/bin/ccx")
STATIONS = (60.0, 80.0, 95.0)
M0 = 1.0e6  # N m, the unit couple (linear problems)
N_CHORD = 4.0  # solid segment length in chords (the validated W1 setting)
N_SHELL_CHORD = 8.0  # shell segment length in chords: the converged setting of the T9c study (converged_setting, odd/tasks/l5-one-way.md T9c)
SHELL_WINDOW = "half"  # shell fit window: half width L/4 (= 2 chords at 8 chords), the T9c converged setting
TOL_PSI = 0.5  # deg, convergence criterion of the shell study: a step to the next segment length moves psi by <= this
TOL_STIFF = 0.02  # relative, same for GJ, EI_flap, EI_edge (bias-corrected); the deficits to resolve are >= 6 %
END_ZONE = 0.05  # shell rigid / clamped zone, fraction of L
WINDOWS = {"third": 1.0 / 6.0, "half": 1.0 / 4.0}  # half widths, fraction of L
LEVELS = {  # solid mesh levels, relative to the local chord: A = the validated W1 mesh A (h 0.1 at chord 2.77), B = half h
    "A": {"h_rel": 0.036, "tmax": 0.03, "hz_rel": 0.108},
    "B": {"h_rel": 0.018, "tmax": 0.02, "hz_rel": 0.108},
}


# --- pure helpers (tests/software/test_l5_t9_arbiter.py) -------------------------------------------------------------


def rot_z3(deg):
    """Right-handed 3x3 rotation about +z."""
    c, s = np.cos(np.deg2rad(deg)), np.sin(np.deg2rad(deg))
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])


def compliance_from_responses(responses, loads):
    """C (3x3): column j = response vector (kappa_flap, kappa_edge, twist rate) under load j, divided by that load."""
    return np.column_stack([np.asarray(r, float) / float(m) for r, m in zip(responses, loads, strict=True)])


def rotate_compliance(C, R):
    """Compliance of the same section in a frame whose axes are the columns of ``R`` ... i.e. v' = R v: C' = R C R^T."""
    R = np.asarray(R, float)
    return R @ np.asarray(C, float) @ R.T


def principal_angle_deg(C):
    """Angle (deg, (-90, 90]) from +x toward +y of the largest-compliance axis of sym(C[:2, :2])."""
    C = np.asarray(C, float)[:2, :2]
    w, v = np.linalg.eigh(0.5 * (C + C.T))
    a = np.degrees(np.arctan2(v[1, np.argmax(w)], v[0, np.argmax(w)]))
    a = (a + 90.0) % 180.0 - 90.0
    return float(90.0 if a == -90.0 else a)


def section_metrics(C):
    """Everything the three-way table prints, from one 3x3 compliance (definitions in the module docstring)."""
    C = np.asarray(C, float)
    B = 0.5 * (C[:2, :2] + C[:2, :2].T)
    ev = np.sort(np.linalg.eigvalsh(np.linalg.inv(B)))
    d = np.sqrt(np.outer(np.abs(np.diag(C)), np.abs(np.diag(C))))
    return {"C": C.tolist(), "psi_deg": principal_angle_deg(C), "ratio_fe_ff": float(C[0, 1] / C[0, 0]),
            "ratio_ef_ee": float(C[1, 0] / C[1, 1]), "GJ": float(1.0 / C[2, 2]), "EI_flap": float(1.0 / C[0, 0]),
            "EI_edge": float(1.0 / C[1, 1]), "EI_min": float(ev[0]), "EI_max": float(ev[1]),
            "sym_err": float(np.max(np.abs(C - C.T) / d))}


def deck_compliance_mesh(K6, R_bd, tz=1.0):
    """3x3 (flap, edge, torque) compliance in the mesh frame from a section-frame BeamDyn 6x6 ``K6`` (order Fx Fy Fz Mx My Mz).

    S = (K^-1)[3:, 3:] is the rotation per unit couple with every force condensed out; R_bd (columns = section axes in the
    blade frame, e.g. Rz(-initial_twist) after the tilt) puts it in the blade frame, MESH_TO_BD (v_mesh = v_bd @ MESH_TO_BD)
    in the mesh frame; /tz converts per arc length to per unit z.
    """
    S = np.linalg.inv(np.asarray(K6, float))[3:, 3:]
    R = np.asarray(R_bd, float)
    P = tp.MESH_TO_BD
    return (P.T @ (R @ S @ R.T) @ P) / tz


def fit_section_kinematics(X, U, z0, half, weights=None, order_t=3, order_r=2):
    """Rotation (theta) and rotation rate (kappa) about mesh x, y, z at ``z0`` from the nodes with |z - z0| <= ``half``.

    X (n, 3) coordinates, U (n, 3) displacements.  Rigid-section model, polynomial in dz = z - z0 (translation degree
    ``order_t``, rotation degree ``order_r``), (x', y') relative to the window centroid:
      u_x = a(dz) - theta_z y' + e_xx x' + (g/2) y',  u_y = b(dz) + theta_z x' + e_yy y' + (g/2) x',  u_z = c(dz) + theta_x y' - theta_y x'
    (the in-plane affine strain e_xx, e_yy, g is free, so theta_z is the antisymmetric part of the in-plane gradient, the estimator
    of tools/l3_solid_section_ccx.py fit_plane; without it Poisson / ovalisation distortion biases the torsion rate).
    Returns kappa = d theta / d dz at dz = 0 (the strains), theta, rms residual, distortion = rms / std(rigid part).
    """
    X, U = np.asarray(X, float), np.asarray(U, float)
    m = np.abs(X[:, 2] - z0) <= half + 1e-9
    nt, nr = order_t + 1, order_r + 1
    nunk = 3 * nt + 6 * nr  # + in-plane affine strain (e_xx, e_yy, gamma_xy) polynomials: theta_z is the ANTISYMMETRIC part
    if m.sum() < 12:
        raise ValueError(f"slab at z={z0} +- {half} has only {int(m.sum())} nodes")
    xp, yp = X[m, 0] - X[m, 0].mean(), X[m, 1] - X[m, 1].mean()
    dz = (X[m, 2] - z0) / half  # scaled abscissa (conditioning); rates are divided back below
    pt = np.column_stack([dz**k for k in range(nt)])
    pr = np.column_stack([dz**k for k in range(nr)])
    n = int(m.sum())
    A = np.zeros((3 * n, nunk))
    A[0:n, 0:nt] = pt
    A[n:2 * n, nt:2 * nt] = pt
    A[2 * n:, 2 * nt:3 * nt] = pt
    o = 3 * nt  # rotations: theta_x | theta_y | theta_z
    A[2 * n:, o:o + nr] = pr * yp[:, None]
    A[2 * n:, o + nr:o + 2 * nr] = -pr * xp[:, None]
    A[0:n, o + 2 * nr:o + 3 * nr] = -pr * yp[:, None]
    A[n:2 * n, o + 2 * nr:o + 3 * nr] = pr * xp[:, None]
    q = o + 3 * nr  # u_x += e_xx x' + (gamma/2) y',  u_y += e_yy y' + (gamma/2) x'
    A[0:n, q:q + nr] = pr * xp[:, None]
    A[n:2 * n, q + nr:q + 2 * nr] = pr * yp[:, None]
    A[0:n, q + 2 * nr:q + 3 * nr] = 0.5 * pr * yp[:, None]
    A[n:2 * n, q + 2 * nr:q + 3 * nr] = 0.5 * pr * xp[:, None]
    y = np.concatenate([U[m, 0], U[m, 1], U[m, 2]])
    w = np.ones(n) if weights is None else np.asarray(weights, float)[m]
    sw = np.sqrt(np.concatenate([w, w, w]))
    if np.linalg.matrix_rank(A * sw[:, None]) < nunk:
        raise ValueError(f"slab at z={z0} +- {half}: rank-deficient section fit ({n} nodes)")
    coef, *_ = np.linalg.lstsq(A * sw[:, None], y * sw, rcond=None)
    res = A @ coef - y
    rigid = A[:, o:o + 3 * nr] @ coef[o:o + 3 * nr]
    rms = float(np.sqrt(np.mean(res**2)))
    return {"kappa": np.array([coef[o + 1], coef[o + nr + 1], coef[o + 2 * nr + 1]]) / half,
            "theta": np.array([coef[o], coef[o + nr], coef[o + 2 * nr]]), "n": n, "rms": rms,
            "distortion": float(rms / max(np.std(rigid), 1e-300))}


def rigid_end_map(points, ref):
    """Sparse (6 n, 6): the 6 dofs of each tied node from the reference node's [u_ref, theta_ref]: u = u_ref + theta x r, theta = theta_ref."""
    from scipy.sparse import coo_matrix

    pts = np.asarray(points, float)
    rows, cols, vals = [], [], []
    for i, p in enumerate(pts):
        r = p - np.asarray(ref, float)
        sk = np.array([[0.0, -r[2], r[1]], [r[2], 0.0, -r[0]], [-r[1], r[0], 0.0]])
        blk = np.block([[np.eye(3), -sk], [np.zeros((3, 3)), np.eye(3)]])
        for a in range(6):
            for b in range(6):
                if blk[a, b] != 0.0:
                    rows.append(6 * i + a)
                    cols.append(b)
                    vals.append(blk[a, b])
    return coo_matrix((vals, (rows, cols)), shape=(6 * len(pts), 6)).tocsr()


def eta_of_z(z, z_values, eta_grid):
    """windIO eta of the reference-axis span coordinate z (the mesh z)."""
    return float(np.interp(z, z_values, eta_grid))


# --- shared inputs -----------------------------------------------------------------------------------------------------


def yaml_axis():
    import yaml

    b = yaml.safe_load(YAML.read_text())["components"]["blade"]["reference_axis"]["z"]
    return np.array(b["values"], float), np.array(b["grid"], float)


def eta_for(z):
    zv, eg = yaml_axis()
    return eta_of_z(z, zv, eg)


def _write(path, doc):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(doc, indent=1, default=float))


def _entry(z, C, extra=None, windows=None):
    out = {"z": float(z), "eta": eta_for(z), **section_metrics(C)}
    if extra:
        out.update(extra)
    if windows:
        out["windows"] = windows
    return out


# --- (b) deck ---------------------------------------------------------------------------------------------------------


def load_deck():
    """BeamDyn deck (key points, 6x6 table) with the twist in BeamDyn's sense (-initial_twist about +z, BeamDyn_IO.f90:898)."""
    primary = (DECK_DIR / "IEA-15-240-RWT_BeamDyn.dat").read_text()
    kpt = tp.read_keypoints(primary)
    blade = DECK_DIR / "IEA-15-240-RWT_BeamDyn_blade.dat"
    return tp.Deck(kpt[:, :3], blade, -kpt[:, 3]), kpt, blade


def deck_C_of_z(deck):
    """z -> (3x3 mesh-frame compliance, K6, tz): the deck 6x6 (linear in eta, BeamDyn's interpolation) and its frame at span z."""
    import l5_section_coupling as sc

    def f(z):
        K, R, tz = sc.deck_at(deck, z)
        return deck_compliance_mesh(K, R, tz), K, tz

    return f


def run_deck(zs, out_dir):
    deck, kpt, blade = load_deck()
    f = deck_C_of_z(deck)
    rows = []
    for z in zs:
        C, K, tz = f(z)
        rows.append(_entry(z, C, {"kp_twist_deg": float(np.interp(z, kpt[:, 2], kpt[:, 3])), "tz": tz, "K66": float(K[5, 5])}))
    doc = {"model": "deck", "stations": rows, "blade_file_sha256": tp.sha256(blade)}
    _write(Path(out_dir) / "deck.json", doc)
    return doc


def beam_displacements(X, C_of_z, zlo, zhi, load, ref_xy, n=801):
    """Nodal displacements of a rigid-section BEAM with compliance C(z) under a uniform couple ``load`` (3-vector, N m), clamped at ``zlo``.

    kappa(z) = C(z) load; theta = int kappa dz; (u_x, u_y)' = (theta_y, -theta_x); a node at (x, y) moves by u + theta x r with
    r = (x, y) - ref_xy.  Used with ``fit_section_kinematics`` on the SAME nodes as the shell, it gives the value the shell fit
    would read from a section whose compliance is exactly the deck's: the taper / twist bias of the segment is then common.
    """
    X = np.asarray(X, float)
    zg = np.linspace(zlo, zhi, n)
    kap = np.array([C_of_z(z) @ np.asarray(load, float) for z in zg])

    def cum(f):
        return np.concatenate([[0.0], np.cumsum(0.5 * (f[1:] + f[:-1]) * np.diff(zg))])

    th = np.column_stack([cum(kap[:, i]) for i in range(3)])
    ux, uy = cum(th[:, 1]), cum(-th[:, 0])
    z = X[:, 2]
    t = [np.interp(z, zg, th[:, i]) for i in range(3)]
    rx, ry = X[:, 0] - ref_xy[0], X[:, 1] - ref_xy[1]
    return np.column_stack([np.interp(z, zg, ux) - t[2] * ry, np.interp(z, zg, uy) + t[2] * rx, t[0] * ry - t[1] * rx])


# --- (a) shell segment ------------------------------------------------------------------------------------------------


def shell_segment(mesh, props, z0, n_chord=N_CHORD):
    """Sub-mesh of the production mesh: elements whose centroid z lies within +-L/2 of z0, L = n_chord * local chord."""
    from aeroelast.core.mesh.entities import ElementSet
    from l5_twist_sense_chain import section_from_ring

    X = np.array([n.coords for n in mesh.nodes], float)
    sl = X[np.abs(X[:, 2] - z0) < 0.12]
    chord = float(section_from_ring(sl[:, :2])["chord"])
    L = n_chord * chord
    sel = {e for e in mesh.elements if abs(np.mean([n.coords[2] for n in e.nodes]) - z0) <= L / 2}
    name = f"T9SEG_{z0:g}_{n_chord:g}"
    if name not in mesh.element_sets:
        mesh.add_element_set(ElementSet(name=name, elements=sel))
    sub = mesh.extract_submesh(name)
    return sub, {k: v for k, v in props.items() if k in sub.element_sets}, chord, L


def solve_shell_couples(sub, props, M=M0):
    """Clamp end A, rigid end B with the reference-node couple (flap, edge, torque); returns (X, [U_flap, U_edge, U_torque], tributary
    mid-surface area per node, info)."""
    import scipy.sparse as sp
    from aeroelast.core.assembler import MeshAssembler
    from aeroelast.elements import ElementFamily
    from scipy.sparse.linalg import splu

    model = {"elements": {"element_family": ElementFamily.SHELL, "span_direction": (0.0, 0.0, 1.0), "properties": props}}
    asm = MeshAssembler(sub, model)
    dpn, n_dof = asm.dofs_per_node, asm.dofs_count
    assert dpn == 6
    X = np.array([n.coords for n in sub.nodes], float)
    r, c, v = asm._rust.assemble_k()
    K = sp.coo_matrix((np.asarray(v), (np.asarray(r), np.asarray(c))), shape=(n_dof, n_dof)).tocsr()
    zlo, zhi = X[:, 2].min(), X[:, 2].max()
    w = END_ZONE * (zhi - zlo)
    a_nodes = np.flatnonzero(X[:, 2] <= zlo + w)
    b_nodes = np.flatnonzero(X[:, 2] >= zhi - w)
    free_nodes = np.setdiff1d(np.arange(len(X)), np.concatenate([a_nodes, b_nodes]))
    ref = np.array([X[b_nodes, 0].mean(), X[b_nodes, 1].mean(), X[b_nodes, 2].mean()])
    nf = len(free_nodes)
    fdofs = (free_nodes[:, None] * 6 + np.arange(6)).ravel()
    T_free = sp.coo_matrix((np.ones(len(fdofs)), (fdofs, np.arange(len(fdofs)))), shape=(n_dof, len(fdofs) + 6)).tocsr()
    Tb = rigid_end_map(X[b_nodes], ref).tocoo()
    bdofs = (b_nodes[:, None] * 6 + np.arange(6)).ravel()
    T_b = sp.coo_matrix((Tb.data, (bdofs[Tb.row], len(fdofs) + Tb.col)), shape=(n_dof, len(fdofs) + 6)).tocsr()
    T = (T_free + T_b).tocsr()
    Kr = (T.T @ K @ T).tocsc()
    lu = splu(Kr)
    Us = []
    for dof in (3, 4, 5):  # M_flap about x, M_edge about y, T about z, at the reference node
        f = np.zeros(Kr.shape[0])
        f[len(fdofs) + dof] = M
        u = T @ lu.solve(f)
        Us.append(u.reshape(-1, 6)[:, :3])
    idx = {n.id: i for i, n in enumerate(sub.nodes)}
    area = np.zeros(len(X))
    for e in sub.elements:  # mid-surface tributary area, element area split equally among its nodes
        q = np.array([n.coords for n in e.nodes], float)
        a = 0.5 * np.linalg.norm(np.cross(q[2] - q[0], q[3] - q[1])) if len(q) == 4 else 0.5 * np.linalg.norm(np.cross(q[1] - q[0], q[2] - q[0]))
        for n in e.nodes:
            area[idx[n.id]] += a / len(q)
    return X, Us, area, {"n_nodes": len(X), "n_elements": sub.elements_count, "n_free_nodes": nf, "n_clamped": len(a_nodes),
                   "n_rigid": len(b_nodes), "ref": ref.tolist(), "z_range": [float(zlo), float(zhi)]}


def window_compliance(X, Us, zc, L, loads, weights=None, include=None, windows=None):
    """3x3 compliance for each window of ``windows`` (default ``WINDOWS``: name -> half width / L) from nodal displacements of the three load cases."""
    out = {}
    for name, frac in (WINDOWS if windows is None else windows).items():
        half = frac * L
        fits = [fit_section_kinematics(X if include is None else X[include], U if include is None else U[include], zc, half,
                                       weights=None if weights is None else (weights if include is None else weights[include]))
                for U in Us]
        C = compliance_from_responses([f["kappa"] for f in fits], loads)
        out[name] = {"C": C, "n": fits[0]["n"], "rms": [f["rms"] for f in fits], "distortion": [f["distortion"] for f in fits]}
    return out


def shell_profile(X, Us, zc, L, deck_f, half=2.0, step=2.0, weights=None):
    """Sliding-window readings along the segment (centres zc +- k step, |s| <= L/3) with the deck's exact value at each centre.

    The shell section compliance is NOT smooth in z (layup steps, mesh events): a single window value depends on where it sits,
    so the profile is part of the answer, not a diagnostic afterthought.
    """
    rows = []
    for s in np.arange(-np.floor(L / 3 / step) * step, L / 3 + 1e-9, step):
        try:
            fits = [fit_section_kinematics(X, U, zc + s, half, weights=weights) for U in Us]
        except ValueError:
            continue
        m = section_metrics(compliance_from_responses([f["kappa"] for f in fits], [M0] * 3))
        d = section_metrics(deck_f(zc + s)[0])
        rows.append({"z": float(zc + s), **{k: m[k] for k in ("psi_deg", "ratio_fe_ff", "GJ", "EI_flap", "EI_edge", "sym_err")},
                     "deck": {k: d[k] for k in ("psi_deg", "ratio_fe_ff", "GJ", "EI_flap", "EI_edge")}})
    return rows


def shell_study(zs, out_dir, n_chords=(2.0, 4.0, 6.0, 8.0), fracs=(1 / 12, 1 / 8, 1 / 6, 1 / 4, 1 / 3), chord_divisions=20):
    """Convergence study of the shell measurement (L5 T9c): for each station z, segment length (chords) and fit half window (fraction
    of the segment length) the shell metrics, the deck read by the SAME fit on the SAME nodes (``beam_displacements``) and the deck's
    exact value at the segment centre.  ``bias`` = same-fit deck - exact deck is what the taper / twist / window does to a section whose
    compliance is exactly the deck's; ``corrected`` removes it from the shell reading (see ``bias_corrected``)."""
    from tests.support.v2_blade import build_v2_blade

    t0 = time.time()
    mesh, props = build_v2_blade(0.25, skin_reference="mid", mesher="gmsh", es_range=(0.05, 0.25), chord_divisions=chord_divisions)
    deck_f = deck_C_of_z(load_deck()[0])
    wins = {f"{f:.4f}": f for f in fracs}
    rows = []
    for z in zs:
        for nc in n_chords:
            try:
                sub, p, chord, L = shell_segment(mesh, props, z, nc)
                X, Us, area, info = solve_shell_couples(sub, p)
            except Exception as e:  # noqa: BLE001  (a segment that crosses the root / tip or fails to assemble is a result, not a crash)
                rows.append({"z": z, "n_chord": nc, "error": f"{type(e).__name__}: {e}"})
                print(f"study z={z:g} n={nc:g}: {e}", flush=True)
                continue
            zlo, zhi = info["z_range"]
            zc = 0.5 * (zlo + zhi)
            Lz = zhi - zlo
            ref_xy = X[:, :2].mean(axis=0)
            Ub = [beam_displacements(X, lambda zz: deck_f(zz)[0], zlo, zhi, M0 * np.eye(3)[j], ref_xy) for j in range(3)]
            sh = window_compliance(X, Us, zc, Lz, [M0] * 3, weights=area, windows=wins)
            bd = window_compliance(X, Ub, zc, Lz, [M0] * 3, weights=area, windows=wins)
            exact = section_metrics(deck_f(zc)[0])
            for k, frac in wins.items():
                ms, mb = section_metrics(sh[k]["C"]), section_metrics(bd[k]["C"])
                rows.append({"z": z, "n_chord": nc, "frac": frac, "half_m": frac * Lz, "zc": zc, "L": Lz, "chord": chord, "n_fit": sh[k]["n"],
                             "shell": ms, "deck_same_fit": mb, "deck_exact": exact, "corrected": bias_corrected(ms, mb, exact),
                             "distortion": sh[k]["distortion"]})
            print(f"study z={z:g} n_chord={nc:g}: L={Lz:.1f} m, {info['n_nodes']} nodes ({time.time() - t0:.0f} s)", flush=True)
    doc = {"model": "shell_study", "chord_divisions": chord_divisions, "rows": rows, "elapsed_s": time.time() - t0, "revision": tp.revision(REPO)}
    _write(Path(out_dir) / "shell_study.json", doc)
    return doc


def converged_setting(rows, tol_psi=TOL_PSI, tol_stiff=TOL_STIFF, frac=0.25):
    """Written convergence criterion of the shell measurement: per station z the SMALLEST segment length (chords) such that every later
    step to the next studied length, at the fit window ``frac``, moves the bias-corrected psi by <= ``tol_psi`` deg and GJ, EI_flap, EI_edge
    by <= ``tol_stiff`` (relative).  The longest length has no later step and is never accepted; ``None`` = not converged in the studied range.
    ``rows`` = shell_study rows (``corrected`` dict); error rows are skipped."""
    out = {}
    for z in sorted({float(r["z"]) for r in rows}):
        by_n = {float(r["n_chord"]): r["corrected"] for r in rows
                if float(r["z"]) == z and "corrected" in r and abs(r["frac"] - frac) < 1e-3}
        ns = sorted(by_n)
        steps = [abs(by_n[b]["psi_deg"] - by_n[a]["psi_deg"]) <= tol_psi
                 and all(abs(by_n[b][k] / by_n[a][k] - 1.0) <= tol_stiff for k in ("GJ", "EI_flap", "EI_edge"))
                 for a, b in zip(ns[:-1], ns[1:], strict=True)]
        out[z] = next((ns[i] for i in range(len(steps)) if all(steps[i:])), None)
    return out


def bias_corrected(shell, deck_fit, deck_exact):
    """Shell reading with the geometric bias of the fit removed: psi and ratios by difference, stiffnesses by ratio.

    bias = (deck read through the same fit on the same nodes) vs (deck exact at the segment centre): it is what taper, twist and the
    window do to a section whose compliance is exactly the deck's, so it is common to any model on those nodes.
    """
    out = {"psi_deg": shell["psi_deg"] - (deck_fit["psi_deg"] - deck_exact["psi_deg"]),
           "ratio_fe_ff": shell["ratio_fe_ff"] - (deck_fit["ratio_fe_ff"] - deck_exact["ratio_fe_ff"])}
    if all("ratio_ef_ee" in m for m in (shell, deck_fit, deck_exact)):
        out["ratio_ef_ee"] = shell["ratio_ef_ee"] - (deck_fit["ratio_ef_ee"] - deck_exact["ratio_ef_ee"])
    for k in ("GJ", "EI_flap", "EI_edge"):
        out[k] = shell[k] * deck_exact[k] / deck_fit[k]
    return out


def run_shell(zs, out_dir, chord_divisions=20, n_chord=N_SHELL_CHORD):
    from tests.support.v2_blade import build_v2_blade

    t0 = time.time()
    mesh, props = build_v2_blade(0.25, skin_reference="mid", mesher="gmsh", es_range=(0.05, 0.25),
                                 chord_divisions=chord_divisions)
    deck_f = deck_C_of_z(load_deck()[0])
    rows = []
    for z in zs:
        sub, p, chord, L = shell_segment(mesh, props, z, n_chord)
        X, Us, area, info = solve_shell_couples(sub, p)
        zc = 0.5 * (info["z_range"][0] + info["z_range"][1])
        Lz = info["z_range"][1] - info["z_range"][0]
        wins = window_compliance(X, Us, zc, Lz, [M0] * 3, weights=area)  # primary: tributary-area weights (as the solid)
        wins_eq = window_compliance(X, Us, zc, Lz, [M0] * 3)
        main = wins[SHELL_WINDOW]
        zlo, zhi = info["z_range"]
        ref_xy = X[:, :2].mean(axis=0)
        Ub = [beam_displacements(X, lambda zz: deck_f(zz)[0], zlo, zhi, M0 * np.eye(3)[j], ref_xy) for j in range(3)]
        bwins = window_compliance(X, Ub, zc, zhi - zlo, [M0] * 3, weights=area)
        rows.append(_entry(z, main["C"], {"chord": chord, "L": L, "z_center": zc, **info, "fit_n": main["n"], "fit_rms_m": main["rms"],
                                          "fit_distortion": main["distortion"]},
                           {**{f"{k}_area": section_metrics(v["C"]) | {"n": v["n"], "rms": v["rms"]} for k, v in wins.items()},
                            **{f"{k}_equal": section_metrics(v["C"]) | {"n": v["n"], "rms": v["rms"]} for k, v in wins_eq.items()}}))
        rows[-1]["deck_same_fit"] = {k: section_metrics(v["C"]) for k, v in bwins.items()}
        rows[-1]["corrected"] = bias_corrected(rows[-1], rows[-1]["deck_same_fit"][SHELL_WINDOW], section_metrics(deck_f(zc)[0]))  # the value the table prints
        rows[-1]["profile"] = shell_profile(X, Us, zc, zhi - zlo, lambda zz: deck_f(zz), weights=area)
        rows[-1]["deck_exact_at_zc"] = section_metrics(deck_f(zc)[0])
        rows[-1]["window"] = SHELL_WINDOW
        print(f"shell z={z:g}: {info['n_nodes']} nodes, L={L:.2f} m, psi={rows[-1]['psi_deg']:.3f} deg, "
              f"GJ={rows[-1]['GJ']:.4e} ({time.time() - t0:.0f} s)", flush=True)
    doc = {"model": "shell", "chord_divisions": chord_divisions, "n_chord": n_chord, "moment_Nm": M0, "stations": rows,
           "elapsed_s": time.time() - t0, "revision": tp.revision(REPO)}
    _write(Path(out_dir) / "shell.json", doc)
    return doc


# --- (c) solid --------------------------------------------------------------------------------------------------------


def solid_build(z, level, out, n_chord=N_CHORD, te_filler=0.0, nmid=5):
    import l3_solid_section_ccx as S

    eta = eta_for(z)
    sec = S.blade_section(eta)
    chord = float(np.linalg.norm(sec.le_point - sec.te_point))
    lv = LEVELS[level]
    m2 = S.build_2d(sec, lv["h_rel"] * chord, lv["tmax"], "bonded",
                    {"mat": "adhesive", "arc": te_filler, "nthick": 1} if te_filler else None)
    meta = S.write_model(sec, m2, out, n_chord, lv["hz_rel"] * chord, nmid, [0.25, 1 / 3, 0.4, 0.5, 0.6, 2 / 3, 0.75])
    meta["GJ_st2d"] = S.torsion_2d(sec, m2)  # weight-free St-Venant GJ of the same 2D mesh (the 3D estimator's own continuum limit)
    meta.update(z=z, eta=eta, level=level, h_m=lv["h_rel"] * chord, tmax=lv["tmax"], hz_m=lv["hz_rel"] * chord, te_filler=te_filler,
                n_chord=n_chord)
    json.dump(meta, open(Path(out) / "meta.json", "w"))
    print(f"built solid z={z:g} eta={eta:.4f} level {level}: NQ={meta['NQ']} nz={meta['nz']} elements={meta['nelem']} "
          f"nodes={meta['nnode']} L={meta['L']:.2f} chord={chord:.3f}", flush=True)
    return meta


def solid_st2d(z, level, out, te_filler=0.0):
    """Add the 2D St-Venant GJ to an existing solid directory's meta.json (rebuilds the 2D mesh only; seconds)."""
    import l3_solid_section_ccx as S

    sec = S.blade_section(eta_for(z))
    chord = float(np.linalg.norm(sec.le_point - sec.te_point))
    lv = LEVELS[level]
    m2 = S.build_2d(sec, lv["h_rel"] * chord, lv["tmax"], "bonded", {"mat": "adhesive", "arc": te_filler, "nthick": 1} if te_filler else None)
    meta = json.load(open(Path(out) / "meta.json"))
    meta["GJ_st2d"] = S.torsion_2d(sec, m2)
    json.dump(meta, open(Path(out) / "meta.json", "w"))
    return meta["GJ_st2d"]


def solid_run(out):
    """The three load cases (T, MX, MY) as three parallel single-thread ccx processes; the axial case F is not needed."""
    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
    procs = []
    for c in ("T", "MX", "MY"):
        lg = open(Path(out) / f"ccx_{c}.log", "w")
        procs.append((c, subprocess.Popen([CCX, f"case_{c}"], cwd=out, stdout=lg, stderr=subprocess.STDOUT, env=env)))
    bad = [c for c, p in procs if p.wait() != 0]
    if bad:
        raise SystemExit(f"ccx failed for load cases {bad}; see {out}/ccx_*.log")


def solid_analyze(out, frame_deg=0.0):
    import l3_solid_section_ccx as S

    out = Path(out)
    meta = json.load(open(out / "meta.json"))
    N2, L = meta["N2"], meta["L"]
    xy = np.array(meta["corner_xy"])
    wn = np.array(meta["wnode"])
    blocks = {c: S.parse_dat(out / f"case_{c}.dat", 1)[0] for c in ("MX", "MY", "T")}
    levels = meta["levels"]
    zs = dict(zip(levels, meta["zs_levels"], strict=True))
    X = np.vstack([np.column_stack([xy, np.full(N2, zs[k])]) for k in levels])
    ids = [[1 + k * N2 + i for i in range(N2)] for k in levels]
    Us = [np.array([[blocks[c][i] for i in row] for row in ids]).reshape(-1, 3) for c in ("MX", "MY", "T")]
    W = np.tile(wn, len(levels))
    zc = 0.5 * L
    res = {}
    for wname, weights in (("area", W), ("equal", None)):
        wins = window_compliance(X, Us, zc, L, [M0] * 3, weights=weights)
        for k, v in wins.items():
            C = rotate_compliance(v["C"], rot_z3(frame_deg))
            res[f"{k}_{wname}"] = section_metrics(C) | {"n": v["n"], "rms": v["rms"], "distortion": v["distortion"]}
    main = res["third_area"]
    doc = {"model": "solid", "z": meta["z"], "eta": meta["eta"], "level": meta["level"], "te_filler": meta["te_filler"], "n_chord": meta["n_chord"],
           "frame_deg": frame_deg, "GJ_st2d": meta.get("GJ_st2d"), **{k: main[k] for k in main if k not in ("n", "rms", "distortion")}, "windows": res,
           "mesh": {k: meta[k] for k in ("NQ", "nz", "nelem", "nnode", "L", "chord", "h_m", "tmax", "hz_m")}}
    _write(out / "solid.json", doc)
    return doc


# --- table ------------------------------------------------------------------------------------------------------------


def make_table(out_dir, solid_dirs=()):
    out = Path(out_dir)
    deck = json.loads((out / "deck.json").read_text())["stations"] if (out / "deck.json").exists() else []
    shell = json.loads((out / "shell.json").read_text())["stations"] if (out / "shell.json").exists() else []
    solids = [json.loads((Path(d) / "solid.json").read_text()) for d in solid_dirs if (Path(d) / "solid.json").exists()]
    rows = []
    for z in sorted({r["z"] for r in deck + shell + solids}):
        d = next((r for r in deck if r["z"] == z), None)
        sh = [r for r in shell if r["z"] == z]
        for label, src in [("deck", [d] if d else []),
                           ("shell, fit bias removed", [r["corrected"] | {"z": z} for r in sh if "corrected" in r]),
                           ("shell raw", sh),
                           ("deck, same fit on the shell nodes", [r["deck_same_fit"][r.get("window", "third")] | {"z": z} for r in sh if "deck_same_fit" in r]),
                           *[(f"solid {s['level']}{' +fill' if s['te_filler'] else ''}", [s]) for s in solids if s["z"] == z]]:
            for r in src:
                rows.append((z, label, r, d))
    lines = ["| z | model | psi deg | C_fe/C_ff | C_ef/C_ee | GJ N m2 | EI_flap N m2 | EI_edge N m2 | GJ/deck | flap/deck | edge/deck |",
             "|---|---|---|---|---|---|---|---|---|---|---|"]
    for z, label, r, d in rows:
        rt = (lambda k: f"{r[k] / d[k]:.3f}") if d else (lambda k: "-")
        lines.append(f"| {z:g} | {label} | {r['psi_deg']:.3f} | {r['ratio_fe_ff']:+.4f} | {r.get('ratio_ef_ee', float('nan')):+.4f} | {r['GJ']:.4e} | "
                     f"{r['EI_flap']:.4e} | {r['EI_edge']:.4e} | {rt('GJ')} | {rt('EI_flap')} | {rt('EI_edge')} |")
    for sd in solids:
        if sd.get("GJ_st2d"):
            d = next((r for r in deck if r["z"] == sd["z"]), None)
            lines.append(f"| {sd['z']:g} | solid {sd['level']} 2D St-Venant GJ (weight-free) | - | - | - | {sd['GJ_st2d']:.4e} | - | - | "
                         f"{(sd['GJ_st2d'] / d['GJ']) if d else float('nan'):.3f} | - | - |")
    text = "\n".join(lines)
    (out / "table.md").write_text(text + "\n")
    return text


def study_table(out_dir):
    """Markdown tables of shell_study.json: raw and bias-corrected psi, C_fe/C_ff, GJ, EI_flap, EI_edge per (z, chords, window)."""
    out = Path(out_dir)
    rows = json.loads((out / "shell_study.json").read_text())["rows"]
    lines = ["| z | chords | L m | half/L | psi raw | psi corr | deck psi | fe/ff raw | fe/ff corr | deck fe/ff | GJ corr/deck | flap corr/deck | edge corr/deck | GJ raw/deck | edge raw/deck |",
             "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for r in rows:
        if "error" in r:
            lines.append(f"| {r['z']:g} | {r['n_chord']:g} | - | - | {r['error']} |")
            continue
        s, c, d = r["shell"], r["corrected"], r["deck_exact"]
        lines.append(f"| {r['z']:g} | {r['n_chord']:g} | {r['L']:.1f} | {r['frac']:.3f} | {s['psi_deg']:.2f} | {c['psi_deg']:.2f} | {d['psi_deg']:.2f} | "
                     f"{s['ratio_fe_ff']:+.4f} | {c['ratio_fe_ff']:+.4f} | {d['ratio_fe_ff']:+.4f} | {c['GJ'] / d['GJ']:.3f} | {c['EI_flap'] / d['EI_flap']:.3f} | "
                     f"{c['EI_edge'] / d['EI_edge']:.3f} | {s['GJ'] / d['GJ']:.3f} | {s['EI_edge'] / d['EI_edge']:.3f} |")
    text = "\n".join(lines)
    (out / "study_table.md").write_text(text + "\n")
    return text


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("model", choices=["deck", "shell", "solid", "table", "study", "study-table"])
    ap.add_argument("--z", type=float, nargs="+", default=list(STATIONS))
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--cmd", choices=["build", "run", "analyze", "all", "st2d"], default="all", help="solid only")
    ap.add_argument("--level", choices=sorted(LEVELS), default="A", help="solid mesh level (relative to the local chord)")
    ap.add_argument("--nchord", type=float, default=None, help=f"segment length in chords (shell {N_SHELL_CHORD:g}, solid {N_CHORD:g})")
    ap.add_argument("--te-filler", dest="te_filler", type=float, default=0.0, help="solid: TE adhesive filler arc, m (L23 factor test)")
    ap.add_argument("--nchords", type=float, nargs="+", default=[2.0, 4.0, 6.0, 8.0], help="study: segment lengths in chords")
    ap.add_argument("--fracs", type=float, nargs="+", default=[1 / 12, 1 / 8, 1 / 6, 1 / 4, 1 / 3], help="study: fit half window / segment length")
    ap.add_argument("--chord-divisions", type=int, default=20)
    ap.add_argument("--solid-dirs", nargs="*", default=[])
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    if a.model == "deck":
        for r in run_deck(a.z, a.out)["stations"]:
            print(f"deck z={r['z']:g}: psi={r['psi_deg']:.3f} fe/ff={r['ratio_fe_ff']:+.4f} GJ={r['GJ']:.4e}")
    elif a.model == "shell":
        run_shell(a.z, a.out, a.chord_divisions, a.nchord or N_SHELL_CHORD)
    elif a.model == "study":
        shell_study(a.z, a.out, a.nchords, a.fracs, a.chord_divisions)
    elif a.model == "study-table":
        print(study_table(a.out))
    elif a.model == "solid":
        z = a.z[0]
        if a.cmd in ("build", "all"):
            solid_build(z, a.level, a.out, a.nchord or N_CHORD, a.te_filler)
        if a.cmd == "st2d":
            print(f"GJ_st2d = {solid_st2d(z, a.level, a.out, a.te_filler):.5e}")
        if a.cmd in ("run", "all"):
            solid_run(a.out)
        if a.cmd in ("analyze", "all"):
            r = solid_analyze(a.out)
            print(f"solid z={z:g} level {a.level}: psi={r['psi_deg']:.3f} fe/ff={r['ratio_fe_ff']:+.4f} GJ={r['GJ']:.4e} "
                  f"EIf={r['EI_flap']:.4e} EIe={r['EI_edge']:.4e} sym_err={r['sym_err']:.3f} GJ_st2d={r['GJ_st2d']:.4e}")
    else:
        print(make_table(a.out, a.solid_dirs))


if __name__ == "__main__":
    main()
