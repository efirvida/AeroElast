"""P1 twist attribution: pure functions and the CONTROLS stage (J1).

Dispatched from ``tools/diagnose_radii_beamdyn.py --twist-p1 {controls,matrix}`` under the
2026-10-23 re-pre-registration (``odd/tasks/twist-discrepancy-closure.md``, R1-R8).  Everything
above ``# --- BeamDyn harness`` is pure numpy and is unit-tested by ``tests/software/test_twist_p1.py``.

Conventions (R1): angles are right-hand rotations about +z (span), nose-up positive.  BeamDyn's blade
frame is x_BD = mesh +y, y_BD = -mesh x, z_BD = z, a proper rotation about z, so RDzr keeps the
sense of the mesh theta_z.  Primary section angle (R2): phi_z = 4 atan(c_z / 4) of the WM component;
secondary: the twist about the deformed tangent, from the swing-twist split.

Controls (R3), each reading a different row, as registered:

* C1  static (DynamicSolve False), pure tip couple Tz = 1e4 N.m, last row, against
      ``integrate_twist`` (the deck's own 6x6 compliance, couple transported into each section frame).
* C2a dynamic, QuasiStaticInit, Omega = 7.56 rpm about x_BD, GlbPos = (0, 0, 3.97), no load, FIRST row
      (its time is printed): RootFzr against ``centrifugal_root_force``.
* C2b dynamic, QuasiStaticInit, Omega = 0, P_shell load (0.25 m strips) at scale 1, FIRST row,
      against the static solve of the same load: tip phi_z and tip flap within 0.1 %.

Exit status of the controls stage: 0 all passed, 3 C1 failed (fatal for P1), 2 only C2 failed.
"""

from __future__ import annotations

import hashlib
import json
import os
import re
import subprocess
import sys
from pathlib import Path

import numpy as np

# mesh vector components -> BeamDyn blade-frame components (R1)
MESH_TO_BD = np.array([[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])

OMEGA_RATED = 7.56 * 2.0 * np.pi / 60.0  # rad/s
HUB_RADIUS = 3.97  # m, ElastoDyn HubRad (bd_mu30/ed.dat:47)
TIP_COUPLE = 1.0e4  # N.m, C1
FRACTIONS = (0.5, 0.75, 0.9, 1.0)  # C1 stations: 50, 75, 90 % and the tip
C1_SUBDIVISIONS = 1024  # reference points per blade-station interval


# --- rotations ---------------------------------------------------------------------------------


def wm_z_angle(c_z):
    """Section angle from the z Wiener-Milenkovic component: phi = 4 atan(c / 4) (R2 primary)."""
    return 4.0 * np.arctan(np.asarray(c_z, float) / 4.0)


def wm_to_matrix(c):
    """Rotation matrix of the Wiener-Milenkovic vector ``c = 4 tan(phi/4) n``."""
    c = np.asarray(c, float)
    norm = np.linalg.norm(c)
    if norm == 0.0:
        return np.eye(3)
    return rotation_matrix(c / norm, 4.0 * np.arctan(norm / 4.0))


def rotation_matrix(axis, angle):
    """Right-handed Rodrigues rotation about the unit vector ``axis``."""
    a = np.asarray(axis, float)
    k = np.array([[0.0, -a[2], a[1]], [a[2], 0.0, -a[0]], [-a[1], a[0], 0.0]])
    return np.eye(3) + np.sin(angle) * k + (1.0 - np.cos(angle)) * k @ k


def rot_to_quat(R):
    """Unit quaternion (w, x, y, z), w >= 0, of a rotation matrix (Shepperd)."""
    R = np.asarray(R, float)
    tr = np.trace(R)
    if tr > 0.0:
        s = 2.0 * np.sqrt(1.0 + tr)
        q = np.array([0.25 * s, (R[2, 1] - R[1, 2]) / s, (R[0, 2] - R[2, 0]) / s, (R[1, 0] - R[0, 1]) / s])
    else:
        i = int(np.argmax(np.diag(R)))
        j, k = (i + 1) % 3, (i + 2) % 3
        s = 2.0 * np.sqrt(1.0 + R[i, i] - R[j, j] - R[k, k])
        q = np.empty(4)
        q[0] = (R[k, j] - R[j, k]) / s
        q[1 + i] = 0.25 * s
        q[1 + j] = (R[j, i] + R[i, j]) / s
        q[1 + k] = (R[k, i] + R[i, k]) / s
    return q if q[0] >= 0.0 else -q


def minimal_rotation(a, b):
    """The shortest rotation taking direction ``a`` to direction ``b`` (no twist about either)."""
    a = np.asarray(a, float) / np.linalg.norm(a)
    b = np.asarray(b, float) / np.linalg.norm(b)
    cross = np.cross(a, b)
    sin, cos = np.linalg.norm(cross), float(a @ b)
    if sin < 1e-15:
        if cos > 0.0:
            return np.eye(3)
        raise ValueError("antiparallel directions: the minimal rotation is not unique")
    return rotation_matrix(cross / sin, np.arctan2(sin, cos))


def swing_twist(R, axis):
    """Split ``R = R_twist(deformed axis) . R_swing`` (R2 secondary).

    ``axis`` is the undeformed tangent.  Returns ``(twist, swing)`` in radians: the signed twist about
    the deformed tangent ``R @ axis`` and the unsigned swing angle that tilts the tangent.
    """
    t0 = np.asarray(axis, float) / np.linalg.norm(axis)
    t_def = R @ t0
    S = minimal_rotation(t0, t_def)
    q = rot_to_quat(R @ S.T)
    twist = 2.0 * np.arctan2(q[1:] @ t_def, q[0])
    swing = np.arctan2(np.linalg.norm(np.cross(t0, t_def)), float(t0 @ t_def))
    return float(twist), float(swing)


# --- sections: C1 and C2a references -----------------------------------------------------------


def cumtrapz(y, x):
    y = np.asarray(y, float)
    return np.concatenate([[0.0], np.cumsum(0.5 * (y[1:] + y[:-1]) * np.diff(x))])


def section_frames(kp, twist_deg):
    """Section frames (columns = local axes in the blade frame) on the key-point curve.

    The local z axis is the curve tangent; the frame is the minimal tilt from +z followed by the
    initial twist about the tangent (right-handed).
    """
    kp = np.asarray(kp, float)
    arc = np.concatenate([[0.0], np.cumsum(np.linalg.norm(np.diff(kp, axis=0), axis=1))])
    tangent = np.gradient(kp, arc, axis=0)
    z = np.array([0.0, 0.0, 1.0])
    return np.array([
        minimal_rotation(z, t) @ rotation_matrix(z, np.deg2rad(tw)) for t, tw in zip(tangent, twist_deg, strict=True)
    ])


def transport_moment(R, M_global):
    """Components of the constant global couple in each section frame: ``R_i^T M``."""
    return np.einsum("nji,j->ni", np.asarray(R, float), np.asarray(M_global, float))


def integrate_twist(arc, K, R, M_global):
    """C1 reference: ``twist(s) = int_0^s (K(z)^-1 [0, 0, 0, M_local])_6 dz`` in radians.

    ``K`` is (n, 6, 6) in the section frame at arc positions ``arc`` (n,), ``R`` (n, 3, 3) the section
    frames; the couple is constant in the global frame and carries no force.
    """
    local = transport_moment(R, M_global)
    rhs = np.concatenate([np.zeros((len(local), 3)), local], axis=1)[..., None]
    rate = np.linalg.solve(np.asarray(K, float), rhs)[:, 5, 0]
    return cumtrapz(rate, arc)


def rotor_axis_radius(kp, hub_radius):
    """Distance from the rotor axis (x_BD through the origin) of key points of a blade rooted at z = hub."""
    kp = np.asarray(kp, float)
    return np.hypot(kp[:, 1], hub_radius + kp[:, 2])


def centrifugal_root_force(s, m, rho, omega):
    """C2a reference: ``int m(s) Omega^2 rho(s) ds`` along the arc ``s``."""
    return float(np.trapezoid(np.asarray(m, float) * omega**2 * np.asarray(rho, float), s))


def within(value, ref, tol):
    """``value`` within relative ``tol`` of ``ref`` and of the same sign (a zero reference never passes)."""
    return bool(ref != 0.0 and value * ref > 0.0 and abs(value / ref - 1.0) <= tol)


# --- BeamDyn harness -----------------------------------------------------------------------------


def revision(repo):
    try:
        out = subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=repo, capture_output=True,
                             text=True, check=True)
        return out.stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return os.environ.get("REV", "unknown")


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def parse_blade(path):
    """Station fractions (n,), stiffness (n, 6, 6) and mass (n, 6, 6) of a BeamDyn blade file."""
    lines = Path(path).read_text().splitlines()
    k0 = next(i for i, ln in enumerate(lines) if "DISTRIBUTED PROPERTIES" in ln) + 1
    rows = [ln.split() for ln in lines[k0:] if ln.strip()]
    eta, K, M = [], [], []
    i = 0
    while i < len(rows):
        if len(rows[i]) != 1:
            break
        eta.append(float(rows[i][0]))
        K.append([[float(v) for v in r] for r in rows[i + 1:i + 7]])
        M.append([[float(v) for v in r] for r in rows[i + 7:i + 13]])
        i += 13
    return np.array(eta), np.array(K), np.array(M)


class Deck:
    """One converted BeamDyn deck: key points, quadrature points and section properties.

    The quadrature points are the blade stations refined by two (the deck sets ``refine 2``), as in
    ``diagnose_radii_beamdyn.read_beamdyn``, placed on the key-point curve by arc-length fraction.
    """

    def __init__(self, kp, blade_path, twist_deg):
        self.kp = np.asarray(kp, float)
        seg = np.linalg.norm(np.diff(self.kp, axis=0), axis=1)
        self.arc_kp = np.concatenate([[0.0], np.cumsum(seg)])
        self.length = float(self.arc_kp[-1])
        st, K, M = parse_blade(blade_path)
        self.eta = np.append(np.ravel(np.column_stack([st[:-1], 0.5 * (st[:-1] + st[1:])])), st[-1])
        self.K = np.array([_interp_matrix(st, K, e) for e in self.eta])
        self.M = np.array([_interp_matrix(st, M, e) for e in self.eta])
        frac = self.arc_kp / self.length
        self.qp = np.column_stack([np.interp(self.eta, frac, self.kp[:, a]) for a in range(3)])
        self.arc = self.eta * self.length
        frames = section_frames(self.kp, twist_deg)
        tangent = np.gradient(self.kp, self.arc_kp, axis=0)
        # frames and tangents are interpolated per component and re-orthonormalised
        self.tangent = _unit(np.column_stack([np.interp(self.eta, frac, tangent[:, a]) for a in range(3)]))
        self.R = np.array([_orthonormal(np.array([[np.interp(e, frac, frames[:, i, j]) for j in range(3)]
                                                  for i in range(3)])) for e in self.eta])
        self._stations, self._K, self._frac, self._frames = st, K, frac, frames

    def dense(self, n_sub):
        """Section data on a grid with ``n_sub`` points per blade-station interval.

        The 1/K integrand falls off a cliff outboard (K66 drops ~190x in the last 5 %), so the C1
        reference is evaluated here, not on the 51 quadrature points, where trapezoid error alone is 4x.
        K is linear between stations (BeamDyn's own interpolation) and the frames follow the key points.
        """
        st = self._stations
        eta = np.unique(np.concatenate([np.linspace(a, b, n_sub + 1) for a, b in zip(st[:-1], st[1:], strict=True)]))
        K = np.empty((len(eta), 6, 6))
        for i in range(6):
            for j in range(6):
                K[:, i, j] = np.interp(eta, st, self._K[:, i, j])
        R = np.empty((len(eta), 3, 3))
        for i in range(3):
            for j in range(3):
                R[:, i, j] = np.interp(eta, self._frac, self._frames[:, i, j])
        u, _, vt = np.linalg.svd(R)
        return eta, eta * self.length, K, u @ vt


def _interp_matrix(st, mats, e):
    return np.array([[np.interp(e, st, mats[:, i, j]) for j in range(6)] for i in range(6)])


def _unit(v):
    return v / np.linalg.norm(v, axis=1, keepdims=True)


def _orthonormal(A):
    u, _, vt = np.linalg.svd(A)
    return u @ vt


def tune_deck(work, with_root_fz=True):
    """R2: OutFmt ES16.8E3 and the rotation channels in the primary file written by ``write_deck``."""
    path = Path(work) / "bd_primary.dat"
    text = path.read_text()
    text, n = re.subn(r'^"ES10\.3E2"(\s+OutFmt)', r'"ES16.8E3"\1', text, flags=re.M)
    if n != 1:
        raise RuntimeError(f"OutFmt line not found in {path}")
    if "RootFxr, RootFyr, RootFzr" not in text:
        raise RuntimeError(f"root reaction channels missing in {path}")
    for ch in ("RDxr", "RDyr", "RDzr"):
        if f'"{ch}"' not in text:
            text = text.replace('"TDzr"\n', f'"TDzr"\n"{ch}"\n', 1)
    path.write_text(text)


def read_keypoints(primary_text):
    """Key points (50, 4): x, y, z in m and the initial twist in degrees."""
    lines = primary_text.splitlines()
    i0 = next(i for i, ln in enumerate(lines) if "kp_xr" in ln) + 2
    return np.array([[float(v) for v in lines[i].split()[:4]] for i in range(i0, i0 + 50)])


def twist_deg_of(primary_text):
    return read_keypoints(primary_text)[:, 3]


def write_driver(case, deck, dynamic, omega, rows, t_final, dt):
    """Write ``bd_driver.inp``.  ``rows``: (n, 7) array ``eta Fx Fy Fz Mx My Mz`` in the BeamDyn frame."""
    rows = np.zeros((0, 7)) if rows is None else np.asarray(rows, float)
    glb_z = HUB_RADIUS if dynamic else 0.0
    zeros = "\n".join(f"0   {name}" for name in
                      [f"DistrLoad({i})" for i in range(1, 7)] + [f"TipLoad({i})" for i in range(1, 7)])
    table = "".join(f"{r[0]:.8f} " + " ".join(f"{v:.8e}" for v in r[1:]) + "\n" for r in rows)
    (case / "bd_driver.inp").write_text(f"""------- BEAMDYN Driver INPUT FILE --------------------------------
IEA-15 twist P1 control (tools/twist_p1.py)
---------------------- SIMULATION CONTROL --------------------------------------
{'True' if dynamic else 'False'}   DynamicSolve
0   t_initial
{t_final:g}   t_final
{dt:g}   dt
---------------------- GRAVITY PARAMETER --------------------------------------
0   Gx
0   Gy
0   Gz
---------------------- FRAME PARAMETER --------------------------------------
0   GlbPos(1)
0   GlbPos(2)
{glb_z:g}   GlbPos(3)
---The following 3 by 3 matrix is the direction cosine matrix, GlbDCM(3,3),
---relates global frame to the initial blade root frame
1.0 0.0 0.0
0.0 1.0 0.0
0.0 0.0 1.0
True   GlbRotBladeT0
---------------------- ROOT VELOCITY PARAMETER ----------------------------------
{omega:.12g}   RootVel(4)
0   RootVel(5)
0   RootVel(6)
---------------------- APPLIED FORCE ----------------------------------
{zeros}
{len(rows)}   NumPointLoads
eta Fx Fy Fz Mx My Mz
(-) (N) (N) (N) (N-m) (N-m) (N-m)
{table}---------------------- PRIMARY INPUT FILE --------------------------------------
"bd_primary.dat"    InputFile
---------------------- OUTPUTS --------------------------------------
0   WrVTK
1   VTK_fps
""")


def read_out(path):
    """``(names, data)`` of a BeamDyn driver output: ``data`` has one row per time step."""
    lines = Path(path).read_text().splitlines()
    ih = next(i for i, ln in enumerate(lines) if ln.startswith("Time"))
    names = [h.strip() for h in lines[ih].split("\t")]
    data = []
    for ln in lines[ih + 2:]:
        parts = ln.split()
        if len(parts) == len(names):
            data.append([float(p) for p in parts])
    return names, np.array(data)


def prismatic_blade(path):
    """Rewrite a blade file with the root section's K and M at every station (harness-validity diagnostic)."""
    st, K, M = parse_blade(path)
    lines = Path(path).read_text().splitlines()
    k0 = next(i for i, ln in enumerate(lines) if "DISTRIBUTED PROPERTIES" in ln) + 1
    fmt = lambda A: ["\t " + " \t ".join(f"{v:.16e}" for v in r) for r in A]  # noqa: E731
    out = lines[:k0]
    for e in st:
        out += [f"\t {e:.6f}", *fmt(K[0]), "", *fmt(M[0]), ""]
    Path(path).write_text("\n".join(out) + "\n")


def run_case(ctx, name, *, dynamic, omega=0.0, rows=None, t_final=None, dt=0.01, prismatic=False):
    """Run one BeamDyn case under ``<work>/<name>/``; return the loaded deck and the output table."""
    case = ctx["work"] / name
    case.mkdir(parents=True, exist_ok=True)
    kp = ctx["dmr"].write_deck(ctx["deck_dir"], case)
    tune_deck(case)
    if prismatic:
        prismatic_blade(case / "IEA-15-240-RWT_BeamDyn_blade.dat")
    write_driver(case, kp, dynamic, omega, rows, t_final if t_final is not None else (dt if dynamic else 1.0), dt)
    run = subprocess.run([ctx["beamdyn"], "bd_driver.inp"], cwd=case, capture_output=True, text=True)
    if run.returncode != 0 or not (case / "bd_driver.out").exists():
        raise RuntimeError(f"beamdyn_driver failed in {case}:\n{(run.stdout + run.stderr)[-2000:]}")
    names, data = read_out(case / "bd_driver.out")
    twist = twist_deg_of((case / "bd_primary.dat").read_text())
    deck = Deck(kp, case / "IEA-15-240-RWT_BeamDyn_blade.dat", twist)
    ctx["hashes"][name] = {f: sha256(case / f) for f in
                           ("bd_driver.inp", "bd_primary.dat", "IEA-15-240-RWT_BeamDyn_blade.dat")}
    return deck, names, data


def pick_row(data, dynamic):
    """The first row of a dynamic run (the quasi-static initialisation), the last of a static one."""
    return data[0] if dynamic else data[-1]


def node_fields(names, row, n):
    col = {h: i for i, h in enumerate(names)}
    td = np.array([[row[col[f"N{j + 1:03d}_TD{a}r"]] for a in "xyz"] for j in range(n)])
    rd = np.array([[row[col[f"N{j + 1:03d}_RD{a}r"]] for a in "xyz"] for j in range(n)])
    return td, rd


def section_angles(deck, rd):
    """Primary phi_z (R2) and the tangent twist (secondary) at every quadrature point."""
    phi = wm_z_angle(rd[:, 2])
    tang = np.array([swing_twist(wm_to_matrix(c), t)[0] for c, t in zip(rd, deck.tangent, strict=True)])
    return phi, tang


# --- controls ------------------------------------------------------------------------------------


def control_c1(ctx):
    rows = np.array([[1.0, 0, 0, 0, 0, 0, TIP_COUPLE]])
    deck, names, data = run_case(ctx, "C1", dynamic=False, rows=rows)
    row = pick_row(data, False)
    _, rd = node_fields(names, row, len(deck.eta))
    phi, tang = section_angles(deck, rd)
    eta_d, arc_d, K_d, R_d = deck.dense(C1_SUBDIVISIONS)
    ref = integrate_twist(arc_d, K_d, R_d, [0.0, 0.0, TIP_COUPLE])
    _, arc_c, K_c, R_c = deck.dense(C1_SUBDIVISIONS // 4)  # convergence of the reference itself
    tip_coarse = float(integrate_twist(arc_c, K_c, R_c, [0.0, 0.0, TIP_COUPLE])[-1])
    out = []
    for f in FRACTIONS:
        bd, rf = float(np.interp(f, deck.eta, phi)), float(np.interp(f, eta_d, ref))
        out.append({"eta": f, "bd_phi_z": bd, "bd_tangent_twist": float(np.interp(f, deck.eta, tang)),
                        "reference": rf, "naive_T_over_K66": float(np.interp(f, deck.eta, _naive(deck))),
                        "rel_err": bd / rf - 1.0 if rf else None, "ok": within(bd, rf, 0.01)})
    return {"time": float(row[0]), "couple_Nm": TIP_COUPLE, "stations": out, "passed": all(r["ok"] for r in out),
                "reference_tip_at_quarter_subdivisions": tip_coarse, "subdivisions": C1_SUBDIVISIONS}


def diagnostic_c1_prismatic(ctx):
    """NOT a gate.  C1 on a prismatic blade (root K everywhere): separates a wrong harness (frames, WM,
    units) from BeamDyn's single order-10 element failing to resolve the outboard stiffness collapse."""
    rows = np.array([[1.0, 0, 0, 0, 0, 0, TIP_COUPLE]])
    deck, names, data = run_case(ctx, "C1_prismatic", dynamic=False, rows=rows, prismatic=True)
    _, rd = node_fields(names, pick_row(data, False), len(deck.eta))
    phi, _ = section_angles(deck, rd)
    eta_d, arc_d, K_d, R_d = deck.dense(16)
    ref = integrate_twist(arc_d, K_d, R_d, [0.0, 0.0, TIP_COUPLE])
    st = [{"eta": f, "bd_phi_z": float(np.interp(f, deck.eta, phi)), "reference": float(np.interp(f, eta_d, ref))}
          for f in FRACTIONS]
    for r in st:
        r["rel_err"] = r["bd_phi_z"] / r["reference"] - 1.0
    return {"stations": st, "max_abs_rel_err": max(abs(r["rel_err"]) for r in st)}


def _naive(deck):  # T/K66 integrated at the quadrature points: printed for context, not a reference
    return cumtrapz(TIP_COUPLE / deck.K[:, 5, 5], deck.arc)


def control_c2a(ctx):
    deck, names, data = run_case(ctx, "C2a", dynamic=True, omega=OMEGA_RATED)
    row = pick_row(data, True)
    value = float(row[names.index("RootFzr")])
    rho = rotor_axis_radius(deck.qp, HUB_RADIUS)
    ref = centrifugal_root_force(deck.arc, deck.M[:, 0, 0], rho, OMEGA_RATED)
    print(f"C2a first output row at t = {row[0]:g} s (rows: {len(data)}, t_last {data[-1, 0]:g})", flush=True)
    return {"time": float(row[0]), "omega": OMEGA_RATED, "RootFzr": value, "reference": ref,
                "rel_err": abs(value) / ref - 1.0, "sign_of_RootFzr": float(np.sign(value)),
                "passed": bool(abs(abs(value) / ref - 1.0) <= 0.01)}


def shell_load_rows(ctx):
    """P_shell: the production rated rigid-blade nodal load at 0.25 m, binned to BeamDyn point loads."""
    s = shell_state_025(ctx["dmr"])
    kp = read_keypoints((ctx["deck_dir"] / "IEA-15-240-RWT_BeamDyn.dat").read_text())
    arc = np.concatenate([[0.0], np.cumsum(np.linalg.norm(np.diff(kp, axis=0), axis=1))])
    eta = np.interp(s["ref"][:, 2], kp[:, 2], arc / arc[-1])
    F, M = s["F"] @ MESH_TO_BD.T, s["M"] @ MESH_TO_BD.T
    return np.column_stack([eta, F, M]), s


def shell_state_025(dmr):
    """``diagnose_radii_beamdyn.shell_state`` at element_size 0.25 (never its 1.0 m default, AGENTS rule 9).

    Same construction and binning; the linear solve is left out (the controls need only the load).
    """
    dmr.Node._id_counter = 0
    dmr.MeshElement._id_counter = 0
    model = dmr.Blade(str(dmr.fx.YAML), element_size=0.25)
    model.generate_mesh()
    mesh = model.mesh
    X = np.array([[nd.x, nd.y, nd.z] for nd in mesh.nodes])
    part = dmr.BEMFSIParticipant(mesh, dmr.build_blade_aero_from_aerodyn(dmr.fx.AD_PRIMARY), dmr.fx.BEM_CONFIG)
    fa, _ = part._compute_forces(np.zeros((len(X), 3)))
    ref = np.stack([part._ref_coords[i].mean(axis=0) for i in part._strip_node_indices])
    k_of = np.argmin(abs(X[:, 2][:, None] - ref[None, :, 2]), axis=1)
    F, M = np.zeros((len(ref), 3)), np.zeros((len(ref), 3))
    for k in range(len(ref)):
        m = k_of == k
        F[k] = fa[m].sum(0)
        M[k] = np.cross(X[m] - np.array([0.0, ref[k, 1], ref[k, 2]]), fa[m]).sum(0)
    return {"ref": ref, "F": F, "M": M, "n_nodes": len(X)}


def control_c2b(ctx):
    rows, s = shell_load_rows(ctx)
    ctx["load_rows"] = rows
    deck, names, data = run_case(ctx, "C2b_static", dynamic=False, rows=rows)
    row_s = pick_row(data, False)
    td_s, rd_s = node_fields(names, row_s, len(deck.eta))
    _, names_q, data_q = run_case(ctx, "C2b_qs", dynamic=True, omega=0.0, rows=rows)
    row_q = pick_row(data_q, True)
    td_q, rd_q = node_fields(names_q, row_q, len(deck.eta))
    print(f"C2b quasi-static first output row at t = {row_q[0]:g} s", flush=True)
    phi_s, phi_q = float(wm_z_angle(rd_s[-1, 2])), float(wm_z_angle(rd_q[-1, 2]))
    flap_s, flap_q = float(td_s[-1, 0]), float(td_q[-1, 0])
    err_phi = phi_q / phi_s - 1.0 if phi_s else None
    err_flap = flap_q / flap_s - 1.0
    return {"time_qs": float(row_q[0]), "time_static": float(row_s[0]), "n_strips": len(rows),
                "static": {"tip_phi_z": phi_s, "tip_flap": flap_s}, "quasi_static": {"tip_phi_z": phi_q, "tip_flap": flap_q},
                "rel_err_phi": err_phi, "rel_err_flap": err_flap,
                "passed": bool(err_phi is not None and abs(err_phi) <= 1e-3 and abs(err_flap) <= 1e-3)}


def run_controls(args, dmr):
    repo = Path(__file__).resolve().parents[1]
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    work = Path(args.work).resolve()
    work.mkdir(parents=True, exist_ok=True)
    ctx = {"work": work, "deck_dir": Path(args.deck), "beamdyn": args.beamdyn, "dmr": dmr, "hashes": {}}
    src = {f: sha256(Path(args.deck) / f) for f in ("IEA-15-240-RWT_BeamDyn.dat", "IEA-15-240-RWT_BeamDyn_blade.dat")}
    results = {}
    for name, fn in (("C1", control_c1), ("C1_diag_prismatic", diagnostic_c1_prismatic),
                     ("C2a", control_c2a), ("C2b", control_c2b)):
        try:
            results[name] = fn(ctx)
        except Exception as exc:  # a failed control is a result, not a crash
            results[name] = {"passed": False, "error": f"{type(exc).__name__}: {exc}"}
        verdict = "diagnostic" if name.startswith("C1_diag") else ("PASS" if results[name].get("passed") else "FAIL")
        print(f"{name}: {verdict}", flush=True)
        print(json.dumps(results[name], indent=1, default=float)[:1500], flush=True)
    doc = {"revision": revision(repo), "element_size_shell": 0.25, "source_deck_sha256": src,
               "generated_sha256": ctx["hashes"], "controls": results,
               "c1_passed": results["C1"]["passed"],
               "all_passed": all(r["passed"] for k, r in results.items() if not k.startswith("C1_diag"))}
    (out_dir / "controls.json").write_text(json.dumps(doc, indent=1, default=float) + "\n")
    print(f"wrote {out_dir / 'controls.json'}", flush=True)
    if not doc["c1_passed"]:
        return 3
    return 0 if doc["all_passed"] else 2


def main(args, dmr):
    if args.output_dir is None:
        sys.exit("--twist-p1 requires --output-dir")
    if args.work is None:
        sys.exit("--twist-p1 requires --work")
    if args.twist_p1 == "controls":
        return run_controls(args, dmr)
    print("twist-p1 matrix stage: not implemented in this unit (controls only)", file=sys.stderr)
    return 2
