"""BeamDyn reference for the deformed BEM radii (#33 L5 T2, ``radii_linear_kinematics_axial_term``).

Solves the de-loading fixture's state (``tests/validation/blade/test_blade_deloading_vs_reference.py``:
official IEA-15MW shell, ``assemble_k()`` only, rigid-blade rated aero load), bins the nodal aero load
into one BeamDyn point load per BEM strip, and runs a static, geometrically exact BeamDyn solve of the
official IEA-15 BeamDyn deck under it (dead loads, no gravity, no rotation).  It then reports:

1. **The formula on the reference's own field.**  ``lin + nl`` with ``lin = -int (x0' u_x' + y0' u_y')``
   and ``nl = -1/2 int (u_x'^2 + u_y'^2)`` evaluated on BeamDyn's nodal displacements, against
   BeamDyn's own axial displacement.  This separates "is the kinematic formula right" from "shell vs
   beam stiffness".  Truncation bound: the next term of the inextensibility expansion is O(theta^2)
   relative to ``nl`` (theta the deflection slope), printed as ``theta_max^2 / 4``.
2. **The shell against BeamDyn** at the tip-strip centroid span: first-order (old) feed, the
   participant's ``r_def`` increment, and BeamDyn's axial displacement, with the flap ratio.

Frames: BeamDyn blade x = mesh +y (flap, prebend -4 m at the tip in both), BeamDyn y = -mesh x,
z = span in both.  The v4 deck is converted to the v5 input format on the fly (pitch-actuator block
dropped, damping block re-headed, ``damp_flag = 0``: damping is irrelevant to a static solve).

Usage::

    scripts/aeroenv.sh python -W ignore tools/diagnose_radii_beamdyn.py \
        --deck <dir with IEA-15-240-RWT_BeamDyn.dat> --beamdyn <path to beamdyn_driver>

Defaults: ``$IEA15_OPENFAST_DIR`` and ``$BEAMDYN_DRIVER``.  Recorded 2026-10-11 (OpenFAST v5.0.0):
BeamDyn tip axial -1.339 m (at z0 = 116.403 m: -1.3146 m); formula on BeamDyn's field -1.3313 m
(error 0.3 % of nl); shell feed old +1.1056 m, new -1.0478 m; shell flap 17.007 m vs BeamDyn 18.720 m.
"""

from __future__ import annotations

import argparse
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from aeroelast.core.mesh.entities import MeshElement, Node  # noqa: E402
from aeroelast.core.assembler import MeshAssembler  # noqa: E402
from aeroelast.elements import ElementFamily  # noqa: E402
from aeroelast.models.blade.model import Blade  # noqa: E402
from aeroelast.solvers.bem.fsi_participant import BEMFSIParticipant  # noqa: E402
import tests.validation.blade.test_blade_deloading_vs_reference as fx  # noqa: E402
from tests.support.openfast_bem import build_blade_aero_from_aerodyn  # noqa: E402

N_KP = 50  # key points of the official IEA-15 BeamDyn deck


def cumtrapz(y, x):
    return np.concatenate([[0.0], np.cumsum(0.5 * (y[1:] + y[:-1]) * np.diff(x))])


def kinematics(X, U):
    """``lin`` and ``nl`` along a line X (n, 3), span along z, displacements U (n, 3)."""
    z = X[:, 2]
    g = lambda f: np.gradient(f, z)  # noqa: E731
    lin = -cumtrapz(g(X[:, 0]) * g(U[:, 0]) + g(X[:, 1]) * g(U[:, 1]), z)
    nl = -0.5 * cumtrapz(g(U[:, 0]) ** 2 + g(U[:, 1]) ** 2, z)
    return lin, nl


def shell_state():
    Node._id_counter = 0
    MeshElement._id_counter = 0
    model = Blade(str(fx.YAML), element_size=1.0)
    model.generate_mesh()
    mesh = model.mesh
    asm = MeshAssembler(mesh, {"elements": {"element_family": ElementFamily.SHELL,
                                            "span_direction": (0.0, 0.0, 1.0),
                                            "properties": model.get_element_properties()}})
    n, dpn = asm.dofs_count, asm.dofs_per_node
    r, c, v = asm._rust.assemble_k()
    K = coo_matrix((np.asarray(v), (np.asarray(r), np.asarray(c))), shape=(n, n)).tocsr()
    root = {mesh.node_id_to_index[nid] for nid in mesh.get_node_set("RootNodes").node_ids}
    free = np.array(sorted(set(range(n)) - {dpn * i + d for i in root for d in range(dpn)}))
    X = np.array([[nd.x, nd.y, nd.z] for nd in mesh.nodes])
    part = BEMFSIParticipant(mesh, build_blade_aero_from_aerodyn(fx.AD_PRIMARY), fx.BEM_CONFIG)
    fa, _ = part._compute_forces(np.zeros((len(X), 3)))
    f = np.zeros(n)
    for d in range(3):
        f[d::6] = fa[:, d]
    u = np.zeros(n)
    u[free] = spsolve(K[np.ix_(free, free)].tocsc(), f[free])
    U = np.column_stack([u[0::6], u[1::6], u[2::6]])
    ref = np.stack([part._ref_coords[i].mean(axis=0) for i in part._strip_node_indices])
    um = np.stack([U[i].mean(axis=0) for i in part._strip_node_indices])
    r_def, _ = part._compute_deformed_geometry(U)
    # One point load per strip: every node binned to its nearest strip centroid along z,
    # moments about the strip's axis point (x = 0, y = centroid prebend).
    k_of = np.argmin(abs(X[:, 2][:, None] - ref[None, :, 2]), axis=1)
    F, M = np.zeros((len(ref), 3)), np.zeros((len(ref), 3))
    for k in range(len(ref)):
        m = k_of == k
        F[k] = fa[m].sum(0)
        M[k] = np.cross(X[m] - np.array([0.0, ref[k, 1], ref[k, 2]]), fa[m]).sum(0)
    return dict(ref=ref, um=um, F=F, M=M, dr=r_def - part._ref_r)


def write_deck(deck: Path, work: Path):
    """Copy the v4 deck into ``work`` in v5 format; return the key points (n, 3)."""
    prim = (deck / "IEA-15-240-RWT_BeamDyn.dat").read_text()
    prim = re.sub(r"-+ PITCH ACTUATOR PARAMETERS.*?PitchC[^\n]*\n", "", prim, flags=re.S)
    prim = re.sub(r"^False(\s+SumPrint)", r"True \1", prim, flags=re.M)
    prim = prim.replace('"TDyr"\n', '"TDyr"\n"TDzr"\n', 1)
    (work / "bd_primary.dat").write_text(prim)
    blade = (deck / "IEA-15-240-RWT_BeamDyn_blade.dat").read_text().splitlines()
    out = []
    for i, line in enumerate(blade):
        if "damp_type" in line:
            out.append(" 0   damp_flag")
            continue
        out.append(line)
        if i >= 2 and "mu1" in blade[i - 2]:  # the mu values row: v5 adds a modal block
            out += [" ---------------------- MODAL DAMPING ------------------", " 1   n_modes",
                    " 0.0   zeta"]
    (work / "IEA-15-240-RWT_BeamDyn_blade.dat").write_text("\n".join(out) + "\n")
    lines = prim.splitlines()
    i0 = next(i for i, ln in enumerate(lines) if "kp_xr" in ln) + 2
    return np.array([[float(v) for v in lines[i].split()[:3]] for i in range(i0, i0 + N_KP)])


def write_driver(work: Path, kp, s):
    arc = np.concatenate([[0.0], np.cumsum(np.linalg.norm(np.diff(kp, axis=0), axis=1))])
    eta = np.interp(s["ref"][:, 2], kp[:, 2], arc / arc[-1])
    Fb = np.column_stack([s["F"][:, 1], -s["F"][:, 0], s["F"][:, 2]])
    Mb = np.column_stack([s["M"][:, 1], -s["M"][:, 0], s["M"][:, 2]])
    rows = "\n".join(f"{e:.8f} " + " ".join(f"{v:.6e}" for v in (*fb, *mb))
                     for e, fb, mb in zip(eta, Fb, Mb))
    zeros = "\n".join(f"0   {name}" for name in
                      [f"DistrLoad({i})" for i in range(1, 7)] + [f"TipLoad({i})" for i in range(1, 7)])
    (work / "bd_driver.inp").write_text(f"""------- BEAMDYN Driver INPUT FILE --------------------------------
IEA-15 static, de-loading fixture aero load (dead), no gravity, no rotation
---------------------- SIMULATION CONTROL --------------------------------------
False   DynamicSolve
0   t_initial
1   t_final
0.01   dt
---------------------- GRAVITY PARAMETER --------------------------------------
0   Gx
0   Gy
0   Gz
---------------------- FRAME PARAMETER --------------------------------------
0   GlbPos(1)
0   GlbPos(2)
0   GlbPos(3)
---The following 3 by 3 matrix is the direction cosine matrix, GlbDCM(3,3),
---relates global frame to the initial blade root frame
1.0 0.0 0.0
0.0 1.0 0.0
0.0 0.0 1.0
True   GlbRotBladeT0
---------------------- ROOT VELOCITY PARAMETER ----------------------------------
0   RootVel(4)
0   RootVel(5)
0   RootVel(6)
---------------------- APPLIED FORCE ----------------------------------
{zeros}
{len(eta)}   NumPointLoads
eta Fx Fy Fz Mx My Mz
(-) (N) (N) (N) (N-m) (N-m) (N-m)
{rows}
---------------------- PRIMARY INPUT FILE --------------------------------------
"bd_primary.dat"    InputFile
---------------------- OUTPUTS --------------------------------------
0   WrVTK
1   VTK_fps
""")


def read_beamdyn(work: Path, kp):
    """Quadrature-point positions and displacements of the converged static step."""
    # v5 leaves Init_QP_E1 empty in the summary: rebuild the trapezoidal QPs (blade stations
    # refined by 2, as the deck sets) on the key-point curve by arc-length fraction.
    blade = (work / "IEA-15-240-RWT_BeamDyn_blade.dat").read_text().splitlines()
    k0 = next(i for i, ln in enumerate(blade) if "DISTRIBUTED PROPERTIES" in ln) + 1
    st = np.array([float(ln.split()[0]) for ln in blade[k0:] if len(ln.split()) == 1])
    eta = np.append(np.ravel(np.column_stack([st[:-1], 0.5 * (st[:-1] + st[1:])])), st[-1])
    arc = np.concatenate([[0.0], np.cumsum(np.linalg.norm(np.diff(kp, axis=0), axis=1))])
    qp = np.column_stack([np.interp(eta, arc / arc[-1], kp[:, a]) for a in range(3)])
    out = (work / "bd_driver.out").read_text().splitlines()
    hdr = next(ln for ln in out if ln.startswith("Time")).split("\t")
    col = {h.strip(): k for k, h in enumerate(hdr)}
    last = np.array(out[-1].split(), float)
    td = np.array([[last[col[f"N{j + 1:03d}_TD{a}r"]] for a in "xyz"] for j in range(len(qp))])
    return qp, td


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--deck", type=Path, default=os.environ.get("IEA15_OPENFAST_DIR"))
    ap.add_argument("--beamdyn", default=os.environ.get("BEAMDYN_DRIVER", "beamdyn_driver"))
    ap.add_argument("--work", type=Path, default=None, help="keep the BeamDyn run here")
    args = ap.parse_args()
    if args.deck is None:
        ap.error("--deck (or IEA15_OPENFAST_DIR) is required")
    work = args.work or Path(tempfile.mkdtemp(prefix="bd_radii_"))
    work.mkdir(parents=True, exist_ok=True)

    s = shell_state()
    kp = write_deck(args.deck, work)
    write_driver(work, kp, s)
    run = subprocess.run([args.beamdyn, "bd_driver.inp"], cwd=work, capture_output=True, text=True)
    if run.returncode != 0:
        sys.exit(f"beamdyn_driver failed in {work}:\n{run.stdout[-2000:]}")
    qp, td = read_beamdyn(work, kp)

    lin, nl = kinematics(qp, td)
    theta = np.hypot(*np.gradient(td[:, :2], qp[:, 2], axis=0).T).max()
    err = td[-1, 2] - lin[-1] - nl[-1]
    print(f"BeamDyn run in {work}")
    print(f"[BeamDyn own field, tip z={qp[-1, 2]:.3f}] u_flap {td[-1, 0]:+.4f}  u_axial {td[-1, 2]:+.4f}"
          f"  lin {lin[-1]:+.4f}  nl {nl[-1]:+.4f}  lin+nl {lin[-1] + nl[-1]:+.4f}")
    print(f"  formula error {err:+.4f} m = {err / abs(nl[-1]):+.2%} of nl;"
          f" truncation bound theta_max^2/4 = {theta**2 / 4:.2%} (theta_max {theta:.3f} rad)")
    k = int(np.argmax(s["ref"][:, 2]))
    z0 = s["ref"][k, 2]
    bd = np.array([np.interp(z0, qp[:, 2], td[:, a]) for a in range(3)])
    print(f"[shell vs BeamDyn at the tip-strip centroid z0={z0:.3f}]")
    print(f"  axial: shell first order {s['um'][k, 2]:+.4f}  shell r_def increment {s['dr'][k]:+.4f}"
          f"  BeamDyn {bd[2]:+.4f}")
    print(f"  flap: shell {s['um'][k, 1]:.4f}  BeamDyn {bd[0]:.4f}  ({s['um'][k, 1] / bd[0] - 1:+.2%})")


if __name__ == "__main__":
    main()
