"""T5c M2: shell twist vs the BeamDyn deck 6x6 for one TE-gap option (A literal, B te_fill, C gelcoat+adhesive).

One rated solve per case (mesh_transform hook of tools/l3_torque_sensitivity.solve_rated, K_ff factorised once by intercepting
spsolve as logs/check/unit_force_twist/shell2.py does) gives (i) the rated one-way field and (ii) a unit tip torque case on
the same K. Beam side: logs/check/unit_force_twist/beam.py (deck 6x6, section frames) integrated for the same loads:
unit tip couple, and the production nodal rated loads as distributed V(s), M(s) resultants of all nodes outboard of s about the
deck axis point. Both twists are the raw rotation about +z (mesh = BD span), so the ratio is sign-safe. Diagnostic only.

    scripts/aeroenv.sh python tools/diagnose_unit_force_twist_te.py N15 A --out logs/chord_scaled/t5c
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import scipy.sparse.linalg as spl

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / "tools"), str(ROOT / "logs/check/unit_force_twist")]

T_UNIT = 1.0e5  # N m tip couple; results per MN m
STATIONS_R = (60.0, 90.0, 110.0)


def main() -> int:
    import diagnose_gc_twist_chord as D
    import l3_torque_sensitivity as T
    import beam as B
    import twist_p1 as t

    ap = argparse.ArgumentParser()
    ap.add_argument("case")
    ap.add_argument("option", choices=["A", "B", "C"])
    ap.add_argument("--out", default="logs/chord_scaled/t5c")
    a = ap.parse_args()
    nd = float(a.case[1:])
    transform = {"A": None, "B": D.te_fill, "C": D.te_adhesive}[a.option]
    cap = {}

    def grab(A, b):
        cap["lu"] = spl.splu(A.tocsc())
        return cap["lu"].solve(b)

    spl.spsolve = grab
    part, disp = T.solve_rated(0.25, "gmsh", wall_flow=True, chord_divisions=nd, mesh_transform=transform)
    lu = cap["lu"]
    mesh = part._mesh
    n = part._n_nodes * 6
    root = {mesh.node_id_to_index[nid] for nid in mesh.get_node_set("RootNodes").node_ids}
    fixed = {6 * r + d for r in root for d in range(6)}
    free = np.array([i for i in range(n) if i not in fixed])
    xyz = np.asarray(part._ref_coords, float)
    rings = [np.asarray(r) for r in part._strip_ring_indices]
    zr = np.array([xyz[r, 2].mean() if len(r) else np.nan for r in rings])
    rr = np.asarray(part._ref_r, float)
    hub = float(np.nanmean(rr - zr))
    deck, arc, K, R = B.load_deck()
    kp, akp = deck.kp, deck.arc_kp
    z_arc = np.interp(arc, akp, kp[:, 2])
    r_arc = np.column_stack([np.interp(arc, akp, kp[:, k]) for k in range(3)])  # BD frame axis points
    r_mesh = np.column_stack([-r_arc[:, 1], r_arc[:, 0], r_arc[:, 2]])

    def beam_twist(V, M):  # V, M (n, 3) in mesh frame -> theta(s) about +z
        Vb, Mb = (t.MESH_TO_BD @ V.T).T, (t.MESH_TO_BD @ M.T).T
        loc = np.concatenate([np.einsum("nji,nj->ni", R, Vb), np.einsum("nji,nj->ni", R, Mb)], axis=1)[..., None]
        return t.cumtrapz(np.linalg.solve(K, loc)[:, 5, 0], arc)

    # (i) unit tip torque: pure couple on the last ring with >= 4 nodes
    k = len(rings) - 1
    while len(rings[k]) < 4:
        k -= 1
    P = xyz[rings[k]]
    d = P - P.mean(0)
    A = np.zeros((6, 3 * len(P)))
    for j, dj in enumerate(d):
        A[:3, 3 * j:3 * j + 3] = np.eye(3)
        A[3:, 3 * j:3 * j + 3] = np.array([[0, -dj[2], dj[1]], [dj[2], 0, -dj[0]], [-dj[1], dj[0], 0]])
    f = A.T @ np.linalg.solve(A @ A.T, np.array([0, 0, 0, 0, 0, T_UNIT]))
    load = np.zeros(n)
    for j, nd_ in enumerate(rings[k]):
        load[6 * nd_:6 * nd_ + 3] = f[3 * j:3 * j + 3]
    u = np.zeros(n)
    u[free] = lu.solve(load[free])
    disp_u = np.column_stack([u[0::6], u[1::6], u[2::6]])
    V0 = np.zeros((len(arc), 3))
    M0 = np.tile([0.0, 0.0, T_UNIT], (len(arc), 1))
    th_u = beam_twist(V0, M0)

    # (ii) rated one-way strip loads as distributed resultants outboard of s
    forces = np.asarray(part._compute_forces(np.zeros((disp.shape[0], 3)))[0], float)
    zn = xyz[:, 2]
    V = np.zeros((len(arc), 3))
    M = np.zeros((len(arc), 3))
    for i in range(len(arc)):
        sel = zn > z_arc[i]
        V[i] = forces[sel].sum(0)
        M[i] = np.cross(xyz[sel] - r_mesh[i], forces[sel]).sum(0)
    th_r = beam_twist(V, M)

    sc = np.rad2deg(1.0)
    rows = []
    for rs in STATIONS_R:
        ks = int(np.nanargmin(np.abs(rr - rs)))
        while len(rings[ks]) < 4:
            ks -= 1
        ring = rings[ks]
        zc = float(zr[ks])
        wu = part._section_rotation(disp_u, ring)
        wr = part._section_rotation(disp, ring)
        bu = float(np.interp(zc, z_arc, th_u))
        br = float(np.interp(zc, z_arc, th_r))
        rows.append({"r_target": rs, "r": float(rr[ks]), "z": zc, "n_ring": int(len(ring)),
                     "unit_shell_deg_per_MNm": wu * sc * 1e6 / T_UNIT, "unit_beam_deg_per_MNm": bu * sc * 1e6 / T_UNIT,
                     "rated_shell_deg": wr * sc, "rated_beam_deg": br * sc})
        r_ = rows[-1]
        r_["unit_ratio"] = r_["unit_shell_deg_per_MNm"] / r_["unit_beam_deg_per_MNm"]
        r_["rated_ratio"] = r_["rated_shell_deg"] / r_["rated_beam_deg"]
        print(a.case, a.option, {kk: (round(v, 4) if isinstance(v, float) else v) for kk, v in r_.items()}, flush=True)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    (out / f"m2-{a.case}-{a.option}.json").write_text(json.dumps({"case": a.case, "option": a.option, "hub_offset_r_minus_z": hub,
                                                                  "T_unit_Nm": T_UNIT, "tip_ring_k": int(k), "stations": rows}, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
