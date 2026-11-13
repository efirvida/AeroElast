"""T5d: why the rated shell/beam twist ratio is -0.99 at r 59 m and 0.44 at r 90 m while a unit tip couple gives +1.11.

One rated N15 solve (K_ff factorised once). Saves forces/coords/rings to an npz, then
 (1) beam side: per-station twist rate split by 6x6 compliance row 6 (S66*Mz, S65*My, S64*Mx, S61..63*V), reference-axis variants
     (deck axis, shell ring centroid) for the strip moment arm;
 (2) shell side: rated load = torsional part (pure z-couple per strip about the deck axis, zero net force) + force part (rest); two solves.
Diagnostic only.   scripts/aeroenv.sh timeout 1500 python tools/diagnose_rated_twist_sign_t5d.py N15 B
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
OUT = ROOT / "logs/chord_scaled/t5d"
STATIONS_R = (60.0, 90.0, 110.0)


def main() -> int:
    import beam as B
    import diagnose_gc_twist_chord as D
    import l3_torque_sensitivity as T
    import twist_p1 as t

    ap = argparse.ArgumentParser()
    ap.add_argument("case")
    ap.add_argument("option", choices=["A", "B"])
    a = ap.parse_args()
    nd = float(a.case[1:])
    transform = {"A": None, "B": D.te_fill}[a.option]
    cap = {}

    def grab(A, b):
        cap["lu"] = spl.splu(A.tocsc())
        return cap["lu"].solve(b)

    spl.spsolve = grab
    part, disp = T.solve_rated(0.25, "gmsh", wall_flow=True, chord_divisions=nd, mesh_transform=transform)
    lu, mesh = cap["lu"], part._mesh
    n = part._n_nodes * 6
    root = {mesh.node_id_to_index[nid] for nid in mesh.get_node_set("RootNodes").node_ids}
    fixed = {6 * r + d for r in root for d in range(6)}
    free = np.array([i for i in range(n) if i not in fixed])
    xyz = np.asarray(part._ref_coords, float)
    rings = [np.asarray(r) for r in part._strip_ring_indices]
    zr = np.array([xyz[r, 2].mean() if len(r) else np.nan for r in rings])
    rr = np.asarray(part._ref_r, float)
    forces = np.asarray(part._compute_forces(np.zeros((disp.shape[0], 3)))[0], float)
    np.savez_compressed(OUT / f"{a.case}-{a.option}-rated.npz", xyz=xyz, forces=forces, disp=disp, zr=zr, rr=rr)

    deck, arc, K, R = B.load_deck()
    kp, akp = deck.kp, deck.arc_kp
    z_arc = np.interp(arc, akp, kp[:, 2])
    r_arc = np.column_stack([np.interp(arc, akp, kp[:, k]) for k in range(3)])
    r_mesh = np.column_stack([-r_arc[:, 1], r_arc[:, 0], r_arc[:, 2]])
    S = np.linalg.inv(K)  # (n,6,6)

    def resultants(F, axis):
        V = np.zeros((len(arc), 3)); M = np.zeros((len(arc), 3))
        for i in range(len(arc)):
            sel = xyz[:, 2] > z_arc[i]
            V[i] = F[sel].sum(0)
            M[i] = np.cross(xyz[sel] - axis[i], F[sel]).sum(0)
        return V, M

    def loc_of(V, M):
        Vb, Mb = (t.MESH_TO_BD @ V.T).T, (t.MESH_TO_BD @ M.T).T
        return np.concatenate([np.einsum("nji,nj->ni", R, Vb), np.einsum("nji,nj->ni", R, Mb)], axis=1)

    def twist_terms(loc):  # per-component theta(s) contributions, shape (6, n)
        return np.array([t.cumtrapz(S[:, 5, k] * loc[:, k], arc) for k in range(6)])

    sel_rings = []
    for rs in STATIONS_R:
        ks = int(np.nanargmin(np.abs(rr - rs)))
        while len(rings[ks]) < 4:
            ks -= 1
        sel_rings.append(ks)
    deg = np.rad2deg(1.0)
    res = {"case": a.case, "option": a.option, "stations": []}

    # shell ring-centroid axis (mid-surface centroid per strip) interpolated on arc z
    cen = np.array([xyz[r].mean(0) if len(r) else [np.nan] * 3 for r in rings])
    ok = ~np.isnan(zr)
    cen_axis = np.column_stack([np.interp(z_arc, zr[ok], cen[ok, k]) for k in range(3)])

    V, M = resultants(forces, r_mesh)
    terms = twist_terms(loc_of(V, M))
    Vc, Mc = resultants(forces, cen_axis)
    terms_c = twist_terms(loc_of(Vc, Mc))
    # rotating local frame at the station: also print local loc comps and the S row
    loc = loc_of(V, M)

    # shell split: pure z couple per strip about deck axis (zero net force) + rest
    axis_at = lambda z: np.array([np.interp(z, z_arc, r_mesh[:, k]) for k in range(3)])
    FT = np.zeros_like(forces)
    strip_Mz = []
    for k, r in enumerate(rings):
        if len(r) < 4:
            strip_Mz.append(0.0); continue
        p = xyz[r] - axis_at(zr[k])
        Mz = np.sum(p[:, 0] * forces[r, 1] - p[:, 1] * forces[r, 0])
        q = p[:, :2]
        den = np.sum(q ** 2)
        FT[r, 0] = -Mz * q[:, 1] / den
        FT[r, 1] = Mz * q[:, 0] / den
        strip_Mz.append(float(Mz))
    FF = forces - FT
    # beam-equivalent resultant of each strip (force F_k and moments Mx, My about the deck axis, Mz = 0) by the min-norm nodal field
    # (equal weights, same realisation as logs/check/unit_force_twist/shell2.py) and the beam-invisible remainder D (zero force and
    # zero moment per strip about the deck axis)
    FR = np.zeros_like(forces)
    in_ring = np.zeros(len(xyz), bool)
    for k, r in enumerate(rings):
        if len(r) < 4:
            continue
        in_ring[r] = True
        p = xyz[r] - axis_at(zr[k])
        A6 = np.zeros((6, 3 * len(r)))
        for j, dj in enumerate(p):
            A6[:3, 3 * j:3 * j + 3] = np.eye(3)
            A6[3:, 3 * j:3 * j + 3] = np.array([[0, -dj[2], dj[1]], [dj[2], 0, -dj[0]], [-dj[1], dj[0], 0]])
        tgt = A6 @ FF[r].reshape(-1)
        tgt[5] = 0.0  # Mz already carried by FT
        FR[r] = np.linalg.lstsq(A6, tgt, rcond=None)[0].reshape(-1, 3)
    FD = FF - FR
    res["nodes_outside_rings"] = int((~in_ring).sum())
    res["force_outside_rings_N"] = forces[~in_ring].sum(0).tolist()
    res["load_norm_N"] = {"full": float(np.linalg.norm(forces)), "T": float(np.linalg.norm(FT)), "R": float(np.linalg.norm(FR)), "D": float(np.linalg.norm(FD))}
    th = {}
    for name, F in (("full", forces), ("torque_only", FT), ("forces_only", FF), ("resultant_R", FR), ("residual_D", FD)):
        load = np.zeros(n)
        for d in range(3):
            load[d::6] = F[:, d]
        u = np.zeros(n)
        u[free] = lu.solve(load[free])
        dsp = np.column_stack([u[0::6], u[1::6], u[2::6]])
        th[name] = [float(part._section_rotation(dsp, rings[ks])) * deg for ks in sel_rings]
        th[name + "_slopefit"] = [float(T._ring_rotation_with_slope(part, dsp, rings[ks])) * deg for ks in sel_rings]
        if name == "full":
            ztilt = []
            for ks in sel_rings:
                P = xyz[rings[ks]]
                ztilt.append({"z_spread_m": float(np.ptp(P[:, 2])), "chord_x_spread_m": float(np.ptp(P[:, 0])),
                              "flap_slope_dudz_uy": float(np.polyfit(P[:, 2], dsp[rings[ks], 1], 1)[0])})
            res["ring_tilt_and_slope_full"] = ztilt
    # beam: torque-only (V=0, M=(0,0,Mz outboard sum)) vs coupling
    Mz_cum = np.zeros(len(arc))
    zs = np.array([zr[k] if len(rings[k]) >= 4 else np.nan for k in range(len(rings))])
    for i in range(len(arc)):
        Mz_cum[i] = sum(strip_Mz[k] for k in range(len(rings)) if len(rings[k]) >= 4 and zs[k] > z_arc[i])
    Mt = np.zeros((len(arc), 3)); Mt[:, 2] = Mz_cum
    th_bt = t.cumtrapz(np.einsum("nk,nk->n", S[:, 5, :], loc_of(np.zeros_like(Mt), Mt)), arc) * deg

    for ks in sel_rings:
        zc = float(zr[ks])
        i = int(np.argmin(abs(z_arc - zc)))
        row = {"r": float(rr[ks]), "z": zc, "n_ring": int(len(rings[ks])),
               "deck_axis_mesh_xy": r_mesh[i, :2].tolist(), "ring_centroid_xy": cen[ks, :2].tolist(),
               "V_out_N": V[i].tolist(), "M_out_Nm_deckaxis": M[i].tolist(), "M_out_Nm_ringcentroid": Mc[i].tolist(),
               "loc_V_M_deckaxis": loc[i].tolist(),
               "S6row": S[i, 5, :].tolist(),
               "beam_twist_terms_deg_deckaxis(V1,V2,V3,Mx,My,Mz)": (terms[:, i] * deg).tolist(),
               "beam_twist_total_deg_deckaxis": float(terms[:, i].sum() * deg),
               "beam_twist_total_deg_ringcentroid": float(terms_c[:, i].sum() * deg),
               "beam_twist_torque_only_Mz_deg": float(th_bt[i])}
        res["stations"].append(row)
    for name in th:
        for j, ks in enumerate(sel_rings):
            res["stations"][j]["shell_" + name + "_deg"] = th[name][j]
    (OUT / f"t5d-{a.case}-{a.option}.json").write_text(json.dumps(res, indent=1))
    for r in res["stations"]:
        print(json.dumps({k: (np.round(v, 5).tolist() if isinstance(v, (list, float)) else v) for k, v in r.items()}), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
