"""T5b: per-station G-C elastic twist on the chord-scaled gmsh mesh (odd/tasks/gmsh-chord-scaled-mesh.md).

For one case (N<k> or es<h>) solves the rated one-way field (production wall flow) and writes, per BEM strip: radius,
the production ring twist (``BEMFSIParticipant._section_rotation``), the same fit with a per-ring bending-slope column
(``l3_torque_sensitivity._ring_rotation_with_slope``), the ring node count, its span spread (non-planarity) and the
strip-mean flap/edge displacement. Diagnostic only.

    scripts/aeroenv.sh python tools/diagnose_gc_twist_chord.py N20 --out logs/chord_scaled/t5b
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))


def main() -> int:
    import l3_torque_sensitivity as l3

    ap = argparse.ArgumentParser()
    ap.add_argument("case")
    ap.add_argument("--out", default="logs/chord_scaled/t5b")
    a = ap.parse_args()
    if a.case.startswith("N"):
        es, nd = 0.25, float(a.case[1:])
    else:
        es, nd = float(a.case[2:]), None
    p, disp = l3.solve_rated(es, "gmsh", wall_flow=True, chord_divisions=nd)
    span = p._span_dir
    sense = p._twist_mesh_to_bem
    # the applied nodal load, recomputed exactly as solve_rated builds it (undeformed geometry, no nodal moments)
    forces = np.asarray(p._compute_forces(np.zeros((disp.shape[0], 3)))[0], float)
    rows = []
    for k, ring in enumerate(p._strip_ring_indices):
        idx = p._strip_node_indices[k]
        if len(idx) == 0:
            continue
        w0 = p._section_rotation(disp, ring)
        w1 = l3._ring_rotation_with_slope(p, disp, ring)
        zs = p._ref_coords[ring] @ span
        um = disp[idx].mean(axis=0)
        X = p._ref_coords[idx]
        F = forces[idx]
        cen = p._ref_coords[ring].mean(axis=0)
        mt_axis = float(np.cross(X - np.outer(X @ span, span), F).sum(axis=0) @ span)
        mt_cen = float(np.cross(X - cen, F).sum(axis=0) @ span)
        # extraction-independent rotation: the in-plane LE-TE chord line (farthest ring node pair) before/after
        R = p._ref_coords[ring]
        Rp = R - np.outer(R @ span, span)
        dd = np.linalg.norm(Rp[:, None, :] - Rp[None, :, :], axis=2)
        i, j = np.unravel_index(int(np.argmax(dd)), dd.shape)
        c0 = Rp[j] - Rp[i]
        Up = disp[ring] - np.outer(disp[ring] @ span, span)
        c1 = c0 + Up[j] - Up[i]
        chord_rot = float(np.arctan2(np.cross(c0, c1) @ span, c0 @ c1))
        # in-plane affine fit u = A x + t: symmetric part (section distortion) next to the antisymmetric part (rotation)
        e1 = c0 / np.linalg.norm(c0)
        e2 = np.cross(span, e1)
        P2 = np.column_stack([(Rp - Rp.mean(0)) @ e1, (Rp - Rp.mean(0)) @ e2])
        U2 = np.column_stack([Up @ e1, Up @ e2])
        G = np.column_stack([P2, np.ones(len(P2))])
        A = np.linalg.lstsq(G, U2, rcond=None)[0][:2].T
        resid = U2 - G @ np.linalg.lstsq(G, U2, rcond=None)[0]
        rows.append({
            "k": k, "r": float(p._ref_r[k]), "n_ring": int(len(ring)), "ring_z_ptp": float(np.ptp(zs)),
            "twist_prod_deg": None if w0 is None else float(np.rad2deg(sense * w0)),
            "twist_slope_deg": None if w1 is None else float(np.rad2deg(sense * w1)),
            "u_mean": [float(v) for v in um],
            "f_sum": [float(v) for v in F.sum(axis=0)],
            "mt_pitch_axis": float(sense * mt_axis), "mt_ring_centroid": float(sense * mt_cen),
            "ring_centroid": [float(v) for v in cen],
            "chord_len": float(np.linalg.norm(c0)), "chord_line_rot_deg": float(np.rad2deg(sense * chord_rot)),
            "affine_rot_deg": float(np.rad2deg(sense * 0.5 * (A[1, 0] - A[0, 1]))),
            "affine_shear": float(0.5 * (A[1, 0] + A[0, 1])), "affine_strain": [float(A[0, 0]), float(A[1, 1])],
            "fit_resid_rms": float(np.sqrt((resid ** 2).sum(1).mean())),
            "resid_max": float(np.sqrt((resid ** 2).sum(1)).max()),
            "resid_max_xyz": [float(v) for v in R[int(np.argmax((resid ** 2).sum(1)))]],
            "resid_max_u": [float(v) for v in disp[ring][int(np.argmax((resid ** 2).sum(1)))]],
            "resid_max_xi": float((P2[int(np.argmax((resid ** 2).sum(1))), 0] - P2[:, 0].min()) / np.ptp(P2[:, 0])),
            "resid_max_side": float(np.sign(P2[int(np.argmax((resid ** 2).sum(1))), 1])),
            "ring_inplane_force_max_over_mean": float(np.linalg.norm((forces[ring] - np.outer(forces[ring] @ span, span)), axis=1).max()
                                                      / max(np.linalg.norm((forces[ring] - np.outer(forces[ring] @ span, span)), axis=1).mean(), 1e-30)),
            "ring_inplane_force_sum_abs": float(np.linalg.norm((forces[ring] - np.outer(forces[ring] @ span, span)), axis=1).sum()),
            "ring_inplane_force_net": float(np.linalg.norm((forces[ring] - np.outer(forces[ring] @ span, span)).sum(0))),
        })
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    (out / f"{a.case}.json").write_text(json.dumps({"case": a.case, "n_nodes": int(disp.shape[0]), "strips": rows}, indent=1))
    for r in rows[-8:]:
        print(r["k"], round(r["r"], 2), r["n_ring"], round(r["ring_z_ptp"], 4), r["twist_prod_deg"], r["twist_slope_deg"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
