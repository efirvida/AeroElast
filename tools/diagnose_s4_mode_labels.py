"""Diagnose the S-4 mode-label classifier at fine meshes.

``classify_modes`` labels a mode from the translation of the single
highest-z tip node with a 2:1 dominance threshold.  At the 0.125 m mesh the
first-edgewise pick jumps out of the convergence trend of the S-x artifact
(0.6970 -> 0.7019 Hz), so this dumps, per mesh size and parked mode:

  * the tip-node components and the |u_y|/|u_x| ratio,
  * the span-integrated participation ratio (all blade nodes) and the ratio
    restricted to the top 10 % of the span,
  * the legacy label and a robust integrated label.

Usage:
  python tools/diagnose_s4_mode_labels.py --sizes 0.5 0.25 0.125
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "tools"))

DPN = 6


def _participation(mesh, mode, free):
    """Span-integrated absolute participation per translation component."""
    node = free // DPN
    comp = free % DPN
    z = np.asarray([n.coords for n in mesh.nodes])[node, 2]
    ztip = float(z.max())
    a = np.abs(mode)
    out = {}
    for c, name in ((0, "x"), (1, "y"), (2, "z")):
        m = comp == c
        out[name] = float(a[m].sum())
        out["top_" + name] = float(a[m & (z > 0.9 * ztip)].sum())
    return out


def _label(sx, sy, sz, factor=2.0):
    if sy >= factor * sx and sy >= factor * sz:
        return "flap"
    if sx >= factor * sy and sx >= factor * sz:
        return "edge"
    if sz >= factor * sx and sz >= factor * sy:
        return "axial"
    return "mixed"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sizes", type=float, nargs="+", default=[0.5, 0.25, 0.125])
    ap.add_argument("--n-modes", type=int, default=8)
    ap.add_argument(
        "--csv",
        type=Path,
        default=REPO / "docs" / "validation_data" / "generated" / "s4_mode_labels.csv",
    )
    args = ap.parse_args()

    from run_s4_rotating_modal import assemble_rotating, classify_modes
    from run_s5_oneway import build_assembly

    rows = []
    for es in args.sizes:
        asm, mesh = build_assembly(es)
        freqs, modes, free = assemble_rotating(asm, mesh, 0.0, n_modes=args.n_modes)
        labels = classify_modes(mesh, modes, free)

        nid2idx = mesh.node_id_to_index
        xyz = np.asarray([n.coords for n in mesh.nodes])
        tip = int(np.argmax(xyz[:, 2]))
        tdof0 = nid2idx[tip] * DPN

        print(f"\n=== element size {es} m  ({len(mesh.nodes)} nodes) ===")
        for k in np.argsort(freqs):
            u = np.zeros(3)
            for d in range(3):
                dof = tdof0 + d
                pos = np.searchsorted(free, dof)
                if pos < len(free) and free[pos] == dof:
                    u[d] = modes[k, pos]
            p = _participation(mesh, modes[k], free)
            r_tip = abs(u[1]) / max(abs(u[0]), 1e-30)
            r_span = p["y"] / max(p["x"], 1e-30)
            r_top = p["top_y"] / max(p["top_x"], 1e-30)
            robust = _label(p["x"], p["y"], p["z"])
            rows.append(
                (
                    es, len(mesh.nodes), int(k), float(freqs[k]),
                    u[0], u[1], u[2], r_tip, labels[k], r_span, r_top, robust,
                )
            )
            print(
                f"  {k:2d} f={freqs[k]:8.4f}  tip=({u[0]:+.4f},{u[1]:+.4f},{u[2]:+.4f})"
                f"  r_tip={r_tip:7.3f}  legacy={labels[k]:6s}"
                f" | r_span={r_span:7.3f}  r_top10={r_top:7.3f}  robust={robust}"
            )

    args.csv.parent.mkdir(parents=True, exist_ok=True)
    with args.csv.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(
            [
                "element_size_m", "n_nodes", "mode", "freq_Hz",
                "ux_tip", "uy_tip", "uz_tip", "ratio_uy_ux_tip",
                "label_legacy", "ratio_span", "ratio_top10pct", "label_robust",
            ]
        )
        w.writerows(rows)
    print(f"\nwrote {args.csv}")


if __name__ == "__main__":
    main()
