#!/usr/bin/env python
"""S-2 at an interior station: shell arc ring mean vs BeamDyn FE-node translation (l3-closure S12, T10c-station).

Reads the shell ring profiles written by ``tools/s2_ring_metric_blade.py`` (one JSON per element size) and the
BeamDyn node results written by ``tools/bd_static_reference.py post --nodes ...`` (scales 1e-4 and 2e-4), and prints
the signed error of the shell at each BeamDyn output node's initial z: the shell value is the arc-weighted ring mean,
linear in z between the two rings that bracket the node. The predeclared stations are the two nodes nearest the tip
whose tip distance is >= 2 tip chords; the verdict uses the inner one. Gates (predeclaration): BeamDyn scales agree,
the two stations within 0.5 pp, es 0.25 -> 0.125 moves < 0.5 pp, primary vs ``@equal`` load introduction < 0.3 pp.

Usage:  python tools/s2_station_compare.py --blade iea15 \\
            --shell logs/check/l3/t10c-station/iea15-es025.json logs/check/l3/t10c-station/iea15-es0125.json \\
            --out logs/check/l3/t10c-station/iea15-compare.json
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
BD_LOGS = {"iea15": ROOT / "logs/check/bd_static", "iea34": ROOT / "logs/check/bd_static/iea34"}
CASES = {"LC1": "LC1_gravity", "LC2": "LC2_tip_flap", "LC3": "LC3_tip_edge", "LC4": "LC4_distributed"}
SCALES = (1e-4, 2e-4)
MIN_CHORDS = 2.0


def component(case):
    """Dominant shell component: ux for the edgewise tip load, uy otherwise."""
    return 0 if case == "LC3" else 1


def bd_nodes(blade, case, order):
    """Linear BeamDyn node translations (shell axes) per output node, checked across the two load scales."""
    runs = []
    for scale in SCALES:
        path = BD_LOGS[blade] / f"{case}-p{order}-s{scale:g}-nodes" / "result.json"
        runs.append(json.loads(path.read_text())["nodes"])
    k = component(case)
    spread = max(abs(a[f"linear_shell_u{'xy'[k]}_m"] / b[f"linear_shell_u{'xy'[k]}_m"] - 1)
                 for a, b in zip(*runs, strict=True))
    return runs[0], spread


def shell_at(profile, z):
    """Arc ring mean translation linear in z between the two bracketing rings."""
    prof = np.array(sorted(profile, key=lambda r: r[0]))
    if not prof[0, 0] <= z <= prof[-1, 0]:
        raise ValueError(f"z {z:.3f} outside the stored ring profile [{prof[0, 0]:.3f}, {prof[-1, 0]:.3f}]")
    return np.array([np.interp(z, prof[:, 0], prof[:, 1 + k]) for k in range(3)])


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--blade", choices=tuple(BD_LOGS), required=True)
    ap.add_argument("--shell", nargs="+", required=True, help="s2_ring_metric_blade.py JSONs, coarse first")
    ap.add_argument("--order", type=int, default=60)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    shells = [json.loads(Path(p).read_text()) for p in args.shell]
    out = {"blade": args.blade, "shell": args.shell, "cases": {}}
    for case, key in CASES.items():
        nodes, spread = bd_nodes(args.blade, case, args.order)
        k = component(case)
        comp = f"linear_shell_u{'xy'[k]}_m"
        rec = {"bd_scale_spread": spread, "rows": []}
        print(f"\n{args.blade} {case} u{'xy'[k]}  (BeamDyn scales 1e-4 / 2e-4 agree within {spread:.2e})")
        z_tip = shells[0]["tip_z"]
        chord = shells[0]["tip_scale_m"]
        inner = [n for n in nodes if (z_tip - n["z"]) >= MIN_CHORDS * chord]
        stations = sorted(inner, key=lambda n: -n["z"])[:2]
        variants = [key] + ([f"{key}@equal"] if f"{key}@equal" in shells[0]["cases"] else [])
        for n in sorted(nodes, key=lambda n: n["z"]):
            depth = (z_tip - n["z"]) / chord
            row = {"node": n["node"], "z": n["z"], "chords": depth, "bd": n[comp], "err": {}}
            cells = []
            for var in variants:
                for sh in shells:
                    if depth <= 0 or n["z"] < sh["tip_z"] - sh["profile_depth_m"]:
                        continue
                    if var not in sh["cases"]:  # a case-subset run (--cases), e.g. es 0.0625 LC1 / LC4
                        continue
                    u = shell_at(sh["cases"][var]["profile"], n["z"])[k]
                    err = u / n[comp] - 1
                    row["err"][f"{var}|es{sh['element_size']}"] = err
                    cells.append(f"{'eq' if '@' in var else 'pr'} es{sh['element_size']} {err:+.3%}")
            tag = " <- station" if n in stations else ""
            print(f"  node {n['node']:2d} z {n['z']:8.3f} ({depth:5.2f} chords) BD {n[comp]: .6f}  " + "  ".join(cells) + tag)
            row["station"] = n in stations
            rec["rows"].append(row)
        out["cases"][case] = rec
    Path(args.out).write_text(json.dumps(out, indent=2) + "\n")


if __name__ == "__main__":
    main()
