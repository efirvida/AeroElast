#!/usr/bin/env python
"""W8 / L3 T11: 2D diagnostics of the TE-filler section: GJ_st2d vs in-plane refinement (h, filler dx), coincident-node census,
and where the torsion shear flow sits (energy share of the filler / step region). Run: aeroenv.sh env PYTHONPATH=$PWD python tools/l3_w8_diag2d.py ETA ARC H TMAX [DX [NTHICK [ESCALE]]]"""
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
import l3_solid_section_ccx as T  # noqa: E402


def census(m2, tol=1e-7):
    xy = m2["xy"]
    key = np.round(xy / tol).astype(np.int64)
    _, cnt = np.unique(key, axis=0, return_counts=True)
    ids = {}
    for q in m2["quads"]:
        for i in q[0]:
            ids[i] = 1
    return len(xy), int((cnt > 1).sum()), len(ids)


def main():
    eta, arc, h, tmax = float(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4])
    dx = float(sys.argv[5]) if len(sys.argv) > 5 else 0.02
    nth = int(sys.argv[6]) if len(sys.argv) > 6 else 1
    T_dx = T.te_filler_stations.__defaults__
    T.te_filler_stations.__defaults__ = (dx, T_dx[1])
    sec = T.blade_section(eta)
    esc = float(sys.argv[7]) if len(sys.argv) > 7 else 1.0
    sec.materials["adhesive"] = {**sec.materials["adhesive"], "E": sec.materials["adhesive"]["E"] * esc, "G": sec.materials["adhesive"]["G"] * esc}
    t0 = time.time()
    fill = {"mat": "adhesive", "arc": arc, "nthick": nth} if arc > 0 else None
    m2 = T.build_2d(sec, h, tmax, "bonded", fill)
    n, dup, used = census(m2)
    gj = T.torsion_2d(sec, m2)
    sf = T.deck_metrics(T.deck_station(eta))[1]["GJ"]
    tags = {}
    for q in m2["quads"]:
        tags[q[4]] = tags.get(q[4], 0) + 1
    print(f"eta {eta} arc {arc} h {h} tmax {tmax} dx {dx} nthick {nth}: nodes {n} coincident-unmerged {dup} unused {n-used} quads {len(m2['quads'])} {tags}")
    print(f"  GJ_st2d {gj:.5e} = {gj/sf:.4f} of deck GJ_sc ({time.time()-t0:.0f} s)")


if __name__ == "__main__":
    main()
