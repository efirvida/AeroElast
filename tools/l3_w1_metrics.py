#!/usr/bin/env python
"""W1 / L3 T3+T4: one frame-invariant metric (K4 = shear-free (F,T,Mx,My) stiffness, bending = eigenvalues after axial condensation,
GJ = K4[T,T]) for the deck, SONATA 2026, SONATA 2021 and the CalculiX solid. Run: python tools/l3_w1_metrics.py (needs numpy only)."""
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
import l3_solid_section_ccx as T  # noqa: E402

S = Path("/scratch/leahk/eduardo.donestevez/tmp")
SON26 = {0.5: S / "t3/out_r1200/sonata_t3.json", 0.8: S / "t3/out_r1200/sonata_t3.json"}
SON21 = {0.5: S / "t3b/out/lo_r800/sonata_t3b.json", 0.8: S / "t3b/out/hi_r800/sonata_t3b.json"}


def sonata_inv(path, eta):
    r = json.load(open(path))
    i = [k for k, x in enumerate(r["stations"]) if abs(x - eta) < 1e-6][0]
    C = np.linalg.inv(np.array(r["TS"][i]))  # SONATA order: axial 0, shears 1-2, torsion 3, bending 4-5
    K4 = np.linalg.inv(C[np.ix_([0, 3, 4, 5], [0, 3, 4, 5])])
    return {**T.invariants(K4), "GJ_pure": gj_pure(K4)}, np.array(r["TS"][i])


def variants(eta):
    """SONATA 2026 factor variants (tools/l3_w1_sonata_variants.py output)."""
    f = S / f"w1/var26_{eta}.json"
    out = {}
    if f.exists():
        for k, v in json.load(open(f))["variants"].items():
            C = np.linalg.inv(np.array(v["TS"]))
            K4 = np.linalg.inv(C[np.ix_([0, 3, 4, 5], [0, 3, 4, 5])])
            out["S26:" + k] = {**T.invariants(K4), "GJ_pure": gj_pure(K4)}
    return out


def gj_pure(K4):
    """torsion stiffness under a pure torque (axial force and both bending moments zero): 1/C[T,T] of the shear-free (F,T,Mx,My) compliance.
    The deck / SONATA K66 is about the reference axis, so it carries any torsion-bending coupling; this one does not."""
    return 1.0 / np.linalg.inv(K4)[1, 1]


def main():
    rows = []
    for eta in (0.5, 0.8):
        raw, sf = T.deck_metrics(T.deck_station(eta))
        D = T.deck_station(eta)
        print(f"eta {eta}: deck K66 diag GJ {D[5,5]:.4e}  GJ_sc(shear-free) {sf['GJ']:.4e}  ratio K66diag/GJ_sc {D[5,5]/sf['GJ']:.3f}")
        K4d = np.linalg.inv(np.linalg.inv(D)[np.ix_([2, 5, 3, 4], [2, 5, 3, 4])])
        cand = {"deck": {**sf, "GJ_pure": gj_pure(K4d)}}
        for nm, d in (("SONATA2026", SON26), ("SONATA2021", SON21)):
            if d[eta].exists():
                try:
                    cand[nm] = sonata_inv(d[eta], eta)[0]
                except IndexError:
                    pass
        cand.update(variants(eta))
        for nm, tag in (("solid_bonded", "bon_B"), ("solid_through", "thr_B")):
            p = S / f"t4w/e{int(eta*100)}_{tag}/result.json"
            if p.exists():
                m = json.load(open(p))["mid_third"]
                cand[nm] = {**{k: m[k] for k in ("EA", "GJ", "EI_min", "EI_max")}, "GJ_pure": gj_pure(np.array(m["K4"]))}
        for nm, v in cand.items():
            rows.append((eta, nm, *(v[k] / sf[k] for k in ("GJ", "EI_min", "EI_max", "EA")), v["GJ"], v["GJ_pure"] / cand["deck"]["GJ_pure"]))
    print("eta  model          GJ/deck  flap/deck edge/deck EA/deck   GJ [N m2]   GJ_pure/deck_pure")
    for r in rows:
        print(f"{r[0]:.2f} {r[1]:14s} {r[2]:.3f}    {r[3]:.3f}     {r[4]:.3f}     {r[5]:.3f}    {r[6]:.4e}   {r[7]:.3f}")


if __name__ == "__main__":
    main()
