"""W4/T5: tables from logs/chord_scaled/w4/<case>.json: per z x N FE/Bredt (compliance, +-3 m line fit) and distortion; decay length; bands."""

import json
import sys

import numpy as np

cases = sys.argv[1:]  # e.g. n15 n20 n30 [n15z105 ...] (file stems)
D = {c: json.load(open(f"logs/chord_scaled/w4/{c}.json")) for c in cases}
key = "fe_over_bredt_th_W3.0"
zs = sorted({r["z"] for d in D.values() for r in d["rows"] if 88 <= r["z"] <= 116 and key in r})
print(
    "FE local GJ / Bredt (line fit +-3 m) | distortion rel. to rms(r)*theta (rigid fit with z-gradient):"
)
print("  z  " + "  ".join(f"{c:>16s}" for c in cases))
for z in zs:
    cells = []
    for c in cases:
        r = next((r for r in D[c]["rows"] if r["z"] == z), None)
        cells.append(
            f"{r[key]:6.3f} | {r['dist_rel_theta']:6.3f}" if r and key in r else " " * 7 + "-"
        )
    print(f"{z:5.0f} " + "  ".join(f"{x:>16s}" for x in cells))
print("\nband compliance ratio (rigid / affine twist):")
bk = list(next(iter(D.values()))["bands"])
print("band    " + "  ".join(f"{c:>14s}" for c in cases))
for b in bk:
    print(
        f"{b:8s}"
        + "  ".join(
            (
                f"{D[c]['bands'][b]['fe_over_bredt_th']:6.3f}/{D[c]['bands'][b]['fe_over_bredt_th_aff']:6.3f}"
                if b in D[c]["bands"]
                else "       -      "
            ).rjust(14)
            for c in cases
        )
    )
print(
    "\ndecay length from the tip ring (tip z - highest-inboard z where FE/Bredt >= 0.9*mid-span mean of z 60-85, and stays so for the next 2 stations inward... first scanning inward from z 112):"
)
for c in cases:
    rows = [r for r in D[c]["rows"] if key in r]
    mid = np.mean([r[key] for r in rows if 60 <= r["z"] <= 85])
    zt = D[c]["zmax"]
    zz = [r["z"] for r in rows if r["z"] <= 112]
    vv = [r[key] for r in rows if r["z"] <= 112]
    out = None
    for i in range(len(zz) - 1, 1, -1):  # scan inward from the outermost station
        if vv[i] >= 0.9 * mid and vv[i - 1] >= 0.9 * mid and vv[i - 2] >= 0.9 * mid:
            out = zz[i]
            break
    print(
        f"{c:10s} mid={mid:.3f} threshold={0.9 * mid:.3f} outermost station of 3 consecutive >= thr: z={out} -> decay length {None if out is None else round(zt - out, 1)} m (tip z {zt:.1f}); zload={D[c]['zload']}"
    )
