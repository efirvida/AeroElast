#!/usr/bin/env python
"""L3 / T9 (W2): estimate of the deck EA / EI_flap consequence of the SONATA cap-area loss: dEA = E1*dA, dEI = E1*dA*y^2 (y = cap centre
height above the mid-height of the two caps, blended airfoil, chord frame; ignores the shift of the neutral axis and the foam in the wedge).
args: eta dA_total(m2 CarbonUD missing in SONATA) t_cap   (uses l3_w2_cap_hand.py's geometry via import is avoided: re-derives y here)"""

import sys
import numpy as np
import yaml

eta, dA, t = float(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3])
Y = "tests/IEA15MW/reference/upstream/IEA-15-240-RWT.yaml"
D = "tests/IEA15MW/reference/IEA-15-240-RWT_BeamDyn_blade.dat"
y = yaml.safe_load(open(Y))
E1 = [m["E"][0] for m in y["materials"] if m["name"] == "CarbonUD"][0]
# deck: station eta line then 6x6 K, then 6x6 M
L = [ln.split() for ln in open(D).read().splitlines()]
st = []
i = 0
while i < len(L):
    if len(L[i]) == 1:
        try:
            e = float(L[i][0])
            K = np.array([[float(v) for v in L[i + 1 + j]] for j in range(6)])
            st.append((e, K))
            i += 13
            continue
        except (ValueError, IndexError):
            pass
    i += 1
es = np.array([s[0] for s in st])
EA = np.interp(eta, es, [s[1][2, 2] for s in st])
EI44 = np.interp(eta, es, [s[1][3, 3] for s in st])
EI55 = np.interp(eta, es, [s[1][4, 4] for s in st])
# cap y from the blended airfoil (same blend as l3_w2_cap_hand.py)
bl = y["components"]["blade"]
o = bl["outer_shape"]
pos = [(a["name"], a["spanwise_position"]) for a in o["airfoils"]]
afs = {a["name"]: a["coordinates"] for a in y["airfoils"]}


def rs(n):
    c = afs[n]
    xy = np.c_[c["x"], c["y"]]
    xy = np.vstack([xy, xy[:1]]) if np.hypot(*(xy[0] - xy[-1])) > 1e-9 else xy
    s = np.r_[0, np.cumsum(np.hypot(*np.diff(xy, axis=0).T))]
    s /= s[-1]
    g = np.linspace(0, 1, 4001)
    return np.c_[np.interp(g, s, xy[:, 0]), np.interp(g, s, xy[:, 1])]


for (n0, p0), (n1, p1) in zip(pos[:-1], pos[1:], strict=False):
    if p0 <= eta <= p1:
        w = (eta - p0) / (p1 - p0)
        xy = ((1 - w) * rs(n0) + w * rs(n1)) * float(
            np.interp(eta, o["chord"]["grid"], o["chord"]["values"])
        )
        break
anc = {a["name"]: a for a in bl["structure"]["anchors"]}
ndc = lambda name: (
    0.5
    * (
        np.interp(eta, anc[name]["start_nd_arc"]["grid"], anc[name]["start_nd_arc"]["values"])
        + np.interp(eta, anc[name]["end_nd_arc"]["grid"], anc[name]["end_nd_arc"]["values"])
    )
)
seg = np.hypot(*np.diff(xy, axis=0).T)
s = np.r_[0, np.cumsum(seg)]
P = s[-1]
ycap = [np.interp(ndc(n) * P, s, xy[:, 1]) for n in ("Spar_Cap_SS", "Spar_Cap_PS")]
y0 = 0.5 * (ycap[0] + ycap[1])
yc = 0.5 * abs(ycap[0] - ycap[1]) - 0.5 * t  # inward by t/2 (cap centre)
print(
    f"eta {eta}: deck EA {EA:.4e} EI44 {EI44:.4e} EI55 {EI55:.4e}; E1 {E1:.3e}; cap y_c = {yc:.3f} m (OML cap y {ycap[0]:.3f}/{ycap[1]:.3f})"
)
print(
    f"dA {dA:.5f} m2 -> dEA = {E1 * dA:.4e} ({100 * E1 * dA / EA:.2f} % of deck EA); dEI_flap = {E1 * dA * yc**2:.4e} ({100 * E1 * dA * yc**2 / EI55:.2f} % of deck EI55 = EI_min = flap; {100 * E1 * dA * yc**2 / EI44:.2f} % of EI44 = edge)"
)
