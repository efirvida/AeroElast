#!/usr/bin/env python
"""L3 / T9 (W2): independent hand computation of cap / reinforcement layer area from the windIO yaml only.

area_hand = t(eta) * (end_nd_arc - start_nd_arc)(eta) * P_oml(eta)         (flat-strip estimate, arc on the OML)
area_exact = t * L - t^2/2 * dtheta   (OML arc L, tangent turning dtheta; inner arc shorter by t*dtheta, layer grown inward)
Perimeter / curvature from the blended airfoil (linear blend of arc-length-resampled coordinates) scaled by chord. Not the
solid tool's code (tools/l3_solid_section_ccx.py). Usage: python tools/l3_w2_cap_hand.py [yaml] eta ...
"""
import sys
import numpy as np
import yaml

Y = sys.argv[1] if len(sys.argv) > 1 and sys.argv[1].endswith(".yaml") else "tests/IEA15MW/reference/upstream/IEA-15-240-RWT.yaml"
ETAS = [float(a) for a in sys.argv[1:] if not a.endswith(".yaml")] or [0.5, 0.8]
y = yaml.safe_load(open(Y))
bl = y["components"]["blade"]
v2 = "outer_shape" in bl
os_ = bl["outer_shape"] if v2 else bl["outer_shape_bem"]
st = bl["structure"] if v2 else bl["internal_structure_2d_fem"]
afs = {a["name"]: a["coordinates"] for a in y["airfoils"]}
pos = ([(a["name"], a["spanwise_position"]) for a in os_["airfoils"]] if v2
       else list(zip(os_["airfoil_position"]["labels"], os_["airfoil_position"]["grid"])))

ANCH = {a["name"]: a for a in st.get("anchors", [])} if v2 else {}
for l in st["layers"]:
    for k in ("start_nd_arc", "end_nd_arc"):
        if v2 and k in l and "grid" in l[k]:
            ANCH.setdefault(l["name"], {})[k] = l[k]

def f(d, eta):
    return float(np.interp(eta, d["grid"], d["values"]))

def resample(name, n=4001):
    c = afs[name]; xy = np.c_[c["x"], c["y"]]
    if np.hypot(*(xy[0] - xy[-1])) > 1e-9:
        xy = np.vstack([xy, xy[:1]])  # close open TE (arc length then includes the TE edge; only used for blending)
    s = np.r_[0, np.cumsum(np.hypot(*np.diff(xy, axis=0).T))]; s /= s[-1]
    g = np.linspace(0, 1, n)
    return np.c_[np.interp(g, s, xy[:, 0]), np.interp(g, s, xy[:, 1])]

def oml(eta):
    for (n0, p0), (n1, p1) in zip(pos[:-1], pos[1:]):
        if p0 <= eta <= p1:
            w = 0.0 if p1 == p0 else (eta - p0) / (p1 - p0)
            xy = (1 - w) * resample(n0) + w * resample(n1)
            return xy * f(os_["chord"], eta)

def layer_geom(xy, a0, a1):
    """arc length and turning angle of the OML between nd-arc a0..a1 (nd arc measured on the open airfoil, TE->TE)."""
    seg = np.hypot(*np.diff(xy, axis=0).T); s = np.r_[0, np.cumsum(seg)]; P = s[-1]
    g = np.linspace(a0 * P, a1 * P, 800); X = np.c_[np.interp(g, s, xy[:, 0]), np.interp(g, s, xy[:, 1])]
    ang = np.unwrap(np.arctan2(*np.diff(X, axis=0).T[::-1]))
    return (a1 - a0) * P, abs(ang[-5:].mean() - ang[:5].mean()), P

for eta in ETAS:
    xy = oml(eta); print(f"\n== {Y.split('/')[-1]} eta {eta}  chord {f(os_['chord'], eta):.4f} m")
    tot = {}
    for l in st["layers"]:
        if "web" in l or l["material"] in ("Gelcoat",) or "filler" in l["name"] or l["name"].startswith("Shell"):
            continue
        t = f(l["thickness"], eta)
        if t <= 0:
            continue
        def nd(key):
            d = l.get(key)
            if d is None:
                return None
            if "anchor" in d:  # windIO v2: resolve the named anchor (top-level anchors, or the layer's own)
                a = ANCH[d["anchor"]["name"]][d["anchor"]["handle"]]
                return f(a, eta)
            if "grid" in d:
                return f(d, eta)
            return None
        a0, a1 = nd("start_nd_arc"), nd("end_nd_arc")
        if a0 is None: a0 = 0.0  # v1 TE_SS reinforcement starts at TE; end given only by yaml width (see below)
        if a1 is None: a1 = 1.0
        if not v2 and l["name"].startswith("TE_reinforcement") and "width" in l:
            wv = f(l["width"], eta)  # v1: arc extent = yaml width (m) from the TE, measured on the OML
            P0 = layer_geom(xy, 0, 1)[2]
            if "SS" in l["name"]: a1 = wv / P0
            else: a0 = 1 - wv / P0
        L, dth, P = layer_geom(xy, a0, a1)
        flat = t * L; exact = t * L - 0.5 * t * t * dth
        w = f(l["width"], eta) if "width" in l and isinstance(l["width"], dict) and "grid" in l["width"] else float("nan")
        print(f"{l['name']:22s} {l['material']:9s} t={t*1e3:7.3f} mm nd=[{a0:.4f},{a1:.4f}] L_oml={L:.4f} (yaml width {w:.3f}) dtheta={dth:.3f}"
              f" t*L={flat:.5f}  curv-corr={-0.5*t*t*dth:+.5f} -> {exact:.5f} m2")
        tot.setdefault(l["material"], [0, 0]); tot[l["material"]][0] += flat; tot[l["material"]][1] += exact
    print("per material (t*L | with curvature):", {k: (round(v[0], 4), round(v[1], 4)) for k, v in tot.items()})
