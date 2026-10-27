"""Prismatic copy of a windIO v1 blade: the section at eta held constant over the whole span (core-torsion C4).

Held at its value interpolated at eta: chord, pitch_axis, twist (set 0), the airfoil (a blend of the two airfoils
bracketing eta, y-coordinates and relative thickness weighted by the position fraction; airfoil_position becomes
[0, 1] -> the blend) and every layer / web grid (start/end/midpoint_nd_arc, thickness, width, offset_y_pa,
fiber_orientation) as the constant [0, 1] grid. A layer whose thickness grid does not contain eta is removed (with
a web whose layers are all removed); a partial-span layer present at eta becomes full span. reference_axis x = y = 0
(no prebend or sweep), z unchanged. Scratch input for odd/tasks/core-torsion.md (C4), never a test input.

usage: python tools/core_torsion_prismatic_yaml.py <in.yaml> <out.yaml> [--eta 0.8]
"""

import argparse

import numpy as np
import yaml

BLEND = "prismatic_blend"


def at(d, eta):
    return float(np.interp(eta, d["grid"], d["values"]))


def const(v):
    return {"grid": [0.0, 1.0], "values": [float(v), float(v)]}


def resample(x, y, xt):
    """y(xt) of a closed airfoil listed TE -> LE -> TE: each side is interpolated on its own."""
    ile = int(np.argmin(x))
    jle = int(np.argmin(xt))
    up = np.interp(xt[: jle + 1], x[: ile + 1][::-1], y[: ile + 1][::-1])
    lo = np.interp(xt[jle:], x[ile:], y[ile:])
    return np.concatenate([up, lo[1:]])


def blend_airfoil(doc, eta):
    pos = doc["components"]["blade"]["outer_shape_bem"]["airfoil_position"]
    g, lab = pos["grid"], pos["labels"]
    i = int(np.clip(np.searchsorted(g, eta) - 1, 0, len(g) - 2))
    a, b = (next(x for x in doc["airfoils"] if x["name"] == n) for n in lab[i : i + 2])
    w = (eta - g[i]) / (g[i + 1] - g[i]) if g[i + 1] > g[i] else 0.0
    xa, ya = np.asarray(a["coordinates"]["x"]), np.asarray(a["coordinates"]["y"])
    xb, yb = np.asarray(b["coordinates"]["x"]), np.asarray(b["coordinates"]["y"])
    yb = resample(xb, yb, xa)  # b onto the x stations of a, side by side (TE -> LE -> TE)
    blend = {
        "name": BLEND,
        "coordinates": {"x": list(a["coordinates"]["x"]), "y": [float(v) for v in (1 - w) * ya + w * yb]},
        "relative_thickness": float((1 - w) * a["relative_thickness"] + w * b["relative_thickness"]),
        "aerodynamic_center": a["aerodynamic_center"],
        "polars": a["polars"],
    }
    doc["airfoils"].append(blend)
    print(f"airfoil blend {lab[i]} / {lab[i + 1]} w={w:.4f} rthick={blend['relative_thickness']:.4f}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("src")
    ap.add_argument("dst")
    ap.add_argument("--eta", type=float, default=0.8)
    a = ap.parse_args()
    eta = a.eta
    doc = yaml.safe_load(open(a.src))
    blade = doc["components"]["blade"]
    osb = blade["outer_shape_bem"]
    blend_airfoil(doc, eta)
    osb["airfoil_position"] = {"grid": [0.0, 1.0], "labels": [BLEND, BLEND]}
    osb["chord"] = const(at(osb["chord"], eta))
    osb["pitch_axis"] = const(at(osb["pitch_axis"], eta))
    osb["twist"] = const(0.0)
    ra = osb["reference_axis"]
    ra["x"] = {"grid": list(ra["x"]["grid"]), "values": [0.0] * len(ra["x"]["grid"])}
    ra["y"] = {"grid": list(ra["y"]["grid"]), "values": [0.0] * len(ra["y"]["grid"])}

    st = blade["internal_structure_2d_fem"]
    layers = []
    for layer in st["layers"]:
        g = layer["thickness"]["grid"]
        if not g[0] - 1e-12 <= eta <= g[-1] + 1e-12:
            continue
        for k, d in list(layer.items()):
            if isinstance(d, dict) and "grid" in d and "values" in d:
                layer[k] = const(at(d, eta))
        layers.append(layer)
    st["layers"] = layers
    used = {x["web"] for x in layers if "web" in x}
    st["webs"] = [w for w in st["webs"] if w["name"] in used]
    for w in st["webs"]:
        for k, d in list(w.items()):
            if isinstance(d, dict) and "grid" in d and "values" in d:
                w[k] = const(at(d, eta))
    yaml.safe_dump(doc, open(a.dst, "w"), sort_keys=False)

    # check: reload and compare at 0.1 / 0.5 / 0.9
    c = yaml.safe_load(open(a.dst))["components"]["blade"]
    cs = c["outer_shape_bem"]
    for e in (0.1, 0.5, 0.9):
        n = len(c["internal_structure_2d_fem"]["layers"])
        print(
            f"eta {e}: chord {at(cs['chord'], e):.4f} twist {at(cs['twist'], e):.1f} pa {at(cs['pitch_axis'], e):.4f} "
            f"airfoil {cs['airfoil_position']['labels'][0]} layers {n} webs {len(c['internal_structure_2d_fem']['webs'])}"
        )
    print(f"tip z {at(cs['reference_axis']['z'], 1.0):.4f}; layers kept {len(layers)}/{len(blade['internal_structure_2d_fem']['layers'])}")


if __name__ == "__main__":
    main()
