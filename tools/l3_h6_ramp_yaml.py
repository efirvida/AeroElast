"""H6 deck-equivalent copy of a windIO v1 blade: partial-span layers and webs ramp like the 21-station deck.

For each layer (and web) whose thickness grid does not span [0, 1], the thickness across each deck interval
(0.05 step) that contains a grid end becomes the deck's linear interpolation: zero at the outer deck station,
the yaml value at the inner one. Every other key of the layer or web keeps its points and is extended to the
new ends with its end value held. Scratch input for odd/tasks/l3-closure.md (H6), never a test input.

usage: python tools/l3_h6_ramp_yaml.py <in.yaml> <out.yaml>
"""

import sys

import numpy as np
import yaml

DS = np.round(np.arange(0.0, 1.0001, 0.05), 10)
EPS = 1e-9


def ramp_ends(grid):
    g0, g1 = grid[0], grid[-1]
    a0 = DS[DS < g0 - EPS].max() if g0 > EPS else None
    b1 = DS[DS > g1 + EPS].min() if g1 < 1 - EPS else None
    return a0, b1


def ramp_thickness(t, a0, b1):
    g, v = np.asarray(t["grid"], float), np.asarray(t["values"], float)
    ng, nv = list(g), list(v)
    if a0 is not None:  # zero at a0, yaml value at the first deck station inside the grid
        b0 = DS[DS >= g[0] - EPS].min()
        keep = g >= b0 - EPS
        ng = [a0, b0] + [x for x in g[keep] if x > b0 + EPS]
        nv = [0.0, float(np.interp(b0, g, v))] + [
            y for x, y in zip(g[keep], v[keep], strict=True) if x > b0 + EPS
        ]
    g, v = np.asarray(ng), np.asarray(nv)
    if b1 is not None:  # yaml value at the last deck station inside the grid, zero at b1
        a1 = DS[DS <= g[-1] + EPS].max()
        keep = g <= a1 + EPS
        ng = [x for x in g[keep] if x < a1 - EPS] + [a1, b1]
        nv = [y for x, y in zip(g[keep], v[keep], strict=True) if x < a1 - EPS] + [
            float(np.interp(a1, g, v)),
            0.0,
        ]
    return {"grid": [float(x) for x in ng], "values": [float(y) for y in nv]}


def extend(d, lo, hi):
    g, v = list(d["grid"]), list(d["values"])
    if lo is not None and g[0] > lo + EPS:
        g, v = [lo] + g, [v[0]] + v
    if hi is not None and g[-1] < hi - EPS:
        g, v = g + [hi], v + [v[-1]]
    return {**d, "grid": [float(x) for x in g], "values": v}


def main():
    src, dst = sys.argv[1], sys.argv[2]
    doc = yaml.safe_load(open(src))
    st = doc["components"]["blade"]["internal_structure_2d_fem"]
    web_ends = {}
    for layer in st["layers"]:
        a0, b1 = ramp_ends(layer["thickness"]["grid"])
        if a0 is None and b1 is None:
            continue
        layer["thickness"] = ramp_thickness(layer["thickness"], a0, b1)
        for k, d in list(layer.items()):
            if k != "thickness" and isinstance(d, dict) and "grid" in d:
                layer[k] = extend(d, a0, b1)
        if "web" in layer:
            lo, hi = web_ends.get(layer["web"], (a0, b1))
            web_ends[layer["web"]] = (
                min(x for x in (lo, a0) if x is not None) if a0 is not None else lo,
                max(x for x in (hi, b1) if x is not None) if b1 is not None else hi,
            )
        print(
            f"{layer['name']:24s} {layer['thickness']['grid'][0]:.3f}-{layer['thickness']['grid'][-1]:.3f} "
            f"head {np.round(layer['thickness']['values'][:2], 4)} tail {np.round(layer['thickness']['values'][-2:], 4)}"
        )
    for web in st["webs"]:
        lo, hi = web_ends.get(web["name"], (None, None))
        for k, d in list(web.items()):
            if isinstance(d, dict) and "grid" in d:
                web[k] = extend(d, lo, hi)
        print(f"web {web['name']} -> {lo}-{hi}")
    yaml.safe_dump(doc, open(dst, "w"), sort_keys=False)


if __name__ == "__main__":
    main()
