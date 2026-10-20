"""Mesh-free principal EI of the WindIO layup (thin-wall Steiner integral).

Only NumPy/PyYAML/SciPy plus tools/yaml_layup_independent.py (itself mesh-free).
Every layer is a thin strip on the airfoil contour (nd_arc interval, TE-anchored
width, or straight web segment) with Ex(theta) * t * ds. EA-weighted centroid,
then eigenvalues of the 2x2 bending block. Shear lag, warping and the offset of a
layer from the contour are ignored; the 4 % mid-surface vs OML placement effect is
therefore NOT in this number (it is a lower-order thin-wall estimate).

    python -B tools/yaml_ei_independent.py [n_stations]
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from yaml_layup_independent import contour, interp, load, moduli  # noqa: E402


def strips(data, eta, resolution=4001):
    blade = data["components"]["blade"]
    structure = blade["internal_structure_2d_fem"]
    u, xy, perimeter = contour(data, eta, resolution=resolution)
    ds = np.r_[np.diff(u) * perimeter, 0.0]
    mid_xy = xy
    mats = {m["name"]: m for m in data["materials"]}
    webs = {}
    for web in structure["webs"]:
        a, b = [np.array([np.interp(interp(web[k], eta), u, xy[:, c]) for c in (0, 1)])
                for k in ("start_nd_arc", "end_nd_arc")]
        webs[web["name"]] = (a, b)
    out = []  # (EA, x, y) point masses of Ex*t*ds
    for layer in structure["layers"]:
        h = interp(layer["thickness"], eta)
        mat = mats[layer["material"]]
        ex = moduli(mat, interp(layer["fiber_orientation"], eta))[1]
        if "web" in layer:
            a, b = webs[layer["web"]]
            n = 400
            t = (np.arange(n) + .5) / n
            length = np.linalg.norm(b - a)
            for ti in t:
                p = a + ti * (b - a)
                out.append((ex * h * length / n, p[0], p[1]))
            continue
        if "start_nd_arc" in layer and "end_nd_arc" in layer:
            s, e = interp(layer["start_nd_arc"], eta), interp(layer["end_nd_arc"], eta)
        else:
            w = interp(layer["width"], eta) / perimeter
            s, e = (0.0, w) if "start_nd_arc" in layer else (1.0 - w, 1.0)
        sel = (u >= s) & (u <= e)
        for x, y, d in zip(mid_xy[sel, 0], mid_xy[sel, 1], ds[sel]):
            out.append((ex * h * d, x, y))
    return np.array(out)


def principal_ei(data, eta):
    p = strips(data, eta)
    ea = p[:, 0].sum()
    xc, yc = (p[:, 0] * p[:, 1]).sum() / ea, (p[:, 0] * p[:, 2]).sum() / ea
    dx, dy = p[:, 1] - xc, p[:, 2] - yc
    block = np.array([[(p[:, 0] * dy * dy).sum(), -(p[:, 0] * dx * dy).sum()],
                      [-(p[:, 0] * dx * dy).sum(), (p[:, 0] * dx * dx).sum()]])
    w = np.linalg.eigvalsh(block)
    return w[0], w[1], ea  # (small=flap, large=edge: chord >> thickness)


if __name__ == "__main__":
    from beam_reference import load_beamdyn_blade
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    root = Path(__file__).resolve().parents[1]
    data = load()
    bd = load_beamdyn_blade(root / "tests/IEA15MW/reference/IEA-15-240-RWT_BeamDyn_blade.dat")
    span = 117.0
    rf, re = [], []
    for eta in np.linspace(.25, .8, n):
        f, e, ea = principal_ei(data, eta)
        bf, be = np.interp(eta * bd.z.max(), bd.z, bd.ei_flap_ec), np.interp(eta * bd.z.max(), bd.z, bd.ei_edge_ec)
        rf.append(f / bf - 1); re.append(e / be - 1)
        print(f"eta={eta:.3f} EI_flap yaml/6x6={f/bf-1:+.3%} EI_edge={e/be-1:+.3%} EA={ea/np.interp(eta*bd.z.max(), bd.z, bd.ea)-1:+.3%}")
    print(f"MEAN yaml-layup vs 6x6: flap {100*np.mean(rf):+.2f} %  edge {100*np.mean(re):+.2f} %")
