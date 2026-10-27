#!/usr/bin/env python
"""S-2 tip metric on prismatic thin-wall tubes against the Timoshenko cantilever (l3-closure T10a).

The S-2 test judged one tip node, ``argmax(z)``: the first tip-ring node in mesh order. The tip force is split
equally over the tip-ring nodes (``tools/run_s2_static_cases.py``), so the beam-equivalent tip translation is the
ring mean. This tool checks that claim where the answer is known: a clamped prismatic tube under the same load
introduction, against ``P L^3 / (3 EI) + P L / (G A_s)`` with the thin-wall mid-line EI (exact segment integral)
and the shear area of the walls parallel to the load (rect) or A/2 (circle).

Usage:  scripts/aeroenv.sh python tools/s2_ring_metric_tube.py [--length 117] [--force 1e6]
Prediction and verdict: odd/tasks/l3-closure.md "T10a prediction".
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))

import diagnose_s1_gj_tube as tube  # noqa: E402
from run_s2_static_cases import run_shell_case, saint_venant_tip, tip_ring_scale  # noqa: E402

CASES = {  # name: (shape, load axis index in the section plane: 0 = x, 1 = y)
    "rect-flap": ("rect", 1),  # load along the 2 m walls (H)
    "rect-edge": ("rect", 0),  # load along the 4 m walls (B)
    "circle": ("circle", 1),
}


def section_constants(shape, axis):
    """Thin-wall mid-line I about the neutral axis normal to ``axis`` [m^3 * t] and shear area / t [m]."""
    p = np.asarray(tube.ring_points(shape))
    q = np.roll(p, -1, axis=0)
    a, b = p[:, axis], q[:, axis]
    seg = np.linalg.norm(q - p, axis=1)
    centroid = np.sum(seg * (a + b) / 2.0) / np.sum(seg)
    a, b = a - centroid, b - centroid
    i_over_t = np.sum(
        seg * (a * a + a * b + b * b) / 3.0
    )  # exact for a linear coordinate on each segment
    if shape == "circle":
        shear_over_t = 0.5 * np.sum(seg)
    else:  # walls parallel to the load carry the shear
        d = q - p
        shear_over_t = np.sum(seg[np.abs(d[:, axis]) > np.abs(d[:, 1 - axis])])
    return i_over_t, shear_over_t


def run(name, length, force):
    shape, axis = CASES[name]
    n_layers = int(round(length / tube.ES))
    mesh, props, _zl = tube.build_tube(
        np.full(n_layers, tube.G_UNIFORM), shape, tube.T_WALL, length
    )
    f = [0.0, 0.0, 0.0]
    f[axis] = force
    u = np.asarray(run_shell_case(mesh, props, {"type": "tip", "value": tuple(f)})).reshape(-1, 6)

    coords = mesh.coords_array
    ring = np.nonzero(coords[:, 2] >= coords[:, 2].max() - 1e-6)[0]
    node = int(np.argmax(coords[:, 2]))  # the metric the S-2 test used
    u_ring = u[ring, axis]
    ring_mean = float(u_ring.mean())

    e_mod = 2.0 * (1.0 + tube.NU) * tube.G_UNIFORM
    i_over_t, shear_over_t = section_constants(shape, axis)
    ei, gas = e_mod * i_over_t * tube.T_WALL, tube.G_UNIFORM * shear_over_t * tube.T_WALL
    bend, shear = force * length**3 / (3.0 * ei), force * length / gas
    ref = bend + shear

    print(f"== {name}: {len(mesh.nodes)} nodes, ring {len(ring)} nodes, load axis {'xy'[axis]}")
    print(
        f"   closed form {ref:.6f} m (bending {bend:.6f}, shear {shear:.6f}; EI {ei:.4e}, G A_s {gas:.4e})"
    )
    print(f"   ring mean   {ring_mean:.6f} m  ({ring_mean / ref - 1:+.4%} vs closed form)")
    print(
        f"   node {node:6d} {u[node, axis]:.6f} m  ({u[node, axis] / ring_mean - 1:+.4%} vs ring mean)"
    )
    print(f"   ring spread (max-min)/mean {(u_ring.max() - u_ring.min()) / ring_mean:+.4%}")
    # Saint-Venant decay: ring mean inboard of the loaded ring vs the closed-form deflection at that z
    depth = 2.0 * max(abs(p[axis]) for p in tube.ring_points(shape))  # section depth along the load
    z_all = coords[:, 2]
    for d in (0.5, 1.0, 2.0, 4.0):
        z_t = length - d * depth
        z_s = z_all[np.argmin(np.abs(z_all - z_t))]
        sel = np.nonzero(np.abs(z_all - z_s) < 1e-6)[0]
        u_s = float(u[sel, axis].mean())
        ref_s = force * z_s**2 * (3.0 * length - z_s) / (6.0 * ei) + force * z_s / gas
        print(
            f"   ring mean at {d:.1f} depth(s) inboard (z {z_s:.2f}): {u_s:.6f} m  ({u_s / ref_s - 1:+.4%} vs closed form)"
        )
    scale = tip_ring_scale(mesh)
    u_sv = saint_venant_tip(mesh, u.ravel(), scale)[axis]
    print(
        f"   Saint-Venant tip (rings 1-3 x {scale:.2f} m inboard, cubic): {u_sv:.6f} m  ({u_sv / ref - 1:+.4%} vs closed form)"
    )
    return ring_mean / ref - 1, u[node, axis] / ring_mean - 1


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--length", type=float, default=tube.LENGTH)
    ap.add_argument("--force", type=float, default=1.0e6)
    ap.add_argument("--es", type=float, default=tube.ES, help="element size around and along [m]")
    ap.add_argument("--case", nargs="*", default=list(CASES), choices=list(CASES))
    args = ap.parse_args()
    tube.ES = args.es  # ring_points and build_tube read the module constant
    for name in args.case:
        run(name, args.length, args.force)


if __name__ == "__main__":
    main()
