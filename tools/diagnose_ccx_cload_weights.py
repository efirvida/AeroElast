"""Does the nodal weighting of a distributed ``*CLOAD`` on S8R change the CalculiX response?

Clamped strip (L = 1, b = 0.1, h = 1e-3 steel), small dead tip load alpha = P L^2 / EI = 0.05 on the
free edge, linear step. Same total force P for every weighting; the tip-edge mean deflection is
compared with the Euler-Bernoulli ``P L^3 / (3 EI)`` (the strip is a plate: 1-5 % stiffer) and with
the AeroElast MITC4 solution under the consistent load. A load lost by the weighting would show as a
deflection ratio well below 1 against the consistent weighting.

Usage: python tools/diagnose_ccx_cload_weights.py [--n 32] [--nv 2]
"""

from __future__ import annotations

import argparse
import tempfile
from pathlib import Path

import numpy as np

from tests.conftest import ccx_bin_or_skip
from tests.support.ccx_io import fail_ccx, run_ccx
from tests.support.ccx_surface import node_id, read_dat_displacements, write_s8r_surface_deck
from tests.support.strip_cantilever import E, H, L

B, NU, ALPHA = 0.1, 0.3, 0.05
EI = E * B * H**3 / 12.0


def weights(kind: str, nv: int) -> np.ndarray:
    m = 2 * nv + 1
    if kind == "consistent":
        w = np.array([1.0 if j in (0, m - 1) else (4.0 if j % 2 else 2.0) for j in range(m)])
    elif kind == "equal":
        w = np.ones(m)
    elif kind == "corners_only":
        w = np.array([1.0 if j % 2 == 0 else 0.0 for j in range(m)])
        w[1:-1:2] = 0.0
        w[2:-1:2] = 2.0
    else:
        raise ValueError(kind)
    return w / w.sum()


def ccx_tip_w(kind: str, n: int, nv: int, ccx_bin: str, workdir: Path) -> float:
    load = ALPHA * EI / L**2
    w = weights(kind, nv)
    inp = workdir / f"strip_{kind}_{n}_{nv}.inp"
    nodes = [node_id(2 * n, j, n) for j in range(2 * nv + 1)]
    write_s8r_surface_deck(
        inp,
        lambda u, v: (L * u, B * v, 0.0),
        n,
        nv,
        young=E,
        poisson=NU,
        thickness=H,
        boundary=lambda x, y, z: range(1, 7) if abs(x) < 1e-12 else (),
        nodal_loads={nid: (3, load * w[j]) for j, nid in enumerate(nodes)},
        print_nodes=nodes,
    )
    result = run_ccx(inp, ccx_bin)
    if result.returncode != 0:
        fail_ccx(result, inp)
    disp = read_dat_displacements(inp.with_suffix(".dat"))
    return float(np.mean([abs(disp[nid][2]) for nid in nodes]))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n", type=int, default=32)
    ap.add_argument("--nv", type=int, default=2)
    args = ap.parse_args()
    ccx = ccx_bin_or_skip()
    beam = ALPHA * EI / L**2 * L**3 / (3.0 * EI)
    with tempfile.TemporaryDirectory() as tmp:
        res = {k: ccx_tip_w(k, args.n, args.nv, ccx, Path(tmp)) for k in ("consistent", "equal", "corners_only")}
    print(f"beam P L^3/(3EI) = {beam:.6e} m   (n={args.n}, nv={args.nv})")
    for k, v in res.items():
        print(f"{k:13s} w_tip = {v:.6e}  /beam = {v / beam:.4f}  /consistent = {v / res['consistent']:.4f}")


if __name__ == "__main__":
    main()
