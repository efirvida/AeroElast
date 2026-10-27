#!/usr/bin/env python3
"""Group 5 (blade modal mesh convergence vs CalculiX): pair the modes by shape, not by frequency.

``tests/validation/parity/test_blade_iea15mw_mesh_convergence.py`` pairs AeroElast and CCX
modes by the smallest one-to-one relative frequency distance (Hungarian on |f_ae - f_ccx| / f_ccx).
CCX modes 7 and 8 sit within a few percent of each other, so a frequency pairing can swap their
partners between meshes; T29 then read CCX mode 7 as "not converging" (0.89 / 1.01 / 1.01 %).

For every mesh this tool solves both sides as the test does (same mesh, same properties, CCX S8R
of it) and reports, per CCX mode:

* the frequency pairing (the test's) and the MAC pairing (Hungarian on 1 - MAC), with the gap and
  the MAC under each;
* the mode character of the AeroElast partner: share of the squared translations in x / y / z and
  the share of the rotation about the span axis (theta_z) among the rotations.

MAC uses the translations of the mesh (corner) nodes, which keep their index in the quadratic CCX
deck (``_build_quadratic_data`` stacks the corner nodes first); the tool checks the coordinates.

Usage (SLURM, not the login node):

    python tools/diagnose_g5_mode7.py --meshes 0.5 0.35 0.25 --modes 10 --out logs/check/l3/g5/mac.json
"""

from __future__ import annotations

import argparse
import json
import sys
import tempfile
from pathlib import Path

import numpy as np
from scipy.optimize import linear_sum_assignment

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from _aeroelast import PyMeshAssembler, modal_solve_coo  # noqa: E402
from tests.conftest import ccx_bin_or_skip  # noqa: E402
from tests.support.ccx_io import fail_ccx, parse_ccx_frequencies, run_ccx  # noqa: E402
from aeroelast.core.mesh.entities import MeshElement, Node  # noqa: E402
from aeroelast.core.mesh.io.writers import write_ccx_mesh  # noqa: E402
from aeroelast.models.blade.model import Blade  # noqa: E402

import tests.validation.blade.test_blade_iea15mw_validation as blade_val  # noqa: E402


def parse_frd_mode_translations(frd: Path, n_nodes: int, n_modes: int) -> np.ndarray:
    """Return (n_modes, n_nodes, 3) translations of nodes 1..n_nodes from every DISP block."""
    modes: list[np.ndarray] = []
    current: np.ndarray | None = None
    with open(frd, encoding="utf-8", errors="replace") as handle:
        for line in handle:
            stripped = line.strip()
            if stripped.startswith("-4") and "DISP" in stripped.upper():
                current = np.zeros((n_nodes, 3))
                continue
            if current is None:
                continue
            if stripped.startswith("-3"):
                modes.append(current)
                current = None
                if len(modes) == n_modes:
                    break
                continue
            if not stripped.startswith("-1"):
                continue
            node_id = int(line[3:13])
            if 1 <= node_id <= n_nodes:
                current[node_id - 1] = [float(line[13:25]), float(line[25:37]), float(line[37:49])]
    return np.array(modes)


def mac_matrix(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """MAC between the rows of ``a`` (n_a, m) and of ``b`` (n_b, m)."""
    num = np.abs(a @ b.T) ** 2
    den = np.outer(np.einsum("ij,ij->i", a, a), np.einsum("ij,ij->i", b, b))
    return num / np.maximum(den, 1e-300)


def solve_mesh(element_size: float, n_modes: int, workdir: Path, ccx_bin: str) -> dict:
    Node._id_counter = 0
    MeshElement._id_counter = 0
    blade = Blade(str(blade_val.YAML), element_size=element_size)
    blade.generate_mesh()
    mesh = blade.mesh
    props = blade.get_element_properties()
    coords = np.array([[nd.x, nd.y, nd.z] for nd in mesh.nodes])
    n_nodes = len(coords)

    assembler = PyMeshAssembler.from_model(
        blade_val._to_rust_mesh(mesh, props), props, list(blade_val.SPAN_DIRECTION), None
    )
    n = assembler.dofs_count
    k_rows, k_cols, k_vals = assembler.assemble_k()
    m_rows, m_cols, m_vals = assembler.assemble_m()
    root = {mesh.node_id_to_index[nid] for nid in mesh.get_node_set("RootNodes").node_ids}
    fixed = {6 * i + d for i in root for d in range(6)}
    free = np.array([i for i in range(n) if i not in fixed], dtype=np.int64)
    freqs, flat = modal_solve_coo(
        np.asarray(k_rows, dtype=np.int64),
        np.asarray(k_cols, dtype=np.int64),
        np.asarray(k_vals, dtype=np.float64),
        np.asarray(m_rows, dtype=np.int64),
        np.asarray(m_cols, dtype=np.int64),
        np.asarray(m_vals, dtype=np.float64),
        n,
        free,
        n_modes,
    )
    freqs_ae = np.asarray(freqs, dtype=float)
    vecs = np.asarray(flat, dtype=float).reshape(len(freqs_ae), len(free))
    full = np.zeros((len(freqs_ae), n))
    full[:, np.sort(free)] = vecs
    full = full.reshape(len(freqs_ae), n_nodes, 6)
    order = np.argsort(freqs_ae)
    freqs_ae, full = freqs_ae[order], full[order]

    inp = workdir / f"blade_es{element_size:.2f}.inp"
    write_ccx_mesh(
        mesh,
        str(inp),
        properties=props,
        boundary_nodeset="RootNodes",
        solver_type="Modal",
        num_modes=n_modes,
        quadratic=True,
        span_direction=blade_val.SPAN_DIRECTION,
    )
    # Corner nodes keep their index: check the first coordinates written to the deck.
    # The main deck *INCLUDEs the mesh file, which holds the *NODE block.
    with open(inp.with_suffix(".msh"), encoding="utf-8") as handle:
        lines = handle.readlines()
    start = next(i for i, ln in enumerate(lines) if ln.upper().startswith("*NODE,"))
    deck = np.array(
        [[float(v) for v in ln.split(",")[1:4]] for ln in lines[start + 1 : start + 1 + n_nodes]]
    )
    coord_err = float(np.max(np.abs(deck - coords)))
    if coord_err > 1e-4 * float(np.max(np.abs(coords))):
        raise RuntimeError(
            f"corner nodes are not first in the CCX deck (max coord err {coord_err})"
        )

    result = run_ccx(inp, ccx_bin)
    if result.returncode != 0:
        fail_ccx(result, inp)
    freqs_ccx = np.asarray(parse_ccx_frequencies(inp, n_modes=n_modes), dtype=float)
    shapes_ccx = parse_frd_mode_translations(inp.with_suffix(".frd"), n_nodes, n_modes)
    order_c = np.argsort(freqs_ccx)
    freqs_ccx, shapes_ccx = freqs_ccx[order_c], shapes_ccx[order_c]

    ae_t = full[:, :, :3].reshape(len(freqs_ae), -1)
    ccx_t = shapes_ccx.reshape(len(freqs_ccx), -1)
    mac = mac_matrix(ae_t, ccx_t)

    cost_f = np.abs(freqs_ae[:, None] - freqs_ccx[None, :]) / freqs_ccx[None, :]
    rf, cf = linear_sum_assignment(cost_f)
    rm, cm = linear_sum_assignment(1.0 - mac)
    pair_f = {int(c): int(r) for r, c in zip(rf, cf, strict=True)}
    pair_m = {int(c): int(r) for r, c in zip(rm, cm, strict=True)}

    span = int(np.argmax(np.abs(blade_val.SPAN_DIRECTION)))
    rows = []
    for j in range(len(freqs_ccx)):

        def character(i: int) -> dict:
            t = full[i, :, :3]
            r = full[i, :, 3:]
            ts = np.sum(t**2, axis=0)
            rs = np.sum(r**2, axis=0)
            return {
                "ae_mode": i + 1,
                "trans_share_xyz": [round(float(v), 3) for v in ts / ts.sum()],
                "rot_share_span": round(float(rs[span] / rs.sum()), 3),
            }

        i_f, i_m = pair_f[j], pair_m[j]
        rows.append(
            {
                "ccx_mode": j + 1,
                "f_ccx": float(freqs_ccx[j]),
                "by_frequency": {
                    **character(i_f),
                    "f_ae": float(freqs_ae[i_f]),
                    "gap": float((freqs_ae[i_f] - freqs_ccx[j]) / freqs_ccx[j]),
                    "mac": float(mac[i_f, j]),
                },
                "by_mac": {
                    **character(i_m),
                    "f_ae": float(freqs_ae[i_m]),
                    "gap": float((freqs_ae[i_m] - freqs_ccx[j]) / freqs_ccx[j]),
                    "mac": float(mac[i_m, j]),
                },
            }
        )
    return {
        "element_size": element_size,
        "nodes": n_nodes,
        "freqs_ae": freqs_ae.tolist(),
        "freqs_ccx": freqs_ccx.tolist(),
        "mac": np.round(mac, 3).tolist(),
        "rows": rows,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--meshes", type=float, nargs="+", default=[0.5, 0.35, 0.25])
    parser.add_argument("--modes", type=int, default=10)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    ccx_bin = ccx_bin_or_skip()
    report = []
    with tempfile.TemporaryDirectory() as tmp:
        for es in args.meshes:
            res = solve_mesh(es, args.modes, Path(tmp), ccx_bin)
            report.append(res)
            print(f"es={es:.2f} m ({res['nodes']} nodes)")
            print(
                "  CCX  f_ccx   | freq-pair: AE  gap%   MAC  | MAC-pair: AE  gap%   MAC  | char (x,y,z ; rot_span)"
            )
            for r in res["rows"]:
                bf, bm = r["by_frequency"], r["by_mac"]
                print(
                    f"  {r['ccx_mode']:3d}  {r['f_ccx']:.4f} | {bf['ae_mode']:3d} {bf['gap'] * 100:+6.2f} "
                    f"{bf['mac']:.3f} | {bm['ae_mode']:3d} {bm['gap'] * 100:+6.2f} {bm['mac']:.3f} | "
                    f"{bm['trans_share_xyz']} ; {bm['rot_share_span']}"
                )
            args.out.parent.mkdir(parents=True, exist_ok=True)
            args.out.write_text(json.dumps(report, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
