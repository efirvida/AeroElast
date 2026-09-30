#!/usr/bin/env python3
"""Blade AeroElast (linear MITC4) vs CalculiX (S8R): mesh convergence.

The blade parity tests compare two *different* formulations on the same mesh, so a
fixed-mesh gap says nothing on its own -- it could be discretisation or a model
difference.  This instrument answers that by refining the mesh and checking
whether the two sides tend to the same value.  Both sides get the same mesh;
CalculiX gets the quadratic element of it.

Measured 2026-09-30 (loads: 1e5 N flap in +y, 1e5 N edge in +x, 1e6 N axial in +z,
all applied at the tip section, clamped root):

    element_size  nodes   flap          edge          axial
    2.0            1460   10.56%        11.60%        12.43%
    1.0            3043    6.22%         5.15%         7.43%
    0.5            9277    3.68%         1.56%         4.90%

The gap shrinks monotonically for every case, so the two formulations do converge
to the same limit and the ~5-7% at 1.0 m is discretisation, not a model
difference.  A stringer check is still wanted: the node-averaged tip metric used
here moves with the tip node distribution (flap reads 7.55 / 7.94 / 7.12 m across
the three meshes), so the trend of the *gap* is the signal and the individual
values carry a few percent of metric noise.  The mesh-independent version of this
check compares the work done by the load, which is a global scalar.

Usage:

    python tools/blade_ccx_convergence.py [--element-sizes 2.0,1.0,0.5] [--ccx PATH]
"""

from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "tests"))

from _aeroelast import PyMeshAssembler  # noqa: E402

from aeroelast.core.mesh.entities import NodeSet  # noqa: E402
from aeroelast.core.mesh.io.writers import write_ccx_mesh  # noqa: E402
from aeroelast.models.blade.model import Blade  # noqa: E402

import test_blade_ccx_parity as blade_parity  # noqa: E402
from test_blade_iea15mw_validation import SPAN_DIRECTION, _to_rust_mesh  # noqa: E402

BLADE_YAML = REPO / "tests" / "IEA-15-240-RWT.yaml"
CASES = {
    "flap": [0.0, 1e5, 0.0],
    "edge": [1e5, 0.0, 0.0],
    "axial": [0.0, 0.0, 1e6],
}


def _tip_displacement(K, n_dofs, mesh, tip_idx, load) -> float:
    f = np.zeros(n_dofs)
    for i in tip_idx:
        for d in range(3):
            f[i * 6 + d] = load[d] / len(tip_idx)
    clamped = {
        mesh.node_id_to_index[nd.id] * 6 + d
        for nd in mesh.get_node_set("RootNodes").nodes.values()
        for d in range(6)
    }
    mask = np.ones(n_dofs, dtype=bool)
    mask[list(clamped)] = False
    free = np.where(mask)[0]
    u = np.zeros(n_dofs)
    u[free] = spsolve(K[np.ix_(free, free)], f[free])
    u = u.reshape(-1, 6)
    return float(np.mean(np.linalg.norm(u[tip_idx, :3], axis=1)))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--element-sizes", default="2.0,1.0,0.5")
    ap.add_argument("--ccx", default=shutil.which("ccx") or "ccx")
    ap.add_argument("--work-dir", type=Path, default=Path("blade_ccx_convergence"))
    args = ap.parse_args()

    sizes = [float(x) for x in args.element_sizes.split(",")]
    print(f"{'es':>5} {'nodes':>7} {'case':>6} {'AE (m)':>9} {'CCX (m)':>9} {'gap':>8}")
    for es in sizes:
        model = Blade(str(BLADE_YAML), element_size=es)
        model.generate_mesh()
        mesh, props = model.mesh, model.get_element_properties()
        coords = np.asarray([[nd.x, nd.y, nd.z] for nd in mesh.nodes])
        tip_idx = sorted(np.nonzero(coords[:, 2] > coords[:, 2].max() - 1e-6)[0])
        mesh.add_node_set(NodeSet("tip", {mesh.nodes[int(i)] for i in tip_idx}))

        asm = PyMeshAssembler.from_model(
            _to_rust_mesh(mesh, props), props, list(SPAN_DIRECTION), None
        )
        n_dofs = asm.dofs_count
        r, c, v = asm.assemble_k()
        K = coo_matrix((v, (r, c)), shape=(n_dofs, n_dofs)).tocsr()

        for name, load in CASES.items():
            ae = _tip_displacement(K, n_dofs, mesh, tip_idx, load)
            case_dir = args.work_dir / f"es{es}_{name}"
            shutil.rmtree(case_dir, ignore_errors=True)
            case_dir.mkdir(parents=True)
            inp = case_dir / f"blade_{name}.inp"
            write_ccx_mesh(
                mesh,
                str(inp),
                properties=props,
                boundary_nodeset="RootNodes",
                solver_type="LinearStatic",
                load_nodeset="tip",
                load_vector=[float(x) for x in load],
                span_direction=(0.0, 0.0, 1.0),
                quadratic=True,
            )
            proc = blade_parity.run_ccx(inp, args.ccx)
            frd = case_dir / f"blade_{name}.frd"
            if proc.returncode != 0 or not frd.exists():
                print(f"{es:5.2f} {mesh.node_count:7d} {name:>6} {ae:9.4f} {'CCX failed':>9}")
                continue
            ccx_disp = blade_parity._ccx_tip_displacements(frd, mesh, coords)
            ccx, _ = blade_parity._tip_metrics(mesh, ccx_disp, coords, sorted(ccx_disp))
            gap = abs(ae - ccx) / ccx
            print(f"{es:5.2f} {mesh.node_count:7d} {name:>6} {ae:9.4f} {ccx:9.4f} {100 * gap:7.2f}%")
            sys.stdout.flush()


if __name__ == "__main__":
    main()
