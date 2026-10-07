#!/usr/bin/env python3
"""Blade AeroElast (linear MITC4) vs CalculiX (S8R): mesh convergence.

The blade parity tests compare two *different* formulations on the same mesh, so a
fixed-mesh gap says nothing on its own -- it could be discretisation or a model
difference.  This instrument answers that by refining the mesh and checking
whether the two sides tend to the same value.  Both sides get the same mesh;
CalculiX gets the quadratic element of it.

Two scalars are reported per case:

* the **node-mean tip displacement**.  Cheap, but it moves with the tip node
  distribution, so across meshes it carries a few percent of metric noise and only
  its trend is meaningful;
* the **work done by the load**, ``W = 1/2 sum f_i . u_i``.  A global scalar, and
  the one to read: each side uses its own load representation (AeroElast spreads
  the resultant over the tip corner nodes, the writer spreads it over the corner
  nodes *and* the mid-side nodes of the quadratic conversion), and both
  representations converge to the same uniform traction, so the two compliances
  must converge to the same number.

Measured 2026-09-30 (loads: 1e5 N flap in +y, 1e5 N edge in +x, 1e6 N axial in +z,
applied at the tip section, clamped root).  Node-mean gap:

    element_size  nodes   flap      edge      axial
    2.0            1460   10.56%    11.60%    12.43%
    1.0            3043    6.22%     5.15%     7.43%
    0.5            9277    3.68%     1.56%     4.90%

The gap shrinks monotonically for every case, so the two formulations do converge
to the same limit and the ~5-7% at 1.0 m is discretisation.  The same is true of
the work-based gap, which is the one to read:

    element_size  nodes   flap              edge              axial
    2.0            1460   10.56% -> 10.70%  11.60% -> 11.65%  12.43% -> 14.55%
    1.0            3043    6.22% ->  6.38%   5.15% ->  5.17%   7.43% ->  9.52%
    0.5            9277    3.68% ->  3.87%   1.56% ->  1.55%   4.90% ->  7.73%

(the first number of each pair is the node-mean gap, the second the work gap.)
Every column falls with the mesh.  Axial settles slowest, which is expected: the
resultant is spread over the tip section, so its convergence carries the section's
local distortion as well as the global one.

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
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "tests"))

from _aeroelast import PyMeshAssembler  # noqa: E402

from aeroelast.core.mesh.entities import NodeSet  # noqa: E402
from aeroelast.core.mesh.io.writers import write_ccx_mesh  # noqa: E402
from aeroelast.models.blade.model import Blade  # noqa: E402

import test_blade_ccx_parity as blade_parity  # noqa: E402
from tests.validation.blade.test_blade_iea15mw_validation import (  # noqa: E402
    SPAN_DIRECTION,
    _to_rust_mesh,
)

BLADE_YAML = REPO / "tests" / "IEA-15-240-RWT.yaml"
CASES = {
    "flap": [0.0, 1e5, 0.0],
    "edge": [1e5, 0.0, 0.0],
    "axial": [0.0, 0.0, 1e6],
}


def _ae_solve(K, n_dofs, mesh, tip_idx, load):
    """Return (node-mean tip displacement, work done by the load)."""
    f = np.zeros(n_dofs)
    per_node = np.asarray(load, dtype=float) / len(tip_idx)
    for i in tip_idx:
        f[i * 6 : i * 6 + 3] = per_node
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
    mean = float(np.mean(np.linalg.norm(u.reshape(-1, 6)[tip_idx, :3], axis=1)))
    return mean, 0.5 * float(np.dot(f, u))


def _load_nset(nam_path: Path, name: str) -> list[int]:
    """Read the integer labels of one NSET from a CalculiX .nam file."""
    labels: list[int] = []
    collecting = False
    for line in nam_path.read_text().splitlines():
        stripped = line.strip()
        if stripped.upper().startswith("*NSET"):
            collecting = stripped.upper().endswith(f"N{name.upper()}")
            continue
        if stripped.startswith("*"):
            collecting = False
            continue
        if collecting and stripped:
            labels.extend(int(tok) for tok in stripped.split(",") if tok.strip())
    return labels


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--element-sizes", default="2.0,1.0,0.5")
    ap.add_argument("--ccx", default=shutil.which("ccx") or "ccx")
    ap.add_argument("--work-dir", type=Path, default=Path("blade_ccx_convergence"))
    args = ap.parse_args()

    sizes = [float(x) for x in args.element_sizes.split(",")]
    header = (
        f"{'es':>5} {'nodes':>7} {'case':>6} {'AE mean':>9} {'CCX mean':>9} {'gap':>8}"
        f" {'AE work':>12} {'CCX work':>12} {'gapW':>8}"
    )
    print(header)
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
            ae_mean, ae_work = _ae_solve(K, n_dofs, mesh, tip_idx, load)
            case_dir = args.work_dir / f"es{es}_{name}"
            shutil.rmtree(case_dir, ignore_errors=True)
            case_dir.mkdir(parents=True)
            stem = f"blade_{name}"
            inp = case_dir / f"{stem}.inp"
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
            frd = case_dir / f"{stem}.frd"
            nam = case_dir / f"{stem}.nam"
            if proc.returncode != 0 or not frd.exists():
                print(f"{es:5.2f} {mesh.node_count:7d} {name:>6} {ae_mean:9.4f} {'CCX failed':>9}")
                continue

            loaded = _load_nset(nam, "tip")
            disp = blade_parity._ccx_tip_displacements(frd, mesh, coords)
            ccx_mean, _ = blade_parity._tip_metrics(mesh, disp, coords, sorted(disp))
            # The work must integrate over every loaded node: the writer spreads the
            # resultant over the corner nodes AND the mid-side nodes the quadratic
            # conversion adds, so dividing by that count but summing only the
            # corners would count half the force and halve the work.
            from _ccx_io import parse_frd_disp  # noqa: PLC0415

            full = parse_frd_disp(frd, loaded)
            per_node = np.asarray(load, dtype=float) / max(len(loaded), 1)
            ccx_work = 0.5 * sum(float(np.dot(per_node, u)) for u in full.values())
            gap = abs(ae_mean - ccx_mean) / ccx_mean
            gap_w = abs(ae_work - ccx_work) / ccx_work
            print(
                f"{es:5.2f} {mesh.node_count:7d} {name:>6} {ae_mean:9.4f} {ccx_mean:9.4f} "
                f"{100 * gap:7.2f}% {ae_work:12.6e} {ccx_work:12.6e} {100 * gap_w:7.2f}%"
            )
            sys.stdout.flush()


if __name__ == "__main__":
    main()
