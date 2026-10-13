"""Shell vs CalculiX S8R on the Zhou-load blade case (Line A of odd/tasks/shell-twist-mixed-mesh.md).

Called from ``diagnose_zhou_loads_reverse.py`` when ``ZHOU_CCX_DIR`` is set. It receives the mesh,
properties and the nodal force vector the Zhou tool already built, writes the same case as a
CalculiX S8R deck (``load_field``, translations only), runs it, and compares our displacement
field with CalculiX's through the SAME ring decomposition. The gates were committed beforehand
(``bb4cc5d``): tip ``omega`` within 5 %, work ``W = 1/2 f.u`` within 5 %.
"""

from __future__ import annotations

import shutil
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

GATE = 0.05
MATCH_TOL = 2e-3  # [m] the FRD prints 5 significant digits: ~5e-4 m resolution at z ~ 100 m


def _ccx_displacement(frd: Path, coords: np.ndarray) -> np.ndarray:
    """CalculiX translations on OUR corner nodes, as a 6-DOF-strided vector (rotations zero)."""
    from tests.validation.parity.test_orthotropic_shell_parity import (
        _parse_ccx_frd,
        _parse_ccx_frd_coords,
    )

    ccx_xyz = _parse_ccx_frd_coords(frd)
    ids = np.fromiter(ccx_xyz.keys(), dtype=int)
    pts = np.array([ccx_xyz[i] for i in ids])
    disp = _parse_ccx_frd(frd, list(ids))
    dist, nn = cKDTree(pts).query(coords)
    if dist.max() > MATCH_TOL:
        raise RuntimeError(f"{int((dist > MATCH_TOL).sum())} of our nodes have no FRD match (max {dist.max():.2e} m)")
    if len(set(nn.tolist())) != len(nn):
        raise RuntimeError("ambiguous match: two of our nodes map to the same FRD node")
    u = np.zeros(6 * len(coords))
    for i, k in enumerate(nn):
        d = disp.get(int(ids[k]))
        if d is None:
            raise RuntimeError(f"node {i}: no displacement in the FRD")
        u[6 * i : 6 * i + 3] = d[:3]
    return u


def compare(mesh, props, f, u_ours, free, coords, ring_sets, outdir, ring_fn, span_direction):
    """Write, run and compare. ``ring_sets`` maps a label to the node-index array of that ring."""
    from aeroelast.core.mesh.io.writers import write_ccx_mesh
    from tests.support.ccx_io import run_ccx

    ccx_bin = shutil.which("ccx") or shutil.which("CalculiX")
    if ccx_bin is None:
        raise SystemExit("ccx not found in PATH")
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    load_field = {
        mesh.nodes[i].id: [float(f[6 * i]), float(f[6 * i + 1]), float(f[6 * i + 2])]
        for i in range(len(mesh.nodes))
        if np.any(f[6 * i : 6 * i + 3] != 0.0)
    }
    inp = outdir / "zhou_s8r.inp"
    write_ccx_mesh(
        mesh,
        str(inp),
        properties=props,
        boundary_nodeset="RootNodes",
        solver_type="LinearStatic",
        load_field=load_field,
        span_direction=tuple(float(x) for x in span_direction),
        shell_element_type="S8R",
    )

    # Precondition, before any result is read: the deck carries the same load.
    res_ours = f.reshape(-1, 6)[:, :3].sum(axis=0)
    res_deck = np.sum(list(load_field.values()), axis=0)
    rel = float(np.linalg.norm(res_deck - res_ours) / max(np.linalg.norm(res_ours), 1e-30))
    print(f"[ccx] load resultant ours {res_ours} deck {res_deck} rel diff {rel:.2e}")
    if rel > 1e-6:
        raise SystemExit("precondition failed: the CLOAD resultant differs from f")

    proc = run_ccx(inp, ccx_bin)
    print(f"[ccx] returncode {proc.returncode}; stdout tail: {proc.stdout[-300:]!r}")
    frd = inp.with_suffix(".frd")
    if proc.returncode != 0 or not frd.exists():
        raise SystemExit("CalculiX did not produce a result; nothing to compare")

    u_ccx = _ccx_displacement(frd, coords)

    w_ours = 0.5 * float(f[free] @ u_ours[free])
    w_ccx = 0.5 * float(f[free] @ u_ccx[free])
    print(f"[ccx] work W=1/2 f.u  ours {w_ours:.6e}  ccx {w_ccx:.6e}  ratio {w_ours / w_ccx:.4f}")

    print(f"[ccx] {'ring':>8} {'omega ours':>12} {'omega ccx':>12} {'rel diff':>9} {'dist/|om| ours':>15} {'ccx':>8}")
    gates = {}
    for label, ring in ring_sets.items():
        ko = ring_fn(coords, u_ours, ring)
        kc = ring_fn(coords, u_ccx, ring)
        wo, wc = float(np.rad2deg(ko["omega"])), float(np.rad2deg(kc["omega"]))
        so, sc = ko["distortion"] / abs(ko["omega"]), kc["distortion"] / abs(kc["omega"])
        gates[label] = (wo, wc)
        print(f"[ccx] {label:>8} {wo:12.4f} {wc:12.4f} {(wo - wc) / wc:+9.2%} {so:15.3f} {sc:8.3f}")

    wo, wc = gates["tip"]
    d_tip = abs(wo - wc) / abs(wc)
    d_w = abs(w_ours - w_ccx) / abs(w_ccx)
    print(f"[ccx] GATE tip omega {d_tip:.2%} (<= {GATE:.0%}): {'PASS' if d_tip <= GATE else 'FAIL'}")
    print(f"[ccx] GATE work      {d_w:.2%} (<= {GATE:.0%}): {'PASS' if d_w <= GATE else 'FAIL'}")
