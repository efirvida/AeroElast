"""Bounded S-6 instrumentation; no reference export or validation-store writes.

Run via scripts/aeroenv.sh python -u tools/diagnose_s6_divergence.py.
Default: official WindIO es0.25, existing S-5/S-6 CSV artifacts. The accepted
mean/std bounds remain 25%; this diagnostic does not promote a verdict.
Optional --eigenvalues 3 requests smallest algebraic reduced stiffness
values (NOT frequencies: translation and rotation DOFs have mixed units).
ARPACK work is capped by --eigen-maxiter; nonconvergence is reported, not
interpreted as evidence of positive stiffness. Budget a separate blade job.
"""
from __future__ import annotations

import argparse
import hashlib
import time

import numpy as np
from scipy.sparse.linalg import ArpackNoConvergence, eigsh, norm

from run_s6_oneway_dynamic import (
    DT, ETA_K, ETA_M, REPO, RPM_RATED, assemble_keff, build_assembly,
    newmark_tip_series, nodal_load_series,
)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--element-size", type=float, choices=(.15, .2, .25), default=.25,
                    help="official WindIO mesh size in metres (default .25)")
    ap.add_argument("--eigenvalues", type=int, choices=range(0, 7), default=0,
                    help="smallest algebraic K_rot eigenvalues; 0 skips (default)")
    ap.add_argument("--eigen-maxiter", type=int, default=300,
                    help="ARPACK iteration cap, each iteration can be expensive")
    args = ap.parse_args()
    import aeroelast
    import _aeroelast._aeroelast as extension
    from pathlib import Path

    print(f"Python package: {aeroelast.__file__}")
    origin = Path(extension.__file__)
    print(f"extension: {origin}; sha256={hashlib.sha256(origin.read_bytes()).hexdigest()}")
    start = time.monotonic()
    reference = REPO / "tests/IEA15MW/reference"
    loads = np.loadtxt(reference / "s6_ad_blade1_timeseries.csv", delimiter=",", skiprows=1)
    s5 = np.loadtxt(reference / "s5_ad_blade1_loads.csv", delimiter=",", skiprows=1)
    t, ed = loads[:, 0], loads[:, 1]
    if not np.allclose(np.diff(t), DT, rtol=0, atol=1e-9):
        raise ValueError("S-6 artifact must have uniform DT=0.01 s")
    print(f"official WindIO es={args.element_size}; {len(t)} samples {t[0]}..{t[-1]} s")
    asm, mesh = build_assembly(args.element_size)
    components = {}
    k, m, free = assemble_keff(asm, mesh, RPM_RATED * 2 * np.pi / 60,
                              components=components)
    c = ETA_K * k + ETA_M * m
    step_matrix = k + 4 / DT**2 * m + 2 / DT * c
    matrices = dict(components, Krot=k, M=m, C=c, Newmark_step=step_matrix)
    print(f"nodes={mesh.node_count}, free={len(free)}; assembly seconds={time.monotonic()-start:.1f}")
    for name, matrix in matrices.items():
        reduced = matrix[free][:, free].tocsr()
        size = norm(reduced)
        skew = norm(reduced - reduced.T) / max(size, np.finfo(float).tiny)
        print(f"{name}: reduced Frobenius={size:.6e}, relative skew={skew:.3e}, nnz={reduced.nnz}")
    if args.eigenvalues:
        reduced = k[free][:, free].tocsr()
        skew = norm(reduced - reduced.T) / max(norm(reduced), np.finfo(float).tiny)
        if skew > 1e-10:
            print("Krot eigenvalues skipped: nonsymmetric; do not silently symmetrize")
        else:
            try:
                values, vectors = eigsh(reduced, k=args.eigenvalues, which="SA", tol=1e-7,
                                        maxiter=args.eigen_maxiter,
                                        v0=np.ones(len(free)))
            except ArpackNoConvergence as exc:
                values, vectors = exc.eigenvalues, exc.eigenvectors
                print("ARPACK NOT CONVERGED: no certificate of the lowest eigenvalue")
            for i, value in enumerate(values):
                q = vectors[:, i]
                residual = np.linalg.norm(reduced @ q - value * q)
                print(f"Krot eigenvalue={value:.6e}, absolute Ritz residual={residual:.3e}")
                print("  component Rayleigh: " + ", ".join(
                    f"{name}={q @ (matrix[free][:, free] @ q):.6e}"
                    for name, matrix in matrices.items()))
    else:
        print("eigenvalues skipped (use --eigenvalues 3); norms alone do not establish stability")
    f = nodal_load_series(asm, mesh, s5[:, 0], loads[:, 2:52], loads[:, 52:102], s5[:, 3])
    force_norm = np.linalg.norm(f[:, free], axis=1)
    resultant = f.reshape(len(t), -1, asm.dofs_per_node)[:, :, :3].sum(axis=1)
    print(f"prescribed force L2 min/max={force_norm.min():.6e}/{force_norm.max():.6e} N")
    print(f"resultant XYZ min={resultant.min(axis=0)}, max={resultant.max(axis=0)} N")
    print("Loads are precomputed, independent of u/v/a: no displacement feedback in this tool.")
    first = {}
    records = []

    def observe(step, record):
        records.append(record)
        for threshold in (1., 10., 100., 1e3, 1e6):
            if record["translation_max"] > threshold and threshold not in first:
                first[threshold] = step
                print(f"first translation >{threshold:g} m: sample={step}, "
                      f"artifact_time={t[step]:.3f}, elapsed={(step+1)*DT:.3f} s")
        if step < 8 or step == len(t) - 1:
            print(f"step={step}: {record}")

    hist = newmark_tip_series(asm, mesh, k, m, free, f, DT, observer=observe)
    shell = hist[:, 1]
    print("Sample zero follows one step from static mean-load IC; artifact labels are unchanged.")
    print(f"OoP first8={shell[:8]}, last4={shell[-4:]}")
    for label, actual, ref in (("mean magnitude", abs(shell.mean()), abs(ed.mean())),
                               ("std", shell.std(), ed.std())):
        relative = abs(actual - ref) / ref
        print(f"{label}: shell={actual:.8e} m, ElastoDyn={ref:.8e} m, "
              f"relative_error={relative:.6%}, accepted_bound=25%")
    print(f"max balance residual={max(r['relative_balance'] for r in records):.3e}; "
          f"solve reasons={sorted({r['solve_reason'] for r in records})}")
    print(f"total seconds={time.monotonic()-start:.1f}")


if __name__ == "__main__":
    main()
