"""Project blade shell modes on beam kinematics to tell global modes from local ones.

``classify_modes`` (tools/run_s4_rotating_modal.py) labels a mode from the single highest-z
tip node. T12 of odd/tasks/beamdyn-arbiter-midsurface.md found a mode at 1.5325 Hz (oml,
es 0.25, 7.56 rpm) with a small tip amplitude that it labels "2nd flap". This tool says what
such modes are, before any classifier change (odd/tasks/s4-mode-classifier.md, T1).

Per span bin, the nodal translations are fitted by rigid-section motion
u_i = t + theta x (p_i - c), the plane-section kinematics of a beam. Per mode it prints:

  * beam fraction  sum |u_fit|^2 / sum |u|^2 over all nodes (1 = beam-like, 0 = no beam
    content; unweighted by mass),
  * the beam energy shares of flap (t_y), edge (t_x), axial (t_z) and torsion (theta_z x r),
  * where the residual (non-beam) motion peaks: span fraction and chordwise position,
  * the legacy tip-node label (frozen copy of the old rule), the current ``classify_modes``
    label, and the tip-node components.

Modes are solved once and written to an .npz, so the projection can be redone offline:

  python tools/diagnose_mode_beam_projection.py solve --ref oml --es 0.25 --rpm 0 7.56 --out X.npz
  python tools/diagnose_mode_beam_projection.py analyse X.npz [--bins 200]
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from types import SimpleNamespace

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "tools"))

from run_s4_rotating_modal import beam_projection, classify_modes, nodal_translations  # noqa: E402

def legacy_tip_label(xyz, u):
    """The pre-2026-10-18 rule of ``classify_modes``: highest-z node, 2:1 dominance (frozen copy)."""
    ux, uy, uz = np.abs(u[int(np.argmax(xyz[:, 2]))])
    if uy >= 2.0 * ux and uy >= 2.0 * uz:
        return "flap"
    if ux >= 2.0 * uy and ux >= 2.0 * uz:
        return "edge"
    if uz >= 2.0 * ux and uz >= 2.0 * uy:
        return "axial"
    return "mixed"


def solve(args):
    from run_s4_rotating_modal import assemble_rotating, omega_from_rpm

    from aeroelast.core.assembler import MeshAssembler
    from aeroelast.core.mesh.generators import BladeMesh
    from aeroelast.elements import ElementFamily
    from aeroelast.models.blade.model import build_rust_properties

    gen = BladeMesh(
        element_size=args.es,
        yaml_file=str(REPO / "tests/IEA-15-240-RWT.yaml"),
        skin_reference=args.ref,
    )
    mesh = gen.generate(renumber="rcm")
    props = build_rust_properties(gen.numad_mesh_data)
    model = {"elements": {"element_family": ElementFamily.SHELL, "span_direction": (0.0, 0.0, 1.0), "properties": props}}
    asm = MeshAssembler(mesh, model)
    xyz = np.asarray([n.coords for n in mesh.nodes])
    out = {"xyz": xyz, "ref": args.ref, "es": args.es, "rpm": np.asarray(args.rpm)}
    for i, rpm in enumerate(args.rpm):
        freqs, modes, free = assemble_rotating(asm, mesh, omega_from_rpm(rpm), n_modes=args.n_modes)
        out[f"freqs_{i}"] = np.asarray(freqs)
        out[f"modes_{i}"] = np.asarray(modes)
        out[f"free_{i}"] = np.asarray(free)
        print(f"solved ref={args.ref} es={args.es} rpm={rpm}", flush=True)
    np.savez_compressed(args.out, **out)
    print(f"wrote {args.out}", flush=True)
    analyse_file(args.out, args.bins)


def analyse_file(path, n_bins):
    d = np.load(path, allow_pickle=False)
    xyz = d["xyz"]
    span = float(xyz[:, 2].max() - xyz[:, 2].min())
    for i, rpm in enumerate(d["rpm"]):
        freqs, modes, free = d[f"freqs_{i}"], d[f"modes_{i}"], d[f"free_{i}"]
        new = classify_modes(SimpleNamespace(nodes=[SimpleNamespace(coords=p) for p in xyz]), modes, free)
        print(f"\n=== ref={d['ref']} es={float(d['es'])} rpm={float(rpm)} bins={n_bins} ===")
        for k in np.argsort(freqs):
            u = nodal_translations(modes[k], free, len(xyz))
            frac, en, res = beam_projection(xyz, u, n_bins)
            tot = sum(en.values()) or 1.0
            rn = np.sum(res**2, axis=1)
            j = int(np.argmax(rn))
            tip = np.abs(u[int(np.argmax(xyz[:, 2]))])
            print(
                f"BEAM k={int(k):2d} {float(freqs[k]):9.6f} Hz legacy={legacy_tip_label(xyz, u):<5} new={new[k]:<7} frac={frac:6.3f}"
                f" flap={en['flap'] / tot:5.3f} edge={en['edge'] / tot:5.3f} tors={en['torsion'] / tot:5.3f}"
                f" axial={en['axial'] / tot:5.3f} | res peak z/L={(xyz[j, 2] - xyz[:, 2].min()) / span:5.3f}"
                f" x={xyz[j, 0]:+7.3f} y={xyz[j, 1]:+7.3f} | tip |u|=({tip[0]:.2e},{tip[1]:.2e},{tip[2]:.2e})",
                flush=True,
            )


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("solve")
    s.add_argument("--ref", choices=("oml", "mid"), required=True)
    s.add_argument("--es", type=float, required=True)
    s.add_argument("--rpm", type=float, nargs="+", default=[0.0, 7.56])
    s.add_argument("--n-modes", type=int, default=12)
    s.add_argument("--bins", type=int, default=200)
    s.add_argument("--out", type=Path, required=True)
    a = sub.add_parser("analyse")
    a.add_argument("npz", type=Path)
    a.add_argument("--bins", type=int, default=200)
    args = ap.parse_args()
    if args.cmd == "solve":
        solve(args)
    else:
        analyse_file(args.npz, args.bins)


if __name__ == "__main__":
    main()
