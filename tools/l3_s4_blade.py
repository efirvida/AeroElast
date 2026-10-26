"""Parked (or rotating) S-4 blade modes of any windIO v1 yaml, same code path as the S-4 test.

python tools/l3_s4_blade.py --yaml tests/IEA-3.4-130-RWT.yaml --element-size 0.25 --out r.json

Mirrors tests/test_iea15mw_s4_rotating_modal.py: MITC4 shell, skin_reference="mid", clamped
root nodes, consistent mass, ``assemble_rotating`` (K + K_G + K_SP, parked by default) and
``classify_modes`` (beam-kinematics projection; 1F/1E/2F = lowest flap, lowest edge, second
flap). The JSON keeps every mode with its label, beam fraction and per-family amplitudes so
the identification can be audited.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import numpy as np

sys.dont_write_bytecode = True
REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "tools"))

from run_s4_rotating_modal import (  # noqa: E402
    assemble_rotating,
    beam_projection,
    classify_modes,
    nodal_translations,
    omega_from_rpm,
)


def build_assembly(yaml_file, element_size):
    from aeroelast.core.assembler import MeshAssembler
    from aeroelast.core.mesh.generators import BladeMesh
    from aeroelast.elements import ElementFamily
    from aeroelast.models.blade.model import build_rust_properties

    gen = BladeMesh(element_size=element_size, yaml_file=str(yaml_file), skin_reference="mid")
    mesh = gen.generate(renumber="rcm")
    props = build_rust_properties(gen.numad_mesh_data)
    model = {"elements": {"element_family": ElementFamily.SHELL,
                          "span_direction": (0.0, 0.0, 1.0), "properties": props}}
    return MeshAssembler(mesh, model), mesh


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--yaml", required=True)
    ap.add_argument("--element-size", type=float, required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--rpm", type=float, default=0.0)
    ap.add_argument("--n-modes", type=int, default=12, help="12 = the S-4 test")
    args = ap.parse_args()
    if args.element_size >= 0.3:
        ap.error("element_size >= 0.3 is not allowed (docs/blade_mesh_size_rule.md)")

    asm, mesh = build_assembly(Path(args.yaml), args.element_size)
    freqs, modes, free = assemble_rotating(asm, mesh, omega_from_rpm(args.rpm), n_modes=args.n_modes)
    labels = classify_modes(mesh, modes, free)
    xyz = np.asarray([n.coords for n in mesh.nodes], dtype=np.float64)
    rows = []
    for k, (f, lab) in enumerate(zip(freqs, labels, strict=True)):
        frac, energy, _ = beam_projection(xyz, nodal_translations(modes[k], free, len(xyz)))
        rows.append({"index": k, "f_hz": float(f), "label": lab, "beam_fraction": frac,
                     "amplitude": {fam: float(np.sqrt(e)) for fam, e in energy.items()}})
    rows.sort(key=lambda r: r["f_hz"])
    flaps = [r["f_hz"] for r in rows if r["label"] == "flap"]
    edges = [r["f_hz"] for r in rows if r["label"] == "edge"]
    result = {"yaml": str(args.yaml), "element_size": args.element_size, "rpm": args.rpm,
              "n_nodes": len(xyz), "n_elements": len(mesh.elements),
              "span_m": float(xyz[:, 2].max() - xyz[:, 2].min()),
              "1F": flaps[0] if flaps else None, "1E": edges[0] if edges else None,
              "2F": flaps[1] if len(flaps) > 1 else None, "modes": rows}
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(json.dumps(result, indent=2) + "\n")
    for r in rows:
        print(f"  f={r['f_hz']:9.5f} Hz  {r['label']:8s} beam_fraction {r['beam_fraction']:.3f}")
    print(f"1F {result['1F']}  1E {result['1E']}  2F {result['2F']}  (nodes {len(xyz)})")


if __name__ == "__main__":
    main()
