#!/usr/bin/env python
"""S-1 sectional metrics of any windIO v1 blade against its BeamDyn 6x6 (odd/tasks/l3-closure.md, T9c).

Same extraction as ``tests/test_iea15mw_s1_sectional.py`` (mid-surface skin, native yaml bay centres,
tangent-frame EI from two static solves, T40 wall-shear-flow GJ, CLT EA, station-membership mass), with
the yaml, the deck and the blade length as arguments, so the IEA-15 and a second blade run the same code.
Prints, per metric, the mid-span (0.2-0.85 of the station range) signed mean and mean |rel err| vs the deck,
GJ band means, and writes the per-station values to ``--out`` (JSON). Decides nothing.

  python tools/l3_s1_blade.py --yaml tests/IEA-3.4-130-RWT.yaml \\
      --beamdyn tests/IEA34/reference/IEA-3.4-130-RWT_BeamDyn_blade.dat --blade-length 63.0 \\
      --element-size 0.25 --out logs/check/l3/t9/s1-iea34-es025.json
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "tools"))

from beam_reference import load_beamdyn_blade  # noqa: E402
from run_s1_sectional import ea_gj_from_static_response, native_bay_centres  # noqa: E402

GJ_BANDS = ((0.2, 0.3), (0.3, 0.6), (0.6, 0.85))


def _summary(z, ours, ref, lo=0.2, hi=0.85):
    rel = ours / ref - 1.0
    mid = (z > lo * z.max()) & (z < hi * z.max())
    return float(np.nanmean(rel[mid])), float(np.nanmean(np.abs(rel[mid]))), rel


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--yaml", required=True)
    ap.add_argument("--beamdyn", required=True, help="BeamDyn blade file (6x6 per station)")
    ap.add_argument("--blade-length", type=float, required=True,
                    help="length that scales the deck eta (yaml reference_axis z at the tip)")
    ap.add_argument("--element-size", type=float, default=0.25)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    if args.element_size >= 0.3:
        raise SystemExit("element_size >= 0.3 is not converged (docs/blade_mesh_size_rule.md)")

    from aeroelast.core.mesh.generators import BladeMesh
    from aeroelast.models.blade.model import build_rust_properties
    from aeroelast.postprocess.sectional import SectionalExtractor

    bd = load_beamdyn_blade(args.beamdyn, blade_length=args.blade_length)
    gen = BladeMesh(yaml_file=args.yaml, element_size=args.element_size, skin_reference="mid")
    mesh = gen.generate(renumber="rcm")
    props = build_rust_properties(gen.numad_mesh_data)
    ext = SectionalExtractor(mesh, props)
    stations = native_bay_centres(args.yaml)
    sec = ext.extract(stations=stations)
    z, ei_flap, ei_edge = ext.ei_from_static_response(props, frame="tangent", stations=stations)
    z_eg, _ea_resp, gj, ea = ea_gj_from_static_response(
        mesh, props, numad_data=gen.numad_mesh_data, stations=stations)

    ip = lambda x, y: np.interp(x, bd.z, y)  # noqa: E731
    cases = {
        "EI_flap": (z, ei_flap, ip(z, bd.ei_flap_ec)),
        "EI_edge": (z, ei_edge, ip(z, bd.ei_edge_ec)),
        "GJ": (z_eg, gj, ip(z_eg, bd.gj)),
        "EA_CLT": (z_eg, ea, ip(z_eg, bd.ea)),
        "mass": (sec.z, sec.mass_dens, ip(sec.z, bd.mass_dens)),
    }
    out = {"yaml": args.yaml, "beamdyn": args.beamdyn, "blade_length": args.blade_length,
           "element_size": args.element_size, "n_nodes": int(len(mesh.coords_array)), "metrics": {}}
    print(f"S-1 {Path(args.yaml).name} es {args.element_size} nodes {out['n_nodes']}; "
          "mid-span 0.2-0.85 vs deck: signed mean % (mean |rel| %)")
    for name, (zz, ours, ref) in cases.items():
        signed, mabs, rel = _summary(zz, np.asarray(ours, float), np.asarray(ref, float))
        rec = {"signed_pct": 100 * signed, "mean_abs_pct": 100 * mabs,
               "z": np.asarray(zz).tolist(), "ours": np.asarray(ours).tolist(), "ref": np.asarray(ref).tolist()}
        line = f"  {name:8s} {100 * signed:+7.2f} ({100 * mabs:5.2f})"
        if name == "GJ":
            eta = np.asarray(zz) / args.blade_length
            bands = {f"{lo}-{hi}": float(100 * np.nanmean(rel[(eta > lo) & (eta < hi)])) for lo, hi in GJ_BANDS}
            rec["band_signed_pct"] = bands
            line += "  bands " + " ".join(f"{k} {v:+.1f}" for k, v in bands.items())
        out["metrics"][name] = rec
        print(line)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(json.dumps(out, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
