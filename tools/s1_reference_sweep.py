"""S-1 numbers (signed and absolute mid-span error vs the BeamDyn 6x6) for one skin reference.

Same mesh, extractors and window as tests/test_iea15mw_s1_sectional.py; only the
skin reference surface is selectable, so a layup change can be measured before/after
with one command.

    scripts/aeroenv.sh python tools/s1_reference_sweep.py <oml|mid> [element_size]
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))

from beam_reference import load_beamdyn_blade  # noqa: E402
from run_s1_sectional import ea_gj_from_static_response  # noqa: E402
from aeroelast.core.mesh.generators import BladeMesh  # noqa: E402
from aeroelast.models.blade.model import build_rust_properties  # noqa: E402
from aeroelast.postprocess.sectional import SectionalExtractor  # noqa: E402


def main():
    reference = sys.argv[1]
    size = float(sys.argv[2]) if len(sys.argv) > 2 else 0.25
    bd = load_beamdyn_blade(ROOT / "tests/IEA15MW/reference/IEA-15-240-RWT_BeamDyn_blade.dat")
    gen = BladeMesh(yaml_file=str(ROOT / "tests/IEA-15-240-RWT.yaml"),
                    element_size=size, skin_reference=reference)
    mesh = gen.generate(renumber="rcm")
    props = build_rust_properties(gen.numad_mesh_data)
    ext = SectionalExtractor(mesh, props)
    sec = ext.extract()
    z, ei_flap, ei_edge = ext.ei_from_static_response(props)
    z_eg, _, gj, ea = ea_gj_from_static_response(mesh, props, numad_data=gen.numad_mesh_data)
    rows = {
        "EI_flap": (z, ei_flap, np.interp(z, bd.z, bd.ei_flap_ec)),
        "EI_edge": (z, ei_edge, np.interp(z, bd.z, bd.ei_edge_ec)),
        "mass": (sec.z, sec.mass_dens, np.interp(sec.z, bd.z, bd.mass_dens)),
        "GJ": (z_eg, gj, np.interp(z_eg, bd.z, bd.gj)),
        "EA_clt": (z_eg, ea, np.interp(z_eg, bd.z, bd.ea)),
    }
    for name, (zz, value, ref) in rows.items():
        m = (zz > 0.2 * zz.max()) & (zz < 0.85 * zz.max())
        rel = value[m] / ref[m] - 1
        print(f"S1 {reference} es={size} {name:8s} signed={100*np.nanmean(rel):+7.2f} % "
              f"abs={100*np.nanmean(np.abs(rel)):6.2f} %", flush=True)


if __name__ == "__main__":
    main()
