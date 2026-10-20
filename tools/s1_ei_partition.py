"""Per-station EI: static solve vs mesh-CLT Steiner (same mesh), with BeamDyn 6x6 for scale.

    scripts/aeroenv.sh python -B tools/s1_ei_partition.py <oml|mid> [element_size]
"""
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))
from s1_clt_steiner_ei import (BladeMesh, SectionalExtractor, build_rust_properties,  # noqa: E402
                               clt_ei_table, load_beamdyn_blade)

ref = sys.argv[1]
size = float(sys.argv[2]) if len(sys.argv) > 2 else 0.25
gen = BladeMesh(yaml_file=str(ROOT / "tests/IEA-15-240-RWT.yaml"), element_size=size, skin_reference=ref)
mesh = gen.generate(renumber="rcm")
props = build_rust_properties(gen.numad_mesh_data)
ext = SectionalExtractor(mesh, props)
z, flap, edge = ext.ei_from_static_response(props)
tab = clt_ei_table(gen, mesh, ext)
bd = load_beamdyn_blade(ROOT / "tests/IEA15MW/reference/IEA-15-240-RWT_BeamDyn_blade.dat")
print("PART z     solve/CLT flap  solve/CLT edge  solve/6x6 edge  twist_deg(bd)")
for k in np.linspace(.12, .95, 14):
    zz = k * z.max()
    cf, ce = np.interp(zz, tab[:, 0], tab[:, 1]), np.interp(zz, tab[:, 0], tab[:, 2])
    sf, se = np.interp(zz, z, flap), np.interp(zz, z, edge)
    be = np.interp(zz, bd.z, bd.ei_edge_ec)
    print(f"PART {zz:6.1f} {sf/cf-1:+10.2%} {se/ce-1:+14.2%} {se/be-1:+14.2%}", flush=True)
