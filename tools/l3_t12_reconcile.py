"""T12: reconcile S-1 EI_edge CLT (+5.1 %, L21) with the older -0.40 % (beamdyn-arbiter-midsurface.md T8).

Same function as tools/s1_clt_steiner_ei.py (clt_ei_table, deck term = elastic-centre principal EI of the
BeamDyn 6x6, i.e. eigenvalues of K[3:5,3:5]-k k^T/EA), but keeps every station and reports the mean under
different station weightings. Writes per-station rows to the --out csv.
    scripts/aeroenv.sh env PYTHONPATH=$PWD python -B tools/l3_t12_reconcile.py mid 0.25 logs/chord_scaled/t12/rows.csv
"""

import sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))
sys.path.insert(0, str(ROOT))
from s1_clt_steiner_ei import clt_ei_table  # noqa: E402
from beam_reference import load_beamdyn_blade  # noqa: E402
from yaml_ei_independent import principal_ei  # noqa: E402
from yaml_layup_independent import load  # noqa: E402
from aeroelast.core.mesh.generators import BladeMesh  # noqa: E402
from aeroelast.models.blade.model import build_rust_properties  # noqa: E402
from aeroelast.postprocess.sectional import SectionalExtractor  # noqa: E402

ref, size, out = sys.argv[1], float(sys.argv[2]), sys.argv[3]
import os  # noqa: E402

if (
    os.environ.get("T12_A11") == "eff"
):  # hoop-free axial stiffness 1/(A^-1)11 (run_s1_sectional.laminate_a11_eff)
    import s1_clt_steiner_ei as _m
    from run_s1_sectional import laminate_a11_eff

    _m.laminate_a11 = laminate_a11_eff
zmax = None
if Path(out).exists():  # reuse the per-station rows (the mesh+CLT pass takes ~3.5 min)
    r = np.loadtxt(out, delimiter=",")
    zmax = r[:, 0].max()
else:
    gen = BladeMesh(
        yaml_file=str(ROOT / "tests/IEA-15-240-RWT.yaml"), element_size=size, skin_reference=ref
    )
    mesh = gen.generate(renumber="rcm")
    ext = SectionalExtractor(mesh, build_rust_properties(gen.numad_mesh_data))
    tab = clt_ei_table(gen, mesh, ext)
    data = load()
    bd = load_beamdyn_blade(ROOT / "tests/IEA15MW/reference/IEA-15-240-RWT_BeamDyn_blade.dat")
    zmax = tab[:, 0].max()
    rows = []
    for z, f, e in tab:
        eta = z / zmax
        yf, ye, _ = principal_ei(data, eta)
        rows.append(
            (
                z,
                eta,
                f / np.interp(z, bd.z, bd.ei_flap_ec) - 1,
                e / np.interp(z, bd.z, bd.ei_edge_ec) - 1,
                e / np.interp(z, bd.z, bd.ei_edge) - 1,
                e / ye - 1,
            )
        )
    r = np.array(rows)
    np.savetxt(
        out, r, delimiter=",", header="z,eta,flap_vs_ec,edge_vs_ec,edge_vs_K44ref,edge_vs_yaml"
    )
m = (r[:, 1] > 0.2) & (r[:, 1] < 0.85)
s = r[m]
print(f"{ref} es={size}: stations total {len(r)}, in 0.2-0.85 {m.sum()}")
print(
    f"A) mean over all stations in window (tool as run in L21): edge {100 * s[:, 3].mean():+.2f} %  flap {100 * s[:, 2].mean():+.2f} %"
)
zg = np.linspace(0.2 * zmax, 0.85 * zmax, 32)
for name, col in (
    ("edge_vs_ec", 3),
    ("edge_vs_K44ref (raw K44, not elastic-centre)", 4),
    ("edge_vs_yaml", 5),
):
    u = np.interp(zg, s[:, 0], s[:, col])
    print(
        f"B) {name}: 32 uniform-z samples {100 * u.mean():+.2f} %, span-integral (trapz) {100 * np.trapezoid(s[:, col], s[:, 0]) / (s[-1, 0] - s[0, 0]):+.2f} %"
    )
print("C) per eta band, edge_vs_ec (count, mean %):")
for lo in np.arange(0.2, 0.85, 0.1):
    b = s[(s[:, 1] >= lo) & (s[:, 1] < lo + 0.1)]
    if len(b):
        print(f"   eta {lo:.1f}-{lo + 0.1:.1f}: n={len(b):3d} {100 * b[:, 3].mean():+.2f} %")
