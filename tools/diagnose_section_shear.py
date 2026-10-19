"""A55 (what the MITC4 composite shear block consumes) vs the energy-equivalent Cs per blade section.

T7(iii) L17 of odd/tasks/rotating-local-tip-modes.md. IEA-15 official, es 0.25, oml:
  scripts/aeroenv.sh python -I tools/diagnose_section_shear.py "$PWD"
"""
import sys
from pathlib import Path
import numpy as np
REPO = Path(sys.argv[1])
sys.path.insert(0, str(REPO / "src"))
from aeroelast.core.mesh.generators import BladeMesh  # noqa: E402
from aeroelast.core.laminate import Laminate, Ply  # noqa: E402
from aeroelast.core.material import OrthotropicMaterial  # noqa: E402
gen = BladeMesh(element_size=0.25, yaml_file=str(REPO / "tests/IEA-15-240-RWT.yaml"), skin_reference="oml")
gen.generate(renumber="rcm")
nd = gen.numad_mesh_data
mats = {m["name"]: m for m in nd["materials"]}
def om(name):
    el = mats[name]["elastic"]
    E, nu, G = el["E"], el["nu"], el.get("G")
    if np.ndim(E) == 0:
        nu_v = float(nu if np.ndim(nu) == 0 else list(nu)[0])
        E, nu = (float(E),) * 3, (nu_v,) * 3
        G = (E[0] / (2 * (1 + nu_v)),) * 3
    elif G is None:
        G = tuple(E[i] / (2 * (1 + nu[i])) for i in range(3))
    return OrthotropicMaterial(name=name, E=tuple(map(float, E)), G=tuple(map(float, G)),
                               nu=tuple(map(float, nu)), rho=1.0)


rows = []
for sec in nd["sections"]:
    lay = sec["layup"]
    if len(lay) < 2:
        continue
    lam = Laminate([Ply(om(m), float(t), float(a)) for m, t, a in lay])
    a55, cs = lam.shear_stiffness_uncorrected()[0, 0], lam.Cs[0, 0]
    rows.append((sec["elementSet"], lam.total_thickness, a55 / cs, "/".join(m for m, _, _ in lay)))
r = np.array([x[2] for x in rows])
print(f"multi-ply sections {len(rows)}: A55/Cs median {np.median(r):.2f}, max {r.max():.1f}, >2: {(r > 2).sum()}")
for name, h, ratio, mats_ in rows:
    if "TE_PANEL" in name and any(f"_{k}_" in f"_{name}_" or name.startswith(f"09_{k}") for k in range(40, 50)):
        print(f"{name:28s} h={h * 1e3:6.1f} mm  A55/Cs={ratio:6.1f}  {mats_[:70]}")
