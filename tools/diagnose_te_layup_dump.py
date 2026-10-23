"""B3-3 K13.A3: glass_uni thickness of the TE panel sets in the mesh vs the yaml TE_reinforcement thickness.

    scripts/aeroenv.sh python -u -B tools/diagnose_te_layup_dump.py <element_size>
"""
import re
import sys
from pathlib import Path

import numpy as np
import yaml
from scipy.interpolate import PchipInterpolator

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from aeroelast.core.mesh.generators import BladeMesh  # noqa: E402

L = 117.0


def main():
    size = float(sys.argv[1])
    gen = BladeMesh(yaml_file=str(ROOT / "tests/IEA-15-240-RWT.yaml"), element_size=size, skin_reference="mid")
    mesh = gen.generate(renumber="rcm")
    y = yaml.safe_load(open(ROOT / "tests/IEA-15-240-RWT.yaml"))
    layers = {l["name"]: l for l in y["components"]["blade"]["internal_structure_2d_fem"]["layers"]}
    f = {k: PchipInterpolator(layers[k]["thickness"]["grid"], layers[k]["thickness"]["values"]) for k in
         ("TE_reinforcement_SS", "TE_reinforcement_PS")}
    coords = np.asarray([n.coords for n in mesh.nodes])
    idx = mesh.node_id_to_index
    secs = {s["elementSet"]: s for s in gen.numad_mesh_data["sections"]}
    rows = {}
    for name, es in mesh.element_sets.items():
        if name not in secs or not re.search(r"(HP|LP)_TE_(PANEL|FLAT)", name):
            continue
        ids = {idx[n] for e in es.elements for n in e.node_ids}
        z = coords[sorted(ids), 2].mean()
        uni = sum(h for m, h, a in secs[name]["layup"] if m == "glass_uni")
        tot = sum(h for m, h, a in secs[name]["layup"])
        side = "HP" if "HP" in name else "LP"
        rows.setdefault((round(z, 2), side, "PANEL" if "PANEL" in name else "FLAT"), []).append((uni, tot, name))
    print(f"TE-LAYUP mid es={size}: mesh glass_uni thickness [mm] (sum of layup) vs yaml pchip SS/PS [mm] at z/117")
    print("  z[m]  side kind  uni_mm  total_mm  yaml_SS  yaml_PS  set")
    for (z, side, kind), v in sorted(rows.items()):
        if not 50 <= z <= 102:
            continue
        s = z / L
        print(f"  {z:6.2f} {side} {kind:5s} {1e3 * v[0][0]:7.2f} {1e3 * v[0][1]:8.2f} {1e3 * float(f['TE_reinforcement_SS'](s)):8.2f} "
              f"{1e3 * float(f['TE_reinforcement_PS'](s)):8.2f}  {v[0][2]}")


if __name__ == "__main__":
    main()
