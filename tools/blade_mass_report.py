"""Blade mass report of the windIO 2.x layup, mesh-free (independent evaluator).

Run: python -B tools/blade_mass_report.py [--yaml FILE] [--stations N] [--output FILE]

Per layer mass per length at N stations (tools/yaml_layup_independent), the layup mass per length against the
declared elastic_properties mass, and the blade mass: both integrated over the reference-axis arc length.
Prints to stdout; --output also writes the text. (The pre-T18 version summarised a NuMAD-mesh diagnostic JSON under
logs/check/mass; that mesher and its JSON are retired, so this report no longer depends on any mesh artefact.)
"""

import argparse
from pathlib import Path

import numpy as np

from yaml_layup_independent import DEFAULT_YAML, blade_length, blade_mass, load, section


def report(data, stations=12):
    lines = [f"blade length (reference-axis arc) = {blade_length(data):.4f} m",
             "eta layer material thickness[m] width[m] mass[kg/m]"]
    totals = {}
    for eta in np.linspace(0.2, 0.85, stations):
        r = section(data, eta)
        for layer in r["layers"]:
            lines.append(f"{eta:.6f} {layer['name']} {layer['material']} {layer['thickness_m']:.9f} "
                         f"{layer['width_m']:.6f} {layer['mass']:.6f}")
            totals.setdefault(layer["material"], []).append(layer["mass"] / r["beam_mass"])
        lines.append(f"{eta:.6f} TOTAL layup={r['mass']:.5f} declared={r['beam_mass']:.5f} kg/m "
                     f"({100*(r['mass']/r['beam_mass']-1):+.4f} %)")
    for material, v in sorted(totals.items()):
        lines.append(f"MATERIAL_MEAN {material} contribution_pct_of_declared={100*np.sum(v)/stations:.6f}")
    layup, declared = blade_mass(data)
    lines.append(f"BLADE_MASS layup={layup:.1f} kg declared(elastic_properties)={declared:.1f} kg "
                 f"ratio-1={100*(layup/declared-1):+.3f} %")
    return lines


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--yaml", type=Path, default=DEFAULT_YAML)
    parser.add_argument("--stations", type=int, default=12)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    lines = report(load(args.yaml), args.stations)
    print("\n".join(lines))
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
