"""Summarise the mass diagnostic JSON, without importing aeroelast.

Run: python -B tools/blade_mass_report.py
Writes logs/check/mass/layer-report.txt. Delta decomposes exactly into coverage,
thickness and geometry terms; it does not infer SONATA internals from a deck.
"""

import json
from pathlib import Path

import numpy as np

from yaml_layup_independent import contour, interp, load

ROOT = Path(__file__).resolve().parents[1]


def main():
    evidence = json.loads((ROOT / "logs/check/mass/diagnostic.json").read_text())
    lines = ["kg/m; widths are actual mesh area / axial span; density=1940 kg/m^3",
             "Post-fix mesh_h is region-effective thickness (includes partial coverage).",
             "coverage/thickness/geometry is an algebraic decomposition, not three independent physical effects.",
             "reference eta layer yaml_h mesh_h_min mesh_h_max yaml_width mesh_width delta coverage thickness geometry"]
    for reference, result in evidence.items():
        terms = {}
        for row in result["stations"]:
            if not .2 <= row["eta"] <= .85:
                continue
            for layer in row["layers"]:
                if layer["material"] != "glass_triax":
                    continue
                mapped = layer["mapped"]
                h = layer["thickness_m"]
                w = mapped["width_m"]
                # Both affected skins are mapped over the full perimeter.
                # This decomposition is exact for those uniform-thickness layers.
                hm = mapped["mass"]/(1940*w)
                coverage = 1940*h*(row["perimeter_m"]-layer["width_m"])
                thickness = 1940*(hm-h)*w
                geometry = 1940*h*(w-row["perimeter_m"])
                terms.setdefault(layer["name"], []).append([coverage, thickness, geometry])
                lines.append(f"{reference} {row['eta']:.6f} {layer['name']} {h:.9f} "
                             f"{min(mapped['thicknesses_m']):.9f} {max(mapped['thicknesses_m']):.9f} "
                             f"{layer['width_m']:.6f} {w:.6f} {layer['delta_kg_m']:+.6f} "
                             f"{coverage:+.6f} {thickness:+.6f} {geometry:+.6f}")
        for name, values in terms.items():
            mean = np.mean(values, axis=0)
            lines.append(f"NATIVE_MEAN {reference} {name} coverage/thickness/geometry="
                         f"{mean[0]:+.6f}/{mean[1]:+.6f}/{mean[2]:+.6f} kg/m")
        yaml_ratios = [row["mass"]/row["beam_mass"]-1 for row in result["stations"]
                       if .2 <= row["eta"] <= .85]
        lines.append(f"NATIVE_MEAN yaml/beam signed={100*np.mean(yaml_ratios):+.6f}% "
                     "(native midpoints; not the prior 12-station mean)")
    data = load()
    structure = data["components"]["blade"]["internal_structure_2d_fem"]
    le = next(l for l in structure["layers"] if l["name"] == "LE_reinforcement")
    rho = next(m["rho"] for m in data["materials"] if m["name"] == le["material"])
    for eta in (3.5/49, 40.5/49):
        a, b = (interp(le[key], eta) for key in ("start_nd_arc", "end_nd_arc"))
        h = interp(le["thickness"], eta)
        _, _, perimeter = contour(data, eta)
        lines.append(f"YAML_INVALID_CHECK eta={eta:.9f} nd=[{a:.17g},{b:.17g}] "
                     f"h_m={h:.9g} signed_strip_mass_kg_m={rho*h*(b-a)*perimeter:.9g}; "
                     "diagnostic skips the independent comparison row, NOT a production mesh layer")
    original = json.loads((ROOT / "logs/check/layup/independent.json").read_text())
    delta = 100*np.mean([r["mass"]/r["beam_mass"]-1 for r in original])
    lines.append(f"ORIGINAL_12_STATION_MEAN yaml/beam signed={delta:+.6f}% (unmodified independent integral)")
    for material in sorted({l["material"] for r in original for l in r["layers"]}):
        fraction = 100*np.mean([sum(l["mass"] for l in r["layers"]
                                   if l["material"] == material)/r["beam_mass"]
                                for r in original])
        lines.append(f"MATERIAL_REMOVAL_MEAN {material} contribution_pct_of_beam={fraction:.6f} "
                     f"remaining_yaml_minus_beam_pct={delta-fraction:+.6f}; "
                     "sensitivity only, NOT evidence that BeamDyn omits this material")
    path = ROOT / "logs/check/mass/layer-report.txt"
    path.write_text("\n".join(lines)+"\n")
    for line in lines:
        if "MEAN" in line or "0.500000" in line or "YAML_INVALID_CHECK" in line:
            print(line)
    print("WROTE", path)


if __name__ == "__main__":
    main()
