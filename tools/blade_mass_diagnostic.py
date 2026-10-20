"""Localise YAML-to-shell mass differences without static stiffness solves.

Run through SLURM: scripts/aeroenv.sh python -B tools/blade_mass_diagnostic.py
Outputs logs/check/mass/diagnostic.json. OML remains the production default.
Uses the native S-1 station mass path, not a substitute mass definition.
Layer masses use the actual mesh areas and the component identity retained
by NuMAD's stack database. Width is area/span (includes taper and curvature).
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))

from beam_reference import load_beamdyn_blade  # noqa: E402
from yaml_layup_independent import interp, load, section  # noqa: E402
from aeroelast.core.mesh.generators import BladeMesh  # noqa: E402
from aeroelast.models.blade.model import build_rust_properties  # noqa: E402
from aeroelast.postprocess.sectional import SectionalExtractor  # noqa: E402


def main():
    data = load()
    layers = data["components"]["blade"]["internal_structure_2d_fem"]["layers"]
    densities = {m["name"]: float(m["rho"]) for m in data["materials"]}
    bd = load_beamdyn_blade(ROOT / "tests/IEA15MW/reference/IEA-15-240-RWT_BeamDyn_blade.dat")
    output = {}
    for reference in ("oml", "mid"):
        gen = BladeMesh(yaml_file=str(ROOT / "tests/IEA-15-240-RWT.yaml"),
                        element_size=0.25, skin_reference=reference)
        mesh = gen.generate(renumber="rcm")
        ext = SectionalExtractor(mesh, build_rust_properties(gen.numad_mesh_data))
        # Preserve S-1's station acceptance as well as its mass formula:
        # extract() omits stations whose experimental stiffness solve fails.
        native = ext.extract()
        z, mass = native.z, native.mass_dens
        mid = (z > 0.2*z.max()) & (z < 0.85*z.max())
        error = mass[mid]/np.interp(z[mid], bd.z, bd.mass_dens)-1
        print(f"S1_MASS {reference} total_kg={ext.node_mass.sum():.6f} "
              f"signed_pct={100*error.mean():.6f} absolute_pct={100*abs(error).mean():.6f}",
              flush=True)
        blade = gen.numad_blade
        stacks = {s.name: s for s in blade.stackdb.stacks.flat}
        stacks.update({s.name: s for s in blade.stackdb.swstacks.flat})
        by_station = {}
        xyz = np.asarray(mesh.coords_array)
        idx = mesh.node_id_to_index
        total_area_mass = 0.0
        for name, stack in stacks.items():
            if name not in mesh.element_sets:
                continue
            station = stack.indices[0]
            dz = float(blade.ispan[station+1]-blade.ispan[station])
            area = 0.0
            for e in mesh.element_sets[name].elements:
                p = xyz[[idx[n] for n in e.node_ids]]
                area += (0.5*np.linalg.norm(np.cross(p[2]-p[0], p[3]-p[1]))
                         if len(p) == 4 else
                         0.5*np.linalg.norm(np.cross(p[1]-p[0], p[2]-p[0])))
            row = by_station.setdefault(station, {})
            for pg in stack.plygroups:
                h = .001*pg.thickness*pg.nPlies
                contribution = densities[pg.materialid]*h*area
                total_area_mass += contribution
                entry = row.setdefault(pg.component, {"mass": 0., "width_m": 0.,
                                                      "thicknesses_m": [], "regions": []})
                entry["mass"] += contribution/dz
                entry["width_m"] += area/dz
                entry["thicknesses_m"].append(h)
                entry["regions"].append(name)
        rows = []
        for station, mapped in sorted(by_station.items()):
            eta = float((blade.ispan[station]+blade.ispan[station+1])/234.)
            try:
                independent = section(data, eta)
            except ValueError as exc:
                print(f"YAML_INVALID eta={eta:.9f}: {exc}", flush=True)
                continue
            for layer in independent["layers"]:
                spec = next(l for l in layers if l["name"] == layer["name"])
                layer["density_kg_m3"] = densities[layer["material"]]
                layer["yaml_nd"] = {
                    key: interp(spec[key], eta) for key in ("start_nd_arc", "end_nd_arc")
                    if key in spec and "grid" in spec[key]
                }
                layer["mapped"] = mapped.get(layer["name"], {})
                layer["delta_kg_m"] = layer["mapped"].get("mass", 0.)-layer["mass"]
            rows.append(independent)
        print(f"AREA_MASS {reference} kg={total_area_mass:.6f}", flush=True)
        common = []
        for eta in np.linspace(.2, .85, 12):
            independent = section(data, eta)
            differences = {}
            for layer in independent["layers"]:
                name = layer["name"]
                mapped_mass = np.interp(eta, [r["eta"] for r in rows],
                                        [r["layers"][next(i for i, l in enumerate(r["layers"])
                                                          if l["name"] == name)]
                                         ["mapped"].get("mass", 0.) for r in rows])
                differences[name] = float(mapped_mass-layer["mass"])
            common.append({"eta": float(eta), "differences_kg_m": differences})
        for layer in layers:
            if layer["material"] == "glass_triax":
                name = layer["name"]
                delta = np.mean([r["differences_kg_m"][name] for r in common])
                print(f"LAYER {reference} {name} mean_delta_kg_m={delta:+.6f}", flush=True)
        output[reference] = {"total_kg": float(ext.node_mass.sum()),
                             "area_mass_kg": total_area_mass,
                             "s1_signed_pct": float(100*error.mean()),
                             "s1_absolute_pct": float(100*abs(error).mean()),
                             "stations": rows, "common": common}
    path = ROOT / "logs/check/mass/diagnostic.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(output, indent=2)+"\n")
    print("WROTE", path, flush=True)


if __name__ == "__main__":
    main()
