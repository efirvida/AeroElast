# `aeroelast mesh` — Mesh Generation CLI

Generate blade, rotor, nacelle and tower meshes from a WindIO turbine YAML
and write STL / VTK / OBJ / MSH files. The command is a subcommand of the main
`aeroelast` entry point:

```bash
aeroelast mesh <component> [INPUT] [OPTIONS]
```

| Component | Description |
|-----------|-------------|
| `blade` | Single blade shell mesh (with or without shear webs) |
| `rotor` | Multi-blade rotor (blades rotated about the rotor axis) |
| `nacelle` | Nacelle **single body**: hemispherical hub cap + constant-radius cylinder + rear hemisphere |
| `tower` | Tapered tower lofted from the WindIO outer-shape profile |
| `turbine` | Blades + nacelle + tower, one file per component |

For a single-component command the output format is taken from the `--out`
extension; `--format` is the fallback when the path has no extension.
`turbine` writes one file per component into `--out-dir`.

## Blade without shear webs (CFD surface)

`--no-webs` omits the shear webs and returns only the blade outer-mold-line
surface — the wall boundary a CFD case needs. It is the default for
`aeroelast mesh turbine`; use `--webs` there to include the internal structure.

```bash
aeroelast mesh blade IEA-15-240-RWT.yaml --out blade.vtk --element-size 0.1
aeroelast mesh blade IEA-15-240-RWT.yaml --out blade.stl --no-webs
```

## Nacelle (hub and body are one mesh)

The hub cap and the nacelle body are **the same geometry**: a hemispherical hub
cap centred on the rotor origin, a **constant-radius cylinder** extending
backwards along the rotor axis, and a **rear hemisphere** at the tail so the
body stays streamlined instead of ending in a flat disc. The turbine therefore
exports a single `nacelle` mesh (one STL).

```bash
# From the WindIO definition (radius = components.hub.diameter / 2,
# length = nacelle overhang)
aeroelast mesh nacelle IEA-15-240-RWT.yaml --out nacelle.stl

# From parameters
aeroelast mesh nacelle --length 12 --diameter 7.94 --out nacelle.stl
aeroelast mesh nacelle --length 28 --radius 3.97 --axis 1,0,0 --n-circ 48 --out nacelle.stl

# Flat tail instead of the rear hemisphere
aeroelast mesh nacelle --length 12 --diameter 7.94 --flat-tail --out nacelle.stl
```

## Tower

The tower follows the imported `components.tower.outer_shape_bem`
(`reference_axis` + `outer_diameter`); the profile is resampled along `z` using
`--n-axial` or `--element-size`. Without an input file, pass parameters:

```bash
aeroelast mesh tower IEA-15-240-RWT.yaml --out tower.stl
aeroelast mesh tower --height 144 --base-diameter 10 --top-diameter 6.5 \
    --n-axial 30 --n-circ 48 --out tower.stl
```

## Full turbine

```bash
aeroelast mesh turbine IEA-15-240-RWT.yaml \
    --out-dir meshes/ --format stl --element-size 0.5
```

Writes `meshes/blade_1.stl … meshes/blade_N.stl`, `meshes/nacelle.stl` and
`meshes/tower.stl`. Coordinate convention:

- The **rotor is centred at the origin** `(0, 0, 0)`.
- The **rotor axis is configurable** with `--rotor-axis x,y,z`, default
  **`0,1,0`** (the same convention `RotorMesh` uses). The rotor plane is the
  `X-Z` plane by default, so the blade base mesh needs no remap.
- The **tower is the displaced component**: vertical along `+Z`, with its top at
  `--tower-offset x,y,z` relative to the rotor centre. When the YAML carries the
  data the default offset is `+overhang * rotor_axis` (plus an optional
  `-distance-tt-hub * Z` drop), i.e. the tower top sits **behind** the rotor
  plane (`+rotor_axis`) and rises to meet the horizontal nacelle. Override with
  `--tower-offset`, `--overhang`, `--distance-tt-hub` and/or `--tower-base-z`.
- The **nacelle is always horizontal** along the rotor axis. Its length defaults
  to `2 * overhang + radius`, so the body is **symmetric about the tower axis**
  and the tail clears the tower by more than the tower diameter.
- The blades are coned by the hub `cone_angle` (tips away from the tower) and
  distributed azimuthally about the rotor axis.

Example with an explicit frame:

```bash
aeroelast mesh turbine turbine.yaml --out-dir meshes/ \
    --rotor-axis 1,0,0 --tower-offset 12,0,0 --quiet
```

## Options

Shared by the surface components (`nacelle`, `tower`):

| Flag | Description |
|------|-------------|
| `--n-circ N` | Circumferential elements (default 32) |
| `--no-caps` | Leave the open ends uncapped |

Nacelle:

| Flag | Description |
|------|-------------|
| `--length L` | Axial extent from the hub centre to the rear tip |
| `--radius` / `--diameter` | Body radius (constant along the cylinder) |
| `--axis x,y,z` | Body axis (default `0,1,0`) |
| `--n-axial`, `--n-tip`, `--n-tail` | Cylinder, front-cap and rear-cap divisions |
| `--flat-tail` | Flat tail cap instead of the rear hemisphere |

Blade / rotor:

| Flag | Description |
|------|-------------|
| `--element-size E` | Target element size (default 0.1) |
| `--n-samples N` | Airfoil samples (default 300) |
| `--no-webs` | Surface only (CFD outer-mold-line) |
| `--span-grading` | `chord` (default) or `uniform` |
| `--renumber` | `simple`, `rcm`, or omitted |

`turbine` accepts `--n-blades`, `--element-size`, `--webs`, `--rotor-axis`,
`--tower-offset`, `--tower-base-z`, `--overhang`, `--distance-tt-hub`,
`--tower-n-axial`, `--tower-n-circ`, `--nacelle-n-circ`, `--nacelle-n-axial`,
`--nacelle-n-tip`, `--nacelle-n-tail`, `--nacelle-flat-tail`,
`--nacelle-diameter`, `--nacelle-length`, `--prefix` and `--renumber`.

Run `aeroelast mesh <component> --help` for the complete list.

## Python API

```python
from aeroelast.core.mesh import BladeMesh, NacelleMesh, TowerMesh, write_mesh
from aeroelast.core.mesh.turbine import TurbineMesh

blade = BladeMesh("turbine.yaml", element_size=0.1, include_webs=False).generate()
write_mesh(blade, "blade.stl")

tower = TowerMesh.from_windio("turbine.yaml").generate()
write_mesh(tower, "tower.obj")

# hub cap + nacelle + rear hemisphere: one mesh
body = NacelleMesh.from_windio("turbine.yaml", axis=(0, 1, 0)).generate()
write_mesh(body, "nacelle.stl")

result = TurbineMesh("turbine.yaml", element_size=0.5).generate()
result.write("meshes/", format="stl")  # {"blade_1": "...", "nacelle": "...", "tower": "..."}
result.rotor_axis      # (0, 1, 0)
result.tower_top       # tower-top position relative to the rotor centre
result.nacelle_radius  # blade rooting radius
```
