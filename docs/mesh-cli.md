# `aeroelast mesh` — Mesh Generation CLI

Generate blade, rotor, tower and hub/nacelle meshes from a WindIO turbine
YAML and write STL / VTK / OBJ / MSH files. The command is a subcommand of the
main `aeroelast` entry point:

```bash
aeroelast mesh <component> [INPUT] [OPTIONS]
```

| Component | Description |
|-----------|-------------|
| `blade` | Single blade shell mesh (with or without shear webs) |
| `rotor` | Multi-blade rotor (blades rotated about the rotor axis) |
| `hub` | Hub + nacelle **single body**: hemispherical hub at the rotor + cylinder behind |
| `tower` | Tapered tower lofted from the WindIO outer-shape profile |
| `turbine` | Blades + hub/nacelle + tower, one file per component |

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

## Hub + nacelle (one body)

The hub and the nacelle are **the same geometry**: a hemispherical hub centred
on the rotor origin, then a **constant-radius cylinder** extending backwards
along the rotor axis, closed by a flat tail cap. The turbine therefore exports a
single `hub` mesh (one STL).

```bash
# From the WindIO definition (radius = components.hub.diameter / 2,
# length = nacelle overhang)
aeroelast mesh hub IEA-15-240-RWT.yaml --out hub.stl

# From parameters
aeroelast mesh hub --length 12 --diameter 7.94 --out hub.stl
aeroelast mesh hub --length 12 --radius 3.97 --axis 1,0,0 --n-circ 48 --out hub.stl
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

Writes `meshes/blade_1.stl … meshes/blade_N.stl`, `meshes/hub.stl` and
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
- The **nacelle is always horizontal** along the rotor axis, and the hub/nacelle
  body is a hemisphere centred on the rotor origin followed by a cylinder
  extending backwards (`+rotor_axis`).
- The blades are coned by the hub `cone_angle` (tips away from the tower) and
  distributed azimuthally about the rotor axis.

Example with an explicit frame:

```bash
aeroelast mesh turbine turbine.yaml --out-dir meshes/ \
    --rotor-axis 1,0,0 --tower-offset -12,0,-6 --quiet
```

## Options

Shared by the surface components (`hub`, `tower`):

| Flag | Description |
|------|-------------|
| `--n-circ N` | Circumferential elements (default 32) |
| `--no-caps` | Leave the open ends uncapped |

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
`--tower-n-axial`, `--tower-n-circ`, `--hub-n-circ`, `--hub-n-axial`,
`--hub-n-tip`, `--hub-diameter`, `--hub-length`, `--prefix` and `--renumber`.

Run `aeroelast mesh <component> --help` for the complete list.

## Python API

```python
from aeroelast.core.mesh import BladeMesh, HubNacelleMesh, TowerMesh, write_mesh
from aeroelast.core.mesh.turbine import TurbineMesh

blade = BladeMesh("turbine.yaml", element_size=0.1, include_webs=False).generate()
write_mesh(blade, "blade.stl")

tower = TowerMesh.from_windio("turbine.yaml").generate()
write_mesh(tower, "tower.obj")

# hub + nacelle: one body
body = HubNacelleMesh.from_windio("turbine.yaml", axis=(0, 1, 0)).generate()
write_mesh(body, "hub.stl")

result = TurbineMesh("turbine.yaml", element_size=0.5).generate()
result.write("meshes/", format="stl")  # {"blade_1": "...", "hub": "...", "tower": "..."}
result.rotor_axis   # (0, 1, 0)
result.tower_top    # tower-top position relative to the rotor centre
```
