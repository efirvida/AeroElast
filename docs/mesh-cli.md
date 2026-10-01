# `aeroelast mesh` — Mesh Generation CLI

Generate blade, rotor, hub, nacelle and tower meshes from a WindIO turbine
YAML and write STL / VTK / OBJ / MSH files. The command is a subcommand of the
main `aeroelast` entry point:

```bash
aeroelast mesh <component> [INPUT] [OPTIONS]
```

| Component | Description |
|-----------|-------------|
| `blade` | Single blade shell mesh (with or without shear webs) |
| `rotor` | Multi-blade rotor (blades rotated about the rotor axis) |
| `hub` | Hub spheroid |
| `nacelle` | Nacelle capsule (body + tapered nose + flat tail) |
| `tower` | Tapered tower lofted from the WindIO outer-shape profile |
| `turbine` | Blades + hub + nacelle + tower, one file per component |

For a single-component command the output format is taken from the `--out`
extension; `--format` is the fallback when the path has no extension.
`turbine` writes one file per component into `--out-dir`.

## Blade without shear webs (CFD surface)

`--no-webs` omits the shear webs and returns only the blade outer-mold-line
surface — the wall boundary a CFD case needs. It is the default for
`aeroelast mesh turbine`; use `--webs` there to include the internal structure.

```bash
# Structural shell mesh (default: includes shear webs)
aeroelast mesh blade IEA-15-240-RWT.yaml --out blade.vtk --element-size 0.1

# CFD wall surface
aeroelast mesh blade IEA-15-240-RWT.yaml --out blade.stl --no-webs
```

## Hub, nacelle and tower

Without an input file, or when the input file lacks the block, pass explicit
parameters:

```bash
# From the WindIO definition
aeroelast mesh tower IEA-15-240-RWT.yaml --out tower.stl
aeroelast mesh hub   IEA-15-240-RWT.yaml --out hub.stl
aeroelast mesh nacelle IEA-15-240-RWT.yaml --out nacelle.obj

# From parameters (no WindIO tower/hub/nacelle block needed)
aeroelast mesh tower --height 144 --base-diameter 10 --top-diameter 6.5 \
    --n-axial 30 --n-circ 48 --out tower.stl
aeroelast mesh hub --diameter 7.94 --n-merid 24 --n-circ 48 --out hub.stl
aeroelast mesh nacelle --length 12 --body-diameter 3 --nose-diameter 2.2 \
    --out nacelle.stl
```

The tower follows the imported `components.tower.outer_shape_bem`
(`reference_axis` + `outer_diameter`); tower `x/y` centreline offsets are kept
and the profile is resampled along `z` using `--n-axial` or
`--element-size`. Hub uses `components.hub.diameter`; nacelle uses
`components.nacelle.drivetrain.lss_diameter` / `nose_diameter` and `overhang`
as the body length.

## Full turbine

```bash
aeroelast mesh turbine IEA-15-240-RWT.yaml \
    --out-dir meshes/ \
    --format stl \
    --element-size 0.5
```

Writes `meshes/blade_1.stl … meshes/blade_N.stl`, `hub.stl`, `nacelle.stl` and
`tower.stl` in one shared frame:

- tower base at `z = --tower-base-z` (default `0`), tower axis `+Z`;
- hub centre at `(overhang, 0, hub_height)` with the imported `hub_height`;
- rotor axis `+X` tilted by the nacelle `uptilt`;
- blades coned by the hub `cone_angle` and distributed azimuthally.

## Options

Shared by the surface components (`hub`, `nacelle`, `tower`):

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

`turbine` accepts `--n-blades`, `--element-size`, `--webs`,
`--tower-base-z`, `--tower-n-axial`, `--tower-n-circ`, `--hub-n-circ`,
`--hub-n-merid`, `--nacelle-n-circ`, `--prefix` and `--renumber`.

Run `aeroelast mesh <component> --help` for the complete list.

## Python API

The CLI is a thin wrapper; the same generators are importable:

```python
from aeroelast.core.mesh import BladeMesh, HubMesh, NacelleMesh, TowerMesh, write_mesh
from aeroelast.core.mesh.turbine import TurbineMesh

blade = BladeMesh("turbine.yaml", element_size=0.1, include_webs=False).generate()
write_mesh(blade, "blade.stl")

tower = TowerMesh.from_windio("turbine.yaml").generate()
write_mesh(tower, "tower.obj")

result = TurbineMesh("turbine.yaml", element_size=0.5).generate()
result.write("meshes/", format="stl")  # {"blade_1": "...", "tower": "...", ...}
```
