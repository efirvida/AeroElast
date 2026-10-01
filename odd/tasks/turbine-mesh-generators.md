# Feature: turbine mesh generators (surface-only blade, hub/tower/nacelle, TurbineMesh, CLI)

## Objective

Extend the existing blade/rotor mesher so that it can:

1. Emit a blade **without shear webs** — the outer-mold-line surface used directly
   as a CFD wall.
2. Generate **hub, nacelle and tower** meshes from the imported WindIO
   definition, or from explicit parameters when the input file lacks them.
   Exportable to STL / VTK / OBJ (and `.msh`).
3. Provide a **`TurbineMesh`** generator that produces the blade(s), hub/nacelle
   and tower as **separate** meshes in a single turbine coordinate system.
4. Expose all of it through the CLI as `aeroelast mesh <component>`.

## Decisions taken (user)

| Topic | Decision |
|---|---|
| No-web blade | Flag `include_webs=False` on `BladeMesh` / `get_shell_mesh` (same generator) |
| Hub/nacelle/tower fidelity | Fully parametric: tower = lofted tapered cylinder from `outer_diameter` + `reference_axis`; hub = sphere/cylinder from `diameter` + `cone_angle`; nacelle = capsule from `lss_diameter` / `nose_diameter` / `overhang` |
| Turbine output | Separate `MeshModel` per component written to a directory; rotor axis **+X**, tower base at `z=0`, hub at `hub_height`, tilt/cone from WindIO |
| CLI | `aeroelast mesh blade|rotor|hub|tower|turbine ...` |

## Coordinate convention (turbine frame)

- Tower base at `(0, 0, 0)`, tower axis **+Z** (reference-axis x/y offsets kept).
- Hub centre at `(0, 0, hub_height)`, rotor axis **+X**, tilted by
  `nacelle.drivetrain.uptilt` about **Y** and shifted by `overhang` along the
  tilted axis.
- Blade span is the base-blade **+Z** direction; each blade is coned about the
  rotor-plane tangential axis by `hub.cone_angle`, azimuthally rotated about the
  rotor axis by `2*pi*i/n_blades`, then translated to the hub.
- Nacelle axis **+X**, its tail at the tower top and nose toward the hub.

## Work units

- **WU1** — `include_webs` flag through `get_shell_mesh`, `BladeMesh`,
  `RotorMesh`, `Blade`, `Rotor`.
- **WU2** — axisymmetric component builder + `TowerMesh`, `HubMesh`,
  `NacelleMesh`, WindIO component reader.
- **WU3** — `TurbineMesh` assembly and per-component export.
- **WU4** — `aeroelast mesh` subcommands.
- **WU5** — docs (`docs/cli-reference.md`, `README.md`).

## Acceptance criteria

- [ ] `get_shell_mesh(..., include_webs=False)` produces the outer shell with no
      shear-web element/node sets and fewer elements than the with-webs mesh.
- [ ] `BladeMesh`/`RotorMesh`/`Blade`/`Rotor` accept and forward `include_webs`.
- [ ] `TowerMesh` from the IEA-15-240-RWT YAML reproduces the outer diameter
      profile and reference axis (mean radius/height checked against the YAML).
- [ ] `HubMesh` / `NacelleMesh` produce closed surfaces whose bounding box
      matches the requested radius/length.
- [ ] `TurbineMesh` returns separate meshes for every blade, hub/nacelle and
      tower, in one consistent frame.
- [ ] `aeroelast mesh ...` writes STL/VTK/OBJ for each component.
- [ ] Focused pytest for each work unit is green; existing suite not degraded.

## Verification command

```bash
source ~/miniconda3/etc/profile.d/conda.sh && conda activate aeroelast-dev &&
python -m pytest tests/test_blade_mesh.py tests/test_blade_no_webs.py \
  tests/test_axisymmetric_components.py tests/test_turbine_mesh.py \
  tests/test_cli_mesh.py -q
```

## Tasks

- [ ] T1 `include_webs` flag (WU1).
- [ ] T2 Axisymmetric builder + tower/hub/nacelle (WU2).
- [ ] T3 WindIO component reader (WU2).
- [ ] T4 `TurbineMesh` assembly (WU3).
- [ ] T5 `aeroelast mesh` CLI (WU4).
- [ ] T6 Docs + verification (WU5).

## Result

Pending.
