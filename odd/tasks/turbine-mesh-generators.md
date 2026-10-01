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

- [x] T1 `include_webs` flag (WU1).
- [x] T2 Axisymmetric builder + tower/hub/nacelle (WU2).
- [x] T3 WindIO component reader (WU2).
- [x] T4 `TurbineMesh` assembly (WU3).
- [x] T5 `aeroelast mesh` CLI (WU4).
- [x] T6 Docs + verification (WU5).

## Result

Implemented and committed as six work units:

| Commit | Work unit | Files |
|---|---|---|
| `5320032` | `include_webs` for blade/rotor | `mesh_gen/mesh_gen.py`, `generators.py`, `model.py`, `tests/test_blade_no_webs.py` |
| `95a9aea` | tower/hub/nacelle + WindIO reader | `core/mesh/components.py`, `core/mesh/__init__.py`, `tests/test_axisymmetric_components.py` |
| `97fe652` | nacelle length = overhang + typing | `components.py`, test |
| `0ebbf62` | `TurbineMesh` assembly | `core/mesh/turbine.py`, `__init__.py`, `tests/test_turbine_mesh.py` |
| WU4 | `aeroelast mesh` CLI | `cli/mesh.py`, `cli/aeroelast.py`, `tests/test_cli_mesh.py` |
| WU5 | docs | `docs/mesh-cli.md`, `README.md`, `docs/cli-reference.md` |

### Design notes

1. **No webs** is a flag, not a new mesher: `get_shell_mesh` collapses its
   shear-web loop range to `range(rws - 1 if include_webs else 0)`, so the
   outer-shell code path is byte-identical. The web element count is exactly
   the element-count difference between the two modes (asserted).
2. **Components are analytic**, not gmsh-meshed: a meridian profile is
   revolved by `build_revolved_shell`, giving exact, deterministic element
   counts and correct outward orientation (verified by edge-manifold and
   signed-volume checks). gmsh remains available for `.msh`; meshio covers
   STL/VTK/OBJ through the existing `write_mesh` dispatch.
3. **Nacelle length is `overhang`** (yaw axis -> hub), not
   `distance_tt_hub` (tower-top -> hub vertical distance); corrected during
   WU2 before the turbine assembly relied on it.
4. **Turbine frame**: tower base at `tower_base_z` (default 0) via a rigid
   translation, hub at `hub_height + base_offset`, rotor axis `+X` tilted by
   the nacelle uptilt, blades coned by `cone_angle` and azimuthally
   distributed. The base-blade frame (span +Z, chord +X, thickness +Y) is
   remapped with `Rz(-90 deg)` to (span +Z, rotor axis +X, chord -Y).

### Verification evidence

| Suite | Result |
|---|---|
| `tests/test_blade_no_webs.py` | 6 passed |
| `tests/test_axisymmetric_components.py` | 15 passed |
| `tests/test_turbine_mesh.py` | 10 passed |
| `tests/test_cli_mesh.py` | 13 passed |
| `ruff check` + `ruff format --check` (changed files) | clean |
| `pyright` (new/changed files, standalone) | 0 errors, 0 warnings |
| `tests/test_blade_mesh.py` (existing regression) | 1 passed |

### Known environment issue

The lens LSP reports `reportMissingImports` for the three files created during
this session (`core/mesh/components.py`, `core/mesh/turbine.py`,
`cli/mesh.py`) while the standalone `pyright` binary — same engine, same
`pyrightconfig.json`, same workspace root — reports **0 errors** and the
imports resolve at runtime. This is a stale language-server workspace index,
not a code defect.

---

## Corrections — iteration 2 (user feedback)

The first implementation was rejected on three points and reworked:

1. **Rotor axis.** The first pass remapped the blade with `Rz(-90 deg)` to spin
   about **+X**, breaking the `RotorMesh` convention. Now the **rotor axis is
   configurable** (`rotor_axis`, CLI `--rotor-axis x,y,z`) and defaults to
   **+Y**, with **no remap** for the default axis (`_rotation_between(+Y, axis)`
   handles any other axis).
2. **Tower side.** The rotor is now **centred at the origin** and the **tower is
   the displaced component**: vertical `+Z`, top at `tower_offset` (default
   `-overhang * rotor_axis - distance_tt_hub * Z`), so the tower sits **behind**
   the rotor plane. `tower_offset` / `overhang` / `distance_tt_hub` /
   `tower_base_z` are all parameters, so a file with no nacelle data can still
   be placed.
3. **Hub = nacelle.** `HubMesh` (spheroid) and `NacelleMesh` (tapered capsule)
   were removed and replaced by a single **`HubNacelleMesh`**: a
   **constant-radius cylinder + hemispherical tip** (flat tail cap). The
   turbine now exports **one `hub` mesh / one STL**, and the `nacelle` CLI
   subcommand is gone. The radius defaults to `components.hub.diameter / 2`
   and the length to the nacelle overhang.

Also flipped the **coning sign** so the blade tips tilt **away** from the
tower (an upwind turbine cones away from the tower).

### Iteration-2 verification

| Suite | Result |
|---|---|
| `tests/test_axisymmetric_components.py` | 15 passed |
| `tests/test_turbine_mesh.py` | 13 passed |
| `tests/test_cli_mesh.py` | 14 passed |
| `tests/test_blade_no_webs.py` | 6 passed |
| `pyright` (new/changed files, standalone) | 0 errors, 0 warnings |
| `ruff check` | clean |

Geometry is asserted, not assumed: the new tests prove the rotation axis is `+Y`
by checking that a rotation about it preserves every blade's `y` profile, that
the tower plane lies entirely at negative `y` (behind the rotor), that the body
radius is constant up to the round tip, and that the body volume matches the
analytic cylinder+hemisphere. The base blade's **prebend/sweep** (the tip is not
at `y=0`) is what the first iteration's plane test missed.
