#!/usr/bin/env python3
"""CLI entry point for the BEM-FSI fluid participant.

Runs a BEM-based preCICE fluid participant as a standalone process,
equivalent to launching OpenFOAM on the fluid side.  The structural
solver (``aeroelast-fsi``) must be started in a separate process before
or concurrently with this command.

Usage
-----
::

    aeroelast-bem-fsi simulation.yaml
    aeroelast-bem-fsi simulation.yaml --workdir /path/to/case
    aeroelast-bem-fsi --template

Configuration file schema
-------------------------
The YAML file must contain at minimum::

    participant: "Fluid"            # preCICE participant name
    config_file: "precice-config.xml"
    coupling_mesh: "Fluid-Mesh"     # mesh name this participant provides
    blade_file: "path/to/blade.yaml"

    mesh:
      source: "generator"           # "generator" or "file"
      generator:
        type: "BladeMesh"
        params:
          yaml_file: "path/to/blade.yaml"
          element_size: 0.5

    bem:
      wind_speed: 10.59             # v_inf [m/s]
      omega: 7.56                   # angular velocity [rad/s]
      pitch: 0.0                    # collective pitch [deg]
      air_density: 1.225
      dynamic_viscosity: 1.81e-5
      hub_height: 150.0
      shear_exp: 0.2
      n_blades: 3
      hub_radius: 3.0
      default_re: 1.0e7
      neuralfoil_model: "large"
      precone: 0.0
      tilt: 0.0
      span_direction: [0.0, 0.0, 1.0]
      normal_direction: [1.0, 0.0, 0.0]
      tangential_direction: [0.0, 1.0, 0.0]

    output:
      folder: "bem_fsi_results"
      log_interval: 10              # windows between diagnostic lines

Notes
-----
The mesh section accepts any generator type supported by ``FSIRunner``
(``BladeMesh``, ``RotorMesh``, etc.) or ``source: file`` to load an
existing ``.h5`` mesh.  The user decides the mesh parameters;
the BEM-FSI participant does **not** impose any constraints.  preCICE
performs the spatial mapping between the fluid and solid meshes.
"""

import argparse
import logging
import os
import sys
from pathlib import Path
from typing import TYPE_CHECKING

import yaml

if TYPE_CHECKING:
    from aeroelast.core.mesh import MeshModel

_TEMPLATE = """\
# BEM-FSI fluid participant configuration
# ========================================
# Launch with: aeroelast-bem-fsi simulation.yaml

participant: "Fluid"
config_file: "precice-config.xml"
coupling_mesh: "Fluid-Mesh"
blade_file: "path/to/blade.yaml"

mesh:
  source: "generator"
  generator:
    type: "BladeMesh"
    params:
      yaml_file: "path/to/blade.yaml"
      element_size: 0.5
  # Optional: write mesh to file for inspection / preCICE support-radius tuning
  # output_file: "fluid_mesh.vtk"
    # Optional: export separate rotor blades when using RotorMesh
    # export_parts:
    #   enabled: true
    #   format: "stl"
    #   output_dir: "blade_parts"

bem:
  wind_speed: 10.59           # m/s
  omega: 7.56                 # rad/s (7.23 rpm for IEA-15)
  pitch: 0.0                  # deg
    azimuth: 0.0                # deg (blade position around rotor disc)
  air_density: 1.225          # kg/m³
  dynamic_viscosity: 1.81e-5  # Pa·s
  hub_height: 150.0           # m
  shear_exp: 0.2
  n_blades: 3
  hub_radius: 3.0             # m
  default_re: 1.0e7
  neuralfoil_model: "large"
  precone: 0.0                # deg
  tilt: 0.0                   # deg
  span_direction: [0.0, 0.0, 1.0]
  normal_direction: [1.0, 0.0, 0.0]
  tangential_direction: [0.0, 1.0, 0.0]
  # The wall-flow moment realisation is OFF by default: on the coupled gate case it moves the
  # tip section rotation by ~30 deg at an identical per-strip resultant, and that movement is
  # not attributed yet, so the strips keep the minimum-norm distribution unless a case asks.
  wall_flow_realisation: false

output:
  folder: "bem_fsi_results"
  log_interval: 10
"""


def _element_properties_for(mesh: "MeshModel", generator: object) -> dict | None:
    """Build the deck's section-property map for the coupling mesh, or ``None``.

    The wall-flow moment realisation resolves each coupling element through
    ``mesh.element_sets[set_name]``, so only the keys the projection mesh
    actually carries are useful here.  A key the mesh does not have is dead
    weight, and ``_RingSection.from_element_properties`` is set from
    ``element_properties is not None`` rather than from whether a wall resolved
    a property: a partially matching map would claim a physical split while
    some walls silently keep the geometric ``S = 1.0``.

    The map is therefore restricted to ``mesh.element_sets`` and returned only
    when the kept sets cover every coupling element.  When the generator has no
    deck data (a mesh loaded from file, or a generator without NuMAD data),
    when nothing matches, or when coverage is partial, the helper returns
    ``None`` so the caller keeps the minimum-norm ``_distribute`` fallback - the
    honest result beats a half-resolved split.

    Parameters
    ----------
    mesh : MeshModel
        The final (already filtered) coupling mesh.
    generator : object
        The generator instance that produced the mesh.

    Returns
    -------
    dict or None
        ``{element_set_name: property}`` restricted to the mesh's element sets,
        or ``None``.
    """
    numad_data = getattr(generator, "numad_mesh_data", None)
    if not numad_data:
        return None

    from aeroelast.models.blade.model import build_rust_properties  # noqa: PLC0415

    props = build_rust_properties(numad_data)
    props = {name: prop for name, prop in props.items() if name in mesh.element_sets}
    if not props:
        logging.warning(
            "[BEM-FSI] The deck has no section property matching the coupling mesh's "
            "element sets; withholding the map so strips use the minimum-norm distribution"
        )
        return None

    covered: set = set()
    for name in props:
        covered.update(mesh.element_sets[name].element_ids)
    if len(covered) != len(mesh.elements):
        logging.warning(
            "[BEM-FSI] %d of %d coupling elements have no section property; "
            "withholding the map so strips use the minimum-norm distribution",
            len(mesh.elements) - len(covered),
            len(mesh.elements),
        )
        return None
    return props


def _effective_element_properties(cfg: dict, element_properties: dict | None) -> dict | None:
    """The property map to hand to the projector: ``None`` unless the case asks for it.

    The wall-flow moment realisation is **off by default**. It is implemented, guarded and
    exercised by its own group, but on the coupled gate case it moved the tip section rotation
    from ``+2.583`` to ``-27.404`` deg at an identical per-strip resultant, with the flapwise
    deflection unchanged, and that movement is not attributed yet (``docs/validation/gaps.yaml``,
    ``moment_realization_over_delivers``). A case activates it with::

        bem:
          wall_flow_realisation: true
    """
    if element_properties is None:
        return None
    if not bool(cfg.get("bem", {}).get("wall_flow_realisation", False)):
        logging.info(
            "[BEM-FSI] Wall-flow moment realisation off (bem.wall_flow_realisation is not set); "
            "strips use the minimum-norm distribution"
        )
        return None
    logging.info("[BEM-FSI] Wall-flow moment realisation enabled by bem.wall_flow_realisation")
    return element_properties


def _build_mesh(cfg: dict, config_path: Path):
    """Build or load the coupling mesh from the YAML ``mesh`` section.

    Returns
    -------
    mesh : MeshModel
        The coupling mesh (node-set filtered when ``coupling_node_set`` is set).
    viz_mesh : MeshModel or None
        The full, unfiltered mesh kept for surface VTU output; ``None`` when no
        filter was applied.
    element_properties : dict or None
        The deck's section-property map restricted to the coupling mesh's
        element sets, when the generator carries deck data and the sets cover
        every coupling element; ``None`` otherwise (see
        :func:`_element_properties_for`).  This is the **available** map: what is
        actually handed to the projector is gated by ``bem.wall_flow_realisation``
        (off by default) through :func:`_effective_element_properties`.
    """
    from aeroelast.core.config import MeshGeneratorType, MeshSource
    from aeroelast.core.mesh import (
        BladeMesh,
        BoxSurfaceMesh,
        MeshModel,
        MultiFlapMesh,
        RotorMesh,
        SquareShapeMesh,
    )

    mesh_cfg = cfg.get("mesh", {})
    source = mesh_cfg.get("source", "generator")
    # Only the generator branch reads it; the node-set filter below asks for
    # ``coupling_node_set`` only when the source is not a file.
    # ``generator`` is bound only by the BladeMesh/RotorMesh generators; the
    # deck-data lookup below is skipped for every generator without one.
    gen_cfg: dict = {}
    generator = None

    def _resolve(path_str: str) -> str:
        p = Path(path_str)
        if not p.is_absolute():
            return str(config_path.parent / p)
        return path_str

    if source == MeshSource.FILE.value:
        file_cfg = mesh_cfg.get("file", {})
        file_path = _resolve(file_cfg["path"])
        fmt = file_cfg.get("format", "auto")
        logging.info("[BEM-FSI] Loading mesh from %s", file_path)
        mesh = MeshModel.load(file_path, format=fmt)

    else:
        gen_cfg = mesh_cfg.get("generator", {})
        gen_type = gen_cfg.get("type", MeshGeneratorType.BLADE.value)
        params = gen_cfg.get("params", {})

        logging.info("[BEM-FSI] Generating mesh with %s", gen_type)

        if gen_type == MeshGeneratorType.BLADE.value:
            yaml_file = params.get("yaml_file")
            excel_file = params.get("excel_file")
            airfoil_dir = params.get("airfoil_dir")
            if yaml_file:
                yaml_file = _resolve(yaml_file)
            if excel_file:
                excel_file = _resolve(excel_file)
            if airfoil_dir:
                airfoil_dir = _resolve(airfoil_dir)
            generator = BladeMesh(
                yaml_file=yaml_file,
                excel_file=excel_file,
                airfoil_dir=airfoil_dir,
                element_size=params.get("element_size", 0.5),
                n_samples=params.get("n_samples", 300),
                span_grading=params.get("span_grading", "chord"),
                airfoil_spacing=params.get("airfoil_spacing", "cosine"),
            )
            mesh = generator.generate(renumber=None)

        elif gen_type == MeshGeneratorType.ROTOR.value:
            yaml_file = params.get("yaml_file")
            excel_file = params.get("excel_file")
            airfoil_dir = params.get("airfoil_dir")
            hub_diameter = params.get("hub_diameter")
            rotor_diameter = params.get("rotor_diameter")
            if yaml_file:
                yaml_file = _resolve(yaml_file)
            if excel_file:
                excel_file = _resolve(excel_file)
            if airfoil_dir:
                airfoil_dir = _resolve(airfoil_dir)
            generator = RotorMesh(
                yaml_file=yaml_file,
                excel_file=excel_file,
                airfoil_dir=airfoil_dir,
                n_blades=params.get("n_blades", 3),
                hub_radius=params.get("hub_radius"),
                hub_diameter=hub_diameter,
                rotor_diameter=rotor_diameter,
                element_size=params.get("element_size", 0.5),
                n_samples=params.get("n_samples", 300),
                airfoil_spacing=params.get("airfoil_spacing", "cosine"),
            )
            mesh = generator.generate(renumber="rcm")

        elif gen_type == MeshGeneratorType.SQUARE.value:
            mesh = SquareShapeMesh(
                width=params["width"],
                height=params["height"],
                nx=params["nx"],
                ny=params["ny"],
                quadratic=params.get("quadratic", False),
                triangular=params.get("triangular", False),
            ).generate()

        elif gen_type == MeshGeneratorType.BOX.value:
            mesh = BoxSurfaceMesh(
                center=tuple(params["center"]),
                dims=tuple(params["dims"]),
                nx=params["nx"],
                ny=params["ny"],
                nz=params["nz"],
                quadratic=params.get("quadratic", False),
                triangular=params.get("triangular", False),
            ).generate()

        elif gen_type == MeshGeneratorType.MULTIFLAP.value:
            mesh = MultiFlapMesh(
                n_flaps=params["n_flaps"],
                flap_width=params["flap_width"],
                flap_height=params["flap_height"],
                x_spacing=params["x_spacing"],
                base_height=params.get("base_height", 0.05),
                nx_flap=params.get("nx_flap", 4),
                ny_flap=params.get("ny_flap", 20),
                nx_base_segment=params.get("nx_base_segment", 10),
                ny_base=params.get("ny_base", 2),
                quadratic=params.get("quadratic", False),
            ).generate()

        else:
            raise ValueError(f"Unknown mesh generator type: '{gen_type}'")

    # Optional: write full mesh for inspection (before any node-set filter)
    output_file = mesh_cfg.get("output_file")
    if output_file:
        mesh.write_mesh(_resolve(output_file))
        logging.info("[BEM-FSI] Mesh written to %s", output_file)

    export_parts_cfg = mesh_cfg.get("export_parts")
    if export_parts_cfg and export_parts_cfg.get("enabled", False):
        mesh_format = str(export_parts_cfg.get("format", "stl")).lstrip(".").lower()
        output_dir = _resolve(export_parts_cfg.get("output_dir", "mesh_parts"))
        set_names = sorted(
            name
            for name in mesh.element_sets_names
            if name.startswith(RotorMesh.BLADE_PART_SET_PREFIX)
        )
        if not set_names:
            raise ValueError(
                "mesh.export_parts is enabled, but no rotor blade element sets were found. "
                "This option currently supports RotorMesh blade exports."
            )
        written_files = mesh.write_element_sets(output_dir, set_names, mesh_format)
        logging.info(
            "[BEM-FSI] Wrote %d mesh parts to %s",
            len(written_files),
            output_dir,
        )

    # Optional: filter to a specific node set (e.g. "allOuterShellNods" to
    # exclude shear-web nodes from the aerodynamic coupling mesh).
    # The BEM participant needs node coordinates for preCICE registration
    # AND the fully contained elements, because the wall-flow moment
    # realisation walks the section's chordwise wall graph; a nodes-only
    # coupling mesh makes every ring section unrealisable and silently falls
    # back to the minimum-norm distribution (#16).  The node list - and so the
    # preCICE vertex order - is exactly the node set's, unchanged.
    # The full mesh (with elements) is kept as viz_mesh for VTU surface output.
    viz_mesh = None
    coupling_node_set = (
        gen_cfg.get("coupling_node_set") if source != MeshSource.FILE.value else None
    )
    if coupling_node_set:
        try:
            ns = mesh.get_node_set(coupling_node_set)
        except ValueError:
            available = list(mesh.node_sets.keys())
            raise ValueError(
                f"coupling_node_set '{coupling_node_set}' not found. "
                f"Available node sets: {available}"
            )
        viz_mesh = mesh  # full mesh with elements, for surface VTU
        mesh = mesh.subset_to_nodes(list(ns.nodes.keys()))
        logging.info(
            "[BEM-FSI] Filtered to node set '%s': %d nodes, %d elements",
            coupling_node_set,
            len(mesh.nodes),
            len(mesh.elements),
        )

    # The deck's section-property map is built from the very generator that
    # produced the mesh, after the node-set filter, so coverage is measured on
    # the filtered mesh and the keys match its element sets.  A file-sourced
    # mesh has no generator and keeps the minimum-norm distribution.
    element_properties = None
    if generator is not None:
        element_properties = _element_properties_for(mesh, generator)

    logging.info(
        "[BEM-FSI] Mesh ready: %d nodes, %d elements",
        len(mesh.nodes),
        len(mesh.elements),
    )
    return mesh, viz_mesh, element_properties


def main(argv=None) -> int:
    """CLI entry point."""
    parser = argparse.ArgumentParser(
        prog="aeroelast-bem-fsi",
        description=(
            "BEM-based preCICE fluid participant for FSI coupling. "
            "Runs as an independent process alongside the structural solver."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "examples:\n"
            "  aeroelast-bem-fsi simulation.yaml\n"
            "  aeroelast-bem-fsi simulation.yaml --workdir /path/to/case\n"
            "  aeroelast-bem-fsi --template > simulation.yaml\n"
        ),
    )
    parser.add_argument(
        "config",
        nargs="?",
        metavar="CONFIG_FILE",
        help="Path to the YAML configuration file.",
    )
    parser.add_argument(
        "--workdir",
        metavar="DIR",
        help="Working directory (default: directory of CONFIG_FILE).",
    )
    parser.add_argument(
        "--template",
        action="store_true",
        help="Print a template YAML configuration and exit.",
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        help="Logging verbosity (default: INFO).",
    )

    args = parser.parse_args(argv)

    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%H:%M:%S",
    )

    if args.template:
        print(_TEMPLATE, end="")
        return 0

    if args.config is None:
        parser.error("CONFIG_FILE is required (or use --template to generate one).")

    config_path = Path(args.config).resolve()
    if not config_path.exists():
        logging.error("Configuration file not found: %s", config_path)
        return 1

    # Change working directory
    workdir = Path(args.workdir).resolve() if args.workdir else config_path.parent
    logging.info("[BEM-FSI] Working directory: %s", workdir)

    # Load YAML
    with open(config_path) as f:
        cfg = yaml.safe_load(f)

    # Resolve file paths that are relative to the config file's directory BEFORE
    # changing CWD.  The working directory (workdir) controls where preCICE sockets
    # and output files land; it must NOT affect how config-referenced data files are
    # found.  All relative paths in the YAML are anchored at config_path.parent.
    yaml_dir = config_path.parent

    def _resolve(val: str) -> str:
        p = Path(val)
        return str((yaml_dir / p).resolve()) if not p.is_absolute() else val

    if "config_file" in cfg:
        cfg["config_file"] = _resolve(cfg["config_file"])
    if "blade_file" in cfg:
        cfg["blade_file"] = _resolve(cfg["blade_file"])
    mesh_gen = cfg.get("mesh", {}).get("generator", {}).get("params", {})
    if "excel_file" in mesh_gen:
        mesh_gen["excel_file"] = _resolve(mesh_gen["excel_file"])
    if "airfoil_dir" in mesh_gen:
        mesh_gen["airfoil_dir"] = _resolve(mesh_gen["airfoil_dir"])

    os.chdir(workdir)

    # Build the mesh (user's responsibility to configure)
    try:
        mesh, viz_mesh, element_properties = _build_mesh(cfg, config_path)
    except Exception as exc:
        logging.error("Failed to build mesh: %s", exc)
        return 1

    # Build and run the BEM-FSI participant
    from aeroelast.solvers.bem.fsi_participant import build_from_config

    try:
        participant = build_from_config(
            mesh,
            cfg,
            viz_mesh=viz_mesh,
            element_properties=_effective_element_properties(cfg, element_properties),
        )
        participant.run()
    except Exception as exc:
        logging.exception("BEM-FSI participant failed: %s", exc)
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
