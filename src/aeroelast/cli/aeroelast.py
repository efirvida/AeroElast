#!/usr/bin/env python3
"""
aeroelast — Structural FEM solver CLI.

Runs modal, static and dynamic FEM analyses from YAML configuration files.
Assembly and element kernels execute in Rust; mesh generation and
post-processing remain in Python.

Usage:
    aeroelast model_config.yaml
    aeroelast model_config.yaml --workdir /path/to/case
    aeroelast model_config.yaml --preview
    aeroelast --template > config.yaml

Supported solver types (solver.type in YAML):
    Modal           — natural frequencies and mode shapes (SLEPc)
    LinearStatic    — static linear analysis
    LinearDynamic   — transient Newmark-β integration
    AeroFSI         — standalone aerodynamic participant for preCICE coupling
"""

import argparse
import logging
import sys
from pathlib import Path


def setup_logging(verbose: bool = False) -> None:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%H:%M:%S",
    )


def _resolve_cli_path(path_str: str | None, base_dir: Path) -> str | None:
    """Resolve CLI paths relative to the active working directory."""
    if path_str is None:
        return None
    path = Path(path_str)
    if not path.is_absolute():
        path = base_dir / path
    return str(path)


def _build_mesh_from_cli(args: argparse.Namespace, base_dir: Path):
    """Build a BladeMesh or RotorMesh directly from CLI arguments."""
    from aeroelast.core.mesh import BladeMesh, RotorHubMesh, RotorMesh

    if args.mesh_generator is None:
        raise ValueError("--mesh-generator is required with --export-mesh or --export-parts-dir")

    if args.mesh_generator == "BladeMesh":
        if not args.blade_yaml and not args.excel_file:
            raise ValueError("BladeMesh export requires --blade-yaml or --excel-file")
        generator = BladeMesh(
            yaml_file=_resolve_cli_path(args.blade_yaml, base_dir),
            excel_file=_resolve_cli_path(args.excel_file, base_dir),
            airfoil_dir=_resolve_cli_path(args.airfoil_dir, base_dir),
            element_size=args.element_size,
            n_samples=args.n_samples,
            span_grading=args.span_grading,
            airfoil_spacing=args.airfoil_spacing,
        )
    elif args.mesh_generator == "RotorMesh":
        if not args.blade_yaml and not args.excel_file:
            raise ValueError("RotorMesh export requires --blade-yaml or --excel-file")
        generator = RotorMesh(
            yaml_file=_resolve_cli_path(args.blade_yaml, base_dir),
            excel_file=_resolve_cli_path(args.excel_file, base_dir),
            airfoil_dir=_resolve_cli_path(args.airfoil_dir, base_dir),
            n_blades=args.n_blades,
            hub_radius=args.hub_radius,
            hub_diameter=args.hub_diameter,
            rotor_diameter=args.rotor_diameter,
            element_size=args.element_size,
            n_samples=args.n_samples,
            airfoil_spacing=args.airfoil_spacing,
        )
    elif args.mesh_generator == "RotorHubMesh":
        if not args.blade_yaml and not args.excel_file:
            raise ValueError("RotorHubMesh export requires --blade-yaml or --excel-file")
        generator = RotorHubMesh(
            yaml_file=_resolve_cli_path(args.blade_yaml, base_dir),
            excel_file=_resolve_cli_path(args.excel_file, base_dir),
            airfoil_dir=_resolve_cli_path(args.airfoil_dir, base_dir),
            n_blades=args.n_blades,
            hub_radius=args.hub_radius,
            hub_diameter=args.hub_diameter,
            rotor_diameter=args.rotor_diameter,
            element_size=args.element_size,
            n_samples=args.n_samples,
            airfoil_spacing=args.airfoil_spacing,
            hub_length=args.hub_length,
            connector_radius=args.hub_connector_radius,
            nose_radius=args.hub_nose_radius,
        )
    else:
        raise ValueError(f"Unsupported mesh generator for CLI export: {args.mesh_generator}")

    return generator.generate(renumber=args.renumber, verbose=args.verbose)


def _export_mesh_from_cli(args: argparse.Namespace) -> int:
    """Export a mesh directly from CLI arguments and exit."""
    from aeroelast.core.mesh import RotorMesh

    base_dir = Path(args.workdir) if args.workdir else Path.cwd()
    base_dir.mkdir(parents=True, exist_ok=True)

    mesh = _build_mesh_from_cli(args, base_dir)

    if args.export_mesh:
        output_path = Path(args.export_mesh)
        if not output_path.is_absolute():
            output_path = base_dir / output_path
        output_path.parent.mkdir(parents=True, exist_ok=True)
        mesh.write_mesh(str(output_path))
        print(f"Mesh written to: {output_path}")

    if args.export_parts_dir:
        if args.mesh_generator != "RotorMesh":
            raise ValueError(
                "--export-parts-dir currently supports only --mesh-generator RotorMesh"
            )
        output_dir = Path(args.export_parts_dir)
        if not output_dir.is_absolute():
            output_dir = base_dir / output_dir
        set_names = sorted(
            name
            for name in mesh.element_sets_names
            if name.startswith(RotorMesh.BLADE_PART_SET_PREFIX)
        )
        if not set_names:
            raise ValueError("No rotor blade element sets were found for separate export")
        written_files = mesh.write_element_sets(
            output_dir,
            set_names,
            args.export_parts_format,
        )
        print(f"Wrote {len(written_files)} mesh parts to: {output_dir}")

    if not args.export_mesh and not args.export_parts_dir:
        raise ValueError("Specify --export-mesh and/or --export-parts-dir")

    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="aeroelast",
        description="Structural FEM solver (Rust-accelerated assembly).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  aeroelast modal.yaml                  Run modal analysis
  aeroelast modal.yaml --workdir /tmp   Set working directory
  aeroelast modal.yaml --preview        Preview config without running
  aeroelast modal.yaml --validate       Validate config and exit
  aeroelast --template > config.yaml    Generate template config
  aeroelast --export-mesh blade.vtk --mesh-generator BladeMesh --blade-yaml IEA.yaml
    aeroelast --export-mesh rotor.stl --mesh-generator RotorMesh --blade-yaml IEA.yaml --n-blades 3
        aeroelast --export-mesh hub.stl --mesh-generator RotorHubMesh --blade-yaml IEA.yaml --n-blades 3
    aeroelast --export-parts-dir blades --mesh-generator RotorMesh --excel-file blade.xlsx --airfoil-dir airfoils --rotor-diameter 242.23775645 --n-blades 3
        """,
    )

    parser.add_argument(
        "config",
        nargs="?",
        help="Path to YAML configuration file",
    )
    parser.add_argument(
        "--workdir",
        "-w",
        help="Working directory (default: directory of the config file)",
    )
    parser.add_argument(
        "--preview",
        "-p",
        action="store_true",
        help="Print parsed configuration and exit",
    )
    parser.add_argument(
        "--validate",
        action="store_true",
        help="Validate configuration file and exit",
    )
    parser.add_argument(
        "--template",
        "-t",
        action="store_true",
        help="Print a template YAML configuration to stdout and exit",
    )
    parser.add_argument(
        "--verbose",
        "-v",
        action="store_true",
        help="Enable verbose logging",
    )
    parser.add_argument(
        "--export-ccx",
        metavar="FILE",
        help="Export model to CalculiX .inp format and exit",
    )
    parser.add_argument(
        "--export-mesh",
        metavar="FILE",
        help="Generate and export a BladeMesh or RotorMesh directly from CLI options and exit",
    )
    parser.add_argument(
        "--export-parts-dir",
        metavar="DIR",
        help="Export rotor blades as separate mesh files in DIR and exit",
    )
    parser.add_argument(
        "--export-parts-format",
        default="stl",
        help="Format used with --export-parts-dir (default: stl)",
    )
    parser.add_argument(
        "--mesh-generator",
        choices=("BladeMesh", "RotorMesh", "RotorHubMesh"),
        help="Mesh generator to use with --export-mesh or --export-parts-dir",
    )
    parser.add_argument(
        "--blade-yaml",
        help="WindIO/NuMAD YAML blade definition used by direct mesh export",
    )
    parser.add_argument(
        "--excel-file",
        help="NuMAD Excel blade definition used by direct BladeMesh or RotorMesh export",
    )
    parser.add_argument(
        "--airfoil-dir",
        help="Airfoil directory used with --excel-file for direct BladeMesh or RotorMesh export",
    )
    parser.add_argument(
        "--n-blades",
        type=int,
        default=3,
        help="Number of blades used by RotorMesh export (default: 3)",
    )
    parser.add_argument(
        "--hub-radius",
        type=float,
        help="Optional hub radius override for RotorMesh export",
    )
    parser.add_argument(
        "--hub-diameter",
        type=float,
        help="Optional hub diameter used to derive hub radius for RotorMesh export",
    )
    parser.add_argument(
        "--rotor-diameter",
        type=float,
        help="Optional rotor diameter used with blade span to derive hub radius for RotorMesh export",
    )
    parser.add_argument(
        "--hub-length",
        type=float,
        help="Optional axial length override for RotorHubMesh export",
    )
    parser.add_argument(
        "--hub-connector-radius",
        type=float,
        help="Optional connector-cylinder radius override for RotorHubMesh export",
    )
    parser.add_argument(
        "--hub-nose-radius",
        type=float,
        help="Optional spherical nose radius override for RotorHubMesh export",
    )
    parser.add_argument(
        "--element-size",
        type=float,
        default=0.1,
        help="Target mesh element size for direct mesh export (default: 0.1)",
    )
    parser.add_argument(
        "--n-samples",
        type=int,
        default=300,
        help="Airfoil resampling points used by direct mesh export (default: 300)",
    )
    parser.add_argument(
        "--span-grading",
        default="chord",
        help="Spanwise grading strategy for BladeMesh export (default: chord)",
    )
    parser.add_argument(
        "--airfoil-spacing",
        default="cosine",
        choices=("constant", "cosine", "half-cosine", "auto"),
        help="Airfoil sampling distribution used by direct mesh export",
    )
    parser.add_argument(
        "--renumber",
        choices=("simple", "rcm"),
        help="Optional mesh renumbering applied before export",
    )
    parser.add_argument(
        "--ccx-quadratic",
        action="store_true",
        default=False,
        help=(
            "Upgrade mesh to second-order elements (S8R/S6) for CalculiX export "
            "and use *SHELL SECTION, COMPOSITE with per-ply detail. "
            "Default: first-order S4/S3 with orthotropic equivalent material."
        ),
    )

    args = parser.parse_args(argv)

    # --template: no config needed
    if args.template:
        # Reuse the template from run_fsi to keep them in sync
        from aeroelast.cli.run_fsi import TEMPLATE_CONFIG

        print(TEMPLATE_CONFIG)
        return 0

    if args.export_mesh or args.export_parts_dir:
        setup_logging(args.verbose)
        try:
            return _export_mesh_from_cli(args)
        except Exception as e:
            logging.exception("Mesh export failed")
            print(f"\nError: {e}")
            return 1

    if not args.config:
        parser.print_help()
        return 1

    config_path = Path(args.config)
    if not config_path.exists():
        print(f"Error: configuration file not found: {config_path}")
        return 1

    setup_logging(args.verbose)

    # --preview
    if args.preview:
        from aeroelast.core.config import FSISimulationConfig

        cfg = FSISimulationConfig.from_yaml(str(config_path))
        print(cfg)
        return 0

    # --validate
    if args.validate:
        from aeroelast.cli.run_fsi import validate_config

        return 0 if validate_config(str(config_path)) else 1

    # --export-ccx
    if args.export_ccx:
        try:
            from aeroelast.solvers.fsi.runner import FSIRunner

            runner = FSIRunner(str(config_path), args.workdir)
            num_modes = 10
            if runner.config.solver.num_modes:
                num_modes = runner.config.solver.num_modes
            span_direction = None
            elem_cfg = getattr(runner.config, "elements", None)
            if elem_cfg is not None and getattr(elem_cfg, "span_direction", None) is not None:
                span_direction = tuple(float(v) for v in elem_cfg.span_direction)
            runner.export_calculix(
                args.export_ccx,
                num_modes=num_modes,
                span_direction=span_direction,
                quadratic=args.ccx_quadratic,
            )
            return 0
        except Exception as e:
            logging.exception("CalculiX export failed")
            print(f"\nError: {e}")
            return 1

    # Run simulation
    try:
        # Peek at the YAML to detect BEM-FSI configs without importing the
        # full runner (avoids the LinearDynamicFSI default path).
        import yaml as _yaml

        with open(config_path) as _f:
            _raw = _yaml.safe_load(_f)
        _solver_type = _raw.get("solver", {}).get("type")
        _is_aero_fsi = (
            (_solver_type == "AeroFSI")
            or ("aero" in _raw and _solver_type in (None, "AeroFSI"))
            or ("bem" in _raw and _solver_type in (None, "BEMFSI", "AeroFSI"))
        )
        if _is_aero_fsi:
            from aeroelast.cli.run_aero_fsi import main as _aero_main

            _argv = [str(config_path)]
            if args.workdir:
                _argv += ["--workdir", args.workdir]
            return _aero_main(_argv)

        from aeroelast.solvers.fsi.runner import FSIRunner

        runner = FSIRunner(str(config_path), args.workdir)
        runner.run()
        return 0

    except KeyboardInterrupt:
        print("\nInterrupted by user.")
        return 130

    except Exception as e:
        logging.exception("Simulation failed")
        print(f"\nError: {e}")
        return 1


if __name__ == "__main__":
    sys.exit(main())
