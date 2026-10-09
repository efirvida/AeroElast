#!/usr/bin/env python3
"""CLI entry point for generic aerodynamic FSI participants.

This command runs an aerodynamic participant as a standalone preCICE process.
The YAML can use either the new generic ``aero`` section or the legacy
top-level ``bem`` section.
"""

from __future__ import annotations

import argparse
import logging
import os
import sys
from pathlib import Path

import yaml

from aeroelast.cli.run_bem_fsi import _build_mesh
from aeroelast.solvers.aero import build_aero_participant_from_config

logger = logging.getLogger(__name__)

_TEMPLATE = """\
# Generic aerodynamic FSI participant configuration
# =================================================
# Launch with: fem-shell-aero-fsi simulation.yaml

solver:
    type: "AeroFSI"

participant: "Fluid"
config_file: "precice-config.xml"
coupling_mesh: "Fluid-Mesh"

aero:
    backend: "bem"  # or "vlm" or "sharpy"
    blade_file: "path/to/blade.yaml"
    rotor:
        n_blades: 3
        hub_radius: 3.0
        rotation_axis: [0.0, 1.0, 0.0]
        rotation_center: [0.0, 0.0, 0.0]
    backend_config:
        wind_speed: 10.59
        omega: 7.56
        pitch: 0.0
        air_density: 1.225
        dynamic_viscosity: 1.81e-5
        hub_height: 150.0
        shear_exp: 0.2
        default_re: 1.0e7
        neuralfoil_model: "large"
        precone: 0.0
        tilt: 0.0
        span_direction: [0.0, 0.0, 1.0]
        normal_direction: [1.0, 0.0, 0.0]
        tangential_direction: [0.0, 1.0, 0.0]
        # For backend: "vlm"
        # wake_length_scale: 25.0
        # wake_history_steps: 8
        # wake_corrector_iterations: 1
        # diagonal_regularization: 1.0e-8
        # core_radius: 1.0e-8
        # vertical_axis: [0.0, 0.0, 1.0]
        # For backend: "sharpy"
        # model: "uvlm"              # "uvlm" or "vlm"
        # case_adapter: "pkg.module.AdapterClass"
        # solver_settings:
        #   convection_scheme: 3
        #   dt: 0.05
        # Dynamic inflow (linearly interpolated in time)
        # wind:
        #   type: "table"
        #   times: [0.0, 1.0, 2.0]
        #   velocities:
        #     - [10.0, -1.0, 0.0]
        #     - [13.0, -2.0, 0.0]
        #     - [9.0, 0.5, 0.0]
        #   # Or load the table from CSV:
        #   # file: "wind_timeseries.csv"
        #   # time_column: "time"
        #   # velocity_columns: ["u", "v", "w"]

mesh:
    source: "generator"
    generator:
        type: "RotorMesh"
        params:
            yaml_file: "path/to/blade.yaml"
            element_size: 0.5

output:
    folder: "aero_fsi_results"
    log_interval: 10  # windows between summary log lines
    sectional_bins: 8
    # Writes <folder>/aero_report_schema.json describing report files,
    # schema versions, and stable field ordering.
    # Writes <folder>/aero_report.csv with schema version plus rotor-global and
    # blade-wise force, moment, wind, and wake diagnostics.
    # Writes <folder>/aero_spanwise_report.csv with schema version, stable column
    # order, and per-panel normalized radial/span coordinates.
    # Writes <folder>/aero_sectional_report.csv with stable column order and
    # radial-bin aggregates for validation-facing postprocessing.
"""


def _resolve_config_paths(cfg: dict, config_path: Path) -> None:
    """Resolve config file references relative to the YAML location."""
    yaml_dir = config_path.parent

    def _resolve(val: str) -> str:
        p = Path(val)
        return str((yaml_dir / p).resolve()) if not p.is_absolute() else val

    if "config_file" in cfg:
        cfg["config_file"] = _resolve(cfg["config_file"])
    if "blade_file" in cfg:
        cfg["blade_file"] = _resolve(cfg["blade_file"])

    aero_cfg = cfg.get("aero")
    if isinstance(aero_cfg, dict) and "blade_file" in aero_cfg:
        aero_cfg["blade_file"] = _resolve(aero_cfg["blade_file"])
    if isinstance(aero_cfg, dict):
        backend_cfg = aero_cfg.get("backend_config")
        if isinstance(backend_cfg, dict):
            wind_cfg = backend_cfg.get("wind")
            if isinstance(wind_cfg, dict) and "file" in wind_cfg:
                wind_cfg["file"] = _resolve(wind_cfg["file"])

    mesh_gen = cfg.get("mesh", {}).get("generator", {}).get("params", {})
    if "excel_file" in mesh_gen:
        mesh_gen["excel_file"] = _resolve(mesh_gen["excel_file"])
    if "airfoil_dir" in mesh_gen:
        mesh_gen["airfoil_dir"] = _resolve(mesh_gen["airfoil_dir"])


def _resolve_output_path_to_workdir(path_value: str, workdir: Path) -> Path:
    """Resolve output paths relative to the active case working directory."""
    path = Path(path_value)
    return path if path.is_absolute() else (workdir / path).resolve()


def _pin_relative_mesh_output_to_workdir(cfg: dict, workdir: Path) -> None:
    """Rewrite relative mesh.output_file paths to land in the active case dir."""
    mesh_cfg = cfg.get("mesh")
    if not isinstance(mesh_cfg, dict):
        return

    output_file = mesh_cfg.get("output_file")
    if not output_file:
        return

    output_path = Path(output_file)
    if output_path.is_absolute():
        return

    mesh_cfg["output_file"] = str(_resolve_output_path_to_workdir(output_file, workdir))


def _derive_sharpy_coupling_mesh_output(aero_mesh_path: Path) -> Path:
    """Keep the SHARPy aerogrid on the requested path and move the 3D coupling mesh aside."""
    suffix = aero_mesh_path.suffix or ".vtu"
    stem = aero_mesh_path.stem if aero_mesh_path.suffix else aero_mesh_path.name
    return aero_mesh_path.with_name(f"{stem}_coupling{suffix}")


def _configure_sharpy_mesh_exports(
    cfg: dict,
    workdir: Path,
    resolved_backend: str | None,
) -> None:
    """Reserve mesh.output_file for the SHARPy aerogrid and spill the coupling mesh to a sibling file."""
    if str(resolved_backend).strip().lower() != "sharpy":
        return

    mesh_cfg = cfg.get("mesh")
    if not isinstance(mesh_cfg, dict):
        return

    raw_mesh_output = mesh_cfg.get("output_file")
    if not raw_mesh_output:
        return

    output_cfg = cfg.get("output")
    if not isinstance(output_cfg, dict):
        output_cfg = {}
        cfg["output"] = output_cfg

    aero_mesh_path = _resolve_output_path_to_workdir(
        str(output_cfg.get("aero_mesh_file", raw_mesh_output)),
        workdir,
    )

    raw_coupling_output = output_cfg.get("coupling_mesh_file")
    if raw_coupling_output:
        coupling_mesh_path = _resolve_output_path_to_workdir(str(raw_coupling_output), workdir)
    else:
        coupling_mesh_path = _derive_sharpy_coupling_mesh_output(aero_mesh_path)

    if coupling_mesh_path == aero_mesh_path:
        coupling_mesh_path = _derive_sharpy_coupling_mesh_output(aero_mesh_path)

    mesh_cfg["output_file"] = str(coupling_mesh_path)
    output_cfg["aero_mesh_file"] = str(aero_mesh_path)
    output_cfg["coupling_mesh_file"] = str(coupling_mesh_path)

    logger.info(
        "[AERO-FSI] SHARPy mesh export remapped: aerogrid=%s coupling=%s",
        aero_mesh_path,
        coupling_mesh_path,
    )


def main(argv=None, *, default_backend: str | None = None) -> int:
    """CLI entry point."""
    parser = argparse.ArgumentParser(
        prog="fem-shell-aero-fsi",
        description=(
            "Generic aerodynamic preCICE participant for FSI coupling. "
            "Runs as an independent process alongside the structural solver."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "examples:\n"
            "  fem-shell-aero-fsi simulation.yaml\n"
            "  fem-shell-aero-fsi simulation.yaml --workdir /path/to/case\n"
            "  fem-shell-aero-fsi --template > simulation.yaml\n"
        ),
    )
    parser.add_argument(
        "config", nargs="?", metavar="CONFIG_FILE", help="Path to the YAML configuration file."
    )
    parser.add_argument(
        "--workdir", metavar="DIR", help="Working directory (default: directory of CONFIG_FILE)."
    )
    parser.add_argument(
        "--template", action="store_true", help="Print a template YAML configuration and exit."
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

    workdir = Path(args.workdir).resolve() if args.workdir else config_path.parent
    logging.info("[AERO-FSI] Working directory: %s", workdir)

    with open(config_path) as f:
        cfg = yaml.safe_load(f)

    _resolve_config_paths(cfg, config_path)
    aero_cfg = cfg.get("aero") if isinstance(cfg, dict) else None
    resolved_backend = (
        aero_cfg.get("backend", default_backend) if isinstance(aero_cfg, dict) else default_backend
    )
    _pin_relative_mesh_output_to_workdir(cfg, workdir)
    _configure_sharpy_mesh_exports(cfg, workdir, resolved_backend)
    logger.info(
        "[AERO-FSI] Loaded config: %s (backend=%s)",
        config_path,
        resolved_backend if resolved_backend is not None else "auto",
    )
    os.chdir(workdir)

    try:
        logger.info("[AERO-FSI] Building coupling mesh...")
        mesh, viz_mesh, _ = _build_mesh(cfg, config_path)
        logger.info(
            "[AERO-FSI] Mesh build completed: nodes=%d elements=%d",
            len(mesh.nodes),
            len(mesh.elements),
        )
    except Exception as exc:
        logging.error("Failed to build mesh: %s", exc)
        return 1

    try:
        logger.info("[AERO-FSI] Building aerodynamic participant instance...")
        participant = build_aero_participant_from_config(
            mesh,
            cfg,
            viz_mesh=viz_mesh,
            default_backend=default_backend,
        )
        runtime_context = getattr(participant, "runtime_context", None)
        if runtime_context is None:
            runtime_context = getattr(participant, "aero_runtime_context", None)
        if runtime_context is not None:
            logger.info(
                "[AERO-FSI] Participant ready: class=%s backend=%s",
                type(participant).__name__,
                getattr(runtime_context, "backend", "unknown"),
            )
        else:
            logger.info(
                "[AERO-FSI] Participant ready: class=%s",
                type(participant).__name__,
            )
        logger.info("[AERO-FSI] Entering preCICE coupling loop...")
        participant.run()
        logger.info("[AERO-FSI] Coupling loop finished cleanly.")
    except Exception as exc:
        logging.error("AERO-FSI participant failed: %s", exc, exc_info=True)
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
