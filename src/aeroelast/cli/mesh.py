#!/usr/bin/env python3
"""
aeroelast mesh — mesh generation CLI.

Generates blade, rotor, hub, nacelle and tower meshes from a WindIO turbine
YAML (or from explicit parameters) and writes STL / VTK / OBJ / MSH files.

Usage:
    aeroelast mesh blade   <input.yaml> --out blade.stl [--no-webs]
    aeroelast mesh rotor   <input.yaml> --out rotor.vtk --n-blades 3
    aeroelast mesh hub     [input.yaml] --out hub.stl [--diameter 7.94]
    aeroelast mesh nacelle [input.yaml] --out nacelle.obj [--length 12]
    aeroelast mesh tower   [input.yaml] --out tower.stl [--height 144]
    aeroelast mesh turbine <input.yaml> --out-dir meshes/ --format stl

The single-component commands infer the format from ``--out``; ``--format``
is the fallback when the path has no extension.  ``turbine`` writes one file
per component (blade_1..blade_N, hub, nacelle, tower) into ``--out-dir``.
"""

import argparse
import sys
from pathlib import Path


class MeshCliError(Exception):
    """Invalid component arguments; reported as a CLI error, not a crash."""


def _write_mesh(mesh, out: str, quiet: bool) -> None:
    from aeroelast.core.mesh.io.writers import write_mesh

    path = Path(out)
    write_mesh(mesh, str(path))
    if not quiet:
        print(f"wrote {path} ({mesh.node_count} nodes, {mesh.elements_count} elements)")


def _add_output(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--out", "-o", required=True, help="Output mesh file (format = extension)")
    parser.add_argument(
        "--format",
        default="stl",
        help="Format fallback when --out has no extension (default: stl)",
    )


def _add_surface_options(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--n-circ", type=int, default=32, help="Circumferential elements")
    parser.add_argument("--no-caps", action="store_true", help="Leave the open ends uncapped")


# ---------------------------------------------------------------------------
# commands
# ---------------------------------------------------------------------------


def _cmd_blade(args) -> None:
    from aeroelast.core.mesh.generators import BladeMesh

    mesh = BladeMesh(
        yaml_file=args.input,
        element_size=args.element_size,
        n_samples=args.n_samples,
        include_webs=not args.no_webs,
        span_grading=args.span_grading,
    ).generate(renumber=args.renumber, verbose=not args.quiet)
    _write_mesh(mesh, _output_path(args), args.quiet)


def _cmd_rotor(args) -> None:
    from aeroelast.core.mesh.generators import RotorMesh

    mesh = RotorMesh(
        yaml_file=args.input,
        n_blades=args.n_blades,
        hub_radius=args.hub_radius,
        element_size=args.element_size,
        n_samples=args.n_samples,
        include_webs=not args.no_webs,
    ).generate(renumber=args.renumber, verbose=not args.quiet)
    _write_mesh(mesh, _output_path(args), args.quiet)


def _cmd_hub(args) -> None:
    from aeroelast.core.mesh.components import HubMesh

    if args.diameter is not None:
        hub = HubMesh.from_params(args.diameter, n_circ=args.n_circ, n_merid=args.n_merid)
    elif args.input:
        hub = HubMesh.from_windio(args.input, n_circ=args.n_circ, n_merid=args.n_merid)
    else:
        raise MeshCliError("hub needs an INPUT yaml or --diameter")
    _write_mesh(hub.generate(), _output_path(args), args.quiet)


def _cmd_nacelle(args) -> None:
    from aeroelast.core.mesh.components import NacelleMesh

    if args.length is not None or args.body_diameter is not None:
        if args.length is None or args.body_diameter is None:
            raise MeshCliError("--length and --body-diameter must be given together")
        nacelle = NacelleMesh.from_params(
            length=args.length,
            body_diameter=args.body_diameter,
            nose_diameter=args.nose_diameter,
            n_circ=args.n_circ,
            n_axial=args.n_axial,
            n_nose=args.n_nose,
            cap_tail=not args.no_caps,
        )
    elif args.input:
        nacelle = NacelleMesh.from_windio(
            args.input,
            n_circ=args.n_circ,
            n_axial=args.n_axial,
            n_nose=args.n_nose,
            cap_tail=not args.no_caps,
        )
    else:
        raise MeshCliError("nacelle needs an INPUT yaml or --length/--body-diameter")
    _write_mesh(nacelle.generate(), _output_path(args), args.quiet)


def _cmd_tower(args) -> None:
    from aeroelast.core.mesh.components import TowerMesh

    if args.height is not None or args.base_diameter is not None:
        if args.height is None or args.base_diameter is None:
            raise MeshCliError("--height and --base-diameter must be given together")
        tower = TowerMesh.from_params(
            height=args.height,
            base_diameter=args.base_diameter,
            top_diameter=args.top_diameter,
            n_axial=args.n_axial,
            n_circ=args.n_circ,
            cap_start=not args.no_caps,
            cap_end=not args.no_caps,
        )
    elif args.input:
        tower = TowerMesh.from_windio(
            args.input,
            n_circ=args.n_circ,
            n_axial=args.n_axial,
            element_size=args.element_size,
            cap_start=not args.no_caps,
            cap_end=not args.no_caps,
        )
    else:
        raise MeshCliError("tower needs an INPUT yaml or --height/--base-diameter")
    _write_mesh(tower.generate(), _output_path(args), args.quiet)


def _cmd_turbine(args) -> None:
    from aeroelast.core.mesh.turbine import TurbineMesh

    result = TurbineMesh(
        yaml_file=args.input,
        n_blades=args.n_blades,
        element_size=args.element_size,
        n_samples=args.n_samples,
        include_webs=args.webs,
        span_grading=args.span_grading,
        tower_element_size=args.tower_element_size,
        tower_n_axial=args.tower_n_axial,
        tower_n_circ=args.tower_n_circ,
        hub_n_circ=args.hub_n_circ,
        hub_n_merid=args.hub_n_merid,
        nacelle_n_circ=args.nacelle_n_circ,
        tower_base_z=args.tower_base_z,
    ).generate(renumber=args.renumber, verbose=not args.quiet)

    written = result.write(args.out_dir, format=args.format, prefix=args.prefix)
    if not args.quiet:
        for name, path in written.items():
            print(f"wrote {name}: {path}")


def _output_path(args) -> str:
    path = Path(args.out)
    if path.suffix:
        return str(path)
    return f"{path}.{args.format}"


# ---------------------------------------------------------------------------
# parser
# ---------------------------------------------------------------------------


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="aeroelast mesh",
        description="Generate blade, rotor, hub, nacelle and tower meshes.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  aeroelast mesh blade turbine.yaml --out blade.stl --element-size 0.05 --no-webs
  aeroelast mesh rotor turbine.yaml --out rotor.vtk --n-blades 3
  aeroelast mesh hub   turbine.yaml --out hub.stl
  aeroelast mesh tower --height 144 --base-diameter 10 --top-diameter 6.5 --out tower.obj
  aeroelast mesh turbine turbine.yaml --out-dir meshes/ --format stl
        """,
    )
    sub = parser.add_subparsers(dest="component", required=True)

    blade = sub.add_parser("blade", help="Blade shell/surface mesh from a WindIO YAML")
    blade.add_argument("input", help="WindIO turbine YAML")
    blade.add_argument("--element-size", type=float, default=0.1)
    blade.add_argument("--n-samples", type=int, default=300)
    blade.add_argument("--no-webs", action="store_true", help="Surface only (CFD outer mold line)")
    blade.add_argument("--span-grading", choices=["chord", "uniform"], default="chord")
    blade.add_argument("--renumber", choices=["simple", "rcm"], default=None)
    blade.add_argument("--quiet", action="store_true")
    _add_output(blade)
    blade.set_defaults(func=_cmd_blade)

    rotor = sub.add_parser("rotor", help="Multi-blade rotor mesh")
    rotor.add_argument("input", help="WindIO turbine YAML")
    rotor.add_argument("--n-blades", type=int, required=True)
    rotor.add_argument("--hub-radius", type=float, default=None)
    rotor.add_argument("--element-size", type=float, default=0.1)
    rotor.add_argument("--n-samples", type=int, default=300)
    rotor.add_argument("--no-webs", action="store_true", help="Surface only (CFD)")
    rotor.add_argument("--renumber", choices=["simple", "rcm"], default=None)
    rotor.add_argument("--quiet", action="store_true")
    _add_output(rotor)
    rotor.set_defaults(func=_cmd_rotor)

    hub = sub.add_parser("hub", help="Hub spheroid")
    hub.add_argument("input", nargs="?", help="WindIO turbine YAML")
    hub.add_argument("--diameter", type=float, default=None)
    hub.add_argument("--n-merid", type=int, default=16)
    hub.add_argument("--quiet", action="store_true")
    _add_surface_options(hub)
    _add_output(hub)
    hub.set_defaults(func=_cmd_hub)

    nacelle = sub.add_parser("nacelle", help="Nacelle capsule")
    nacelle.add_argument("input", nargs="?", help="WindIO turbine YAML")
    nacelle.add_argument("--length", type=float, default=None)
    nacelle.add_argument("--body-diameter", type=float, default=None)
    nacelle.add_argument("--nose-diameter", type=float, default=None)
    nacelle.add_argument("--n-axial", type=int, default=12)
    nacelle.add_argument("--n-nose", type=int, default=6)
    nacelle.add_argument("--quiet", action="store_true")
    _add_surface_options(nacelle)
    _add_output(nacelle)
    nacelle.set_defaults(func=_cmd_nacelle)

    tower = sub.add_parser("tower", help="Tapered tower")
    tower.add_argument("input", nargs="?", help="WindIO turbine YAML")
    tower.add_argument("--height", type=float, default=None)
    tower.add_argument("--base-diameter", type=float, default=None)
    tower.add_argument("--top-diameter", type=float, default=None)
    tower.add_argument("--n-axial", type=int, default=20)
    tower.add_argument("--element-size", type=float, default=None)
    tower.add_argument("--quiet", action="store_true")
    _add_surface_options(tower)
    _add_output(tower)
    tower.set_defaults(func=_cmd_tower)

    turbine = sub.add_parser("turbine", help="Blades + hub + nacelle + tower, separate files")
    turbine.add_argument("input", help="WindIO turbine YAML")
    turbine.add_argument("--out-dir", "-d", required=True, help="Output directory")
    turbine.add_argument("--format", default="stl", help="Output format (default: stl)")
    turbine.add_argument("--prefix", default="", help="Filename prefix")
    turbine.add_argument("--n-blades", type=int, default=None)
    turbine.add_argument("--element-size", type=float, default=0.5)
    turbine.add_argument("--tower-element-size", type=float, default=None)
    turbine.add_argument("--n-samples", type=int, default=300)
    turbine.add_argument("--span-grading", choices=["chord", "uniform"], default="chord")
    turbine.add_argument("--webs", action="store_true", help="Include shear webs (default: off)")
    turbine.add_argument("--tower-base-z", type=float, default=0.0)
    turbine.add_argument("--tower-n-axial", type=int, default=20)
    turbine.add_argument("--tower-n-circ", type=int, default=32)
    turbine.add_argument("--hub-n-circ", type=int, default=32)
    turbine.add_argument("--hub-n-merid", type=int, default=16)
    turbine.add_argument("--nacelle-n-circ", type=int, default=32)
    turbine.add_argument("--renumber", choices=["simple", "rcm"], default=None)
    turbine.add_argument("--quiet", action="store_true")
    turbine.set_defaults(func=_cmd_turbine)

    return parser


def main(argv: list[str] | None = None) -> int:
    """Entry point for ``aeroelast mesh``."""
    parser = _build_parser()
    args = parser.parse_args(argv)

    try:
        args.func(args)
        return 0
    except SystemExit:
        raise
    except MeshCliError as exc:
        print(f"Error: {exc}")
        return 1
    except KeyboardInterrupt:
        print("\nInterrupted by user.")
        return 130
    except Exception as exc:  # noqa: BLE001 - CLI boundary
        print(f"Error: {exc}")
        return 1


if __name__ == "__main__":
    sys.exit(main())
