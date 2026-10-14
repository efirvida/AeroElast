"""Compare static D-Tube solvers under one fixed global traction, not follower pressure.

Run one mode per process: PYTHONPATH=$SCRATCH/validenv:$PWD scripts/aeroenv.sh python
  tools/diagnose_dtube_static_parity.py {ae-linear,ae-newton,ccx-linear,ccx-nlgeom} --n-z 100 --scale 0.1
The CCX deck and solver logs remain under logs/check/dtube-static-parity/. Compare
at the same n_z, scale, cross-section, clamp, load resultant and mean arc-tip uy.
S8R and MITC4 have different interpolation: a fixed-mesh gap is not an accuracy verdict.
"""

from __future__ import annotations

import argparse
import json
import math
import shutil
import tempfile
from pathlib import Path

import numpy as np

from tests.support.ccx_io import run_ccx
from tests.support.ccx_surface import node_id, read_dat_displacements, write_s8r_surface_deck
from tests.validation.benchmarks import test_dtube_camarena_bending as dtube


def _rule_weight(index: int, last: int, rule: str) -> float:
    """Quadrature weight in units of the grid step: trapezoid, or Simpson (quadratic-consistent)."""
    if rule == "trapezoid":
        return 0.5 if index in (0, last) else 1.0
    if index in (0, last):
        return 1.0 / 3.0
    return 4.0 / 3.0 if index % 2 else 2.0 / 3.0


def fixed_arc_forces(
    *, n_arc: int, n_z: int, subdivisions: int, pressure: float, rule: str = "trapezoid"
) -> list[tuple[int, int, tuple[float, float, float]]]:
    """Reference-configuration integration of the radial dead load (trapezoid or Simpson)."""
    forces = []
    n_theta, n_span = n_arc * subdivisions, n_z * subdivisions
    dtheta, dz = math.pi / n_theta, dtube.L / n_span
    for k in range(n_span + 1):
        wz = dz * _rule_weight(k, n_span, rule)
        for i in range(n_theta + 1):
            theta = math.pi * i / n_theta
            ds = dtube.R * dtheta * _rule_weight(i, n_theta, rule)
            w = pressure * ds * wz
            forces.append((i, k, (-w * math.cos(theta), w * math.sin(theta), 0.0)))
    return forces


def tip_shear_forces(
    *, n_arc: int, n_z: int, subdivisions: int, force: float
) -> list[tuple[int, int, tuple[float, float, float]]]:
    """Control load: a known tip shear spread equally over the arc nodes of the free ring."""
    n_theta = n_arc * subdivisions
    share = force / (n_theta + 1)
    return [(i, n_z * subdivisions, (0.0, share, 0.0)) for i in range(n_theta + 1)]


RHO = 7850.0


def gravity_acceleration(scale: float) -> float:
    """Acceleration whose body load rho*a*A per unit span equals the scaled pressure resultant."""
    return scale * dtube.W_RES / (RHO * dtube.A_SEC)


def gravity_beam_reference(accel: float) -> float:
    """Euler-Bernoulli tip deflection of the cantilever under w = rho a A through the centroid."""
    return RHO * accel * dtube.A_SEC * dtube.L**4 / (8.0 * dtube.E * dtube.I_C)


def add_root_reaction_probe(inp: Path, n_u: int) -> None:
    """Request CalculiX RF on the original clamped ring for a load-transfer audit."""
    text = inp.read_text(encoding="utf-8")
    steps = [line for line in text.splitlines(keepends=True) if line.startswith("*STEP")]
    if len(steps) != 1 or text.count("*END STEP") != 1:
        raise ValueError("expected one static step for root reaction probe")
    ring = "\n".join(
        ", ".join(str(i) for i in range(start, min(start + 16, 2 * n_u + 1)))
        for start in range(1, 2 * n_u + 1, 16)
    )
    text = text.replace(steps[0], f"*NSET, NSET=ROOT_RF\n{ring}\n{steps[0]}", 1)
    text = text.replace("*END STEP", "*NODE PRINT, NSET=ROOT_RF\nRF\n*END STEP", 1)
    inp.write_text(text, encoding="utf-8")


def read_root_reaction(dat: Path) -> tuple[float, float, float]:
    """Sum the ROOT_RF block of a CalculiX ``.dat``: the equilibrium audit of the applied load."""
    total, active = [0.0, 0.0, 0.0], False
    for line in dat.read_text(encoding="utf-8").splitlines():
        if "forces (fx,fy,fz) for set ROOT_RF" in line:
            active = True
        elif active and "total force" in line:
            break
        elif active and len(line.split()) == 4:
            for k, token in enumerate(line.split()[1:]):
                total[k] += float(token)
    return total[0], total[1], total[2]


def _aeroelast(
    mode: str,
    n_z: int,
    pressure: float,
    tip_force: float | None = None,
    accel: float | None = None,
    n_arc: int = dtube.N_ARC,
) -> float:
    from aeroelast.core.bc import DirichletCondition, NodalLoad
    from aeroelast.core.material import IsotropicMaterial
    from aeroelast.core.properties import ShellProperty
    from aeroelast.elements import ElementFamily
    from aeroelast.solvers.elasticity.static_linear import StaticLinearSolver
    from aeroelast.solvers.elasticity.static_nonlinear import StaticNonlinearSolver

    mesh, arc, _ = dtube._build_dtube_mesh(n_arc=n_arc, n_web=n_arc // 2, n_z=n_z)
    prop = ShellProperty(
        material=IsotropicMaterial(name="Steel", E=dtube.E, nu=dtube.NU, rho=RHO),
        thickness=dtube.T,
    )
    cfg = {
        "solver": {"scheme": "newton", "continuation_steps": 8, "continuation_max_steps": 128}
        if mode == "ae-newton"
        else {},
        "elements": {
            "element_family": ElementFamily.SHELL,
            "properties": {"tube": prop},
            "span_direction": (0.0, 0.0, 1.0),
        },
    }
    cls = StaticNonlinearSolver if mode == "ae-newton" else StaticLinearSolver
    solver = cls(mesh, cfg)
    dpn = solver.domain.dofs_per_node
    clamp = sorted(
        mesh.node_id_to_index[node.id] * dpn + d
        for node in mesh.get_node_set("clamped").nodes.values()
        for d in range(dpn)
    )
    solver.add_dirichlet_conditions([DirichletCondition(clamp, 0.0)])
    if accel is not None:
        from aeroelast.core.bc import BodyForce

        solver.add_body_forces([BodyForce([0.0, accel, 0.0])])
    dofs, values = [], []
    for i, k, force in (
        ()
        if accel is not None
        else (
            fixed_arc_forces(n_arc=n_arc, n_z=n_z, subdivisions=1, pressure=pressure)
            if tip_force is None
            else tip_shear_forces(n_arc=n_arc, n_z=n_z, subdivisions=1, force=tip_force)
        )
    ):
        node = arc[k][i]
        assert node is not None
        base = mesh.node_id_to_index[node.id] * dpn
        dofs.extend((base, base + 1))
        values.extend((force[0], force[1]))
    if values:
        solver.add_nodal_loads([NodalLoad(dofs, values)])
    u = np.asarray(solver.solve(), dtype=float)
    tip_nodes = [node for node in arc[n_z] if node is not None]
    return float(np.mean([u[mesh.node_id_to_index[node.id] * dpn + 1] for node in tip_nodes]))


def _calculix(
    mode: str,
    n_z: int,
    pressure: float,
    root: Path,
    tip_force: float | None = None,
    rule: str = "trapezoid",
    accel: float | None = None,
    n_arc: int = dtube.N_ARC,
) -> tuple[float, tuple[float, float, float]]:
    ccx = shutil.which("ccx")
    if ccx is None:
        raise RuntimeError("CalculiX ccx not found")
    n_web = n_arc // 2
    n_u = n_arc + n_web

    def point(u: float, v: float) -> tuple[float, float, float]:
        s = u * n_u
        if s <= n_arc:
            theta = math.pi * s / n_arc
            return (dtube.R * math.cos(theta), -dtube.R * math.sin(theta), v * dtube.L)
        return (-dtube.R + 2 * dtube.R * (s - n_arc) / n_web, 0.0, v * dtube.L)

    if accel is not None:
        fields = []
    elif tip_force is not None:
        fields = tip_shear_forces(n_arc=n_arc, n_z=n_z, subdivisions=2, force=tip_force)
    else:
        fields = fixed_arc_forces(
            n_arc=n_arc, n_z=n_z, subdivisions=2, pressure=pressure, rule=rule
        )
    load_field = (
        None
        if mode == "ccx-pressure-control" or accel is not None
        else {node_id(i, k, n_u, periodic=True): force for i, k, force in fields}
    )
    tip_ids = [node_id(i, 2 * n_z, n_u, periodic=True) for i in range(2 * n_arc + 1)]
    run_dir = Path(tempfile.mkdtemp(prefix=f"{mode}_nz{n_z}_", dir=root))
    inp = run_dir / "dtube.inp"
    write_s8r_surface_deck(
        inp,
        point,
        n_u,
        n_z,
        young=dtube.E,
        poisson=dtube.NU,
        thickness=dtube.T,
        boundary=lambda x, y, z: [1, 2, 3, 4, 5, 6] if abs(z) < 1e-9 else [],
        periodic_u=True,
        load_field=load_field,
        pressure=(lambda ei, ej: ei < n_arc, pressure)
        if load_field is None and accel is None
        else None,
        density=RHO if accel is not None else 0.0,
        gravity=(0.0, accel, 0.0) if accel is not None else None,
        print_nodes=tip_ids,
        nlgeom=8 if mode == "ccx-nlgeom" else None,
    )
    add_root_reaction_probe(inp, n_u)
    result = run_ccx(inp, ccx)
    (run_dir / "ccx.stdout").write_text(result.stdout)
    (run_dir / "ccx.stderr").write_text(result.stderr)
    if result.returncode:
        raise RuntimeError(f"ccx exited {result.returncode}; see {run_dir}")
    disp = read_dat_displacements(inp.with_suffix(".dat"))
    print(f"CCX_DECK={inp}", flush=True)
    tip = float(np.mean([disp[nid][1] for nid in tip_ids]))
    return tip, read_root_reaction(inp.with_suffix(".dat"))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "mode",
        choices=(
            "ae-linear",
            "ae-newton",
            "ccx-linear",
            "ccx-nlgeom",
            "ccx-pressure-control",
            "ae-tipshear",
            "ccx-tipshear",
        ),
    )
    parser.add_argument(
        "--load",
        choices=("pressure", "gravity"),
        default="pressure",
        help="gravity = non-follower body load rho*a*A matched to the scaled pressure resultant",
    )
    parser.add_argument("--n-arc", type=int, default=dtube.N_ARC, help="arc elements; web = n/2")
    parser.add_argument("--rule", choices=("trapezoid", "simpson"), default="trapezoid")
    parser.add_argument("--tip-force", type=float, default=1.0e6, help="N, *-tipshear modes")
    parser.add_argument("--n-z", type=int, default=100, help="common axial element count")
    parser.add_argument("--scale", type=float, default=0.1, help="fraction of 63 kPa nominal load")
    args = parser.parse_args()
    if args.n_z < 2 or not math.isfinite(args.scale) or not 0 < args.scale <= 1:
        parser.error("require n-z >= 2 and 0 < scale <= 1")
    root = Path("logs/check/dtube-static-parity")
    root.mkdir(parents=True, exist_ok=True)
    pressure = dtube.P_TRACTION * args.scale
    accel = gravity_acceleration(args.scale) if args.load == "gravity" else None
    shear = args.tip_force if args.mode.endswith("tipshear") else None
    reaction = None
    if args.mode.startswith("ccx"):
        tip, reaction = _calculix(
            args.mode, args.n_z, pressure, root, shear, args.rule, accel, args.n_arc
        )
    else:
        tip = _aeroelast(args.mode, args.n_z, pressure, shear, accel, args.n_arc)
    loads = (
        []
        if accel is not None
        else fixed_arc_forces(n_arc=args.n_arc, n_z=args.n_z, subdivisions=1, pressure=pressure)
        if shear is None
        else tip_shear_forces(n_arc=args.n_arc, n_z=args.n_z, subdivisions=1, force=shear)
    )
    if accel is not None:
        beam = gravity_beam_reference(accel)
    elif shear is not None:
        beam = shear * dtube.L**3 / (3.0 * dtube.E * dtube.I_C)
    else:
        beam = dtube.DELTA_REF * args.scale
    print(
        json.dumps(
            {
                "load": args.load,
                "mode": args.mode,
                "n_arc": args.n_arc,
                "n_z": args.n_z,
                "scale": args.scale,
                "tip_uy_m": tip,
                "dead_load_y_N": RHO * accel * dtube.A_SEC * dtube.L
                if accel is not None
                else sum(f[1] for _, _, f in loads),
                "beam_linear_m": beam,
                "root_reaction_N": reaction,
            },
            sort_keys=True,
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
