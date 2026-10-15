"""Shell geometric-nonlinear increment of the official blade tip under gravity vs a beam reference.

Runs the linear and one nonlinear scheme on the same mesh and prints the tip axial increment
``uz(nonlinear) - uz(linear)``, next to the second-order beam foreshortening of
``blade_gravity_beam_shortening`` on that mesh. Run one (element size, scheme) per process and
compare the trend over element sizes, not one mesh.

Usage: python tools/blade_gravity_nl_increment.py --element-size 0.5 --scheme newton
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(REPO_ROOT / "tools"))

from beam_reference import BeamReference, load_beamdyn_blade, load_elastodyn_blade  # noqa: E402
from blade_gravity_beam_shortening import foreshortening  # noqa: E402
from run_s2_static_cases import shell_tip_displacement  # noqa: E402

GRAVITY = (0.0, -9.81, 0.0)
STEPS = {"newton": 8, "ul": 64}


def tip_ring_mean(mesh, u) -> np.ndarray:
    """Mean translation of every node on the tip ring (z at its maximum): less mesh-sensitive
    than the single max-z node, whose offset from the beam axis couples to the section slope."""
    coords = mesh.coords_array
    dpn = len(u) // coords.shape[0]
    ring = np.nonzero(coords[:, 2] >= coords[:, 2].max() - 1e-6)[0]
    return np.mean([np.asarray(u[i * dpn : i * dpn + 3]) for i in ring], axis=0)


def _solve(cls, mesh, props, solver_cfg):
    from aeroelast.core.bc import BodyForce, DirichletCondition
    from aeroelast.elements import ElementFamily

    cfg = {
        "solver": solver_cfg,
        "elements": {
            "element_family": ElementFamily.SHELL,
            "span_direction": (0.0, 0.0, 1.0),
            "properties": props,
        },
    }
    solver = cls(mesh, cfg)
    dpn = solver.domain.dofs_per_node
    root = sorted(
        mesh.node_id_to_index[n] * dpn + d
        for n in mesh.get_node_set("RootNodes").nodes
        for d in range(dpn)
    )
    solver.add_dirichlet_conditions([DirichletCondition(root, 0.0)])
    solver.add_body_forces([BodyForce(list(GRAVITY))])
    start = time.time()
    u = solver.solve()
    return (shell_tip_displacement(mesh, u), tip_ring_mean(mesh, u)), time.time() - start


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--element-size", type=float, required=True)
    ap.add_argument("--scheme", choices=sorted(STEPS), required=True)
    ap.add_argument("--n-slices", type=int, default=60)
    args = ap.parse_args()

    from aeroelast.core.mesh.generators import BladeMesh
    from aeroelast.models.blade.model import build_rust_properties
    from aeroelast.solvers.elasticity.static_linear import StaticLinearSolver
    from aeroelast.solvers.elasticity.static_nonlinear import StaticNonlinearSolver

    gen = BladeMesh(
        excel_file=None,
        yaml_file=str(REPO_ROOT / "tests" / "IEA-15-240-RWT.yaml"),
        airfoil_dir=str(REPO_ROOT / "tests" / "airfoils"),
        element_size=args.element_size,
    )
    mesh = gen.generate(renumber="rcm")
    props = build_rust_properties(gen.numad_mesh_data)
    ref = REPO_ROOT / "tests" / "IEA15MW" / "reference"
    beam = BeamReference.from_mesh(
        mesh,
        load_elastodyn_blade(ref / "IEA-15-240-RWT_ElastoDyn_blade.dat"),
        load_beamdyn_blade(ref / "IEA-15-240-RWT_BeamDyn_blade.dat"),
        n_slices=args.n_slices,
    )
    ub = beam.solve_gravity(GRAVITY).reshape(-1, 6)
    dz_ref = foreshortening(beam.z, ub[:, 0], ub[:, 1])

    (lin, lin_ring), t_lin = _solve(StaticLinearSolver, mesh, props, {})
    steps = STEPS[args.scheme]
    (nl, nl_ring), t_nl = _solve(
        StaticNonlinearSolver,
        mesh,
        props,
        {"scheme": args.scheme, "continuation_steps": steps, "continuation_max_steps": 16 * steps},
    )
    dz = float(nl[2] - lin[2])
    dz_ring = float(nl_ring[2] - lin_ring[2])
    print(
        json.dumps(
            {
                "element_size": args.element_size,
                "n_nodes": int(mesh.coords_array.shape[0]),
                "scheme": args.scheme,
                "tip_linear": [float(x) for x in lin],
                "tip_nonlinear": [float(x) for x in nl],
                "dz_shell_nl_minus_lin": dz,
                "dz_beam_reference": dz_ref,
                "dz_relative_error": dz / dz_ref - 1.0,
                "dz_ring_mean": dz_ring,
                "dz_ring_mean_relative_error": dz_ring / dz_ref - 1.0,
                "tip_ring_mean_linear": [float(x) for x in lin_ring],
                "beam_tip_uy": float(ub[-1, 1]),
                "seconds_linear": t_lin,
                "seconds_nonlinear": t_nl,
            },
            sort_keys=True,
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
