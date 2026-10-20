"""V-02 mesh study: official IEA 15 MW blade, clamped-root shell modal frequencies vs element size.

Same construction as tests/validation/blade/test_iea15mw_v02_natural_frequencies.py (BladeMesh +
build_rust_properties + ModalSolver, root fixed). Prints the closest-to-reference frequency for the
first flap/edge modes and the 2nd flap/edge, per element size. Run through scripts/aeroenv.sh on a
compute node; es=0.125 needs ~12 GB.
"""
import sys
from pathlib import Path

import numpy as np

REF = {"1F": 0.5585, "1E": 0.6406, "2F": 1.659, "2E": 2.167}
YAML = Path("tests/IEA-15-240-RWT.yaml")


def modal(es: float) -> np.ndarray:
    from aeroelast.core.bc import DirichletCondition
    from aeroelast.core.mesh.generators import BladeMesh
    from aeroelast.elements import ElementFamily
    from aeroelast.models.blade.model import build_rust_properties
    from aeroelast.solvers.modal import ModalSolver

    gen = BladeMesh(yaml_file=str(YAML), element_size=es)
    mesh = gen.generate(renumber="rcm")
    props = build_rust_properties(gen.numad_mesh_data)
    cfg = {
        "solver": {"num_modes": 10},
        "elements": {"element_family": ElementFamily.SHELL, "span_direction": (0.0, 0.0, 1.0), "properties": props},
    }
    solver = ModalSolver(mesh, cfg)
    dpn = solver.domain.dofs_per_node
    ids = sorted(mesh.get_node_set("RootNodes").nodes.keys())
    idx = mesh.node_id_to_index
    solver.add_dirichlet_conditions([DirichletCondition(sorted(idx[n] * dpn + d for n in ids for d in range(dpn)), 0.0)])
    f, _ = solver.solve()
    print(f"es={es} nodes={len(mesh.nodes) if hasattr(mesh, 'nodes') else '?'} f[Hz]=" + " ".join(f"{x:.5f}" for x in np.asarray(f)), flush=True)
    return np.asarray(f)


if __name__ == "__main__":
    for es in (float(x) for x in sys.argv[1:]):
        f = modal(es)
        for k, r in REF.items():
            v = float(f[np.argmin(np.abs(f - r))])
            print(f"RESULT es={es} {k} {v:.6f} Hz vs {r} err {100 * (v - r) / r:+.4f} %", flush=True)
