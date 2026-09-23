"""Quick LC3 edgewise convergence probe: element_size 0.25 vs 0.5."""
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tools"))

from tools.beam_reference import load_elastodyn_blade, load_beamdyn_blade, BeamReference  # noqa: E402
from tools.run_s2_static_cases import run_shell_case, shell_tip_displacement  # noqa: E402
from aeroelast.core.mesh.generators import BladeMesh  # noqa: E402
from aeroelast.models.blade.model import build_rust_properties  # noqa: E402

ED = ROOT / "tests/IEA15MW/reference/IEA-15-240-RWT_ElastoDyn_blade.dat"
BD = ROOT / "tests/IEA15MW/reference/IEA-15-240-RWT_BeamDyn_blade.dat"

for h in (0.5, 0.25):
    gen = BladeMesh(excel_file=str(ROOT / "tests/NuMAD_utd_iea15mw.xlsx"),
                    airfoil_dir=str(ROOT / "tests/airfoils"), element_size=h)
    mesh = gen.generate(renumber="rcm")
    props = build_rust_properties(gen.numad_mesh_data)
    ed = load_elastodyn_blade(ED)
    bd = load_beamdyn_blade(BD)
    beam = BeamReference.from_mesh(mesh, ed, bd)
    load = {"type": "tip", "value": (1e6, 0.0, 0.0)}
    u = run_shell_case(mesh, props, load)
    d = shell_tip_displacement(mesh, u)
    db = np.asarray(beam.solve_tip_load(load["value"])[-6:-3])
    print(f"h={h}: shell ux={d[0]:.3f} uy={d[1]:.3f} uz={d[2]:.3f} | "
          f"beam ux={db[0]:.3f} | ratio={d[0]/db[0]:.3f}", flush=True)
