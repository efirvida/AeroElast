"""S-4 mesh study: official IEA 15 MW blade, rotating modal (K + K_G + K_SP) vs element size.

Same construction as tests/test_iea15mw_s4_rotating_modal.py (BladeMesh + build_rust_properties +
assemble_rotating from tools/run_s4_rotating_modal.py, clamped root, mode families by beam projection, classify_modes),
but with the element size as an argument so the 0.5 m default can be shown to be unconverged.
Prints 1st flap / 1st edge / 2nd flap per element size and rpm against OpenFAST v5.0.0 MBC3,
preceded by one MODE line per computed mode (frequency, family label, tip translations).
Usage: scripts/aeroenv.sh python tools/s4_modal_convergence.py 0.5 0.25 [0.125] [--rpm=7.56]
es=0.125 needs --ntasks=12 on SLURM (OOM with one task).
"""
import sys
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "tools"))

OF = {
    0.0: {"1st flap": 0.544, "1st edge": 0.739, "2nd flap": 1.610},
    7.56: {"1st flap": 0.5666, "1st edge": 0.7461, "2nd flap": 1.6372},
}
FAMILIES = [("1st flap", "flap", 1), ("1st edge", "edge", 1), ("2nd flap", "flap", 2)]


def run(es: float, rpm: float) -> dict:
    from run_s4_rotating_modal import assemble_rotating, classify_modes, omega_from_rpm

    from aeroelast.core.assembler import MeshAssembler
    from aeroelast.core.mesh.generators import BladeMesh
    from aeroelast.elements import ElementFamily
    from aeroelast.models.blade.model import build_rust_properties

    gen = BladeMesh(element_size=es, yaml_file=str(REPO_ROOT / "tests/IEA-15-240-RWT.yaml"))
    mesh = gen.generate(renumber="rcm")
    props = build_rust_properties(gen.numad_mesh_data)
    model = {"elements": {"element_family": ElementFamily.SHELL, "span_direction": (0.0, 0.0, 1.0), "properties": props}}
    asm = MeshAssembler(mesh, model)
    freqs, modes, free = assemble_rotating(asm, mesh, omega_from_rpm(rpm), n_modes=12)
    labels = classify_modes(mesh, modes, free)
    # The whole spectrum with the labels and the tip components, so a mode that jumps family
    # between meshes can be traced (S-4 oml 2F, es 0.25 -> 0.125).
    xyz = np.asarray([n.coords for n in mesh.nodes])
    tip_dof0 = mesh.node_id_to_index[int(np.argmax(xyz[:, 2]))] * 6
    for k, (f, lab) in enumerate(zip(freqs, labels, strict=True)):
        u = [0.0, 0.0, 0.0]
        for d in range(3):
            pos = np.searchsorted(free, tip_dof0 + d)
            if pos < len(free) and free[pos] == tip_dof0 + d:
                u[d] = abs(float(modes[k, pos]))
        print(f"MODE es={es} rpm={rpm} k={k} {float(f):.6f} Hz {lab:<5} tip |ux|={u[0]:.3e} |uy|={u[1]:.3e} |uz|={u[2]:.3e}", flush=True)
    out = {}
    for name, direction, order in FAMILIES:
        cand = sorted(f for f, lab in zip(freqs, labels) if lab == direction)
        out[name] = float(cand[order - 1]) if len(cand) >= order else float("nan")
    return out


if __name__ == "__main__":
    # Optional --rpm=<value> restricts the run to one speed (es=0.125 at both speeds exceeds 90 min).
    rpms = tuple(float(a.split("=", 1)[1]) for a in sys.argv[1:] if a.startswith("--rpm=")) or (0.0, 7.56)
    for es in (float(x) for x in sys.argv[1:] if not x.startswith("--")):
        for rpm in rpms:
            res = run(es, rpm)
            for name, f in res.items():
                ref = OF[rpm][name]
                print(f"RESULT es={es} rpm={rpm} {name} {f:.6f} Hz vs {ref} err {100 * (f - ref) / ref:+.4f} %", flush=True)
