"""T18 follow-up: is the production dump's nodal F_AERO the minimum-norm (wall-flow off) or the wall-flow distribution?

One dump of the T6 run, compared node by node with the production participant (rebuilt branch, theta = the dump's)
built with ``element_properties=None`` (what ``_effective_element_properties`` returns for fluid.yaml) and with the
deck property map (wall-flow realisation on).  Output: ``<out>/nodal.json``.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "tools"))

import l5_oneway_bd_parity as par  # noqa: E402
import l5_t18_loads as lo  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--dump", type=Path, required=True, help="results/<t> directory")
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    import meshio

    from aeroelast.solvers.bem.fsi_participant import build_from_config
    from aeroelast.solvers.fsi.corotational import CoordinateTransforms

    mesh, props, asm = par.build_shell(20.0)
    m = meshio.read(args.dump / "fields.vtu")
    th = float(np.load(args.dump / "state.npz")["theta"])
    fa = np.asarray(m.point_data["F_AERO"])[:, :3]
    ct = CoordinateTransforms(rotation_axis=par.AXIS, rotation_center=np.array([0.0, 0.0, -3.97]))
    fa_r = ct.to_rotating(fa, th)
    out = {"dump": str(args.dump), "theta_deg": float(np.degrees(th) % 360), "dump_sum": fa_r.sum(0).tolist(),
           "dump_n_nonzero": int(np.count_nonzero(np.abs(fa_r).sum(1))), "variants": {}}
    u = np.zeros((mesh.node_count, 3))
    u[:, 0] = 1e-6
    for label, ep in (("wall_off_production", None), ("wall_on", props)):
        p = build_from_config(mesh, lo.production_cfg(), element_properties=ep)
        p._partner_theta = th
        f, _ = p._compute_forces(u)
        d = f - fa_r
        out["variants"][label] = {"sum": f.sum(0).tolist(), "n_nonzero": int(np.count_nonzero(np.abs(f).sum(1))),
                                  "rel_norm": float(np.linalg.norm(d) / np.linalg.norm(fa_r)),
                                  "rel_max_node": float(np.abs(d).max() / np.abs(fa_r).max()),
                                  "abs_max_node_N": float(np.abs(d).max()),
                                  "max_dump_node_N": float(np.abs(fa_r).max())}
        print(label, out["variants"][label], flush=True)
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "nodal.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
