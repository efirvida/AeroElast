"""T18: the production coupling mesh (``coupling_node_set: allOuterShellNods``) vs the full shell mesh the T3 tool loads.

The production Fluid participant runs on the runner's own mesh (``_build_mesh``: node-set filtered); its forces reach
the Solid by nearest-neighbour mapping, so web nodes receive nothing.  This tool builds that participant exactly as
``run_bem_fsi.main`` does, takes its azimuth-mean forces on the rebuilt branch, maps them onto the full mesh by
coordinates, and solves the same K_eff: it must reproduce the T6 ``aero_t6`` row (-3.204 deg) from the production
config alone.  Output ``<dir>/subset.json`` (+ ``subset.npz``).
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "tools"))

import l5_oneway_bd_parity as par  # noqa: E402
import l5_t18_loads as lo  # noqa: E402
import twist_p1 as tp  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--dir", type=Path, required=True)
    args = ap.parse_args()
    t0 = time.time()
    from scipy.spatial import cKDTree

    from aeroelast.cli.run_bem_fsi import _build_mesh, _effective_element_properties
    from aeroelast.solvers.bem.fsi_participant import BEMFSIParticipant, build_from_config
    from run_s5_oneway import solve_static
    from tests.support.openfast_bem import build_blade_aero_from_aerodyn

    cfg = lo.production_cfg()
    cmesh, viz, ep = _build_mesh(cfg, lo.FLUID_YAML)
    prod = build_from_config(cmesh, cfg, viz_mesh=viz, element_properties=_effective_element_properties(cfg, ep))
    mesh, props, asm = par.build_shell(20.0)
    dpn = asm.dofs_per_node
    Xc, Xf = np.asarray(prod._ref_coords), np.asarray(mesh.coords_array)
    dist, idx = cKDTree(Xf).query(Xc)
    info = {"coupling_nodes": int(len(Xc)), "full_mesh_nodes": int(len(Xf)), "max_match_dist_m": float(dist.max()),
            "unique_targets": int(len(np.unique(idx)))}
    print(info, flush=True)
    u = np.zeros((len(Xc), 3))
    u[:, 0] = 1e-6
    az = np.radians(np.arange(lo.N_AZ) * 360.0 / lo.N_AZ)
    acc = np.zeros((len(Xc), 3))
    for a in az:
        prod._partner_theta = float(a)
        acc += prod._compute_forces(u)[0] / len(az)
    f3 = np.zeros((len(Xf), 3))
    f3[idx] = acc  # nearest-neighbour conservative write: 1:1 on identical vertices
    print(f"forces done ({time.time() - t0:.0f} s) sum {f3.sum(0)}", flush=True)
    f = np.zeros(asm.dofs_count)
    for d in range(3):
        f[d::dpn] = f3[:, d]
    K, free = par.rotating_keff(asm, mesh, par.OMEGA, 3.97)
    uu = solve_static(asm, K, free, f)
    U = np.column_stack([uu[d::dpn] for d in range(3)])
    part = BEMFSIParticipant(mesh, build_blade_aero_from_aerodyn(par.AD_PRIMARY), par.bem_config(-6.0, -4.0), element_properties=props)
    kp_bd = tp.read_keypoints((par.DECK_DIR / "IEA-15-240-RWT_BeamDyn.dat").read_text())[:, :3]
    zq = np.array(json.loads((par.REPO / "logs/check/l5_oneway/t11/bem.json").read_text())["table"]["z"])
    sh = par.shell_sections(mesh, part, U, kp_bd)
    r = {"flap": np.interp(zq, sh["z"], sh["u"][:, 1]).tolist(), "edge": np.interp(zq, sh["z"], sh["u"][:, 0]).tolist(),
         "twist_eq_deg": np.rad2deg(np.interp(zq, sh["z"], sh["twist_eq"])).tolist(),
         "twist_skin_deg": np.rad2deg(np.interp(zq, sh["z"], sh["twist_skin"])).tolist(), "F": f3.sum(0).tolist()}
    print(f"coupling-subset azmean rebuilt: flap {r['flap'][-1]:.4f} edge {r['edge'][-1]:.4f} tw_eq {r['twist_eq_deg'][-1]:.3f} tw_skin {r['twist_skin_deg'][-1]:.3f}")
    (args.dir / "subset.json").write_text(json.dumps({"info": info, "z": zq.tolist(), "case": r, "revision": tp.revision(par.REPO)}, indent=1))
    np.savez_compressed(args.dir / "subset.npz", f3=f3)


if __name__ == "__main__":
    main()
