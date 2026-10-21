"""B3-3 M1: frame-invariant bending compliance of the shell vs BeamDyn (diagnostic, no verdict).

S-1 extracts EI from one displacement component per load (global +x for edge, +y for flap) and
compares it with the PRINCIPAL EIs of the 6x6 deck; S-2 LC3 compares the tip ux under a global +x
load. Both read one global axis. This script keeps both components of both loads, so the 2x2
compliance it reports can be rotated freely and its eigenvalues do not depend on the frame
(odd/tasks/validation-l3-b3-audit.md, "B3-3 rejected-rows analysis", K12.M1).

Four linear static solves on the mid-surface official blade (clamped root, 1 MN dead loads):
  S-2 style: load split equally over the max-z ring nodes (tools/run_s2_static_cases.py), +x, +y;
  S-1 style: load split over the tip-STATION nodes (SectionalExtractor.ei_from_static_response).
Prints
  (1) tip compliance 2x2 [m/MN] (tip node and ring mean, S-2 load) vs the BeamDyn v5 driver p60
      2x2 (logs/check/bd_static/LC{2,3}-p60-s0.0001, shell frame), with reciprocity residual,
      eigenvalues and principal angle; tip-ring best-fit rotation about z;
  (2) per station (S-1 stations and load): the single-axis S-1 EIs (must reproduce the stored S-1
      rows) and the eigenvalue EIs of the symmetrised sectional compliance S = kappa / (F (L - z)),
      both against the deck principal EIs, mean signed and mean |.| over 0.2-0.85 R.

    scripts/aeroenv.sh python -B tools/diagnose_bending_compliance_2x2.py <element_size>
"""

import json
import sys
from pathlib import Path

import numpy as np
from scipy.interpolate import UnivariateSpline

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))
sys.path.insert(0, str(ROOT))

from beam_reference import load_beamdyn_blade  # noqa: E402
from aeroelast.core.bc import DirichletCondition, NodalLoad  # noqa: E402
from aeroelast.core.mesh.generators import BladeMesh  # noqa: E402
from aeroelast.elements import ElementFamily  # noqa: E402
from aeroelast.models.blade.model import build_rust_properties  # noqa: E402
from aeroelast.postprocess.sectional import SectionalExtractor  # noqa: E402
from aeroelast.solvers.elasticity.static_linear import StaticLinearSolver  # noqa: E402

F = 1.0e6
BD_DIR = ROOT / "logs" / "check" / "bd_static"


def solve(mesh, props, nodes, comp):
    """Linear static solve with F along global `comp`, split equally over `nodes`."""
    cfg = {"solver": {}, "elements": {"element_family": ElementFamily.SHELL,
                                      "span_direction": (0.0, 0.0, 1.0), "properties": props}}
    s = StaticLinearSolver(mesh, cfg)
    dpn = s.domain.dofs_per_node
    idx = mesh.node_id_to_index
    root = sorted(mesh.get_node_set("RootNodes").nodes.keys())
    s.add_dirichlet_conditions([DirichletCondition(sorted(idx[n] * dpn + d for n in root
                                                          for d in range(dpn)), 0.0)])
    s.add_nodal_loads([NodalLoad([int(i) * dpn + comp for i in nodes], [F / len(nodes)] * len(nodes))])
    return np.asarray(s.solve()).reshape(-1, dpn)


def bd_tip_compliance():
    """BeamDyn p60 linear tip 2x2 in the shell frame [m/MN]: column j = load along shell axis j."""
    c = np.zeros((2, 2))
    for case, j in (("LC3", 0), ("LC2", 1)):
        r = json.loads((BD_DIR / f"{case}-p60-s0.0001" / "result.json").read_text())
        s, t = r["scale"], r["tip_bd_m"]
        # BD (x,y,z) -> shell (-y,x,z): shell ux = -TipTDyr, uy = TipTDxr. The LC3 BD load is
        # -y (shell +x) and LC2 is BD +x (shell +y), both 1 MN * scale.
        c[:, j] = np.array([-t["TipTDyr"], t["TipTDxr"]]) / s
    return c


def describe(c):
    ev, vec = np.linalg.eigh(0.5 * (c + c.T))
    ang = np.degrees(np.arctan2(vec[1, 0], vec[0, 0]))
    ang = (ang + 90.0) % 180.0 - 90.0
    recip = abs(c[0, 1] - c[1, 0]) / np.abs(c).max()
    return (f"[[{c[0, 0]:9.4f} {c[0, 1]:8.4f}] [{c[1, 0]:8.4f} {c[1, 1]:9.4f}]]  recip {recip:.2%}  "
            f"eig {ev[0]:.4f} / {ev[1]:.4f}  stiff-axis angle vs x {ang:+.2f} deg")


def ring_rotation_z(xy, uxy):
    r = xy - xy.mean(axis=0)
    return float(np.sum(r[:, 0] * uxy[:, 1] - r[:, 1] * uxy[:, 0]) / np.sum(r ** 2))


def main():
    es = float(sys.argv[1])
    gen = BladeMesh(yaml_file=str(ROOT / "tests/IEA-15-240-RWT.yaml"), airfoil_dir=str(ROOT / "tests/airfoils"),
                    element_size=es, skin_reference="mid")
    mesh = gen.generate(renumber="rcm")
    props = build_rust_properties(gen.numad_mesh_data)
    coords = mesh.coords_array
    print(f"es={es} nodes={len(coords)} skin=mid")

    # (1) tip compliance, S-2 load
    ring = np.nonzero(coords[:, 2] >= coords[:, 2].max() - 1e-6)[0]
    tip = int(np.argmax(coords[:, 2]))
    u2 = {comp: solve(mesh, props, ring, comp) for comp in (0, 1)}
    c_node = np.column_stack([u2[j][tip, :2] for j in (0, 1)]) / (F / 1e6)
    c_ring = np.column_stack([u2[j][ring, :2].mean(axis=0) for j in (0, 1)]) / (F / 1e6)
    c_bd = bd_tip_compliance()
    print(f"TIP shell node  {describe(c_node)}")
    print(f"TIP shell ring  {describe(c_ring)}")
    print(f"TIP BeamDyn p60 {describe(c_bd)}")
    for j, name in ((0, "LC3 +x"), (1, "LC2 +y")):
        print(f"TIP ring {name}: ux {c_ring[0, j]:.4f} uy {c_ring[1, j]:.4f} uz "
              f"{u2[j][ring, 2].mean():.4f} m, rot_z {np.degrees(ring_rotation_z(coords[ring, :2], u2[j][ring, :2])):+.4f} deg")
    print(f"TIP ring/BD - 1: xx {c_ring[0, 0] / c_bd[0, 0] - 1:+.2%}  yy {c_ring[1, 1] / c_bd[1, 1] - 1:+.2%}  "
          f"eig-min {np.linalg.eigvalsh(c_ring)[0] / np.linalg.eigvalsh(c_bd)[0] - 1:+.2%}  "
          f"eig-max {np.linalg.eigvalsh(c_ring)[1] / np.linalg.eigvalsh(c_bd)[1] - 1:+.2%}")

    # (2) sectional 2x2, S-1 load and stations
    ext = SectionalExtractor(mesh, props)
    mem = ext.station_membership()
    sts = sorted(mem)
    zs = np.array([mem[st]["z"] for st in sts])
    L = zs.max()
    u1 = {comp: solve(mesh, props, mem[sts[-1]]["nodes"], comp) for comp in (0, 1)}
    curv = np.zeros((len(sts), 2, 2))  # [station, displacement comp, load comp]
    for j in (0, 1):
        for i in (0, 1):
            prof = np.array([u1[j][mem[st]["nodes"], i].mean() for st in sts])
            curv[:, i, j] = UnivariateSpline(zs, prof, s=1e-2).derivative(2)(zs)
    bd = load_beamdyn_blade(ROOT / "tests/IEA15MW/reference/IEA-15-240-RWT_BeamDyn_blade.dat")
    ref_f, ref_e = np.interp(zs, bd.z, bd.ei_flap_ec), np.interp(zs, bd.z, bd.ei_edge_ec)
    mid = (zs > 0.2 * L) & (zs < 0.85 * L)
    rows = []
    for k in np.nonzero(mid)[0]:
        m = F * (L - zs[k])
        single_e, single_f = m / abs(curv[k, 0, 0]), m / abs(curv[k, 1, 1])
        s = curv[k] / m
        ev, vec = np.linalg.eigh(0.5 * (s + s.T))
        ei_lo, ei_hi = 1.0 / ev[1], 1.0 / ev[0]  # larger compliance = flap (smaller EI)
        ang = (np.degrees(np.arctan2(vec[1, 0], vec[0, 0])) + 90.0) % 180.0 - 90.0
        rows.append((zs[k], single_f / ref_f[k] - 1, single_e / ref_e[k] - 1, ei_lo / ref_f[k] - 1,
                     ei_hi / ref_e[k] - 1, abs(s[0, 1] - s[1, 0]) / np.abs(s).max(), ang))
    r = np.array(rows)
    print(f"SECTION n={len(r)} (0.2-0.85 R)          flap signed / |.|     edge signed / |.|")
    for label, cf, ce in (("single-axis (S-1 method)", 1, 2), ("2x2 eigenvalue", 3, 4)):
        print(f"  {label:26s} {100 * np.mean(r[:, cf]):+7.2f} / {100 * np.mean(np.abs(r[:, cf])):5.2f} %   "
              f"{100 * np.mean(r[:, ce]):+7.2f} / {100 * np.mean(np.abs(r[:, ce])):5.2f} %")
    print(f"  sectional reciprocity residual median {np.median(r[:, 5]):.2%} max {np.max(r[:, 5]):.2%}")
    print("  z [m]   single f/e [%]     2x2 f/e [%]    stiff (edge) axis angle vs x [deg]")
    for z, sf, se, ef, ee, _, ang in r:
        print(f"  {z:6.2f} {100 * sf:+8.2f} {100 * se:+8.2f}   {100 * ef:+8.2f} {100 * ee:+8.2f}   {ang:+7.2f}")


if __name__ == "__main__":
    main()
