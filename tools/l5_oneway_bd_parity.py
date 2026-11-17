"""L5 one-way parity: the shell and BeamDyn under the same aerodynamic station resultants.

``odd/tasks/l5-one-way.md`` T1 (``--loads ad``: the S-5 AeroDyn azimuth-mean blade-1 loads) and T3
(``--loads prod``: the production BEM participant, built as ``run_bem_fsi.main`` builds it from
``cases/l5_oneway_rated/fluid/fluid.yaml``, on the rigid blade; ``--loads bem``: the pre-T18 T3 load, kept
for the earlier records).  Design: L3 of that document; the T3 load definition: L37 and L39.

1. Shell: reference gmsh chord-scaled N20 AR 1 mid-surface mesh; rotating linear stiffness
   K + K_G + K_SP at 7.56 rpm about mesh +y through (0, 0, -hub_radius); u = K_eff^-1 f.
2. Same loads: the shell's nodal force vector f is binned into 0.5 m span strips; each strip's
   resultant force and its moment about the BeamDyn reference axis (key-point curve, mapped to the
   mesh frame) become one BeamDyn point load.  Both models therefore carry identical station
   resultants whatever the shell's nodal distribution.
3. BeamDyn: ``beamdyn_driver`` rotating at the same speed (GlbPos z = hub radius), QuasiStaticInit,
   first row; delta = run(load * s) - run(no load), s = 1 (``nl``) and s = 0.01 rescaled (``lin``).
4. Compared along span: deflection at the reference axis (flap = mesh y, edge = mesh x) and the
   elastic twist about +z, nose-up positive (BeamDyn phi_z = 4 atan(c_z / 4), R1 of
   ``tools/twist_p1.py``; shell: production equal-weight ring estimator and the skin-weighted one).

Gravity is off and the structure has no precone on both sides.  Output: one JSON.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "tools"))

import twist_p1 as tp  # noqa: E402

RPM_RATED = 7.56
OMEGA = RPM_RATED * 2.0 * np.pi / 60.0
AXIS = np.array([0.0, 1.0, 0.0])  # rotor axis in the mesh frame (flap direction), as run_s6_oneway_dynamic
DECK_DIR = Path("/scratch/leahk/eduardo.donestevez/ofruns/OpenFAST/IEA-15-240-RWT")
BEAMDYN = os.environ.get("BEAMDYN_DRIVER", "/scratch/leahk/eduardo.donestevez/conda-envs/openfast/bin/beamdyn_driver")
AD_LOADS = REPO / "tests" / "IEA15MW" / "reference" / "s5_ad_blade1_loads.csv"
AD_PRIMARY = REPO / "tests" / "reference" / "iea15mw_openfast" / "case" / "IEA-15-240-RWT_AeroDyn15.dat"
STRIP_M = 0.5
FLUID_YAML = REPO / "cases/l5_oneway_rated/fluid/fluid.yaml"
N_AZ = 36
SPAN_FRACTIONS = (0.3, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 1.0)


def bem_config(tilt, precone):
    return {
        "wind_speed": 10.59, "omega": OMEGA, "pitch": 0.0, "azimuth": 0.0, "air_density": 1.225,
        "dynamic_viscosity": 1.81206e-5, "hub_height": 150.0, "shear_exp": 0.0, "tilt": tilt, "precone": precone,
    }


def build_shell(chord_divisions):
    from aeroelast.core.assembler import MeshAssembler
    from aeroelast.elements import ElementFamily
    from tests.support.v2_blade import build_v2_blade

    mesh, props = build_v2_blade(0.25, skin_reference="mid", mesher="gmsh", es_range=(0.05, 0.25),
                                 chord_divisions=chord_divisions)
    model = {"elements": {"element_family": ElementFamily.SHELL, "span_direction": (0.0, 0.0, 1.0), "properties": props}}
    return mesh, props, MeshAssembler(mesh, model)


def rotating_keff(asm, mesh, omega, hub_radius):
    """K + K_G + K_SP (csr) and the clamped-root free DOFs; ``run_s6_oneway_dynamic.assemble_keff`` with the hub offset."""
    from scipy.sparse import coo_matrix, csr_matrix

    n_dof, dpn = asm.dofs_count, asm.dofs_per_node
    r, c, v = asm._rust.assemble_k()
    K = coo_matrix((np.asarray(v), (np.asarray(r), np.asarray(c))), shape=(n_dof, n_dof)).tocsr()
    root = {mesh.node_id_to_index[n] for n in mesh.get_node_set("RootNodes").node_ids}
    free = np.array(sorted(set(range(n_dof)) - {i * dpn + d for i in root for d in range(dpn)}), dtype=np.int64)
    kg = asm.assemble_geometric_stiffness(omega=omega, rotation_axis=AXIS,
                                          rotation_center=np.array([0.0, 0.0, -hub_radius]), free_dofs=free)
    ai, aj, av = kg.getValuesCSR()
    Kg = csr_matrix((av, aj, ai), shape=(n_dof, n_dof))
    mr, _, mv = asm._rust.assemble_m()
    nodal_m = np.zeros(mesh.node_count)
    for row, val in zip(np.asarray(mr), np.asarray(mv), strict=True):
        if row % dpn < 3:
            nodal_m[row // dpn] += val / 3.0
    blk = -omega**2 * (np.eye(3) - np.outer(AXIS, AXIS))
    idx = np.arange(mesh.node_count) * dpn
    rows = np.concatenate([idx + a for a in range(3) for b in range(3)])
    cols = np.concatenate([idx + b for a in range(3) for b in range(3)])
    vals = np.concatenate([nodal_m * blk[a, b] for a in range(3) for b in range(3)])
    Ksp = coo_matrix((vals, (rows, cols)), shape=(n_dof, n_dof)).tocsr()
    return (K + Kg + Ksp).tocsr(), free


def production_cfg():
    """The L5 rated Fluid config with its blade and BEM files resolved to absolute paths."""
    import yaml

    cfg = yaml.safe_load(FLUID_YAML.read_text())
    for k in ("blade_file", "config_file"):
        cfg[k] = str((FLUID_YAML.parent / cfg[k]).resolve())
    return cfg


def production_loads(mesh, tilt, precone):
    """(n_nodes, 3) production aero load on the full shell mesh: the T6 load (L37), from the config alone.

    The Fluid participant is built exactly as ``run_bem_fsi.main`` builds it (node-set filtered coupling mesh,
    windIO polars, the runner's own wall-flow decision), evaluated on the rebuilt branch (max|u| > 1e-12, the branch
    production always takes) and averaged over ``N_AZ`` azimuths; its nodal forces are written onto the full mesh
    by coordinate match, so web nodes carry nothing, as in the coupled run.  ``tilt``/``precone`` must equal the
    config's values: the config, not the caller, defines this load.
    """
    from scipy.spatial import cKDTree

    from aeroelast.cli.run_bem_fsi import _build_mesh, _effective_element_properties
    from aeroelast.solvers.bem.fsi_participant import build_from_config

    cfg = production_cfg()
    if (float(cfg["bem"]["tilt"]), float(cfg["bem"]["precone"])) != (float(tilt), float(precone)):
        raise SystemExit(f"--loads prod: tilt/precone {tilt}/{precone} != config "
                         f"{cfg['bem']['tilt']}/{cfg['bem']['precone']} ({FLUID_YAML})")
    cmesh, viz, ep = _build_mesh(cfg, FLUID_YAML)
    prod = build_from_config(cmesh, cfg, viz_mesh=viz, element_properties=_effective_element_properties(cfg, ep))
    Xc, Xf = np.asarray(prod._ref_coords), np.asarray(mesh.coords_array)
    dist, idx = cKDTree(Xf).query(Xc)
    if dist.max() > 1e-6 or len(np.unique(idx)) != len(Xc):
        raise SystemExit(f"coupling mesh is not a vertex subset of the shell mesh: max dist {dist.max():.3e} m, "
                         f"{len(np.unique(idx))} targets for {len(Xc)} nodes")
    u = np.zeros((len(Xc), 3))
    u[:, 0] = 1e-6
    acc = np.zeros((len(Xc), 3))
    for a in np.radians(np.arange(N_AZ) * 360.0 / N_AZ):
        prod._partner_theta = float(a)
        acc += prod._compute_forces(u)[0] / N_AZ
    f3 = np.zeros((len(Xf), 3))
    f3[idx] = acc
    return f3


def shell_loads(kind, asm, mesh, props, tilt, precone):
    """Nodal force vector (n_dof,) and the participant (rings, production twist estimator).

    ``prod``: the production load (``production_loads``); ``bem``: the pre-T18 T3 load (AeroDyn polars, wall-flow
    realisation on, azimuth 0, every node); ``ad``: the S-5 AeroDyn azimuth-mean loads.  The returned participant is
    always the full-mesh one, so the rings and the twist estimator do not depend on ``kind``.
    """
    from aeroelast.solvers.bem.fsi_participant import BEMFSIParticipant
    from tests.support.openfast_bem import build_blade_aero_from_aerodyn

    part = BEMFSIParticipant(mesh, build_blade_aero_from_aerodyn(AD_PRIMARY), bem_config(tilt, precone),
                             element_properties=props)
    if kind == "prod":
        fa = production_loads(mesh, tilt, precone)
        f = np.zeros(asm.dofs_count)
        for d in range(3):
            f[d::asm.dofs_per_node] = fa[:, d]
    elif kind == "bem":
        fa, _ = part._compute_forces(np.zeros((mesh.node_count, 3)))
        f = np.zeros(asm.dofs_count)
        for d in range(3):
            f[d::asm.dofs_per_node] = fa[:, d]
    else:
        from run_s5_oneway import apply_ad_loads

        r_m, fn, ft, twist = np.loadtxt(AD_LOADS, delimiter=",", skiprows=1).T
        f = apply_ad_loads(asm, mesh, r_m, fn, ft, twist)
    return f, part


def ref_axis_mesh(kp_bd, z):
    """Reference-axis point(s) in the mesh frame at span z (key points are BeamDyn-frame, v_mesh = v_bd @ MESH_TO_BD)."""
    kp = kp_bd @ tp.MESH_TO_BD
    return np.column_stack([np.interp(z, kp[:, 2], kp[:, a]) for a in range(3)])


def bin_rows(X, f3, kp_bd):
    """BeamDyn point-load rows (eta F M, BD frame) from the nodal forces, and their mesh-frame resultants."""
    k = np.floor(X[:, 2] / STRIP_M).astype(int)
    arc = np.concatenate([[0.0], np.cumsum(np.linalg.norm(np.diff(kp_bd, axis=0), axis=1))])
    out, F_all, M_all, z_all = [], [], [], []
    for b in np.unique(k):
        m = k == b
        if not np.any(f3[m]):
            continue
        z = float(X[m, 2].mean())
        r0 = ref_axis_mesh(kp_bd, z)[0]
        F = f3[m].sum(0)
        M = np.cross(X[m] - r0, f3[m]).sum(0)
        eta = float(np.interp(z, kp_bd[:, 2], arc / arc[-1]))
        out.append(np.concatenate([[eta], F @ tp.MESH_TO_BD.T, M @ tp.MESH_TO_BD.T]))
        F_all.append(F), M_all.append(M), z_all.append(z)
    return np.array(out), np.array(F_all), np.array(M_all), np.array(z_all)


def shell_sections(mesh, part, U, kp_bd):
    """Per production ring: z, deflection at the reference axis (affine fit), twist (equal and skin weights)."""
    from aeroelast.solvers.bem.section_contour import ring_section_rotation, ring_tributary_lengths
    from tests.validation.blade.test_blade_section_torsion_stiffness import _skin_wall_edges

    X = part._ref_coords
    edges = _skin_wall_edges(mesh)
    rows = []
    for ring in part._strip_ring_indices:
        ring = np.asarray(ring)
        if len(ring) < 4:
            continue
        z = float(X[ring, 2].mean())
        c = X[ring, :2].mean(0)
        A = np.column_stack([X[ring, :2] - c, np.ones(len(ring))])
        coef, *_ = np.linalg.lstsq(A, U[ring, :2], rcond=None)
        r0 = ref_axis_mesh(kp_bd, z)[0]
        u_axis = np.concatenate([np.array([*(r0[:2] - c), 1.0]) @ coef, [U[ring, 2].mean()]])
        eq = part._section_rotation(U, ring)
        sk = ring_section_rotation(X[ring], U[ring], np.array([0.0, 0.0, 1.0]),
                                   ring_tributary_lengths(ring, edges, X))
        rows.append([z, *u_axis, np.nan if eq is None else eq, np.nan if sk is None else sk])
    a = np.array(sorted(rows))
    return {"z": a[:, 0], "u": a[:, 1:4], "twist_eq": a[:, 4], "twist_skin": a[:, 5]}


def run_beamdyn(work, rows, hub_radius, scales, omega=OMEGA):
    import diagnose_radii_beamdyn as dmr

    tp.HUB_RADIUS = hub_radius
    ctx = {"work": work, "deck_dir": DECK_DIR, "beamdyn": BEAMDYN, "dmr": dmr, "hashes": {}}

    def state(name, r):
        deck, names, data = tp.run_case(ctx, name, dynamic=True, omega=omega, rows=r)
        td, rd = tp.node_fields(names, tp.pick_row(data, True), len(deck.eta))
        return deck, td @ tp.MESH_TO_BD, tp.wm_z_angle(rd[:, 2])

    deck, td0, phi0 = state("noload", None)
    out = {"z": deck.qp[:, 2]}
    for s in scales:
        r = rows.copy()
        r[:, 1:] *= s
        _, td, phi = state(f"load_x{s:g}", r)
        out[s] = {"u": (td - td0) / s, "twist": (phi - phi0) / s}
    return out, ctx["hashes"]


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--loads", choices=("ad", "bem", "prod"), required=True)
    ap.add_argument("--out", type=Path, required=True, help="result JSON")
    ap.add_argument("--work", type=Path, required=True, help="BeamDyn case directory")
    ap.add_argument("--chord-divisions", type=float, default=20.0)
    ap.add_argument("--hub-radius", type=float, default=tp.HUB_RADIUS)
    ap.add_argument("--tilt", type=float, default=-6.0, help="BEM shaft tilt [deg] (deck)")
    ap.add_argument("--precone", type=float, default=-4.0, help="BEM precone [deg] (deck)")
    ap.add_argument("--keep", default="xyz", help="diagnostic: mesh force components kept (e.g. x = edgewise only)")
    ap.add_argument("--rpm", type=float, default=RPM_RATED, help="structural rotor speed; 0 = parked diagnostic")
    ap.add_argument("--scales", type=float, nargs="+", default=[1.0, 0.01])
    args = ap.parse_args()
    t0 = time.time()

    from run_s5_oneway import solve_static  # PETSc LU (MUMPS) on the free DOFs

    mesh, props, asm = build_shell(args.chord_divisions)
    dpn = asm.dofs_per_node
    f, part = shell_loads(args.loads, asm, mesh, props, args.tilt, args.precone)
    for d, a in enumerate("xyz"):
        if a not in args.keep:
            f[d::dpn] = 0.0
    omega = args.rpm * 2.0 * np.pi / 60.0
    K, free = rotating_keff(asm, mesh, omega, args.hub_radius)
    u = solve_static(asm, K, free, f)
    U = np.column_stack([u[d::dpn] for d in range(3)])
    print(f"shell solved ({mesh.node_count} nodes, {time.time() - t0:.0f} s)", flush=True)

    kp_bd = tp.read_keypoints((DECK_DIR / "IEA-15-240-RWT_BeamDyn.dat").read_text())[:, :3]
    X = np.asarray(part._ref_coords, float)
    f3 = np.column_stack([f[d::dpn] for d in range(3)])
    rows, F, M, zs = bin_rows(X, f3, kp_bd)
    root = ref_axis_mesh(kp_bd, np.array([0.0]))[0]
    conservation = {
        "F_nodal": f3.sum(0).tolist(), "F_rows": F.sum(0).tolist(),
        "Mroot_nodal": np.cross(X - root, f3).sum(0).tolist(),
        "Mroot_rows": (np.cross(ref_axis_mesh(kp_bd, zs) - root, F) + M).sum(0).tolist(),
    }
    sh = shell_sections(mesh, part, U, kp_bd)
    bd, hashes = run_beamdyn(args.work, rows, args.hub_radius, args.scales, omega)
    print(f"beamdyn done ({time.time() - t0:.0f} s)", flush=True)

    L = float(kp_bd[-1, 2])
    zq = np.array([fr * L for fr in SPAN_FRACTIONS])
    zq[-1] = min(zq[-1], sh["z"].max(), bd["z"].max())

    def at(z, y):
        return np.interp(zq, z, y)

    table = {"z": zq.tolist(), "shell": {
        "flap": at(sh["z"], sh["u"][:, 1]).tolist(), "edge": at(sh["z"], sh["u"][:, 0]).tolist(),
        "twist_eq_deg": np.rad2deg(at(sh["z"], sh["twist_eq"])).tolist(),
        "twist_skin_deg": np.rad2deg(at(sh["z"], sh["twist_skin"])).tolist()}}
    for s in args.scales:
        b = bd[s]
        table[f"beamdyn_x{s:g}"] = {"flap": at(bd["z"], b["u"][:, 1]).tolist(), "edge": at(bd["z"], b["u"][:, 0]).tolist(),
                                    "twist_deg": np.rad2deg(at(bd["z"], b["twist"])).tolist()}
    W_shell = 0.5 * float(f @ u)
    doc = {"loads": args.loads, "chord_divisions": args.chord_divisions, "hub_radius": args.hub_radius,
           "tilt": args.tilt, "precone": args.precone, "rpm": args.rpm, "keep": args.keep, "n_nodes": mesh.node_count,
           "n_point_loads": len(rows), "W_shell": W_shell, "conservation": conservation, "table": table,
           "revision": tp.revision(REPO), "beamdyn": BEAMDYN, "generated_sha256": hashes,
           "elapsed_s": time.time() - t0}
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(doc, indent=1) + "\n")
    print(json.dumps(table, indent=1))


if __name__ == "__main__":
    main()
