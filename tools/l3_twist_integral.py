#!/usr/bin/env python
"""Derivative-free torsion check: FE twist increment vs the integral of T / GJ (odd/tasks/l3-closure.md, T9).

Decides whether the S-1 FE GJ deficit of a blade is a soft FE (element) or a wrong per-station derivative
(extraction). Mesh, tip shear-flow torque (T40 ``section_torque_wall_forces``), boundary conditions and the
rigid-fit twist theta(z) projected on the station tangent are those of the GJ branch of
``tools/run_s1_sectional.py:ea_gj_from_static_response``, but on EVERY mesh station and without any
derivative. Per band [z1, z2]: FE delta-theta = theta(z2) - theta(z1) against the integral of
T / (GJ tz) dz (the extractor takes the rate per arc length, ds = dz / tz) with GJ = the mesh Bredt of
``tools/l3_bredt_mesh.py`` on a dense z grid, and with the deck GJ_sc. T is the delivered moment about
the station tangent. Ratio ~1: the FE is torsionally right (extraction wrong); ~1.2-1.3: the FE is soft.

  python tools/l3_twist_integral.py --yaml tests/IEA-3.4-130-RWT.yaml \\
      --beamdyn tests/IEA34/reference/IEA-3.4-130-RWT_BeamDyn_blade.dat --blade-length 63.0 \\
      --element-size 0.15 --out logs/check/l3/t9/twist/iea34.json
  python tools/l3_twist_integral.py --selftest
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from l3_bredt_mesh import laminate_a66, mesh_section  # noqa: E402,F401
from l3_bredt_yaml import condense_shear, interp_deck, load_beamdyn_6x6  # noqa: E402

BANDS = ((0.3, 0.6), (0.6, 0.85), (0.7, 0.85), (0.4, 0.6))
SLICE_ETAS = tuple(round(0.05 * k, 2) for k in range(1, 20))


def band_integral(z, f, z1, z2) -> float:
    """Trapezoid integral of the samples (z, f) over [z1, z2], endpoints linearly interpolated."""
    z, f = np.asarray(z, float), np.asarray(f, float)
    zz = np.concatenate([[z1], z[(z > z1) & (z < z2)], [z2]])
    return float(np.trapezoid(np.interp(zz, z, f), zz))


def selftest() -> int:
    z = np.linspace(0.0, 10.0, 41)
    t, gj = 1.0e6, 4.0e8
    got = band_integral(z, np.full_like(z, t / gj), 1.3, 8.7)
    ref = t * (8.7 - 1.3) / gj
    print(f"selftest constant GJ: {got:.6e} vs T L / GJ {ref:.6e} ({100 * (got / ref - 1):+.2e} %)")
    ok = abs(got / ref - 1) < 1e-12
    print("selftest", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def fe_twist(gen, mesh, torque=1e6):
    """Per-station z, theta (rigid fit on the tangent), tangent, delivered moment about it; all stations."""
    from aeroelast.core.bc import DirichletCondition, NodalLoad
    from aeroelast.elements import ElementFamily
    from aeroelast.models.blade.model import build_rust_properties
    from aeroelast.postprocess.sectional import SectionalExtractor, rigid_rotation_fit
    from aeroelast.solvers.bem.force_projection import section_torque_wall_forces
    from aeroelast.solvers.elasticity.static_linear import StaticLinearSolver

    props = build_rust_properties(gen.numad_mesh_data)
    model = {"solver": {}, "elements": {"element_family": ElementFamily.SHELL,
                                        "span_direction": (0.0, 0.0, 1.0), "properties": props}}
    membership = SectionalExtractor(mesh, props).station_membership()
    sts = sorted(membership)
    zs = np.array([membership[st]["z"] for st in sts])
    tip_nodes = membership[sts[-1]]["nodes"]
    coords = np.asarray([n.coords for n in mesh.nodes])
    dpn = mesh.dofs_per_node if hasattr(mesh, "dofs_per_node") else 6
    idx = mesh.node_id_to_index
    root = sorted(mesh.get_node_set("RootNodes").nodes.keys())
    rootd = sorted(idx[n] * dpn + d for n in root for d in range(dpn))
    cents = np.array([coords[membership[st]["nodes"]].mean(axis=0) for st in sts])
    tangs = np.gradient(cents, zs, axis=0)
    tangs /= np.linalg.norm(tangs, axis=1)[:, None]

    prev_nodes = set(membership[sts[-2]]["nodes"].tolist())
    end_ring = np.array([n for n in tip_nodes.tolist() if n not in prev_nodes], dtype=np.intp)
    ring_xyz = coords[end_ring]
    ring_c = ring_xyz.mean(axis=0)
    normal = np.linalg.svd(ring_xyz - ring_c)[2][-1]
    if normal @ tangs[-1] < 0.0:
        normal = -normal
    ring, f_ring = section_torque_wall_forces(mesh, end_ring, normal, props, torque, coords=coords)
    m_delivered = np.cross(coords[ring] - ring_c, f_ring).sum(axis=0)

    s = StaticLinearSolver(mesh, model)
    s.add_dirichlet_conditions([DirichletCondition(rootd, 0.0)])
    for t, f in zip(ring, f_ring, strict=True):
        s.add_nodal_loads([NodalLoad([int(t) * dpn + d], [float(f[d])]) for d in range(3)])
    u = s.solve()
    theta = np.empty(len(sts))
    for i, st in enumerate(sts):
        nodes = membership[st]["nodes"]
        rot = rigid_rotation_fit(coords[nodes], coords[nodes].mean(axis=0), u.reshape(-1, dpn)[nodes, :3])
        theta[i] = float(rot @ tangs[i])
    return zs, theta, tangs, tangs @ m_delivered, m_delivered, len(end_ring)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--yaml")
    ap.add_argument("--beamdyn")
    ap.add_argument("--blade-length", type=float)
    ap.add_argument("--element-size", type=float, default=0.15)
    ap.add_argument("--torque", type=float, default=1e6)
    ap.add_argument("--bands", type=float, nargs="+", help="eta pairs z1 z2 [z1 z2 ...] (default: spec bands)")
    ap.add_argument("--bredt-step", type=float, default=0.5, help="[m] Bredt z grid spacing")
    ap.add_argument("--out")
    args = ap.parse_args(argv)
    if args.selftest:
        return selftest()
    if not all((args.yaml, args.beamdyn, args.blade_length, args.out)):
        ap.error("--yaml, --beamdyn, --blade-length and --out are required")
    if args.element_size >= 0.3:
        raise SystemExit("element_size >= 0.3 is not converged (docs/blade_mesh_size_rule.md)")
    bands = tuple(zip(args.bands[::2], args.bands[1::2], strict=True)) if args.bands else BANDS
    length = args.blade_length

    from aeroelast.core.mesh.generators import BladeMesh

    gen = BladeMesh(yaml_file=args.yaml, element_size=args.element_size, skin_reference="mid")
    mesh = gen.generate(renumber="rcm")
    data = gen.numad_mesh_data
    sections = {s["elementSet"]: s for s in data["sections"]}
    mats = {m["name"]: m for m in data["materials"]}
    elem_section = {el.id: name for name, eset in mesh.element_sets.items() if name in sections
                    for el in eset.elements}
    coords, idx = mesh.coords_array, mesh.node_id_to_index
    elements = []
    for el in mesh.elements:
        p = coords[[idx[n] for n in el.node_ids]]
        elements.append((elem_section[el.id], p, float(p[:, 2].mean())))
    rows = np.unique(np.round([zm for _, _, zm in elements], 6))

    zs, theta, tangs, m_t, m_del, n_ring = fe_twist(gen, mesh, args.torque)
    print(f"{Path(args.yaml).name} es {args.element_size}: {len(elements)} elements, {len(zs)} stations, "
          f"tip ring {n_ring} nodes; requested T {args.torque:.3e}; delivered T.t min {m_t.min():.4e} "
          f"max {m_t.max():.4e} (ratio {m_t.mean() / args.torque:.4f}); |M| {np.linalg.norm(m_del):.4e}")

    # Bredt z grid: every bredt-step over eta 0.25-0.95 plus the band and slice etas
    lo, hi = 0.25 * length, 0.95 * length
    targets = set(np.arange(lo, hi, args.bredt_step).tolist())
    targets |= {e * length for e in SLICE_ETAS if lo <= e * length <= hi}
    targets |= {e * length for b in bands for e in b}
    zb, gjb, secs = [], [], []
    for zt in sorted(targets):
        sec = mesh_section(gen, sections, mats, elements, rows, zt)
        if sec["n_open"] or zb and abs(sec["z_row"] - zb[-1]) < 1e-9:
            continue
        zb.append(sec["z_row"])
        gjb.append(sec["gj"])
        secs.append({"z": sec["z_row"], "gj": sec["gj"], "n_cells": sec["n_cells"]})
    zb, gjb = np.array(zb), np.array(gjb)
    etas, ks = load_beamdyn_6x6(args.beamdyn)
    gjsc = np.array([condense_shear(interp_deck(etas, ks, z / length)) for z in zb])
    tz = np.interp(zb, zs, tangs[:, 2])
    tb = np.interp(zb, zs, m_t)
    print(f"Bredt grid: {len(zb)} z in [{zb[0]:.2f}, {zb[-1]:.2f}] m, step {args.bredt_step} m")

    def compare(e1, e2):
        z1, z2 = e1 * length, e2 * length
        # band ends snap to the Bredt rows already (targets include them); use the actual row z
        i1, i2 = np.argmin(abs(zb - z1)), np.argmin(abs(zb - z2))
        z1, z2 = zb[i1], zb[i2]
        d_fe = float(np.interp(z2, zs, theta) - np.interp(z1, zs, theta))
        i_b = band_integral(zb, tb / (gjb * tz), z1, z2)
        i_d = band_integral(zb, tb / (gjsc * tz), z1, z2)
        return {"eta1": e1, "eta2": e2, "z1": z1, "z2": z2, "fe_dtheta": d_fe, "int_mesh_bredt": i_b,
                "int_deck_gjsc": i_d, "fe_over_mesh_bredt": d_fe / i_b, "fe_over_deck": d_fe / i_d}

    res = [compare(*b) for b in bands]
    sl = [compare(a, b) for a, b in zip(SLICE_ETAS[:-1], SLICE_ETAS[1:], strict=True)
          if lo <= a * length and b * length <= hi]
    print(f"{'eta band':>11} {'z1':>7} {'z2':>7} {'FE dtheta':>11} {'int Bredt':>11} {'int deck':>11} "
          f"{'FE/Bredt':>9} {'FE/deck':>8}")
    for r in res + sl:
        print(f"{r['eta1']:5.2f}-{r['eta2']:5.2f} {r['z1']:7.2f} {r['z2']:7.2f} {r['fe_dtheta']:11.4e} "
              f"{r['int_mesh_bredt']:11.4e} {r['int_deck_gjsc']:11.4e} {r['fe_over_mesh_bredt']:9.4f} "
              f"{r['fe_over_deck']:8.4f}" + ("   <- slice" if r in sl and r is sl[0] else ""))
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(json.dumps({
        "yaml": args.yaml, "beamdyn": args.beamdyn, "blade_length": length, "element_size": args.element_size,
        "torque_requested": args.torque, "torque_delivered_t": m_t.tolist(), "m_delivered": m_del.tolist(),
        "stations_z": zs.tolist(), "theta": theta.tolist(), "bredt": secs, "gj_sc": gjsc.tolist(),
        "bands": res, "slices": sl}, indent=1, default=float))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
