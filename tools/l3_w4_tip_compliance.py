"""W4/T5 (L3 D3): is the FE local GJ drop below Bredt from z~100 load introduction / section distortion (b) or mesh/extraction (c)?

usage: python tools/l3_w4_tip_compliance.py CASE OUT.json      CASE = gmsh<N>
Same case as the S-7 test (mid-surface, tip torque 1000 N.m as wall shear flow on TipNodes, root clamped).  Variant ``--zload Z``:
the same torque is applied on the ring of nodes nearest z = Z (full wall-flow, nothing outboard) - a load-introduction probe: if the drop
is tip introduction it moves with the load plane; if it is the section it stays at the same z.
Per station z: rings = node slabs; rigid-fit rotation (S-7 estimator) and affine-fit antisymmetric rotation; FE local GJ = T / theta' (central
difference over +-H m); mesh-cut Bredt GJ (tools/bredt_batho_bands.py machinery, same mid-surface cut); in-plane distortion = rms residual
of the rigid fit and of the affine fit of the slab, absolute and relative to rms(r)*theta.
"""

import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__))))
from tests.support.paths import REPO_ROOT  # noqa: E402
from tests.support.v2_blade import build_v2_blade, load_raw  # noqa: E402

sys.path.insert(0, str(REPO_ROOT / "tools"))
sys.path.insert(0, str(REPO_ROOT / "tests/validation/blade/v2"))

TORQUE = 1000.0
ZS = [90, 92, 94, 96, 98, 100, 102, 104, 106, 108, 110, 112, 114]
ZREF = [60, 70, 80]  # mid-span reference


def parse():
    a = sys.argv[1:]
    case, out = a[0], a[1]
    zload = float(a[a.index("--zload") + 1]) if "--zload" in a else None
    return case, out, zload


def mesh_kw(case):
    if not case.startswith("gmsh"):
        raise SystemExit(f"CASE {case!r}: only gmsh<N> exists (the windio2 mesher is retired, gmsh-official-mesher S9)")
    return {"mesher": "gmsh", "chord_divisions": int(case[4:])}


def solve(mesh, props, zload):
    from aeroelast.core.bc import DirichletCondition, NodalLoad
    from aeroelast.elements import ElementFamily
    from aeroelast.solvers.bem.force_projection import section_torque_wall_forces
    from aeroelast.solvers.elasticity.static_linear import StaticLinearSolver

    cfg = {
        "solver": {},
        "elements": {
            "element_family": ElementFamily.SHELL,
            "span_direction": (0.0, 0.0, 1.0),
            "properties": props,
        },
    }
    solver = StaticLinearSolver(mesh, cfg)
    dpn = solver.domain.dofs_per_node
    idx_of = mesh.node_id_to_index
    root_ids = sorted(mesh.get_node_set("RootNodes").nodes.keys())
    solver.add_dirichlet_conditions(
        [DirichletCondition(sorted(idx_of[n] * dpn + d for n in root_ids for d in range(dpn)), 0.0)]
    )
    coords = mesh.coords_array
    if zload is None:
        ring = np.asarray(
            sorted(idx_of[n] for n in mesh.get_node_set("TipNodes").nodes), dtype=np.intp
        )
        nodes, f = section_torque_wall_forces(
            mesh, ring, (0.0, 0.0, 1.0), props, TORQUE, coords=coords
        )
    else:  # one clean span station near zload: z levels = clusters with gaps > 5 mm; take the nearest cluster whose node count equals the
        # median of its 4 neighbours (a transition station, with a split or half-filled level, would load an incomplete ring)
        zl = coords[:, 2]
        order = np.argsort(zl)
        zs_sorted = zl[order]
        cut = np.nonzero(np.diff(zs_sorted) > 0.005)[0] + 1
        groups = np.split(order, cut)
        cnt = np.array([len(g) for g in groups])
        zmean = np.array([zl[g].mean() for g in groups])
        good = [
            i
            for i in range(2, len(groups) - 2)
            if cnt[i] == int(np.median(np.r_[cnt[i - 2 : i], cnt[i + 1 : i + 3]]))
            and np.ptp(zl[groups[i]]) < 0.08
        ]
        i = min(good, key=lambda k: abs(zmean[k] - zload))
        ring = np.sort(groups[i])
        nodes, f = section_torque_wall_forces(
            mesh, ring, (0.0, 0.0, 1.0), props, TORQUE, coords=coords
        )
        print(
            f"load ring: z={zmean[i]:.3f} nodes={len(ring)} z-range=[{zl[ring].min():.3f},{zl[ring].max():.3f}]",
            flush=True,
        )
    dofs, vals = [], []
    for t, ff in zip(nodes, f, strict=True):
        for d in range(3):
            dofs.append(int(t) * dpn + d)
            vals.append(float(ff[d]))
    solver.add_nodal_loads([NodalLoad(dofs, vals)])
    return np.asarray(solver.solve(), dtype=np.float64), len(ring)


def slab_stats(coords, u, dpn, zc, half):
    """Rigid fit with a linear-in-z rotation, u = t + (th0 + th1 dz) e_z x r, over the nodes of the slab, plus the affine antisymmetric rotation.

    The z-gradient term removes the twist variation inside the slab, so the residual is the in-plane distortion (not the twist increment)."""
    idx = np.nonzero(abs(coords[:, 2] - zc) <= half)[0]
    if len(idx) < 8:
        return None
    xy = coords[idx, :2]
    cen = xy.mean(axis=0)
    r = xy - cen
    dz = coords[idx, 2] - zc
    ux = u[idx * dpn + 0]
    uy = u[idx * dpn + 1]
    n = len(idx)
    A = np.zeros((2 * n, 4))
    A[:n, 0] = 1.0
    A[n:, 1] = 1.0
    A[:n, 2] = -r[:, 1]
    A[n:, 2] = r[:, 0]
    A[:n, 3] = -r[:, 1] * dz
    A[n:, 3] = r[:, 0] * dz
    b = np.concatenate([ux, uy])
    c = np.linalg.lstsq(A, b, rcond=None)[0]
    res_r = b - A @ c
    D = np.column_stack([r, np.ones(n)])
    cx = np.linalg.lstsq(D, ux, rcond=None)[0]
    cy = np.linalg.lstsq(D, uy, rcond=None)[0]
    rrms = float(np.sqrt((r**2).sum(axis=1).mean()))
    return {
        "n": n,
        "z": float(coords[idx, 2].mean()),
        "th": float(c[2]),
        "th_aff": float(0.5 * (cy[0] - cx[1])),
        "rrms": rrms,
        "res_rigid": float(np.sqrt((res_r**2).sum() / n)),
    }


def main():
    case, out, zload = parse()
    t0 = time.time()
    kw = mesh_kw(case)
    mesh, props = build_v2_blade(0.25, skin_reference="mid", **kw)
    u, nring = solve(mesh, props, zload)
    t_solve = time.time() - t0
    coords = mesh.coords_array
    dpn = 6
    zmax = float(coords[:, 2].max())
    # Bredt on the same cut (mid-surface inset of the raw mesh of the case)
    from l3_w4_bredt import make

    raw_kw = (
        {"chord_divisions": kw["chord_divisions"]}
        if "chord_divisions" in kw
        else {"es_range": (0.05, 0.25)}
    )
    raw = load_raw(0.25, mesher=kw["mesher"], generate=False, **raw_kw)
    gj_at, eta_deck, K66_deck = make(
        raw, str(REPO_ROOT / "tests/IEA15MW/reference/IEA-15-240-RWT_BeamDyn_blade.dat")
    )
    K66_deck, eta_deck = np.asarray(K66_deck), np.asarray(eta_deck)

    # slabs every 0.5 m (half width 0.25 m); local GJ = T / (d theta/dz), theta'(z) from a line fit over the slabs within +-W
    zhi = zmax if zload is None else zload
    sl = [
        s
        for s in (slab_stats(coords, u, dpn, zc, 0.25) for zc in np.arange(58.0, zmax, 0.5))
        if s is not None
    ]
    zsl = np.array([s["z"] for s in sl])
    ths = {k: np.array([s[k] for s in sl]) for k in ("th", "th_aff")}
    rows = []
    zlist = ZREF + ZS if zload is None else ZREF + list(range(90, int(zload) + 1, 2))
    for z in sorted(set(zlist) | set(range(60, int(zhi) + 1, 2))):
        row = {"z": float(z)}
        k0 = int(np.argmin(abs(zsl - z)))
        sc = sl[k0]
        gb, ncell, _ = gj_at(z + 0.0137)
        kd = float(np.interp(z / 117.0, eta_deck, K66_deck))
        row.update(
            nodes=sc["n"],
            gj_bredt=gb,
            cells=ncell,
            k66=kd,
            bredt_over_k66=gb / kd,
            theta_deg=float(np.degrees(sc["th"])),
        )
        for W in (1.5, 3.0):
            m = (abs(zsl - z) <= W) & (zsl <= zhi - 0.25)
            if m.sum() < 4 or z + W > zhi + 1e-9:
                continue
            for k in ("th", "th_aff"):
                slope = np.polyfit(zsl[m], ths[k][m], 1)[0]
                row[f"gj_{k}_W{W}"] = TORQUE / slope
                row[f"fe_over_bredt_{k}_W{W}"] = TORQUE / slope / gb
            slope = np.polyfit(zsl[m], ths["th"][m], 1)[0]
            row[f"dist_per_m_W{W}"] = float(
                np.mean(
                    [
                        s["res_rigid"] / (s["rrms"] * abs(slope))
                        for s, mm in zip(sl, m, strict=True)
                        if mm
                    ]
                )
            )
        row["dist_rel_theta"] = (
            sc["res_rigid"] / (sc["rrms"] * abs(sc["th"])) if sc["th"] else float("nan")
        )
        row["dist_abs_m"] = sc["res_rigid"]
        rows.append(row)

    # band compliance: (theta(z1)-theta(z0))/T against the integral of 1/GJ_Bredt over the band (Bredt every 1 m, trapezoid);
    # theta(z) from the +-1.5 m line fit of the slab rotations evaluated at the band ends.
    def theta_at(z, key="th", W=1.5):
        m = (abs(zsl - z) <= W) & (zsl <= zhi - 0.25)
        pf = np.polyfit(zsl[m], ths[key][m], 1)
        return float(np.polyval(pf, z))

    zb = np.arange(88.0, zhi + 1e-9, 1.0)
    gbz = np.array([gj_at(z + 0.0137)[0] for z in zb])
    bands = {}
    for z0, z1 in (
        (60, 70),
        (70, 80),
        (80, 90),
        (90, 96),
        (94, 100),
        (96, 102),
        (100, 106),
        (102, 108),
        (106, 112),
        (108, 114),
        (99, 105),
        (102, 105),
        (98, 102),
        (100, 104),
        (104, 108),
    ):
        if z1 > zhi - 1.0 or z0 < zb[0] and z0 < 88 and False:
            continue
        zz = np.arange(max(z0, 60), z1 + 1e-9, 1.0)
        if z0 < 88:
            gz = np.array([gj_at(z + 0.0137)[0] for z in zz])
        else:
            gz = np.interp(zz, zb, gbz)
        comp_bredt = float(np.trapezoid(1.0 / gz, zz))
        d = {}
        for k in ("th", "th_aff"):
            d[f"fe_over_bredt_{k}"] = float(
                comp_bredt / ((theta_at(z1, k) - theta_at(z0, k)) / TORQUE)
            )
        bands[f"{z0}-{z1}"] = d
    res = {
        "case": case,
        "zload": zload,
        "n_nodes": int(len(coords)),
        "zmax": zmax,
        "load_ring_nodes": int(nring),
        "torque": TORQUE,
        "solve_s": round(t_solve, 1),
        "total_s": round(time.time() - t0, 1),
        "rows": rows,
        "bands": bands,
    }
    json.dump(res, open(out, "w"), indent=1)
    print(
        " z    FE/Bredt(W1.5) (W3.0)  aff(W1.5)  Bredt/K66 cells  dist/(rrms*th')[W1.5]  dist_rel_th  nodes"
    )
    for r in rows:
        if "fe_over_bredt_th_W1.5" not in r:
            continue
        print(
            "%5.1f  %6.3f  %6.3f  %6.3f  %6.3f  %d  %9.3f  %9.4f  %d"
            % (
                r["z"],
                r["fe_over_bredt_th_W1.5"],
                r.get("fe_over_bredt_th_W3.0", np.nan),
                r["fe_over_bredt_th_aff_W1.5"],
                r["bredt_over_k66"],
                r["cells"],
                r["dist_per_m_W1.5"],
                r["dist_rel_theta"],
                r["nodes"],
            )
        )
    print(
        "band FE/Bredt (compliance ratio) rigid | affine:",
        {
            k: (round(v["fe_over_bredt_th"], 3), round(v["fe_over_bredt_th_aff"], 3))
            for k, v in bands.items()
        },
    )
    print("zmax", zmax, "load_ring", nring, "solve_s", res["solve_s"], "total_s", res["total_s"])


if __name__ == "__main__":
    main()
