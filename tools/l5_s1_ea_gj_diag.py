"""T14b: is the S-1 v2 GJ / EA(CLT) deficit physics or extraction artifact?

  scripts/aeroenv.sh python tools/l5_s1_ea_gj_diag.py extract OUT.json [--strips 0.25 0.5 ...]   (V2_BLADE_* env as the S-1 test)
  scripts/aeroenv.sh python tools/l5_s1_ea_gj_diag.py indep   OUT.json                           (no solves, ~1 min)
  python tools/l5_s1_ea_gj_diag.py report OUT1.json [OUT2.json ...]

extract: the S-1 channels (compute_s1_v2) per z-strip width (GMSH_STRIP_M; windio2 ignores it).
indep:   per station z = const plane cut of the mid-surface mesh (raw builder arrays, no S-1 helper):
           EA   = sum_cut-segments A11_eff(layup) * L,   A11_eff = 1 / (A^-1)[0,0] of the CLT A matrix (own code),
                  L = length of the segment the plane cuts out of the element (skin + web);
           GJ_B = multi-cell Bredt 4 A^2 / oint ds/(G t), G t = A66 of the element layup (tools/l3_bredt_mesh.graph_bredt).
"""
import argparse, json, sys
from pathlib import Path
import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO)); sys.path.insert(0, str(REPO / "tools"))


def mid(z, n=None):
    return (z > 0.2 * z.max()) & (z < 0.85 * z.max())


def cmd_extract(a):
    import tests.support.v2_sectional as vs
    from tests.validation.blade.v2.test_v2_iea15mw_s1_sectional import compute_s1_v2, mid_errors
    out = {}
    for w in a.strips:
        vs.GMSH_STRIP_M = w
        r = compute_s1_v2()
        e = mid_errors(r)
        out[str(w)] = {"z": r["z"].tolist(), "ea": r["ea"].tolist(), "ea_bd": r["ea_bd"].tolist(),
                       "gj": r["gj"].tolist(), "gj_bd": r["gj_bd"].tolist(),
                       "ei_edge": r["ei_edge"].tolist(), "ei_edge_bd": r["ei_edge_bd"].tolist(),
                       "err": {k: v for k, v in e.items()}}
        print(w, {k: (round(100 * v[0], 2), round(100 * v[1], 2)) for k, v in e.items()}, flush=True)
        Path(a.out).write_text(json.dumps(out))


def a_matrix(layup, mats):
    A = np.zeros((3, 3))
    for name, t, ang in layup:
        m = mats[name]
        f3 = lambda v: list(v) if isinstance(v, (list, tuple)) else [v] * 3
        E1, E2 = f3(m["E"])[:2]; n12 = f3(m["nu"])[0]; G = f3(m["G"])[0] or E1 / (2 * (1 + n12))
        d = 1 - n12 * n12 * E2 / E1
        Q = np.array([[E1 / d, n12 * E2 / d, 0], [n12 * E2 / d, E2 / d, 0], [0, 0, G]])
        c, s = np.cos(np.radians(ang)), np.sin(np.radians(ang))
        T = np.array([[c * c, s * s, 2 * c * s], [s * s, c * c, -2 * c * s], [-c * s, c * s, c * c - s * s]])
        R = np.diag([1, 1, 2.0])  # Reuter: engineering shear strain
        A += (np.linalg.inv(T) @ Q @ R @ T @ np.linalg.inv(R)) * t
    return A


def cmd_indep(a):
    from tests.support.v2_blade import build_v2_blade
    from tests.support import v2_blade as vb
    from tests.support.v2_sectional import bay_centres_v2
    from l3_bredt_mesh import graph_bredt
    from beam_reference import load_beamdyn_blade
    from l3_bredt_yaml import load_beamdyn_6x6, condense_shear
    bdp = REPO / "tests/IEA15MW/reference/IEA-15-240-RWT_BeamDyn_blade.dat"
    bd = load_beamdyn_blade(bdp)
    etas, ks = load_beamdyn_6x6(bdp)
    mesh, props = build_v2_blade(0.25, skin_reference="mid")
    raw = vb._MEMO[mesh.v2_info["cache_key"]]
    tab = json.loads(bytes(raw["tables_json"]).decode())
    mats = tab["materials"]
    lay = tab["layups"]
    ea_l = np.array([1 / np.linalg.inv(a_matrix(L, mats))[0, 0] for L in lay])
    # shear stiffness G t: A66 of the laminate (own Qbar66 via the transformed Q)
    g_l = np.array([a_matrix(L, mats)[2, 2] for L in lay])
    x = vb._inset(raw)
    el, lid = raw["elements"], raw["layup_id"]
    zs = bay_centres_v2()
    zmin = np.array([x[e[e >= 0], 2].min() for e in el]); zmax = np.array([x[e[e >= 0], 2].max() for e in el])
    rows = []
    for zc in zs:
        ea = 0.0; segs = []; gts = []; skipped = 0
        for k in np.nonzero((zmin < zc) & (zmax > zc))[0]:
            p = x[el[k][el[k] >= 0]]; hits = []
            for u, v in zip(p, np.roll(p, -1, axis=0)):
                if (u[2] <= zc) != (v[2] <= zc):
                    hits.append(u + (v - u) * (zc - u[2]) / (v[2] - u[2]))
            if len(hits) != 2:
                skipped += 1; continue
            L = np.linalg.norm(hits[0] - hits[1])
            ea += ea_l[lid[k]] * L
            segs.append((hits[0][:2], hits[1][:2])); gts.append(g_l[lid[k]])
        pts = np.array([s for sg in segs for s in sg])
        keys, inv = np.unique(np.round(pts / 1e-6).astype(np.int64), axis=0, return_inverse=True)
        inv = inv.ravel()
        edges, gt = [], []
        for i, g in enumerate(gts):
            u, v = int(inv[2 * i]), int(inv[2 * i + 1])
            if u != v:
                edges.append((u, v)); gt.append(g)
        try:
            b = graph_bredt(keys * 1e-6, edges, gt); gjb, nc, nopen = b["gj"], b["n_cells"], b["n_open"]
        except Exception as ex:
            gjb, nc, nopen = float("nan"), 0, -1
        k6 = np.array([[np.interp(zc / bd.z.max(), etas, ks[:, r, c]) for c in range(6)] for r in range(6)])
        rows.append({"z": float(zc), "ea": ea, "ea_bd": float(np.interp(zc, bd.z, bd.ea)), "k33": float(k6[2, 2]),
                     "gj_bredt": gjb, "gj_bd": float(np.interp(zc, bd.z, bd.gj)), "gj_sc": condense_shear(k6),
                     "n_cells": nc, "n_open": nopen, "n_cut": len(segs), "skipped": skipped})
        print(rows[-1], flush=True)
    Path(a.out).write_text(json.dumps(rows))


def cmd_clt(a):
    """No-solve replica of the CLT EA loop of tools/run_s1_sectional.ea_gj_from_static_response (the part that fails on
    gmsh), with its three ingredients split: elements counted (all nodes in the station's node set vs assigned to the
    bay), and the dz divisor (node-mean z difference vs nominal strip pitch)."""
    import tests.support.v2_sectional as vs
    from tests.support.v2_blade import build_v2_blade
    from tests.support.v2_sectional import V2SectionalExtractor, ring_ids_for_v2_mesh, bay_centres_v2
    from run_s1_sectional import laminate_a11_eff
    mesh, props = build_v2_blade(0.25, skin_reference="mid")
    out = {}
    mats = {m["name"]: m for m in mesh.numad_mesh_data["materials"]}
    secs = {s["elementSet"]: s for s in mesh.numad_mesh_data["sections"]}
    e2s = {}
    for en, es in mesh.element_sets.items():
        for e in es.elements:
            e2s.setdefault(e.id, en)
    a11 = {i: laminate_a11_eff(secs[n]["layup"], mats) for i, n in e2s.items() if n in secs}
    idx = mesh.node_id_to_index
    coords = np.asarray([n.coords for n in mesh.nodes])
    area = np.zeros(len(mesh.elements)); nodes = []
    for i, e in enumerate(mesh.elements):
        ii = [idx[n] for n in e.node_ids]; nodes.append(ii); c = coords[ii]
        area[i] = 0.5 * np.linalg.norm(np.cross(c[2] - c[0], c[3] - c[1])) if len(c) == 4 else 0.5 * np.linalg.norm(np.cross(c[1] - c[0], c[2] - c[0]))
    for w in a.strips:
        vs.GMSH_STRIP_M = w
        ring = ring_ids_for_v2_mesh(mesh); ext = V2SectionalExtractor(mesh, props, ring_ids=ring)
        mem = ext.station_membership()
        # element -> bay actually assigned by station_membership
        assigned = {}
        for st, m in mem.items():
            pass
        from aeroelast.postprocess.sectional import select_stations
        sts = sorted(mem); sel = [sts[i] for i in select_stations([mem[s]["z"] for s in sts], bay_centres_v2())]
        rows = []
        for st in sel:
            ns = set(mem[st]["nodes"].tolist()); dz = max(mem[st]["z"] - (mem[st - 1]["z"] if st - 1 in mem else mem[st]["z"]), 1e-9)
            tot = sum(a11[e.id] * area[i] for i, e in enumerate(mesh.elements) if e.id in a11 and all(n in ns for n in nodes[i]))
            cnt = sum(1 for i, e in enumerate(mesh.elements) if e.id in a11 and all(n in ns for n in nodes[i]))
            rows.append({"z": mem[st]["z"], "ea_helper": tot / dz, "dz": dz, "n_el": cnt, "sum_a11A": tot})
        out[str(w)] = rows
        print(w, "median dz", np.median([r["dz"] for r in rows]), "median n_el", np.median([r["n_el"] for r in rows]), flush=True)
    Path(a.out).write_text(json.dumps(out))


def cmd_report(a):
    for f in a.files:
        d = json.loads(Path(f).read_text())
        print("==", f)
        if isinstance(d, list):
            z = np.array([r["z"] for r in d]); m = mid(z)
            for q, ref in (("ea", "ea_bd"), ("gj_bredt", "gj_bd"), ("gj_bredt", "gj_sc")):
                r = np.array([x[q] / x[ref] - 1 for x in d])
                print(f"indep {q} vs {ref}: mid mean |rel| {np.nanmean(abs(r[m]))*100:.2f}% signed {np.nanmean(r[m])*100:+.2f}%")
        else:
            for w, v in d.items():
                print(w, {k: tuple(round(100 * t, 2) for t in e) for k, e in v["err"].items()})


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); sp = ap.add_subparsers(dest="c", required=True)
    e = sp.add_parser("extract"); e.add_argument("out"); e.add_argument("--strips", type=float, nargs="+", default=[0.25]); e.set_defaults(f=cmd_extract)
    i = sp.add_parser("indep"); i.add_argument("out"); i.set_defaults(f=cmd_indep)
    c = sp.add_parser("clt"); c.add_argument("out"); c.add_argument("--strips", type=float, nargs="+", default=[0.25]); c.set_defaults(f=cmd_clt)
    r = sp.add_parser("report"); r.add_argument("files", nargs="+"); r.set_defaults(f=cmd_report)
    a = ap.parse_args(); a.f(a)
