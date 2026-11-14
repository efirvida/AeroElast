"""T4 (odd/tasks/gmsh-chord-scaled-mesh.md): mesh quality of the uniform vs chord-scaled gmsh IEA-15 meshes. No FEM solve.

    scripts/aeroenv.sh python logs/chord_scaled/t4/analyze.py run <label>     # label: uniform | n10 | n20  (writes <label>.json, <label>.vtu)
    scripts/aeroenv.sh python logs/chord_scaled/t4/analyze.py merge            # summary.json from the three
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
sys.path.insert(0, str(REPO / "tests/support"))
YAML = REPO / "tests/IEA15MW/reference/upstream/IEA-15-240-RWT.yaml"  # windIO v2 (tests/IEA-15-240-RWT.yaml is v1: no windio2 topology, the gmsh mesher cannot read it)
CASES = {
    "uniform": dict(element_size=0.25),
    "n10": dict(element_size=0.25, chord_divisions=10, span_aspect=1.0, min_band_elements=2),
    "n20": dict(element_size=0.25, chord_divisions=20, span_aspect=1.0, min_band_elements=2),
}
TIP = (109.0, 117.0)
STATIONS_FIXED = {"root": 0.5, "mid_span": 58.0, "z100": 100.0, "z115": 115.0}


def pct(a, qs=(50, 95, 100)):
    a = np.asarray(a, float)
    return {f"p{q}": float(np.percentile(a, q)) for q in qs} if len(a) else {}


def pct_low(a):  # lower-is-worse quantities: p50 / p5 / min
    a = np.asarray(a, float)
    return {"p50": float(np.percentile(a, 50)), "p5": float(np.percentile(a, 5)), "min": float(a.min())}


def physical_rings(raw, mesh_mid):
    """Rings as the BEM projector defines them (physical_ring_labels, RING_GAP_FRACTION * span, in-ring edges of slope < 0.3),
    on the mid-surface mesh; returns (list of RAW node index arrays, raw->mesh position map).  Node sets include web-only nodes."""
    from aeroelast.solvers.bem.force_projection import RING_GAP_FRACTION, physical_ring_labels
    from scipy.spatial import cKDTree

    mid = np.asarray(__import__("v2_blade")._inset(raw))
    coords = mesh_mid.coords_array
    _, pos = cKDTree(coords).query(mid)
    inv = np.empty(len(pos), int)
    inv[pos] = np.arange(len(pos))
    span = coords[:, 2]
    lab = physical_ring_labels(mesh_mid, span, RING_GAP_FRACTION * float(span.max() - span.min()))
    lab_raw = lab[pos]
    order = np.argsort(lab_raw, kind="stable")
    cuts = np.nonzero(np.diff(lab_raw[order]))[0] + 1
    return [g for g in np.split(order, cuts)], pos


def ring_stats(x, rings, skin_nodes):
    out = []
    for r in rings:
        rs = r[np.isin(r, skin_nodes)]
        if len(rs) < 3:
            continue
        z = x[rs, 2]
        p = x[rs]
        # best-fit plane z = a + b x + c y: residual separates a tilted-but-flat ring from a zig-zag one
        A = np.c_[np.ones(len(rs)), p[:, 0], p[:, 1]]
        coef, *_ = np.linalg.lstsq(A, z, rcond=None)
        res = z - A @ coef
        d = np.linalg.norm(p[:, None, :2] - p[None, :, :2], axis=-1)
        out.append({"n": len(rs), "zmean": float(z.mean()), "spread": float(np.abs(z - z.mean()).max()),
                    "plane_resid": float(np.abs(res).max()), "chord": float(d.max())})
    return out


def summarise_rings(st, sel=None):
    s = [r for r in st if sel is None or sel(r)]
    if not s:
        return {"n_rings": 0}
    sp = np.array([r["spread"] for r in s])
    worst = sorted(s, key=lambda r: -r["spread"])[:5]
    pr = np.array([r["plane_resid"] for r in s])
    wp = sorted(s, key=lambda r: -r["plane_resid"])[:5]
    return {"n_rings": len(s), "spread_m": pct(sp), "plane_resid_m": pct(pr), "n_plane_resid_over_1mm": int((pr > 1e-3).sum()),
            "worst5_plane": [{"z": round(r["zmean"], 3), "resid_m": r["plane_resid"], "nodes": r["n"]} for r in wp], "n_over_1mm": int((sp > 1e-3).sum()), "n_over_10mm": int((sp > 1e-2).sum()),
            "worst5": [{"z": round(r["zmean"], 3), "spread_m": r["spread"], "nodes": r["n"]} for r in worst]}


def run(label):
    from aeroelast.models.blade.numad.mesh_gen.gmsh_skin import DEFAULT_LIMITS, element_metrics, quality_gate
    import v2_blade

    t0 = time.time()
    raw = v2_blade.load_raw(mesher="gmsh", yaml_file=YAML, **CASES[label])
    tgen = time.time() - t0
    x = np.asarray(raw["nodes"], float)
    el = np.asarray(raw["elements"])
    surf = np.asarray(raw["surfaces"])
    res = {"label": label, "args": CASES[label], "cache_key": str(raw["key"]), "wall_s_incl_cache": tgen}

    skin = surf == "skin"
    res["counts"] = {"nodes": len(x), "elements": len(el), "skin_elements": int(skin.sum()),
                     "web_elements": int((~skin).sum()), "triangles": int((el[:, 3] < 0).sum())}

    # (b) element quality on the OML mesh the mesher gated (skin + webs)
    m = element_metrics(x, el)
    rep = quality_gate(m, raise_on_fail=False)
    res["quality"] = {"aspect": pct(m["aspect"]), "min_angle_deg": pct_low(m["min_angle"]),
                      "max_angle_deg": pct(m["max_angle"]), "jacobian_min_sine": pct_low(m["jacobian"]),
                      "gate_limits": rep["limits"], "gate_violations": rep["violations"], "gate_passed": rep["passed"]}
    zc = x[el.clip(0)][:, :, 2].mean(1)
    tipel = (zc >= TIP[0]) & (zc <= TIP[1])
    res["quality_tip"] = {"elements": int(tipel.sum()), "aspect": pct(m["aspect"][tipel]),
                          "min_angle_deg": pct_low(m["min_angle"][tipel]), "jacobian_min_sine": pct_low(m["jacobian"][tipel])}

    # (c)(d) ring planarity of the skin rows (OML nodes)
    mesh_mid, props = v2_blade.mesh_from_raw(raw, "mid")
    rings, pos = physical_rings(raw, mesh_mid)
    skin_nodes = np.unique(el[skin][el[skin] >= 0])
    st = ring_stats(x, rings, skin_nodes)
    res["rings_all"] = summarise_rings(st)
    res["rings_tip"] = summarise_rings(st, lambda r: TIP[0] <= r["zmean"] <= TIP[1])
    full = [r for r in st if r["n"] >= 8]  # stub groups (<8 skin nodes) carry no section
    ring_z = np.array([r["zmean"] for r in full])
    ring_ch = np.array([r["chord"] for r in full])
    res["rings_all"]["n_stub_groups_lt8"] = len(st) - len(full)

    # (e) chordwise skin elements at stations
    skin_idx = np.nonzero(skin)[0]
    zmin = np.array([x[el[k][el[k] >= 0], 2].min() for k in skin_idx])
    zmax = np.array([x[el[k][el[k] >= 0], 2].max() for k in skin_idx])
    stations = dict(STATIONS_FIXED)
    stations["max_chord"] = float(ring_z[np.argmax(ring_ch)])
    sec = {}
    for name, z in sorted(stations.items(), key=lambda t: t[1]):
        i = int(np.argmin(np.abs(ring_z - z)))
        zz = float(ring_z[i])  # use a ring-centred cut slightly above the ring so one element row is cut
        zq = zz + 1e-6 if name != "root" else zz + 1e-6
        cut = (zmin < zq) & (zmax > zq)
        n = int(cut.sum())
        c = float(ring_ch[i])
        sec[name] = {"z_target": z, "z_ring": zz, "chord_m": c, "skin_elements_around": n, "elements_per_chord": n / 2,
                     "mean_chordwise_size_m": float(2 * c / max(n, 1))}
    res["sections"] = sec

    # (f) unrealisable rings through the BEM projector helpers (cheap: ring_section per ring on the mid mesh)
    res["realisable"] = realisable(raw, rings, pos, mesh_mid, props)

    # VTU
    import meshio

    cells = {"quad": [], "triangle": []}
    cid = {"quad": [], "triangle": []}
    for k, e in enumerate(el):
        t = "quad" if e[3] >= 0 else "triangle"
        cells[t].append(e[e >= 0])
        cid[t].append(k)
    cd, blocks = {}, []
    for t in ("quad", "triangle"):
        if cells[t]:
            blocks.append((t, np.array(cells[t])))
            ks = np.array(cid[t])
            for name, arr in (("aspect", m["aspect"]), ("min_angle", m["min_angle"]), ("jacobian", m["jacobian"]),
                              ("is_web", (~skin).astype(float))):
                cd.setdefault(name, []).append(arr[ks])
    # node data: ring z deviation from ring mean (skin nodes)
    dz = np.zeros(len(x))
    for r in rings:
        dz[r] = x[r, 2] - x[r, 2].mean()  # ring_dz: z deviation from the physical ring mean (OML)
    meshio.Mesh(x, blocks, point_data={"ring_dz": dz}, cell_data=cd).write(HERE / f"{label}.vtu")
    res["wall_s_total"] = time.time() - t0
    (HERE / f"{label}.json").write_text(json.dumps(res, indent=1))
    print(label, "done", res["counts"], f"{res['wall_s_total']:.0f}s")


def realisable(raw, rings, pos, mesh, props):
    """Rings of the projector (physical_ring_labels groups, >= 3 nodes) whose ring_section is not realisable
    (no closed cell -> ValueError/None, or a non-positive wall stiffness) on the mid mesh with the real laminates."""
    try:
        from aeroelast.solvers.bem.force_projection import _RingSection, _ring_section_is_realisable, ring_section

        coords = mesh.coords_array
        span = np.array([0.0, 0.0, 1.0])
        bad, kinds, tot = [], {"no_cell": 0, "zero_stiffness": 0}, 0
        for r in rings:
            if len(r) < 3:
                continue
            tot += 1
            zm = float(coords[pos[r], 2].mean())
            try:
                cells, adj, el_len, stiff = ring_section(mesh, pos[r], span, props, coords=coords)
            except ValueError:
                bad.append(zm); kinds["no_cell"] += 1
                continue
            sec = _RingSection(cells, adj, el_len, stiff, np.arange(len(r)), from_element_properties=True)
            if not _ring_section_is_realisable(sec):
                bad.append(zm); kinds["zero_stiffness"] += 1
        bz = np.array(bad)
        return {"rings_checked": tot, "unrealisable": len(bad), "kinds": kinds,
                "z_min_max": [float(bz.min()), float(bz.max())] if len(bad) else None,
                "in_tip_109_117": int(((bz >= TIP[0]) & (bz <= TIP[1])).sum()) if len(bad) else 0,
                "z_sorted_first40": [round(b, 2) for b in sorted(bad)[:40]],
                "note": "rings = projector's physical_ring_labels groups on the whole mid mesh (skin+web), span +z; ring_section with real laminates; no BEM strips involved"}
    except Exception as exc:
        return {"error": f"{type(exc).__name__}: {exc}"}


def merge():
    out = {l: json.loads((HERE / f"{l}.json").read_text()) for l in CASES if (HERE / f"{l}.json").exists()}
    (HERE / "metrics.json").write_text(json.dumps(out, indent=1))
    print("wrote metrics.json:", list(out))


if __name__ == "__main__":
    if sys.argv[1] == "run":
        run(sys.argv[2])
    else:
        merge()
