"""L3 T7 (W6): outboard shell-stiffness scaling for the torque/power sensitivity (odd/tasks/l3-outboard-section-stiffness.md).

Scalings (eta = span fraction; the factor is 1 up to 0.60 and ramps linearly to its value at 0.65, then piecewise linear on ``ETA_NODES``; targets capped at +25 % over the baseline):
  V1 (GJ):   in-plane shear modulus G12 of every orthotropic ply (skin and webs; isotropic gelcoat/foam untouched) times fG(eta);
  V2 (edge): fibre modulus E1 of every orthotropic ply of the skin elements in the edge-lever zone (chordwise offset from the
             section mean > ZONE_FRAC of the chord, i.e. the LE and TE panels) times fE(eta);
  V3: both.
The factors come from a closed-form thin-wall calibration on the same gmsh chord-scaled N20 AR 1 mesh cut (mid-surface):
Bredt multi-cell GJ (``tools/l3_w4_bredt.make``) against the deck K66 and the EA-weighted section EI_edge (pattern of
``tools/l3_w3_ei_edge.py``) against max(K44, K55); a station already at or above the deck keeps factor 1.

    scripts/aeroenv.sh env PYTHONPATH=$PWD python tools/l3_w6_scale.py calibrate logs/chord_scaled/w6/factors.json
"""

from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tools"))

DECK = ROOT / "tests/IEA15MW/reference/IEA-15-240-RWT_BeamDyn_blade.dat"
ETA_NODES = [0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00]
ZONE_FRAC = 0.30
N_CHORD = 20
CAP = 1.25


def raw_case():
    os.environ.setdefault("V2_BLADE_MESHER", "gmsh")
    from tests.support.v2_blade import load_raw

    return load_raw(0.25, mesher="gmsh", chord_divisions=N_CHORD, es_range=(0.05, 0.25))


def factor_at(eta, nodes, values):
    """Piecewise-linear factor; 1 for eta < 0.6 (node 0.55 carries 1), clamp beyond the last node."""
    return np.interp(eta, nodes, values)


def element_geometry(centroids_xyz, zmin, zmax, is_skin):
    """eta and the chordwise zone flag per element (same formula for raw arrays and for the MeshModel)."""
    z = centroids_xyz[:, 2]
    eta = (z - zmin) / (zmax - zmin)
    zone = np.zeros(len(z), bool)
    bins = np.floor((z - zmin) / 1.0).astype(int)  # 1 m span bins
    for b in np.unique(bins):
        sel = bins == b
        sk = sel & is_skin
        if sk.sum() < 4:
            continue
        x = centroids_xyz[sk, 0]
        ref, chord = x.mean(), np.ptp(x)
        zone[sk] = np.abs(centroids_xyz[sk, 0] - ref) > ZONE_FRAC * chord
    return eta, zone


def scale_material(m, fg, fe):
    """Copy of a numad/tables material with G12 x fg and E1 x fe (orthotropic only)."""
    out = dict(m)
    if np.ndim(m["E"]) == 0:
        return out
    E, G = list(map(float, m["E"])), list(map(float, m["G"]))
    E[0] *= fe
    G[0] *= fg
    out["E"], out["G"] = E, G
    return out


def key(name, fg, fe):
    return name if (fg == 1.0 and fe == 1.0) else f"{name}|{fg:.5f}|{fe:.5f}"


def expand_tables(tables, layup_id, fg_el, fe_el):
    """Per-element-group tables with scaled materials: returns (tables2, layup_id2) for the thin-wall tools."""
    combos, new_lay, new_id = {}, [], np.empty(len(layup_id), int)
    mats = dict(tables["materials"])
    for k, (l, fg, fe) in enumerate(zip(layup_id, fg_el, fe_el, strict=True)):
        c = (int(l), round(float(fg), 5), round(float(fe), 5))
        if c not in combos:
            lay = []
            for n, th, ang in tables["layups"][c[0]]:
                nn = key(n, c[1], c[2])
                if nn not in mats:
                    mats[nn] = scale_material(tables["materials"][n], c[1], c[2])
                lay.append([nn, th, ang])
            combos[c] = len(new_lay)
            new_lay.append(lay)
        new_id[k] = combos[c]
    return {"materials": mats, "layups": new_lay}, new_id


def read_deck(path=DECK):
    d = open(path).read().split("\n")
    i = [k for k, ln in enumerate(d) if "DISTRIBUTED PROPERTIES" in ln][0] + 1
    n = int([ln for ln in d if "station_total" in ln][0].split()[0])
    eta, K = [], []
    while i < len(d) and len(eta) < n:
        t = d[i].split()
        if len(t) == 1:
            eta.append(float(t[0]))
            K.append(np.array([[float(x) for x in d[i + 1 + r].split()] for r in range(6)]))
            i += 14
        else:
            i += 1
    return np.array(eta), np.array(K)


def _g(v, i):
    return float(np.atleast_1d(v)[min(i, np.size(v) - 1)])


def q11bar0(m):  # 0 deg ply, all layups here are at 0 deg
    E1, E2, nu12 = _g(m["E"], 0), _g(m["E"], 1), _g(m["nu"], 0)
    nu21 = nu12 * E2 / E1
    return E1 / (1 - nu12 * nu21)


def thin_wall(raw, tables2, lid2, zmin, zmax, etas, X):
    """EA, EI_edge (max principal), EI_flap (min), Bredt GJ per eta on the mid-surface cut (copy of l3_w3_ei_edge.section)."""
    import l3_w4_bredt

    E = np.asarray(raw["elements"])
    a11 = np.array(
        [
            sum(q11bar0(tables2["materials"][n]) * th for n, th, a in lay)
            for lay in tables2["layups"]
        ]
    )[lid2]
    zz = X[:, 2]
    emin = np.array([zz[e[e >= 0]].min() for e in E])
    emax = np.array([zz[e[e >= 0]].max() for e in E])
    raw2 = dict(raw)
    raw2["tables_json"] = json.dumps(tables2).encode()
    raw2["layup_id"] = lid2
    gj_at, _, _ = l3_w4_bredt.make(raw2, str(DECK))
    out = {}
    L = zmax - zmin
    for eta in etas:
        zc = zmin + eta * L
        mids, ln, tv, w = [], [], [], []
        for k in np.nonzero((emin < zc) & (emax > zc))[0]:
            ids = E[k][E[k] >= 0]
            pts = []
            for a, b in zip(ids, np.roll(ids, -1), strict=True):
                za, zb = zz[a] - zc, zz[b] - zc
                if za * zb < 0:
                    pts.append(X[a, :2] + za / (za - zb) * (X[b, :2] - X[a, :2]))
            if len(pts) == 2:
                d = pts[1] - pts[0]
                length = float(np.linalg.norm(d))
                if length > 0:
                    mids.append(0.5 * (pts[0] + pts[1]))
                    ln.append(length)
                    tv.append(d / length)
                    w.append(a11[k] * length)
        mids, ln, tv, w = map(np.array, (mids, ln, tv, w))
        ea = w.sum()
        c = (w[:, None] * mids).sum(0) / ea
        d = mids - c
        own = w * ln**2 / 12.0
        ixx = (w * d[:, 1] ** 2).sum() + (own * tv[:, 1] ** 2).sum()
        iyy = (w * d[:, 0] ** 2).sum() + (own * tv[:, 0] ** 2).sum()
        ixy = (w * d[:, 0] * d[:, 1]).sum() + (own * tv[:, 0] * tv[:, 1]).sum()
        ev = np.linalg.eigvalsh(np.array([[ixx, -ixy], [-ixy, iyy]]))
        gj = gj_at(zc)[0]
        out[eta] = {
            "EA": float(ea),
            "EI_edge": float(ev[1]),
            "EI_flap": float(ev[0]),
            "GJ": float(gj),
        }
    return out


def calibrate(out_path):
    from tests.support.v2_blade import _inset

    raw = raw_case()
    tables = json.loads(bytes(raw["tables_json"]).decode())
    E = np.asarray(raw["elements"])
    Xo = np.asarray(raw["nodes"], float)
    X = _inset(raw)
    zmin, zmax = float(Xo[:, 2].min()), float(Xo[:, 2].max())
    cen = np.array([X[e[e >= 0]].mean(axis=0) for e in E])
    is_skin = np.array([str(s) == "skin" for s in raw["surfaces"]])
    eta_el, zone = element_geometry(cen, zmin, zmax, is_skin)
    lid = np.asarray(raw["layup_id"])
    deck_eta, K = read_deck()
    eta_cal = [0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95]
    deck_gj = {e: float(np.interp(e, deck_eta, K[:, 5, 5])) for e in eta_cal}
    deck_ed = {
        e: float(max(np.interp(e, deck_eta, K[:, 3, 3]), np.interp(e, deck_eta, K[:, 4, 4])))
        for e in eta_cal
    }

    def run(fG, fE):
        fg_el = np.where(eta_el >= 0.6, factor_at(eta_el, ETA_NODES, fG), 1.0)
        fe_el = np.where((eta_el >= 0.6) & zone, factor_at(eta_el, ETA_NODES, fE), 1.0)
        t2, l2 = expand_tables(tables, lid, fg_el, fe_el)
        return thin_wall(raw, t2, l2, zmin, zmax, eta_cal, X)

    one = [1.0] * len(ETA_NODES)
    base = run(one, one)
    # target = deck, but never more than +CAP over the baseline (the gap quoted by the spec is +15-25 %; beyond it the factor
    # would be a fit of the thin-wall model, not a stiffness the input can carry)
    tgt_gj = {e: min(deck_gj[e], CAP * base[e]["GJ"]) for e in eta_cal}
    tgt_ed = {e: min(deck_ed[e], CAP * base[e]["EI_edge"]) for e in eta_cal}
    fG, fE = list(one), list(one)
    idx = {e: ETA_NODES.index(e) for e in eta_cal}
    # first guess from the thin-wall ratios, then two fixed-point corrections
    for e in eta_cal:
        fG[idx[e]] = max(1.0, tgt_gj[e] / base[e]["GJ"])
        fE[idx[e]] = max(1.0, tgt_ed[e] / base[e]["EI_edge"])
    for it in range(6):
        cur = run(fG, fE)
        print(
            "iter",
            it,
            {
                e: (round(cur[e]["GJ"] / tgt_gj[e], 3), round(cur[e]["EI_edge"] / tgt_ed[e], 3))
                for e in eta_cal
            },
            flush=True,
        )
        for e in eta_cal:
            if base[e]["GJ"] < tgt_gj[e]:
                fG[idx[e]] *= tgt_gj[e] / cur[e]["GJ"]
            if base[e]["EI_edge"] < tgt_ed[e]:
                fE[idx[e]] *= tgt_ed[e] / cur[e]["EI_edge"]
    fG[0], fE[0], fG[1], fE[1] = 1.0, 1.0, 1.0, 1.0  # nodes 0.55 and 0.60: ramp 0.60 -> 0.65
    fG[-1], fE[-1] = fG[-2], fE[-2]  # node 1.00 (clamp the tip)
    v1 = run(fG, one)
    v2 = run(one, fE)
    v3 = run(fG, fE)
    res = {
        "eta_nodes": ETA_NODES,
        "fG": fG,
        "fE": fE,
        "zone_frac": ZONE_FRAC,
        "n_chord": N_CHORD,
        "deck": {"GJ": deck_gj, "edge": deck_ed},
        "target": {"GJ": tgt_gj, "edge": tgt_ed},
        "cap": CAP,
        "zone_elements_outboard": int((zone & (eta_el >= 0.6)).sum()),
        "thin_wall": dict(zip(("base", "V1", "V2", "V3"), (base, v1, v2, v3), strict=True)),
    }
    Path(out_path).write_text(json.dumps(res, indent=1, default=float))
    print("eta  | fG fE | GJ/deck b,V1,V3 | edge/deck b,V2,V3 | flap b,V2,V3 /base | EA V2/base")
    for e in eta_cal:
        b = base[e]
        print(
            f"{e:.2f} | {fG[idx[e]]:.3f} {fE[idx[e]]:.3f} | "
            + " ".join(f"{x[e]['GJ'] / deck_gj[e]:.3f}" for x in (base, v1, v3))
            + " | "
            + " ".join(f"{x[e]['EI_edge'] / deck_ed[e]:.3f}" for x in (base, v2, v3))
            + " | "
            + " ".join(f"{x[e]['EI_flap'] / b['EI_flap']:.3f}" for x in (base, v2, v3))
            + f" | {v2[e]['EA'] / b['EA']:.3f}"
            + f" | V1 edge/flap/EA vs base {v1[e]['EI_edge'] / b['EI_edge']:.3f} {v1[e]['EI_flap'] / b['EI_flap']:.3f} {v1[e]['EA'] / b['EA']:.3f}"
            + f" | V2 GJ vs base {v2[e]['GJ'] / b['GJ']:.3f}"
        )
    return 0


def make_transform(factors_path, do_g, do_e):
    """``mesh_transform`` for ``l3_torque_sensitivity.solve_rated``: rebuild the property of every layup set with scaled materials.

    Layup sets are station-specific (no set outboard of eta 0.6 is shared with an inboard one, checked), but the zone flag splits an
    element set in two (``<set>__E``); both get their own property."""
    f = json.loads(Path(factors_path).read_text())
    nodes, fG, fE = f["eta_nodes"], f["fG"], f["fE"]

    def transform(mesh, props):
        from aeroelast.core.mesh.entities import ElementSet
        from aeroelast.models.blade.model import build_rust_properties

        numad = mesh.numad_mesh_data
        mats = {m["name"]: m for m in numad["materials"]}
        layups = {s["elementSet"]: s["layup"] for s in numad["sections"]}
        pos = {n.id: n.coordinates if hasattr(n, "coordinates") else n.coords for n in mesh.nodes}
        elems = {}
        for name, es in mesh.element_sets.items():
            if name in layups:
                elems[name] = list(es.elements)
        allels = [e for v in elems.values() for e in v]
        cen = np.array([np.mean([pos[i] for i in e.node_ids], axis=0) for e in allels])
        zs = np.array([pos[i][2] for e in allels for i in e.node_ids])
        zmin, zmax = float(zs.min()), float(zs.max())
        skin = {e.id for e in mesh.get_element_set("allOuterShellEls").elements}
        is_skin = np.array([e.id in skin for e in allels])
        eta, zone = element_geometry(cen, zmin, zmax, is_skin)
        info = {e.id: (eta[i], bool(zone[i])) for i, e in enumerate(allels)}
        new_props, new_mats, secs, nmoved = dict(props), {}, [], 0
        for name, els in elems.items():
            groups = {False: [], True: []}
            for e in els:
                groups[info[e.id][1] and do_e].append(e)
            for flag, g in groups.items():
                if not g:
                    continue
                eta_g = float(np.mean([info[e.id][0] for e in g]))
                if eta_g < 0.6:
                    continue
                fg = float(factor_at(eta_g, nodes, fG)) if do_g else 1.0
                fe = float(factor_at(eta_g, nodes, fE)) if flag else 1.0
                fg, fe = round(fg, 5), round(fe, 5)
                sname = name
                if flag:
                    sname = name + "__E"
                    mesh.element_sets[name].elements.difference_update(g)
                    mesh.add_element_set(ElementSet(name=sname, elements=set(g)))
                lay = []
                for n, th, a in layups[name]:
                    nn = key(n, fg, fe)
                    new_mats.setdefault(
                        nn,
                        {
                            "name": nn,
                            "density": mats[n]["density"],
                            "elastic": scale_material(mats[n]["elastic"], fg, fe),
                        },
                    )
                    lay.append([nn, th, a])
                secs.append({"elementSet": sname, "layup": lay})
                nmoved += len(g)
        built = build_rust_properties({"materials": list(new_mats.values()), "sections": secs})
        new_props.update(built)
        for name in [n for n in elems if not mesh.element_sets[n].elements]:  # emptied sets
            new_props.pop(name, None)
        print(
            f"w6 transform: {len(secs)} property sets rebuilt for {nmoved} elements (do_g={do_g}, do_e={do_e})",
            flush=True,
        )
        return mesh, new_props

    return transform


if __name__ == "__main__":
    if len(sys.argv) >= 3 and sys.argv[1] == "calibrate":
        raise SystemExit(calibrate(sys.argv[2]))
    raise SystemExit(__doc__)
