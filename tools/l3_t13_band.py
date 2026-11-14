"""T13 (L3 D3): does the z 90-92 in-plane distortion move local GJ / twist, and which walls carry it.

Band G12 variant (same device as tools/l3_w7_tipstiff.py, but a finite band): G12 of orthotropic plies x F for element centroids with
Z0 <= z <= Z1 (0.5 m bins, one-bin ramps on both sides).  Spec ``band:<F>:<Z0>:<Z1>``.

    python tools/l3_t13_band.py N OUT.json [--variant band:F:Z0:Z1]       (S-7 case: 1000 N m wall-flow torque on TipNodes, root clamped)
Output: slab rotation/residual (0.5 m slabs), mesh-cut Bredt on a 0.5 m grid (as tools/l3_w7_tip.py), and per 0.5 m z bin from 86 to
97 the share of the squared residual carried by node classes (node class = median thickness of the elements touching the node:
skin < 3 mm, panel 3-15 mm, cap > 15 mm), against the node share of the class.
"""

import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import l3_w4_tip_compliance as w4  # noqa: E402
import l3_w6_scale as w6  # noqa: E402
from l3_w4_bredt import make  # noqa: E402
from tests.support.paths import REPO_ROOT  # noqa: E402
from tests.support.v2_blade import build_v2_blade, load_raw  # noqa: E402


def make_transform(factor, z0, z1):
    def fz(z):
        up = np.clip((z - (z0 - 0.5)) / 0.5, 0, 1)
        dn = np.clip(((z1 + 0.5) - z) / 0.5, 0, 1)
        return 1.0 + (factor - 1.0) * float(up * dn)

    def transform(mesh, props):
        from aeroelast.core.mesh.entities import ElementSet
        from aeroelast.models.blade.model import build_rust_properties

        numad = mesh.numad_mesh_data
        mats = {m["name"]: m for m in numad["materials"]}
        layups = {s["elementSet"]: s["layup"] for s in numad["sections"]}
        pos = {n.id: n.coordinates if hasattr(n, "coordinates") else n.coords for n in mesh.nodes}
        new_props, new_mats, secs = dict(props), {}, []
        for name in [n for n in mesh.element_sets if n in layups]:
            groups: dict[float, list] = {}
            for e in mesh.element_sets[name].elements:
                zc = float(np.mean([pos[i][2] for i in e.node_ids]))
                zb = (np.floor(zc / 0.5) + 0.5) * 0.5
                groups.setdefault(round(fz(zb), 4), []).append(e)
            for f, g in groups.items():
                if f == 1.0:
                    continue
                sname = f"{name}__B{f:.4f}"
                mesh.element_sets[name].elements.difference_update(g)
                mesh.add_element_set(ElementSet(name=sname, elements=set(g)))
                lay = []
                for n, th, a in layups[name]:
                    nn = w6.key(n, f, 1.0)
                    new_mats.setdefault(
                        nn,
                        {
                            "name": nn,
                            "density": mats[n]["density"],
                            "elastic": w6.scale_material(mats[n]["elastic"], f, 1.0),
                        },
                    )
                    lay.append([nn, th, a])
                secs.append({"elementSet": sname, "layup": lay})
        new_props.update(
            build_rust_properties({"materials": list(new_mats.values()), "sections": secs})
        )
        for name in [
            n for n in layups if n in mesh.element_sets and not mesh.element_sets[n].elements
        ]:
            new_props.pop(name, None)
        return mesh, new_props

    return transform


def from_spec(spec):
    _, f, z0, z1 = spec.split(":")
    return make_transform(float(f), float(z0), float(z1))


def node_classes(mesh, props):
    idx = mesh.node_id_to_index
    acc = [[] for _ in range(mesh.node_count)]
    for name, es in mesh.element_sets.items():
        if name not in props:
            continue
        t = float(props[name].total_thickness)
        for e in es.elements:
            for i in e.node_ids:
                acc[idx[i]].append(t)
    tm = np.array([np.median(a) if a else np.nan for a in acc])
    return np.where(tm < 3e-3, 0, np.where(tm < 15e-3, 1, 2))  # 0 skin, 1 panel, 2 cap


def class_shares(coords, u, cls, zc, half=0.25):
    idx = np.nonzero(abs(coords[:, 2] - zc) <= half)[0]
    xy = coords[idx, :2]
    r = xy - xy.mean(axis=0)
    dz = coords[idx, 2] - zc
    n = len(idx)
    A = np.zeros((2 * n, 4))
    A[:n, 0] = 1
    A[n:, 1] = 1
    A[:n, 2] = -r[:, 1]
    A[n:, 2] = r[:, 0]
    A[:n, 3] = -r[:, 1] * dz
    A[n:, 3] = r[:, 0] * dz
    b = np.concatenate([u[idx * 6], u[idx * 6 + 1]])
    c = np.linalg.lstsq(A, b, rcond=None)[0]
    res = b - A @ c
    e = res[:n] ** 2 + res[n:] ** 2
    out = {}
    for k, nm in enumerate(("skin", "panel", "cap")):
        m = cls[idx] == k
        out[nm] = {"node_share": float(m.mean()), "res2_share": float(e[m].sum() / e.sum())}
    return out


def main():
    a = sys.argv[1:]
    n, out = int(a[0]), a[1]
    variant = a[a.index("--variant") + 1] if "--variant" in a else None
    t0 = time.time()
    mesh, props = build_v2_blade(0.25, skin_reference="mid", mesher="gmsh", chord_divisions=n)
    if variant:
        mesh, props = from_spec(variant)(mesh, props)
    cls = node_classes(mesh, props)
    u, nring = w4.solve(mesh, props, None)
    coords = mesh.coords_array
    zmax = float(coords[:, 2].max())
    sl = [
        s
        for s in (w4.slab_stats(coords, u, 6, zc, 0.25) for zc in np.arange(58.0, zmax, 0.5))
        if s is not None
    ]
    raw = load_raw(0.25, mesher="gmsh", generate=False, chord_divisions=n)
    deck = str(REPO_ROOT / "tests/IEA15MW/reference/IEA-15-240-RWT_BeamDyn_blade.dat")
    gj_at, _, _ = make(raw, deck)
    zb = np.arange(58.0, 116.01, 0.5)
    shares = {f"{z:.2f}": class_shares(coords, u, cls, z) for z in np.arange(86.0, 97.01, 0.5)}
    res = {
        "N": n,
        "variant": variant,
        "n_nodes": int(len(coords)),
        "zmax": zmax,
        "ring_nodes": int(nring),
        "torque": w4.TORQUE,
        "slab_z": [s["z"] for s in sl],
        "slab_th": [s["th"] for s in sl],
        "slab_th_aff": [s["th_aff"] for s in sl],
        "slab_res_rigid": [s["res_rigid"] for s in sl],
        "slab_rrms": [s["rrms"] for s in sl],
        "bredt_z": zb.tolist(),
        "bredt_off0137": [gj_at(z + 0.0137)[0] for z in zb],
        "class_shares": shares,
        "total_s": round(time.time() - t0, 1),
    }
    json.dump(res, open(out, "w"))
    print("done", res["total_s"], flush=True)


if __name__ == "__main__":
    main()
