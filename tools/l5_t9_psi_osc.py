"""L5 T9d: cause of the +-1.5 deg oscillation of the shell principal angle psi along z (wavelength ~6 m).

One factor at a time on the same shell segment solve: fit window half width, mesh N (chord_divisions), station / segment length.
For each factor level: sliding fine profile (step 0.25 m) of the bias-corrected psi (shell - (deck same fit - deck exact)), the
deck same-fit psi (fit-only content), the detrended shell oscillation, its dominant wavelength (zero-padded FFT of a cubic-detrended
series) and its correlation with the mean layup thickness / ply count per z bin of the segment.
  run:      python tools/l5_t9_psi_osc.py run --z 80 --nchord 8 --cd 20 --halves 1 2 3 4 --out DIR
  analyze:  python tools/l5_t9_psi_osc.py analyze DIR [DIR ...]
"""
import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import l5_t9_arbiter as A  # noqa: E402

STEP = 0.25


def run(z, nchord, cd, halves, out):
    from tests.support.v2_blade import build_v2_blade

    t0 = time.time()
    mesh, props = build_v2_blade(0.25, skin_reference="mid", mesher="gmsh", es_range=(0.05, 0.25), chord_divisions=cd)
    deck_f = A.deck_C_of_z(A.load_deck()[0])
    sub, p, chord, L = A.shell_segment(mesh, props, z, nchord)
    X, Us, area, info = A.solve_shell_couples(sub, p)
    zlo, zhi = info["z_range"]
    ref_xy = X[:, :2].mean(axis=0)
    Ub = [A.beam_displacements(X, lambda zz: deck_f(zz)[0], zlo, zhi, A.M0 * np.eye(3)[j], ref_xy) for j in range(3)]
    series = {}
    for half in halves:
        rows = []
        for zc in np.arange(zlo + half, zhi - half + 1e-9, STEP):
            try:
                fs = [A.fit_section_kinematics(X, U, zc, half, weights=area) for U in Us]
                fb = [A.fit_section_kinematics(X, U, zc, half, weights=area) for U in Ub]
            except ValueError:
                continue
            ms = A.section_metrics(A.compliance_from_responses([f["kappa"] for f in fs], [A.M0] * 3))
            mb = A.section_metrics(A.compliance_from_responses([f["kappa"] for f in fb], [A.M0] * 3))
            me = A.section_metrics(deck_f(zc)[0])
            rows.append([zc, ms["psi_deg"], mb["psi_deg"], me["psi_deg"], ms["GJ"], mb["GJ"], me["GJ"]])
        series[f"{half:g}"] = rows
    ze, he, ne, ex, ey, ea, ek = [], [], [], [], [], [], []
    name_of = {}
    for s, es in sub.element_sets.items():
        if s in p:
            for e in es.elements:
                name_of[e.id] = s
    edge = []
    for e in sub.elements:
        q = np.array([n.coords for n in e.nodes], float)
        ze.append(q[:, 2].mean())
        ex.append(q[:, 0].mean())
        ey.append(q[:, 1].mean())
        ea.append(0.5 * float(np.linalg.norm(np.cross(q[2] - q[0], q[3] - q[1]))) if len(q) == 4 else
                  0.5 * float(np.linalg.norm(np.cross(q[1] - q[0], q[2] - q[0]))))
        lam = p.get(name_of.get(e.id))
        he.append(float(lam.total_thickness) if lam is not None else np.nan)
        ne.append(int(lam.n_plies) if lam is not None else -1)
        a = np.asarray(lam.a_matrix(), float) if lam is not None else np.zeros((3, 3))
        ek.append(float(0.5 * (a[0, 0] + a[1, 1])))  # membrane stiffness proxy (N/m), orientation ignored
        edge.append(float(np.linalg.norm(q[1] - q[0])))
    node_z = np.sort(np.unique(np.round(X[:, 2], 4)))
    doc = {"z": z, "nchord": nchord, "cd": cd, "chord": chord, "L": L, "z_range": [zlo, zhi], "halves": list(halves), "series": series,
           "elem_z": ze, "elem_h": he, "elem_nply": ne, "elem_x": ex, "elem_y": ey, "elem_area": ea, "elem_kmem": ek, "elem_edge_median": float(np.median(edge)), "n_nodes": len(X),
           "n_elements": len(ze), "elapsed_s": time.time() - t0}
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    json.dump(doc, open(out / "psi_osc.json", "w"))
    print(f"done z={z:g} n={nchord:g} cd={cd} nodes={len(X)} L={zhi - zlo:.1f} median edge={doc['elem_edge_median']:.3f} ({doc['elapsed_s']:.0f} s)")


def dominant(zs, y, deg=3, lo=2.0, hi=20.0):
    """Dominant wavelength (m) of y(z) after a polynomial detrend: peak of the zero-padded FFT amplitude within wavelengths [lo, hi]."""
    zs, y = np.asarray(zs), np.asarray(y)
    r = y - np.polyval(np.polyfit(zs - zs.mean(), y, deg), zs - zs.mean())
    r = r * np.hanning(len(r))
    n = 1 << 14
    f = np.fft.rfftfreq(n, STEP)[1:]
    a = np.abs(np.fft.rfft(r, n))[1:]
    m = (1 / f >= lo) & (1 / f <= hi)
    return float(1 / f[m][np.argmax(a[m])]), float(np.std(r) * np.sqrt(2)), r


def analyze(dirs):
    print("| z | nchord | cd | half | n pts | corrected-psi osc wavelength m | amp deg (sqrt2 std) | deck-same-fit osc wavelength m | amp deg | shell GJ osc wavelength m |")
    print("|---|---|---|---|---|---|---|---|---|---|")
    for d in dirs:
        o = json.load(open(Path(d) / "psi_osc.json"))
        for h, rows in o["series"].items():
            a = np.array(rows)
            if len(a) < 12:
                continue
            corr = a[:, 1] - (a[:, 2] - a[:, 3])
            wl, amp, _ = dominant(a[:, 0], corr)
            wb, ab, _ = dominant(a[:, 0], a[:, 2] - a[:, 3])
            wg, _, _ = dominant(a[:, 0], np.log(a[:, 4] / a[:, 5]))
            print(f"| {o['z']:g} | {o['nchord']:g} | {o['cd']} | {h} | {len(a)} | {wl:.1f} | {amp:.2f} | {wb:.1f} | {ab:.2f} | {wg:.1f} |")
        ze, he, ne = (np.array(o[k], float) for k in ("elem_z", "elem_h", "elem_nply"))
        bins = np.arange(o["z_range"][0], o["z_range"][1] + 1e-9, STEP)
        idx = np.digitize(ze, bins)
        hm = np.array([np.nanmean(he[idx == i]) if (idx == i).any() else np.nan for i in range(1, len(bins))])
        zb = 0.5 * (bins[1:] + bins[:-1])
        ok = np.isfinite(hm)
        w, a_, _ = dominant(zb[ok], hm[ok] / np.nanmean(hm))
        print(f"  [{d}] mean layup thickness per 0.25 m bin: dominant wavelength {w:.1f} m, rel amp {a_:.3f}; median element edge {o['elem_edge_median']:.3f} m")
        # correlation of the corrected-psi oscillation (half 2) with mean thickness oscillation
        if "2" in o["series"]:
            a = np.array(o["series"]["2"])
            corr = a[:, 1] - (a[:, 2] - a[:, 3])
            _, _, r = dominant(a[:, 0], corr)
            hz = np.interp(a[:, 0], zb[ok], hm[ok])
            hr = hz - np.polyval(np.polyfit(a[:, 0] - a[:, 0].mean(), hz, 3), a[:, 0] - a[:, 0].mean())
            print(f"  [{d}] corr(psi osc, thickness osc) at half 2 = {np.corrcoef(r, hr)[0, 1]:+.2f}")


def yaml_events(zlo=64.0, zhi=94.0, L=117.0):
    """z of every grid point of the official windIO layers / webs at which a field (thickness, arc start / end, orientation, offset)
    changes slope, plus start / end of its support, inside [zlo, zhi] (L = blade length; eta = z / L as in ``A.eta_for``)."""
    import yaml

    d = yaml.safe_load(open(A.REPO / "tests" / "IEA15MW" / "reference" / "upstream" / "IEA-15-240-RWT.yaml"))
    s = d["components"]["blade"]["structure"]
    out = []
    for item in list(s["webs"]) + list(s["layers"]):
        for k, f in item.items():
            if not (isinstance(f, dict) and "grid" in f and "values" in f):
                continue
            g, v = np.array(f["grid"], float), np.array(f["values"], float)
            for i in range(1, len(g) - 1):
                s1, s2 = (v[i] - v[i - 1]) / (g[i] - g[i - 1] + 1e-12), (v[i + 1] - v[i]) / (g[i + 1] - g[i] + 1e-12)
                if zlo <= g[i] * L <= zhi and abs(s1 - s2) > 1e-6 * max(abs(s1), abs(s2), 1e-9):
                    out.append((round(float(g[i] * L), 1), f"{item['name']}:{k}"))
    return sorted(out)


def geom_psi(o, half, zs):
    """Principal angle (deg, from the mesh x axis) of the second moments of area about the weighted centroid of the shell elements in
    |z - zc| <= half, for two weights: thickness*area ('h') and membrane-stiffness proxy (a11+a22)/2 * area ('k').  Pure layup/geometry
    read: no FE solve, no fit."""
    ez, ex, ey, ea = (np.array(o[k], float) for k in ("elem_z", "elem_x", "elem_y", "elem_area"))
    eh, ek = np.array(o["elem_h"], float), np.array(o["elem_kmem"], float)
    res = {"h": [], "k": []}
    for zc in zs:
        m = np.abs(ez - zc) <= half
        for nm, w in (("h", eh * ea), ("k", ek * ea)):
            ww = w[m]
            cx, cy = np.sum(ww * ex[m]) / ww.sum(), np.sum(ww * ey[m]) / ww.sum()
            ixx, iyy, ixy = np.sum(ww * (ey[m] - cy) ** 2), np.sum(ww * (ex[m] - cx) ** 2), -np.sum(ww * (ex[m] - cx) * (ey[m] - cy))
            res[nm].append(np.arctan2(-2 * ixy, ixx - iyy))
    return {k: np.degrees(0.5 * np.unwrap(np.array(v))) for k, v in res.items()}


def geom(d, half=2.0):
    """Does the same oscillation appear in the layup/geometry second moments (no FE)?  Prints wavelengths and correlation of
    detrended series with the FE corrected psi at the same window."""
    o = json.load(open(Path(d) / "psi_osc.json"))
    a = np.array(o["series"][f"{half:g}"])
    zs = a[:, 0]
    corr = a[:, 1] - (a[:, 2] - a[:, 3])
    g = geom_psi(o, half, zs)
    wl, amp, rfe = dominant(zs, corr)
    print(f"[{d}] FE corrected psi, half {half:g}: wavelength {wl:.1f} m amp {amp:.2f} deg")
    for k, v in g.items():
        w, am, r = dominant(zs, v)
        print(f"   geometry psi ({k} weight): wavelength {w:.1f} m amp {am:.2f} deg; corr with FE detrended = {np.corrcoef(rfe, r)[0, 1]:+.2f}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["run", "analyze", "geom"])
    ap.add_argument("dirs", nargs="*")
    ap.add_argument("--z", type=float, default=80.0)
    ap.add_argument("--nchord", type=float, default=8.0)
    ap.add_argument("--cd", type=int, default=20)
    ap.add_argument("--halves", type=float, nargs="+", default=[1.0, 2.0, 3.0, 4.0])
    ap.add_argument("--out", type=Path)
    a = ap.parse_args()
    if a.cmd == "run":
        run(a.z, a.nchord, a.cd, a.halves, a.out)
    elif a.cmd == "geom":
        for d in a.dirs:
            geom(d)
    else:
        analyze(a.dirs)
