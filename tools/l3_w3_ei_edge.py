"""L3 T6 (W3): thin-wall EI_edge / EI_flap / EA of the gmsh chord-scaled mesh cut (mid-surface and OML) vs the BeamDyn deck.

Pattern of tools/bredt_batho_bands.py: every shell element cut by the plane z = const gives one wall segment (length ds, in-plane
CLT A11 along the span, from the layup table with the ply angle measured from the span). Segment axial stiffness A11*ds sits at the
segment midpoint; EA = sum, centroid = EA-weighted mean, EI tensor about that centroid = sum A11 ds (d_x^2, d_y^2, d_x d_y) plus the
segment's own bending A11 ds^3/12 about its own axis. Principal values: max = edge (chordwise), min = flap. The deck K44/K55/K33 are
printed next to them. The deck K44/K55 are about the reference axis (K34 = -EA x_ec != 0): the ratios use the elastic-centre values
K44 - K34^2/K33, K55 - K35^2/K33 (the thin-wall tensor is about the centroid); the reference-axis ratio is stored too.
usage: python tools/l3_w3_ei_edge.py OUT.json [N=20]
"""

import json
import os
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
os.environ["V2_BLADE_MESHER"] = "gmsh"
from tests.support.v2_blade import _inset, load_raw  # noqa: E402

DECK = ROOT / "tests/IEA15MW/reference/IEA-15-240-RWT_BeamDyn_blade.dat"


def read_deck(path):
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


def q11bar(m, deg):
    g = lambda v, i: float(np.atleast_1d(v)[min(i, np.size(v) - 1)])  # noqa: E731
    E1, E2, G12, nu12 = g(m["E"], 0), g(m["E"], 1), g(m["G"], 0), g(m["nu"], 0)
    nu21 = nu12 * E2 / E1
    d = 1 - nu12 * nu21
    q11, q22, q12 = E1 / d, E2 / d, nu12 * E2 / d
    c, s = np.cos(np.radians(deg)), np.sin(np.radians(deg))
    return q11 * c**4 + q22 * s**4 + 2 * (q12 + 2 * G12) * s * s * c * c


def qbar(m, deg):
    """(Qbar11, Qbar22, Qbar12) of a ply, angle from the span."""
    g = lambda v, i: float(np.atleast_1d(v)[min(i, np.size(v) - 1)])  # noqa: E731
    E1, E2, G12, nu12 = g(m["E"], 0), g(m["E"], 1), g(m["G"], 0), g(m["nu"], 0)
    d = 1 - nu12 * nu12 * E2 / E1
    q11, q22, q12, q66 = E1 / d, E2 / d, nu12 * E2 / d, G12
    c, s = np.cos(np.radians(deg)), np.sin(np.radians(deg))
    return (
        q11 * c**4 + q22 * s**4 + 2 * (q12 + 2 * q66) * s * s * c * c,
        q11 * s**4 + q22 * c**4 + 2 * (q12 + 2 * q66) * s * s * c * c,
        (q11 + q22 - 4 * q66) * s * s * c * c + q12 * (s**4 + c**4),
    )


def main():
    out = sys.argv[1]
    N = int(sys.argv[2]) if len(sys.argv) > 2 else 20
    raw = load_raw(0.25, mesher="gmsh", generate=True, chord_divisions=N)
    E = np.asarray(raw["elements"])
    tables = json.loads(bytes(raw["tables_json"]).decode())
    a11 = np.array(
        [
            sum(q11bar(tables["materials"][n], a) * th for n, th, a in lay)
            for lay in tables["layups"]
        ]
    )[raw["layup_id"]]
    # reduced axial stiffness A11* = A11 - A12^2/A22 (free transverse contour strain; A16/A26 ignored): the thin-wall beam value
    ap = np.array(
        [
            [sum(qbar(tables["materials"][n], a)[k] * th for n, th, a in lay) for k in range(3)]
            for lay in tables["layups"]
        ]
    )[raw["layup_id"]]
    a11r = ap[:, 0] - ap[:, 2] ** 2 / ap[:, 1]
    Xs = {"oml": np.asarray(raw["nodes"], float), "mid": _inset(raw)}
    L = float(Xs["oml"][:, 2].max() - Xs["oml"][:, 2].min())
    deck_eta, K = read_deck(DECK)
    zz = Xs["oml"][:, 2]
    zmin = np.array([zz[e[e >= 0]].min() for e in E])
    zmax = np.array([zz[e[e >= 0]].max() for e in E])

    def section(X, zc, a11=a11):
        mids, ln, tv, w = [], [], [], []
        for k in np.nonzero((zmin < zc) & (zmax > zc))[0]:
            ids = E[k][E[k] >= 0]
            pts = []
            for a, b in zip(ids, np.roll(ids, -1), strict=True):
                za, zb = zz[a] - zc, zz[b] - zc
                if za * zb < 0:
                    pts.append(X[a, :2] + za / (za - zb) * (X[b, :2] - X[a, :2]))
            if len(pts) == 2:
                d = pts[1] - pts[0]
                l = float(np.linalg.norm(d))
                if l > 0:
                    mids.append(0.5 * (pts[0] + pts[1]))
                    ln.append(l)
                    tv.append(d / l)
                    w.append(a11[k] * l)
        mids, ln, tv, w = map(np.array, (mids, ln, tv, w))
        ea = w.sum()
        c = (w[:, None] * mids).sum(0) / ea
        d = mids - c
        own = w * ln**2 / 12.0
        ixx = (w * d[:, 1] ** 2).sum() + (
            own * tv[:, 1] ** 2
        ).sum()  # bending about x: curvature in y
        iyy = (w * d[:, 0] ** 2).sum() + (own * tv[:, 0] ** 2).sum()
        ixy = (w * d[:, 0] * d[:, 1]).sum() + (own * tv[:, 0] * tv[:, 1]).sum()
        ev, evec = np.linalg.eigh(np.array([[ixx, -ixy], [-ixy, iyy]]))
        ang = np.degrees(
            np.arctan2(evec[1, 1], evec[0, 1])
        )  # direction (x,y) of the max-EI (edge) axis: bending about the axis perpendicular
        return ea, ev[0], ev[1], ixx, iyy, ang

    etas = [round(0.30 + 0.05 * i, 2) for i in range(14)]  # 0.30 .. 0.95
    rows = []
    hdr = "eta   z[m]   | deck: EA  K44(ref axis)  K55 | mid: EA  EIedge  EIflap | oml: EA EIedge EIflap | ratios mid/deck(ec): EA edge flap | oml/deck: EA edge flap | principal angle(deg, mid)"
    print(hdr)
    for eta in etas:
        z = eta * L + 0.0137
        kd = [np.interp(eta, deck_eta, K[:, i, i]) for i in (2, 3, 4)]
        r = {"eta": eta, "z": eta * L, "deck_EA": kd[0], "deck_K44": kd[1], "deck_K55": kd[2]}
        for ref in ("mid", "oml"):
            ea, fl, ed, ixx, iyy, ang = section(Xs[ref], z)
            r[ref] = {"EA": ea, "EI_edge": ed, "EI_flap": fl, "Ixx": ixx, "Iyy": iyy, "ang": ang}
        k33, k34, k35 = (np.interp(eta, deck_eta, K[:, 2, i]) for i in (2, 3, 4))
        # deck K44/K55 are about the reference axis; the thin-wall tensor is about the section centroid -> compare with the
        # elastic-centre values K44 - K34^2/K33 (BeamDyn: K34 = -EA * x_ec), the frame trap of this task
        ec44, ec55 = kd[1] - k34**2 / k33, kd[2] - k35**2 / k33
        dmax, dmin = max(ec44, ec55), min(ec44, ec55)
        for ref in ("mid", "oml"):
            r[ref + "_red"] = dict(zip(("EA", "EI_flap", "EI_edge"), section(Xs[ref], z, a11r)[:3], strict=False))
        r["deck_edge"], r["deck_flap"] = dmax, dmin
        r["deck_edge_refaxis"], r["deck_flap_refaxis"] = max(kd[1], kd[2]), min(kd[1], kd[2])
        rows.append(r)
        print(
            "%.2f %6.1f | %.3e %.3e %.3e | %.3e %.3e %.3e | %.3e %.3e %.3e | %.3f %.3f %.3f | %.3f %.3f %.3f | %.1f"
            % (
                eta,
                eta * L,
                kd[0],
                kd[1],
                kd[2],
                r["mid"]["EA"],
                r["mid"]["EI_edge"],
                r["mid"]["EI_flap"],
                r["oml"]["EA"],
                r["oml"]["EI_edge"],
                r["oml"]["EI_flap"],
                r["mid"]["EA"] / kd[0],
                r["mid"]["EI_edge"] / dmax,
                r["mid"]["EI_flap"] / dmin,
                r["oml"]["EA"] / kd[0],
                r["oml"]["EI_edge"] / dmax,
                r["oml"]["EI_flap"] / dmin,
                r["mid"]["ang"],
            )
        )
    json.dump({"L": L, "N": N, "rows": rows}, open(out, "w"), indent=1)


if __name__ == "__main__":
    main()
