"""H6 prior: how much of the IEA-3.4 parked-mode deficit do the S-1 station ratios explain, and how much
of it lives in the bands where H6 acts (deck smear of the partial-span layer ends)?

Decoupled Euler-Bernoulli cantilever (flap: EI_flap, edge: EI_edge, mass per length), Hermite elements on
the deck eta scaled to the blade length. Base = deck. Perturbed = deck x (ours/ref) of the S-1 json per
station (EI_flap / EI_edge / mass), linearly interpolated on the station centres, unity outside a mask.
Scratch estimate for the written prediction (Handoff 2 item 1); not a reference.

usage: python tools/l3_beam_station_ratio.py <beamdyn.dat> <s1.json> <blade_length> [bridged]
"""

import json
import sys

import numpy as np
from scipy.linalg import eigh

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from beam_reference import load_beamdyn_blade  # noqa: E402


def cantilever_freqs(z, ei, m, n_modes=3):
    n = len(z)
    K = np.zeros((2 * n, 2 * n))
    M = np.zeros((2 * n, 2 * n))
    for e in range(n - 1):
        L = z[e + 1] - z[e]
        EI = 0.5 * (ei[e] + ei[e + 1])
        rho = 0.5 * (m[e] + m[e + 1])
        k = (
            EI
            / L**3
            * np.array(
                [
                    [12, 6 * L, -12, 6 * L],
                    [6 * L, 4 * L * L, -6 * L, 2 * L * L],
                    [-12, -6 * L, 12, -6 * L],
                    [6 * L, 2 * L * L, -6 * L, 4 * L * L],
                ]
            )
        )
        mm = (
            rho
            * L
            / 420
            * np.array(
                [
                    [156, 22 * L, 54, -13 * L],
                    [22 * L, 4 * L * L, 13 * L, -3 * L * L],
                    [54, 13 * L, 156, -22 * L],
                    [-13 * L, -3 * L * L, -22 * L, 4 * L * L],
                ]
            )
        )
        idx = slice(2 * e, 2 * e + 4)
        K[idx, idx] += k
        M[idx, idx] += mm
    w2 = eigh(K[2:, 2:], M[2:, 2:], eigvals_only=True)[:n_modes]
    return np.sqrt(w2) / (2 * np.pi)


def main():
    bd = load_beamdyn_blade(sys.argv[1], blade_length=float(sys.argv[3]))
    s1 = json.load(open(sys.argv[2]))["metrics"]
    L = float(sys.argv[3])
    zf = np.linspace(0.0, L, 841)  # fine beam grid
    zd = np.asarray(bd.z, float)
    fields = {"flap": np.asarray(bd.ei_flap, float), "edge": np.asarray(bd.ei_edge, float)}
    mass = np.asarray(bd.mass_dens, float)
    base_m = np.interp(zf, zd, mass)
    base = {k: np.interp(zf, zd, v) for k, v in fields.items()}

    def ratio(metric, mask):
        zs = np.asarray(s1[metric]["z"], float)
        r = np.asarray(s1[metric]["ours"], float) / np.asarray(s1[metric]["ref"], float)
        r = np.where(mask(zs / L), r, 1.0)
        return np.interp(zf, zs, r)

    f0 = {k: cantilever_freqs(zf, base[k], base_m, 2 if k == "flap" else 1) for k in base}
    print(f"deck beam Hz: 1F {f0['flap'][0]:.4f} 2F {f0['flap'][1]:.4f} 1E {f0['edge'][0]:.4f}")
    masks = {
        "all <=0.955 (tip extraction excluded)": lambda e: e <= 0.955,
        "all incl. tip >0.955": lambda e: e >= 0.0,
        "root 0-0.20": lambda e: e <= 0.20,
        "mid 0.20-0.85": lambda e: (e > 0.20) & (e <= 0.85),
        "0.85-0.955": lambda e: (e > 0.85) & (e <= 0.955),
        "H6 bands 0.05-0.10, 0.80-0.85, 0.90-0.955": lambda e: (
            ((e >= 0.05) & (e <= 0.10)) | ((e >= 0.80) & (e <= 0.85)) | ((e >= 0.90) & (e <= 0.955))
        ),
        "0-0.05 (no partial layer, not H6)": lambda e: e < 0.05,
    }
    for name, mk in masks.items():
        rm = ratio("mass", mk)
        out = []
        for k, metric in (("flap", "EI_flap"), ("edge", "EI_edge")):
            f = cantilever_freqs(
                zf, base[k] * ratio(metric, mk), base_m * rm, 2 if k == "flap" else 1
            )
            out.append(100 * (f / f0[k] - 1))
        print(f"{name:45s} 1F {out[0][0]:+6.2f} 1E {out[1][0]:+6.2f} 2F {out[0][1]:+6.2f} %")


if __name__ == "__main__":
    main()


def bridged(argv):
    """H6-realistic: inside each H6 band the ratio is bridged linearly between the stations just outside it
    (H6 can only remove the band's excess over its neighbours; the 0.04 and 0.15 deficits sit where no layer
    ends in either model). Prints the frequency change from 'measured ratios' to 'bridged ratios'."""
    bd = load_beamdyn_blade(argv[1], blade_length=float(argv[3]))
    s1 = json.load(open(argv[2]))["metrics"]
    L = float(argv[3])
    zf = np.linspace(0.0, L, 841)
    zd = np.asarray(bd.z, float)
    base_m = np.interp(zf, zd, np.asarray(bd.mass_dens, float))
    bands = [(0.05, 0.10), (0.80, 0.85), (0.90, 0.955)]

    def ratio(metric, bridge):
        zs = np.asarray(s1[metric]["z"], float)
        r = np.asarray(s1[metric]["ours"], float) / np.asarray(s1[metric]["ref"], float)
        keep = zs / L <= 0.955
        if bridge:
            for lo, hi in bands:
                keep &= ~((zs / L >= lo) & (zs / L <= hi))
        return np.interp(zf, zs[keep], r[keep])

    res = {}
    for bridge in (False, True):
        rm = ratio("mass", bridge)
        f = {}
        for k, metric, nm in (("flap", "EI_flap", 2), ("edge", "EI_edge", 1)):
            ei = np.interp(zf, zd, np.asarray(getattr(bd, f"ei_{k}"), float))
            f[k] = cantilever_freqs(zf, ei * ratio(metric, bridge), base_m * rm, nm)
        res[bridge] = f
    d = {
        m: 100 * (res[True][k][i] / res[False][k][i] - 1)
        for m, k, i in (("1F", "flap", 0), ("1E", "edge", 0), ("2F", "flap", 1))
    }
    print(
        "H6-realistic shift (bridged bands vs measured):",
        " ".join(f"{m} {v:+.2f}" for m, v in d.items()),
        "%",
    )


if __name__ == "__main__" and len(sys.argv) > 4 and sys.argv[4] == "bridged":
    bridged(sys.argv)
