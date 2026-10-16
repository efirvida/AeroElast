"""V-05 deck consistency: ElastoDyn FlpStff/EdgStff vs BeamDyn 6x6 bending block, same eta.

Both official IEA-15-240-RWT decks describe the same blade; no aeroelast code is involved.
BeamDyn is reduced to the elastic centre (Schur complement of K33 into the K44/K55 block) and
to principal axes (eigenvalues), which is what ElastoDyn tabulates; ElastoDyn is linearly
interpolated onto the 26 BeamDyn stations (ElastoDyn has 50 stations, so no extrapolation).
Mass: ElastoDyn BMassDen vs BeamDyn M11.

Usage: python3 tools/v05_deck_consistency.py   (numpy only, no package import)
"""

from pathlib import Path

import numpy as np

REF = Path(__file__).resolve().parents[1] / "tests" / "IEA15MW" / "reference"


def load_beamdyn(path: Path) -> list[tuple[float, np.ndarray, np.ndarray]]:
    text = path.read_text().split("DISTRIBUTED PROPERTIES")[1].split("\n", 1)[1]
    v = [float(t) for t in text.split()]
    n = int(path.read_text().split("station_total")[0].split()[-1])
    return [(v[i], np.array(v[i + 1:i + 37]).reshape(6, 6), np.array(v[i + 37:i + 73]).reshape(6, 6))
            for i in range(0, 73 * n, 73)]


def main() -> None:
    bd = load_beamdyn(REF / "IEA-15-240-RWT_BeamDyn_blade.dat")
    ed = np.loadtxt(REF / "IEA-15-240-RWT_ElastoDyn_blade.dat", skiprows=16, max_rows=50)
    print(" eta   BD_EIflap  BD_EIedge |  ED_Flp     ED_Edg   | Flp/BD-1  Edg/BD-1  mED/mBD-1")
    for eta, K, M in bd:
        k = K[2, 3:5]
        lo, hi = np.sort(np.linalg.eigvalsh(K[3:5, 3:5] - np.outer(k, k) / K[2, 2]))
        flp, edg, med = (np.interp(eta, ed[:, 0], ed[:, c]) for c in (4, 5, 3))
        print(f"{eta:.3f} {lo:.3e} {hi:.3e} | {flp:.3e} {edg:.3e} | {100 * (flp / lo - 1):+7.2f} % "
              f"{100 * (edg / hi - 1):+7.2f} % {100 * (med / M[0, 0] - 1):+7.2f} %")


if __name__ == "__main__":
    main()
