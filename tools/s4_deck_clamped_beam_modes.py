"""Clamped-root Euler-Bernoulli modes of the official IEA-15 decks, parked (0 rpm).

Independent check of what the S-4 OpenFAST arbiter measures: one 1D beam per deck, flap and
edge uncoupled in principal axes (structural twist ignored), cubic Hermite elements on a
uniform grid over the 117 m flexible length, properties linearly interpolated per element
midpoint. ElastoDyn: FlpStff/EdgStff/BMassDen. BeamDyn: principal values of the elastic-centre
bending block (Schur complement of K33) and M11.

Usage: python3 tools/s4_deck_clamped_beam_modes.py [n_elements]   (numpy/scipy only)
"""

import sys
from pathlib import Path

import numpy as np
from scipy.linalg import eigh

sys.path.insert(0, str(Path(__file__).resolve().parent))
from v05_deck_consistency import REF, load_beamdyn  # noqa: E402

LENGTH = 117.0  # m, TipRad 120 - HubRad 3


def beam_freqs(eta, ei, m, n_el: int, n_modes: int = 3) -> np.ndarray:
    x = np.linspace(0.0, LENGTH, n_el + 1)
    nd = 2 * (n_el + 1)
    K = np.zeros((nd, nd))
    M = np.zeros((nd, nd))
    for e in range(n_el):
        h = x[e + 1] - x[e]
        s = 0.5 * (x[e] + x[e + 1]) / LENGTH
        EI, mu = np.interp(s, eta, ei), np.interp(s, eta, m)
        ke = EI / h**3 * np.array([[12, 6 * h, -12, 6 * h], [6 * h, 4 * h * h, -6 * h, 2 * h * h],
                                   [-12, -6 * h, 12, -6 * h], [6 * h, 2 * h * h, -6 * h, 4 * h * h]])
        me = mu * h / 420 * np.array([[156, 22 * h, 54, -13 * h], [22 * h, 4 * h * h, 13 * h, -3 * h * h],
                                      [54, 13 * h, 156, -22 * h], [-13 * h, -3 * h * h, -22 * h, 4 * h * h]])
        idx = slice(2 * e, 2 * e + 4)
        K[idx, idx] += ke
        M[idx, idx] += me
    w2 = eigh(K[2:, 2:], M[2:, 2:], eigvals_only=True, subset_by_index=[0, n_modes - 1])
    return np.sqrt(w2) / (2 * np.pi)


def main() -> None:
    n_el = int(sys.argv[1]) if len(sys.argv) > 1 else 400
    ed = np.loadtxt(REF / "IEA-15-240-RWT_ElastoDyn_blade.dat", skiprows=16, max_rows=50)
    bd = load_beamdyn(REF / "IEA-15-240-RWT_BeamDyn_blade.dat")
    eta_bd = np.array([b[0] for b in bd])
    lam = []
    for _, K, _ in bd:
        k = K[2, 3:5]
        lam.append(np.sort(np.linalg.eigvalsh(K[3:5, 3:5] - np.outer(k, k) / K[2, 2])))
    lam = np.array(lam)
    m_bd = np.array([b[2][0, 0] for b in bd])
    for name, eta, flp, edg, m in (("ElastoDyn", ed[:, 0], ed[:, 4], ed[:, 5], ed[:, 3]),
                                   ("BeamDyn", eta_bd, lam[:, 0], lam[:, 1], m_bd)):
        f_flp = beam_freqs(eta, flp, m, n_el)
        f_edg = beam_freqs(eta, edg, m, n_el)
        mass = np.trapezoid(np.interp(np.linspace(0, 1, 2001), eta, m), dx=LENGTH / 2000)
        print(f"RESULT {name:9s} n_el={n_el} mass={mass:.0f} kg  1F {f_flp[0]:.4f}  1E {f_edg[0]:.4f}  "
              f"2F {f_flp[1]:.4f}  2E {f_edg[1]:.4f} Hz")


if __name__ == "__main__":
    main()
