"""B3-3 M4: is the shell's per-station edge softness (M1) consistent with its 1E? (no FE, numpy only)

Reads the per-station table of tools/diagnose_bending_compliance_2x2.py, scales the BeamDyn deck
principal EIs by the 2x2-eigenvalue shell/deck ratio, and predicts with a clamped Euler-Bernoulli
beam (tools/s4_deck_clamped_beam_modes.beam_freqs) the change in tip compliance and in 1F/1E
(odd/tasks/validation-l3-b3-audit.md, K12.M4). Variant A: ratio 1 outside the 0.2-0.85 R window;
variant B: nearest window ratio held to root and tip.

    python3 tools/diagnose_edge_static_modal_consistency.py logs/check/b3-3/m1-025-<job>.log
"""

import re
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from s4_deck_clamped_beam_modes import LENGTH, beam_freqs  # noqa: E402
from v05_deck_consistency import REF, load_beamdyn  # noqa: E402


def tip_compliance(eta, ei):
    z = np.linspace(0.0, LENGTH, 4001)
    return np.trapezoid((LENGTH - z) ** 2 / np.interp(z / LENGTH, eta, ei), z)


def main():
    rows = [tuple(map(float, l.split())) for l in Path(sys.argv[1]).read_text().splitlines()
            if re.match(r"^\s+\d+\.\d+\s+[+-]", l)]
    t = np.array(rows)
    zr, rf, re_ = t[:, 0] / LENGTH, t[:, 3] / 100, t[:, 4] / 100
    bd = load_beamdyn(REF / "IEA-15-240-RWT_BeamDyn_blade.dat")
    eta = np.array([b[0] for b in bd])
    lam = np.array([np.sort(np.linalg.eigvalsh(K[3:5, 3:5] - np.outer(K[2, 3:5], K[2, 3:5]) / K[2, 2]))
                    for _, K, _ in bd])
    m = np.array([b[2][0, 0] for b in bd])
    f0, e0 = beam_freqs(eta, lam[:, 0], m, 400)[0], beam_freqs(eta, lam[:, 1], m, 400)[0]
    c0f, c0e = tip_compliance(eta, lam[:, 0]), tip_compliance(eta, lam[:, 1])
    print(f"deck EB: 1F {f0:.4f} 1E {e0:.4f} Hz; tip compliance flap {c0f / 1e-6:.3f} edge {c0e / 1e-6:.3f} m/MN"
          f" (n={len(t)} window stations)")
    inside = (eta >= zr.min()) & (eta <= zr.max())
    for name in ("A", "B"):
        sf = np.interp(eta, zr, rf) if name == "B" else np.where(inside, np.interp(eta, zr, rf), 0.0)
        se = np.interp(eta, zr, re_) if name == "B" else np.where(inside, np.interp(eta, zr, re_), 0.0)
        ef, ee = lam[:, 0] * (1 + sf), lam[:, 1] * (1 + se)
        f1, e1 = beam_freqs(eta, ef, m, 400)[0], beam_freqs(eta, ee, m, 400)[0]
        print(f"variant {name}: 1F {f1 / f0 - 1:+.2%}  1E {e1 / e0 - 1:+.2%}   tip compliance flap "
              f"{tip_compliance(eta, ef) / c0f - 1:+.2%}  edge {tip_compliance(eta, ee) / c0e - 1:+.2%}")


if __name__ == "__main__":
    main()
