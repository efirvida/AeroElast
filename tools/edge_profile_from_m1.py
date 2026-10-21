"""B3-3 K13.C: edge-EI ratio profile (eta, shell/deck) from a tools/diagnose_bending_compliance_2x2.py log.

Variant A: ratio 1 outside the 0.2-0.85 R window; variant B: nearest window ratio held to root and tip
(the same two variants as tools/diagnose_edge_static_modal_consistency.py, K12.M4).

    python3 tools/edge_profile_from_m1.py <m1 log> <A|B> <out.txt>
"""
import re
import sys
from pathlib import Path

import numpy as np

LENGTH = 117.0  # tools/s4_deck_clamped_beam_modes.LENGTH


def main():
    log, variant, out = Path(sys.argv[1]), sys.argv[2], Path(sys.argv[3])
    rows = np.array([tuple(map(float, l.split())) for l in log.read_text().splitlines()
                     if re.match(r"^\s+\d+\.\d+\s+[+-]", l)])
    zr, ratio = rows[:, 0] / LENGTH, 1.0 + rows[:, 4] / 100
    eta = np.linspace(0.0, 1.0, 101)
    held = np.interp(eta, zr, ratio)
    inside = (eta >= zr.min()) & (eta <= zr.max())
    prof = held if variant == "B" else np.where(inside, held, 1.0)
    out.write_text("".join(f"{e:.6f} {p:.8f}\n" for e, p in zip(eta, prof)))
    print(f"variant {variant}: {len(rows)} window stations, ratio min {prof.min():.4f} max {prof.max():.4f}")


if __name__ == "__main__":
    main()
