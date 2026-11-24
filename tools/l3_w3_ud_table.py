"""Parse logs/chord_scaled/w3/ud-*.log RESULT lines (SONATA 2021 raw TS diag: EA, shear, shear, GJ, flap, edge) -> ratio to the deck.
usage: python tools/l3_w3_ud_table.py LOG [LOG ...]"""

import re
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from l3_w3_ei_edge import DECK, read_deck  # noqa: E402

eta_d, K = read_deck(DECK)
for log in sys.argv[1:]:
    tag = ""
    for ln in open(log):
        if ln.startswith("==="):
            tag = ln.strip()
            print(tag)
        m = re.match(r"RESULT eta (\S+) TS diag (.*)", ln)
        if m:
            e = float(m.group(1))
            v = [float(x) for x in m.group(2).split()]
            d = lambda i: np.interp(e, eta_d, K[:, i, i])  # noqa: E731
            print(
                "  eta %.2f  EA %.3f  GJ %.3f  flap %.3f  edge %.3f"
                % (e, v[0] / d(2), v[3] / d(5), v[4] / d(4), v[5] / d(3))
            )
        elif ln.startswith("RESULT") and "FAILED" in ln:
            print("  eta %s FAILED" % ln.split()[2])
