"""L5 T9c verdict table: per station, shell (fit bias removed) and deck 6x6 against the independent layered solid (CalculiX, level B).

Agreement tolerance on the principal bending angle: |psi_model - psi_solid| <= TOL_PSI_AGREE = 1.0 deg, fixed as twice the 0.5 deg step
tolerance of the shell measurement (tools/l5_t9_arbiter.py TOL_PSI); the solid A -> B change is < 0.15 deg.  The cross-term difference
C_fe/C_ff is printed next to it; the stiffness ratios are model / solid.  Nothing here widens a band: a row outside it is "no".

  scripts/aeroenv.sh env PYTHONPATH=$PWD python tools/l5_t9_verdict.py logs/check/l5_oneway/t9c [--solid z45=solid/z45_B_n4 z80=... z95=...]
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

TOL_PSI_AGREE = 1.0  # deg


def compare(model, ref):
    """Differences / ratios of one metrics dict ``model`` against ``ref`` (keys psi_deg, ratio_fe_ff, GJ, EI_flap, EI_edge)."""
    return {"dpsi": model["psi_deg"] - ref["psi_deg"], "dfe": model["ratio_fe_ff"] - ref["ratio_fe_ff"], "GJ": model["GJ"] / ref["GJ"],
            "EI_flap": model["EI_flap"] / ref["EI_flap"], "EI_edge": model["EI_edge"] / ref["EI_edge"],
            "agrees_psi": abs(model["psi_deg"] - ref["psi_deg"]) <= TOL_PSI_AGREE}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("root", type=Path)
    ap.add_argument("--solid", nargs="+", default=["45=solid/z45_B_n4", "80=solid/z80_B_n4", "95=solid/z95_B_n4"])
    a = ap.parse_args()
    deck = {r["z"]: r for r in json.loads((a.root / "stations/deck.json").read_text())["stations"]}
    shell = {r["z"]: r for r in json.loads((a.root / "stations/shell.json").read_text())["stations"]}
    lines = ["| z | model | psi deg | dpsi vs solid | agrees (<= 1 deg) | C_fe/C_ff | dC_fe/C_ff | GJ/solid | EI_flap/solid | EI_edge/solid |", "|" + "---|" * 10]
    for item in a.solid:
        z, path = item.split("=")
        z = float(z)
        s = json.loads((a.root / path / "solid.json").read_text())
        lines.append(f"| {z:g} | solid B | {s['psi_deg']:.2f} | - | - | {s['ratio_fe_ff']:+.4f} | - | 1 | 1 | 1 |")
        for label, m in (("shell (bias removed)", shell[z]["corrected"]), ("deck 6x6", deck[z])):
            c = compare(m, s)
            lines.append(f"| {z:g} | {label} | {m['psi_deg']:.2f} | {c['dpsi']:+.2f} | {'yes' if c['agrees_psi'] else 'no'} | {m['ratio_fe_ff']:+.4f} | "
                         f"{c['dfe']:+.4f} | {c['GJ']:.3f} | {c['EI_flap']:.3f} | {c['EI_edge']:.3f} |")
    text = "\n".join(lines)
    (a.root / "verdict.md").write_text(text + "\n")
    print(text)


if __name__ == "__main__":
    main()
