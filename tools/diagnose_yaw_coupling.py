"""Rung 3 of the yaw diagnosis: the aeroelastic coupling, isolated.

Both sides use OUR aero model (CCBlade + NeuralFoil), so the only difference
is the coupling:

  R1 rigid  : undeformed blade   (docs/validation_data/generated/rigid_bem_check.csv)
  coupled   : the FSI campaign   (frontiersin_results_corotational_100s_mitc3fix)

  ratio = P_fsi / P_rigid -> the coupling's de-loading at each yaw.

Reference: Zhou 2025 (yaw 0) and Ma 2025 (yaw sweep), both beam FSI, i.e.
aeroelastic too — so the comparison is like-for-like in kind (coupled vs
coupled) even though the structural class differs (shell vs beam).

Usage:
  python tools/diagnose_yaw_coupling.py [--csv <out>]
"""

from __future__ import annotations

import argparse
import csv
import pathlib
import sys

import numpy as np

REPO = pathlib.Path(__file__).resolve().parents[1]
RIGID_CSV = REPO / "docs" / "validation_data" / "generated" / "rigid_bem_check.csv"
CAMPAIGN = pathlib.Path(
    "/scratch/leahk/eduardo.donestevez/frontiersin_results_corotational_100s_mitc3fix"
)
LIT_P = {0: 14.76, 10: 14.5, 20: 13.6, 30: 11.9, 40: 9.7}   # Zhou(0) + Ma(10-40)
LIT_T = {0: 2.20}
YAWS = (0, 10, 20, 30, 40)
T0, T1 = 40.0, 100.0


def rigid_table() -> dict:
    """Rigid (uncoupled) rotor values from the aero ladder artifact."""
    rows = list(csv.DictReader(RIGID_CSV.open()))
    out = {}
    for r in rows:
        if r["sweep"] != "yaw":
            continue
        out[int(float(r["yaw_deg"]))] = (float(r["power_MW"]), float(r["thrust_MN"]))
    return out


def fsi_table() -> dict:
    out = {}
    for y in YAWS:
        f = CAMPAIGN / f"yaw_{y}" / "fluid" / "bem_report.csv"
        if not f.exists():
            continue
        d = np.loadtxt(f, delimiter=",", skiprows=1)
        m = (d[:, 0] >= T0) & (d[:, 0] <= T1)
        out[y] = (float(d[m, 4].mean() / 1e6), float(d[m, 2].mean() / 1e6))
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", type=pathlib.Path, default=None)
    args = ap.parse_args()

    rigid = rigid_table()
    fsi = fsi_table()
    print("=" * 100)
    print("RUNG 3 — COUPLING ISOLATION (both sides = our CCBlade + NeuralFoil)")
    print("=" * 100)
    print(f"{'yaw':>4} | {'P_rigid':>8} {'P_fsi':>8} {'ratio':>6} {'de-load%':>9} | "
          f"{'T_rigid':>8} {'T_fsi':>8} {'ratio':>6} | {'P_lit':>6} {'fsi/lit':>8}")
    rows = []
    for y in YAWS:
        r = rigid.get(y)
        f = fsi.get(y)
        if r is None or f is None:
            print(f"{y:4d} |  (faltan datos)"); continue
        pr, tr = r
        pf, tf = f
        ratio = pf / pr
        rT = tf / tr
        pl = LIT_P.get(y, float("nan"))
        print(f"{y:4d} | {pr:8.3f} {pf:8.3f} {ratio:6.3f} {100*(1-ratio):9.2f} | "
              f"{tr:8.3f} {tf:8.3f} {rT:6.3f} | {pl:6.2f} {pf/pl:8.3f}")
        rows.append([y, pr, pf, ratio, tr, tf, rT, pl])

    if rows:
        d0 = 1 - rows[0][3]                                  # ours: 1 - P_fsi/P_rigid
        dl_own = 1 - LIT_P[0] / rows[0][1]                   # literature vs OUR rigid
        AD_RIGID = 15.815                                    # AeroDyn rigid (fair base)
        dl_ad = 1 - LIT_P[0] / AD_RIGID                      # literature vs AeroDyn rigid
        ours_ad = 1 - rows[0][2] / AD_RIGID                  # ours vs AeroDyn rigid
        print(f"yaw 0 power de-loading (vs OUR rigid 16.511):  ours {100*d0:+.1f}%  "
              f"lit {100*dl_own:+.1f}%  -> x{d0/dl_own:.1f}")
        print(f"yaw 0 power de-loading (vs AeroDyn rigid 15.815, fair base):  "
              f"ours {100*ours_ad:+.1f}%  lit {100*dl_ad:+.1f}%  -> x{ours_ad/dl_ad:.1f}")
        t0 = 1 - fsi[0][1] / rigid[0][1]
        tl = 1 - LIT_T[0] / rigid[0][1]
        print(f"yaw 0 thrust de-loading (vs our rigid):  ours {100*t0:+.1f}%  "
              f"lit {100*tl:+.1f}%  -> x{t0/tl:.1f}")

    if args.csv:
        with args.csv.open("w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["yaw_deg", "P_rigid_MW", "P_fsi_MW", "ratio", "T_rigid_MN",
                        "T_fsi_MN", "T_ratio", "P_lit_MW"])
            w.writerows(rows)
        print(f"\nwrote {args.csv}")


if __name__ == "__main__":
    main()
