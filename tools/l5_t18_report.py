"""T18 report: per-strip (0.5 m, ``bin_rows``) T6-mean vs T3 resultants, banded summary, Cm/arm split.  Reads the T18 jsons.

Writes ``<dir>/strip_table.csv`` (every strip: z, Fx edge, Fy flap, Mz for t6 and t3, their Cm-free Mz) and prints/writes
``<dir>/report.md``.  Mz = strip moment about the BeamDyn reference axis (+z, nose-up positive); ``T`` = section torsion
accumulated from the tip (strip Mz + the arm of the outer strips' F about the local axis).
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

BANDS = [(0, 20), (20, 40), (40, 60), (60, 80), (80, 100), (100, 117.5)]


def grid(tab, z):
    return {k: np.interp(z, tab["z"], tab[k], left=0.0, right=0.0) if k in tab else np.zeros(len(z)) for k in ("Fx_edge", "Fy_flap", "Mz_torsion", "T_cum_from_tip")}


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--dir", type=Path, required=True)
    args = ap.parse_args()
    L = json.loads((args.dir / "loads.json").read_text())["tables"]
    R = json.loads((args.dir / "rebuild.json").read_text())["tables"]
    z = np.array(L["t6"]["z"])
    sets = {"t6": L["t6"], "t3": L["aerodyn_wall_az0"], "t3_nocm": L["aerodyn_wall_az0_nocm"], "t6_like_nocm": L["windio_minnorm_azmean_nocm"],
            "t6_like": R["prod_rebuilt_dumptheta_mean"]}
    G = {k: grid(v, z) for k, v in sets.items()}
    # same-bin check for the zero-padded interpolation
    same = {k: len(v["z"]) == len(z) and np.allclose(v["z"], z) for k, v in sets.items()}
    lines = ["# T18 per-strip report", "", f"strip grids identical to t6's (0.5 m bins): {same}", ""]
    hdr = "| band z [m] | dFy flap (t6-t3) [kN] | dFx edge [kN] | Fy t6 | Fy t3 | Fx t6 | Fx t3 | Mz t6 | Mz t3 | dMz | of which Cm part | arm part |"
    lines += [hdr, "|" + "---|" * 12]
    for lo, hi in BANDS:
        m = (z >= lo) & (z < hi)
        s = lambda k, f: G[k][f][m].sum() / 1e3  # noqa: E731
        dmz = s("t6", "Mz_torsion") - s("t3", "Mz_torsion")
        cm_t6 = (s("t6", "Mz_torsion") - s("t6_like_nocm", "Mz_torsion"))  # t6 total - (windio, nocm, azmean): Cm part of t6
        cm_t3 = (s("t3", "Mz_torsion") - s("t3_nocm", "Mz_torsion"))
        arm = s("t6_like_nocm", "Mz_torsion") - s("t3_nocm", "Mz_torsion")
        lines.append(f"| {lo}-{hi} | {s('t6','Fy_flap')-s('t3','Fy_flap'):+.1f} | {s('t6','Fx_edge')-s('t3','Fx_edge'):+.1f} | {s('t6','Fy_flap'):.1f} | "
                     f"{s('t3','Fy_flap'):.1f} | {s('t6','Fx_edge'):.1f} | {s('t3','Fx_edge'):.1f} | {s('t6','Mz_torsion'):.1f} | {s('t3','Mz_torsion'):.1f} | "
                     f"{dmz:+.1f} | {cm_t6 - cm_t3:+.1f} | {arm:+.1f} |")
    tot = lambda k, f: G[k][f].sum() / 1e3  # noqa: E731
    lines += ["", f"totals [kN, kN*m]: Fy t6 {tot('t6','Fy_flap'):.1f} t3 {tot('t3','Fy_flap'):.1f}; Fx t6 {tot('t6','Fx_edge'):.1f} t3 {tot('t3','Fx_edge'):.1f}; "
              f"Mz t6 {tot('t6','Mz_torsion'):.1f} t3 {tot('t3','Mz_torsion'):.1f}; T_root t6 {G['t6']['T_cum_from_tip'][0]/1e3:.1f} t3 {G['t3']['T_cum_from_tip'][0]/1e3:.1f}",
              "(Cm part = strip Mz minus the same load with Mp = None; 'arm part' = difference of the Mp = None resultants, t6-like windIO azimuth mean minus T3)", ""]
    # section torsion at spans
    lines += ["| z [m] | T(z) t6 [kN*m] | T(z) t3 | ratio |", "|---|---|---|---|"]
    for zz in (20, 40, 60, 80, 100, 110):
        a, b = np.interp(zz, z, G["t6"]["T_cum_from_tip"]) / 1e3, np.interp(zz, z, G["t3"]["T_cum_from_tip"]) / 1e3
        lines.append(f"| {zz} | {a:.1f} | {b:.1f} | {a / b:.2f} |")
    (args.dir / "report.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    with open(args.dir / "strip_table.csv", "w") as fh:
        fh.write("z,Fx_t6,Fy_t6,Mz_t6,Fx_t3,Fy_t3,Mz_t3,Mz_t3_nocm,Mz_t6like_nocm,T_t6,T_t3\n")
        for i, zz in enumerate(z):
            fh.write(",".join(f"{v:.6g}" for v in (zz, G["t6"]["Fx_edge"][i], G["t6"]["Fy_flap"][i], G["t6"]["Mz_torsion"][i], G["t3"]["Fx_edge"][i],
                                                G["t3"]["Fy_flap"][i], G["t3"]["Mz_torsion"][i], G["t3_nocm"]["Mz_torsion"][i], G["t6_like_nocm"]["Mz_torsion"][i],
                                                G["t6"]["T_cum_from_tip"][i], G["t3"]["T_cum_from_tip"][i])) + "\n")


if __name__ == "__main__":
    main()
