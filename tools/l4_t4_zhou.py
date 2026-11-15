"""L4 T4: Zhou et al. 2025 (Energy 336:138488, Table 6) rigid means 16.11 MW / 2.53 MN versus our BEM
and AeroDyn at Zhou's stated operating point: 10.59 m/s, 7.55 rpm, cone and shaft tilt zero
(odd/tasks/composite-bend-twist-verdict.md, quoting the paper); pitch, density, tower, shear not stated.
Measurement instrument only. Output JSON under logs/l4_t4t5/.
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from l4_t5_options import ad_run, ours, p, ob  # noqa: E402

ZHOU_P, ZHOU_T = 16.11e6, 2.53e6
V, RPM = 10.59, 7.55


def e(x, ref):
    return (x - ref) / ref


def row(d):
    return dict(
        P_MW=d["power"] / 1e6,
        T_MN=d["thrust"] / 1e6,
        Q_MNm=d["torque"] / 1e6,
        eP_vs_Zhou=e(d["power"], ZHOU_P),
        eT_vs_Zhou=e(d["thrust"], ZHOU_T),
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", type=Path, required=True)
    a = ap.parse_args()
    out = {"zhou": dict(P_W=ZHOU_P, T_N=ZHOU_T, V=V, rpm=RPM, cone=0, tilt=0)}
    yb = p.load_blade_aero(str(p.YAML_BLADE))
    deck_tilt, deck_cone = ob.IEA15MW_SHAFT_TILT_DEG, ob.IEA15MW_PRECONE_DEG
    for tag, rpm, tilt, cone in (
        ("zhou_point", RPM, 0.0, 0.0),
        ("zhou_rpm_deck_geometry", RPM, deck_tilt, deck_cone),
        ("rated_7.56_zero_geometry", 7.56, 0.0, 0.0),
    ):
        case = (V, rpm, 0.0)
        om = rpm * 2 * np.pi / 60
        blade, ref = ad_run(case, tilt=tilt, precone=cone)
        ad = dict(thrust=ref.thrust, torque=ref.torque, power=ref.torque * om)
        o_deck = ours(blade, case, tilt=tilt, precone=cone)
        o_yaml = ours(yb, case, tilt=tilt, precone=cone)
        out[tag] = dict(
            AeroDyn=row(ad),
            ours_deck_polars=row(o_deck),
            ours_yaml_production=row(o_yaml),
            ours_deck_vs_AeroDyn={k: e(o_deck[k], ad[k]) for k in ad},
            ours_yaml_vs_AeroDyn={k: e(o_yaml[k], ad[k]) for k in ad},
        )
    # pitch is not stated by Zhou: sensitivity of ours (yaml, zero geometry)
    out["pitch_sensitivity_ours_yaml"] = {
        f"{pt:+.1f}": row(ours(yb, (V, RPM, pt), tilt=0.0, precone=0.0))
        for pt in (-1.0, -0.5, 0.0, 0.5, 1.0)
    }
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(out, indent=2, default=float) + "\n")
    for k, v in out.items():
        print(k, json.dumps(v, default=lambda x: round(x, 5)))


if __name__ == "__main__":
    main()
