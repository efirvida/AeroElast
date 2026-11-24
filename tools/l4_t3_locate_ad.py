"""L4 T3 locate (step 1): AeroDyn driver run with extended nodal outputs (Theta, Phi, Curve, Vx, Vy, Cl, Cd, Re, ...).
Saves last-revolution means per node to logs/l4_t3locate/ad_<key>.json. Diagnostic only."""

import json
import re
import sys
import tempfile
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import l4_t5_options as t5  # noqa: E402

ob, p = t5.ob, t5.p
EXTRA = (
    sys.argv[1].split(",")
    if len(sys.argv) > 1
    else [
        "Theta",
        "Phi",
        "Curve",
        "Vx",
        "Vy",
        "Cl",
        "Cd",
        "Re",
        "Chi",
        "Fl",
        "Fd",
    ]
)
CASES = [(15.471, 7.518, 12.19), (10.659, 7.518, 0.0)]
OUT = Path("logs/l4_t3locate")


def run(case, extra=EXTRA, tilt=None, precone=None, keep=None):
    tilt = ob.IEA15MW_SHAFT_TILT_DEG if tilt is None else tilt
    precone = ob.IEA15MW_PRECONE_DEG if precone is None else precone
    driver = p._openfast_driver_or_skip()
    with tempfile.TemporaryDirectory(prefix="l4loc-") as td:
        td = Path(td)
        primary = ob.write_bem_primary(p.PRIMARY, td / "AD_clean.dat")
        txt = primary.read_text()
        txt = re.sub(r"^Vrel\s*$", "Vrel\n" + "\n".join(extra), txt, count=1, flags=re.M)
        primary.write_text(txt)
        dvr = ob.write_aerodyn_dvr(
            td / "t.dvr",
            primary,
            [case],
            dt=0.25,
            tmax=40.0,
            precone_deg=precone,
            shaft_tilt_deg=tilt,
        )
        outs = ob.run_aerodyn_driver(driver, dvr, td)
        if keep:
            import shutil

            shutil.copy(outs[0], keep)
        lines = Path(outs[0]).read_text(errors="replace").splitlines()
    hi = next(i for i, l in enumerate(lines) if l.startswith(" Time") or "\tTime" in l)
    hd = lines[hi].split()
    rows = []
    for ln in lines[hi + 2 :]:
        f = ln.split()
        if len(f) == len(hd):
            try:
                rows.append([float(x) for x in f])
            except ValueError:
                pass
    d = np.array(rows)
    col = {n: j for j, n in enumerate(hd)}
    az = d[:, col["Azimuth"]]
    w = np.where(np.diff(az) < 0)[0]
    seg = d[w[-2] + 1 :] if len(w) >= 2 else d
    res = {}
    for n in hd:
        m = re.match(r"AB1N(\d{3})(\w+)", n)
        if m:
            res.setdefault(m.group(2), {})[int(m.group(1))] = float(seg[:, col[n]].mean())
    out = {k: [v[i] for i in sorted(v)] for k, v in res.items()}
    out["_thrust"] = float(seg[:, col["RtFldFxh"]].mean())
    out["_torque"] = float(seg[:, col["RtFldMxh"]].mean())
    return out, hd


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    for case in CASES:
        res, hd = run(case, keep=OUT / f"ad_{case[0]}.out")
        key = f"{case[0]}_{case[1]}_{case[2]}"
        (OUT / f"ad_{key}.json").write_text(json.dumps(res))
        print(key, sorted(res))
