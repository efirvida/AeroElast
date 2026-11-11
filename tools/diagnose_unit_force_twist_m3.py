"""T5c M3: S-4 modal / S-2+S-7 static of tools/chord_scaled_convergence.py with the TE-gap option B or C applied to the mesh.

Wraps ``chord_scaled_convergence.build`` so every group sees the transformed (mesh, props) (te_fill / te_adhesive of
tools/diagnose_gc_twist_chord.py). Option A is the existing logs/chord_scaled/t5/<case>/{modal,static}.json.

    scripts/aeroenv.sh python tools/diagnose_unit_force_twist_m3.py modal N15 B --out-dir logs/chord_scaled/t5c
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / "tools")]


def main() -> int:
    import chord_scaled_convergence as C
    import diagnose_gc_twist_chord as D

    ap = argparse.ArgumentParser()
    ap.add_argument("group", choices=["modal", "static"])
    ap.add_argument("case")
    ap.add_argument("option", choices=["B", "C"])
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--rpm", type=float, default=C.RPM_RATED)
    ap.add_argument("--n-modes", type=int, default=20)
    a = ap.parse_args()
    fn = {"B": D.te_fill, "C": D.te_adhesive}[a.option]
    orig = C.build

    def build(case):
        mesh, props, t = orig(case)
        mesh, props = fn(mesh, props)
        return mesh, props, t

    C.build = build
    res = C.GROUPS[a.group](a.case, a)
    a.out_dir.mkdir(parents=True, exist_ok=True)
    (a.out_dir / f"m3-{a.case}-{a.option}-{a.group}.json").write_text(json.dumps({"case": a.case, "option": a.option, "group": a.group, **res}, indent=1, default=float))
    print("done", a.group, a.case, a.option)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
