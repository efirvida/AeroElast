"""W1: factor-by-factor GJ (2D St-Venant on the solid's own 2D mesh) at one eta: TE filler extent x web junction mode.
  python tools/l3_w1_factors.py ETA   (aeroenv, PYTHONPATH=$PWD). Ratios are to the deck GJ_sc (shear-free, tools/l3_solid_section_ccx.deck_metrics)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import l3_solid_section_ccx as T  # noqa: E402

eta = float(sys.argv[1])
deck = T.deck_metrics(T.deck_station(eta))[1]["GJ"]
sec = T.blade_section(eta)
base = None
print(f"eta {eta}: deck GJ_sc {deck:.4e}")
print("web_mode  te_filler_arc  filler_area   GJ_st2d      /deck   delta vs baseline(bonded, none)")
for mode in ("bonded", "through"):
    for arc in (None, 0.04, 0.06, 0.08, 0.10):
        try:
            m2 = T.build_2d(sec, 0.1, 0.03, mode, None if arc is None else {"mat": "adhesive", "arc": arc})
        except (AssertionError, KeyError) as e:
            print(f"{mode:8s} {arc} FAILED {e!r:.60}")
            continue
        A = sum(q[3] for q in m2["quads"] if q[1] == "adhesive")
        g = T.torsion_2d(sec, m2)
        base = base or g
        print(f"{mode:8s}  {arc!s:12s}  {A:.3e}  {g:.5e}  {g/deck:.3f}  {100*(g/base-1):+.1f} %")
