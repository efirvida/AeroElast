"""W1: material intervals along horizontal cuts (flap direction) through the webs: solid vs a registered SONATA export.
  python tools/l3_w1_webcut.py ETA SONATA.json GEOMPREFIX.json   (GEOM json from l3_w1_geom.py has the rigid fit)"""
import json
import sys
from pathlib import Path

from shapely import affinity
from shapely.geometry import LineString

sys.path.insert(0, str(Path(__file__).parent))
import l3_w1_geom as G  # noqa: E402


def cut(polys, y, x0=-0.9, x1=0.9):
    line = LineString([(x0, y), (x1, y)])
    out = []
    for m, g in polys:
        i = g.intersection(line)
        for q in (i.geoms if hasattr(i, "geoms") else [i]):
            if q.length > 1e-6:
                out.append((q.bounds[0], q.bounds[2], m))
    out.sort()
    merged = []
    for a, b, m in out:  # merge neighbours of the same material
        if merged and merged[-1][2] == m and a - merged[-1][1] < 1e-4:
            merged[-1][1] = b
        else:
            merged.append([a, b, m])
    return merged


eta, jf, gf = float(sys.argv[1]), sys.argv[2], sys.argv[3]
sec, m2, sol = G.solid_polys(eta)
d, son = G.sonata_polys(jf)
p = json.load(open(gf))
son = [(m, affinity.translate(affinity.rotate(g, p["rot"], origin=(0, 0)), *p["shift"])) for m, g in son]
for y in (-1.70, -1.80, -1.92, -2.05, -2.15):
    print(f"--- cut y={y}: material intervals x0..x1 (width mm)")
    for nm, pl in (("solid  ", sol), ("SONATA ", son)):
        iv = [(a, b, m) for a, b, m in cut(pl, y) if (abs(a + 0.35) < 0.08 or abs(b + 0.35) < 0.08 or abs(a - 0.35) < 0.08 or abs(b - 0.35) < 0.08)]
        print(nm, "; ".join(f"{m[:6]} [{a:.4f},{b:.4f}] {1e3*(b-a):.1f}" for a, b, m in iv))
