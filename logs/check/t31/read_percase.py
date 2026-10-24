"""usage: read_percase.py [logs/check/t31]  -> one table from pc-*.log"""
import re, sys, glob
d = sys.argv[1] if len(sys.argv) > 1 else "logs/check/t31"
rows = []
for f in sorted(glob.glob(d + "/pc-*.log")):
    t = open(f, errors="replace").read()
    h = re.search(r"PERCASE tree=(\S+) case=(\S+) es=(\S+)", t)
    if not h: rows.append((f.split("/")[-1], "no result (running/failed)")); continue
    g = re.search(r"PERCASE_GAP\S* edge=(\S+) flap=(\S+) axial=(\S+) twist=(\S+) W=(\S+)", t)
    a = re.search(r"PERCASE_AE .*W=(\S+)", t); c = re.search(r"PERCASE_CCX .*W=(\S+)", t)
    rows.append((*h.groups(), *g.groups(), a.group(1), c.group(1)))
print("| tree | case | es | gap edge % | gap flap % | gap axial % | gap twist % | gap W % | W_AE | W_CCX |\n|---|---|---|---|---|---|---|---|---|---|")
for r in rows: print("| " + " | ".join(r) + " |" if len(r) == 10 else f"| {r[0]} | {r[1]} |")
