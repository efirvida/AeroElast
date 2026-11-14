import sys, collections, numpy as np
sys.path.insert(0, "tools")
import l3_solid_section_ccx as T
for eta in map(float, sys.argv[1].split(",")):
    sec = T.blade_section(eta)
    for mode in ("bonded", "through"):
        m = T.build_2d(sec, 0.1, 0.03, mode)
        xy = np.asarray(m["xy"]); A = collections.defaultdict(float)
        for q in m["quads"]:
            p = xy[q[0][:4], :2]; x, y = p[:, 0], p[:, 1]
            A[str(q[1])] += 0.5 * abs(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))
        print(f"eta {eta} {mode} total {sum(A.values()):.5e}")
        for k in sorted(A): print(f"   {k:28s} {A[k]:.5e}")
