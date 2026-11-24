"""W4/T5: build (and cache) the gmsh chord-scaled v2 blade mesh for N = argv[1], AR 1 (sbatch entry)."""

import sys
import time

from tests.support.v2_blade import build_v2_blade

n = int(sys.argv[1])
t0 = time.time()
mesh, props = build_v2_blade(0.25, skin_reference="mid", mesher="gmsh", chord_divisions=n)
print(
    f"N={n} nodes={mesh.node_count} elements={mesh.elements_count} cache={mesh.v2_info['cache_key'][:24]} t={time.time() - t0:.0f}s"
)
