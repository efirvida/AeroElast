"""T14 diagnostic: S-1 v2 station membership / mass on the gmsh mesh + independent strip-mass cross-check.

usage: scripts/aeroenv.sh python tools/l5_s1_extractor_diag.py   (V2_BLADE_* env as the S-1 test; the windio2 mesher is retired, gmsh-official-mesher S9)
"""
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "tools"))
from tests.support.paths import REPO_ROOT  # noqa: E402
from tests.support.v2_blade import build_v2_blade  # noqa: E402
from tests.support.v2_sectional import V2SectionalExtractor, bay_centres_v2, ring_ids_for_v2_mesh  # noqa: E402
from beam_reference import load_beamdyn_blade  # noqa: E402

bd = load_beamdyn_blade(REPO_ROOT / "tests/IEA15MW/reference/IEA-15-240-RWT_BeamDyn_blade.dat")
mesh, props = build_v2_blade(0.25)
info = mesh.v2_info
print("mesh", sys.argv[1:], "nodes", mesh.node_count, "elements", mesh.elements_count, "key", info["cache_key"])
ring = ring_ids_for_v2_mesh(mesh)
idx = mesh.node_id_to_index
nr = int(ring.max()) + 1
print("n_rings", nr, "ring node counts min/max", np.bincount(ring)[np.bincount(ring)>0].min(), np.bincount(ring).max())
span = np.array([np.ptp(ring[[idx[n] for n in e.node_ids]]) for e in mesh.elements])
print("element ring-span histogram", dict(zip(*np.unique(span, return_counts=True))))
ext = V2SectionalExtractor(mesh, props, ring_ids=ring)
mem = ext.station_membership()
stations = bay_centres_v2()
sec = ext.extract(stations=stations)
print("n stations", len(mem), "selected", len(sec.z))
m = (sec.z > 0.2 * sec.z.max()) & (sec.z < 0.85 * sec.z.max())
mb = np.interp(sec.z, bd.z, bd.mass_dens)
print("mass mid mean|rel|", np.mean(np.abs(sec.mass_dens[m] / mb[m] - 1)))
# independent: element mass from layup rho*t*area, strip by element centroid z
raw = None
from tests.support import v2_blade as vb
raw = vb._MEMO[info["cache_key"]]
import json
tab = json.loads(bytes(raw["tables_json"]).decode())
lm = np.array([sum(t * tab["materials"][n]["rho"] for n, t, _ in L) for L in tab["layups"]])  # kg/m2
xr = np.asarray(raw["nodes"], float)
el = raw["elements"]
emass = np.empty(len(el)); cz = np.empty(len(el))
for k, e in enumerate(el):
    ids = e[e >= 0]
    cz[k] = xr[ids, 2].mean()
    emass[k] = lm[raw["layup_id"][k]] * raw["areas"][k]
print("total shell mass (raw areas*rho*t) kg", emass.sum(), " BeamDyn int m dz", np.trapezoid(bd.mass_dens, bd.z))
print(" z      m_ext   m_strip  m_bd   ratio_ext/bd ratio_strip/bd  nnodes elems")
edges = np.linspace(0, xr[:, 2].max(), 0)  # placeholder
for i in np.linspace(0, len(sec.z) - 1, 10).astype(int):
    z0 = sec.z[i]; half = 0.5 * np.median(np.diff(stations))
    sel = (cz > z0 - half) & (cz < z0 + half)
    ms = emass[sel].sum() / (2 * half)
    st = min(mem, key=lambda s: abs(mem[s]["z"] - z0))
    print(f"{z0:7.2f} {sec.mass_dens[i]:9.1f} {ms:9.1f} {mb[i]:9.1f} {sec.mass_dens[i]/mb[i]:8.3f} {ms/mb[i]:8.3f} {len(mem[st]['nodes']):6d} {sel.sum()}")
