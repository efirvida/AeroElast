"""S-7's beam reference: staircase integral against the integral to z (E8, 2026-10-11).

`run_s7_torsion.beam_twist_profile_deg` integrates 1/GJ only up to the last BeamDyn station at or
below each query z and drops the partial segment, so it under-states the beam twist (mean 0.84,
minimum 0.63 over the 0.3-0.9 band) and inflates the shell/beam ratio. This script runs the S-7
production case unchanged and prints the band ratio against the staircase reference, the integral
to z with K66, and the integral to z with the effective torsional stiffness 1/C66 (C = K^-1, which
frees the shear-torsion coupling of the 6x6).

    scripts/aeroenv.sh python tools/diagnose_s7_beam_reference.py <element_size>

Recorded in odd/tasks/s7-torsion-ratio.md section E8.
"""

import sys
import numpy as np
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "tools"))
sys.path.insert(0, str(REPO_ROOT))
import run_s7_torsion as s7  # noqa: E402

es = float(sys.argv[1])
T = 1e6
L = 117.0
path = s7.BEAMDYN_FILE
lines = [l.strip() for l in open(path)]
i = next(k for k, l in enumerate(lines) if "DISTRIBUTED PROPERTIES" in l) + 1
zb, k66, c66 = [], [], []
while i < len(lines):
    p = lines[i].split()
    if len(p) == 1:
        try:
            eta = float(p[0])
        except ValueError:
            i += 1
            continue
        j = i + 1
        while not lines[j]:
            j += 1
        K = np.array([[float(x) for x in lines[j + k].split()] for k in range(6)])
        zb.append(eta * L)
        k66.append(K[5, 5])
        c66.append(1 / np.linalg.inv(K)[5, 5])
        i = j + 6
    else:
        i += 1
zb, k66, c66 = map(np.asarray, (zb, k66, c66))


def integral_to(zq, gj):
    v = zb <= 111.2
    z, inv = zb[v], 1 / gj[v]
    zq = min(zq, 111.2)
    m = z <= zq
    zz = np.append(z[m], zq) if zq > z[m][-1] else z[m]
    return np.rad2deg(T * np.trapezoid(np.interp(zz, z, inv), zz))


mesh, props = s7.build_mesh_and_properties(es)
u = s7.run_torsion_case(mesh, props, T)
z7, th7 = s7.section_twists_deg(mesh, u)
stair = s7.beam_twist_profile_deg(T, z7)
full_k = np.array([integral_to(z, k66) for z in z7])
full_c = np.array([integral_to(z, c66) for z in z7])
m = (z7 >= 0.3 * L) & (z7 <= 0.9 * L)
print(f"es={es} nodes={len(mesh.coords_array)}")
print(f"  ratio staircase K66 (S-7 as published): {np.mean(th7[m] / stair[m]):.4f}")
print(f"  ratio integral-to-z K66:                {np.mean(th7[m] / full_k[m]):.4f}")
print(f"  ratio integral-to-z 1/C66:              {np.mean(th7[m] / full_c[m]):.4f}")
print(
    f"  staircase/full beam twist over band: min {np.min(stair[m] / full_k[m]):.3f} mean {np.mean(stair[m] / full_k[m]):.3f}"
)
