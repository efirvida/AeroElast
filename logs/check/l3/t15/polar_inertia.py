"""T15 follow-up: section polar mass inertia about the mass centroid, shell vs the BeamDyn deck mass matrix.

1T ~ sqrt(GJ / I_p): a shell GJ about 10 % below the deck with 1T 1.6 % above it is consistent only if the shell I_p is
about 13-15 % below the deck. No solve: mesh + lumped nodal mass (the S-1 mass path) only.
Deck: BeamDyn mass matrix about the reference axis; xc = M[1,5]/m, yc = -M[0,5]/m, I_p,cg = M[5,5] - m (xc^2 + yc^2)
(rotation-invariant, so the deck/shell axis mapping does not enter).
"""
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / "tools"))
from run_s1_sectional import native_bay_centres  # noqa: E402

from aeroelast.core.mesh.generators import BladeMesh  # noqa: E402
from aeroelast.models.blade.model import build_rust_properties  # noqa: E402
from aeroelast.postprocess.sectional import SectionalExtractor  # noqa: E402


def deck_mass(path):
    """(eta, m, I_p about the mass centroid) per station of a BeamDyn blade file."""
    lines = [ln.split() for ln in Path(path).read_text().splitlines()]
    eta, m, ip = [], [], []
    i = 0
    while i < len(lines):
        row = lines[i]
        if len(row) == 1:
            try:
                e = float(row[0])
            except ValueError:
                i += 1
                continue
            num = []
            j = i + 1
            while len(num) < 12 and j < len(lines):
                try:
                    vals = [float(x) for x in lines[j]]
                except ValueError:
                    vals = []
                if len(vals) == 6:
                    num.append(vals)
                j += 1
            if len(num) == 12:
                M = np.array(num[6:])
                mm = M[0, 0]
                xc, yc = M[1, 5] / mm, -M[0, 5] / mm
                eta.append(e)
                m.append(mm)
                ip.append(M[5, 5] - mm * (xc ** 2 + yc ** 2))
                i = j
                continue
        i += 1
    return np.array(eta), np.array(m), np.array(ip)


def shell_mass(yaml):
    gen = BladeMesh(yaml_file=str(yaml), element_size=0.25, skin_reference="mid")
    mesh = gen.generate(renumber="rcm")
    ext = SectionalExtractor(mesh, build_rust_properties(gen.numad_mesh_data))
    memb = ext.station_membership()
    z, m, ip = [], [], []
    for st in sorted(memb):
        sl = ext._station_slice(st, memb)
        w = np.concatenate([0.5 * ext.node_mass[sl["lower"]], 0.5 * ext.node_mass[sl["upper"]],
                            ext.node_mass[sl["interior"]]])
        xyz = np.concatenate([ext.coords[sl["lower"]], ext.coords[sl["upper"]], ext.coords[sl["interior"]]])
        dz = max(sl["z1"] - sl["z0"], 1e-9)
        c = (w[:, None] * xyz).sum(axis=0) / w.sum()
        r2 = ((xyz[:, :2] - c[:2]) ** 2).sum(axis=1)
        z.append(0.5 * (sl["z0"] + sl["z1"]))
        m.append(w.sum() / dz)
        ip.append((w * r2).sum() / dz)
    return np.array(z), np.array(m), np.array(ip)


for yaml, deck in (("tests/IEA-15-240-RWT.yaml", "tests/IEA15MW/reference/IEA-15-240-RWT_BeamDyn_blade.dat"),
                   ("tests/IEA-3.4-130-RWT.yaml", sys.argv[1] if len(sys.argv) > 1 else "")):
    if not deck:
        continue
    z, m, ip = shell_mass(ROOT / yaml)
    L = z.max()
    eta_d, m_d, ip_d = deck_mass(ROOT / deck)
    eta = z / L
    md, ipd = np.interp(eta, eta_d, m_d), np.interp(eta, eta_d, ip_d)
    print(f"{yaml}: {len(z)} stations, L {L:.2f}; deck {deck} ({len(eta_d)} stations)")
    for a, b in ((0.2, 0.3), (0.3, 0.6), (0.6, 0.85), (0.2, 0.85), (0.6, 1.0)):
        k = (eta >= a) & (eta < b)
        print(f"  eta {a:.2f}-{b:.2f}: mass {100 * np.mean(m[k] / md[k] - 1):+6.2f} %  "
              f"I_p,cg {100 * np.mean(ip[k] / ipd[k] - 1):+6.2f} %  (n {int(k.sum())})")
