"""B3-3 K13.D: S-1 sectional GJ at mid, per station and as a window integral (extractor dispersion?).

Runs tools/run_s1_sectional.ea_gj_from_static_response as the S-1 test does (mid skin), captures the
twist profile theta(z) that feeds its smoothing spline, and compares: (1) the per-station GJ against
K66 and against 1/C66 of the BeamDyn deck; (2) the spline-free window estimate
GJ_win = T dz / (theta_b - theta_a) against the harmonic mean of each reference over the same window
(odd/tasks/validation-l3-b3-audit.md, K13.D).

    scripts/aeroenv.sh python -u -B tools/diagnose_s1_gj_station.py <element_size>
"""
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))
sys.path.insert(0, str(ROOT))

import run_s1_sectional as rs  # noqa: E402
from beam_reference import load_beamdyn_blade  # noqa: E402
from aeroelast.core.mesh.generators import BladeMesh  # noqa: E402
from aeroelast.models.blade.model import build_rust_properties  # noqa: E402

TORQUE = 1e6


def deck_gj(path):
    """(z, K66, 1/C66) from the BeamDyn blade file, as tools/diagnose_s7_beam_reference.py."""
    from v05_deck_consistency import load_beamdyn
    st = load_beamdyn(path)
    length = 117.0
    return (np.array([e * length for e, _, _ in st]), np.array([k[5, 5] for _, k, _ in st]),
            np.array([1 / np.linalg.inv(k)[5, 5] for _, k, _ in st]))


def stats(x):
    return f"signed {100 * np.mean(x):+6.2f}  abs {100 * np.mean(np.abs(x)):5.2f}  std {100 * np.std(x):5.2f}  frac>0 {np.mean(x > 0):.2f}"


def main():
    size = float(sys.argv[1])
    captured = []
    original = rs.UnivariateSpline

    def spy(x, y, *a, **k):
        captured.append((np.asarray(x), np.asarray(y), k.get("s")))
        return original(x, y, *a, **k)

    rs.UnivariateSpline = spy
    gen = BladeMesh(yaml_file=str(ROOT / "tests/IEA-15-240-RWT.yaml"), element_size=size, skin_reference="mid")
    mesh = gen.generate(renumber="rcm")
    props = build_rust_properties(gen.numad_mesh_data)
    z, _, gj, _ = rs.ea_gj_from_static_response(mesh, props, torque=TORQUE, numad_data=gen.numad_mesh_data)
    theta = next(y for x, y, s in captured if s == 1e-3)
    zb, k66, c66 = deck_gj(ROOT / "tests/IEA15MW/reference/IEA-15-240-RWT_BeamDyn_blade.dat")
    bd = load_beamdyn_blade(ROOT / "tests/IEA15MW/reference/IEA-15-240-RWT_BeamDyn_blade.dat")
    m = (z > 0.2 * z.max()) & (z < 0.85 * z.max())
    print(f"S1-GJ-STATION mid es={size} stations in window {m.sum()}")
    for label, ref in (("K66", k66), ("1/C66", c66)):
        e = gj[m] / np.interp(z[m], zb, ref) - 1
        print(f"  per-station vs {label:6s} {stats(e)}")
    e = gj[m] / np.interp(z[m], bd.z, bd.gj) - 1
    print(f"  per-station vs bd.gj (test ref) {stats(e)}   [test criterion: abs < 5 %]")
    idx = np.nonzero(m)[0]
    a, b = idx[0], idx[-1]
    zz = np.linspace(z[a], z[b], 2001)
    for label, ref in (("K66", k66), ("1/C66", c66)):
        gj_ref = (z[b] - z[a]) / np.trapezoid(1 / np.interp(zz, zb, ref), zz)
        gj_win = TORQUE * (z[b] - z[a]) / abs(theta[b] - theta[a])
        print(f"  window z {z[a]:.2f}-{z[b]:.2f} m: GJ_win {gj_win:.4e}  harmonic {label} {gj_ref:.4e}  ratio-1 {100 * (gj_win / gj_ref - 1):+6.2f} %")
    for step in (4, 8):
        loc = [(0.5 * (z[i] + z[i + step]), TORQUE * (z[i + step] - z[i]) / abs(theta[i + step] - theta[i])) for i in idx if i + step <= b]
        r = np.array([g / np.interp(zc, zb, k66) - 1 for zc, g in loc])
        print(f"  chord over {step} stations vs K66 (local pointwise): {stats(r)}")
    print("  z[m]  GJ/K66-1[%]  GJ/(1/C66)-1[%]")
    for i in idx:
        print(f"  {z[i]:6.2f} {100 * (gj[i] / np.interp(z[i], zb, k66) - 1):+9.2f} {100 * (gj[i] / np.interp(z[i], zb, c66) - 1):+12.2f}")


if __name__ == "__main__":
    main()
