"""Second-order beam reference for the geometric-nonlinear increment of the blade under gravity.

For a slender cantilever the first geometric nonlinearity of the tip is the axial foreshortening
``dz = -1/2 * integral((dux/dz)^2 + (duy/dz)^2) dz``, computed from the linear
Euler-Bernoulli ``BeamReference`` deflection (official ElastoDyn/BeamDyn properties). It is the
independent reference for ``tip_z(nonlinear) - tip_z(linear)`` of the shell solvers, valid while the
tip slope is small (here ~0.02 rad, so the neglected terms are O(slope^2) ~ 0.05 %).

Usage: python tools/blade_gravity_beam_shortening.py [--element-size 1.0]
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "tools"))

from beam_reference import BeamReference, load_beamdyn_blade, load_elastodyn_blade  # noqa: E402


def foreshortening(z: np.ndarray, ux: np.ndarray, uy: np.ndarray) -> float:
    """Tip axial shortening of a deflected inextensible line, to second order in the slope."""
    slope_x, slope_y = np.gradient(ux, z), np.gradient(uy, z)
    return float(-0.5 * np.trapezoid(slope_x**2 + slope_y**2, z))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--element-size", type=float, default=1.0)
    ap.add_argument("--n-slices", type=int, default=60)
    args = ap.parse_args()

    from aeroelast.core.mesh.generators import BladeMesh

    ref = REPO_ROOT / "tests" / "IEA15MW" / "reference"
    ed = load_elastodyn_blade(ref / "IEA-15-240-RWT_ElastoDyn_blade.dat")
    bd = load_beamdyn_blade(ref / "IEA-15-240-RWT_BeamDyn_blade.dat")
    gen = BladeMesh(
        excel_file=None,
        yaml_file=str(REPO_ROOT / "tests" / "IEA-15-240-RWT.yaml"),
        airfoil_dir=str(REPO_ROOT / "tests" / "airfoils"),
        element_size=args.element_size,
    )
    mesh = gen.generate(renumber="rcm")
    beam = BeamReference.from_mesh(mesh, ed, bd, n_slices=args.n_slices)
    u = beam.solve_gravity((0.0, -9.81, 0.0)).reshape(-1, 6)
    z = beam.z
    tip = u[-1, :3]
    print(f"beam linear tip  ux={tip[0]:+.4f} uy={tip[1]:+.4f} uz={tip[2]:+.4f} m")
    print(f"beam foreshortening (second order) dz = {foreshortening(z, u[:, 0], u[:, 1]):+.5f} m")


if __name__ == "__main__":
    main()
