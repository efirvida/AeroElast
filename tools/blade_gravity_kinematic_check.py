"""Second-order kinematic identity for the blade tip under gravity, from the shell's OWN linear field.

For a slender cantilever the geometric-nonlinear axial increment of the tip is the foreshortening
``-1/2 * integral((dux/dz)^2 + (duy/dz)^2) dz`` of the deflection it actually has. Evaluating it on
the shell's linear solution removes every beam-vs-shell difference in mass and stiffness (the beam
reference is not converged in its slice count: mass and tip deflection move with ``n_slices``), so
it checks the solver's geometric nonlinearity alone: ``tip_z(nonlinear) - tip_z(linear)`` of
``blade_gravity_nl_increment`` must equal this number.

Usage: python tools/blade_gravity_kinematic_check.py --element-size 0.5
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(REPO_ROOT / "tools"))

from blade_gravity_beam_shortening import foreshortening  # noqa: E402
from blade_gravity_nl_increment import GRAVITY, _solve  # noqa: E402


def slab_means(mesh, u, n_bins: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Span stations and the mean x/y translation of the nodes in each equal-width z slab."""
    coords = mesh.coords_array
    trans = np.asarray(u).reshape(coords.shape[0], -1)[:, :3]
    edges = np.linspace(coords[:, 2].min(), coords[:, 2].max() + 1e-9, n_bins + 1)
    slab = np.digitize(coords[:, 2], edges) - 1
    z, ux, uy = [], [], []
    for k in range(n_bins):
        inside = slab == k
        if inside.any():
            z.append(coords[inside, 2].mean())
            ux.append(trans[inside, 0].mean())
            uy.append(trans[inside, 1].mean())
    return np.array(z), np.array(ux), np.array(uy)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--element-size", type=float, required=True)
    ap.add_argument("--bins", type=int, nargs="+", default=[30, 60, 120])
    args = ap.parse_args()

    from aeroelast.core.mesh.generators import BladeMesh
    from aeroelast.models.blade.model import build_rust_properties
    from aeroelast.solvers.elasticity.static_linear import StaticLinearSolver

    gen = BladeMesh(
        excel_file=None,
        yaml_file=str(REPO_ROOT / "tests" / "IEA-15-240-RWT.yaml"),
        airfoil_dir=str(REPO_ROOT / "tests" / "airfoils"),
        element_size=args.element_size,
    )
    mesh = gen.generate(renumber="rcm")
    props = build_rust_properties(gen.numad_mesh_data)
    captured = {}

    class _Capture(StaticLinearSolver):
        def solve(self, *a, **k):
            captured["u"] = np.asarray(super().solve(*a, **k))
            return captured["u"]

    _solve(_Capture, mesh, props, {})
    result = {"element_size": args.element_size, "gravity": list(GRAVITY)}
    for n_bins in args.bins:
        z, ux, uy = slab_means(mesh, captured["u"], n_bins)
        result[f"dz_kinematic_bins{n_bins}"] = foreshortening(z, ux, uy)
    print(json.dumps(result, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
