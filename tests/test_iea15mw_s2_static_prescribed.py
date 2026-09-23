"""S-2 — IEA 15 MW blade: static prescribed-load validation (shell vs 1D beam).

Compares AeroElast (MITC3/4 shell) against a 3D Euler-Bernoulli beam reference
built from the SAME mesh geometry (section centroids + measured twist) with
stiffness/mass from the official ElastoDyn + BeamDyn distributed properties
(tests/IEA15MW/reference/).

Load cases (blade clamped at root, span along +Z):
  LC1 gravity          (0, -9.81, 0)
  LC2 tip load 1 MN    (0, +1, 0)    flapwise
  LC3 tip load 1 MN    (+1, 0, 0)    edgewise-ish (twisted blade)
  LC4 distributed 5 kN/m (0, +1, 0)

Two blade inputs are compared:

  official (WindIO, IEA-15-240-RWT.yaml) — the official layup from the IEA
    15 MW definition (TE reinforcement ≈ 30 mm mid-span).  All four cases
    close within ~8% of the beam references (2026-09-08):

        LC1 uy +0.8% | LC2 uy +1.4% | LC3 ux +8.1% | LC4 uy -2.5%

  utd (NuMAD_utd_iea15mw.xlsx, UTD AIAA 2023-2093) — independent re-modeling.
    Its TE reinforcement is 30-75% thinner than the official layup between
    10-75 m span, which reduces the outboard edgewise stiffness: LC3 measures
    +40-55% vs the beam references.  Kept as an input-sensitivity record;
    the UTD model is NOT used for the closing claims.

Marks: slow (blade mesh generation + 4 static solves per blade).
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

petsc4py = pytest.importorskip("petsc4py", reason="PETSc not available")

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "tools"))

from beam_reference import load_beamdyn_blade, load_elastodyn_blade  # noqa: E402

GRAVITY = (0.0, -9.81, 0.0)
TIP_FORCE_N = 1.0e6
DISTRIB_Q_NM = 5.0e3

# Per-source expectations (measured 2026-09-08, element_size 0.5).
_BLADES = {
    "official": {"rel_tol_dominant": 0.05, "lc3_ratio": (0.9, 1.2)},
    "utd": {"rel_tol_dominant": 0.05, "lc3_ratio": (1.0, 2.2)},
}

_DOMINANT = {
    "LC1_gravity": {"component": "uy", "expect_negative": True},
    "LC2_tip_flap": {"component": "uy"},
    "LC4_distributed": {"component": "uy"},
}

# Cache of expensive per-blade results shared across the test session.
_S2_CACHE: dict = {}


@pytest.mark.slow
@pytest.mark.parametrize("blade_name", ["official", "utd"])
class TestS2StaticPrescribed:
    @pytest.fixture
    def s2_results(self, blade_name, iea_blade_yaml, iea_blade_xlsx):
        if blade_name in _S2_CACHE:
            return _S2_CACHE[blade_name]

        from beam_reference import BeamReference  # noqa: E402
        from run_s2_static_cases import (  # noqa: E402
            beam_tip_displacement,
            run_shell_case,
            shell_tip_displacement,
        )
        from aeroelast.core.mesh.generators import BladeMesh
        from aeroelast.models.blade.model import build_rust_properties

        ref_dir = Path(__file__).resolve().parents[1] / "tests" / "IEA15MW" / "reference"
        ed = load_elastodyn_blade(ref_dir / "IEA-15-240-RWT_ElastoDyn_blade.dat")
        bd = load_beamdyn_blade(ref_dir / "IEA-15-240-RWT_BeamDyn_blade.dat")

        generator = BladeMesh(
            excel_file=None if blade_name == "official" else str(iea_blade_xlsx),
            yaml_file=iea_blade_yaml if blade_name == "official" else None,
            airfoil_dir=str(iea_blade_xlsx.parent / "airfoils"),
            element_size=0.5,
        )
        mesh = generator.generate(renumber="rcm")
        properties = build_rust_properties(generator.numad_mesh_data)
        beam = BeamReference.from_mesh(mesh, ed, bd)

        cases = [
            ("LC1_gravity", {"type": "gravity", "value": GRAVITY}),
            ("LC2_tip_flap", {"type": "tip", "value": (0.0, TIP_FORCE_N, 0.0)}),
            ("LC3_tip_edge", {"type": "tip", "value": (TIP_FORCE_N, 0.0, 0.0)}),
            ("LC4_distributed", {"type": "distributed",
                                 "value": (0.0, DISTRIB_Q_NM, 0.0), "n_slices": 60}),
        ]
        results = {}
        for name, load in cases:
            u_shell = run_shell_case(mesh, properties, load)
            d_shell = shell_tip_displacement(mesh, u_shell)
            d_beam = beam_tip_displacement(beam, load)
            results[name] = {"shell": d_shell, "beam": d_beam}
        _S2_CACHE[blade_name] = results
        return results

    @pytest.mark.parametrize("case_name", ["LC1_gravity", "LC2_tip_flap", "LC4_distributed"])
    def test_dominant_component_matches_beam(self, s2_results, blade_name, case_name):
        spec = _DOMINANT[case_name]
        comp_idx = {"ux": 0, "uy": 1, "uz": 2}[spec["component"]]
        d_shell = s2_results[case_name]["shell"][comp_idx]
        d_beam = s2_results[case_name]["beam"][comp_idx]
        rel = abs(d_shell - d_beam) / abs(d_beam)
        tol = _BLADES[blade_name]["rel_tol_dominant"]
        assert rel < tol, (
            f"[{blade_name}] {case_name} {spec['component']}: shell={d_shell:.4f} "
            f"beam={d_beam:.4f} (rel={rel:.2%}, tol={tol:.0%})"
        )
        if spec.get("expect_negative"):
            assert d_shell < 0, f"[{blade_name}] {case_name}: expected downward deflection"

    def test_edgewise_ratio_within_band(self, s2_results, blade_name):
        """LC3 edgewise: official blade closes at ~8%; UTD shows +40-55%."""
        lo, hi = _BLADES[blade_name]["lc3_ratio"]
        d_shell = s2_results["LC3_tip_edge"]["shell"][0]
        d_beam = s2_results["LC3_tip_edge"]["beam"][0]
        ratio = d_shell / d_beam
        assert lo < ratio < hi, (
            f"[{blade_name}] LC3 edgewise ratio shell/beam = {ratio:.2f} "
            f"outside [{lo}, {hi}]"
        )
