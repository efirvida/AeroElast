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

Arbiter for the official blade (odd/tasks/validation-l3-b3-audit.md K9, 2026-10-21): OpenFAST v5.0.0
BeamDyn driver, single member p60, linear static, official v1.1 deck (T3, d19ef69, report
logs/check/bd_static/report.md, p40-p60 within 0.13 %). Tip uy [m]: LC1 -2.231483, LC2 82.273278,
LC4 9.062631. The in-house EB beam above is printed as a comparison only: its geometry is sliced from
the shell mesh, so it moves with the mesh it judges (es 0.125: -13.3 % on LC1). Shell vs BeamDyn
(mid-surface skin): es 0.25 / 0.125 LC1 -7.41 / -7.20 %, LC2 -3.68 / -4.05 %, LC4 -6.22 / -6.24 %.
LC2 is printed, not asserted: node and ring tip uy differ by 9-11 % from local load introduction
(T4 report), so it is not a clean metric. The K9 choice was made after these numbers were known.

LC3 (K11.3, D7): the shell tip ux is judged against the same BeamDyn driver, LC3 case (shell tip force
(+1 MN, 0, 0) = BD (0, -1 MN, 0)), p60 linear tip ux 34.469441 m (p50-p60 0.047 %, scale 1e-4 vs 2e-4
0.0005 %; logs/check/bd_static/LC3-p*). The 0.9-1.2 band on the shell/EB ratio is gone; that ratio
(1.25, known before K11) is printed only. The tip-ring mean ux is printed next to the node for the
K11.3 clean-metric gate (node vs ring within 0.5 %, as K9.3). Measured at es 0.25 (SLURM 11621347):
shell 40.310002 m, +16.94 % vs BeamDyn, ring mean 39.837109 m, node vs ring +1.19 % > 0.5 %, so by
K11.3 LC3 is printed, not asserted (`undecided`), like LC2.

Marks: slow (blade mesh generation + 4 static solves per blade).
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

petsc4py = pytest.importorskip("petsc4py", reason="PETSc not available")

REPO_ROOT = Path('/scratch/leahk/eduardo.donestevez/fem-shell')
sys.path.insert(0, str(REPO_ROOT / "tools"))

from beam_reference import load_beamdyn_blade, load_elastodyn_blade  # noqa: E402

GRAVITY = (0.0, -9.81, 0.0)
TIP_FORCE_N = 1.0e6
DISTRIB_Q_NM = 5.0e3

# Per-source expectations (measured 2026-09-08 at element_size 0.5 (stale: unconverged mesh, re-measure at 0.25)).
# The UTD LC3 band stays an input-sensitivity record against the EB beam; the official blade has none (K11.3).
_BLADES = {
    "official": {"rel_tol_dominant": 0.05},
    "utd": {"rel_tol_dominant": 0.05, "lc3_ratio": (1.0, 2.2)},
}

_DOMINANT = {
    "LC1_gravity": {"component": "uy", "expect_negative": True},
    "LC2_tip_flap": {"component": "uy"},
    "LC4_distributed": {"component": "uy"},
}

# OpenFAST v5.0.0 BeamDyn driver, p60, linear static, official deck: tip uy [m] (K9).
_BD_STATIC_UY = {"LC1_gravity": -2.231483, "LC2_tip_flap": 82.273278, "LC4_distributed": 9.062631}
_NOT_ASSERTED_VS_BD = {"LC2_tip_flap"}
# Same driver, LC3 case (K11.3): linear tip ux [m] in shell axes, -TipTDyr / scale.
_BD_STATIC_LC3_UX = 34.469441
# K11.3 clean-metric gate failed at es 0.25 (node vs ring +1.19 % > 0.5 %): printed, not asserted.
_LC3_ASSERTED_VS_BD = False

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

        ref_dir = Path('/scratch/leahk/eduardo.donestevez/fem-shell') / "tests" / "IEA15MW" / "reference"
        ed = load_elastodyn_blade(ref_dir / "IEA-15-240-RWT_ElastoDyn_blade.dat")
        bd = load_beamdyn_blade(ref_dir / "IEA-15-240-RWT_BeamDyn_blade.dat")

        generator = BladeMesh(
            excel_file=None if blade_name == "official" else str(iea_blade_xlsx),
            yaml_file=iea_blade_yaml if blade_name == "official" else None,
            airfoil_dir=str(iea_blade_xlsx.parent / "airfoils"),
            element_size=0.125,
            skin_reference="mid",  # blade skin on the mid-surface (odd/tasks/beamdyn-arbiter-midsurface.md, T4)
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
            ring = np.nonzero(mesh.coords_array[:, 2] >= mesh.coords_array[:, 2].max() - 1e-6)[0]
            u_ring = np.asarray(u_shell).reshape(-1, 6)[ring, :3].mean(axis=0)
            results[name] = {"shell": d_shell, "beam": d_beam, "ring": u_ring}
        _S2_CACHE[blade_name] = results
        return results

    @pytest.mark.parametrize("case_name", ["LC1_gravity", "LC2_tip_flap", "LC4_distributed"])
    def test_dominant_component_matches_beam(self, s2_results, blade_name, case_name):
        spec = _DOMINANT[case_name]
        comp_idx = {"ux": 0, "uy": 1, "uz": 2}[spec["component"]]
        d_shell = s2_results[case_name]["shell"][comp_idx]
        d_beam = s2_results[case_name]["beam"][comp_idx]
        tol = _BLADES[blade_name]["rel_tol_dominant"]
        if blade_name == "official":
            d_bd = _BD_STATIC_UY[case_name]
            print(
                f"S-2 {case_name} uy [m]: shell {d_shell:.6f} vs OpenFAST BeamDyn p60 {d_bd:.6f} "
                f"(signed error: {d_shell / d_bd - 1:+.4%}, bound: {tol:.4%}); in-house EB beam "
                f"{d_beam:.6f} ({d_shell / d_beam - 1:+.4%}, comparison, not asserted)"
            )
            if case_name in _NOT_ASSERTED_VS_BD:
                return
            d_ref, ref_name = d_bd, "BeamDyn p60"
        else:
            d_ref, ref_name = d_beam, "beam"
        rel = abs(d_shell - d_ref) / abs(d_ref)
        assert rel < tol, (
            f"[{blade_name}] {case_name} {spec['component']}: shell={d_shell:.4f} "
            f"{ref_name}={d_ref:.4f} (rel={rel:.2%}, tol={tol:.0%})"
        )
        if spec.get("expect_negative"):
            assert d_shell < 0, f"[{blade_name}] {case_name}: expected downward deflection"

    def test_edgewise_ratio_within_band(self, s2_results, blade_name):
        """LC3 edgewise tip ux: official vs BeamDyn p60 (K11.3); UTD vs the EB beam band (+40-55 %)."""
        d_shell = s2_results["LC3_tip_edge"]["shell"][0]
        d_beam = s2_results["LC3_tip_edge"]["beam"][0]
        ratio = d_shell / d_beam
        if blade_name == "official":
            d_ring = s2_results["LC3_tip_edge"]["ring"][0]
            tol = _BLADES[blade_name]["rel_tol_dominant"]
            d_bd = _BD_STATIC_LC3_UX
            print(
                f"S-2 LC3_tip_edge ux [m]: shell {d_shell:.6f} vs OpenFAST BeamDyn p60 {d_bd:.6f} "
                f"(signed error: {d_shell / d_bd - 1:+.4%}, bound: {tol:.4%}); tip ring mean {d_ring:.6f} "
                f"(node vs ring {d_shell / d_ring - 1:+.4%}); in-house EB beam {d_beam:.6f} "
                f"(ratio {ratio:.4f}, comparison, not asserted)"
            )
            if not _LC3_ASSERTED_VS_BD:
                return
            rel = abs(d_shell - d_bd) / abs(d_bd)
            assert rel < tol, (
                f"[official] LC3_tip_edge ux: shell={d_shell:.4f} BeamDyn p60={d_bd:.4f} "
                f"(rel={rel:.2%}, tol={tol:.0%})"
            )
            return
        lo, hi = _BLADES[blade_name]["lc3_ratio"]
        assert lo < ratio < hi, (
            f"[{blade_name}] LC3 edgewise ratio shell/beam = {ratio:.2f} "
            f"outside [{lo}, {hi}]"
        )
