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

Interior station (l3-closure T10, S12; `TestS2InteriorStation`, official blade only). The tip metric moved
three times and each move is a finding, so the LC1-LC4 claim now lives away from the free tip:

  1. tip node (`argmax z`, the first tip-ring node in mesh order; the K9/K11 metric above): an arbitrary material
     point of a ring whose load is split equally, 9-11 % from the ring mean on LC2/LC3 (T4 report, T10b).
  2. tip-ring mean (L12(4)): node-count dependent where the ring spacing is non-uniform, so the load
     introduction was made defined (tip `arc` = uniform line traction, LC4 `area` = tributary-area traction;
     T10c). On the closed-form tube that did not converge either: an edge load not tangent to the wall has a
     transverse-shear edge layer of order t that es 0.25 / 0.125 do not resolve (rect flap tip ring +1.24 ->
     +1.54 % vs Timoshenko, T10c RESULT), while rings >= 2 section depths inboard converge to the closed form.
  3. interior station (S12): BeamDyn p60 FE-node translation (Gauss quadrature, so output nodes are the GLL FE
     nodes) vs the shell arc-weighted ring mean at the node's z, linearly interpolated between the bracketing
     rings. Predeclared node rule: inner of the two nodes nearest the tip that are >= 2 tip chords inboard; z is
     read from the BeamDyn reference output, not the GLL estimate. Dominant component uy (LC1/LC2/LC4) or ux
     (LC3: BeamDyn -N#TDyr). BeamDyn is normalised by its load scale (1e-4 and 2e-4 must agree), the shell
     response is linear in the load. Tolerance K3, 5 %. The tip-node tests above stay as they were (printed or
     asserted per K9/K11) and the `utd` param is untouched.

  4. what is asserted (S13, after T10c-station refuted the predeclared convergence gate P2): the distributed
     loads LC1 (gravity) and LC4 pass every gate on both blades (es 0.25 -> 0.125 <= 0.32 pp, stations <= 0.21 pp,
     arc vs area <= 0.18 pp) and are asserted at K3. The tip point loads LC2 / LC3 move 0.5-1 pp with refinement
     at every station and stay station dependent 2-4 chords inboard: the shell has no rigid tip coupling (MPC),
     so BeamDyn's axis point load is not reproduced, and they are printed as diagnostics, not asserted.

Marks: slow (blade mesh generation + 4 static solves per blade, +4 for the station loads).
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

petsc4py = pytest.importorskip("petsc4py", reason="PETSc not available")

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "tools"))

sys.path.insert(0, str(REPO_ROOT / "tests"))

from beam_reference import load_beamdyn_blade, load_elastodyn_blade  # noqa: E402
from support.assertions import assert_residual_below  # noqa: E402
from support.s2_station import (  # noqa: E402
    bd_station_nodes,
    component_index,
    pick_station_node,
    ring_mean_at_z,
)

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

# Interior station (S12): BeamDyn p60 reference nodes (tools/bd_static_reference.py post --nodes 53,...,61).
_BD_STATIC_DIR = REPO_ROOT / "tests" / "IEA15MW" / "bd_static"
_STATION_TOL = 0.05  # K3
_BD_SCALE_SPREAD = 1.0e-3  # K9 linearity gate: scales 1e-4 / 2e-4 agree within 0.1 %
# Load introduction as predeclared: tip `arc`, LC4 `area` (defined tractions, independent of the ring node count).
# S13: only the distributed loads are a validation claim; LC2 / LC3 are printed (see the module docstring, item 4).
_STATION_ASSERTED = ("LC1", "LC4")
_STATION_CASES = {
    "LC1": ("LC1_gravity", {"type": "gravity", "value": GRAVITY}),
    "LC2": ("LC2_tip_flap", {"type": "tip", "value": (0.0, TIP_FORCE_N, 0.0), "distribution": "arc"}),
    "LC3": ("LC3_tip_edge", {"type": "tip", "value": (TIP_FORCE_N, 0.0, 0.0), "distribution": "arc"}),
    "LC4": ("LC4_distributed", {"type": "distributed", "value": (0.0, DISTRIB_Q_NM, 0.0),
                                "n_slices": 60, "distribution": "area"}),
}

# Cache of expensive per-blade results shared across the test session.
_S2_CACHE: dict = {}


def _blade_mesh(blade_name, iea_blade_yaml, iea_blade_xlsx):
    """(mesh, properties) of a blade at element_size 0.25, mid-surface skin; cached per blade."""
    key = ("mesh", blade_name)
    if key not in _S2_CACHE:
        from aeroelast.core.mesh.generators import BladeMesh
        from aeroelast.models.blade.model import build_rust_properties

        generator = BladeMesh(
            excel_file=None if blade_name == "official" else str(iea_blade_xlsx),
            yaml_file=iea_blade_yaml if blade_name == "official" else None,
            airfoil_dir=str(iea_blade_xlsx.parent / "airfoils"),
            element_size=0.25,
            skin_reference="mid",  # blade skin on the mid-surface (odd/tasks/beamdyn-arbiter-midsurface.md, T4)
        )
        mesh = generator.generate(renumber="rcm")
        _S2_CACHE[key] = (mesh, build_rust_properties(generator.numad_mesh_data))
    return _S2_CACHE[key]


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
        ref_dir = Path(__file__).resolve().parents[1] / "tests" / "IEA15MW" / "reference"
        ed = load_elastodyn_blade(ref_dir / "IEA-15-240-RWT_ElastoDyn_blade.dat")
        bd = load_beamdyn_blade(ref_dir / "IEA-15-240-RWT_BeamDyn_blade.dat")

        mesh, properties = _blade_mesh(blade_name, iea_blade_yaml, iea_blade_xlsx)
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


@pytest.mark.slow
class TestS2InteriorStation:
    """S12: BeamDyn p60 node displacement vs the shell arc ring mean at the node's z (official blade)."""

    @pytest.fixture
    def station_profiles(self, iea_blade_yaml, iea_blade_xlsx):
        """Per case: (ring profile near the station, tip z, tip chord); cached, loads as predeclared."""
        if "station" in _S2_CACHE:
            return _S2_CACHE["station"]

        from run_s2_static_cases import run_shell_case, tip_ring_scale  # noqa: E402
        from s2_ring_metric_blade import profile_means, ring_profile  # noqa: E402

        mesh, properties = _blade_mesh("official", iea_blade_yaml, iea_blade_xlsx)
        z_tip = float(mesh.coords_array[:, 2].max())
        chord = tip_ring_scale(mesh)
        # the deepest station any case can pick is bounded by the BeamDyn nodes; keep the rings of the last 8 m
        rings = ring_profile(mesh, z_tip - 8.0)
        out = {}
        for case, (_name, load) in _STATION_CASES.items():
            u = np.asarray(run_shell_case(mesh, properties, load))
            out[case] = profile_means(rings, u)
        _S2_CACHE["station"] = (out, z_tip, chord)
        return _S2_CACHE["station"]

    @pytest.mark.parametrize("case", list(_STATION_CASES))
    def test_station_matches_beamdyn_node(self, station_profiles, case):
        profiles, z_tip, chord = station_profiles
        nodes, spread = bd_station_nodes(_BD_STATIC_DIR, case)
        assert spread < _BD_SCALE_SPREAD, f"{case}: BeamDyn scales 1e-4 / 2e-4 differ by {spread:.3e} (gate {_BD_SCALE_SPREAD})"
        node = pick_station_node(nodes, z_tip=z_tip, chord=chord)
        k = component_index(case)
        comp = "xy"[k]
        d_bd = node[f"linear_shell_u{comp}_m"]
        d_shell = ring_mean_at_z(profiles[case], node["z"])[k]
        rel = abs(d_shell - d_bd) / abs(d_bd)
        if case not in _STATION_ASSERTED:
            print(f"S-2 {case} station u{comp} [m] diagnostic (S13, not asserted: tip point load without MPC): "
                  f"shell {d_shell:.6f} vs BeamDyn p60 node {node['node']} {d_bd:.6f}, signed {d_shell / d_bd - 1:+.4%}")
            return
        assert_residual_below(
            rel,
            tol=_STATION_TOL,
            reference_name=f"OpenFAST v5.0.0 BeamDyn p60 node {node['node']} (z {node['z']:.3f} m, {(z_tip - node['z']) / chord:.2f} chords)",
            kind="code",
            what=f"S-2 {case} station u{comp} [m] shell arc ring mean {d_shell:.6f} (BeamDyn {d_bd:.6f}, signed {d_shell / d_bd - 1:+.4%})",
        )
