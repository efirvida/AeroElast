"""S-1 — IEA 15 MW blade: sectional property validation (shell vs beam targets).

The WindIO yaml publishes BOTH structural representations: ``elastic_properties_mb.six_x_six``
(the beam 6x6 matrices, identical to the BeamDyn blade file) and ``internal_structure_2d_fem``
(the 2D layup the shell is built from). The 6x6 is the sectional arbiter (decision (b),
odd/tasks/validation-l3-b3-audit.md); EI is compared with its elastic-centre principal values
(Schur complement of K33, eigenvalues of the bending block), not the reference-axis K44/K55.

Two independent extractions from the assembled shell model:

  1. Mass distribution m(r) via station membership (element sets
     {component}_{station}_{...}) + lumped mass partition.
  2. EI_flap(r) / EI_edge(r) via the classical curvature method on two
     full-mesh static tip-load solves (EI = M(r)/kappa(r), smoothed with a
     spline over the station-mean transverse deflection).

Targets: the official BeamDyn 6x6 section matrices
(tests/IEA15MW/reference/IEA-15-240-RWT_BeamDyn_blade.dat), which are the
PreComp/VABS/BECAS results published for the IEA 15 MW.

The anchor builds the skin on the mid-surface (BladeMesh(skin_reference="mid"), 7324508).
Measured 2026-10-17 after the layup mass-mapping fix (df446c3), element_size 0.25, mid-span
(0.2-0.85R) mean |rel err| vs BeamDyn, signed in brackets (SLURM 11617746 OML / 11617752 mid):

  EI_flap 21.63 % / 5.18 % (+21.63 / +3.33)   EI_edge 11.96 % / 12.70 % (-11.11 / -12.10)
  GJ 10.20 % / 6.13 % (+7.37 / -4.30)         EA(CLT) 8.82 % / 7.59 % (+8.82 / +7.57)
  mass signed +1.79 % / +0.73 %

Anchor move vs 2026-10-16 (SLURM 11617384): the mapping fix removed surplus skin glass, so mass
went green and EI_edge went red (4.7 % -> 12.7 % mid); EI_edge is the S-1 extraction, not the layup
(odd/tasks/beamdyn-arbiter-midsurface.md, "T8 result"). Record: docs/validation_closures.md,
"Superficie media en las anclas de pala".

The flap and GJ excess is the skin sitting on the OML with the laminate stacked outward (E8,
tools/s1_midsurface_counterfactual.py). EA and mass are layup integrals independent of that: an
independent integral of the yaml layup is +9.6 % (EA) and +2.8 % (mass) above the 6x6 itself
(tools/yaml_layup_independent.py), so part of the residual is inside the yaml. Bounds are 5 %
(maintainer rule); the red ones stay red.

EI is read in the section tangent frame (``ei_from_static_response(frame="tangent")``, L28 in
odd/tasks/l3-physical-validation.md): the default global-X/Y reading mixes the prebend-tilted torsion into the
edge moment (1/EI_g = cos^2(phi)/EI_loc + sin^2(phi)/GJ). Measured 2026-10-22, es 0.25, mid, 0.2-0.85R, signed
(mean |rel|) vs the 6x6, global -> tangent (SLURM 11622464): EI_edge -11.45 (11.45) -> -7.51 (8.82) %,
EI_flap +3.93 (5.17) -> +3.79 (5.70) %; GJ -2.76 (5.05) %, EA(CLT) +7.57 (7.59) % are unchanged by the
frame. Bounds unchanged; EI_edge and EI_flap stay red (SLURM 11622473). On the prebend-removed blade the two
frames agree (EI_edge -7.88 global / -8.17 tangent, SLURM 11622465).

Marks: slow (mesh generation + 2 static solves + extraction).
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

petsc4py = pytest.importorskip("petsc4py", reason="PETSc not available")

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "tools"))

from beam_reference import load_beamdyn_blade  # noqa: E402
from run_s1_sectional import ea_gj_from_static_response as _ea_gj_static_response  # noqa: E402


def _interp(x_src, y_src, x_dst):
    return np.interp(x_dst, x_src, y_src)


@pytest.mark.slow
class TestS1SectionalProperties:
    @pytest.fixture(scope="class")
    def s1_results(self, iea_blade_yaml):
        from aeroelast.core.mesh.generators import BladeMesh
        from aeroelast.models.blade.model import build_rust_properties
        from aeroelast.postprocess.sectional import SectionalExtractor

        ref_dir = Path(__file__).resolve().parents[1] / "tests" / "IEA15MW" / "reference"
        bd = load_beamdyn_blade(ref_dir / "IEA-15-240-RWT_BeamDyn_blade.dat")

        generator = BladeMesh(
            yaml_file=str(iea_blade_yaml),
            element_size=0.25,
            skin_reference="mid",  # blade skin on the mid-surface (odd/tasks/beamdyn-arbiter-midsurface.md, T4)
        )
        mesh = generator.generate(renumber="rcm")
        props = build_rust_properties(generator.numad_mesh_data)
        ext = SectionalExtractor(mesh, props)

        sec = ext.extract()  # mass path
        z_static, ei_flap, ei_edge = ext.ei_from_static_response(props, frame="tangent")
        z_eg, ea_resp, gj_static, ea_clt = _ea_gj_static_response(
            mesh, props, numad_data=generator.numad_mesh_data)

        mid = (z_static > 0.2 * z_static.max()) & (z_static < 0.85 * z_static.max())
        return {
            "z": z_static, "mid": mid,
            "ei_flap": ei_flap, "ei_edge": ei_edge,
            "ei_flap_bd": _interp(bd.z, bd.ei_flap_ec, z_static),
            "ei_edge_bd": _interp(bd.z, bd.ei_edge_ec, z_static),
            "z_mass": sec.z, "mass": sec.mass_dens,
            "mass_bd": _interp(bd.z, bd.mass_dens, sec.z),
            "gj": gj_static, "gj_bd": _interp(bd.z, bd.gj, z_eg),
            "ea": ea_clt, "ea_bd": _interp(bd.z, bd.ea, z_eg),
        }

    def test_edgewise_ei_matches_beamdyn(self, s1_results):
        """EI_edge mid-span mean |rel err| vs the BeamDyn elastic-centre EI below 5 %."""
        m = s1_results["mid"]
        rel = np.abs(s1_results["ei_edge"][m] / s1_results["ei_edge_bd"][m] - 1)
        err = float(np.nanmean(rel))
        assert err < 0.05, (
            f"EI_edge mid-span mean |rel err| = {err:.2%} (tol 5%)"
        )

    def test_flapwise_ei_within_band(self, s1_results):
        """EI_flap mid-span mean |rel err| vs the BeamDyn elastic-centre EI below 5 %."""
        m = s1_results["mid"]
        rel = np.abs(s1_results["ei_flap"][m] / s1_results["ei_flap_bd"][m] - 1)
        err = float(np.nanmean(rel))
        assert err < 0.05, (
            f"EI_flap mid-span mean |rel err| = {err:.2%} (tol 5%)"
        )

    def test_mass_distribution_matches_beamdyn(self, s1_results):
        """Mass per unit length mid-span mean |rel err| vs BeamDyn below 5 %."""
        z, m, mb = s1_results["z_mass"], s1_results["mass"], s1_results["mass_bd"]
        mid = (z > 0.2 * z.max()) & (z < 0.85 * z.max())
        rel = np.abs(m[mid] / mb[mid] - 1)
        err = float(np.nanmean(rel))
        assert err < 0.05, (
            f"m(r) mid-span mean |rel err| = {err:.2%} (tol 5%)"
        )

    def test_gj_distribution_matches_beamdyn(self, s1_results):
        """GJ(r) mid-span mean |rel err| vs BeamDyn below 5 %.

        Static-response torsion test: tip torque couple about the local pre-bent span tangent.
        """
        z, gj, gjb = s1_results["z"], s1_results["gj"], s1_results["gj_bd"]
        mid = (z > 0.2 * z.max()) & (z < 0.85 * z.max())
        rel = np.abs(gj[mid] / gjb[mid] - 1)
        err = float(np.nanmean(rel))
        assert err < 0.05, (
            f"GJ(r) mid-span mean |rel err| = {err:.2%} (tol 5%)"
        )

    def test_ea_clt_distribution_matches_beamdyn(self, s1_results):
        """EA(r) from the direct CLT section integral, mid-span mean |rel err| vs BeamDyn below 5 %.

        The CLT path (sum of A11(layup) per element, reciprocal-Poisson reduced stiffness) is
        immune to the pre-bend bending contamination of the global-response extraction.
        """
        z, ea, eab = s1_results["z"], s1_results["ea"], s1_results["ea_bd"]
        mid = (z > 0.2 * z.max()) & (z < 0.85 * z.max())
        rel = np.abs(ea[mid] / eab[mid] - 1)
        err = float(np.nanmean(rel))
        assert err < 0.05, (
            f"EA(r) mid-span mean |rel err| = {err:.2%} (tol 5%)"
        )
