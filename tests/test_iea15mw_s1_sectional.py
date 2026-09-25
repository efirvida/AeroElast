"""S-1 — IEA 15 MW blade: sectional property validation (shell vs beam targets).

Input-chain verification (2026-09-10): the WindIO yaml publishes BOTH
structural representations side by side — ``elastic_properties_mb.six_x_six``
(the beam 6x6 matrices, identical to machine precision with the BeamDyn
blade file) and ``internal_structure_2d_fem`` (the 2D layup the shell is
built from).  The measured shell-vs-beam bands (GJ ~+17%, EI 10-25%,
mass +4-8%) are the disagreement between the yaml's own two
representations, not a shell error — the shell is CLT-consistent with its
layup.  The tests below are regression guards around those documented
bands, NOT claims of exact agreement with the beam tables.

Two independent extractions from the assembled shell model:

  1. Mass distribution m(r) via station membership (element sets
     {component}_{station}_{...}) + lumped mass partition.
  2. EI_flap(r) / EI_edge(r) via the classical curvature method on two
     full-mesh static tip-load solves (EI = M(r)/kappa(r), smoothed with a
     spline over the station-mean transverse deflection).

Targets: the official BeamDyn 6x6 section matrices
(tests/IEA15MW/reference/IEA-15-240-RWT_BeamDyn_blade.dat), which are the
PreComp/VABS/BECAS results published for the IEA 15 MW.

Measured (2026-09-09, element_size 0.5, official WindIO blade):
  mid-span (0.2-0.85R) mean |rel err| vs BeamDyn:
    EI_edge ≈ 10%, EI_flap ≈ 20-25% (systematically high — the global load
    direction mixes flap/edge through the blade twist, contaminating the
    flap extraction upward), m(r) ≈ 1-10%.

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

        generator = BladeMesh(yaml_file=str(iea_blade_yaml), element_size=0.5)
        mesh = generator.generate(renumber="rcm")
        props = build_rust_properties(generator.numad_mesh_data)
        ext = SectionalExtractor(mesh, props)

        sec = ext.extract()  # mass path
        z_static, ei_flap, ei_edge = ext.ei_from_static_response(props)
        z_eg, ea_resp, gj_static, ea_clt = _ea_gj_static_response(
            mesh, props, numad_data=generator.numad_mesh_data)

        mid = (z_static > 0.2 * z_static.max()) & (z_static < 0.85 * z_static.max())

        # Compare by SPAN FRACTION.  Both datums are valid but different: the
        # mesh now lives in the rotor frame (3.97-120.97 m, matching the
        # official ElastoDyn HubRad/TipRad), while the BeamDyn reference file
        # is blade-local.  The fraction is frame-free — and it is how the IEA
        # tabulates the reference anyway.  It must be
        # (x - x.min())/(x.max() - x.min()): with different origins,
        # x/x.max() is NOT the span fraction.
        def _frac(x):
            x = np.asarray(x, dtype=float)
            return (x - x.min()) / (x.max() - x.min())

        zn = _frac(z_static)
        z_mass_n = _frac(sec.z)
        z_eg_n = _frac(z_eg)
        bdn = _frac(bd.z)
        return {
            "z": z_static, "mid": mid,
            "ei_flap": ei_flap, "ei_edge": ei_edge,
            "ei_flap_bd": _interp(bdn, bd.ei_flap, zn),
            "ei_edge_bd": _interp(bdn, bd.ei_edge, zn),
            "z_mass": sec.z, "mass": sec.mass_dens,
            "mass_bd": _interp(bdn, bd.mass_dens, z_mass_n),
            "gj": gj_static, "gj_bd": _interp(bdn, bd.gj, z_eg_n),
            "ea": ea_clt, "ea_bd": _interp(bdn, bd.ea, z_eg_n),
        }

    def test_edgewise_ei_matches_beamdyn(self, s1_results):
        """EI_edge mid-span within 25% of the BeamDyn targets."""
        m = s1_results["mid"]
        rel = np.abs(s1_results["ei_edge"][m] / s1_results["ei_edge_bd"][m] - 1)
        err = float(np.nanmean(rel))
        assert err < 0.25, (
            f"EI_edge mid-span mean |rel err| = {err:.2%} (tol 25%)"
        )

    def test_flapwise_ei_within_band(self, s1_results):
        """EI_flap mid-span within 40% (global-load direction mixing biases
        the flap extraction upward; see module docstring)."""
        m = s1_results["mid"]
        rel = np.abs(s1_results["ei_flap"][m] / s1_results["ei_flap_bd"][m] - 1)
        err = float(np.nanmean(rel))
        assert err < 0.40, (
            f"EI_flap mid-span mean |rel err| = {err:.2%} (tol 40%)"
        )

    def test_mass_distribution_matches_beamdyn(self, s1_results):
        """Mass per unit length mid-span within 15% of the BeamDyn targets."""
        z, m, mb = s1_results["z_mass"], s1_results["mass"], s1_results["mass_bd"]
        mid = (z > 0.2 * z.max()) & (z < 0.85 * z.max())
        rel = np.abs(m[mid] / mb[mid] - 1)
        err = float(np.nanmean(rel))
        assert err < 0.15, (
            f"m(r) mid-span mean |rel err| = {err:.2%} (tol 15%)"
        )

    def test_gj_distribution_matches_beamdyn(self, s1_results):
        """GJ(r) mid-span within 30% of the BeamDyn targets.

        Static-response torsion test (tip torque couple -> section twist
        rate).  Post twist-sign fix (2026-09-15), the couple is applied
        about the LOCAL pre-bent span tangent and theta is projected on it:
        a global-z couple on the swept blade mixes bending into the twist
        (same artefact as EI_flap/EA).  Measured 23.9% mean |rel| — the
        yaml's two structural representations (2D layup vs the 6x6),
        consistent with the CLT/Bredt layup judge (~20%); see
        docs/shell_vs_beam_sectional_validation.md section 6.
        """
        z, gj, gjb = s1_results["z"], s1_results["gj"], s1_results["gj_bd"]
        mid = (z > 0.2 * z.max()) & (z < 0.85 * z.max())
        rel = np.abs(gj[mid] / gjb[mid] - 1)
        err = float(np.nanmean(rel))
        assert err < 0.30, (
            f"GJ(r) mid-span mean |rel err| = {err:.2%} (tol 30%)"
        )

    def test_ea_clt_distribution_matches_beamdyn(self, s1_results):
        """EA(r) from the direct CLT section integral within 30% of BeamDyn.

        The CLT path (sum of A11(layup) per element) is immune to the
        pre-bend bending contamination that biased the global-response
        extraction.  Measured ~20% (2026-09-10) — same representation band
        as the GJ.
        """
        z, ea, eab = s1_results["z"], s1_results["ea"], s1_results["ea_bd"]
        mid = (z > 0.2 * z.max()) & (z < 0.85 * z.max())
        rel = np.abs(ea[mid] / eab[mid] - 1)
        err = float(np.nanmean(rel))
        assert err < 0.30, (
            f"EA(r) mid-span mean |rel err| = {err:.2%} (tol 30%)"
        )
