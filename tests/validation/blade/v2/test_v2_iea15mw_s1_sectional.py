"""S-1 on the windIO v2 blade: sectional EI / mass / GJ / EA vs BeamDyn (v2 counterpart of group 38).

Same reference (official BeamDyn 6x6 deck, elastic-centre principal EI), same metrics (mid-span 0.2-0.85 of the
tip z, mean |rel err| of EI_edge, EI_flap, m(r), GJ, EA(CLT)), same tolerance (5 %, ``_TOL`` as in
``tests/test_iea15mw_s1_sectional.py``), same extractions (tangent-frame curvature EI on two tip-load solves, lumped
mass partition, static torsion GJ with the end-ring wall shear flow, direct CLT EA).  Only the mesh and the station
membership differ: ``build_v2_blade(0.25, skin_reference="mid")`` with the v2 extractor of
``tests/support/v2_sectional.py`` (stations = ring-to-ring bays; evaluated on the 49 bays of the reference-axis grid,
as v1 uses ``native_bay_centres``).  See docs/validation/_pending/63-v2_blade_sectional_distribution/NOTES.md.

Marks: slow (4 static solves on the 28 k-node mesh; run it as a SLURM job, never on a login node).
"""

from __future__ import annotations

import sys

import numpy as np
import pytest

pytest.importorskip("petsc4py", reason="PETSc not available")

from tests.support.assertions import assert_residual_below  # noqa: E402
from tests.support.paths import REPO_ROOT  # noqa: E402
from tests.support.v2_blade import build_v2_blade  # noqa: E402
from tests.support.v2_sectional import (  # noqa: E402
    V2SectionalExtractor,
    bay_centres_v2,
    ea_gj_from_static_response_v2,
    ring_ids_for_v2_mesh,
)

sys.path.insert(0, str(REPO_ROOT / "tools"))

from beam_reference import load_beamdyn_blade  # noqa: E402

ELEMENT_SIZE = 0.25
SKIN_REFERENCE = "mid"
#: Maintainer rule 2026-10-16 (every L3 tolerance <= 5 %), identical to the v1 test.
_TOL = 0.05
_REF = ("BeamDyn blade deck of the official IEA 15 MW reference turbine (IEA-15-240-RWT_BeamDyn_blade.dat, "
        "6x6 section matrices; EI at the elastic centre, principal axes)")


def _interp(x_src, y_src, x_dst):
    return np.interp(x_dst, x_src, y_src)


def compute_s1_v2(element_size=ELEMENT_SIZE, skin_reference=SKIN_REFERENCE):
    """All S-1 channels of the v2 blade (shared by the test fixture and the sbatch driver)."""
    bd = load_beamdyn_blade(REPO_ROOT / "tests/IEA15MW/reference/IEA-15-240-RWT_BeamDyn_blade.dat")
    mesh, props = build_v2_blade(element_size, skin_reference=skin_reference)
    ring = ring_ids_for_v2_mesh(mesh)
    ext = V2SectionalExtractor(mesh, props, ring_ids=ring)
    stations = bay_centres_v2()
    sec = ext.extract(stations=stations)  # mass path
    z_static, ei_flap, ei_edge = ext.ei_from_static_response(props, frame="tangent", stations=stations)
    z_eg, _ea_resp, gj_static, ea_clt = ea_gj_from_static_response_v2(
        mesh, props, numad_data=mesh.numad_mesh_data, stations=stations, ring_ids=ring)
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


def mid_errors(r):
    """(mean |rel err|, signed mean rel err) per channel over 0.2-0.85 of the tip z."""
    out = {}
    for name, (z, a, b) in {
        "edgewise_ei": (r["z"], r["ei_edge"], r["ei_edge_bd"]),
        "flapwise_ei": (r["z"], r["ei_flap"], r["ei_flap_bd"]),
        "mass": (r["z_mass"], r["mass"], r["mass_bd"]),
        "gj": (r["z"], r["gj"], r["gj_bd"]),
        "ea": (r["z"], r["ea"], r["ea_bd"]),
    }.items():
        m = r["mid"] if name in ("edgewise_ei", "flapwise_ei") else (z > 0.2 * z.max()) & (z < 0.85 * z.max())
        rel = a[m] / b[m] - 1
        out[name] = (float(np.nanmean(np.abs(rel))), float(np.nanmean(rel)))
    return out


@pytest.mark.slow
class TestV2S1SectionalProperties:
    @pytest.fixture(scope="class")
    def s1_results(self):
        return compute_s1_v2()

    def test_edgewise_ei_matches_beamdyn(self, s1_results):
        """EI_edge mid-span mean |rel err| vs the BeamDyn elastic-centre EI below 5 %."""
        err = mid_errors(s1_results)["edgewise_ei"][0]
        assert_residual_below(err, tol=_TOL, kind="code", reference_name=_REF,
                              what="S-1 EI_edge mid-span mean |rel err| (windIO v2 mesh)")

    def test_flapwise_ei_within_band(self, s1_results):
        """EI_flap mid-span mean |rel err| vs the BeamDyn elastic-centre EI below 5 %."""
        err = mid_errors(s1_results)["flapwise_ei"][0]
        assert_residual_below(err, tol=_TOL, kind="code", reference_name=_REF,
                              what="S-1 EI_flap mid-span mean |rel err| (windIO v2 mesh)")

    def test_mass_distribution_matches_beamdyn(self, s1_results):
        """Mass per unit length mid-span mean |rel err| vs BeamDyn below 5 %."""
        err = mid_errors(s1_results)["mass"][0]
        assert_residual_below(err, tol=_TOL, kind="code", reference_name=_REF,
                              what="S-1 m(r) mid-span mean |rel err| (windIO v2 mesh)")

    def test_gj_distribution_matches_beamdyn(self, s1_results):
        """GJ(r) mid-span mean |rel err| vs BeamDyn below 5 % (static torsion, end-ring wall shear flow)."""
        err = mid_errors(s1_results)["gj"][0]
        assert_residual_below(err, tol=_TOL, kind="code", reference_name=_REF,
                              what="S-1 GJ(r) mid-span mean |rel err| (windIO v2 mesh)")

    def test_ea_clt_distribution_matches_beamdyn(self, s1_results):
        """EA(r) from the direct CLT section integral, mid-span mean |rel err| vs BeamDyn below 5 %."""
        err = mid_errors(s1_results)["ea"][0]
        assert_residual_below(err, tol=_TOL, kind="code", reference_name=_REF,
                              what="S-1 EA(r) CLT mid-span mean |rel err| (windIO v2 mesh)")
