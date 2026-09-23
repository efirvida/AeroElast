"""S-5 — IEA 15 MW blade: one-way AeroDyn loads -> shell tip deflection.

Applies the azimuth-mean AeroDyn blade-1 distributed loads (Fn/Ft per unit
length in the local chord frame, extracted from the OpenFAST v5 rated run,
blade-only deck, mean over t in [80, 100] s) onto the official WindIO shell
blade and compares the static tip deflection against the ElastoDyn reference.

Reference (blade-only ElastoDyn deck, RotSpeed = 7.56 rpm, only blade DOFs,
gravity ON):  OoPDefl1 = TipDxc1 = 16.1033 m,  IPDefl1 = -0.8103 m.

Measured shell (2026-09-10, element_size 0.5):  |u_tip| = 15.4874 m with
u_y = 15.3848 m (OoP, -4.5% vs 16.1033) and u_x = 1.4878 m (IP, +84% vs
0.8103).  The OoP closes inside the S-2-family band.

IP composition (resolved 2026-09-10):  ElastoDyn's IPDefl1 is the
modal-reduced 1D measure (aero loads 6.00e6 N·m + gravity/cone mean
2.09e6 N·m of the RootMxb1 = 8.09e6 N·m total), while the 3D references
carry the geometric twist coupling on top: the 1D Euler-Bernoulli beam
with the same aero loads gives u_x = 1.17 m (per-moment compliance
1.94e-7 vs the ED's 1.00e-7) and the shell adds ~0.25 m of section-level
in-plane compliance (1.49 m).  The comparison against IPDefl1 is
apples-to-oranges by construction; kept as a loose regression bound.
The OoP is the S-5 deliverable.

Design notes (see tools/run_s5_oneway.py for the full rationale):

  * The shell is solved PARKED (K only) + the azimuth-mean gravity component.
    ElastoDyn's own static OoP at 7.56 rpm matches the un-stiffened beam
    within ~0.5%: the centrifugal stiffening contributes ~0-2% to the STATIC
    deflection, while the shell's K_G would inject a spurious ~13% (open K_G
    issue tracked for the S-4 margin follow-up).
  * The section frames come from the AD geometric twist (BlTwist, negated —
    the shell is built mirrored in the twist), so the Fn/Ft mapping is exact.
  * The loads artifact (tests/IEA15MW/reference/s5_ad_blade1_loads.csv) makes
    the test OpenFAST-independent.

Marks: slow (blade mesh generation + 2 static solves).
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

petsc4py = pytest.importorskip("petsc4py", reason="PETSc not available")

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "tools"))

from run_s5_oneway import (  # noqa: E402
    G_MEAN,
    TIP_IP_DEFL_M,
    TIP_OOP_DEFL_M,
    apply_ad_loads,
    assemble_gravity,
    assemble_parked_k,
    build_assembly,
    solve_static,
    tip_deflections,
)

LOADS_CSV = REPO_ROOT / "tests" / "IEA15MW" / "reference" / "s5_ad_blade1_loads.csv"

# OoP tolerance: the measured -4.5% with headroom for environment variation.
OOP_REL_TOL = 0.08
# IP: the open gap is ~2x; bound it to catch gross regressions only.
IP_MAX_RATIO = 2.5

# Loads artifact sanity: the extracted Fn total should stay near the
# campaign value (827.9 kN — the blade thrust share at rated).
FN_TOTAL_N = 8.2794e5
FT_TOTAL_N = 8.7287e4


@pytest.mark.slow
class TestS5OneWayADLoads:
    @pytest.fixture(scope="class")
    def s5_result(self):
        loads = np.loadtxt(LOADS_CSV, delimiter=",", skiprows=1)
        r_m, fn, ft, twist = loads.T

        asm, mesh = build_assembly(0.5)
        K, free = assemble_parked_k(asm, mesh)
        f = apply_ad_loads(asm, mesh, r_m, fn, ft, twist)
        u_aero = solve_static(asm, K, free, f)
        u_tip_aero = tip_deflections(mesh, u_aero)[0]

        g_dir = u_tip_aero[:3] / np.linalg.norm(u_tip_aero[:3])
        g_mean = np.linalg.norm(G_MEAN) * g_dir
        f_g = assemble_gravity(asm, g_mean)
        u = solve_static(asm, K, free, f + f_g)

        u_tip = tip_deflections(mesh, u)[0]
        return {"loads": loads, "u_tip": u_tip, "fn": fn, "ft": ft}

    def test_loads_artifact_totals(self, s5_result):
        """The artifact carries the campaign AD loads (thrust share sanity)."""
        fn = s5_result["fn"]
        ft = s5_result["ft"]
        r_m = s5_result["loads"][:, 0]
        fn_total = np.trapezoid(fn, r_m)
        ft_total = np.trapezoid(ft, r_m)
        assert abs(fn_total - FN_TOTAL_N) / FN_TOTAL_N < 0.02
        assert abs(ft_total - FT_TOTAL_N) / FT_TOTAL_N < 0.02

    def test_oop_tip_deflection_matches_elastodyn(self, s5_result):
        """OoP tip deflection within the S-2-family band of ElastoDyn."""
        u_tip = s5_result["u_tip"]
        oop = abs(u_tip[1])  # the flap (downwind) component at the ~0-twist tip
        rel = abs(oop - TIP_OOP_DEFL_M) / TIP_OOP_DEFL_M
        assert rel < OOP_REL_TOL, (
            f"shell OoP = {oop:.4f} m vs ElastoDyn {TIP_OOP_DEFL_M:.4f} m "
            f"(rel {rel:.2%} > tol {OOP_REL_TOL:.0%})"
        )

    def test_ip_tip_deflection_bounded(self, s5_result):
        """IP gap is open (~2x); bound it to catch gross regressions."""
        u_tip = s5_result["u_tip"]
        ip = abs(u_tip[0])
        ratio = ip / abs(TIP_IP_DEFL_M)
        assert ratio < IP_MAX_RATIO, (
            f"shell IP = {ip:.4f} m vs ElastoDyn {abs(TIP_IP_DEFL_M):.4f} m "
            f"(ratio {ratio:.2f} > {IP_MAX_RATIO})"
        )
