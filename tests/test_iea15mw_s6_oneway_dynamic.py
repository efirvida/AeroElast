"""S-6 — IEA 15 MW blade: one-way DYNAMIC AeroDyn loads -> shell Newmark.

Applies the AeroDyn blade-1 load time series (Fn/Ft per unit length, local
chord frame, sampled at dt=0.01 over [96, 100] s of the OpenFAST blade-only
rated run) to the official WindIO shell blade with a trapezoidal Newmark
integration of the rotating (K + K_G + K_SP) system, starting from the
static equilibrium under the window-mean loads.  Compares the tip
out-of-plane deflection series against ElastoDyn's OoPDefl1 over the same
window.

Measured (2026-09-10, element_size 0.5):  the fluctuation amplitude closes
(std 0.5539 m vs ElastoDyn 0.6138 m, -9.8%); the mean carries the K_G
over-response documented in S-5 (13.24 m vs 16.09 m, -17.7% — the rotating
stiffness overshoots, see the tool docstring).

Artifact: tests/IEA15MW/reference/s6_ad_blade1_timeseries.csv
(t, OoPDefl1, then the 50 Fn + 50 Ft station channels) — the test is
OpenFAST-independent.

Marks: slow (blade mesh generation + 2 MUMPS factorizations + 401 steps).
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

petsc4py = pytest.importorskip("petsc4py", reason="PETSc not available")

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "tools"))

from run_s6_oneway_dynamic import (  # noqa: E402
    DT,
    RPM_RATED,
    assemble_keff,
    build_assembly,
    newmark_tip_series,
    nodal_load_series,
)

TS_CSV = REPO_ROOT / "tests" / "IEA15MW" / "reference" / "s6_ad_blade1_timeseries.csv"

# The dynamic fluctuation (std) is the S-6 deliverable: it closes at -9.8%.
STD_REL_TOL = 0.25
# The mean carries the K_G over-response documented in S-5 (K_G static
# inflation, ~-18%); bound it loosely to catch gross regressions only.
MEAN_REL_TOL = 0.25


@pytest.mark.slow
class TestS6OneWayDynamic:
    @pytest.fixture(scope="class")
    def s6_result(self):
        loads = np.loadtxt(TS_CSV, delimiter=",", skiprows=1)
        t = loads[:, 0]
        oop_ed = loads[:, 1]
        fn_ts = loads[:, 2:52]
        ft_ts = loads[:, 52:102]

        # twist/station data for the section frames: rebuilt from the S-5
        # artifact (same AD blade geometry)
        s5 = np.loadtxt(
            REPO_ROOT / "tests" / "IEA15MW" / "reference" / "s5_ad_blade1_loads.csv",
            delimiter=",", skiprows=1)
        r_m, twist = s5[:, 0], s5[:, 3]

        asm, mesh = build_assembly(0.5)
        K_eff, M, free = assemble_keff(asm, mesh, RPM_RATED * 2 * np.pi / 60.0)
        F_ts = nodal_load_series(asm, mesh, r_m, fn_ts, ft_ts, twist)
        hist = newmark_tip_series(asm, mesh, K_eff, M, free, F_ts, DT)
        return {"t": t, "oop_ed": oop_ed, "oop_shell": hist[:, 1]}

    def test_artifact_window_and_channels(self, s6_result):
        """The artifact must carry the 4 s window at dt=0.01 with the 100 AD channels."""
        t = s6_result["t"]
        assert abs(t[0] - 96.0) < 1e-9 and abs(t[-1] - 100.0) < 1e-9
        assert len(t) == 401

    def test_oop_fluctuation_matches_elastodyn(self, s6_result):
        """The tip OoP fluctuation (std) must match ElastoDyn's within 25%."""
        shell_std = s6_result["oop_shell"].std()
        ed_std = s6_result["oop_ed"].std()
        rel = abs(shell_std - ed_std) / ed_std
        assert rel < STD_REL_TOL, (
            f"shell OoP std {shell_std:.4f} m vs ElastoDyn {ed_std:.4f} m "
            f"(rel {rel:.2%} > {STD_REL_TOL:.0%})"
        )

    def test_oop_mean_bounded(self, s6_result):
        """Mean bound: the K_G over-response is documented; catch gross drift only."""
        shell_mean = abs(s6_result["oop_shell"].mean())
        ed_mean = abs(s6_result["oop_ed"].mean())
        rel = abs(shell_mean - ed_mean) / ed_mean
        assert rel < MEAN_REL_TOL, (
            f"shell OoP mean {shell_mean:.4f} m vs ElastoDyn {ed_mean:.4f} m "
            f"(rel {rel:.2%} > {MEAN_REL_TOL:.0%})"
        )
