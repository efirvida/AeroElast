"""Campbell cross-validation: shell vs OpenFAST MBC3 at the sweep speeds.

The MBC3 collective members live between the backward/forward whirl pairs;
the single-blade shell's frequencies must fall close to the OF frequency
NEAREST to them in each family band (flap / edge).  The comparison
tolerates 10%: the shell's edgewise sits ~-5% below the MBC3 and the flap
shows the documented K_G over-response that grows with Omega (see
docs/shell_vs_beam_sectional_validation.md and the S-4 notes).

Consumes:
  * docs/validation_data/generated/sx_campbell.csv (the shell sweep)
  * docs/validation_data/generated/sx_campbell_openfast.csv (the MBC3 sweep)
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
SHELL_CSV = REPO_ROOT / "docs" / "validation_data" / "generated" / "sx_campbell.csv"
OF_CSV = REPO_ROOT / "docs" / "validation_data" / "generated" / "sx_campbell_openfast.csv"

REL_TOL = 0.10


def _shell_rows():
    if not SHELL_CSV.exists():
        pytest.skip(f"missing {SHELL_CSV}")
    d = np.loadtxt(SHELL_CSV, delimiter=",", skiprows=1)
    if d.ndim == 1:
        d = d[None, :]
    return {"rpm": d[:, 0], "f1f": d[:, 1], "f2f": d[:, 2], "f1e": d[:, 3]}


def _of_rows():
    if not OF_CSV.exists():
        pytest.skip(f"missing {OF_CSV} — run tools/run_sx_of_campbell.py")
    d = np.loadtxt(OF_CSV, delimiter=",", skiprows=1)
    if d.ndim == 1:
        d = d[None, :]
    return {"rpm": d[:, 0], "freq": d[:, 1]}


@pytest.mark.slow
class TestCampbellCrossValidation:
    @pytest.fixture(scope="class")
    def data(self):
        return {"shell": _shell_rows(), "of": _of_rows()}

    def _match(self, rpm, of_rpm_arr, of_freq, target):
        sel = np.abs(of_rpm_arr - rpm) < 1e-9
        cands = of_freq[sel]
        near = cands[np.argmin(np.abs(cands - target))]
        return near

    def test_flap_matches_mbc3_across_sweep(self, data):
        sh, of = data["shell"], data["of"]
        for rpm in np.unique(of["rpm"]):
            f_shell = np.interp(rpm, sh["rpm"], sh["f1f"])
            near = self._match(rpm, of["rpm"], of["freq"], f_shell)
            rel = abs(f_shell - near) / near
            assert rel < REL_TOL, (
                f"rpm={rpm:g}: shell 1st flap {f_shell:.4f} Hz vs MBC3 nearest "
                f"{near:.4f} Hz (rel {rel:.2%} > {REL_TOL:.0%})"
            )

    def test_edge_matches_mbc3_across_sweep(self, data):
        sh, of = data["shell"], data["of"]
        for rpm in np.unique(of["rpm"]):
            f_shell = np.interp(rpm, sh["rpm"], sh["f1e"])
            near = self._match(rpm, of["rpm"], of["freq"], f_shell)
            rel = abs(f_shell - near) / near
            assert rel < REL_TOL, (
                f"rpm={rpm:g}: shell 1st edge {f_shell:.4f} Hz vs MBC3 nearest "
                f"{near:.4f} Hz (rel {rel:.2%} > {REL_TOL:.0%})"
            )

    def test_sweeps_share_speed_points(self, data):
        """Both sweeps must cover the same speed set."""
        sh_rpm = np.round(data["shell"]["rpm"], 4)
        of_rpm = np.round(data["of"]["rpm"], 4)
        for rpm in sorted(set(of_rpm)):
            assert rpm in sh_rpm, f"shell sweep missing rpm={rpm:g}"
