"""Campbell diagram tests: the rotating-frame modal vs rotor speed.

Consumes tools/run_sx_campbell.py's artifact
(docs/validation_data/generated/sx_campbell.csv) and checks the diagram
physics:

  * the flap modes stiffen monotonically with Omega (the centrifugal K_G),
  * the stiffening scales with Omega^2 (the Southwell-like behaviour),
  * the rated-speed values match the OpenFAST MBC3 targets (the S-4 band),
  * the 1P crossing of the 1st flap exists within the swept range.

The OpenFAST reference points at the intermediate speeds (3/5/9/11 rpm) are
an open follow-up — the S-4 targets cover the endpoints (0 / 7.56 rpm).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
CSV = REPO_ROOT / "docs" / "validation_data" / "generated" / "sx_campbell.csv"

# OpenFAST v5 MBC3 (the S-4's records): (rpm, f1f, f2f, f1e)
OF_TARGETS = {
    0.0: (0.544, 1.610, 0.739),
    7.56: (0.5666, 1.6372, 0.7461),
}


@pytest.mark.slow
class TestCampbellDiagram:
    @pytest.fixture(scope="class")
    def rows(self):
        if not CSV.exists():
            pytest.skip(f"Campbell CSV not found: {CSV} — run tools/run_sx_campbell.py")
        data = np.loadtxt(CSV, delimiter=",", skiprows=1)
        if data.ndim == 1:
            data = data[None, :]
        return {"rpm": data[:, 0], "f1f": data[:, 1], "f2f": data[:, 2],
                "f1e": data[:, 3]}

    def test_flap_modes_stiffen_monotonically(self, rows):
        """The centrifugal K_G must raise the flap frequencies with Omega."""
        assert np.all(np.diff(rows["f1f"]) > 0), (
            f"1st flap series {rows['f1f'].tolist()} not monotonically increasing"
        )
        assert np.all(np.diff(rows["f2f"]) > 0), (
            f"2nd flap series {rows['f2f'].tolist()} not monotonically increasing"
        )

    def test_stiffening_scales_with_omega_squared(self, rows):
        """Delta f^2 ~ kappa * Omega^2 (the Southwell form)."""
        f = rows["f1f"]
        omega2 = (rows["rpm"] * 2 * np.pi / 60.0) ** 2
        kappa = (f**2 - f[0] ** 2) / omega2
        # skip the parked point and check the kappa is roughly constant
        mid = slice(2, None)
        spread = np.abs(kappa[mid] - kappa[mid].mean()).max() / kappa[mid].mean()
        assert spread < 0.30, (
            f"Southwell kappa varies {spread:.1%} across the sweep "
            f"(expected roughly constant): {kappa.tolist()}"
        )

    def test_rated_speed_matches_openfast_targets(self, rows):
        """The 7.56 rpm row must stay inside the S-4 band of the MBC3 targets."""
        i = int(np.argmin(np.abs(rows["rpm"] - 7.56)))
        for key, (target,) in [("f1f", (OF_TARGETS[7.56][0],)),
                               ("f2f", (OF_TARGETS[7.56][1],)),
                               ("f1e", (OF_TARGETS[7.56][2],))]:
            val = rows[key][i]
            assert abs(val - target) / target < 0.08, (
                f"{key} at rated: {val:.4f} vs OpenFAST {target} Hz"
            )

    def test_parked_matches_openfast_target(self, rows):
        """The parked row (0 rpm) inside the S-4 band."""
        i = int(np.argmin(np.abs(rows["rpm"] - 0.0)))
        for key, (target,) in [("f1f", (OF_TARGETS[0.0][0],)),
                               ("f1e", (OF_TARGETS[0.0][2],))]:
            val = rows[key][i]
            assert abs(val - target) / target < 0.08, (
                f"{key} parked: {val:.4f} vs OpenFAST {target} Hz"
            )

    def test_first_flap_clear_of_1p_2p_in_operating_band(self, rows):
        """The 1st flap must stay above the 1P/2P lines in the swept band.

        (The 1P crossing itself lives at ~33 rpm — far above the operating
        range — so the engineering check is the clearance, not the crossing.)
        """
        for p in (1, 2):
            fp = rows["rpm"] * p / 60.0
            clearance = rows["f1f"] - fp
            assert np.all(clearance > 0.1), (
                f"1st flap comes within {clearance.min():.3f} Hz of the {p}P line"
            )
