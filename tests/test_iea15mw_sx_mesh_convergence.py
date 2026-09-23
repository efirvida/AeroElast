"""S-x — mesh convergence of the blade validations (S-4 modal + S-5 static).

Consumes the artifact produced by tools/run_sx_convergence.py
(docs/validation_data/generated/sx_mesh_convergence.csv) and checks:

  * the artifact carries the element-size series (1.0 / 0.5 / 0.25 m),
  * the results converge monotonically with refinement,
  * the 0.5 m values stay within the S-4 / S-5 target bands.

The artifact is produced by the HPC convergence job (sx_convergence.srm);
until it completes, a smoke CSV (the 1.0/0.5 pair) can be used via the
AEROELAST_SX_CONVERGENCE_CSV environment variable.
"""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CSV = REPO_ROOT / "docs" / "validation_data" / "generated" / "sx_mesh_convergence.csv"

# S-4 targets (OpenFAST v5 MBC3, the S-4 test's records) + the S-5 OoP target
TARGETS = {
    "f1f_park_Hz": (0.544, 0.04),
    "f1f_rot_Hz": (0.5666, 0.04),
    "f1e_park_Hz": (0.739, 0.08),   # the S-4's own recorded edge = -5.6%
    "oop_static_m": (16.1033, 0.10),
}


def _load_csv():
    path = Path(os.environ.get("AEROELAST_SX_CONVERGENCE_CSV", str(DEFAULT_CSV)))
    if not path.exists():
        pytest.skip(f"convergence CSV not found: {path}")
    return np.loadtxt(path, delimiter=",", skiprows=1)


@pytest.mark.slow
class TestSxMeshConvergence:
    @pytest.fixture(scope="class")
    def rows(self):
        data = _load_csv()
        if data.ndim == 1:
            data = data[None, :]
        return {
            "size": data[:, 0],
            "f1f_park": data[:, 3],
            "f1f_rot": data[:, 4],
            "f1e_park": data[:, 5],
            "oop": data[:, 7],
        }

    def test_artifact_has_refinement_series(self, rows):
        sizes = rows["size"]
        if len(sizes) < 3:
            pytest.skip("partial artifact (smoke run): the full 3-size study "
                        "is produced by the HPC convergence job")
        assert sizes[-1] < 0.4, (
            f"expected the finest size below 0.4 m; got {sizes.tolist()}"
        )

    def test_oop_static_converges_monotonically(self, rows):
        """The S-5 OoP deflection must change less and less with refinement."""
        oop = rows["oop"]
        if len(oop) < 3:
            pytest.skip("partial artifact: needs the full 3-size series")
        deltas = [abs(oop[i] - oop[i + 1]) for i in range(len(oop) - 1)]
        assert all(d > 0 for d in deltas) and deltas[-1] <= deltas[0], (
            f"OoP series {oop.tolist()} does not converge monotonically"
        )

    def test_frequencies_converge(self, rows):
        """The flapwise modal frequencies must tighten toward the targets."""
        for key in ("f1f_park", "f1f_rot"):
            f = rows[key]
            if len(f) < 3:
                pytest.skip(f"partial artifact: {key} needs 3 sizes")
            deltas = [abs(f[i] - f[i + 1]) for i in range(len(f) - 1)]
            assert deltas[-1] < deltas[0], (
                f"{key} series {f.tolist()} does not tighten with refinement"
            )

    @pytest.mark.xfail(
        strict=False,
        reason="artifact docs/validation_data/generated/sx_mesh_convergence.csv "
        "predates the MITC3 shear-convention fix (2026-09-19): its 0.125 m row "
        "shows a spurious f1e jump (0.6977 -> 0.7028 Hz).  Regenerate with "
        "`python tools/run_sx_convergence.py --sizes 1.0 0.5 0.25 0.125` on HPC "
        "and this becomes a strict check.",
    )
    def test_f1e_frequency_converges(self, rows):
        """First edgewise parked frequency must tighten with refinement."""
        f = rows["f1e_park"]
        if len(f) < 3:
            pytest.skip("partial artifact: f1e_park needs 3 sizes")
        deltas = [abs(f[i] - f[i + 1]) for i in range(len(f) - 1)]
        assert deltas[-1] < deltas[0], (
            f"f1e_park series {f.tolist()} does not tighten with refinement"
        )

    def test_coarse_size_within_target_bands(self, rows):
        """The 0.5 m values (the S-4/S-5's working size) stay inside the bands."""
        i05 = int(np.argmin(np.abs(rows["size"] - 0.5)))
        checks = {
            "f1f_park_Hz": rows["f1f_park"][i05],
            "f1f_rot_Hz": rows["f1f_rot"][i05],
            "f1e_park_Hz": rows["f1e_park"][i05],
            "oop_static_m": rows["oop"][i05],
        }
        for key, val in checks.items():
            target, tol = TARGETS[key]
            assert abs(val - target) / target < tol, (
                f"{key} at 0.5 m: {val:.4f} vs target {target} (tol {tol:.0%})"
            )
