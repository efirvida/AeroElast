"""V-05 on the windIO v2 blade: assembled shell mass vs the BeamDyn deck (v2 counterpart of group 42).

Same reference, metric and tolerance as ``tests/test_iea15mw_v05_structural_properties.py``:

* mass: row-sum of the assembled M (x-DOF rows) vs the integral of M11 of the official v1.1 BeamDyn deck over the
  reference-axis arc length (66,996.9 kg), ``rel_err < 0.05`` (L3 ceiling);
* guards (not store rows in v1): K and M assemble, K size is a multiple of 6, n_nodes >= 50.

Only the mesh differs: ``build_v2_blade(0.25, skin_reference="mid")`` instead of
``BladeMesh(...).generate(renumber="rcm")`` + ``build_rust_properties``. The v1 class
``TestBeamDynK55K66Consistency`` is deck-vs-deck (no aeroelast code, mesher independent) and is not repeated.
See docs/validation/_pending/67-v2_blade_structural_properties/NOTES.md.
"""

from __future__ import annotations

import re
import sys

import pytest

pytest.importorskip("petsc4py", reason="PETSc not available")

from tests.support.assertions import assert_relative_error  # noqa: E402
from tests.support.paths import REPO_ROOT  # noqa: E402
from tests.support.v2_blade import build_v2_blade  # noqa: E402

sys.path.insert(0, str(REPO_ROOT / "tools"))

_BD_PROPERTIES = REPO_ROOT / "tests/IEA15MW/reference/IEA-15-240-RWT_BeamDyn_blade.dat"
# The published v1.1 reference axis (byte-identical to the deck the BeamDyn driver runs).
_BD_GEOMETRY = (
    REPO_ROOT / "tests/IEA15MW/validation papers/github-IEA-15-240-RWT/OpenFAST/IEA-15-240-RWT/IEA-15-240-RWT_BeamDyn.dat"
)
ELEMENT_SIZE = 0.25
SKIN_REFERENCE = "mid"


@pytest.mark.slow
class TestAssembledStiffnessPropertiesV2:
    """V-05c on the v2 blade: assembled FEM model matrix-level checks."""

    @pytest.fixture(scope="class")
    def assembler(self):
        from aeroelast.core.assembler import MeshAssembler
        from aeroelast.elements import ElementFamily

        mesh, properties = build_v2_blade(ELEMENT_SIZE, skin_reference=SKIN_REFERENCE)
        model = {
            "elements": {
                "element_family": ElementFamily.SHELL,
                "span_direction": (0.0, 0.0, 1.0),
                "properties": properties,
            },
        }
        return MeshAssembler(mesh, model)

    def test_assembler_creates_valid_matrices(self, assembler):
        assert assembler.assemble_stiffness_matrix() is not None
        assert assembler.assemble_mass_matrix() is not None

    def test_stiffness_matrix_size_consistent_with_dof_count(self, assembler):
        n_dof = assembler.assemble_stiffness_matrix().getSize()[0]
        assert n_dof % 6 == 0, f"K matrix size {n_dof} is not a multiple of 6"
        n_nodes = n_dof // 6
        assert n_nodes >= 50, f"Only {n_nodes} nodes — blade mesh seems too coarse"

    def test_mass_matrix_total_consistent_with_reference(self, assembler):
        """Total mass from M within 5 % of the official v1.1 BeamDyn deck mass (66,996.9 kg)."""
        import numpy as np
        from beam_reference import beamdyn_deck_mass, load_beamdyn_blade, reference_axis_length

        # Σ_j M_ij over x-DOF rows = m_node (partition of unity), exact for consistent mass.
        m_rows, _, m_vals = assembler._rust.assemble_m()
        m_rows = np.asarray(m_rows)
        total_mass = float(np.asarray(m_vals)[m_rows % 6 == 0].sum())

        lines = _BD_GEOMETRY.read_text().splitlines()
        first = next(k for k, line in enumerate(lines) if "kp_xr" in line) + 2
        n_kp = int(re.search(r"(\d+)\s+kp_total", "\n".join(lines)).group(1))
        kp = np.array([[float(v) for v in lines[first + k].split()[:3]] for k in range(n_kp)])
        deck = load_beamdyn_blade(_BD_PROPERTIES, blade_length=1.0)  # z = eta
        ref_mass = beamdyn_deck_mass(deck.z, deck.mass_dens, reference_axis_length(kp))
        assert_relative_error(
            total_mass,
            ref_mass,
            tol=0.05,
            reference_name="BeamDyn v1.1 deck mass (integral of M11 over the 117.149 m reference axis)",
            kind="code",
            what="V-05 v2 assembled shell mass [kg]",
        )
