"""V-05 — IEA 15 MW RWT: distributed structural properties validation.

Validates that the stiffness and mass properties of the assembled blade model
agree with the reference data from the IEA 15 MW definition files.

These tests are purely parametric — they do NOT run a full solve.  They check
the assembled matrices K and M against known diagonal bounds at specific
spanwise stations.

Reference:
  [R3] IEA-15-240-RWT_tabular.xlsx, sheet "Blade Structural Properties"
       (BeamDyn 6×6 matrices, 26 spanwise stations)
  [R4] IEA-15-240-RWT_ElastoDyn_blade.dat
       (1D distributed properties: BMassDen, FlpStff, EdgStff)
  [R6] IEA-15-240-RWT_BeamDyn_blade.dat
       (full Timoshenko 6×6 K and M per station)

Convention mapping:
  BeamDyn K_55 = flapwise bending stiffness (EI_flap)  [R6]
  BeamDyn K_66 = edgewise bending stiffness (EI_edge)  [R6]
  BeamDyn M_11 = mass per unit length [kg/m]           [R6]
  aeroelast DOFs: [ux, uy, uz, rx, ry, rz] per node   [AGENTS.md]
"""

from __future__ import annotations

import pytest

# ---------------------------------------------------------------------------
# Reference: flapwise bending stiffness EI_flap at selected r/R stations
# Source: [R4] IEA-15-240-RWT_ElastoDyn_blade.dat, column FlpStff [N·m²]
# Values at r/R = 0.00, 0.10, 0.20, 0.30, 0.50, 0.70, 0.90, 1.00
# ---------------------------------------------------------------------------
_EI_FLAP_REF = {
    # r/R: EI_flap (N·m²)  [R4]
    0.00: 1.525e11,
    0.10: 5.671e10,
    0.20: 2.462e10,
    0.30: 1.408e10,
    0.50: 4.925e9,
    0.70: 1.226e9,
    0.90: 1.179e8,
}

# Reference: BeamDyn K_55 (flapwise, Timoshenko) at key stations [R3][R6]
_K55_REF = {
    0.00: 1.497e11,
    0.10: 5.226e10,
    0.50: 4.893e9,
    0.90: 1.133e8,
}


class TestBeamDynK55K66Consistency:
    """V-05b — reference-deck consistency record: BeamDyn K_55 vs ElastoDyn FlpStff.

    NOT a validation of aeroelast: no aeroelast code runs here, so the class is outside the
    L3 contract (validation store group 42 lists it under ``non_validation_tests``). It records
    that the two official decks of the same blade disagree. The hard-coded ElastoDyn values
    labelled r/R 0.10/0.50/0.90 are the rows at BlFract 0.102/0.510/0.898 while the BeamDyn
    values are at eta 0.10/0.50/0.90, so stations differ; compared at the same eta the decks
    differ by up to 24 % (flap) and 39 % (edge)
    (tools/v05_deck_consistency.py, logs/check/v05/deck_consistency.txt;
    odd/tasks/validation-l3-b3-audit.md, "B3-2 finding"). The 5 % assertion is kept, not
    widened; the r/R=0.1 case fails (7.8 %), the others pass only by station choice.
    """

    @pytest.mark.xfail(
        strict=False,
        reason="reference decks disagree by more than 5 % (BeamDyn vs ElastoDyn, same blade); "
        "documented deck spread, not an aeroelast result",
    )
    @pytest.mark.parametrize("rR, k55", sorted(_K55_REF.items()))
    def test_k55_consistent_with_ei_flap(self, rR, k55):
        """BeamDyn K_55 within 5 % of ElastoDyn EI_flap [R3][R4][R6] (deck-vs-deck record)."""
        ei_flap = _EI_FLAP_REF.get(rR)
        if ei_flap is None:
            pytest.skip(f"No ElastoDyn reference at r/R={rR}")

        rel_diff = abs(k55 - ei_flap) / ei_flap
        assert rel_diff < 0.05, (
            f"r/R={rR}: BeamDyn K_55={k55:.3e} vs ElastoDyn EI_flap={ei_flap:.3e} "
            f"(diff={rel_diff:.1%}) [R3][R4][R6]"
        )


@pytest.mark.slow
class TestAssembledStiffnessProperties:
    """V-05c — assembled FEM model matrix-level checks [R3][R4]."""

    @pytest.fixture(scope="class")
    def assembler(self, iea_blade_yaml):
        pytest.importorskip("petsc4py")
        from aeroelast.core.assembler import MeshAssembler
        from aeroelast.elements import ElementFamily
        from aeroelast.core.mesh.generators import BladeMesh
        from aeroelast.models.blade.model import build_rust_properties

        generator = BladeMesh(
            yaml_file=iea_blade_yaml,
            element_size=0.25,
            skin_reference="mid",  # blade skin on the mid-surface (odd/tasks/beamdyn-arbiter-midsurface.md, T4)
        )
        mesh = generator.generate(renumber="rcm")
        properties = build_rust_properties(generator.numad_mesh_data)
        model = {
            "elements": {
                "element_family": ElementFamily.SHELL,
                "span_direction": (0.0, 0.0, 1.0),
                "properties": properties,
            },
        }
        return MeshAssembler(mesh, model)

    def test_assembler_creates_valid_matrices(self, assembler):
        """K and M must be assembled without errors [R3][R4]."""
        K = assembler.assemble_stiffness_matrix()
        M = assembler.assemble_mass_matrix()
        assert K is not None, "Stiffness matrix assembly returned None"
        assert M is not None, "Mass matrix assembly returned None"

    def test_stiffness_matrix_size_consistent_with_dof_count(self, assembler):
        """K matrix size must match n_nodes × 6 DOFs per node [AGENTS.md]."""
        K = assembler.assemble_stiffness_matrix()
        size = K.getSize()
        n_dof = size[0]
        assert n_dof % 6 == 0, (
            f"K matrix size {n_dof} is not a multiple of 6 — expected 6 DOFs/node [AGENTS.md]"
        )
        n_nodes = n_dof // 6
        assert n_nodes >= 50, f"Only {n_nodes} nodes — blade mesh seems too coarse"

    def test_mass_matrix_total_consistent_with_reference(self, assembler):
        """Total mass from M within the documented +4-8% band of 67,921 kg.

        The assembled (layup-based) shell mass sits +4-8% above the beam
        tables: the WindIO yaml publishes two structural representations
        (internal_structure_2d_fem layup vs elastic_properties_mb 6x6) and
        the layup is the heavier one by construction.  Measured +4.1%.
        See docs/shell_vs_beam_sectional_validation.md sections 2 and 5.
        """
        # Total mass = sum of all row entries for x-DOFs (partition of unity identity).
        # For consistent mass matrices: Σ_j M_ij = m_node for x-translation DOFs.
        # This is exact; diag/3 only works for lumped mass matrices.
        m_rows, m_cols, m_vals = assembler._rust.assemble_m()
        dpn = 6
        x_row_mask = m_rows % dpn == 0
        total_mass = float(m_vals[x_row_mask].sum())
        ref_mass = 67_921.0  # kg [R3][R6] beam-table representation
        rel_err = abs(total_mass - ref_mass) / ref_mass
        assert rel_err < 0.05, (
            f"Assembled mass {total_mass:.0f} kg vs beam-table ref {ref_mass:.0f} kg "
            f"(err={rel_err:.1%}) — outside the documented +4-8% layup-vs-6x6 band "
            f"[R3][R6]"
        )
