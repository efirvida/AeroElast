"""S-4 (windIO v2 counterpart) — IEA 15 MW blade rotating modal vs OpenFAST BeamDyn p20.

Same physics, references, metric and tolerance as ``tests/test_iea15mw_s4_rotating_modal.py`` (v1, group 40;
read its docstring for the K + K_G + K_SP formulation, the beam-projection mode classifier and the arbiter
choice): first flap / first edge / second flap at 0 rpm and 7.56 rpm against OpenFAST v5.0.0 BeamDyn single member
p20 on the official v1.1 deck (rated 1E from the BeamDyn driver), ``REL_TOL = 0.05``. ElastoDyn is printed only.

What differs from v1 (odd/tasks/v2-validation-store.md, T6): the mesh and properties come from
``tests.support.v2_blade.build_v2_blade(0.25, skin_reference="mid")`` (windIO v2 mesher, official v2 YAML) instead of
``BladeMesh(...).generate(renumber="rcm")`` + ``build_rust_properties``. The v2 mesh provides ``RootNodes`` (the clamp),
``coords_array`` and ``mesh.nodes``, which is all ``assemble_rotating`` / ``classify_modes`` read, so the v1 tool
helpers are reused unchanged.

Marks: slow (blade eigen-solves, 24-32 min per run in v1). Not run as part of T6 (framework validation postponed);
the cheap ``test_v2_s4_assembly_smoke`` exercises builder + assembly + K only.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

petsc4py = pytest.importorskip("petsc4py", reason="PETSc not available")

REPO_ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(REPO_ROOT / "tools"))
sys.path.insert(0, str(REPO_ROOT))

from tests.support.assertions import assert_relative_error  # noqa: E402

from run_s4_rotating_modal import (  # noqa: E402
    assemble_rotating,
    classify_modes,
    omega_from_rpm,
)

RPM_RATED = 7.56
ELEMENT_SIZE = 0.25

# Arbiter: OpenFAST v5.0.0 BeamDyn single member p20, official v1.1 deck, MBC3 collective (identical to v1);
# rated 1E from the BeamDyn driver (K10).
_BD_TARGETS = {
    "parked": {"1st flap": 0.507082680, "1st edge": 0.693376286, "2nd flap": 1.480900161},
    "rotating": {"1st flap": 0.532003414, "1st edge": 0.699413448, "2nd flap": 1.508057474},
}

# Comparison only: OpenFAST v5.0.0 ElastoDyn MBC3 (collective members).
_ED_COMPARISON = {
    "parked": {"1st flap": 0.544, "1st edge": 0.739, "2nd flap": 1.610},
    "rotating": {"1st flap": 0.5666, "1st edge": 0.7461, "2nd flap": 1.6372},
}

REL_TOL = 0.05  # same as v1 (tests/test_iea15mw_s4_rotating_modal.py)


def _build_assembly(element_size=ELEMENT_SIZE):
    from aeroelast.core.assembler import MeshAssembler
    from aeroelast.elements import ElementFamily
    from tests.support.v2_blade import build_v2_blade

    mesh, props = build_v2_blade(element_size, skin_reference="mid")
    model = {
        "elements": {
            "element_family": ElementFamily.SHELL,
            "span_direction": (0.0, 0.0, 1.0),
            "properties": props,
        }
    }
    return MeshAssembler(mesh, model), mesh


def _family_frequency(freqs, modes, free, mesh, direction, order):
    """Return the n-th lowest mode frequency whose beam-projected motion is `direction`."""
    labels = classify_modes(mesh, modes, free)
    candidates = sorted(f for f, lab in zip(freqs, labels, strict=False) if lab == direction)
    if len(candidates) < order:
        raise AssertionError(f"only {len(candidates)} {direction} mode(s) found; need {order}")
    return float(candidates[order - 1])


def test_v2_s4_assembly_smoke():
    """Cheap fixture check (no eigen-solve): v2 builder + MeshAssembler + K, clamp set and mesh orientation."""
    import numpy as np

    asm, mesh = _build_assembly()
    assert asm.dofs_per_node == 6
    root = mesh.get_node_set("RootNodes").nodes
    assert len(root) > 0
    z = mesh.coords_array[:, 2]
    # the root ring sits at z ~ -0.017 m (mid-surface inset of the v2 skin), tip at ~ +117 m
    assert z.min() == pytest.approx(0.0, abs=0.05), "span must run along +z from the root"
    assert z.max() == pytest.approx(117.0, abs=1.0)
    rows, cols, vals = asm._rust.assemble_k()
    assert len(vals) > 0 and np.all(np.isfinite(np.asarray(vals)))
    assert omega_from_rpm(RPM_RATED) == pytest.approx(RPM_RATED * 2 * np.pi / 60)


def _measure_families(rpm):
    """Shell 1F/1E/2F [Hz] at `rpm`; every family and its comparison numbers are printed before any assertion,
    so a red family does not hide the ones after it (v1 kept the same order of evidence)."""
    key = "parked" if rpm == 0.0 else "rotating"
    asm, mesh = _build_assembly()
    freqs, modes, free = assemble_rotating(asm, mesh, omega_from_rpm(rpm), n_modes=12)
    out = {}
    for family, direction, order in [
        ("1st flap", "flap", 1),
        ("1st edge", "edge", 1),
        ("2nd flap", "flap", 2),
    ]:
        f_shell = _family_frequency(freqs, modes, free, mesh, direction, order)
        f_ed = _ED_COMPARISON[key][family]
        f_bd = _BD_TARGETS[key][family]
        print(
            f"S-4 v2 {key} {family} [Hz]: shell {f_shell:.6f} | BeamDyn p20 {f_bd:.6f} "
            f"({f_shell / f_bd - 1:+.4%}, bound {REL_TOL:.0%}) | ElastoDyn comparison only {f_ed:.6g} "
            f"({f_shell / f_ed - 1:+.4%})"
        )
        out[family] = f_shell
    return out


@pytest.mark.slow
def test_v2_s4_parked_modal_vs_openfast():
    """Parked (0 rpm) first flap, first edge and second flap vs OpenFAST BeamDyn p20 (ElastoDyn printed)."""
    f = _measure_families(0.0)
    assert_relative_error(
        f["1st flap"], _BD_TARGETS["parked"]["1st flap"], tol=REL_TOL,
        reference_name="OpenFAST v5.0.0 BeamDyn p20 (official v1.1 deck) parked first flap, Hz",
        kind="code", what="S-4 v2 parked 1F",
    )
    assert_relative_error(
        f["1st edge"], _BD_TARGETS["parked"]["1st edge"], tol=REL_TOL,
        reference_name="OpenFAST v5.0.0 BeamDyn p20 (official v1.1 deck) parked first edge, Hz",
        kind="code", what="S-4 v2 parked 1E",
    )
    assert_relative_error(
        f["2nd flap"], _BD_TARGETS["parked"]["2nd flap"], tol=REL_TOL,
        reference_name="OpenFAST v5.0.0 BeamDyn p20 (official v1.1 deck) parked second flap, Hz",
        kind="code", what="S-4 v2 parked 2F",
    )


@pytest.mark.slow
def test_v2_s4_rated_modal_vs_openfast():
    """Rated (7.56 rpm) first flap, first edge and second flap vs OpenFAST BeamDyn p20 / BeamDyn driver."""
    f = _measure_families(RPM_RATED)
    assert_relative_error(
        f["1st flap"], _BD_TARGETS["rotating"]["1st flap"], tol=REL_TOL,
        reference_name="OpenFAST v5.0.0 BeamDyn p20 (official v1.1 deck) rated first flap at 7.56 rpm, Hz",
        kind="code", what="S-4 v2 rated 1F",
    )
    assert_relative_error(
        f["1st edge"], _BD_TARGETS["rotating"]["1st edge"], tol=REL_TOL,
        reference_name="OpenFAST v5.0.0 BeamDyn driver, dynamic, p20 (official v1.1 deck, hub 0) rated first edge at 7.56 rpm, Hz",
        kind="code", what="S-4 v2 rated 1E",
    )
    assert_relative_error(
        f["2nd flap"], _BD_TARGETS["rotating"]["2nd flap"], tol=REL_TOL,
        reference_name="OpenFAST v5.0.0 BeamDyn p20 (official v1.1 deck) rated second flap at 7.56 rpm, Hz",
        kind="code", what="S-4 v2 rated 2F",
    )


@pytest.mark.slow
def test_v2_s4_centrifugal_stiffening_increases_first_flap():
    """The 1st flap frequency must rise from parked to rated rotor speed."""
    asm, mesh = _build_assembly()
    f_parked, modes_p, free_p = assemble_rotating(asm, mesh, 0.0, n_modes=8)
    f_rot, modes_r, free_r = assemble_rotating(asm, mesh, omega_from_rpm(RPM_RATED), n_modes=8)
    flap_p = _family_frequency(f_parked, modes_p, free_p, mesh, "flap", 1)
    flap_r = _family_frequency(f_rot, modes_r, free_r, mesh, "flap", 1)
    assert flap_r > flap_p, (
        f"centrifugal stiffening missing: 1st flap {flap_p:.4f} Hz (parked) -> {flap_r:.4f} Hz (7.56 rpm)"
    )
