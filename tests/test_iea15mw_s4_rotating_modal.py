"""S-4 — IEA 15 MW blade: rotating modal (K + K_G + K_SP) vs OpenFAST targets.

Computes the root-clamped rotating-frame modal of the official WindIO blade
(tests/IEA-15-240-RWT.yaml, MITC3 shell, element_size 0.5, consistent mass)
at the rated rotor speed (7.56 rpm, rotation axis +y, span +z) with:

  K_eff = K + K_G + K_SP
  K_G   = geometric stiffness from the centrifugal pre-stress static solve
          (K u = f_cf, root-clamped; membrane stresses recovered from u)
  K_SP  = spin softening -w^2 m_i (I - n n^T) on translational DOFs

Modes are classified by the tip-node translation direction: flap = u_y
(out-of-plane, chord lies along x in the rotor plane), edge = u_x.

Reference targets: OpenFAST v5.0.0, ElastoDyn beam model, blade DOFs only
(no aero, no tower/platform DOFs), MBC3 rotating frame, collective members
(single-blade comparison; the shell model has no whirl pairs):

  mode        parked (0 rpm)      rotating (7.56 rpm)
  1st flap    0.544 Hz            0.5666 Hz
  1st edge    0.739 Hz            0.7461 Hz
  2nd flap    1.610 Hz            1.6372 Hz

Measured agreement (2026-09-09, this test):

  rotating:  1F +0.6% | 1E -4.5% | 2F -1.5%
  parked:    1F -1.5% | 1E -5.6% | 2F -2.3%

Marks: slow (blade mesh generation + static solve + modal solve).
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

petsc4py = pytest.importorskip("petsc4py", reason="PETSc not available")

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "tools"))

from run_s4_rotating_modal import (  # noqa: E402
    assemble_rotating,
    classify_modes,
    omega_from_rpm,
)

RPM_RATED = 7.56

# OpenFAST v5.0.0 MBC3 targets (collective members) — measured 2026-09-09.
_OF_TARGETS = {
    "parked": {"1st flap": 0.544, "1st edge": 0.739, "2nd flap": 1.610},
    "rotating": {"1st flap": 0.5666, "1st edge": 0.7461, "2nd flap": 1.6372},
}

REL_TOL = 0.08


def _build_assembly(element_size=0.5):
    from aeroelast.core.assembler import MeshAssembler
    from aeroelast.core.mesh.generators import BladeMesh
    from aeroelast.elements import ElementFamily
    from aeroelast.models.blade.model import build_rust_properties

    gen = BladeMesh(
        element_size=element_size,
        yaml_file=str(REPO_ROOT / "tests/IEA-15-240-RWT.yaml"),
    )
    mesh = gen.generate(renumber="rcm")
    props = build_rust_properties(gen.numad_mesh_data)
    model = {
        "elements": {
            "element_family": ElementFamily.SHELL,
            "span_direction": (0.0, 0.0, 1.0),
            "properties": props,
        }
    }
    return MeshAssembler(mesh, model), mesh


def _family_frequency(freqs, modes, free, mesh, direction, order):
    """Return the n-th lowest mode frequency whose tip motion is `direction`."""
    labels = classify_modes(mesh, modes, free)
    candidates = sorted(f for f, lab in zip(freqs, labels) if lab == direction)
    if len(candidates) < order:
        raise AssertionError(
            f"only {len(candidates)} {direction} mode(s) found; need {order}"
        )
    return float(candidates[order - 1])


@pytest.mark.slow
@pytest.mark.parametrize("rpm", [0.0, RPM_RATED])
def test_s4_rotating_modal_vs_openfast(rpm):
    """First flap/edge and second flap rotating modal vs OpenFAST MBC3 targets."""
    key = "parked" if rpm == 0.0 else "rotating"
    asm, mesh = _build_assembly()
    freqs, modes, free = assemble_rotating(asm, mesh, omega_from_rpm(rpm), n_modes=12)

    for family, direction, order in [
        ("1st flap", "flap", 1),
        ("1st edge", "edge", 1),
        ("2nd flap", "flap", 2),
    ]:
        f_shell = _family_frequency(freqs, modes, free, mesh, direction, order)
        f_of = _OF_TARGETS[key][family]
        rel = abs(f_shell - f_of) / f_of
        assert rel < REL_TOL, (
            f"{family} ({key}, {rpm} rpm): shell {f_shell:.4f} Hz vs "
            f"OpenFAST {f_of:.4f} Hz — rel {rel:.1%} > {REL_TOL:.0%}"
        )


@pytest.mark.slow
def test_s4_centrifugal_stiffening_increases_first_flap():
    """The 1st flap frequency must rise from parked to rated rotor speed."""
    asm, mesh = _build_assembly()
    f_parked, modes_p, free_p = assemble_rotating(asm, mesh, 0.0, n_modes=8)
    f_rot, modes_r, free_r = assemble_rotating(
        asm, mesh, omega_from_rpm(RPM_RATED), n_modes=8
    )
    flap_p = _family_frequency(f_parked, modes_p, free_p, mesh, "flap", 1)
    flap_r = _family_frequency(f_rot, modes_r, free_r, mesh, "flap", 1)
    assert flap_r > flap_p, (
        f"centrifugal stiffening missing: 1st flap {flap_p:.4f} Hz (parked) -> "
        f"{flap_r:.4f} Hz (7.56 rpm)"
    )
