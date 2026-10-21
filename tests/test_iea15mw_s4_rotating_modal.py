"""S-4 — IEA 15 MW blade: rotating modal (K + K_G + K_SP) vs OpenFAST targets.

Computes the root-clamped rotating-frame modal of the official WindIO blade
(tests/IEA-15-240-RWT.yaml, MITC4 shell with a few MITC3 (oml: 33721 quads, 92 triangles),
element_size 0.25, consistent mass)
at the rated rotor speed (7.56 rpm, rotation axis +y, span +z) with:

  K_eff = K + K_G + K_SP
  K_G   = geometric stiffness from the centrifugal pre-stress static solve
          (K u = f_cf, root-clamped; membrane stresses recovered from u)
  K_SP  = spin softening -w^2 m_i (I - n n^T) on translational DOFs

Modes are classified by projecting the whole mode on beam kinematics (rigid cross-sections
per span bin, ``classify_modes``): flap = u_y (out-of-plane, chord lies along x in the rotor
plane), edge = u_x; a mode the projection explains below 60 % is "local" and never enters a
family. Until 2026-10-18 the label came from the single tip node, which took a local mode
(oml, es 0.25, 7.56 rpm, 1.5325 Hz, beam fraction 0.15) as 2nd flap
(odd/tasks/s4-mode-classifier.md).

Arbiter (odd/tasks/validation-l3-b3-audit.md K8, 2026-10-20): OpenFAST v5.0.0 BeamDyn, single
member, element order 20, on the official v1.1 deck (gravity 0, prescribed rotor speed, MBC3
collective members; odd/tasks/openfast-beamdyn-modal-arbiter.md, SLURM 11617531/11617534).
Converged in order: p15 -> p20 moves 1F/1E/2F by -0.03 / +0.13 / +0.07 % parked and
-0.03 / +0.12 / +0.07 % rated. BeamDyn solves the geometrically exact Timoshenko beam on the
full 6x6 deck; ElastoDyn's blade modes come from fitted mode-shape polynomials on an
Euler-Bernoulli beam, so ElastoDyn is printed as a comparison and decides nothing.

  mode        parked (0 rpm)            rotating (7.56 rpm)
              BeamDyn p20 | ElastoDyn   BeamDyn p20 | ElastoDyn
  1st flap    0.507083    | 0.544       0.532003    | 0.5666
  1st edge    0.693376    | 0.739       0.699413*   | 0.7461
  2nd flap    1.480900    | 1.610       1.508057    | 1.6372

* Rotating 1E does not come from the MBC3 linearization: its rotating/parked rise (+2.66 %)
and edge Southwell coefficient (1.64, equal to flap) there are an artefact of the BeamDyn
linearization (odd/tasks/beamdyn-arbiter-midsurface.md T1/T7). It comes from the OpenFAST
v5.0.0 BeamDyn driver, dynamic, p20, official v1.1 deck, hub 0, 7.56 rpm, free vibration
fitted on the tip response (tools/bd_rotating_freevib.py; validation-l3-b3-audit.md K10,
SLURM 11621099). The method passes the K10.2 gate on a uniform and a published box beam
(parked <= 0.01 %, Southwell K <= 2.8 % from the analytic); on this blade its parked 1F/1E
reproduce BeamDyn p20 MBC3 within 0.03 %, and p15 -> p20 moves rated 1E +0.12 %,
dt -> dt/2 +0.004 %, halving the load -0.005 %. Hub radius 3.97 m (the MBC3 deck's) would
move it +0.12 %; the shell modal rotates about the blade root.

Measured agreement vs ElastoDyn, kept as history (2026-10-18, mid-surface skin, es 0.25,
tools/diagnose_mode_beam_projection.py, SLURM 11618419). The beam-projection labels pick the same three modes the
tip-node labels picked on this mesh, so the numbers did not move with the classifier:

  rotating:  1F -3.40 | 1E -6.09 | 2F -6.13 %     (es 0.125: -3.48 | -6.10 | -6.46 %, T12)
  parked:    1F -3.56 | 1E -5.94 | 2F -6.34 %

After the K_G prestress fix 41e9a16 ((A eps + B kappa)/h; SLURM 11619165, mid es 0.25) the
rotating row reads 1F -3.41 | 1E -6.10 | 2F -6.37 % (2F moved -0.25 % in frequency, the rest
< 0.02 %); parked is unchanged by construction (no K_G at 0 rpm); the es 0.125 values are
pre-fix. The rotating-only tip local modes (2.17-3.61 Hz) are gone with the fix
(odd/tasks/rotating-local-tip-modes.md L39).

Against ElastoDyn 1E and 2F exceed the 5 % bound at both rotor speeds. Against BeamDyn p20
(parked es 0.25, 19ff1a3): 1F +3.42 | 1E +0.20 | 2F +1.80 %; es 0.125 (SLURM 11619955):
1F +3.31 | 1E +0.29 | 2F +1.66 %. The 6-8 % ElastoDyn/BeamDyn disagreement is a finding about
the references, recorded and not used to move the bound (docs/validation_closures.md,
"Superficie media en las anclas de pala"). The K8 choice was made after these numbers were
known; it rests on the model physics above, not on which side passes.

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

# Arbiter: OpenFAST v5.0.0 BeamDyn single member p20, official v1.1 deck, MBC3 collective
# (SLURM 11617531 parked, 11617534 rated); rated 1E from the BeamDyn driver (see above).
_BD_TARGETS = {
    "parked": {"1st flap": 0.507082680, "1st edge": 0.693376286, "2nd flap": 1.480900161},
    "rotating": {"1st flap": 0.532003414, "1st edge": 0.699413448, "2nd flap": 1.508057474},
}
_BD_DRIVER = {("rotating", "1st edge")}

# Comparison only: OpenFAST v5.0.0 ElastoDyn MBC3 (collective members), measured 2026-09-09.
_ED_COMPARISON = {
    "parked": {"1st flap": 0.544, "1st edge": 0.739, "2nd flap": 1.610},
    "rotating": {"1st flap": 0.5666, "1st edge": 0.7461, "2nd flap": 1.6372},
}

REL_TOL = 0.05


def _build_assembly(element_size=0.25):
    from aeroelast.core.assembler import MeshAssembler
    from aeroelast.core.mesh.generators import BladeMesh
    from aeroelast.elements import ElementFamily
    from aeroelast.models.blade.model import build_rust_properties

    gen = BladeMesh(
        element_size=element_size,
        yaml_file=str(REPO_ROOT / "tests/IEA-15-240-RWT.yaml"),
        skin_reference="mid",  # blade skin on the mid-surface (odd/tasks/beamdyn-arbiter-midsurface.md, T4)
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
    """First flap/edge and second flap modal vs OpenFAST BeamDyn p20 (ElastoDyn printed)."""
    key = "parked" if rpm == 0.0 else "rotating"
    asm, mesh = _build_assembly()
    freqs, modes, free = assemble_rotating(asm, mesh, omega_from_rpm(rpm), n_modes=12)

    # Every family is computed and printed before anything is asserted, so a red family
    # (1E) does not hide the ones after it (2F was never measured in a red run before
    # B3-2 batch 1, odd/tasks/validation-l3-b3-audit.md K6).
    failures = []
    for family, direction, order in [
        ("1st flap", "flap", 1),
        ("1st edge", "edge", 1),
        ("2nd flap", "flap", 2),
    ]:
        f_shell = _family_frequency(freqs, modes, free, mesh, direction, order)
        f_ed = _ED_COMPARISON[key][family]
        print(
            f"S-4 {key} {family} [Hz] vs OpenFAST v5.0.0 ElastoDyn MBC3 (comparison, not "
            f"asserted): {f_shell:.6f} vs {f_ed:.6g} (signed error: {f_shell / f_ed - 1:+.4%})"
        )
        f_bd = _BD_TARGETS[key][family]
        if f_bd is None:
            print(f"S-4 {key} {family}: no admissible BeamDyn p20 number, not asserted")
            continue
        rel = abs(f_shell - f_bd) / f_bd
        source = "BeamDyn driver p20" if (key, family) in _BD_DRIVER else "BeamDyn p20"
        print(
            f"S-4 {key} {family} [Hz] vs OpenFAST v5.0.0 {source}: {f_shell:.6f} vs "
            f"{f_bd:.6f} (signed error: {f_shell / f_bd - 1:+.4%}, bound: {REL_TOL:.4%})"
        )
        if rel >= REL_TOL:
            failures.append(
                f"{family} ({key}, {rpm} rpm): shell {f_shell:.4f} Hz vs "
                f"{source} {f_bd:.4f} Hz — rel {rel:.1%} > {REL_TOL:.0%}"
            )
    assert not failures, "; ".join(failures)


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
