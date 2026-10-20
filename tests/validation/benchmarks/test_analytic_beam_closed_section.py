"""T6/S9/S8: independent non-blade diagnostics, never xfailed.

Run: scripts/aeroenv.sh python -m pytest -q <this file> -s
Three mesh levels are retained in each test's diagnostic output. No source
change is part of this task: a converged red assertion is an implementation
finding, not permission to widen a tolerance. Explicit OML geometry is a model
experiment, not a test of the BladeMesh reference option (already covered).

References and conventions live beside the independent formulas in
analytic_beam_helpers.py. Static refinement 160/320/640 (four elements/wall)
shows quadratic span error: the last EI step is 0.38%, giving a Richardson
tail of 0.13%. A 0.5% bound includes that tail plus thin-wall/local-wall and
end-restraint corrections (extrapolated residual ~0.23%). EA's error is
0.00033% and mesh movement <1e-10: the 0.001% bound covers the residual
end restraint. Mass is an exact polynomial integral, roundoff only.

Modal meshes refine BOTH directions: (span, wall)=(160,8),(240,12),(320,16).
Final parked bending errors vs Timoshenko are 0.40-0.48%, last changes <=0.27%,
with quadratic extrapolation tails <=0.35%; torsion error is 0.0173% and its
last change 0.0050%. These trends license the pre-set 2% frequency bound,
without widening it. Timoshenko separates shear+rotary inertia from EB:
first-mode shifts ~0.2%, second-mode shifts 1.27-1.44%.

The native LU graph-ordering hang is bypassed ONLY by a scoped RCM solver
option (no changed matrices). The Southwell tests were RED until 43dce4d:
isotropic centrifugal loading received rho rather than rho*t, giving exactly
20x the analytic resultant for t=.05 (K_flap 23.1293 vs 1.19103, K_lag 21.6663
vs 0.188469). With rho*t they pass the strict pre-set 3% bound: K_flap 0.19 %,
K_lag 0.60 % from Wright et al. 1982 (measured 2026-10-19, SLURM 11619328).

The shared modal convergence fixture takes ~551s on sequana_cpu, so this file
belongs in SLURM, not a login-node smoke run. Static/mass checks remain fast.
"""

import time

import numpy as np
import pytest

pytest.importorskip("petsc4py", reason="PETSc not available")
pytest.importorskip("_aeroelast", reason="Rust backend not available")

from tests.validation.benchmarks.analytic_beam_helpers import (  # noqa: E402
    E, HEIGHT, LENGTH, RHO, THICKNESS, WIDTH,
    assembled_mass, build_box, composite_reference, modal_families,
    parked_reference, rotating_beam_ritz, section, static_rigidities, timoshenko_reference,
)

from tests.support.assertions import assert_residual_below  # noqa: E402

MESHES = [(40, 2), (80, 4), (160, 8)]
STATIC_MESHES = [(160, 4), (320, 4), (640, 4)]
MODAL_MESHES = [(160, 8), (240, 12), (320, 16)]


def report(name, values, reference, meshes=MESHES):
    values, reference = np.array(values), np.asarray(reference)
    print(f"\n{name}: meshes={meshes}; values={values.tolist()}; ref={reference.tolist()}; "
          f"relative_error={(values/reference-1).tolist()}", flush=True)
    return abs(values[-1]/reference-1)


@pytest.fixture(scope="module")
def modal_study():
    start = time.monotonic()
    families, masses, rotating = [], [], []
    reference = parked_reference()
    # Omega/omega_0=0.4 for flap; edge is softer, giving ratio ~0.60.
    omega = .4 * 2*np.pi*reference["flap"][0]
    for span, wall in MODAL_MESHES:
        print(f"\nmesh {span}x{wall}: building", flush=True)
        case = build_box(span, wall, triangles=True)
        masses.append(assembled_mass(case))
        print("parked solve", flush=True)
        families.append(modal_families(case))
        print(f"parked={families[-1]}; rotating solve", flush=True)
        rotating.append(modal_families(case, omega))
        print(f"rotating={rotating[-1]}", flush=True)
    print(f"\nmodal study runtime={time.monotonic()-start:.3f}s", flush=True)
    return families, masses, rotating, omega


@pytest.mark.parametrize("caller_ordering", [None, "natural"])
def test_modal_ordering_option_is_restored(caller_ordering):
    """RCM is scoped: a caller's explicit ordering (or absence) survives."""
    from petsc4py import PETSc

    options = PETSc.Options()
    key = "st_pc_factor_mat_ordering_type"
    original = options.getString(key) if options.hasName(key) else None
    try:
        if caller_ordering is None:
            if options.hasName(key):
                del options[key]
        else:
            options[key] = caller_ordering
        modal_families(build_box(4, 1, triangles=True))
        observed = options.getString(key) if options.hasName(key) else None
        assert observed == caller_ordering
    finally:
        if original is None:
            if options.hasName(key):
                del options[key]
        else:
            options[key] = original


def test_independent_beam_references_converge():
    """Reference self-convergence precedes any shell comparison (no fitted K1)."""
    for family, values in parked_reference().items():
        if family != "torsion":
            np.testing.assert_allclose(timoshenko_reference()[family],
                                       timoshenko_reference(20)[family], rtol=1e-6, atol=0)
            assert np.all(timoshenko_reference()[family] < values)
    for ratio in [.4, .4*np.sqrt(section()[1]/section()[2])]:
        k, low = rotating_beam_ritz(ratio)
        refined, _ = rotating_beam_ritz(ratio, degree=16)
        assert abs(k-refined) < 1e-6
        assert abs(low-1.193) < .0005  # literature rounded to three decimals


def test_exact_closed_section_mass():
    """Rigid x/y/z translation each measures rho*perimeter*t*L, not trace(M)."""
    start = time.monotonic()
    masses = [assembled_mass(build_box(s, w, triangles=True)) for s, w in MESHES]
    print(f"mass study runtime={time.monotonic()-start:.3f}s")
    expected = RHO * 2*(WIDTH+HEIGHT)*THICKNESS*LENGTH
    errors = report("total mass kg (x,y,z)", masses, expected)
    assert_residual_below(np.max(errors), tol=1e-10,
                          reference_name="exact thin-wall mass rho * perimeter * t * L",
                          kind="analytical", what="finest-mesh rigid-translation mass")
    # Strict < replaces inclusive allclose; explicit atol=0 drops no absolute allowance.
    assert_residual_below(np.max(np.abs(np.asarray(masses) / expected - 1.0)), tol=1e-10,
                          reference_name="exact thin-wall mass rho * perimeter * t * L",
                          kind="analytical", what="all-mesh rigid-translation masses")


@pytest.mark.slow
@pytest.mark.parametrize("family,order", [("flap", 0), ("flap", 1),
                                         ("edge", 0), ("edge", 1), ("torsion", 0)])
def test_parked_closed_section_modes(modal_study, family, order):
    """EB/Bredt wave references; bending shear/rotary correction is reported."""
    parked, _, _, _ = modal_study
    values = [m[family][order] for m in parked]
    reference = parked_reference()[family][order]
    if family != "torsion":
        corrected = timoshenko_reference()[family][order]
        refined = timoshenko_reference(degree=20)[family][order]
        assert abs(corrected/refined-1) < 1e-6
        print(f"EB={reference:.9f}; Timoshenko={corrected:.9f}; "
              f"shear+rotary shift={corrected/reference-1:.6%}")
        reference = corrected
    error = report(f"parked {family} {order+1} Hz", values, reference, MODAL_MESHES)
    assert abs(values[-1]/values[-2]-1) < .005, "mesh not converged within 0.5%"
    assert_residual_below(error, tol=.02,
                          reference_name="cantilever Timoshenko shear/rotary bending or Saint-Venant-Bredt torsional wave",
                          kind="analytical", what=f"parked {family} mode {order+1}")


@pytest.mark.slow
@pytest.mark.parametrize("family,subtract", [("flap", 0), ("edge", 1)])
def test_rotating_southwell_same_s4_path(modal_study, family, subtract):
    """Measure K, not just frequency: a missing K_SP or K_G cannot hide."""
    parked, _, rotating, omega = modal_study
    reference_parked = parked_reference()[family][0]
    ratio = omega/(2*np.pi*reference_parked)
    k1, low_speed = rotating_beam_ritz(ratio)
    k1_finer, _ = rotating_beam_ritz(ratio, degree=16)
    assert abs(k1-k1_finer) < 1e-6
    # Literature cross-check is a rounded range, not a fitted shell target.
    assert 1.17 < low_speed < 1.20
    values = [((2*np.pi*r[family][0])**2-(2*np.pi*p[family][0])**2)/omega**2
              for p, r in zip(parked, rotating, strict=True)]
    error = report(f"Southwell {family}, ratio={ratio:.5f}, K1(0)={low_speed:.8f}",
                   values, k1-subtract, MODAL_MESHES)
    reference_rotating = np.sqrt((2*np.pi*reference_parked)**2 +
                                 (k1-subtract)*omega**2)/(2*np.pi)
    report(f"rotating {family} 1 Hz", [r[family][0] for r in rotating],
           reference_rotating, MODAL_MESHES)
    assert abs(values[-1]-values[-2]) < .005, "coefficient not mesh converged"
    assert_residual_below(error, tol=.03,
                          reference_name="rotating Euler-Bernoulli ODE finite-speed Southwell coefficient; Wright et al. 1982 low-speed check",
                          kind="analytical", what=f"rotating {family} Southwell coefficient")


@pytest.mark.parametrize("surface", ["mid", "outer"])
def test_explicit_reference_surface_inertia(surface):
    """Same t=.05 walls: outer placement shifts planes outward by t/2.

    The mid dimensions are (.6,1); OML planes are (.65,1.05). Each reference
    integrates its own geometry, including local t^3/12. This is NOT the solid
    rectangle-minus-hole formula with overlapping corner material.
    """
    start = time.monotonic()
    shift = THICKNESS if surface == "outer" else 0
    width, height = WIDTH+shift, HEIGHT+shift
    expected = E * section(width, height)[1]
    values = [static_rigidities(build_box(s, w, width, height))[1] for s, w in STATIC_MESHES]
    error = report(f"{surface} EI N m2 runtime={time.monotonic()-start:.3f}s", values, expected,
                   STATIC_MESHES)
    assert abs(values[-1]/values[-2]-1) < .005
    assert_residual_below(error, tol=.005,
                          reference_name="Euler-Bernoulli EI from four wall-strip inertias including local t^3/12",
                          kind="analytical", what=f"{surface} surface EI")


@pytest.mark.parametrize("quantity,index,tolerance", [("EA", 0, 1e-5), ("EI", 1, .005)])
def test_composite_closed_box_effective_clt(quantity, index, tolerance):
    """Symmetric balanced glass [0/+45/-45/90]s, free transverse contraction."""
    start = time.monotonic()
    values = [static_rigidities(build_box(s, w, thickness=.01, composite=True))[index]
              for s, w in STATIC_MESHES]
    expected = composite_reference()[index]
    error = report(f"glass {quantity} runtime={time.monotonic()-start:.3f}s", values, expected,
                   STATIC_MESHES)
    assert abs(values[-1]/values[-2]-1) < .005
    assert_residual_below(error, tol=tolerance,
                          reference_name="Jones hand CLT A_eff = 1/a11, integrated over closed-box walls",
                          kind="analytical", what=f"glass closed-box {quantity}")
