"""The thin-walled-tube twist measured through the production ``ForceProjector.project()`` (P1a).

The tube validation in ``test_thin_walled_tube_torsion.py`` realised its load with a construction
written *inside that test* (``_ring_traction`` plus ``_distributed_pattern_load``), so it validated
**that pattern**, not the load path the blade campaigns actually use. This module measures the same
quantity through the production entry point: the distributed station load is handed to
``aeroelast.solvers.bem.force_projection.ForceProjector`` and ``.project()`` - the code the rotor
and blade paths call - and the resulting nodal forces are solved with the same clamped-free MITC4
tube and the same converged metric (``theta_z``) as the reference module. Everything else is reused
from ``tests.validation.parity.test_thin_walled_tube_torsion`` (mesh, properties, Bredt ``GJ``, the
solvers and the metrics), so the only new object under test is the projection.

The adapter is honest about its **artificial inputs**. ``ForceProjector`` expects a ``BladeAero``
(chords, aerofoils, a hub radius); a rectangular tube is dressed as one:

* the chord is the tube's own width ``B`` (with ``aerodynamic_center = 0.5``), so the section frame
  the projector derives is the tube's real section, not an invented aerofoil;
* each station carries a **flat fake polar** - an ``AirfoilAero`` whose ``polars`` list is empty -
  because only the station geometry (radii, chord, aerodynamic centre) reaches the projection; the
  BEM result is supplied directly and no polar is ever evaluated.

Finding (measured before the bound was chosen, on the diagnostic ``probe.py``; see
``odd/tasks/production-path-and-independent-arbiters.md``). Against the continuum closed form
``theta(L) = e f L^2 / (2 GJ)`` the production path is ~3.11% high (10 mm wall) / ~3.24% high
(0.5 mm), and the excess is **entirely the projector's band quadrature**: ``ForceProjector``
extends the first and last strip bands half a spacing past the stations, so the outer bands carry
``DZ`` each and ``sum_k dr_k = L + DZ = 6.2 m`` here, and the projected load's own torque ruler
reads exactly ``e f (L + DZ)`` rather than ``e f L``. Against the exact **discrete** Saint-Venant
response to the projector's own station load - the couples ``e f dr_k`` at ``z_k``, whose tip
rotation is ``(e f / GJ) sum_k dr_k (L - z_k)`` - the same path is only 0.22% low (10 mm) / 0.09%
low (0.5 mm). The continuum comparison carries that quadrature caveat; the discrete one is the sharp
attribution.
"""

from __future__ import annotations

import numpy as np
import pytest

pytest.importorskip("_aeroelast", reason="Rust backend not available")

from aeroelast.core.mesh.entities import ElementSet, ElementType, MeshElement, Node  # noqa: E402
from aeroelast.core.mesh.model import MeshModel  # noqa: E402
from aeroelast.models.blade.aerodynamics import AeroStation, AirfoilAero, BladeAero  # noqa: E402
from aeroelast.solvers.bem.engine import BEMResult  # noqa: E402
from aeroelast.solvers.bem.force_projection import ForceProjector  # noqa: E402

import tests.validation.parity.test_thin_walled_tube_torsion as tube  # noqa: E402
from tests.support.assertions import assert_relative_error  # noqa: E402

# The suite's 5% rule, written as a module-level literal because the store's extractor resolves a
# bound declared in this module and not an alias imported from another one. Same value as
# tests/validation/parity/test_thin_walled_tube_torsion.py::TOL.
TOL = 0.05
FORCE_PER_LENGTH = 1.0e5  # [N/m] uniform transverse load along the span
OFFSET = 0.2  # [m] where the load's line of action sits, from the section centre


def _tube_mesh_model(coords: np.ndarray, conn: list) -> MeshModel:
    """The tube as a ``MeshModel``, with every quad in one ``shell`` element set.

    ``ForceProjector`` reads the element set's property map (``G t`` per wall) to split a strip's
    torsion across its rings, so the set name has to be carried explicitly.
    """
    nodes = [Node(c) for c in coords]
    elems = [MeshElement([nodes[int(i)] for i in row], ElementType.quad) for row in conn]
    mesh = MeshModel(nodes=nodes, elements=elems)
    mesh.add_element_set(ElementSet(name="shell", elements=set(elems)))
    return mesh


def _tube_blade_aero(coords: np.ndarray, rings: list) -> BladeAero:
    """The tube dressed as a ``BladeAero``: one station per ring, chord ``B``, a flat fake polar.

    The artificial inputs are these and only these: the section is given the tube's own width ``B``
    as its chord and a two-point flat "aerofoil" with an empty polar list. They are never used as
    aerodynamic data - ``BEMResult`` is supplied directly and no polar is evaluated - only as the
    station frame the projector derives its chord axis and aerodynamic centre from.
    """
    airfoil = AirfoilAero(
        name="tube-section",
        coordinates=np.array([[0.0, 0.0], [tube.B, 0.0]]),
        relative_thickness=tube.H / tube.B,
        aerodynamic_center=0.5,
        polars=[],
    )
    z = np.asarray([coords[r[0], 2] for r in rings])
    stations = [
        AeroStation(
            span_fraction=float(zi / z[-1]),
            r=float(zi),
            chord=float(tube.B),
            twist=0.0,
            pitch_axis=0.5,
            airfoil=airfoil,
        )
        for zi in z
    ]
    return BladeAero(
        airfoils=[airfoil],
        stations=stations,
        blade_length=float(z[-1]),
        hub_radius=0.0,
        rotor_radius=float(z[-1]),
        n_blades=1,
    )


def _tube_bem_result(
    n_stations: int, force_per_length: float, offset: float, r: np.ndarray
) -> BEMResult:
    """A station load with only ``Np`` and ``Mp`` non-zero: normal force plus a transfer moment."""
    zeros = np.zeros(n_stations)
    return BEMResult(
        r=np.asarray(r, dtype=float),
        Np=np.full(n_stations, force_per_length),
        Tp=zeros.copy(),
        alpha=zeros.copy(),
        cl=zeros.copy(),
        cd=zeros.copy(),
        a=zeros.copy(),
        ap=zeros.copy(),
        thrust=0.0,
        torque=0.0,
        power=0.0,
        Mp=np.full(n_stations, offset * force_per_length),
    )


def _project_tube(
    coords: np.ndarray, conn: list, rings: list, thickness: float, offset: float
) -> tuple[np.ndarray, dict]:
    """Project the tube's station load through production and embed it in the 6-DOF layout.

    Returns the assembled nodal force vector and the projector's own ``verify`` report.
    """
    mesh = _tube_mesh_model(coords, conn)
    blade_aero = _tube_blade_aero(coords, rings)
    projector = ForceProjector(
        mesh,
        blade_aero,
        hub_radius=0.0,
        element_properties={"shell": tube._iso_prop(thickness)},
    )
    z = np.asarray([coords[r[0], 2] for r in rings])
    bem = _tube_bem_result(len(rings), FORCE_PER_LENGTH, offset, z)
    forces = projector.project(bem)
    verify = projector.verify(bem, forces)

    f = np.zeros(6 * len(coords))
    f[0::6] = forces[:, 0]
    f[1::6] = forces[:, 1]
    f[2::6] = forces[:, 2]
    return f, verify


def test_production_projector_reproduces_the_exact_twist():
    """Clamped-free closed tube under the production station load, against exact comparators.

    The load is the one ``ForceProjector.project()`` produces: ``Np = 1e5 N/m``, ``Tp = 0`` and a
    transfer moment ``Mp = 0.2 * 1e5 N.m/m`` per station. The invariant block pins the load's own
    identity (force conservation, the torque ruler, the twist-torque sign pair); the four asserted
    comparisons then read ``theta_z`` against the continuum closed form and against the discrete
    Saint-Venant response to the projector's own station load. The control runs the same projector
    with the transfer moment removed: a load through the section centre must not twist the tube.
    """
    results: dict[float, dict[str, float]] = {}
    print("\n========== production ForceProjector on the clamped-free closed tube ==========")
    for thickness in (tube.THICKNESS, tube.THIN_WALL):
        coords, conn, rings, _ = tube._tube_mesh()
        _, K = tube._assemble(coords, conn, 4, tube._iso_prop(thickness))
        clamped = tube._clamped_dofs(rings)

        f, verify = _project_tube(coords, conn, rings, thickness, OFFSET)
        theta_z = tube._theta_z(tube._solve_clamped(K, f, clamped), rings[-1])
        torque = tube._realised_torque(coords, f)

        gj, _ = tube._bredt_isotropic(thickness)
        theta_closed = OFFSET * FORCE_PER_LENGTH * tube.L**2 / (2.0 * gj)
        z_k = np.asarray([coords[r[0], 2] for r in rings])
        theta_discrete = OFFSET * FORCE_PER_LENGTH * tube.DZ * float(np.sum(tube.L - z_k)) / gj

        # Load identity: no independent reference, so none of these is a store comparison.
        assert verify["force_error"] < 1e-8, (
            f"the projection does not conserve the strip force: {verify['force_error']:.3e} N"
        )
        expected_torque = -(OFFSET * FORCE_PER_LENGTH * tube.DZ * (tube.N_Z + 1))
        assert abs(torque - expected_torque) / abs(expected_torque) < 1e-9, (
            f"the torque ruler reads {torque:.6e} N.m against {expected_torque:.6e}: the outer "
            f"strip bands do not carry DZ each, so sum_k dr_k is not L + DZ"
        )
        assert torque * theta_z > 0.0, (
            f"the twist must follow the applied torque: T = {torque:.6e}, theta_z = {theta_z:.6e}"
        )

        # Control: the same projector with the transfer moment removed (force through the centre).
        f_control, _ = _project_tube(coords, conn, rings, thickness, 0.0)
        theta_control = tube._theta_z(tube._solve_clamped(K, f_control, clamped), rings[-1])
        assert abs(theta_control) < 0.02 * theta_closed, (
            f"a load through the section centre twisted the tube by {theta_control:.3e} rad "
            f"against {theta_closed:.3e}"
        )

        # sum_k dr_k read from the projector's own output: |T| = offset * FORCE_PER_LENGTH * sum dr_k.
        sum_dr = abs(torque) / (OFFSET * FORCE_PER_LENGTH)
        print(
            f"  t={thickness * 1e3:5.1f} mm  theta_z={theta_z:+.6e} rad  "
            f"continuum ratio={theta_z / -theta_closed:.6f}  "
            f"discrete ratio={theta_z / -theta_discrete:.6f}  "
            f"sum dr_k={sum_dr:.4f} vs L + DZ={tube.L + tube.DZ:.4f}"
        )
        results[thickness] = {
            "theta_z": theta_z,
            "theta_closed": theta_closed,
            "theta_discrete": theta_discrete,
        }
    print("===============================================================================")

    # Four comparisons, one per wall per comparator, each at its own call site so the store can map
    # them. The closed forms are magnitudes for a *positive* transfer moment; this projection's
    # station load realises a torque of the opposite sense about +z (the ruler assertion and the
    # sign invariant above pin that), so each signed theta_z is compared with the closed form the
    # applied load carries and the sign is not thrown away with an abs(). The reference_name strings
    # are literals at the call sites on purpose: the store's extractor only reads an ast.Constant.
    assert_relative_error(
        results[tube.THICKNESS]["theta_z"],
        -results[tube.THICKNESS]["theta_closed"],
        tol=TOL,
        kind="analytical",
        reference_name=(
            "the closed form e f L^2 / (2 GJ) for a uniform transverse load at offset e on a "
            "clamped-free closed tube, GJ from Bredt"
        ),
        what="production tip theta_z, 10 mm wall, against the continuum closed form",
    )
    assert_relative_error(
        results[tube.THIN_WALL]["theta_z"],
        -results[tube.THIN_WALL]["theta_closed"],
        tol=TOL,
        kind="analytical",
        reference_name=(
            "the closed form e f L^2 / (2 GJ) for a uniform transverse load at offset e on a "
            "clamped-free closed tube, GJ from Bredt"
        ),
        what="production tip theta_z, 0.5 mm wall, against the continuum closed form",
    )
    assert_relative_error(
        results[tube.THICKNESS]["theta_z"],
        -results[tube.THICKNESS]["theta_discrete"],
        tol=TOL,
        kind="analytical",
        reference_name=(
            "the discrete Saint-Venant response (e f / GJ) sum_k dr_k (L - z_k) to the "
            "projector's station load, dr_k the band width"
        ),
        what="production tip theta_z, 10 mm wall, against the discrete Saint-Venant response",
    )
    assert_relative_error(
        results[tube.THIN_WALL]["theta_z"],
        -results[tube.THIN_WALL]["theta_discrete"],
        tol=TOL,
        kind="analytical",
        reference_name=(
            "the discrete Saint-Venant response (e f / GJ) sum_k dr_k (L - z_k) to the "
            "projector's station load, dr_k the band width"
        ),
        what="production tip theta_z, 0.5 mm wall, against the discrete Saint-Venant response",
    )
