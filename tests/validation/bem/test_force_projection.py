"""Tests for conservative force projection (BEM → shell mesh).

No external BEM dependencies required — uses synthetic BEM results and
a simple rectangular mesh to verify:
- Total force conservation: Σ f_nodes == F_BEM_integrated
- Correct handling of single-node strips (M_strip is dropped; a WARNING
  is emitted but force equilibrium is still preserved)
"""

import numpy as np
import pytest

from aeroelast.core.mesh.entities import ElementType, MeshElement, Node
from aeroelast.core.mesh.model import MeshModel
from aeroelast.models.blade.aerodynamics import AeroStation, AirfoilAero, BladeAero, PolarData
from aeroelast.solvers.bem.engine import BEMResult
from aeroelast.solvers.bem.force_projection import ForceProjector

# =====================================================================
# Synthetic-plate geometry (every expectation below is derived from this,
# never from ForceProjector's own per-strip frame):
#
#   * the strip outline extends along **global x**: the chordwise node index
#     i runs x = 0 .. chord_length, and the plate is flat (y = 0);
#   * the strips are stacked along **global z**: the spanwise node index j
#     runs z = hub_radius .. hub_radius + span_length, so span_dir = z_hat;
#   * the section chord therefore runs along x_hat and the section normal is
#     ``chord_hat x span_dir = x_hat x z_hat = -y_hat``.  Under the corrected
#     physics *Np rides y* (negative sense) and *Tp rides x* - the old
#     expectation "Np on global x, Tp on global y" encoded the defect;
#   * the integrated magnitude of a uniform load is ``sum_k Np_k * dr_k``,
#     with the strip widths dr_k rebuilt here from the BEM station grid and
#     the mesh span datum (:func:`_strip_widths`), so the expected value is
#     never taken from ``ForceProjector.verify()`` - whose ``force_bem`` is
#     assembled from the projector's own frames, the suite-audit circularity
#     (verdict section 22.1).
# =====================================================================


# =====================================================================
# Helpers
# =====================================================================


def _make_rectangular_mesh(n_span: int, n_chord: int, span_length: float, chord_length: float):
    """Create a flat rectangular quad mesh in the XZ plane.

    Span direction: Z (from hub_radius to hub_radius + span_length).
    Chord direction: X.
    Y = 0 (flat plate).

    Returns (mesh, hub_radius).
    """
    hub_radius = 3.0  # arbitrary (the BEM stations' rotor-centre origin)
    Node._id_counter = 0
    MeshElement._id_counter = 0

    nodes = []
    for j in range(n_span):
        # blade-local span (root at 0) — matches the ForceProjector's
        # convention (no hub offset in the mesh coordinates)
        z = j * span_length / (n_span - 1)
        for i in range(n_chord):
            x = i * chord_length / (n_chord - 1)
            nodes.append(Node([x, 0.0, z]))

    elements = []
    for j in range(n_span - 1):
        for i in range(n_chord - 1):
            n0 = j * n_chord + i
            n1 = n0 + 1
            n2 = n0 + n_chord + 1
            n3 = n0 + n_chord
            elem_nodes = [nodes[n0], nodes[n1], nodes[n2], nodes[n3]]
            elements.append(MeshElement(elem_nodes, ElementType.quad))

    mesh = MeshModel(nodes=nodes, elements=elements)
    return mesh, hub_radius


def _strip_widths(blade_aero: BladeAero, mesh: MeshModel, span_dir) -> np.ndarray:
    """Strip widths ``dr_k`` from the BEM station grid and the mesh span datum.

    The widths are the *discretisation* of the BEM stations onto the mesh
    span (station midpoints, with the first and last strips extended half a
    spacing to the mesh ends).  They are geometry, independent of the load
    frame: reading them back from ``ForceProjector.verify()`` would rebuild
    the expectation from the code under test (verdict section 22.1).
    """
    span_dir = np.asarray(span_dir, dtype=float)
    span_dir = span_dir / np.linalg.norm(span_dir)
    r = np.asarray(blade_aero.r, dtype=float)
    span_coords = mesh.coords_array @ span_dir
    r = r + (span_coords.min() - r[0])
    r_mid = 0.5 * (r[:-1] + r[1:])
    edges = np.concatenate(
        [[2.0 * r[0] - r_mid[0]], r_mid, [2.0 * r[-1] - r_mid[-1]]]
    )
    return np.diff(edges)


def _make_simple_blade_aero(n_stations: int, hub_radius: float, span_length: float):
    """Create a BladeAero with uniform stations along the span."""
    alpha = np.linspace(-np.pi, np.pi, 37)
    polar = PolarData(
        alpha=alpha,
        cl=2 * np.pi * np.sin(alpha),
        cd=np.full(37, 0.01),
        cm=np.zeros(37),
        re=1e6,
    )
    airfoil = AirfoilAero(
        name="flat_plate",
        coordinates=np.array([[1, 0], [0, 0]]),
        relative_thickness=0.0,
        aerodynamic_center=0.25,
        polars=[polar],
    )

    stations = []
    for i in range(n_stations):
        eta = i / (n_stations - 1)
        r = hub_radius + eta * span_length
        stations.append(
            AeroStation(
                span_fraction=eta,
                r=r,
                chord=1.0,
                twist=0.0,
                pitch_axis=0.25,
                airfoil=airfoil,
            )
        )

    return BladeAero(
        airfoils=[airfoil],
        stations=stations,
        blade_length=span_length,
        hub_radius=hub_radius,
        rotor_radius=hub_radius + span_length,
        n_blades=3,
    )


def _make_uniform_bem_result(blade_aero: BladeAero, Np_val: float, Tp_val: float):
    """Create a BEMResult with uniform Np and Tp along the span."""
    n = len(blade_aero.stations)
    return BEMResult(
        r=blade_aero.r,
        Np=np.full(n, Np_val),
        Tp=np.full(n, Tp_val),
        alpha=np.zeros(n),
        cl=np.zeros(n),
        cd=np.zeros(n),
        a=np.zeros(n),
        ap=np.zeros(n),
        thrust=0.0,  # not used by projector
        torque=0.0,
        power=0.0,
    )


# =====================================================================
# Tests
# =====================================================================


class TestForceProjectorConstruction:
    """Basic construction tests."""

    def test_creates_strips(self):
        mesh, hub_r = _make_rectangular_mesh(
            n_span=10, n_chord=5, span_length=20.0, chord_length=1.0
        )
        blade_aero = _make_simple_blade_aero(n_stations=10, hub_radius=hub_r, span_length=20.0)
        projector = ForceProjector(mesh, blade_aero)
        # Should have as many strips as BEM stations
        assert len(projector._strips) == 10

    def test_all_nodes_assigned(self):
        """Every mesh node should belong to some strip."""
        mesh, hub_r = _make_rectangular_mesh(
            n_span=10, n_chord=5, span_length=20.0, chord_length=1.0
        )
        blade_aero = _make_simple_blade_aero(n_stations=10, hub_radius=hub_r, span_length=20.0)
        projector = ForceProjector(mesh, blade_aero)

        assigned = set()
        for strip in projector._strips:
            assigned.update(strip.node_indices.tolist())
        assert len(assigned) == len(mesh.nodes)


class TestForceConservation:
    """Force conservation: sum of projected nodal forces == BEM integrated force."""

    @pytest.fixture
    def setup(self):
        n_span, n_chord = 11, 5
        span_length = 20.0
        mesh, hub_r = _make_rectangular_mesh(n_span, n_chord, span_length, chord_length=1.0)
        blade_aero = _make_simple_blade_aero(
            n_stations=n_span, hub_radius=hub_r, span_length=span_length
        )
        return mesh, blade_aero

    def test_uniform_Np_conservation(self, setup):
        """Uniform Np=1000 N/m: total normal force = Np * L, riding -y."""
        mesh, blade_aero = setup
        Np_val = 1000.0
        bem_result = _make_uniform_bem_result(blade_aero, Np_val=Np_val, Tp_val=0.0)

        projector = ForceProjector(mesh, blade_aero, span_direction=[0, 0, 1])
        forces = projector.project(bem_result)

        # Geometry: chord along x, span along z -> the raw section normal is
        # x_hat x z_hat = -y_hat, and the global sign rule then flips it to the
        # +Y half-space declared by normal_direction (the fluid / downwind
        # direction), so Np rides +y with |total| = sum_k Np_k dr_k (dr_k
        # derived in the test from the station grid, never read back from verify()).
        expected = Np_val * _strip_widths(blade_aero, mesh, [0, 0, 1]).sum()
        # Perpendicular components (x and z) must be ~0 ...
        np.testing.assert_allclose(forces[:, 0].sum(), 0.0, atol=1e-6)
        np.testing.assert_allclose(forces[:, 2].sum(), 0.0, atol=1e-6)
        # ... and the full normal force must ride +y.  This replaces the old
        # "y == 0" check, which encoded the global-axis defect.
        np.testing.assert_allclose(forces[:, 1].sum(), expected, rtol=1e-9)

        # The projector's own force bookkeeping is kept as a consistency
        # check; the physics expectation is the geometry-derived line above.
        verification = projector.verify(bem_result, forces)
        assert verification["force_error"] < 1e-6, (
            f"Force error = {verification['force_error']:.2e} N"
        )

    def test_uniform_Tp_conservation(self, setup):
        """Uniform Tp=500 N/m: total tangential force = Tp * L, riding +x."""
        mesh, blade_aero = setup
        Tp_val = 500.0
        bem_result = _make_uniform_bem_result(blade_aero, Np_val=0.0, Tp_val=Tp_val)

        projector = ForceProjector(mesh, blade_aero, span_direction=[0, 0, 1])
        forces = projector.project(bem_result)

        # Geometry: Tp rides the chord, which runs along +x here, with
        # |total| = sum_k Tp_k dr_k.
        expected = Tp_val * _strip_widths(blade_aero, mesh, [0, 0, 1]).sum()
        # Perpendicular components (y and z) must be ~0 ...
        np.testing.assert_allclose(forces[:, 1].sum(), 0.0, atol=1e-6)
        np.testing.assert_allclose(forces[:, 2].sum(), 0.0, atol=1e-6)
        # ... and the full tangential force must ride +x.  This replaces the
        # old "x == 0" check, which encoded the global-axis defect.
        np.testing.assert_allclose(forces[:, 0].sum(), expected, rtol=1e-9)

        verification = projector.verify(bem_result, forces)
        assert verification["force_error"] < 1e-6

    def test_combined_Np_Tp_conservation(self, setup):
        """Combined Np + Tp should conserve both components."""
        mesh, blade_aero = setup
        bem_result = _make_uniform_bem_result(blade_aero, Np_val=800.0, Tp_val=300.0)

        projector = ForceProjector(mesh, blade_aero, span_direction=[0, 0, 1])
        forces = projector.project(bem_result)

        verification = projector.verify(bem_result, forces)
        assert verification["force_error"] < 1e-6

    def test_varying_Np_conservation(self, setup):
        """Linearly varying Np (root=2000, tip=0): force conservation via trapz."""
        mesh, blade_aero = setup
        n = len(blade_aero.stations)
        Np_varying = np.linspace(2000.0, 0.0, n)

        bem_result = BEMResult(
            r=blade_aero.r,
            Np=Np_varying,
            Tp=np.zeros(n),
            alpha=np.zeros(n),
            cl=np.zeros(n),
            cd=np.zeros(n),
            a=np.zeros(n),
            ap=np.zeros(n),
            thrust=0.0,
            torque=0.0,
            power=0.0,
        )

        projector = ForceProjector(mesh, blade_aero, span_direction=[0, 0, 1])
        forces = projector.project(bem_result)

        verification = projector.verify(bem_result, forces)
        # Slightly more relaxed for varying loads — strip discretisation
        # introduces small mismatch vs. continuous trapz integral
        assert verification["force_error"] < 50.0, (
            f"Force error = {verification['force_error']:.2e} N (expected < 50 N for a 20 m blade)"
        )


class TestForceProjectionOutput:
    """Output shape and direction tests."""

    def test_output_shape(self):
        mesh, hub_r = _make_rectangular_mesh(
            n_span=6, n_chord=4, span_length=10.0, chord_length=1.0
        )
        blade_aero = _make_simple_blade_aero(n_stations=6, hub_radius=hub_r, span_length=10.0)
        bem_result = _make_uniform_bem_result(blade_aero, Np_val=100.0, Tp_val=0.0)

        projector = ForceProjector(mesh, blade_aero)
        forces = projector.project(bem_result)

        assert forces.shape == (len(mesh.nodes), 3)

    def test_zero_load_gives_zero_forces(self):
        mesh, hub_r = _make_rectangular_mesh(
            n_span=6, n_chord=4, span_length=10.0, chord_length=1.0
        )
        blade_aero = _make_simple_blade_aero(n_stations=6, hub_radius=hub_r, span_length=10.0)
        bem_result = _make_uniform_bem_result(blade_aero, Np_val=0.0, Tp_val=0.0)

        projector = ForceProjector(mesh, blade_aero)
        forces = projector.project(bem_result)

        np.testing.assert_allclose(forces, 0.0, atol=1e-12)

    def test_forces_only_in_load_direction(self):
        """Np-only load must ride the section normal (+y on this flat plate)."""
        mesh, hub_r = _make_rectangular_mesh(
            n_span=6, n_chord=4, span_length=10.0, chord_length=1.0
        )
        blade_aero = _make_simple_blade_aero(n_stations=6, hub_radius=hub_r, span_length=10.0)
        bem_result = _make_uniform_bem_result(blade_aero, Np_val=500.0, Tp_val=0.0)

        # Production defaults.  This test used to pass normal_direction=[1,0,0] /
        # tangential_direction=[0,1,0] explicitly - the pre-P5 global pair, which is
        # (near-)orthogonal to the section normal, so its sign dot product is ~0 and the
        # resolved sense is round-off.  Pinning a sign through an ill-conditioned
        # reference is exactly the defect that made the aero load upwind, so the test
        # now uses the configured defaults and asserts the sign against geometry.
        projector = ForceProjector(mesh, blade_aero)
        forces = projector.project(bem_result)

        # Geometry: chord along x, span along z -> section normal
        # x_hat x z_hat = -y_hat, and the global sign rule resolves that pair to
        # the +Y half-space (the fluid direction), so an Np-only load rides +y,
        # not the configured global x (see the module-level geometry block).
        # Every nodal force must then be pure +y: the perpendicular x and z
        # components are ~0 and the y component carries the signed total
        # sum_k Np_k dr_k, so a dropped/mis-summed strip fails here too.
        expected_total = 500.0 * _strip_widths(blade_aero, mesh, [0, 0, 1]).sum()
        np.testing.assert_allclose(forces[:, 0], 0.0, atol=1e-8)
        np.testing.assert_allclose(forces[:, 2], 0.0, atol=1e-8)
        np.testing.assert_allclose(forces[:, 1].sum(), expected_total, rtol=1e-9)
        # The sign is asserted explicitly: the geometry-derived total is
        # positive (thrust along +Y = downwind), so a sign inversion fails here.
        assert forces[:, 1].sum() > 0, (
            f"Projected normal force has wrong sign: Σfy={forces[:, 1].sum():.3e}"
        )


class TestSingleNodeStrip:
    """Edge case: strips with a single node should receive the full strip force."""

    def test_single_node_per_strip(self):
        """Mesh with 1 chordwise node per span station."""
        Node._id_counter = 0
        MeshElement._id_counter = 0

        hub_r = 3.0
        span_length = 10.0
        n_stations = 5

        # Create nodes along span, one per strip (blade-local, root at 0)
        nodes = []
        for j in range(n_stations):
            z = j * span_length / (n_stations - 1)
            nodes.append(Node([0.5, 0.0, z]))

        # Minimal triangle elements (not used by projector, but mesh needs them)
        elements = []
        for j in range(n_stations - 1):
            # Degenerate triangles just for mesh validity
            elem_nodes = [nodes[j], nodes[j + 1], nodes[j]]
            elements.append(MeshElement(elem_nodes, ElementType.triangle))

        mesh = MeshModel(nodes=nodes, elements=elements)
        blade_aero = _make_simple_blade_aero(
            n_stations=n_stations, hub_radius=hub_r, span_length=span_length
        )
        bem_result = _make_uniform_bem_result(blade_aero, Np_val=1000.0, Tp_val=0.0)

        projector = ForceProjector(mesh, blade_aero, span_direction=[0, 0, 1])
        forces = projector.project(bem_result)

        # Each node should get the full strip force (F = Np * dr)
        assert forces.shape == (n_stations, 3)
        # A single-node strip has no chordwise extent, so the section frame
        # falls back to the configured sense vectors: normal = +Y (the fluid
        # direction), chord = +X.  The normal force per node is therefore
        # +Np * dr on y and ~0 elsewhere.  The sign is asserted, not
        # |force| > 0, so a sign error fails here.
        dr = span_length / (n_stations - 1)
        assert np.all(forces[:, 1] > 0), f"normal force must be positive, got {forces[:, 1]}"
        assert np.allclose(forces[:, 1], 1000.0 * dr, rtol=1e-6)
        assert np.allclose(forces[:, [0, 2]], 0.0, atol=1e-9)
        # Total force conservation
        verification = projector.verify(bem_result, forces)
        assert verification["force_error"] < 1.0  # relaxed for coarse discretisation
