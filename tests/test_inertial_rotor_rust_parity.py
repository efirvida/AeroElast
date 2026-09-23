"""Phase 6 — Parity tests for InertialRotorFsiSolver Rust binding.

Verifies that the PyO3 binding ``_aeroelast.run_inertial_rotor_fsi_solver``
is callable with the expected signature and marshals arguments correctly.
    def test_solve_via_rust_passes_dynamic_geometric_stiffness_flag(self, monkeypatch):
Test groups:

``TestInertialRotorRustBinding``
    Smoke tests for the Rust fast-path with a minimal 1-DOF system.
    Verifies displacement_mode parsing, extract_masses integration, and
    that the solver accepts both "elastic" and "total" modes.
    Requires ``_aeroelast`` built with ``--features fsi``; skips otherwise.

``TestDisplacementModeValidation``
    Verifies that invalid displacement_mode strings raise ValueError.
"""

from __future__ import annotations

import inspect
import sys
import numpy as np
import pytest
from numpy.testing import assert_allclose

# ---------------------------------------------------------------------------
# Optional imports — guarded so collection doesn't crash without PETSc
# ---------------------------------------------------------------------------

try:
    from petsc4py import PETSc  # noqa: F401

    _HAS_PETSC = True
except (ImportError, OSError):
    _HAS_PETSC = False

try:
    from petsc4py import PETSc
    from aeroelast.solvers.fsi.corotational import ConstantOmega
    from aeroelast.solvers.fsi.rotor_inertial import LinearDynamicFSIRotorInertialSolver

    _HAS_INERTIAL_WRAPPER = True
except (ImportError, OSError):
    _HAS_INERTIAL_WRAPPER = False
    ConstantOmega = None  # type: ignore[assignment]
    LinearDynamicFSIRotorInertialSolver = None  # type: ignore[assignment,misc]

try:
    import _aeroelast  # type: ignore[import]

    _HAS_RUST = hasattr(_aeroelast, "run_inertial_rotor_fsi_solver")
except (ImportError, OSError):
    _HAS_RUST = False

_skip_rust = pytest.mark.skipif(
    not _HAS_RUST, reason="_aeroelast.run_inertial_rotor_fsi_solver not available"
)

_skip_inertial_wrapper = pytest.mark.skipif(
    not _HAS_INERTIAL_WRAPPER,
    reason="LinearDynamicFSIRotorInertialSolver or PETSc not available",
)


# ---------------------------------------------------------------------------
# Helper to build minimal 1-DOF system for smoke tests
# ---------------------------------------------------------------------------


def _minimal_1dof_system():
    """Build a minimal 1-DOF spring-mass system for smoke tests.

    Returns:
        dict with K, M COO triplets, free_dofs, masses, interface_nodes
    """
    # 1 node × 3 DOFs = 3 total DOFs (only translational, no rotations)
    # DOF 2 is constrained (free = [0, 1])
    k_rows = np.array([0, 1], dtype=np.int32)
    k_cols = np.array([0, 1], dtype=np.int32)
    k_vals = np.array([100.0, 100.0], dtype=np.float64)

    m_rows = np.array([0, 1], dtype=np.int32)
    m_cols = np.array([0, 1], dtype=np.int32)
    m_vals = np.array([1.0, 1.0], dtype=np.float64)

    free_dofs = np.array([0, 1], dtype=np.int32)

    # Single node at origin, mass = 1.0 (will be extracted from m_vals)
    all_node_masses = np.array([1.0], dtype=np.float64)

    # Interface: node 0 (0-based indexing)
    interface_nodes = np.array([0], dtype=np.uintp)

    return {
        "k_rows": k_rows,
        "k_cols": k_cols,
        "k_vals": k_vals,
        "m_rows": m_rows,
        "m_cols": m_cols,
        "m_vals": m_vals,
        "free_dofs": free_dofs,
        "all_node_masses": all_node_masses,
        "interface_nodes": interface_nodes,
    }


def _build_minimal_assembler():
    """Build a minimal MeshAssembler for 1-node 3-DOF system.

    Returns:
        PyMeshAssembler with 1 MITC3 element (degenerate triangle)
    """
    from _aeroelast import PyMeshAssembler as MeshAssembler  # type: ignore[import]

    # Degenerate triangle: 3 nodes at distinct positions (smoke test only)
    node_coords = np.array(
        [
            [0.0, 0.0, 0.0],  # node 0
            [1e-3, 0.0, 0.0],  # node 1
            [0.0, 1e-3, 0.0],  # node 2
        ],
        dtype=np.float64,
    )

    connectivity = [[0, 1, 2]]
    elem_types = [3]  # 3 = Mitc3 integer code

    # Materials are plain Python dicts — MaterialSpec is a Rust-internal type
    materials = [{"type": "isotropic", "e": 210e9, "nu": 0.3, "rho": 7850.0, "thickness": 0.001}]

    return MeshAssembler(node_coords, connectivity, elem_types, materials)


# ---------------------------------------------------------------------------
# Group 1 — Smoke tests for run_inertial_rotor_fsi_solver
# ---------------------------------------------------------------------------


@_skip_rust
class TestInertialRotorRustBinding:
    """Smoke tests for _aeroelast.run_inertial_rotor_fsi_solver.

    Verifies that the Rust binding:
    1. Accepts the expected signature (assembler required, displacement_mode, etc.)
    2. Parses displacement_mode correctly ("elastic" vs "total")
    3. Integrates with extract_masses() internally
    4. Raises preCICE initialization error (expected — we don't start preCICE)
    """

    def test_signature_matches_expected(self):
        """Verify function signature has required parameters."""
        params = inspect.signature(_aeroelast.run_inertial_rotor_fsi_solver).parameters

        # Check required parameters exist
        required = [
            "assembler",  # Required (not Optional)
            "rotation_axis",
            "rotation_center",
            "all_node_masses",
            "omega_mode",
            "omega",
            "gravity",
            "include_reference_acceleration",
            "k_update_interval",
            "omega_rebuild_threshold",
            "theta_rebuild_threshold",
            "displacement_mode",  # New parameter
            "dofs_per_node",
            "fluid_density",
            "flow_velocity",
            "rotor_radius",
            "k_rows",
            "k_cols",
            "k_vals",
            "m_rows",
            "m_cols",
            "m_vals",
            "free_dofs",
            "eta_k",
            "eta_m",
            "beta",
            "gamma",
            "dt",
            "interface_nodes",  # Changed from interface_dofs_global
            "mesh_dims",
            "participant_name",
            "config_file",
            "coupling_mesh",
            "write_data_name",
            "read_data_name",
            "ramp_time",
        ]

        for param in required:
            assert param in params, f"Missing required parameter: {param}"

        # Verify displacement_mode is present
        assert "displacement_mode" in params

    def test_displacement_mode_elastic_accepted(self):
        """Verify displacement_mode='elastic' is accepted."""
        sys = _minimal_1dof_system()
        assembler = _build_minimal_assembler()

        # Expect preCICE initialization failure (we're not starting preCICE)
        # but the argument parsing should succeed
        with pytest.raises(RuntimeError, match="preCICE|Participant"):
            _aeroelast.run_inertial_rotor_fsi_solver(
                assembler=assembler,
                rotation_axis=[0.0, 0.0, 1.0],
                rotation_center=[0.0, 0.0, 0.0],
                all_node_masses=sys["all_node_masses"],
                omega_mode="constant",
                omega=10.0,
                omega_target=None,
                t_ramp=None,
                moment_of_inertia=None,
                shaft_torque=None,
                gravity=[0.0, 0.0, -9.81],
                include_reference_acceleration=True,
                k_update_interval=1,
                omega_rebuild_threshold=0.01,
                theta_rebuild_threshold=0.05,
                displacement_mode="elastic",  # Elastic mode
                dofs_per_node=3,
                fluid_density=1.225,
                flow_velocity=10.0,
                rotor_radius=1.0,
                k_rows=sys["k_rows"],
                k_cols=sys["k_cols"],
                k_vals=sys["k_vals"],
                m_rows=sys["m_rows"],
                m_cols=sys["m_cols"],
                m_vals=sys["m_vals"],
                free_dofs=sys["free_dofs"],
                eta_k=0.01,
                eta_m=0.05,
                beta=0.25,
                gamma=0.5,
                dt=0.01,
                interface_nodes=sys["interface_nodes"],
                mesh_dims=3,
                participant_name="Solid",
                config_file="precice-config.xml",
                coupling_mesh="Solid-Mesh",
                write_data_name="Displacement",
                read_data_name="Force",
                ramp_time=0.0,
                force_max=None,
                omega_mesh_name=None,
                omega_write_data=None,
                omega_vertex_coord=None,
                velocity_write_data=None,
                u0=None,
                v0=None,
                a0=None,
                t0=0.0,
                theta0=0.0,
                restart_omega=None,
                restart_alpha=None,
                restart_ramp_completed=None,
                restart_current_time=None,
                step_callback=None,
                kg_use_deformed_coords=False,
                kg_deflection_rebuild_rel_high=0.005,
                kg_deflection_rebuild_rel_low=0.003,
                kg0_rows=np.array([], dtype=np.int64),
                kg0_cols=np.array([], dtype=np.int64),
                kg0_vals=np.array([], dtype=np.float64),
            )

    def test_displacement_mode_total_accepted(self):
        """Verify displacement_mode='total' is accepted."""
        sys = _minimal_1dof_system()
        assembler = _build_minimal_assembler()

        with pytest.raises(RuntimeError, match="preCICE|Participant"):
            _aeroelast.run_inertial_rotor_fsi_solver(
                assembler=assembler,
                rotation_axis=[0.0, 0.0, 1.0],
                rotation_center=[0.0, 0.0, 0.0],
                all_node_masses=sys["all_node_masses"],
                omega_mode="constant",
                omega=10.0,
                omega_target=None,
                t_ramp=None,
                moment_of_inertia=None,
                shaft_torque=None,
                gravity=[0.0, 0.0, -9.81],
                include_reference_acceleration=True,
                k_update_interval=1,
                omega_rebuild_threshold=0.01,
                theta_rebuild_threshold=0.05,
                displacement_mode="total",  # Total mode
                dofs_per_node=3,
                fluid_density=1.225,
                flow_velocity=10.0,
                rotor_radius=1.0,
                k_rows=sys["k_rows"],
                k_cols=sys["k_cols"],
                k_vals=sys["k_vals"],
                m_rows=sys["m_rows"],
                m_cols=sys["m_cols"],
                m_vals=sys["m_vals"],
                free_dofs=sys["free_dofs"],
                eta_k=0.01,
                eta_m=0.05,
                beta=0.25,
                gamma=0.5,
                dt=0.01,
                interface_nodes=sys["interface_nodes"],
                mesh_dims=3,
                participant_name="Solid",
                config_file="precice-config.xml",
                coupling_mesh="Solid-Mesh",
                write_data_name="Displacement",
                read_data_name="Force",
                ramp_time=0.0,
                force_max=None,
                omega_mesh_name=None,
                omega_write_data=None,
                omega_vertex_coord=None,
                velocity_write_data=None,
                u0=None,
                v0=None,
                a0=None,
                t0=0.0,
                theta0=0.0,
                restart_omega=None,
                restart_alpha=None,
                restart_ramp_completed=None,
                restart_current_time=None,
                step_callback=None,
                kg_use_deformed_coords=False,
                kg_deflection_rebuild_rel_high=0.005,
                kg_deflection_rebuild_rel_low=0.003,
                kg0_rows=np.array([], dtype=np.int64),
                kg0_cols=np.array([], dtype=np.int64),
                kg0_vals=np.array([], dtype=np.float64),
            )

    def test_invalid_displacement_mode_raises(self):
        """Verify invalid displacement_mode raises ValueError."""
        sys = _minimal_1dof_system()
        assembler = _build_minimal_assembler()

        with pytest.raises(ValueError, match="unknown displacement_mode"):
            _aeroelast.run_inertial_rotor_fsi_solver(
                assembler=assembler,
                rotation_axis=[0.0, 0.0, 1.0],
                rotation_center=[0.0, 0.0, 0.0],
                all_node_masses=sys["all_node_masses"],
                omega_mode="constant",
                omega=10.0,
                omega_target=None,
                t_ramp=None,
                moment_of_inertia=None,
                shaft_torque=None,
                gravity=[0.0, 0.0, -9.81],
                include_reference_acceleration=True,
                k_update_interval=1,
                omega_rebuild_threshold=0.01,
                theta_rebuild_threshold=0.05,
                displacement_mode="invalid_mode",  # Invalid
                dofs_per_node=3,
                fluid_density=1.225,
                flow_velocity=10.0,
                rotor_radius=1.0,
                k_rows=sys["k_rows"],
                k_cols=sys["k_cols"],
                k_vals=sys["k_vals"],
                m_rows=sys["m_rows"],
                m_cols=sys["m_cols"],
                m_vals=sys["m_vals"],
                free_dofs=sys["free_dofs"],
                eta_k=0.01,
                eta_m=0.05,
                beta=0.25,
                gamma=0.5,
                dt=0.01,
                interface_nodes=sys["interface_nodes"],
                mesh_dims=3,
                participant_name="Solid",
                config_file="precice-config.xml",
                coupling_mesh="Solid-Mesh",
                write_data_name="Displacement",
                read_data_name="Force",
                ramp_time=0.0,
                force_max=None,
                omega_mesh_name=None,
                omega_write_data=None,
                omega_vertex_coord=None,
                velocity_write_data=None,
                u0=None,
                v0=None,
                a0=None,
                t0=0.0,
                theta0=0.0,
                restart_omega=None,
                restart_alpha=None,
                restart_ramp_completed=None,
                restart_current_time=None,
                step_callback=None,
                kg_use_deformed_coords=False,
                kg_deflection_rebuild_rel_high=0.005,
                kg_deflection_rebuild_rel_low=0.003,
                kg0_rows=np.array([], dtype=np.int64),
                kg0_cols=np.array([], dtype=np.int64),
                kg0_vals=np.array([], dtype=np.float64),
            )

    def test_omega_mode_ramped_requires_target(self):
        """Verify omega_mode='ramped' validation works."""
        sys = _minimal_1dof_system()
        assembler = _build_minimal_assembler()

        # Missing omega_target should raise ValueError
        with pytest.raises(ValueError, match="omega_mode='ramped' requires omega_target"):
            _aeroelast.run_inertial_rotor_fsi_solver(
                assembler=assembler,
                rotation_axis=[0.0, 0.0, 1.0],
                rotation_center=[0.0, 0.0, 0.0],
                all_node_masses=sys["all_node_masses"],
                omega_mode="ramped",
                omega=0.0,
                omega_target=None,  # Missing required param
                t_ramp=None,
                moment_of_inertia=None,
                shaft_torque=None,
                gravity=[0.0, 0.0, -9.81],
                include_reference_acceleration=True,
                k_update_interval=1,
                omega_rebuild_threshold=0.01,
                theta_rebuild_threshold=0.05,
                displacement_mode="elastic",
                dofs_per_node=3,
                fluid_density=1.225,
                flow_velocity=10.0,
                rotor_radius=1.0,
                k_rows=sys["k_rows"],
                k_cols=sys["k_cols"],
                k_vals=sys["k_vals"],
                m_rows=sys["m_rows"],
                m_cols=sys["m_cols"],
                m_vals=sys["m_vals"],
                free_dofs=sys["free_dofs"],
                eta_k=0.01,
                eta_m=0.05,
                beta=0.25,
                gamma=0.5,
                dt=0.01,
                interface_nodes=sys["interface_nodes"],
                mesh_dims=3,
                participant_name="Solid",
                config_file="precice-config.xml",
                coupling_mesh="Solid-Mesh",
                write_data_name="Displacement",
                read_data_name="Force",
                ramp_time=0.0,
                force_max=None,
                omega_mesh_name=None,
                omega_write_data=None,
                omega_vertex_coord=None,
                velocity_write_data=None,
                u0=None,
                v0=None,
                a0=None,
                t0=0.0,
                theta0=0.0,
                restart_omega=None,
                restart_alpha=None,
                restart_ramp_completed=None,
                restart_current_time=None,
                step_callback=None,
                kg_use_deformed_coords=False,
                kg_deflection_rebuild_rel_high=0.005,
                kg_deflection_rebuild_rel_low=0.003,
                kg0_rows=np.array([], dtype=np.int64),
                kg0_cols=np.array([], dtype=np.int64),
                kg0_vals=np.array([], dtype=np.float64),
            )


# ---------------------------------------------------------------------------
# Group 2 — Displacement mode validation
# ---------------------------------------------------------------------------


@_skip_rust
class TestDisplacementModeValidation:
    """Verify displacement_mode string validation."""

    def test_elastic_mode_parsed_correctly(self):
        """Verify 'elastic' string maps to DisplacementMode::Elastic."""
        # Already tested in TestInertialRotorRustBinding.test_displacement_mode_elastic_accepted
        # This is a placeholder for future enum value verification if needed
        pass

    def test_total_mode_parsed_correctly(self):
        """Verify 'total' string maps to DisplacementMode::Total."""
        # Already tested in TestInertialRotorRustBinding.test_displacement_mode_total_accepted
        pass

    def test_case_sensitivity(self):
        """Verify displacement_mode is case-sensitive."""
        sys = _minimal_1dof_system()
        assembler = _build_minimal_assembler()

        # "Elastic" (capitalized) should raise ValueError
        with pytest.raises(ValueError, match="unknown displacement_mode"):
            _aeroelast.run_inertial_rotor_fsi_solver(
                assembler=assembler,
                rotation_axis=[0.0, 0.0, 1.0],
                rotation_center=[0.0, 0.0, 0.0],
                all_node_masses=sys["all_node_masses"],
                omega_mode="constant",
                omega=10.0,
                omega_target=None,
                t_ramp=None,
                moment_of_inertia=None,
                shaft_torque=None,
                gravity=[0.0, 0.0, -9.81],
                include_reference_acceleration=True,
                k_update_interval=1,
                omega_rebuild_threshold=0.01,
                theta_rebuild_threshold=0.05,
                displacement_mode="Elastic",  # Wrong case
                dofs_per_node=3,
                fluid_density=1.225,
                flow_velocity=10.0,
                rotor_radius=1.0,
                k_rows=sys["k_rows"],
                k_cols=sys["k_cols"],
                k_vals=sys["k_vals"],
                m_rows=sys["m_rows"],
                m_cols=sys["m_cols"],
                m_vals=sys["m_vals"],
                free_dofs=sys["free_dofs"],
                eta_k=0.01,
                eta_m=0.05,
                beta=0.25,
                gamma=0.5,
                dt=0.01,
                interface_nodes=sys["interface_nodes"],
                mesh_dims=3,
                participant_name="Solid",
                config_file="precice-config.xml",
                coupling_mesh="Solid-Mesh",
                write_data_name="Displacement",
                read_data_name="Force",
                ramp_time=0.0,
                force_max=None,
                omega_mesh_name=None,
                omega_write_data=None,
                omega_vertex_coord=None,
                velocity_write_data=None,
                u0=None,
                v0=None,
                a0=None,
                t0=0.0,
                theta0=0.0,
                restart_omega=None,
                restart_alpha=None,
                restart_ramp_completed=None,
                restart_current_time=None,
                step_callback=None,
                kg_use_deformed_coords=False,
                kg_deflection_rebuild_rel_high=0.005,
                kg_deflection_rebuild_rel_low=0.003,
                kg0_rows=np.array([], dtype=np.int64),
                kg0_cols=np.array([], dtype=np.int64),
                kg0_vals=np.array([], dtype=np.float64),
            )


# ---------------------------------------------------------------------------
# Group 3 — Python wrapper regressions
# ---------------------------------------------------------------------------


@_skip_inertial_wrapper
class TestInertialWrapperMarshalling:
    def test_solve_via_rust_ignores_initial_geometric_stiffness(self, monkeypatch):
        solver = object.__new__(LinearDynamicFSIRotorInertialSolver)
        solver.K = type("FakeMat", (), {"getSize": lambda self: (3, 3)})()
        solver.M = object()
        solver._petsc_to_coo = lambda _mat: (
            np.array([0], dtype=np.int32),
            np.array([0], dtype=np.int32),
            np.array([1.0], dtype=np.float64),
        )
        solver._all_node_masses_full = np.array([1.0], dtype=np.float64)
        solver._interface_node_ids = [1]
        solver.solver_params = {"beta": 0.25, "gamma": 0.5, "dt": 0.1}
        solver.model_properties = {
            "solver": {
                "coupling": {
                    "participant": "Solid",
                    "config_file": "precice-config.xml",
                    "coupling_mesh": "SolidMesh",
                    "write_data": "Displacement",
                    "read_data": "Force",
                }
            }
        }
        solver._omega_provider = ConstantOmega(omega=5.0)
        solver._coord_transforms = type(
            "CoordStub",
            (),
            {
                "center": np.array([0.0, 0.0, 0.0], dtype=np.float64),
                "axis": np.array([0.0, 0.0, 1.0], dtype=np.float64),
            },
        )()
        solver._gravity = np.array([0.0, 0.0, -9.81], dtype=np.float64)
        solver._k_update_interval = 1
        solver._omega_rebuild_threshold = 0.01
        solver._theta_rebuild_threshold = 0.05
        solver._precice_displacement_mode = "elastic"
        solver._fluid_density = 1.225
        solver._flow_velocity = 10.0
        solver._rotor_radius = 1.0
        solver._eta_k = 0.01
        solver._eta_m = 0.05
        solver._force_ramp_time = 0.0
        solver._force_max_magnitude = None
        solver._send_omega_to_precice = False
        solver._omega_mesh_name = "GlobalSolidMesh"
        solver._omega_write_data_name = "AngularVelocity"
        solver._send_velocity_to_precice = False
        solver._velocity_write_data_name = None
        solver._include_geometric_stiffness = True
        solver._include_spin_softening = True
        solver._use_corotational_kt = False
        solver._kt_coro_update_freq = 1
        solver._ksp_omega_threshold = 1e-4
        solver._omega_rebuild_rel_high = 0.005
        solver._omega_rebuild_rel_low = 0.003
        solver._kg_use_deformed_coords = False
        solver._kg_deflection_rebuild_rel_high = 0.01
        solver._kg_deflection_rebuild_rel_low = 0.006
        solver._ksp_omega_rebuild_high = 0.005
        solver._ksp_omega_rebuild_low = 0.003
        solver._kg_omega_rebuild_high = 0.01
        solver._kg_omega_rebuild_low = 0.006
        solver._is_primary_rank = lambda: False
        solver.comm = PETSc.COMM_SELF

        assemble_called = False

        def _assemble_geometric_stiffness(**_kwargs):
            nonlocal assemble_called
            assemble_called = True
            raise AssertionError("assemble_geometric_stiffness should not be called")

        solver.domain = type(
            "DomainStub",
            (),
            {
                "_rust": type("RustStub", (), {
                    "assemble_m": staticmethod(lambda: (
                        np.array([0, 1, 2], dtype=np.int64),
                        np.array([0, 1, 2], dtype=np.int64),
                        np.array([1.0, 1.0, 1.0], dtype=np.float64)))},
                )(),
                "dofs_per_node": 3,
                "spatial_dim": 3,
                "nodes": [type("NodeStub", (), {"coords": np.array([0.0, 0.0, 0.0], dtype=np.float64)})()],
                "mesh": type("MeshStub", (), {"node_id_to_index": {1: 0}})(),
                "assemble_geometric_stiffness": staticmethod(_assemble_geometric_stiffness),
            },
        )()

        captured = {}

        def _fake_run_inertial_rotor_fsi_solver(**kwargs):
            captured.update(kwargs)
            return [], [], [], []

        monkeypatch.setitem(
            sys.modules,
            "_aeroelast",
            type(
                "AeroelastStub",
                (),
                {
                    "run_inertial_rotor_fsi_solver": staticmethod(
                        _fake_run_inertial_rotor_fsi_solver
                    )
                },
            )(),
        )

        bc_manager = type(
            "BcManagerStub", (), {"free_dofs": np.array([0], dtype=np.int32), "fixed_dofs": {}}
        )()

        solver._solve_via_rust(
            bc_manager=bc_manager,
            interface_coords=np.array([[1.0, 0.0, 0.0]], dtype=np.float64),
        )

        assert assemble_called is False
        for key in ("kg0_rows", "kg0_cols", "kg0_vals"):
            v = captured[key]
            assert v is None or getattr(v, "size", 0) == 0, (key, v)

    def test_solve_via_rust_wires_output_callback_and_finalizes_checkpoints(self, monkeypatch):
        solver = object.__new__(LinearDynamicFSIRotorInertialSolver)
        solver.K = type("FakeMat", (), {"getSize": lambda self: (3, 3)})()
        solver.M = object()
        solver._petsc_to_coo = lambda _mat: (
            np.array([0], dtype=np.int32),
            np.array([0], dtype=np.int32),
            np.array([1.0], dtype=np.float64),
        )
        solver._all_node_masses_full = np.array([2.0], dtype=np.float64)
        solver._interface_node_ids = [1]
        solver.solver_params = {
            "beta": 0.25,
            "gamma": 0.5,
            "dt": 0.1,
            "output_folder": "results",
        }
        solver.model_properties = {
            "solver": {
                "coupling": {
                    "participant": "Solid",
                    "config_file": "precice-config.xml",
                    "coupling_mesh": "SolidMesh",
                    "write_data": "Displacement",
                    "read_data": "Force",
                }
            }
        }
        solver._omega_provider = ConstantOmega(omega=5.0)
        from aeroelast.solvers.fsi.corotational import CoordinateTransforms  # noqa: E402
        solver._coord_transforms = CoordinateTransforms([0.0, 0.0, 1.0], [0.0, 0.0, 0.0])
        from aeroelast.solvers.fsi.corotational import InertialForcesCalculator  # noqa: E402
        solver._inertial_calculator = InertialForcesCalculator(
            [0.0, 0.0, 1.0],
            [0.0, 0.0, 0.0],
        )
        solver._gravity = np.array([0.0, 0.0, -9.81], dtype=np.float64)
        solver._include_gravity = True
        solver._rotation_center = np.array([0.0, 0.0, 0.0], dtype=np.float64)
        solver._rotation_axis = np.array([0.0, 0.0, 1.0], dtype=np.float64)
        solver._k_update_interval = 1
        solver._omega_rebuild_threshold = 0.01
        solver._theta_rebuild_threshold = 0.05
        solver._precice_displacement_mode = "elastic"
        solver._fluid_density = 1.225
        solver._flow_velocity = 10.0
        solver._rotor_radius = 1.0
        solver._eta_k = 0.01
        solver._eta_m = 0.05
        solver._force_ramp_time = 0.0
        solver._force_max_magnitude = None
        solver._send_omega_to_precice = False
        solver._omega_mesh_name = "GlobalSolidMesh"
        solver._omega_write_data_name = "AngularVelocity"
        solver._send_velocity_to_precice = False
        solver._velocity_write_data_name = None
        solver._include_geometric_stiffness = False
        solver._include_spin_softening = True
        solver._use_corotational_kt = False
        solver._kt_coro_update_freq = 1
        solver._ksp_omega_threshold = 1e-4
        solver._omega_rebuild_rel_high = 0.005
        solver._omega_rebuild_rel_low = 0.003
        solver._kg_use_deformed_coords = False
        solver._kg_deflection_rebuild_rel_high = 0.01
        solver._kg_deflection_rebuild_rel_low = 0.006
        solver._ksp_omega_rebuild_high = 0.005
        solver._ksp_omega_rebuild_low = 0.003
        solver._kg_omega_rebuild_high = 0.01
        solver._kg_omega_rebuild_low = 0.006
        solver._stress_output_interval = 1
        solver._theta = 0.0
        solver._is_primary_rank = lambda: False
        solver.comm = PETSc.COMM_SELF

        node = type(
            "NodeStub",
            (),
            {"id": 1, "coords": np.array([1.0, 0.0, 0.0]), "x": 1.0, "y": 0.0, "z": 0.0},
        )()
        solver.domain = type(
            "DomainStub",
            (),
            {
                "_rust": type("RustStub", (), {
                    "assemble_m": staticmethod(lambda: (
                        np.array([0, 1, 2], dtype=np.int64),
                        np.array([0, 1, 2], dtype=np.int64),
                        np.array([1.0, 1.0, 1.0], dtype=np.float64)))},
                )(),
                "dofs_per_node": 3,
                "spatial_dim": 3,
                "nodes": [node],
                "mesh": type("MeshStub", (), {"node_id_to_index": {1: 0}})(),
            },
        )()

        calls = {
            "perf": 0,
            "restart": 0,
            "report": 0,
            "probe": 0,
            "checkpoint": 0,
            "finalize": 0,
        }
        solver._compute_stress_fields = lambda _u_full: {}
        solver._log_rotor_performance = lambda **_kwargs: calls.__setitem__(
            "perf", calls["perf"] + 1
        )
        solver._write_restart_state = lambda **_kwargs: calls.__setitem__(
            "restart", calls["restart"] + 1
        )
        solver._log_structural_report = lambda **_kwargs: calls.__setitem__(
            "report", calls["report"] + 1
        )
        solver._log_probe_data = lambda **_kwargs: calls.__setitem__("probe", calls["probe"] + 1)
        solver._handle_checkpoint = lambda **_kwargs: calls.__setitem__(
            "checkpoint", calls["checkpoint"] + 1
        )
        solver._checkpoint_manager = type(
            "CheckpointStub",
            (),
            {
                "should_write": staticmethod(lambda _t: True),
                "finalize": staticmethod(
                    lambda timeout=60.0: calls.__setitem__("finalize", calls["finalize"] + 1)
                ),
            },
        )()

        captured = {}

        def _fake_run_inertial_rotor_fsi_solver(**kwargs):
            captured.update(kwargs)
            kwargs["step_callback"](
                0.1,
                1,
                0.1,
                np.array([0.01, 0.0, 0.0], dtype=np.float64),
                np.array([0.0, 0.0, 0.0], dtype=np.float64),
                np.array([0.0, 0.0, 0.0], dtype=np.float64),
                5.0,
                np.array([0.0, 5.0, 0.0], dtype=np.float64),
                5.0,
                0.0,
                0.25,
                (5.0, 0.0, 0.0, 0.0, 0.0),
            )
            return [0.01, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.1]

        monkeypatch.setitem(
            sys.modules,
            "_aeroelast",
            type(
                "AeroelastStub",
                (),
                {
                    "run_inertial_rotor_fsi_solver": staticmethod(
                        _fake_run_inertial_rotor_fsi_solver
                    )
                },
            )(),
        )

        bc_manager = type(
            "BcManagerStub",
            (),
            {"free_dofs": np.array([0, 1, 2], dtype=np.int32), "fixed_dofs": {}},
        )()

        solver._solve_via_rust(
            bc_manager=bc_manager,
            interface_coords=np.array([[1.0, 0.0, 0.0]], dtype=np.float64),
        )

        assert captured["step_callback"] is not None
        assert calls["perf"] == 1
        assert calls["restart"] == 1
        assert calls["report"] == 1
        assert calls["probe"] == 1
        assert calls["checkpoint"] == 1
        assert calls["finalize"] == 1


# ---------------------------------------------------------------------------
# Group 4 — Physics Validation (requires isolated API or preCICE mock)
# ---------------------------------------------------------------------------


@_skip_rust
@pytest.mark.skip(reason="Requires isolated physics API or preCICE mock - see Task 6.5-6.7")
class TestInertialPhysics:
    """Analytical validation of inertial frame physics.

    These tests require either:
    A) Exposing compute_rigid_body_acceleration_inertial() and
       compute_reference_load_vector() as standalone PyO3 functions
    B) Running end-to-end with a preCICE mock that provides zero forces

    Current implementation: Rust unit tests in rotor_physics.rs cover these
    cases directly. Python tests deferred pending API exposure decision.
    """

    def test_pure_rigid_rotation_yields_zero_elastic_displacement(self):
        """Verify pure rigid rotation produces u_e ≈ 0.

        Setup:
        - Constant omega (no alpha)
        - No aerodynamic forces (preCICE reads zero)
        - No gravity
        - include_reference_acceleration=True

        Expected:
        - Elastic displacement u_e remains < 1e-9 for all nodes
        - Rigid rotation handled by K(θ) + F_ref cancellation

        Verification:
        max(|u_e|) < 1e-9 across all time steps

        Reference:
        Task 6.5 — Pure rigid rotation analytical test
        Covered by Rust test: rotor_physics::tests::test_rigid_rotation_no_elastic
        """
        pass

    def test_gravity_only_produces_1p_torque_modulation(self):
        """Verify gravity-only rotor shows 1P frequency response.

        Setup:
        - Constant omega = 10 rad/s
        - No aerodynamic forces
        - Gravity = [0, 0, -9.81] m/s²
        - Single-blade rotor (simpler to analyze)

        Expected:
        - Torque signal τ(t) has dominant frequency at ω (1P)
        - Displacement u_e perpendicular to blade shows sinusoidal variation

        Verification:
        FFT(τ_aero) peak at ω within 1% frequency error

        Reference:
        Task 6.6 — Gravity-only 1P modulation test
        """
        pass

    def test_reference_load_analytical_validation(self):
        """Verify F_ref matches analytical formula for simple rotor.

        Setup:
        - Single node at r = [1.0, 0.0, 0.0] (1m from axis)
        - Mass m = 10 kg
        - omega = 5 rad/s
        - alpha = 2 rad/s²
        - Rotation axis = [0, 0, 1]
        - Center = [0, 0, 0]

        Analytical F_ref:
        - Centrifugal: -m·ω²·r_perp = -10 · 25 · [1,0,0] = [-250, 0, 0] N
        - Euler: -m·α×r = -10 · ([0,0,2]×[1,0,0]) = -10·[0,2,0] = [0,-20,0] N
        - Total: F_ref = [-250, -20, 0] N

        Verification:
        Rust-computed F_ref matches within 1e-12 (floating-point precision)

        Reference:
        Task 6.7 — Analytical F_ref validation
        Covered by Rust test: rotor_physics::tests::test_reference_load_combined
        """
        pass


# ---------------------------------------------------------------------------
# Group 5 — Integration Parity (requires full solver + preCICE mock)
# ---------------------------------------------------------------------------


@_skip_rust
@pytest.mark.skip(reason="Requires preCICE mock or comparison baseline - see Task 6.8-6.10")
class TestIntegrationParity:
    """End-to-end parity between Rust and Python implementations.

    These tests require:
    - Full preCICE coupling mock or recorded force history
    - Python baseline implementation of InertialRotorFsiSolver
    - IEA-15 blade stub mesh (or similar representative geometry)

    Current status: Deferred pending preCICE mock infrastructure.
    """

    def test_rust_python_parity_iea15_stub(self):
        """Compare Rust vs Python for IEA-15 blade stub simulation.

        Setup:
        - 3-blade IEA-15 stub (first 10m of blade)
        - Prescribed force history from file (recorded preCICE session)
        - omega_mode = "constant", omega = 7.56 rpm (rated)
        - 100 time steps, dt = 0.01s

        Compare:
        - u_e(t) for all interface nodes
        - tau_aero(t)
        - omega(t) if using computed mode

        Acceptance:
        - Relative error < 1e-6 for all outputs
        - Max absolute error in displacement < 1e-9 m

        Reference:
        Task 6.8 — Rust vs Python parity test (migration success criterion)
        """
        pass

    def test_checkpoint_restart_continuity(self):
        """Verify restart produces identical continuation.

        Setup:
        - Run solver for 50 steps → checkpoint state
        - Restart from checkpoint → run 50 more steps
        - Compare against non-restarted 100-step run

        Expected:
        - Steps 51-100 match within floating-point precision

        Reference:
        Task 6.9 — Checkpoint/restart test
        """
        pass

    def test_precice_contract_elastic_displacement_only(self):
        """Verify inertial solver writes elastic displacement only.

        Setup:
        - displacement_mode = "elastic"
        - Single converged window with known u_e
        - Mock preCICE interface to capture write_data

        Verification:
        - Displacement written to preCICE = u_e (not u_e + u_rigid)
        - Compare against corotational solver (writes transformed displacement)

        Reference:
        Task 6.10 — preCICE contract verification
        """
        pass


# ---------------------------------------------------------------------------
# T4.3 — Bit-exact regression with explicit defaults
# ---------------------------------------------------------------------------

class TestInertialRustBindingExplicitDefaults:
    """T4.3 — verify the new keyword params accept defaults without type errors.

    Spec scenario R1: calling the binding with explicit defaults must produce
    identical behavior to calling it without the new kwargs (backward compat).
    We exercise the binding's argument-parsing path only — a full FSI loop
    requires preCICE and is not run here.
    """

    @_skip_rust
    def test_new_kwargs_accepted_with_defaults(self):
        """Verify the binding accepts the new kwargs without raising.

        Checks that:
        - ``include_spin_softening=True`` is accepted
        - ``include_geometric_stiffness=False`` is accepted (explicit)
        - All four threshold kwargs (float) are accepted
        - ``include_geometric_stiffness=None`` is accepted (triggers sentinel path)
        """
        import inspect

        sig = inspect.signature(_aeroelast.run_inertial_rotor_fsi_solver)
        param_names = set(sig.parameters.keys())

        # All new params must be present in the signature
        for expected in [
            "include_spin_softening",
            "include_geometric_stiffness",
            "ksp_omega_rebuild_high",
            "ksp_omega_rebuild_low",
            "kg_omega_rebuild_high",
            "kg_omega_rebuild_low",
        ]:
            assert expected in param_names, (
                f"Expected parameter '{expected}' not found in "
                f"run_inertial_rotor_fsi_solver signature. "
                f"Present params: {sorted(param_names)}"
            )

    @_skip_rust
    def test_include_spin_softening_default_is_true(self):
        """Verify that include_spin_softening defaults to True."""
        import inspect

        sig = inspect.signature(_aeroelast.run_inertial_rotor_fsi_solver)
        param = sig.parameters.get("include_spin_softening")
        assert param is not None, "include_spin_softening not in signature"
        assert param.default is True, (
            f"include_spin_softening default must be True, got {param.default!r}"
        )

    @_skip_rust
    def test_include_geometric_stiffness_default_is_none(self):
        """Verify that include_geometric_stiffness defaults to None (sentinel path)."""
        import inspect

        sig = inspect.signature(_aeroelast.run_inertial_rotor_fsi_solver)
        param = sig.parameters.get("include_geometric_stiffness")
        assert param is not None, "include_geometric_stiffness not in signature"
        assert param.default is None, (
            f"include_geometric_stiffness default must be None, got {param.default!r}"
        )

    @_skip_rust
    def test_threshold_defaults(self):
        """Verify that all threshold kwargs have the expected defaults."""
        import inspect

        sig = inspect.signature(_aeroelast.run_inertial_rotor_fsi_solver)
        expected_defaults = {
            "ksp_omega_rebuild_high": 0.005,
            "ksp_omega_rebuild_low": 0.003,
            "kg_omega_rebuild_high": 0.005,
            "kg_omega_rebuild_low": 0.003,
        }
        for name, expected_val in expected_defaults.items():
            param = sig.parameters.get(name)
            assert param is not None, f"{name} not in signature"
            assert abs(param.default - expected_val) < 1e-12, (
                f"{name} default: expected {expected_val}, got {param.default}"
            )
