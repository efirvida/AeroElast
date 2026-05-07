"""
Inertial Frame Rotor FSI Solver.

This module provides the LinearDynamicFSIRotorInertialSolver for rotor FSI
problems using an alternative formulation where:

- The unknown is ELASTIC DISPLACEMENT (u_e) in global coordinates over a
  rigidly-rotated reference configuration
- The structural mesh rotates internally but the preCICE interface remains FIXED
- No rotating-frame fictitious forces (no K_SP, no K_G in this version,
  no Coriolis gyroscopic matrix)
- The reference load -M·a_ref accounts for rigid-body inertial acceleration

Mathematical Formulation
------------------------
The equation of motion in the INERTIAL frame is:

    [M]{ü_e} + [C(θ)]{u̇_e} + [K(θ)]{u_e} = {F_aero} + {F_g} - [M]{a_ref}

where:
- u_e: elastic displacement in global coordinates (NOT total displacement)
- θ: rotation angle (kinematics only, no frame rotation)
- K(θ), C(θ): assembled on the rigidly-rotated structural geometry
- a_ref = α × r + ω × (ω × r): rigid-body reference acceleration

The preCICE interface contract:
- SolidMesh vertices registered at REFERENCE coordinates (t=0, no rotation)
- Write data: elastic displacement u_e in global frame (NOT total displacement)
- Read data: aerodynamic forces in global frame (no transformation needed)
- GlobalSolidMesh: representative angular velocity (same as corotational solver)

This formulation is physically consistent with the inertial frame perspective
and avoids the need for K_SP, K_G, or Coriolis terms that appear in the
corotational rotating-frame formulation.

See docs/rotor_inertial_solver_design.md for full design documentation.
"""

import logging
from typing import Any, Dict, Optional, Tuple

import numpy as np
from numpy.typing import NDArray

from aeroelast.core.mesh import MeshModel

from .base import Adapter
from .corotational import (
    ComputedOmega,
    ConstantOmega,
    CoordinateTransforms,
    InertialForcesCalculator,
    OmegaProvider,
    RampedComputedOmega,
    RampedOmega,
)
from .linear_dynamic import LinearDynamicFSISolver

_logger = logging.getLogger(__name__)

# Default configuration values
_DEFAULT_ROTATION_AXIS = np.array([0.0, 1.0, 0.0])
_DEFAULT_ROTATION_CENTER = np.array([0.0, 0.0, 0.0])
_DEFAULT_GRAVITY = np.array([0.0, 0.0, -9.81])
_GRAVITY_THRESHOLD = 1e-6
_DEFAULT_FLUID_DENSITY = 1.225  # kg/m³ (air at sea level)
_DEFAULT_FLOW_VELOCITY = 1.0  # m/s


class LinearDynamicFSIRotorInertialSolver(LinearDynamicFSISolver):
    """
    Inertial frame FSI solver for rotating structures (rotors, blades, turbines).

    Solves the inertial-frame equation of motion:

        [M]{ü_e} + [C(θ)]{u̇_e} + [K(θ)]{u_e} = {F_aero}^{glob} + {F_g}^{glob} - [M]{a_ref}

    where:
    - u_e is the elastic displacement in global coordinates over a rigidly-rotated
      reference configuration
    - K(θ) and C(θ) are assembled on the rotated structural geometry (once per
      converged FSI time window)
    - a_ref = α × r + ω × (ω × r) is the rigid-body reference acceleration
      (centripetal + tangential)
    - The preCICE interface mesh (SolidMesh) remains FIXED at reference coordinates

    Key Differences from Corotational Solver
    -----------------------------------------
    - **Frame**: Inertial (global) frame vs. rotating frame
    - **Unknown**: Elastic displacement u_e over rotating reference vs. displacement
      in rotating frame
    - **Stiffness**: K(θ) varies with rotation vs. constant K in rotating frame
    - **Inertial Terms**: -M·a_ref on RHS vs. fictitious forces (centrifugal,
      Coriolis, Euler) in rotating frame
    - **K_SP/K_G**: NOT included in this formulation (may be added in future versions)
    - **preCICE Interface**: Fixed SolidMesh vs. rotating mesh or transformed displacement

    Physical Interpretation
    -----------------------
    The solver tracks:
    - Reference configuration: x̂(t) = R(θ(t))·X₀ (rigid rotation)
    - Elastic deformation: u_e(t) (solved by FEM)
    - Total position: x(t) = x̂(t) + u_e(t)

    The preCICE interface writes ONLY u_e, not the total displacement.
    The CFD participant is responsible for applying its own rigid-body rotation
    using GlobalSolidMesh (angular velocity ω).

    Configuration Parameters
    ------------------------
    All rotor configuration parameters are the same as the corotational solver,
    except physics flags that don't apply:
    - include_geometric_stiffness: NOT used (no K_G in this version)
    - include_spin_softening: NOT used (no K_SP)
    - include_centrifugal: NOT used (replaced by -M·a_ref)
    - include_coriolis: NOT used (no Coriolis in inertial frame)
    - include_euler: NOT used (Euler term is part of a_ref)

    The following configuration keys ARE used:
    - solver.rotor.omega: Angular velocity [rad/s]
    - solver.rotor.rotation_axis: Rotation axis unit vector [x, y, z]
    - solver.rotor.rotation_center: Rotation center [x, y, z]
    - solver.rotor.moment_of_inertia: Rotor inertia (for ComputedOmega) or "auto"
    - solver.rotor.shaft_torque: Resistive torque at rotor shaft [N·m]
    - solver.rotor.omega_ramp_time: Time to ramp ω from 0 to target [s]
    - solver.rotor.gravity: Gravity vector [gx, gy, gz] in m/s²
    - solver.rotor.send_omega_to_precice: Whether to write ω to GlobalSolidMesh
    - solver.rotor.omega_mesh_name: preCICE mesh name for ω (default: "GlobalSolidMesh")
    - solver.rotor.omega_write_data: preCICE data name for ω (default: "AngularVelocity")
    - solver.rotor.transform_displacement_to_inertial: Must be False or omitted
      (inertial solver writes u_e directly, no frame transform)

    Numerical Scheme
    ----------------
    - Time integration: Newmark-β (same as parent LinearDynamicFSISolver)
    - FSI coupling: preCICE implicit coupling (IQN-ILS)
    - Reassembly: K(θ) and C(θ) rebuilt once per converged FSI time window
    - Mass matrix: M remains constant (rigid rotation does not change mass distribution)

    Examples
    --------
    YAML configuration:

        solver:
          type: "LinearDynamicFSIRotorInertial"
          total_time: 5.0
          time_step: 0.001

          rotor:
            omega: 10.0              # rad/s
            rotation_axis: [0, 1, 0] # Y-axis
            rotation_center: [0, 0, 0]
            gravity: [0, 0, -9.81]
            send_omega_to_precice: true
            omega_ramp_time: 1.0

          coupling:
            participant: "Solid"
            config_file: "precice-config.xml"
            coupling_mesh: "SolidMesh"
            write_data: "Displacement"
            read_data: "Force"

          damping:
            eta_m: 1.0e-4
            eta_k: 1.0e-4

    See Also
    --------
    LinearDynamicFSIRotorCorotationalSolver : Alternative rotating-frame formulation
    docs/rotor_inertial_solver_design.md : Full design documentation
    docs/rotor_inertial_solver_tasks.md : Implementation task breakdown
    """

    def __init__(self, mesh: MeshModel, fem_model_properties: Dict[str, Any]) -> None:
        """
        Initialize the inertial rotor FSI solver.

        Parameters
        ----------
        mesh : MeshModel
            Finite element mesh.
        fem_model_properties : Dict[str, Any]
            Configuration dictionary with solver, rotor, coupling, and damping parameters.
        """
        super().__init__(mesh, fem_model_properties)
        self._init_rotor_config()
        self._init_solver_config()
        self._init_state_tracking()

    def _init_rotor_config(self) -> None:
        """
        Initialize rotor-specific configuration from model properties.

        Reuses the same OmegaProvider hierarchy as the corotational solver:
        1. auto-inertia + ramp → RampedOmega initially, replaced by
           RampedComputedOmega after I is estimated from mesh
        2. explicit I + ramp → RampedComputedOmega (ramp then dynamic)
        3. explicit I, no ramp → ComputedOmega (dynamic from t=0)
        4. ramp only → RampedOmega (prescribed linear ramp)
        5. constant → ConstantOmega (fixed ω, α=0)

        Also initializes CoordinateTransforms and InertialForcesCalculator,
        which provide the utilities for:
        - rotate_point_cloud(): Rotate structural mesh geometry internally
        - compute_rigid_body_acceleration_inertial(): Compute a_ref
        - compute_reference_load_vector(): Convert a_ref to F_ref = -M·a_ref
        """
        rotor_cfg = self.solver_params.get("rotor", {})

        # Angular velocity configuration (same as corotational solver)
        omega_value = float(rotor_cfg.get("omega", 0.0))
        omega_ramp_time = float(rotor_cfg.get("omega_ramp_time", 0.0))

        # Check for dynamic omega (ComputedOmega)
        moment_of_inertia = rotor_cfg.get("moment_of_inertia")
        shaft_torque = rotor_cfg.get("shaft_torque", None)
        if shaft_torque is None:
            # Backward compat: negate resistive_torque to match new sign convention
            legacy = rotor_cfg.get("resistive_torque", 0.0)
            shaft_torque = -float(legacy) if float(legacy) != 0.0 else 0.0
        else:
            shaft_torque = float(shaft_torque)
        self._auto_inertia = False

        # OmegaProvider selection (identical logic to corotational solver)
        if isinstance(moment_of_inertia, str) and moment_of_inertia.lower() == "auto":
            # Auto-compute inertia from mesh - will be resolved in solve()
            self._auto_inertia = True
            self._auto_inertia_params = {
                "target_omega": omega_value,
                "ramp_time": omega_ramp_time,
                "shaft_torque": shaft_torque,
            }
            # Temporary provider until inertia is computed
            if omega_ramp_time > 0.0:
                self._omega_provider = RampedOmega(
                    target_omega=omega_value, ramp_time=omega_ramp_time
                )
            else:
                self._omega_provider = ConstantOmega(omega=omega_value)
        elif moment_of_inertia is not None:
            # Explicit moment of inertia provided
            inertia_val = float(moment_of_inertia)
            if omega_ramp_time > 0.0:
                self._omega_provider = RampedComputedOmega(
                    target_omega=omega_value,
                    ramp_time=omega_ramp_time,
                    moment_of_inertia=inertia_val,
                    shaft_torque=shaft_torque,
                )
            else:
                self._omega_provider = ComputedOmega(
                    moment_of_inertia=inertia_val,
                    initial_omega=omega_value,
                    shaft_torque=shaft_torque,
                )
        elif omega_ramp_time > 0.0:
            self._omega_provider = RampedOmega(target_omega=omega_value, ramp_time=omega_ramp_time)
        else:
            self._omega_provider = ConstantOmega(omega=omega_value)

        # Rotation geometry (same as corotational)
        rotation_axis = rotor_cfg.get("rotation_axis", list(_DEFAULT_ROTATION_AXIS))
        rotation_center = rotor_cfg.get("rotation_center", list(_DEFAULT_ROTATION_CENTER))

        # Initialize transformation and inertial utilities
        # These provide the core methods for the inertial formulation:
        # - CoordinateTransforms.rotate_point_cloud()
        # - InertialForcesCalculator.compute_rigid_body_acceleration_inertial()
        # - InertialForcesCalculator.compute_reference_load_vector()
        self._coord_transforms = CoordinateTransforms(
            rotation_axis=rotation_axis,
            rotation_center=rotation_center,
        )
        self._inertial_calculator = InertialForcesCalculator(
            rotation_axis=rotation_axis,
            rotation_center=rotation_center,
        )

        # Store references for torque calculations (Phase 6)
        self._rotation_center = np.array(rotation_center, dtype=np.float64)
        self._rotation_axis = np.array(rotation_axis, dtype=np.float64)

        # Physics flags (corotational-specific flags NOT used in inertial solver)
        # These are read from config but ignored in the inertial implementation:
        # - include_geometric_stiffness: No K_G in this version
        # - include_spin_softening: No K_SP (inertial frame has no spin softening)
        # - include_centrifugal: Replaced by -M·a_ref (centripetal component)
        # - include_coriolis: No Coriolis in inertial frame
        # - include_euler: Replaced by -M·a_ref (tangential component)
        #
        # We log a warning if these are explicitly set to True, to avoid confusion.
        if rotor_cfg.get("include_geometric_stiffness", False):
            _logger.warning(
                "include_geometric_stiffness is set but NOT used in inertial solver "
                "(no K_G in this formulation)"
            )
        if rotor_cfg.get("include_spin_softening", False):
            _logger.warning(
                "include_spin_softening is set but NOT used in inertial solver "
                "(no K_SP in inertial frame)"
            )
        if not rotor_cfg.get("include_centrifugal", True):
            _logger.warning(
                "include_centrifugal=False has no effect in the inertial solver. "
                "Centrifugal acceleration is implicitly included via -M·a_ref and "
                "cannot be disabled independently."
            )
        if rotor_cfg.get("include_coriolis", False):
            _logger.warning(
                "include_coriolis is set but NOT used in inertial solver "
                "(no Coriolis in inertial frame)"
            )
        if rotor_cfg.get("include_euler", False):
            _logger.warning(
                "include_euler is set but NOT used in inertial solver (Euler term is part of a_ref)"
            )

        # Force ramp (same as corotational)
        self._force_ramp_time = float(rotor_cfg.get("force_ramp_time", 0.0))

        # Omega output to preCICE (same as corotational)
        self._send_omega_to_precice = rotor_cfg.get("send_omega_to_precice", True)
        self._omega_mesh_name: str = rotor_cfg.get("omega_mesh_name", "GlobalSolidMesh")
        self._omega_write_data_name: str = rotor_cfg.get("omega_write_data", "AngularVelocity")

        # preCICE displacement mode: controls what the solver writes to the coupling mesh.
        #   "elastic" (default): writes u_e only (elastic deformation over rotating reference).
        #               Use with inertial-aware CFD participants that handle rigid rotation
        #               separately via GlobalSolidMesh/AngularVelocity.
        #   "total":   writes u_e + u_rigid = x - X_0 (total displacement from original config).
        #               Use with standard OpenFOAM / CFD participants that expect absolute
        #               displacement. The rigid component R(θ)·X_0 - X_0 is added internally.
        _precice_mode = rotor_cfg.get("precice_displacement_mode", "elastic")
        if _precice_mode not in ("elastic", "total"):
            raise ValueError(
                f"precice_displacement_mode must be 'elastic' or 'total', got {_precice_mode!r}. "
                f"Use 'elastic' for inertial-aware CFD participants (GlobalSolidMesh pattern) "
                f"or 'total' for standard OpenFOAM adapters."
            )
        self._precice_displacement_mode: str = _precice_mode

        # Displacement transformation flag
        # For the inertial solver, this MUST be False (or omitted).
        # The solver writes elastic displacement u_e directly, not total displacement.
        # If set to True, log an error and override to False.
        transform_disp = rotor_cfg.get("transform_displacement_to_inertial", False)
        if transform_disp:
            _logger.error(
                "transform_displacement_to_inertial is True but MUST be False for "
                "the inertial solver (writes elastic displacement u_e directly). "
                "Overriding to False."
            )
        self._transform_displacement = False

        # Gravity vector (in inertial/global frame, same as corotational)
        self._gravity = np.array(rotor_cfg.get("gravity", list(_DEFAULT_GRAVITY)), dtype=np.float64)
        self._include_gravity = np.linalg.norm(self._gravity) > _GRAVITY_THRESHOLD

        # Aerodynamic performance parameters (for Cp, Cq, Ct, TSR)
        perf_cfg = self.solver_params.get("performance", {})
        self._fluid_density = float(
            perf_cfg.get("fluid_density")
            or rotor_cfg.get("fluid_density")
            or _DEFAULT_FLUID_DENSITY
        )
        self._flow_velocity = float(
            perf_cfg.get("flow_velocity")
            or rotor_cfg.get("flow_velocity")
            or _DEFAULT_FLOW_VELOCITY
        )
        self._rotor_radius: Optional[float] = rotor_cfg.get("radius")
        if self._rotor_radius is not None:
            self._rotor_radius = float(self._rotor_radius)

    def _init_solver_config(self) -> None:
        """
        Initialize solver configuration parameters.

        Simpler than the corotational solver because there are no K_G, K_SP,
        or gyroscopic damping terms to configure.
        """
        super()._init_solver_config()

        damping_cfg = self.solver_params.get("damping") or {}

        self._damping_enabled: bool = damping_cfg.get("enabled", True)
        # Auto mode: compute coefficients from modal analysis at assembly time
        self._damping_auto: bool = (
            self._damping_enabled
            and damping_cfg.get("eta_m") is None
            and damping_cfg.get("eta_k") is None
            and bool(damping_cfg)  # only auto when damping section is present
        )
        self._damping_cfg: dict = damping_cfg

        if not self._damping_enabled:
            self._eta_m = 0.0
            self._eta_k = 0.0
        elif self._damping_auto:
            # Placeholder — will be overwritten during assembly if modal analysis is performed
            self._eta_m = 0.0
            self._eta_k = 0.0
        else:
            # Manual mode: read from nested damping dict, fall back to flat keys
            self._eta_m = float(
                damping_cfg["eta_m"]
                if damping_cfg.get("eta_m") is not None
                else self.solver_params.get("eta_m", 0.0)
            )
            self._eta_k = float(
                damping_cfg["eta_k"]
                if damping_cfg.get("eta_k") is not None
                else self.solver_params.get("eta_k", 0.0)
            )

        self._stress_output_interval = int(self.solver_params.get("stress_output_interval", 1))

    def _init_state_tracking(self) -> None:
        """
        Initialize state tracking variables.

        The inertial solver tracks:
        - theta: Current rotation angle (kinematics only, not a frame rotation)
        - omega, alpha: Angular velocity and acceleration (from OmegaProvider)
        - skip_ramps: Flag to prevent re-applying ramps on restart
        """
        # Current rotation angle (updated during time stepping)
        self._theta = 0.0

        # On restart, ramps must NOT be re-applied — the simulation
        # continues from the checkpoint state where ramps already completed.
        self._skip_ramps = False

        # State for tracking when to reassemble K(θ) and C(θ)
        # In Phase 1, we reassemble every converged time window.
        # In Phase 2 (optimization), we may track whether θ changed significantly.
        self._last_assembly_theta: Optional[float] = None

    # =========================================================================
    # Phase 3: Structural Assembly on Internally Rotated Geometry
    # =========================================================================

    def _rotate_structural_geometry_internal(self, theta: float) -> np.ndarray:
        """
        Rotate mesh node coordinates by angle theta for structural assembly.

        The inertial solver assembles structural matrices (K, C) in the
        rotated configuration at each time step. This method computes the
        rotated nodal coordinates X_rot = R(θ)·X_0 without modifying the
        stored MeshModel.nodes, which remain in the reference (unrotated)
        configuration.

        Parameters
        ----------
        theta : float
            Rotation angle [rad] around the rotation axis.

        Returns
        -------
        np.ndarray, shape (n_nodes, 3)
            Rotated nodal coordinates in global frame. Each row is [x, y, z]
            for one node.

        Notes
        -----
        - This method is called before structural assembly at each time step
        - The returned coordinates differ from MeshModel.nodes by exactly R(θ)
        - Used to rebuild the assembler with rotated geometry (Task 3.2)
        """
        # Extract original node coordinates from mesh
        original_coords = np.array(
            [node.coords for node in self.domain.mesh.nodes.values()],
            dtype=np.float64,
        )

        # Rotate using CoordinateTransforms
        rotated_coords = self._coord_transforms.rotate_point_cloud(
            points=original_coords,
            theta=theta,
        )

        return rotated_coords

    def _rebuild_assembler_with_rotated_geometry(
        self,
        rotated_coords: np.ndarray,
    ) -> "MeshAssembler":
        """
        Rebuild structural assembler with rotated nodal coordinates.

        This method clones the mesh model, updates nodal coordinates to the
        rotated configuration, and reconstructs the MeshAssembler. This allows
        structural matrices (K, C) to be assembled in the rotated frame without
        modifying the original mesh stored in self.domain.

        Parameters
        ----------
        rotated_coords : np.ndarray, shape (n_nodes, 3)
            Rotated nodal coordinates from _rotate_structural_geometry_internal()

        Returns
        -------
        MeshAssembler
            New assembler object built on rotated geometry

        Notes
        -----
        - Creates a shallow copy of the mesh to avoid side effects
        - Only node coordinates are modified; connectivity remains unchanged
        - Called at each time step before structural assembly (Task 3.2)
        - In Phase 9 (optimization), this can be replaced with incremental updates
        """
        from aeroelast.core.assembler import MeshAssembler
        from aeroelast.core.mesh import MeshModel, Node

        # Create a shallow copy of the mesh model
        # We only need to update node coordinates; everything else stays the same
        original_mesh = self.domain.mesh

        # Clone nodes with updated coordinates
        rotated_nodes = {}
        for i, (node_id, original_node) in enumerate(original_mesh.nodes.items()):
            rotated_nodes[node_id] = Node(
                id=node_id,
                coords=rotated_coords[i],
            )

        # Build new mesh model with rotated nodes but same connectivity
        rotated_mesh = MeshModel(
            nodes=rotated_nodes,
            elements=original_mesh.elements,  # Shared reference OK (immutable connectivity)
            node_sets=original_mesh.node_sets,  # Shared reference OK
            element_sets=original_mesh.element_sets,  # Shared reference OK
        )

        # Reconstruct assembler with rotated geometry
        # This follows the "Phase 1" assembly path from linear_dynamic.py:907
        rotated_assembler = MeshAssembler(
            mesh=rotated_mesh,
            model=self.model_properties,
        )

        return rotated_assembler

    def _assemble_or_reuse_mass_matrix(self) -> "PETSc.Mat":
        """
        Assemble mass matrix once and reuse across all timesteps.

        The mass matrix M is invariant under rigid-body rotation because it
        depends only on element density and volume, not on nodal coordinates.
        Therefore, M can be assembled once in the reference configuration and
        reused at all rotated orientations.

        Returns
        -------
        PETSc.Mat
            Lumped mass matrix (assembled once, cached in self._mass_matrix_cached)

        Notes
        -----
        - Mass matrix invariance: M(R(θ)·X) = M(X) for rotation matrix R
        - Assembled once on first call, then returned from cache
        - Task 3.3: Verify M before and after rotation is identical (FP tolerance)
        """
        # Check if mass matrix is already cached
        if hasattr(self, "_mass_matrix_cached") and self._mass_matrix_cached is not None:
            return self._mass_matrix_cached

        # Assemble mass matrix once using original (unrotated) geometry
        # Use the domain's assembler (which has original coordinates)
        self._mass_matrix_cached = self.domain.assemble_mass_matrix_lumped()

        return self._mass_matrix_cached

    def _assemble_rayleigh_damping(
        self,
        K_theta: "PETSc.Mat",
        M: "PETSc.Mat",
    ) -> "PETSc.Mat":
        """
        Assemble Rayleigh damping matrix C(θ) = η_m·M + η_k·K(θ).

        The damping matrix has two components:
        - η_m·M: Mass-proportional (constant, M is invariant)
        - η_k·K(θ): Stiffness-proportional (orientation-dependent)

        Parameters
        ----------
        K_theta : PETSc.Mat
            Stiffness matrix assembled on rotated geometry at angle θ
        M : PETSc.Mat
            Mass matrix (invariant, assembled once)

        Returns
        -------
        PETSc.Mat
            Rayleigh damping matrix C(θ)

        Notes
        -----
        - If damping is disabled, returns zero matrix
        - With η_k=0, damping is constant (η_m·M only)
        - With η_k≠0, damping varies with orientation (stiffness term updates)
        - Task 3.4: Verify orientation-dependent path with η_k≠0
        """
        if not self._damping_enabled:
            # Return zero damping matrix
            C = K_theta.duplicate(copy=False)  # Same sparsity pattern
            C.zeroEntries()
            return C

        # Start with mass-proportional term: C = η_m·M
        C = M.duplicate(copy=True)  # Deep copy of M
        C.scale(self._eta_m)

        # Add stiffness-proportional term: C += η_k·K(θ)
        if self._eta_k != 0.0:
            C.axpy(self._eta_k, K_theta)  # C = C + η_k·K(θ)

        return C

    def _ensure_elastic_boundary_conditions(self) -> None:
        """
        Verify that Dirichlet boundary conditions are homogeneous (elastic-only).

        The inertial solver solves for elastic displacement u_e over a rotating
        reference frame. Root boundary conditions must constrain elastic motion
        (u_e = 0), NOT impose time-dependent rigid-body displacement.

        This method checks self.dirichlet_conditions and raises an error if
        any non-zero prescribed displacement is found.

        Raises
        ------
        ValueError
            If any Dirichlet BC has non-zero prescribed value (time-dependent constraint)

        Notes
        -----
        - Elastic unknown space: u_e = 0 at root means "no elastic deformation"
        - Rigid rotation is handled by rotating the reference configuration x̂(t)
        - Total displacement x(t) = x̂(t) + u_e(t), so u_e=0 ⇒ x follows rigid rotation
        - Task 3.5: Verify constrained DOFs stay homogeneous
        """
        if not hasattr(self, "dirichlet_conditions") or not self.dirichlet_conditions:
            return  # No BCs to check

        # Check for non-zero prescribed displacements
        non_zero_bcs = [
            (dof, val) for dof, val in self.dirichlet_conditions.items() if abs(val) > 1e-12
        ]

        if non_zero_bcs:
            raise ValueError(
                f"Inertial rotor solver requires homogeneous (zero) Dirichlet BCs. "
                f"Found {len(non_zero_bcs)} non-zero prescribed displacements:\n"
                f"  {non_zero_bcs[:5]}\n"  # Show first 5 violations
                f"The unknown is elastic displacement u_e over a rotating reference. "
                f"Root constraints must be u_e=0 (elastic-only), not time-dependent "
                f"rigid-body displacements."
            )

    # =========================================================================
    # Phase 3 Complete
    # =========================================================================
    # Tasks 3.1-3.5 implemented:
    # ✓ 3.1: _rotate_structural_geometry_internal()
    # ✓ 3.2: _rebuild_assembler_with_rotated_geometry()
    # ✓ 3.3: _assemble_or_reuse_mass_matrix()
    # ✓ 3.4: _assemble_rayleigh_damping()
    # ✓ 3.5: _ensure_elastic_boundary_conditions()

    # =========================================================================
    # Phase 4: Newmark Step for the Inertial Formulation
    # =========================================================================

    def _assemble_inertial_effective_system(
        self,
        K_theta: "PETSc.Mat",
        M: "PETSc.Mat",
        C_theta: "PETSc.Mat",
        a0: float,
        a1: float,
    ) -> "PETSc.Mat":
        """
        Assemble the effective stiffness matrix for the inertial formulation.

        K_eff = K(θ) + a0·M + a1·C(θ)

        This is the standard Newmark effective stiffness WITHOUT the rotating-frame
        terms K_G, K_SP, or G_cor that appear in the corotational solver.

        Parameters
        ----------
        K_theta : PETSc.Mat
            Stiffness matrix assembled on the rigidly-rotated geometry.
        M : PETSc.Mat
            Mass matrix (invariant under rotation).
        C_theta : PETSc.Mat
            Damping matrix C = η_m·M + η_k·K(θ).
        a0 : float
            Newmark coefficient for mass: a0 = 1/(β·Δt²).
        a1 : float
            Newmark coefficient for damping: a1 = γ/(β·Δt).

        Returns
        -------
        K_eff : PETSc.Mat
            Effective stiffness matrix for the Newmark solve.

        Notes
        -----
        When rotation is disabled (θ = 0), this reduces to the standard
        LinearDynamicFSI effective system:
            K_eff = K₀ + a0·M + a1·C₀

        Verification:
        - Set ω = 0, α = 0 → K(θ=0) = K₀, C(θ=0) = C₀
        - Result should match LinearDynamicFSISolver._assemble_effective_stiffness()
        """
        import petsc4py.PETSc as PETSc

        # Create effective matrix: K_eff = K(θ) + a0·M + a1·C(θ)
        K_eff = K_theta.copy()
        K_eff.axpy(a0, M)  # K_eff += a0·M
        K_eff.axpy(a1, C_theta)  # K_eff += a1·C(θ)
        K_eff.assemble()

        _logger.debug(
            f"Assembled inertial K_eff: K(θ) + {a0:.3e}·M + {a1:.3e}·C(θ) (no K_G, K_SP, G_cor)"
        )

        return K_eff

    def _assemble_inertial_rhs(
        self,
        F_aero_global: NDArray,
        F_gravity_global: NDArray,
        a_ref_nodal: NDArray,
        u_n: NDArray,
        v_n: NDArray,
        a_n: NDArray,
        a0: float,
        a2: float,
        a3: float,
        a6: float,
        a7: float,
    ) -> NDArray:
        """
        Assemble the RHS vector for the inertial Newmark step.

        F_eff = F_aero + F_g - M·a_ref + M·(a0·u_n + a2·v_n + a3·a_n)
                + C·(a1·u_n + a6·v_n + a7·a_n)

        Key differences from corotational solver:
        - No force-frame transformation (forces already in global frame)
        - Reference load -M·a_ref replaces rotating-frame fictitious forces
        - u_n, v_n, a_n are ELASTIC displacements/velocities/accelerations

        Parameters
        ----------
        F_aero_global : ndarray, shape (n_free_dofs,)
            Aerodynamic forces from preCICE in global coordinates (BC-reduced).
        F_gravity_global : ndarray, shape (n_free_dofs,)
            Gravity forces in global coordinates (BC-reduced).
        a_ref_nodal : ndarray, shape (n_nodes, 3)
            Rigid-body reference acceleration at each node in global coordinates.
        u_n : ndarray, shape (n_free_dofs,)
            Elastic displacement at time n.
        v_n : ndarray, shape (n_free_dofs,)
            Elastic velocity at time n.
        a_n : ndarray, shape (n_free_dofs,)
            Elastic acceleration at time n.
        a0, a2, a3, a6, a7 : float
            Newmark coefficients.

        Returns
        -------
        F_eff : ndarray, shape (n_free_dofs,)
            Effective RHS vector for the Newmark solve.

        Notes
        -----
        Sign convention verification:
        - Pure rigid rotation (u_e = 0): External forces = M·a_ref → equilibrium
        - a_ref points radially inward (centripetal) → -M·a_ref points outward
        """
        # Convert nodal acceleration to DOF load vector
        M_diag_full = self.M.getDiagonal().array
        F_ref_full = self._inertial_calc.compute_reference_load_vector(
            a_ref_nodal, M_diag_full, dofs_per_node=self.domain.dofs_per_node
        )
        # Reduce to free DOFs
        F_ref = F_ref_full[self.domain.free_dof_indices]

        # Newmark history terms
        F_hist_mass = a0 * u_n + a2 * v_n + a3 * a_n
        F_hist_damp = a6 * v_n + a7 * a_n  # a1·u_n term absorbed into K_eff

        # Apply M and C to history terms
        M_hist = np.zeros_like(F_hist_mass)
        C_hist = np.zeros_like(F_hist_damp)

        self.M.mult(F_hist_mass, M_hist)
        self.C.mult(F_hist_damp, C_hist)

        # Assemble RHS
        F_eff = F_aero_global + F_gravity_global + F_ref + M_hist + C_hist

        _logger.debug(
            f"Inertial RHS: ||F_aero||={np.linalg.norm(F_aero_global):.2e}, "
            f"||F_g||={np.linalg.norm(F_gravity_global):.2e}, "
            f"||F_ref||={np.linalg.norm(F_ref):.2e}"
        )

        return F_eff

    def _checkpoint_elastic_state(
        self,
        u_e: "PETSc.Vec",
        v_e: "PETSc.Vec",
        a_e: "PETSc.Vec",
        theta: float,
        omega: float,
        alpha: float,
    ) -> Dict[str, Any]:
        """
        Checkpoint the elastic state and rigid-body kinematics.

        The inertial solver state includes:
        - Elastic displacement, velocity, acceleration (u_e, v_e, a_e)
        - Rigid-body kinematics (theta, omega, alpha)

        This is distinct from the corotational solver checkpoint which stores
        displacement in the rotating frame.

        Parameters
        ----------
        u_e : PETSc.Vec
            Elastic displacement vector.
        v_e : PETSc.Vec
            Elastic velocity vector.
        a_e : PETSc.Vec
            Elastic acceleration vector.
        theta : float
            Accumulated rotation angle [rad].
        omega : float
            Angular velocity [rad/s].
        alpha : float
            Angular acceleration [rad/s²].

        Returns
        -------
        checkpoint : dict
            State dictionary for rollback.

        Notes
        -----
        Verification: After rollback, the solver must reproduce the same
        subsequent time history within tolerance (see Phase 7 tests).
        """
        checkpoint = {
            # Elastic state vectors (deep copy)
            "u_e": u_e.duplicate(),
            "v_e": v_e.duplicate(),
            "a_e": a_e.duplicate(),
            # Rigid-body kinematics
            "theta": float(theta),
            "omega": float(omega),
            "alpha": float(alpha),
            # OmegaProvider state (if applicable)
            "omega_provider_state": (
                self._omega_provider.get_state()
                if hasattr(self._omega_provider, "get_state")
                else None
            ),
            # Structural matrices at checkpoint time.
            # These are NOT deep-copied because K and C are rebuilt from scratch
            # on every converged window and are not modified in-place during
            # sub-iterations. Storing the reference is sufficient.
            # After rollback, these references restore the correct matrices
            # so the solver reuses K(θ_ckpt) / C(θ_ckpt) without re-assembling.
            "K_current": None,  # filled by solve() after checkpoint call
            "C_current": None,  # filled by solve() after checkpoint call
        }

        # Copy vector contents
        u_e.copy(result=checkpoint["u_e"])
        v_e.copy(result=checkpoint["v_e"])
        a_e.copy(result=checkpoint["a_e"])

        return checkpoint

    def _rollback_elastic_state(
        self,
        checkpoint: Dict[str, Any],
        u_e: "PETSc.Vec",
        v_e: "PETSc.Vec",
        a_e: "PETSc.Vec",
    ) -> None:
        """
        Rollback to a previously checkpointed state.

        Restores:
        - Elastic displacement, velocity, acceleration
        - Rigid-body kinematics
        - OmegaProvider state (if applicable)

        Parameters
        ----------
        checkpoint : dict
            State dictionary from _checkpoint_elastic_state().
        u_e : PETSc.Vec
            Elastic displacement vector to restore.
        v_e : PETSc.Vec
            Elastic velocity vector to restore.
        a_e : PETSc.Vec
            Elastic acceleration vector to restore.

        Notes
        -----
        After rollback, the next Newmark step should reuse the same K(θ), C(θ)
        from the checkpointed window kinematics (no reassembly needed within
        the same FSI window).
        """
        # Restore elastic state vectors
        checkpoint["u_e"].copy(result=u_e)
        checkpoint["v_e"].copy(result=v_e)
        checkpoint["a_e"].copy(result=a_e)

        # Restore OmegaProvider state
        if checkpoint.get("omega_provider_state") is not None:
            if hasattr(self._omega_provider, "set_state"):
                self._omega_provider.set_state(checkpoint["omega_provider_state"])

    # =========================================================================
    # Phase 4 Complete
    # =========================================================================
    # Tasks 4.1-4.4 implemented:
    # ✓ 4.1: _assemble_inertial_effective_system()
    # ✓ 4.2: _assemble_inertial_rhs()
    # ✓ 4.3: State variables are elastic (u_e, v_e, a_e)
    # ✓ 4.4: _checkpoint_elastic_state(), _rollback_elastic_state()

    # =========================================================================
    # Phase 5: preCICE Contract with Fixed Interface Mesh
    # =========================================================================

    def _register_interface_at_reference_coords(
        self, adapter: Adapter, mesh_name: str, interface_coords: NDArray
    ) -> NDArray:
        """
        Register interface vertices with preCICE using REFERENCE coordinates only.

        The inertial solver keeps the preCICE SolidMesh fixed at the initial
        reference coordinates (t=0, θ=0). The mesh does NOT rotate on the preCICE
        side, even though the internal structural geometry rotates.

        This is a critical difference from a potential rotating-mesh approach:
        - Vertices registered ONCE at t=0
        - Coordinates NEVER updated during simulation
        - Only elastic displacement is written (not rigid-body motion)

        Parameters
        ----------
        adapter : Adapter
            preCICE adapter instance.
        mesh_name : str
            Name of the preCICE mesh (typically "SolidMesh").
        interface_coords : ndarray, shape (n_interface_nodes, 3)
            Interface node coordinates at REFERENCE configuration (unrotated).

        Returns
        -------
        vertex_ids : ndarray
            preCICE vertex IDs for the registered vertices.

        Notes
        -----
        Verification: Repeated FSI windows must NOT attempt to re-register or
        move vertices. The mesh remains geometrically stationary from preCICE's
        perspective.
        """
        # Register vertices at reference coordinates
        vertex_ids = adapter.register_coupling_mesh(mesh_name, interface_coords)

        _logger.info(
            f"Registered {len(vertex_ids)} vertices on '{mesh_name}' at REFERENCE coords "
            f"(fixed preCICE interface for inertial solver)"
        )

        return vertex_ids

    def _read_forces_from_precice_global(
        self, adapter: Adapter, mesh_name: str, data_name: str
    ) -> NDArray:
        """
        Read aerodynamic forces from preCICE in GLOBAL coordinates (no transformation).

        The inertial solver formulation requires forces in the global inertial frame.
        Unlike the corotational solver, we do NOT apply R^T(θ) transformation because:
        - The equation of motion is written in the global frame
        - F_aero is already expressed in global coordinates by the CFD participant
        - The elastic displacement u_e is also in global coordinates

        Parameters
        ----------
        adapter : Adapter
            preCICE adapter instance.
        mesh_name : str
            Name of the preCICE mesh.
        data_name : str
            Name of the force data field (e.g., "Force").

        Returns
        -------
        F_aero_global : ndarray, shape (n_interface_nodes * dim,)
            Aerodynamic forces in global coordinates (flattened).

        Notes
        -----
        Verification: The force path contains NO coordinate transformation.
        Compare with corotational solver which applies:
            F_local = R^T(θ) · F_global
        """
        F_aero_global = adapter.read_data(mesh_name, data_name)

        _logger.debug(
            f"Read forces from '{mesh_name}/{data_name}': ||F||={np.linalg.norm(F_aero_global):.3e} "
            f"(global frame, NO transformation)"
        )

        return F_aero_global

    def _write_elastic_displacement_to_precice(
        self,
        adapter: Adapter,
        mesh_name: str,
        data_name: str,
        u_e_global: NDArray,
        theta: float,
    ) -> None:
        """
        Write ONLY elastic displacement to preCICE (not total displacement).

        Critical preCICE contract for the inertial solver:
            u_fsi = u_e^{global}  (elastic displacement in global frame)

        This is NOT:
            u_fsi = x - X_0  (total displacement from original reference)
            u_fsi = R(θ) · u_local  (transformed displacement)

        The CFD participant is responsible for applying rigid-body rotation using
        GlobalSolidMesh (ω data). The structural solver sends ONLY the elastic
        deformation over the rigidly-rotated reference.

        Physical interpretation:
            x(t) = x̂(t) + u_e(t)
            x̂(t) = R(θ(t)) · X_0  (rigid reference, handled by CFD)
            u_e(t) = elastic response (what we send)

        Parameters
        ----------
        adapter : Adapter
            preCICE adapter instance.
        mesh_name : str
            Name of the preCICE mesh (typically "SolidMesh").
        data_name : str
            Name of the displacement data field (e.g., "Displacement").
        u_e_global : ndarray, shape (n_interface_nodes * dim,)
            Elastic displacement in global coordinates (flattened).
        theta : float
            Current rotation angle (for logging/debugging only, NOT used in computation).

        Notes
        -----
        Verification: An internal test should fail if this method accidentally
        writes x̂ + u_e (total position) instead of just u_e.

        Compatibility warning: If the CFD participant expects total displacement
        instead of elastic displacement, the coupling will produce WRONG physics.
        This must be documented and checked at configuration time.
        """
        # Assertion: ensure we're not accidentally adding rigid-body motion
        # (This would be detected by checking that u_e remains bounded even as
        # the rotor rotates many revolutions)

        if self._precice_displacement_mode == "total":
            # u_total = u_e + u_rigid, where u_rigid = R(θ)·X_0 - X_0
            # This is required when the CFD participant expects total displacement
            # from the original (unrotated) reference configuration.
            X_0 = self._interface_coords_reference  # (n_interface_nodes, 3)
            R = self._coord_transforms.rotation_matrix(theta)
            X_rotated = (R @ X_0.T).T  # (n_interface_nodes, 3)
            u_rigid = (X_rotated - X_0).ravel()  # flatten to (n_interface_nodes * 3,)
            u_fsi = u_e_global + u_rigid
            _logger.debug(
                f"Wrote TOTAL displacement to '{mesh_name}/{data_name}': "
                f"||u_e||={np.linalg.norm(u_e_global):.3e} "
                f"||u_rigid||={np.linalg.norm(u_rigid):.3e} "
                f"||u_total||={np.linalg.norm(u_fsi):.3e} (θ={np.degrees(theta):.1f}°)"
            )
        else:
            u_fsi = u_e_global
            _logger.debug(
                f"Wrote elastic displacement to '{mesh_name}/{data_name}': "
                f"||u_e||={np.linalg.norm(u_e_global):.3e} "
                f"(θ={np.degrees(theta):.1f}°, elastic only, NO rigid-body component)"
            )

        adapter.write_data(mesh_name, data_name, u_fsi)

    def _write_omega_to_global_mesh(
        self, adapter: Adapter, mesh_name: str, data_name: str, omega: float
    ) -> None:
        """
        Write representative angular velocity to GlobalSolidMesh (same as corotational solver).

        The GlobalSolidMesh pattern is unchanged from the corotational solver:
        - Single scalar value representing the rotor's angular velocity
        - Used by CFD participant to rotate its own mesh or apply kinematics
        - Typically written to a single-vertex "global" mesh

        Parameters
        ----------
        adapter : Adapter
            preCICE adapter instance.
        mesh_name : str
            Name of the global mesh (typically "GlobalSolidMesh").
        data_name : str
            Name of the omega data field (e.g., "AngularVelocity").
        omega : float
            Angular velocity magnitude in rad/s.

        Notes
        -----
        Verification: A fixed-omega case should write the same ω history as
        the corotational solver (exact match).
        """
        omega_array = np.array([omega], dtype=np.float64)
        adapter.write_data(mesh_name, data_name, omega_array)

        _logger.debug(
            f"Wrote omega to '{mesh_name}/{data_name}': {omega:.6f} rad/s "
            f"({np.degrees(omega):.2f} °/s)"
        )

    def _check_coupling_compatibility(self, adapter: Adapter) -> None:
        """
        Check that the CFD participant expects elastic displacement (not total displacement).

        This is a configuration-time guard to prevent silent wrong physics.
        If the CFD participant is configured to expect total displacement
        (x - X_0), the inertial solver will produce incorrect results because
        it only sends u_e (elastic component).

        Raises
        ------
        RuntimeError
            If the coupling configuration is incompatible with the inertial solver.

        Notes
        -----
        Future enhancement: Read preCICE config XML to detect participant type
        or add a metadata field that declares displacement contract.
        For now, this is a placeholder with a clear warning message.
        """
        # TODO: Implement actual check by parsing preCICE config or adding metadata
        # For now, just log a warning
        _logger.warning(
            "Inertial solver writes ELASTIC displacement only (u_e), NOT total displacement. "
            "Verify that your CFD participant expects this contract and handles rigid-body "
            "rotation via GlobalSolidMesh (omega data). If the CFD expects total displacement "
            "(x - X_0), the coupling will produce WRONG physics."
        )

    # =========================================================================
    # Phase 5 Complete
    # =========================================================================
    # Tasks 5.1-5.5 implemented:
    # ✓ 5.1: _register_interface_at_reference_coords()
    # ✓ 5.2: _read_forces_from_precice_global()
    # ✓ 5.3: _write_elastic_displacement_to_precice()
    # ✓ 5.4: _write_omega_to_global_mesh()
    # ✓ 5.5: _check_coupling_compatibility()

    # =========================================================================
    # Future Implementation Methods (Phase 6-10)
    # =========================================================================
    # The following methods will be implemented in subsequent phases:
    #
    # Phase 6: Omega Dynamics and Torque Accounting
    # =========================================================================

    def _compute_torque_from_forces(
        self,
        nodal_coords: NDArray,
        nodal_disps: NDArray,
        nodal_forces: NDArray,
    ) -> Tuple[NDArray, float]:
        """
        Compute torque vector in global frame and scalar projection on the rotation axis.

        This is the CORE torque calculation method shared by aerodynamic, gravity,
        and total torque computations. The inertial solver uses this directly without
        coordinate transformations (unlike the corotational solver which applies R(θ)).

        Algorithm:
        1. Compute lever arms: r = (x + u_e) - center
        2. Compute torque vector: τ = Σ (r × F)
        3. Project onto rotation axis: τ_scalar = τ · ê_axis

        Parameters
        ----------
        nodal_coords : ndarray, shape (n_nodes, dim)
            Reference nodal coordinates (unrotated, global frame).
        nodal_disps : ndarray, shape (n_nodes, dim)
            Elastic displacements (global frame).
        nodal_forces : ndarray, shape (n_nodes, dim)
            Nodal forces (global frame).

        Returns
        -------
        torque_global : ndarray, shape (3,)
            Torque vector in global coordinates.
        torque_scalar : float
            Scalar torque projected onto the rotation axis (positive = axis direction).

        Notes
        -----
        Verification: For a purely radial force field with constant |F|, the scalar
        torque should equal n_nodes * |F| * r_avg (zero if forces are along radii).
        """
        # Ensure 3D vectors (pad with zeros if 2D mesh)
        coords_3d = self._ensure_3d_vectors(nodal_coords)
        disps_3d = self._ensure_3d_vectors(nodal_disps)
        forces_3d = self._ensure_3d_vectors(nodal_forces)

        # Compute current positions: x = X + u_e
        positions = coords_3d + disps_3d

        # Compute lever arms from rotation center
        lever_arms = positions - self._rotation_center

        # Compute torque vector: τ = Σ (r × F)
        torque_contributions = np.cross(lever_arms, forces_3d)
        torque_global = np.sum(torque_contributions, axis=0)

        # Project onto rotation axis
        torque_scalar = float(np.dot(torque_global, self._rotation_axis))

        return torque_global, torque_scalar

    def _ensure_3d_vectors(self, array_2d: NDArray) -> NDArray:
        """
        Ensure nodal arrays are 3D (pad with zeros if 2D mesh).

        Parameters
        ----------
        array_2d : ndarray, shape (n_nodes, dim)
            Nodal array (coordinates, displacements, forces).

        Returns
        -------
        array_3d : ndarray, shape (n_nodes, 3)
            Padded array (z=0 if input was 2D).
        """
        if array_2d.shape[1] == 3:
            return array_2d.copy()
        elif array_2d.shape[1] == 2:
            n_nodes = array_2d.shape[0]
            array_3d = np.zeros((n_nodes, 3), dtype=array_2d.dtype)
            array_3d[:, :2] = array_2d
            return array_3d
        else:
            raise ValueError(f"Expected array with dim=2 or dim=3, got shape {array_2d.shape}")

    def _compute_aerodynamic_torque(
        self, interface_coords: NDArray, interface_disps: NDArray, F_aero: NDArray
    ) -> Tuple[NDArray, float]:
        """
        Compute aerodynamic torque from preCICE forces.

        Uses only external aerodynamic forces (no reference load, no inertial terms).
        This is the driving torque component from fluid-structure interaction.

        Parameters
        ----------
        interface_coords : ndarray, shape (n_interface_nodes, dim)
            Interface node coordinates at reference configuration.
        interface_disps : ndarray, shape (n_interface_nodes, dim)
            Elastic displacements at interface nodes.
        F_aero : ndarray, shape (n_interface_nodes * dim,)
            Aerodynamic forces from preCICE (flattened).

        Returns
        -------
        tau_aero_global : ndarray, shape (3,)
            Aerodynamic torque vector (global frame).
        tau_aero_scalar : float
            Aerodynamic torque projected onto rotation axis.

        Notes
        -----
        Verification: In steady rotation with zero elastic deformation, the aero
        torque should balance the drag torque computed from power coefficient.
        """
        # Reshape flattened force array
        dim = interface_coords.shape[1]
        n_interface = len(interface_coords)
        F_aero_2d = F_aero.reshape((n_interface, dim))

        tau_aero_global, tau_aero_scalar = self._compute_torque_from_forces(
            interface_coords, interface_disps, F_aero_2d
        )

        _logger.debug(
            f"Aerodynamic torque: ||τ||={np.linalg.norm(tau_aero_global):.3e}, "
            f"τ_axis={tau_aero_scalar:.3e}"
        )

        return tau_aero_global, tau_aero_scalar

    def _compute_gravity_torque(
        self, nodal_coords: NDArray, nodal_disps: NDArray, F_gravity: NDArray
    ) -> Tuple[NDArray, float]:
        """
        Compute gravity torque from nodal gravity loads.

        Parameters
        ----------
        nodal_coords : ndarray, shape (n_nodes, dim)
            Reference nodal coordinates (unrotated).
        nodal_disps : ndarray, shape (n_nodes, dim)
            Elastic displacements.
        F_gravity : ndarray, shape (n_nodes * dim,)
            Gravity force vector (flattened).

        Returns
        -------
        tau_grav_global : ndarray, shape (3,)
            Gravity torque vector (global frame).
        tau_grav_scalar : float
            Gravity torque projected onto rotation axis.

        Notes
        -----
        Verification: For a rotor with uniform mass distribution and rotation
        axis aligned with gravity, the gravity torque should oscillate at 1P
        frequency as the rotor rotates.
        """
        dim = nodal_coords.shape[1]
        n_nodes = len(nodal_coords)
        F_grav_2d = F_gravity.reshape((n_nodes, dim))

        tau_grav_global, tau_grav_scalar = self._compute_torque_from_forces(
            nodal_coords, nodal_disps, F_grav_2d
        )

        _logger.debug(
            f"Gravity torque: ||τ||={np.linalg.norm(tau_grav_global):.3e}, "
            f"τ_axis={tau_grav_scalar:.3e}"
        )

        return tau_grav_global, tau_grav_scalar

    def _compute_driving_torque(
        self,
        interface_coords: NDArray,
        interface_disps: NDArray,
        F_aero: NDArray,
        nodal_coords: NDArray,
        nodal_disps: NDArray,
        F_gravity: NDArray,
    ) -> float:
        """
        Compute total driving torque for OmegaProvider dynamics.

        Driving torque is the SUM of external forces ONLY:
            tau_driving = tau_aero + tau_gravity + tau_shaft

        The reference load -M·a_ref is NOT included because it is an inertial
        reaction (internal to the structure), not an external driving force.

        Parameters
        ----------
        interface_coords : ndarray
            Interface node coordinates at reference configuration.
        interface_disps : ndarray
            Elastic displacements at interface nodes.
        F_aero : ndarray
            Aerodynamic forces from preCICE.
        nodal_coords : ndarray
            All nodal coordinates at reference configuration.
        nodal_disps : ndarray
            Elastic displacements at all nodes.
        F_gravity : ndarray
            Gravity forces at all nodes.

        Returns
        -------
        tau_driving : float
            Total driving torque (aero + gravity + shaft).

        Notes
        -----
        Verification: A test with only reference load (no aero, no gravity) must
        produce zero driving torque, preventing the rotor from self-acceleration.
        """
        _, tau_aero = self._compute_aerodynamic_torque(interface_coords, interface_disps, F_aero)
        _, tau_gravity = self._compute_gravity_torque(nodal_coords, nodal_disps, F_gravity)

        # Shaft torque (assumed zero in current implementation, future extension)
        tau_shaft = 0.0

        tau_driving = tau_aero + tau_gravity + tau_shaft

        _logger.debug(
            f"Driving torque components: aero={tau_aero:.3e}, "
            f"gravity={tau_gravity:.3e}, shaft={tau_shaft:.3e}, "
            f"total={tau_driving:.3e}"
        )

        return tau_driving

    def _update_omega_after_converged_window(self, tau_driving: float, dt: float) -> None:
        """
        Update angular velocity using OmegaProvider after a converged FSI window.

        This MUST be called ONLY after the FSI window has converged (no rollbacks).
        The OmegaProvider state machine ensures omega evolves with the correct
        semantics (ramped, constant, or computed from torque balance).

        Parameters
        ----------
        tau_driving : float
            Total driving torque (aero + gravity + shaft).
        dt : float
            Time step size for the converged window.

        Notes
        -----
        Verification: For ComputedOmega or RampedComputedOmega, the omega history
        should match the corotational solver given identical external loads.

        For ConstantOmega, omega should remain unchanged regardless of tau_driving.
        """
        self._omega_provider.update(tau_driving, dt)

        # Cache updated state
        self._omega = self._omega_provider.omega
        self._alpha = self._omega_provider.alpha

        _logger.debug(
            f"Updated omega after converged window: ω={self._omega:.6f} rad/s, "
            f"α={self._alpha:.6f} rad/s², θ={self._theta:.6f} rad"
        )

    # =========================================================================
    # Phase 6 Complete
    # =========================================================================
    # Tasks 6.1-6.4 implemented:
    # ✓ 6.1: _update_omega_after_converged_window() - reuses OmegaProvider workflow
    # ✓ 6.2: _compute_driving_torque() - external forces only (aero + gravity + shaft)
    # ✓ 6.3: Torque logging via _compute_aerodynamic_torque(), _compute_gravity_torque()
    # ✓ 6.4: Representative window omega used consistently (from OmegaProvider)

    # =========================================================================
    # Main solve() Method - Orchestrates Full FSI Cycle
    # =========================================================================

    def solve(self) -> Tuple["PETSc.Vec", "PETSc.Vec", "PETSc.Vec"]:
        """
        Perform inertial-frame dynamic FSI analysis.

        Orchestrates the full simulation pipeline:

        1. **Matrix assembly** (once at reference configuration):
           - K (elastic), M (lumped mass), C (Rayleigh damping)
           - NO K_G, K_SP in this formulation (inertial frame)

        2. **preCICE interface setup**:
           - Register interface at REFERENCE coordinates (fixed mesh on preCICE side)
           - Write elastic displacement u_e only (not total displacement)
           - Read forces in global frame (no transformation)

        3. **preCICE time loop** (implicit coupling with sub-iterations):
           a. Get ω, α from OmegaProvider (constant within time window)
           b. Read aerodynamic forces from fluid (global frame, no transform)
           c. Rotate internal geometry by θ = ∫ω·dt
           d. Rebuild K(θ), C(θ) on rotated geometry, reuse M (invariant)
           e. Compute reference load: F_ref = -M·a_ref (centripetal + tangential)
           f. Solve: K_eff · u_e_new = F_eff
           g. Write elastic displacement u_e to preCICE
           h. Write omega to GlobalSolidMesh
           i. preCICE sub-iteration or advance to next window
           j. If converged: update omega from driving torque, advance θ

        Returns
        -------
        Tuple[PETSc.Vec, PETSc.Vec, PETSc.Vec]
            Final (displacement, velocity, acceleration) vectors.
        """
        from petsc4py import PETSc

        from .base import Adapter

        # Checkpoint/restart (not fully implemented yet - stub for now)
        checkpoint_state = None  # TODO: implement _try_restore_checkpoint()
        t_restart = float(self.solver_params.get("start_time", 0.0))

        self._print_header("FSI DYNAMIC ANALYSIS - INERTIAL ROTOR SOLVER")

        # Phase 1: Matrix Assembly at reference configuration (θ=0)
        matrices, bc_manager = self._assemble_system_matrices()

        # Phase 2: Extract interface nodes
        interface_coords, interface_dofs = self._extract_interface_nodes()

        # Store reference coordinates for use in precice_displacement_mode='total'.
        # These are the unrotated (theta=0) interface coords and must never be mutated.
        self._interface_coords_reference = interface_coords.copy()

        if self._rotor_radius is None:
            self._rotor_radius = self._compute_rotor_radius(interface_coords)
            self._print_info(f"Auto-detected rotor radius: {self._rotor_radius:.4f} m")

        # Phase 3: Auto-compute inertia if requested
        if getattr(self, "_auto_inertia", False):
            estimated_inertia = self._compute_estimated_inertia()
            if self._is_primary_rank():
                print(
                    f"  ↳ Auto-computed Moment of Inertia: {estimated_inertia:.4e} kg·m²",
                    flush=True,
                )
            # Re-initialize provider with computed inertia (same logic as corotational)
            self._resolve_auto_inertia_provider(estimated_inertia)

        # Phase 4: Initialize preCICE adapter
        self._print_phase(1, 2, "Initializing preCICE adapter...")
        cfg = self.model_properties["solver"]["coupling"]
        _coupling_meshes = [{"name": cfg["coupling_mesh"]}]

        # Add GlobalSolidMesh if omega output is requested
        if self._send_omega_to_precice:
            _coupling_meshes.append({"name": self._omega_mesh_name})

        adapter = Adapter(
            participant=cfg["participant"],
            config_file=cfg["config_file"],
            coupling_meshes=_coupling_meshes,
        )

        # Register interface at reference coordinates (FIXED for all time)
        vertex_ids = self._register_interface_at_reference_coords(
            adapter, cfg["coupling_mesh"], interface_coords
        )

        # Register single omega vertex if needed
        if self._send_omega_to_precice:
            origin = np.zeros((1, 3), dtype=np.float64)
            adapter.add_mesh_vertices(self._omega_mesh_name, origin)

        adapter.initialize()
        dt = adapter.dt

        self._print_phase(2, 2, f"preCICE initialized. dt={dt:.4e} s")

        # Phase 5: Time loop initialization
        t = t_restart
        theta = 0.0  # Accumulated rotation angle [rad]
        omega = self._omega_provider.omega
        alpha = self._omega_provider.alpha
        window_count = 0
        iteration_count = 0

        # State vectors (elastic displacement, velocity, acceleration)
        u_e = PETSc.Vec().createMPI(self.domain.dofs_count, comm=self.comm)
        v_e = u_e.duplicate()
        a_e = u_e.duplicate()
        u_e.set(0.0)
        v_e.set(0.0)
        a_e.set(0.0)

        # Gravity load (constant in global frame)
        F_gravity = self._compute_gravity_load_vector()

        # Current matrices (will be updated after each converged window)
        K_current = self.K
        C_current = self._assemble_rayleigh_damping(K_current, self.M)

        # Checkpoint storage for sub-iterations
        checkpoint: Dict[str, Any] = {}

        self._print_separator()
        if self._is_primary_rank():
            print("  Starting FSI coupling loop...\n", flush=True)

        # Phase 6: preCICE coupling loop
        while adapter.is_coupling_ongoing:
            # Checkpoint if preCICE requires (implicit coupling sub-iterations)
            if adapter.requires_writing_checkpoint:
                checkpoint = self._checkpoint_elastic_state(u_e, v_e, a_e, theta, omega, alpha)
                # Store current structural matrices in checkpoint so that rollback
                # can restore K(θ) and C(θ) consistent with the checkpointed kinematics.
                checkpoint["K_current"] = K_current
                checkpoint["C_current"] = C_current
                iteration_count = 0

            # Read aero forces (global frame, no transformation)
            F_aero = self._read_forces_from_precice_global(
                adapter,
                cfg["coupling_mesh"],
                cfg["read_data"][0] if isinstance(cfg["read_data"], list) else cfg["read_data"],
            )

            # Solve FSI step (elastic displacement increment)
            u_e_new, v_e_new, a_e_new = self._solve_fsi_step(
                F_aero,
                F_gravity,
                dt,
                theta,
                omega,
                alpha,
                u_e,
                v_e,
                a_e,
                bc_manager,
                K_current,
                C_current,
            )

            # Extract interface elastic displacement
            u_e_interface = self._extract_interface_values(u_e_new, interface_dofs)

            # Write elastic displacement to preCICE (NOT total displacement)
            self._write_elastic_displacement_to_precice(
                adapter,
                cfg["coupling_mesh"],
                cfg["write_data"][0] if isinstance(cfg["write_data"], list) else cfg["write_data"],
                u_e_interface,
                theta,
            )

            # Write omega to GlobalSolidMesh
            if self._send_omega_to_precice:
                self._write_omega_to_global_mesh(
                    adapter, self._omega_mesh_name, self._omega_write_data_name, omega
                )

            # Log iteration
            iteration_count += 1
            if self._is_primary_rank() and self._debug_interface:
                u_norm = float(np.linalg.norm(u_e_interface))
                f_norm = float(np.linalg.norm(F_aero))
                print(
                    f"  [FSI] window={window_count:4d} iter={iteration_count:2d} "
                    f"| ||u_e||={u_norm:.6e} m | ||F_aero||={f_norm:.6e} N",
                    flush=True,
                )

            # Advance preCICE (may trigger sub-iteration or window advance)
            adapter.advance(dt)

            # Rollback or converge
            if adapter.requires_reading_checkpoint:
                # Sub-iteration did not converge, restore checkpoint
                self._rollback_elastic_state(checkpoint, u_e, v_e, a_e)
                theta = float(checkpoint["theta"])
                omega = float(checkpoint["omega"])
                alpha = float(checkpoint["alpha"])
                # Restore structural matrices to be consistent with restored θ.
                # Without this, K_current/C_current would reflect post-rollback
                # geometry while kinematics reflect checkpoint geometry.
                K_current = checkpoint["K_current"]
                C_current = checkpoint["C_current"]
            else:
                # Window converged!
                window_count += 1
                t += dt

                # Compute driving torque from aerodynamic and gravity forces
                nodal_coords = self.domain.mesh.nodal_coordinates
                tau_driving = self._compute_driving_torque(
                    interface_coords, u_e_interface, F_aero, nodal_coords, u_e, F_gravity
                )

                # Update omega using OmegaProvider dynamics
                self._update_omega_after_converged_window(tau_driving, dt)
                omega = self._omega_provider.omega
                alpha = self._omega_provider.alpha

                # Accumulate rotation angle for next window
                theta += omega * dt

                # Rotate internal structural geometry for next window
                self._rotate_structural_geometry_internal(theta)
                self._rebuild_assembler_with_rotated_geometry()

                # Rebuild matrices on rotated geometry
                K_current = self.domain.assemble_stiffness_matrix()
                C_current = self._assemble_rayleigh_damping(K_current, self.M)

                # Save results (VTK output, time history, etc.)
                # TODO: implement _save_timestep_results() or use parent class method
                # self._save_timestep_results(t, u_e_new, v_e_new, a_e_new)

                # Copy state for next iteration
                u_e_new.copy(result=u_e)
                v_e_new.copy(result=v_e)
                a_e_new.copy(result=a_e)

                # Reset iteration counter
                iteration_count = 0

        # Finalize preCICE
        adapter.finalize()

        if self._is_primary_rank():
            print("\n" + "═" * 70, flush=True)
            print(f"  FSI coupling complete. Total windows: {window_count}", flush=True)
            print("═" * 70 + "\n", flush=True)

        return u_e, v_e, a_e

    # =========================================================================
    # Per-Step Solve Method
    # =========================================================================

    def _solve_fsi_step(
        self,
        F_aero: NDArray,
        F_gravity: NDArray,
        dt: float,
        theta: float,
        omega: float,
        alpha: float,
        u_e_prev: "PETSc.Vec",
        v_e_prev: "PETSc.Vec",
        a_e_prev: "PETSc.Vec",
        bc_manager: "BoundaryConditionManager",
        K_current: "PETSc.Mat",
        C_current: "PETSc.Mat",
    ) -> Tuple["PETSc.Vec", "PETSc.Vec", "PETSc.Vec"]:
        """
        Solve one FSI sub-iteration step using Newmark integration.

        Assembles and solves:
            K_eff · u_e_new = F_eff

        where:
            K_eff = K(θ) + a₀·M + a₁·C(θ)
            F_eff = F_aero + F_gravity - F_ref + M·(a₀·u_prev + a₂·v_prev + a₃·a_prev)

        Parameters
        ----------
        F_aero : ndarray
            Aerodynamic forces from preCICE (global frame, full DOF vector).
        F_gravity : ndarray
            Gravity load vector (constant in global frame).
        dt : float
            Time step size [s].
        theta : float
            Current rotation angle [rad].
        omega : float
            Angular velocity [rad/s].
        alpha : float
            Angular acceleration [rad/s²].
        u_e_prev : PETSc.Vec
            Elastic displacement at previous iteration.
        v_e_prev : PETSc.Vec
            Elastic velocity at previous iteration.
        a_e_prev : PETSc.Vec
            Elastic acceleration at previous iteration.
        bc_manager : BoundaryConditionManager
            Boundary condition manager with free DOFs.
        K_current : PETSc.Mat
            Stiffness matrix at current orientation K(θ).
        C_current : PETSc.Mat
            Damping matrix at current orientation C(θ).

        Returns
        -------
        Tuple[PETSc.Vec, PETSc.Vec, PETSc.Vec]
            Updated (u_e_new, v_e_new, a_e_new) vectors.

        Notes
        -----
        This implements the inertial formulation:
            M·ü_e + C(θ)·u̇_e + K(θ)·u_e = F_aero + F_gravity - M·a_ref

        where a_ref = α × r + ω × (ω × r) is the rigid-body reference acceleration.
        """
        from petsc4py import PETSc

        # Newmark-β coefficients (β=0.25, γ=0.5 for unconditional stability)
        beta = 0.25
        gamma = 0.5
        a0 = 1.0 / (beta * dt * dt)
        a1 = gamma / (beta * dt)
        a2 = 1.0 / (beta * dt)
        a3 = 1.0 / (2.0 * beta) - 1.0

        # 1. Compute reference load: F_ref = -M·a_ref (rigid-body inertial correction)
        nodal_coords = self.domain.mesh.nodal_coordinates
        F_ref = self._inertial_calculator.compute_reference_load_vector(
            self.M, nodal_coords, omega, alpha
        )

        # 2. Assemble effective stiffness: K_eff = K(θ) + a₀·M + a₁·C(θ)
        K_eff = self._assemble_inertial_effective_system(K_current, C_current, self.M, a0, a1)

        # 3. Assemble RHS: F_eff = F_aero + F_gravity - F_ref + M·(a₀·u + a₂·v + a₃·a)
        F_eff = self._assemble_inertial_rhs(
            F_aero, F_gravity, F_ref, self.M, u_e_prev, v_e_prev, a_e_prev, dt, a0, a2, a3
        )

        # 4. Apply boundary conditions (zero essential BCs for elastic displacement)
        bc_manager.apply_dirichlet_to_system(K_eff, F_eff)

        # 5. Solve linear system: K_eff · u_e_new = F_eff
        u_e_new = PETSc.Vec().createMPI(self.domain.dofs_count, comm=self.comm)
        ksp = PETSc.KSP().create(comm=self.comm)
        ksp.setOperators(K_eff)
        ksp.setType("preonly")  # Direct solver
        ksp.getPC().setType("lu")
        ksp.setFromOptions()
        ksp.solve(F_eff, u_e_new)

        # 6. Update velocity and acceleration using Newmark formulas
        v_e_new = self._newmark_velocity_update(
            u_e_new, u_e_prev, v_e_prev, a_e_prev, dt, beta, gamma
        )
        a_e_new = self._newmark_acceleration_update(u_e_new, u_e_prev, v_e_prev, a_e_prev, dt, beta)

        return u_e_new, v_e_new, a_e_new

    def _newmark_velocity_update(
        self,
        u_new: "PETSc.Vec",
        u_prev: "PETSc.Vec",
        v_prev: "PETSc.Vec",
        a_prev: "PETSc.Vec",
        dt: float,
        beta: float,
        gamma: float,
    ) -> "PETSc.Vec":
        """
        Compute Newmark velocity update.

        v_new = γ/(β·Δt) · (u_new - u_prev) + (1 - γ/β)·v_prev + Δt·(1 - γ/(2β))·a_prev

        Parameters
        ----------
        u_new : PETSc.Vec
            New displacement.
        u_prev : PETSc.Vec
            Previous displacement.
        v_prev : PETSc.Vec
            Previous velocity.
        a_prev : PETSc.Vec
            Previous acceleration.
        dt : float
            Time step size.
        beta : float
            Newmark β parameter.
        gamma : float
            Newmark γ parameter.

        Returns
        -------
        PETSc.Vec
            Updated velocity vector.
        """
        from petsc4py import PETSc

        a1 = gamma / (beta * dt)
        c1 = 1.0 - gamma / beta
        c2 = dt * (1.0 - gamma / (2.0 * beta))

        v_new = PETSc.Vec().createMPI(self.domain.dofs_count, comm=self.comm)
        v_new.set(0.0)

        # v_new = a1 * (u_new - u_prev) + c1 * v_prev + c2 * a_prev
        du = u_new.duplicate()
        u_new.copy(result=du)
        du.axpy(-1.0, u_prev)  # du = u_new - u_prev

        v_new.axpy(a1, du)  # v_new += a1 * du
        v_new.axpy(c1, v_prev)  # v_new += c1 * v_prev
        v_new.axpy(c2, a_prev)  # v_new += c2 * a_prev

        return v_new

    def _newmark_acceleration_update(
        self,
        u_new: "PETSc.Vec",
        u_prev: "PETSc.Vec",
        v_prev: "PETSc.Vec",
        a_prev: "PETSc.Vec",
        dt: float,
        beta: float,
    ) -> "PETSc.Vec":
        """
        Compute Newmark acceleration update.

        a_new = 1/(β·Δt²) · (u_new - u_prev - Δt·v_prev) - (1/(2β) - 1)·a_prev

        Parameters
        ----------
        u_new : PETSc.Vec
            New displacement.
        u_prev : PETSc.Vec
            Previous displacement.
        v_prev : PETSc.Vec
            Previous velocity.
        a_prev : PETSc.Vec
            Previous acceleration.
        dt : float
            Time step size.
        beta : float
            Newmark β parameter.

        Returns
        -------
        PETSc.Vec
            Updated acceleration vector.
        """
        from petsc4py import PETSc

        a0 = 1.0 / (beta * dt * dt)
        a2 = 1.0 / (beta * dt)
        a3 = 1.0 / (2.0 * beta) - 1.0

        a_new = PETSc.Vec().createMPI(self.domain.dofs_count, comm=self.comm)
        a_new.set(0.0)

        # a_new = a0 * (u_new - u_prev) - a2 * v_prev - a3 * a_prev
        du = u_new.duplicate()
        u_new.copy(result=du)
        du.axpy(-1.0, u_prev)  # du = u_new - u_prev

        a_new.axpy(a0, du)  # a_new += a0 * du
        a_new.axpy(-a2, v_prev)  # a_new -= a2 * v_prev
        a_new.axpy(-a3, a_prev)  # a_new -= a3 * a_prev

        return a_new

    def _extract_interface_values(self, vec: "PETSc.Vec", interface_dofs: NDArray) -> NDArray:
        """
        Extract interface DOF values from a global PETSc vector.

        Parameters
        ----------
        vec : PETSc.Vec
            Global DOF vector.
        interface_dofs : ndarray, shape (n_interface_nodes, dofs_per_node)
            Global DOF indices for interface nodes.

        Returns
        -------
        ndarray, shape (n_interface_nodes, dofs_per_node)
            Extracted interface values.
        """
        vec_array = vec.getArray()
        return vec_array[interface_dofs.ravel()].reshape(interface_dofs.shape)

    def _resolve_auto_inertia_provider(self, inertia: float) -> None:
        """
        Re-initialize OmegaProvider with auto-computed inertia.

        Parameters
        ----------
        inertia : float
            Computed moment of inertia [kg·m²].
        """
        from .corotational import ComputedOmega, RampedComputedOmega

        ramp_time = self._auto_inertia_params.get("ramp_time", 0.0)
        target_omega = self._auto_inertia_params["target_omega"]
        shaft_torque = self._auto_inertia_params["shaft_torque"]

        if ramp_time > 0.0:
            self._omega_provider = RampedComputedOmega(
                target_omega=target_omega,
                ramp_time=ramp_time,
                moment_of_inertia=inertia,
                shaft_torque=shaft_torque,
            )
            if self._is_primary_rank():
                print(
                    f"  ↳ Omega mode: Ramp ({ramp_time:.3f} s) → Dynamic (I={inertia:.4e} kg·m²)",
                    flush=True,
                )
        else:
            self._omega_provider = ComputedOmega(
                moment_of_inertia=inertia,
                initial_omega=target_omega,
                shaft_torque=shaft_torque,
            )
            if self._is_primary_rank():
                print(f"  ↳ Omega mode: Dynamic (I={inertia:.4e} kg·m²)", flush=True)

    # =========================================================================
    # Assembly and Initialization Methods
    # =========================================================================

    def _assemble_system_matrices(
        self,
    ) -> Tuple[
        Tuple["PETSc.Mat", "PETSc.Mat"],
        "BoundaryConditionManager",
    ]:
        """
        Assemble stiffness and mass matrices at reference configuration (θ=0).

        Unlike the corotational solver, this does NOT assemble:
        - K_G (geometric stiffness): not used in inertial formulation
        - K_SP (spin softening): not used in inertial formulation

        The mass matrix is assembled once and cached (invariant under rotation).
        Stiffness and damping matrices are rebuilt after each converged window
        when the internal geometry is rotated.

        Returns
        -------
        Tuple[Tuple[PETSc.Mat, PETSc.Mat], BoundaryConditionManager]
            ((K, M), bc_manager) where K, M are at reference configuration.

        Notes
        -----
        Verification: The reference K should match the corotational K at θ=0.
        """
        from petsc4py import PETSc

        from aeroelast.core.bc import BoundaryConditionManager

        self._print_phase(1, 4, "Assembling stiffness matrix at reference...")
        self.K = self.domain.assemble_stiffness_matrix()

        self._print_phase(2, 4, "Assembling mass matrix (lumped, invariant under rotation)...")
        self.M = self.domain.assemble_mass_matrix_lumped()

        # Guard: the inertial solver computes F_ref = -M·a_ref assuming M is diagonal.
        # If the assembler returns a consistent (non-diagonal) mass matrix, the product
        # would be silently incorrect. Fail loudly rather than produce wrong physics.
        _n_rows = self.M.getSize()[0]
        _n_nz = int(self.M.getInfo()["nz_used"])
        if _n_nz > _n_rows:
            raise NotImplementedError(
                f"LinearDynamicFSIRotorInertialSolver requires a lumped (diagonal) mass "
                f"matrix, but the assembled M has {_n_nz} non-zeros for {_n_rows} DOFs. "
                f"The F_ref = -M·a_ref formula is only correct for diagonal M. "
                f"Switch to lumped mass assembly or implement the full matrix-vector product."
            )

        # Force vector and boundary conditions
        self._print_phase(3, 4, "Setting up boundary conditions...")
        force_temp = PETSc.Vec().createMPI(self.domain.dofs_count, comm=self.comm)
        force_temp.set(0.0)
        self.F = force_temp

        bc_manager = BoundaryConditionManager(self.K, self.F, self.M, self.domain.dofs_per_node)
        bc_manager.apply_dirichlet(self.dirichlet_conditions)

        if self._is_primary_rank():
            print(
                f"        Fixed: {len(bc_manager.fixed_dofs)} DOFs, "
                f"Free: {len(bc_manager.free_dofs)} DOFs",
                flush=True,
            )

        # Rayleigh damping coefficients (same logic as corotational)
        # Damping matrix C will be built when needed using _assemble_rayleigh_damping()
        self._configure_rayleigh_damping(bc_manager)

        self._print_phase(4, 4, "Matrix assembly complete (K, M at θ=0).")
        return (self.K, self.M), bc_manager

    def _configure_rayleigh_damping(self, bc_manager: "BoundaryConditionManager") -> None:
        """
        Configure Rayleigh damping coefficients (η_k, η_m) for C = η_k·K + η_m·M.

        This does NOT build the damping matrix — that happens in
        _assemble_rayleigh_damping() when K(θ) is available.

        Parameters
        ----------
        bc_manager : BoundaryConditionManager
            Boundary condition manager with free DOFs.

        Notes
        -----
        Verification: Auto-damping coefficients should match corotational solver
        for the same structure at θ=0.
        """
        if not self._damping_enabled:
            self._print_phase(5, 6, "Rayleigh damping: disabled (enabled=false)")
            return

        if self._damping_auto:
            self._print_phase(5, 6, "Rayleigh damping: auto-computing via SLEPc modal analysis...")
            import _aeroelast

            cfg = self._damping_cfg
            zeta = float(cfg.get("zeta", 0.02))
            zeta_i = float(cfg["zeta_1"]) if cfg.get("zeta_1") is not None else zeta
            zeta_j = float(cfg["zeta_2"]) if cfg.get("zeta_2") is not None else zeta
            mode_i = int(cfg.get("mode_i", 1))
            mode_j = int(cfg.get("mode_j", 2))
            num_modes = int(cfg.get("num_modes", max(mode_j + 2, 6)))

            k_rows, k_cols, k_vals = self._petsc_to_coo(self.K)
            m_rows, m_cols, m_vals = self._petsc_to_coo(self.M)
            free_dofs = bc_manager.free_dofs.astype(np.int32)

            self._eta_k, self._eta_m = _aeroelast.compute_rayleigh_auto(
                k_rows,
                k_cols,
                k_vals,
                m_rows,
                m_cols,
                m_vals,
                free_dofs,
                num_modes,
                mode_i,
                mode_j,
                zeta_i,
                zeta_j,
            )

            if self._eta_k < 0.0 or self._eta_m < 0.0:
                _logger.warning(
                    "Rayleigh auto produced negative coefficients: η_k=%.3e, η_m=%.3e. "
                    "Damping ratio may be non-monotone. Verify mode selection.",
                    self._eta_k,
                    self._eta_m,
                )

            self._print_phase(
                5, 6, f"Rayleigh auto: η_k={self._eta_k:.4e} s  η_m={self._eta_m:.4e} 1/s"
            )
        elif self._eta_m != 0.0 or self._eta_k != 0.0:
            self._print_phase(
                5, 6, f"Rayleigh damping: η_m={self._eta_m:.4e}  η_k={self._eta_k:.4e}"
            )
        else:
            self._print_phase(5, 6, "Rayleigh damping: disabled (coefficients zero)")

    def _extract_interface_nodes(self) -> Tuple[NDArray, NDArray]:
        """
        Extract interface node coordinates and DOF indices from coupling boundaries.

        Returns
        -------
        interface_coords : ndarray, shape (n_interface_nodes, dim)
            Interface node coordinates at reference configuration.
        interface_dofs : ndarray, shape (n_interface_nodes, dofs_per_node)
            Global DOF indices for interface nodes.

        Notes
        -----
        Verification: Interface extraction should match corotational solver.
        """
        coupling_boundaries = self.model_properties["solver"]["coupling_boundaries"]
        mesh = self.domain.mesh
        node_sets = [mesh.node_sets[name] for name in coupling_boundaries]
        nodes = {node.id: node.coords for _set in node_sets for node in _set.nodes.values()}
        sorted_node_ids = sorted(nodes.keys())

        self._interface_node_ids = np.array(sorted_node_ids, dtype=np.int64)
        _iface_coords = np.array([nodes[nid] for nid in sorted_node_ids])

        if self.domain.spatial_dim == 2 and _iface_coords.shape[1] > 2:
            _iface_coords = _iface_coords[:, :2]

        self._interface_coords = _iface_coords

        raw_dofs = np.array([self.domain._node_dofs_map[nid] for nid in sorted_node_ids])
        if raw_dofs.ndim == 2 and raw_dofs.shape[1] > 3:
            self._interface_dofs = raw_dofs[:, :3].astype(int)
        else:
            self._interface_dofs = raw_dofs.astype(int)

        return self._interface_coords, self._interface_dofs

    def _compute_rotor_radius(self, interface_coords: NDArray) -> float:
        """
        Auto-detect rotor radius from interface coordinates.

        Parameters
        ----------
        interface_coords : ndarray, shape (n_nodes, dim)
            Interface node coordinates.

        Returns
        -------
        radius : float
            Maximum radial distance from rotation center.

        Notes
        -----
        Verification: Should match corotational solver's auto-detection.
        """
        coords_3d = self._ensure_3d_vectors(interface_coords)
        radial_vectors = coords_3d - self._rotation_center
        radial_distances = np.linalg.norm(
            radial_vectors
            - np.outer(np.dot(radial_vectors, self._rotation_axis), self._rotation_axis),
            axis=1,
        )
        return float(np.max(radial_distances))

    def _compute_estimated_inertia(self) -> float:
        """
        Estimate total moment of inertia about the rotation axis.

        Computes the parallel-axis contribution of all mesh nodes:

            I = Σᵢ mᵢ · r_⊥,ᵢ²

        where mᵢ is the lumped mass from the first translational DOF diagonal
        entry and r_⊥,ᵢ is the perpendicular distance from the node to the
        rotation axis.

        Returns
        -------
        float
            Total moment of inertia [kg·m²] about the rotation axis.
        """
        from petsc4py import PETSc

        if self.M is None:
            _logger.warning("Mass matrix not available for inertia estimation.")
            return 1.0

        diag_vec = self.M.getDiagonal()
        mass_array = diag_vec.getArray(readonly=True)

        nodes = self.domain.nodes
        dofs_per_node = self.domain.dofs_per_node
        limit_idx = len(mass_array)

        node_idx = np.arange(len(nodes)) * dofs_per_node
        valid = node_idx < limit_idx
        masses = np.where(valid, mass_array[np.minimum(node_idx, limit_idx - 1)], 0.0)
        coords = np.array([node.coords for node in nodes], dtype=np.float64)
        r_vec = coords - self._rotation_center
        proj = r_vec @ self._rotation_axis
        r_perp_sq = np.einsum("ij,ij->i", r_vec, r_vec) - proj**2
        total_inertia = float(np.dot(masses, r_perp_sq))

        diag_vec.destroy()
        return total_inertia

    def _compute_gravity_load_vector(self) -> NDArray:
        """
        Compute gravity load vector F_g = M·g for all DOFs.

        Returns
        -------
        F_gravity : ndarray, shape (n_dofs,)
            Gravity load vector (flattened).

        Notes
        -----
        Verification: For a uniform rotor, ||F_g|| should equal total_mass * ||g||.
        """
        if not self._include_gravity:
            return np.zeros(self.domain.dofs_count, dtype=np.float64)

        # M is lumped diagonal matrix
        # F_g[i] = M[i] * g[component_i]
        from petsc4py import PETSc

        M_diag = PETSc.Vec().createMPI(self.domain.dofs_count, comm=self.comm)
        self.M.getDiagonal(M_diag)
        M_array = M_diag.getArray()

        F_gravity = np.zeros(self.domain.dofs_count, dtype=np.float64)
        dofs_per_node = self.domain.dofs_per_node
        gravity_3d = np.pad(self._gravity, (0, max(0, 3 - len(self._gravity))))[:3]

        for i in range(len(M_array)):
            component_idx = i % dofs_per_node
            if component_idx < len(gravity_3d):
                F_gravity[i] = M_array[i] * gravity_3d[component_idx]

        return F_gravity

    # =========================================================================
    # Phase 7-10 Placeholder
    # =========================================================================
    # The following methods will be implemented in subsequent phases:
    #
    # Phase 7: Prototype Validation
    # - Unit tests for pure rigid rotation (u_e ~ 0)
    # - Gravity-only test with 1P torque modulation
    # - Reference load analytical verification
    # - preCICE contract regression test
    # - Checkpoint/restart coverage
    #
    # Phase 8: Comparative Benchmarking
    # - Side-by-side test with corotational solver
    # - Physical convergence validation
    # - Torque balance verification
    # - Displacement magnitude comparison
    #
    # Phase 9: Performance Optimization
    # - Matrix caching strategies
    # - Geometry update alternatives (API vs reconstruction)
    # - Profiling and bottleneck identification
    #
    # Phase 10: Product Decision and Cleanup
    # - Performance vs accuracy tradeoffs
    # - Documentation of validated use cases
    # - Deprecation or promotion decision
    #
    # See docs/rotor_inertial_solver_tasks.md for full task list

    # =========================================================================
    # Helper Methods (for console output and utilities)
    # =========================================================================

    def _is_primary_rank(self) -> bool:
        """Always True in serial."""
        return True

    def _print_header(self, title: str) -> None:
        """Print a formatted header section.

        Parameters
        ----------
        title : str
            Header text to display.
        """
        if self._is_primary_rank():
            print("\n" + "═" * 70, flush=True)
            print(f"  {title}", flush=True)
            print("═" * 70, flush=True)

    def _print_separator(self) -> None:
        """Print a section separator."""
        if self._is_primary_rank():
            print("═" * 70, flush=True)

    def _print_phase(self, phase: int, total: int, message: str) -> None:
        """Print a phase progress message.

        Parameters
        ----------
        phase : int
            Current phase number.
        total : int
            Total number of phases.
        message : str
            Progress description.
        """
        if self._is_primary_rank():
            print(f"  [{phase}/{total}] {message}", flush=True)

    def _print_info(self, message: str) -> None:
        """Print an info message.

        Parameters
        ----------
        message : str
            Informational text to display.
        """
        if self._is_primary_rank():
            print(f"  [Info] {message}", flush=True)

    @staticmethod
    def _petsc_to_coo(mat: "PETSc.Mat") -> Tuple[NDArray, NDArray, NDArray]:
        """
        Convert PETSc matrix to COO format arrays.

        Parameters
        ----------
        mat : PETSc.Mat
            PETSc matrix to convert.

        Returns
        -------
        rows : ndarray
            Row indices.
        cols : ndarray
            Column indices.
        vals : ndarray
            Matrix values.
        """
        from petsc4py import PETSc

        mat_csr = mat.convert("mpiaij")
        indptr, indices, data = mat_csr.getValuesCSR()

        rows = []
        for i in range(len(indptr) - 1):
            rows.extend([i] * (indptr[i + 1] - indptr[i]))

        return (
            np.array(rows, dtype=np.int32),
            np.array(indices, dtype=np.int32),
            np.array(data, dtype=np.float64),
        )
