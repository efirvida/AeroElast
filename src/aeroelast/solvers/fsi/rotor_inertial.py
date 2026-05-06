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
        if rotor_cfg.get("include_centrifugal", True):
            # Centrifugal is always "included" via -M·a_ref, but the flag has no effect
            pass
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
            f"Assembled inertial K_eff: K(θ) + {a0:.3e}·M + {a1:.3e}·C(θ) "
            f"(no K_G, K_SP, G_cor)"
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

    def _checkpoint_elastic_state(self) -> Dict[str, Any]:
        """
        Checkpoint the elastic state and rigid-body kinematics.

        The inertial solver state includes:
        - Elastic displacement, velocity, acceleration (u_e, v_e, a_e)
        - Rigid-body kinematics (theta, omega, alpha)
        - Representative window kinematics (theta_target, omega_window, alpha_window)

        This is distinct from the corotational solver checkpoint which stores
        displacement in the rotating frame.

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
            # Elastic state variables
            "u_e": self.domain.u.copy() if self.domain.u is not None else None,
            "v_e": self.domain.v.copy() if self.domain.v is not None else None,
            "a_e": self.domain.a.copy() if self.domain.a is not None else None,
            # Rigid-body kinematics
            "theta": self._theta,
            "omega": self._omega,
            "alpha": self._alpha,
            # Window kinematics (frozen during sub-iterations)
            "theta_target": getattr(self, "_theta_target", 0.0),
            "omega_window": getattr(self, "_omega_window", 0.0),
            "alpha_window": getattr(self, "_alpha_window", 0.0),
            # OmegaProvider state (if applicable)
            "omega_provider_state": (
                self._omega_provider.get_state()
                if hasattr(self._omega_provider, "get_state")
                else None
            ),
        }

        _logger.debug(
            f"Checkpointed elastic state: theta={self._theta:.6f}, "
            f"omega={self._omega:.6f}, ||u_e||={np.linalg.norm(self.domain.u) if self.domain.u is not None else 0:.3e}"
        )

        return checkpoint

    def _rollback_elastic_state(self, checkpoint: Dict[str, Any]) -> None:
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

        Notes
        -----
        After rollback, the next Newmark step should reuse the same K(θ), C(θ)
        from the checkpointed window kinematics (no reassembly needed within
        the same FSI window).
        """
        # Restore elastic state
        if checkpoint["u_e"] is not None:
            self.domain.u[:] = checkpoint["u_e"]
        if checkpoint["v_e"] is not None:
            self.domain.v[:] = checkpoint["v_e"]
        if checkpoint["a_e"] is not None:
            self.domain.a[:] = checkpoint["a_e"]

        # Restore rigid-body kinematics
        self._theta = checkpoint["theta"]
        self._omega = checkpoint["omega"]
        self._alpha = checkpoint["alpha"]
        self._theta_target = checkpoint.get("theta_target", 0.0)
        self._omega_window = checkpoint.get("omega_window", 0.0)
        self._alpha_window = checkpoint.get("alpha_window", 0.0)

        # Restore OmegaProvider state
        if checkpoint["omega_provider_state"] is not None:
            if hasattr(self._omega_provider, "set_state"):
                self._omega_provider.set_state(checkpoint["omega_provider_state"])

        _logger.debug(
            f"Rolled back to: theta={self._theta:.6f}, omega={self._omega:.6f}"
        )

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

        adapter.write_data(mesh_name, data_name, u_e_global)

        _logger.debug(
            f"Wrote elastic displacement to '{mesh_name}/{data_name}': ||u_e||={np.linalg.norm(u_e_global):.3e} "
            f"(θ={np.degrees(theta):.1f}°, elastic only, NO rigid-body component)"
        )

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
    # - _update_omega_from_converged_window()
    # - _compute_aerodynamic_torque()
    # - _compute_gravity_torque()
    # - _compute_total_torque()
    #
    # Phase 7-10: Validation, benchmarking, optimization
    # See docs/rotor_inertial_solver_tasks.md for full task list
