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
_DEFAULT_ROTATION_AXIS = np.array([0.0, 0.0, 1.0])
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
            rotation_axis: [0, 0, 1] # Z-axis
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
                "include_euler is set but NOT used in inertial solver "
                "(Euler term is part of a_ref)"
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
            perf_cfg.get("flow_velocity") or rotor_cfg.get("flow_velocity") or _DEFAULT_FLOW_VELOCITY
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
    # Future Implementation Methods (Phase 3-10)
    # =========================================================================
    # The following methods will be implemented in subsequent phases:
    #
    # Phase 3: Structural Assembly on Internally Rotated Geometry
    # - _rotate_structural_geometry_internal()
    # - _assemble_stiffness_on_rotated_geometry()
    # - _assemble_damping_from_rayleigh()
    #
    # Phase 4: Newmark Step for the Inertial Formulation
    # - _assemble_inertial_effective_system()
    # - _assemble_inertial_rhs()
    #
    # Phase 5: preCICE Contract with Fixed Interface Mesh
    # - _write_elastic_displacement_to_precice()
    # - _read_forces_from_precice()
    #
    # Phase 6: Omega Dynamics and Torque Accounting
    # - _update_omega_from_converged_window()
    # - _compute_aerodynamic_torque()
    #
    # Phase 7-10: Validation, benchmarking, optimization
    # See docs/rotor_inertial_solver_tasks.md for full task list
