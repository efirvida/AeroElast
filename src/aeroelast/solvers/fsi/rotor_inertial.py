"""
Inertial Frame Rotor FSI Solver — Python configuration wrapper.

The FSI solve loop is implemented in Rust:
    crates/aeroelast-solvers/src/petsc/fsi/rotor_inertial.rs
    (InertialRotorFsiSolver)

This module provides LinearDynamicFSIRotorInertialSolver, a thin Python
configuration and assembly wrapper that:
  1. Parses solver parameters from fem_model_properties.
  2. Assembles K, M, C matrices via PETSc and converts them to COO format.
  3. Extracts interface nodes and computes rotor geometry (radius, inertia).
  4. Delegates the full preCICE coupling loop to the Rust solver via
     ``_aeroelast.run_inertial_rotor_fsi_solver(...)``.

Physical formulation (solved in Rust)
--------------------------------------
The equation of motion in the INERTIAL frame is:

    [M]{ü_e} + [C(θ)]{u̇_e} + [K(θ)]{u_e} = {F_aero} + {F_g} - [M]{a_ref}

where:
- u_e: elastic displacement in global coordinates (NOT total displacement)
- θ: rotation angle (kinematics only, no frame rotation)
- K(θ), C(θ): assembled on the rigidly-rotated structural geometry
- a_ref = α × r + ω × (ω × r): rigid-body reference acceleration

preCICE interface contract:
- SolidMesh vertices at REFERENCE coordinates (t=0, no rotation)
- Write data: elastic displacement u_e in global frame
- Read data: aerodynamic forces in global frame (no transformation)
- GlobalSolidMesh: representative angular velocity

See docs/rotor_inertial_solver_design.md for full design documentation.
"""

import logging
import os
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
from .rotor import LinearDynamicFSIRotorCorotationalSolver

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
        - include_geometric_stiffness: SUPPORTED — enables dynamic K_G(θ, ω)
            assembly in Rust on the rotated geometry used for K(θ).
    - include_spin_softening: SUPPORTED — K_SP diagonal assembled in Rust and
      updated when ω changes (see reassemble_k_if_needed).
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

        YAML keys for the new spin-softening and K_G ω-rebuild features:

        .. code-block:: yaml

            solver:
              rotor:
                # Spin-softening K_SP = -ω²·M·(I - n̂⊗n̂).
                # Default: true. Set false only for debugging (reduces accuracy).
                include_spin_softening: true

                # Geometric stiffness (centrifugal prestress) K_G(θ, ω).
                # Default: false.
                include_geometric_stiffness: false

                # Per-matrix Δω hysteresis thresholds.
                # Rebuild K_SP when |Δ(ω²)|/ω² exceeds ksp_omega_rebuild_high
                # (not recently rebuilt) or ksp_omega_rebuild_low (recently rebuilt).
                ksp_omega_rebuild_high: 0.005   # 0.5% change in ω²
                ksp_omega_rebuild_low: 0.003    # 0.3% change in ω²
                kg_omega_rebuild_high: 0.005
                kg_omega_rebuild_low: 0.003
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

        # Physics flags
        # include_geometric_stiffness: enables dynamic K_G(θ, ω) assembly in
        # the Rust inertial runtime. The wrapper only forwards the flag; the
        # actual K_G is assembled on the current rotated geometry in Rust.
        # include_spin_softening: K_SP = -ω²·M·(I - n̂⊗n̂) updated in Rust when ω
        # changes beyond the configured hysteresis thresholds. Default: True
        # (matches HEAD unconditional behavior — do NOT change the default).
        self._include_geometric_stiffness: bool = bool(
            rotor_cfg.get("include_geometric_stiffness", False)
        )
        self._include_spin_softening: bool = bool(
            rotor_cfg.get("include_spin_softening", True)
        )
        # Per-matrix Δω hysteresis thresholds.
        # K_SP thresholds: rebuild K_SP when |Δ(ω²)|/ω² exceeds these.
        _ksp_high = rotor_cfg.get("ksp_omega_rebuild_high", 0.005)
        _ksp_low = rotor_cfg.get("ksp_omega_rebuild_low", 0.003)
        if not isinstance(_ksp_high, (int, float)):
            raise ValueError(
                f"ksp_omega_rebuild_high must be a float, got {type(_ksp_high).__name__!r}"
            )
        if not isinstance(_ksp_low, (int, float)):
            raise ValueError(
                f"ksp_omega_rebuild_low must be a float, got {type(_ksp_low).__name__!r}"
            )
        self._ksp_omega_rebuild_high: float = float(_ksp_high)
        self._ksp_omega_rebuild_low: float = float(_ksp_low)
        # K_G thresholds: rebuild K_G when |Δ(ω²)|/ω² exceeds these.
        _kg_high = rotor_cfg.get("kg_omega_rebuild_high", 0.005)
        _kg_low = rotor_cfg.get("kg_omega_rebuild_low", 0.003)
        if not isinstance(_kg_high, (int, float)):
            raise ValueError(
                f"kg_omega_rebuild_high must be a float, got {type(_kg_high).__name__!r}"
            )
        if not isinstance(_kg_low, (int, float)):
            raise ValueError(
                f"kg_omega_rebuild_low must be a float, got {type(_kg_low).__name__!r}"
            )
        self._kg_omega_rebuild_high: float = float(_kg_high)
        self._kg_omega_rebuild_low: float = float(_kg_low)

        if self._include_geometric_stiffness:
            _logger.info(
                "include_geometric_stiffness=True: enabling dynamic K_G(theta, omega) "
                "assembly in the inertial Rust runtime.",
            )
        if not self._include_spin_softening:
            _logger.warning(
                "include_spin_softening=False: K_SP is disabled. This reduces physical "
                "accuracy at operating speed (flapwise eigenfrequencies over-estimated)."
            )
        else:
            _logger.debug(
                "include_spin_softening=True: K_SP assembled in Rust "
                "(reassemble_k_if_needed). Δω thresholds: ksp_high=%.4f ksp_low=%.4f, "
                "kg_high=%.4f kg_low=%.4f",
                self._ksp_omega_rebuild_high,
                self._ksp_omega_rebuild_low,
                self._kg_omega_rebuild_high,
                self._kg_omega_rebuild_low,
            )
        if not rotor_cfg.get("include_centrifugal", True):
            _logger.warning(
                "include_centrifugal=False has no effect in the inertial solver. "
                "Centrifugal acceleration is implicitly included via -M·a_ref and "
                "cannot be disabled independently."
            )
        if rotor_cfg.get("include_coriolis", False):
            _logger.warning("(no Coriolis in inertial frame)")
        if rotor_cfg.get("include_euler", False):
            _logger.warning(
                "include_euler is set but NOT used in inertial solver (Euler term is part of a_ref)"
            )

        # Force ramp and magnitude cap (same as corotational)
        self._force_ramp_time = float(rotor_cfg.get("force_ramp_time", 0.0))
        _fmax = rotor_cfg.get("force_max_magnitude", None)
        self._force_max_magnitude: Optional[float] = float(_fmax) if _fmax is not None else None

        # K(θ) reassembly cadence. Each rebuild forces a KSP refactorization
        # (LU/Cholesky), so cadence is a real performance knob. Three triggers
        # are OR-combined inside the Rust solver:
        #   - k_update_interval: forced rebuild every N converged windows
        #   - omega_rebuild_threshold: |Δω²|/ω² (only useful for ramped/computed ω)
        #   - theta_rebuild_threshold: |Δθ| since last rebuild (rad) — physically
        #     correct primary trigger (K depends on θ, not ω)
        # Defaults amortize cost while keeping geometry error bounded.
        self._k_update_interval = int(rotor_cfg.get("k_update_interval", 20))
        if self._k_update_interval < 1:
            self._k_update_interval = 1
        self._omega_rebuild_threshold = float(rotor_cfg.get("omega_rebuild_threshold", 0.01))
        self._theta_rebuild_threshold = float(rotor_cfg.get("theta_rebuild_threshold", 0.05))

        # K_G(u) — Nivel 2 deformed-coords K_G assembly. When enabled, the
        # Rust runtime reassembles K_G at X_rotated + u whenever the radial
        # deflection ratio drift exceeds the high-band threshold (with the
        # low-band reserved for hysteresis). Defaults preserve the legacy
        # explicit-in-geometry baseline.
        self._kg_use_deformed_coords: bool = bool(
            rotor_cfg.get("kg_use_deformed_coords", False)
        )
        self._kg_deflection_rebuild_rel_high = float(
            rotor_cfg.get("kg_deflection_rebuild_rel_high", 0.01)
        )
        self._kg_deflection_rebuild_rel_low = float(
            rotor_cfg.get("kg_deflection_rebuild_rel_low", 0.005)
        )
        if self._kg_use_deformed_coords:
            _logger.info(
                "kg_use_deformed_coords=True: K_G will be reassembled at "
                "X_rotated + u (Nivel 2 foreshortening). High/low rebuild "
                "thresholds = %.4f / %.4f",
                self._kg_deflection_rebuild_rel_high,
                self._kg_deflection_rebuild_rel_low,
            )

        # Omega output to preCICE (same as corotational)
        self._send_omega_to_precice = rotor_cfg.get("send_omega_to_precice", True)
        self._omega_mesh_name: str = rotor_cfg.get("omega_mesh_name", "GlobalSolidMesh")
        self._omega_write_data_name: str = rotor_cfg.get("omega_write_data", "AngularVelocity")

        # Nodal velocity output to preCICE. In Rust this follows the same
        # kinematic contract as displacement_mode: elastic-only in "elastic"
        # mode, or elastic + rigid-body velocity in "total" mode.
        self._send_velocity_to_precice: bool = bool(
            rotor_cfg.get("send_velocity_to_precice", False)
        )
        self._velocity_write_data_name: str = rotor_cfg.get("velocity_write_data", "Velocity")

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
        transform_disp = rotor_cfg.get("transform_displacement_to_inertial")
        if transform_disp is True:
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

        petsc_cfg = self.solver_params.get("petsc") or {}
        self._petsc_factor_options_prefix = str(
            petsc_cfg.get("factor_options_prefix", "rotor_inertial_lu_")
        )
        self._petsc_factor_reuse_ordering = bool(petsc_cfg.get("factor_reuse_ordering", True))
        self._petsc_factor_reuse_fill = bool(petsc_cfg.get("factor_reuse_fill", True))
        self._petsc_factor_mat_ordering_type = petsc_cfg.get("factor_mat_ordering_type")
        self._petsc_factorization_type = (
            str(petsc_cfg.get("factorization_type", "cholesky")).strip().lower()
        )
        if self._petsc_factorization_type not in {"cholesky", "lu"}:
            _logger.warning(
                "Unknown inertial PETSc factorization_type=%r; falling back to 'cholesky'",
                self._petsc_factorization_type,
            )
            self._petsc_factorization_type = "cholesky"

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

        # Reference interface coordinates (set in solve(); used when
        # displacement_mode='total' is configured in the Rust solver).
        self._interface_coords_reference: Optional["NDArray"] = None

        # All-node masses (flat, one scalar per node) — extracted during assembly
        # before BCs reduce the matrix size. Set in _assemble_system_matrices().
        self._all_node_masses_full: Optional["NDArray"] = None

        # Lumped mass diagonal is invariant under rigid rotation; cache it once
        # to avoid recreating PETSc diagonal vectors on every FSI sub-iteration.
        self._mass_diagonal_array: Optional[NDArray] = None

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

        return self._solve_via_rust(
            bc_manager=bc_manager,
            interface_coords=interface_coords,
        )

    # =========================================================================
    # Rust fast-path delegate
    # =========================================================================

    def _solve_via_rust(
        self,
        bc_manager: "BoundaryConditionManager",
        interface_coords: "NDArray",
    ) -> Tuple["PETSc.Vec", "PETSc.Vec", "PETSc.Vec"]:
        """Delegate the FSI loop to the compiled Rust solver.

        Mirrors the pattern used by
        ``LinearDynamicFSIRotorCorotationalSolver._solve_via_rust``.  The
        Python solver remains the authoritative reference / fallback.

        Parameters
        ----------
        bc_manager :
            Boundary condition manager with ``free_dofs`` populated.
        interface_coords :
            (n_iface_nodes × 3) array of interface node coordinates in the
            reference (θ=0) frame.
        """
        import _aeroelast  # noqa: PLC0415
        import numpy as np  # noqa: PLC0415

        from petsc4py import PETSc  # noqa: PLC0415

        rust_asm = getattr(self.domain, "_rust", None)
        if rust_asm is None:
            raise RuntimeError(
                "Rust assembler (domain._rust) is not available; "
                "cannot delegate to run_inertial_rotor_fsi_solver."
            )

        k_rows, k_cols, k_vals = self._petsc_to_coo(self.K)
        m_rows, m_cols, m_vals = self._petsc_to_coo(self.M)

        free_dofs = bc_manager.free_dofs.astype(np.int32)
        n_total: int = int(self.K.getSize()[0])
        _free_dofs_arr = np.asarray(free_dofs, dtype=np.int32)

        # Fixed DOF values (Dirichlet)
        _fixed_dof_vals: dict = dict(bc_manager.fixed_dofs)

        # All-node masses (flat, one scalar per node) — extracted before BC application
        # to include masses for ALL nodes (both free and fixed).
        all_node_masses_np = self._all_node_masses_full

        # Interface node indices (0-based) — same sorted node IDs as _extract_interface_nodes
        node_id_to_idx = self.domain.mesh.node_id_to_index
        iface_node_indices = np.array(
            [node_id_to_idx[nid] for nid in self._interface_node_ids],
            dtype=np.uintp,
        )

        # Newmark parameters
        beta = float(self.solver_params.get("beta", 0.25))
        gamma = float(self.solver_params.get("gamma", 0.5))
        dt_hint = float(self.solver_params.get("dt", 0.0))

        # preCICE config
        cfg = self.model_properties["solver"]["coupling"]
        mesh_name = cfg["coupling_mesh"]
        write_data = (
            cfg["write_data"] if isinstance(cfg["write_data"], str) else cfg["write_data"][0]
        )
        read_data = cfg["read_data"] if isinstance(cfg["read_data"], str) else cfg["read_data"][0]

        # Omega provider → Rust-compatible params
        from .corotational import (  # noqa: PLC0415
            ConstantOmega,
            RampedOmega,
            ComputedOmega,
            RampedComputedOmega,
        )

        p = self._omega_provider
        if isinstance(p, ConstantOmega):
            omega_mode, omega_val, omega_target, t_ramp, moi, shaft_tau = (
                "constant",
                float(p._omega),
                None,
                None,
                None,
                None,
            )
        elif isinstance(p, RampedOmega):
            omega_mode, omega_val, omega_target, t_ramp, moi, shaft_tau = (
                "ramped",
                0.0,
                float(p._target_omega),
                float(p._ramp_time),
                None,
                None,
            )
        elif isinstance(p, ComputedOmega):
            omega_mode, omega_val, omega_target, t_ramp, moi, shaft_tau = (
                "computed",
                float(p._omega),
                None,
                None,
                float(p._I),
                float(p._tau_shaft),
            )
        elif isinstance(p, RampedComputedOmega):
            omega_mode, omega_val, omega_target, t_ramp, moi, shaft_tau = (
                "ramped_computed",
                0.0,
                float(p._target_omega),
                float(p._ramp_time),
                float(p._I),
                float(p._tau_shaft),
            )
        else:
            omega_val_0, _ = p.get_omega(0.0)
            omega_mode, omega_val, omega_target, t_ramp, moi, shaft_tau = (
                "constant",
                float(omega_val_0),
                None,
                None,
                None,
                None,
            )

        # Omega GlobalSolidMesh (optional)
        _omega_mesh = self._omega_mesh_name if self._send_omega_to_precice else None
        _omega_data = self._omega_write_data_name if self._send_omega_to_precice else None
        _omega_coord = list(self._coord_transforms.center) if self._send_omega_to_precice else None

        if self._is_primary_rank():
            print(
                f"  [Rust] Delegating to run_inertial_rotor_fsi_solver "
                f"(mode={omega_mode}, disp={self._precice_displacement_mode})",
                flush=True,
            )

        # ── Geometric stiffness: use explicit kwarg (deprecate sentinel) ─────
        # Pass include_geometric_stiffness explicitly; keep kg0_* as None so
        # callers still using the legacy sentinel pattern are not broken by the
        # presence of the shim in the Rust binding.
        kg0_rows = kg0_cols = kg0_vals = None  # sentinel kept as None (deprecated path)

        all_node_coords_nodes = self._ensure_3d_vectors(
            np.array([n.coords for n in self.domain.nodes], dtype=np.float64)
        )
        interface_coords_nodes = self._ensure_3d_vectors(np.asarray(interface_coords))
        step_callback = self._build_output_step_callback(
            n_total=n_total,
            free_dofs=free_dofs,
            fixed_dof_vals=_fixed_dof_vals,
            all_node_coords_nodes=all_node_coords_nodes,
            interface_coords_nodes=interface_coords_nodes,
            interface_node_indices=iface_node_indices,
        )

        u_final_red, v_final_red, a_final_red, times = _aeroelast.run_inertial_rotor_fsi_solver(
            assembler=rust_asm,
            rotation_axis=list(self._coord_transforms.axis),
            rotation_center=list(self._coord_transforms.center),
            all_node_masses=all_node_masses_np,
            omega_mode=omega_mode,
            omega=omega_val,
            omega_target=omega_target,
            t_ramp=t_ramp,
            moment_of_inertia=moi,
            shaft_torque=shaft_tau,
            gravity=list(self._gravity),
            include_reference_acceleration=True,
            k_update_interval=self._k_update_interval,
            omega_rebuild_threshold=self._omega_rebuild_threshold,
            theta_rebuild_threshold=self._theta_rebuild_threshold,
            kg_use_deformed_coords=self._kg_use_deformed_coords,
            kg_deflection_rebuild_rel_high=self._kg_deflection_rebuild_rel_high,
            kg_deflection_rebuild_rel_low=self._kg_deflection_rebuild_rel_low,
            displacement_mode=self._precice_displacement_mode,
            dofs_per_node=self.domain.dofs_per_node,
            fluid_density=self._fluid_density,
            flow_velocity=self._flow_velocity,
            rotor_radius=float(self._rotor_radius),
            k_rows=k_rows,
            k_cols=k_cols,
            k_vals=k_vals,
            m_rows=m_rows,
            m_cols=m_cols,
            m_vals=m_vals,
            free_dofs=free_dofs,
            eta_k=self._eta_k,
            eta_m=self._eta_m,
            beta=beta,
            gamma=gamma,
            dt=dt_hint,
            interface_nodes=iface_node_indices,
            mesh_dims=self.domain.spatial_dim,
            participant_name=cfg["participant"],
            config_file=cfg["config_file"],
            coupling_mesh=mesh_name,
            write_data_name=write_data,
            read_data_name=read_data,
            ramp_time=float(self._force_ramp_time),
            force_max=self._force_max_magnitude,
            omega_mesh_name=_omega_mesh,
            omega_write_data=_omega_data,
            omega_vertex_coord=_omega_coord,
            velocity_write_data=(
                self._velocity_write_data_name if self._send_velocity_to_precice else None
            ),
            u0=None,
            v0=None,
            a0=None,
            t0=0.0,
            theta0=0.0,
            restart_omega=None,
            restart_alpha=None,
            restart_ramp_completed=None,
            restart_current_time=None,
            kg0_rows=kg0_rows,
            kg0_cols=kg0_cols,
            kg0_vals=kg0_vals,
            step_callback=step_callback,
            include_geometric_stiffness=self._include_geometric_stiffness,
            include_spin_softening=self._include_spin_softening,
            ksp_omega_rebuild_high=self._ksp_omega_rebuild_high,
            ksp_omega_rebuild_low=self._ksp_omega_rebuild_low,
            kg_omega_rebuild_high=self._kg_omega_rebuild_high,
            kg_omega_rebuild_low=self._kg_omega_rebuild_low,
        )

        n_steps = len(times)
        if n_steps > 0 and self._is_primary_rank():
            print(
                f"  ✓ Inertial FSI loop complete: "
                f"{n_steps} converged steps, t_final={times[-1]:.4f} s",
                flush=True,
            )
        elif self._is_primary_rank():
            print("  ⚠️ Inertial FSI loop returned 0 converged steps.", flush=True)

        checkpoint_manager = getattr(self, "_checkpoint_manager", None)
        if checkpoint_manager is not None:
            checkpoint_manager.finalize(timeout=60.0)

        # Expand reduced DOF arrays back to full-DOF PETSc vectors
        u_final_full = np.zeros(n_total, dtype=np.float64)
        v_final_full = np.zeros(n_total, dtype=np.float64)
        a_final_full = np.zeros(n_total, dtype=np.float64)

        if n_steps > 0:
            u_final_full[_free_dofs_arr] = np.asarray(u_final_red, dtype=np.float64)
            for dof, val in _fixed_dof_vals.items():
                u_final_full[dof] = val
            v_final_full[_free_dofs_arr] = np.asarray(v_final_red, dtype=np.float64)
            a_final_full[_free_dofs_arr] = np.asarray(a_final_red, dtype=np.float64)

        u_vec = PETSc.Vec().createWithArray(u_final_full, comm=self.comm)
        v_vec = PETSc.Vec().createWithArray(v_final_full, comm=self.comm)
        a_vec = PETSc.Vec().createWithArray(a_final_full, comm=self.comm)

        self.u = u_vec
        self.v = v_vec
        self.a = a_vec

        return u_vec, v_vec, a_vec

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

        # Extract all-node masses BEFORE applying BCs (for F_ref = -M·a_ref computation)
        # The diagonal of M has size n_dofs_total; we take the x-DOF mass for each node.
        _m_diag_full = self.M.getDiagonal()
        _m_diag_full_arr = _m_diag_full.getArray(readonly=True).copy()
        _m_diag_full.destroy()
        self._all_node_masses_full = _m_diag_full_arr.reshape(-1, self.domain.dofs_per_node)[
            :, 0
        ].copy()

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

            _logger.warning(
                "Rayleigh auto coefficients computed at θ=0 (reference geometry). "
                "They are applied as C(θ)=η_m·M+η_k·K(θ) at all orientations. "
                "For composite shells, K(θ) eigenvalues may differ from K(0) when "
                "material axes are not aligned with the rotation axis — effective ζ "
                "may deviate from the design value at large rotation angles."
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

    def _compute_rotor_radius(
        self,
        interface_coords: NDArray,
        interface_disps: Optional[NDArray] = None,
    ) -> float:
        """Compute rotor radius as max perpendicular distance from rotation axis."""
        coords_3d = self._ensure_3d_vectors(interface_coords)
        if interface_disps is not None:
            coords_3d = coords_3d + self._ensure_3d_vectors(interface_disps)
        radial_vectors = coords_3d - self._rotation_center
        radial_distances = np.linalg.norm(
            radial_vectors
            - np.outer(np.dot(radial_vectors, self._rotation_axis), self._rotation_axis),
            axis=1,
        )
        return float(np.max(radial_distances))

    _write_restart_state = LinearDynamicFSIRotorCorotationalSolver._write_restart_state
    _log_rotor_performance = LinearDynamicFSIRotorCorotationalSolver._log_rotor_performance
    _build_omega_checkpoint_state = (
        LinearDynamicFSIRotorCorotationalSolver._build_omega_checkpoint_state
    )
    _extract_nodal_translation_field = (
        LinearDynamicFSIRotorCorotationalSolver._extract_nodal_translation_field
    )
    _scatter_nodal_translation_field = (
        LinearDynamicFSIRotorCorotationalSolver._scatter_nodal_translation_field
    )
    _compute_performance_coefficients = (
        LinearDynamicFSIRotorCorotationalSolver._compute_performance_coefficients
    )

    def _compute_global_axis_torque(
        self,
        nodal_coords: NDArray,
        nodal_disps: NDArray,
        nodal_forces: NDArray,
    ) -> Tuple[NDArray, float]:
        """Compute torque directly in the inertial frame."""
        coords_3d = self._ensure_3d_vectors(nodal_coords)
        disps_3d = self._ensure_3d_vectors(nodal_disps)
        forces_3d = self._ensure_3d_vectors(nodal_forces)

        rel_pos = coords_3d + disps_3d - self._coord_transforms.center
        torque_global = np.sum(np.cross(rel_pos, forces_3d), axis=0)
        torque_scalar = float(np.dot(torque_global, self._coord_transforms.axis))
        return torque_global, torque_scalar

    def _build_output_step_callback(
        self,
        *,
        n_total: int,
        free_dofs: NDArray,
        fixed_dof_vals: Dict[int, float],
        all_node_coords_nodes: NDArray,
        interface_coords_nodes: NDArray,
        interface_node_indices: NDArray,
    ):
        """Build the per-window callback that writes CSV, probes, and checkpoints."""
        free_dofs_arr = np.asarray(free_dofs, dtype=np.int64)
        all_node_masses = np.asarray(self._all_node_masses_full, dtype=np.float64)

        def _step_cb(
            t,
            time_step,
            dt,
            u_red,
            v_red,
            a_red,
            force_mag,
            forces_iface,
            omega,
            alpha,
            theta,
            rotor_perf_tuple,
        ):
            checkpoint_manager = getattr(self, "_checkpoint_manager", None)
            tau_aero_rust = float(rotor_perf_tuple[0]) if len(rotor_perf_tuple) > 0 else 0.0
            omega_window = float(omega)
            alpha_window = float(alpha)
            self._theta = float(theta)

            u_full = np.zeros(n_total, dtype=np.float64)
            v_full = np.zeros(n_total, dtype=np.float64)
            a_full = np.zeros(n_total, dtype=np.float64)

            u_full[free_dofs_arr] = np.asarray(u_red, dtype=np.float64)
            v_full[free_dofs_arr] = np.asarray(v_red, dtype=np.float64)
            a_full[free_dofs_arr] = np.asarray(a_red, dtype=np.float64)
            for dof, val in fixed_dof_vals.items():
                u_full[dof] = val

            u_nodes_global = self._extract_nodal_translation_field(u_full)

            needs_stress = (
                self._stress_output_interval <= 1
                or (time_step % self._stress_output_interval == 0)
                or (checkpoint_manager is not None and checkpoint_manager.should_write(t))
                or bool(getattr(self, "_probe_node_ids", []))
            )
            stress_fields = self._compute_stress_fields(u_full) if needs_stress else {}

            coords_rotated = self._coord_transforms.rotate_point_cloud(all_node_coords_nodes, theta)
            iface_coords_rotated = self._coord_transforms.rotate_point_cloud(
                interface_coords_nodes, theta
            )

            iface_force_global = np.zeros((len(interface_node_indices), 3), dtype=np.float64)
            if forces_iface is not None and len(forces_iface) > 0:
                iface_force_global = self._ensure_3d_vectors(
                    np.asarray(forces_iface, dtype=np.float64).reshape(-1, self.domain.spatial_dim)
                )

            iface_u_global = u_nodes_global[np.asarray(interface_node_indices, dtype=np.int64)]

            f_aero_nodes_global = np.zeros_like(all_node_coords_nodes)
            for local_idx, node_idx in enumerate(interface_node_indices):
                f_aero_nodes_global[int(node_idx)] += iface_force_global[local_idx]

            a_ref_global = self._inertial_calculator.compute_rigid_body_acceleration_inertial(
                coords_rotated,
                omega_window,
                alpha_window,
            )
            f_inertial_nodes_global = -all_node_masses[:, np.newaxis] * a_ref_global

            f_gravity_nodes_global = np.zeros_like(all_node_coords_nodes)
            if self._include_gravity:
                f_gravity_nodes_global = all_node_masses[:, np.newaxis] * self._gravity

            f_total_nodes_global = (
                f_aero_nodes_global + f_inertial_nodes_global + f_gravity_nodes_global
            )

            tau_aero_global, tau_aero = self._compute_global_axis_torque(
                iface_coords_rotated,
                iface_u_global,
                iface_force_global,
            )
            if abs(tau_aero) <= 1.0e-14 and abs(tau_aero_rust) > 1.0e-14:
                tau_aero = tau_aero_rust
                tau_aero_global = self._coord_transforms.axis * tau_aero_rust

            _, tau_inertial = self._compute_global_axis_torque(
                coords_rotated,
                u_nodes_global,
                f_inertial_nodes_global,
            )
            _, tau_gravity = self._compute_global_axis_torque(
                coords_rotated,
                u_nodes_global,
                f_gravity_nodes_global,
            )
            tau_total_global, tau_total = self._compute_global_axis_torque(
                coords_rotated,
                u_nodes_global,
                f_total_nodes_global,
            )

            thrust = float(np.dot(np.sum(iface_force_global, axis=0), self._coord_transforms.axis))
            deformed_radius = self._compute_rotor_radius(iface_coords_rotated, iface_u_global)
            ct, cp, cq, tsr = self._compute_performance_coefficients(
                thrust,
                tau_aero,
                omega_window,
                deformed_radius,
            )

            torque_non_aero = tau_total - tau_aero
            power_aero = tau_aero * omega_window
            power_total = tau_total * omega_window
            structural_efficiency = (
                float(np.clip(-torque_non_aero / tau_aero, 0.0, 1.0))
                if abs(tau_aero) > 1.0e-14
                else 0.0
            )
            max_displacement = (
                float(np.max(np.linalg.norm(u_nodes_global, axis=1)))
                if len(u_nodes_global) > 0
                else 0.0
            )

            force_fields = {
                "F_AERO_RAW": self._scatter_nodal_translation_field(f_aero_nodes_global, n_total),
                "F_AERO": self._scatter_nodal_translation_field(f_aero_nodes_global, n_total),
                "F_INERTIAL": self._scatter_nodal_translation_field(
                    f_inertial_nodes_global,
                    n_total,
                ),
                "F_GRAVITY": self._scatter_nodal_translation_field(
                    f_gravity_nodes_global,
                    n_total,
                ),
                "F_TOTAL": self._scatter_nodal_translation_field(f_total_nodes_global, n_total),
                "OMEGA": omega_window,
                "ALPHA": alpha_window,
                "OMEGA_STATE": omega,
                "ALPHA_STATE": alpha,
                "THETA": theta,
                "TAU_AERO": tau_aero,
                "TAU_INERTIAL": tau_inertial,
                "TAU_GRAVITY": tau_gravity,
                "TAU_TOTAL": tau_total,
                "THRUST": thrust,
                "DEFORMED_RADIUS": deformed_radius,
                "CT": ct,
                "CP": cp,
                "CQ": cq,
                "TSR": tsr,
            }
            force_fields.update(stress_fields)

            self._log_rotor_performance(
                t=t,
                omega_rpm=omega_window * 60.0 / (2.0 * np.pi),
                omega_rad=omega_window,
                alpha=alpha_window,
                angle_deg=np.degrees(theta),
                thrust=thrust,
                torque_aero=tau_aero,
                torque_non_aero=torque_non_aero,
                torque_inertial=tau_inertial,
                torque_gravity=tau_gravity,
                torque_total=tau_total,
                power_aero=power_aero,
                power_total=power_total,
                structural_efficiency=structural_efficiency,
                cp=cp,
                cq=cq,
                ct=ct,
                tsr=tsr,
                torque_aero_global=tau_aero_global,
                torque_total_global=tau_total_global,
                max_displacement=max_displacement,
                deformed_radius=deformed_radius,
            )
            self._write_restart_state(
                t=t,
                theta=theta,
                omega=omega,
                alpha=alpha,
            )
            self._log_structural_report(
                t=t,
                time_step=time_step,
                u_full=u_full,
                v_full=v_full,
                a_full=a_full,
                stress_fields=force_fields,
                applied_force_mag=force_mag,
            )
            self._log_probe_data(
                t=t,
                time_step=time_step,
                u_full=u_full,
                v_full=v_full,
                stress_fields=force_fields,
            )
            checkpoint_kwargs = {
                "theta": theta,
                **self._build_omega_checkpoint_state(t=t, omega=omega, alpha=alpha),
            }
            self._handle_checkpoint(
                t=t,
                time_step=time_step,
                dt=dt,
                u_red=np.asarray(u_red, dtype=np.float64),
                v_red=np.asarray(v_red, dtype=np.float64),
                a_red=np.asarray(a_red, dtype=np.float64),
                u_full=u_full,
                v_full=v_full,
                a_full=a_full,
                extra_fields=force_fields,
                **checkpoint_kwargs,
            )

        return _step_cb

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

    def _get_mass_diagonal_array(self) -> NDArray:
        """Return the lumped mass diagonal as a cached NumPy array."""
        if self._mass_diagonal_array is None:
            diag_vec = self.M.getDiagonal()
            self._mass_diagonal_array = diag_vec.getArray(readonly=True).copy()
            diag_vec.destroy()
        return self._mass_diagonal_array

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
