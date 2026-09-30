from aeroelast.postprocess.stress_recovery import (
    StrainResult,
    StressLocation,
    StressRecovery,
    StressResult,
    StressType,
)

from .aero import (
    SUPPORTED_AERO_FSI_BACKENDS,
    AeroFSIRuntimeContext,
    BladeAeroSurface,
    BladeAeroSurfaceState,
    FullRotorAeroGeometry,
    FullRotorAeroSurface,
    FullRotorAeroSurfaceState,
    QuasiSteadyVLMBackend,
    RotorAerodynamicBackend,
    RotorAeroFSIParticipant,
    RotorAeroLoads,
    RotorAeroState,
    RotorBladeAeroGeometry,
    build_aero_participant_from_config,
    build_aero_runtime_context,
    build_blade_aero_surface,
    build_full_rotor_aero_geometry,
    build_full_rotor_aero_surface,
    build_rotor_aero_participant,
    build_vlm_backend,
    discover_rotor_blade_set_names,
    normalize_aero_fsi_config,
    resolve_aero_backend,
)

_linear_names = []
try:
    # Import PETSc-backed structural solvers only when the dependency is
    # available, so importing the solvers package does not break CLI preview/view.
    from .elasticity import DynamicNewmarkSolver, StaticLinearSolver, StaticNonlinearSolver
    from .modal import ModalSolver as ModalSolver

    # Legacy aliases kept for backward compatibility
    LinearStaticSolver = StaticLinearSolver
    LinearDynamicSolver = DynamicNewmarkSolver
    NonlinearStaticSolver = StaticNonlinearSolver

    _linear_names = [
        "StaticLinearSolver",
        "StaticNonlinearSolver",
        "DynamicNewmarkSolver",
        "LinearStaticSolver",
        "LinearDynamicSolver",
        "NonlinearStaticSolver",
        "ModalSolver",
    ]
except ImportError:
    # PETSc not available
    pass

try:
    # FSI classes remain optional because preCICE is not part of the minimal
    # installation used for config inspection and mesh preprocessing.
    from .fsi import (
        Adapter as Adapter,
        ConstantOmega as ConstantOmega,
        CoordinateTransforms as CoordinateTransforms,
        ForceClipper as ForceClipper,
        FSIRunner as FSIRunner,
        FunctionOmega as FunctionOmega,
        InertialForcesCalculator as InertialForcesCalculator,
        LinearDynamicFSIRotorSolver as LinearDynamicFSIRotorSolver,
        LinearDynamicFSISolver as LinearDynamicFSISolver,
        OmegaProvider as OmegaProvider,
        TableOmega as TableOmega,
        run_from_yaml as run_from_yaml,
    )
except ImportError:
    # preCICE not available
    pass

__all__ = [
    "StrainResult",
    "StressLocation",
    "StressRecovery",
    "StressResult",
    "StressType",
    "AeroFSIRuntimeContext",
    "BladeAeroSurface",
    "BladeAeroSurfaceState",
    "FullRotorAeroGeometry",
    "FullRotorAeroSurface",
    "FullRotorAeroSurfaceState",
    "QuasiSteadyVLMBackend",
    "RotorAeroFSIParticipant",
    "RotorAeroLoads",
    "RotorAerodynamicBackend",
    "RotorAeroState",
    "RotorBladeAeroGeometry",
    "SUPPORTED_AERO_FSI_BACKENDS",
    "build_aero_runtime_context",
    "build_aero_participant_from_config",
    "build_blade_aero_surface",
    "build_full_rotor_aero_geometry",
    "build_full_rotor_aero_surface",
    "build_rotor_aero_participant",
    "build_vlm_backend",
    "discover_rotor_blade_set_names",
    "normalize_aero_fsi_config",
    "resolve_aero_backend",
] + _linear_names
