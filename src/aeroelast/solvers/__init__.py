from aeroelast.postprocess.stress_recovery import (
    StrainResult,
    StressLocation,
    StressRecovery,
    StressResult,
    StressType,
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
] + _linear_names
