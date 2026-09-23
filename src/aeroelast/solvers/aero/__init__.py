"""Aerodynamic-solver foundations for external rotor participants.

This package contains the geometry and typing primitives used by future
low- and medium-fidelity aerodynamic participants that couple to the
structural solver through preCICE.
"""

from .participant import (
    SUPPORTED_AERO_FSI_BACKENDS,
    build_aero_participant_from_config,
    build_aero_runtime_context,
    normalize_aero_fsi_config,
    resolve_aero_backend,
)
from .report_schema import (
    DEFAULT_AERO_REPORT_SCHEMA_MANIFEST_FILENAME,
    AeroReportDefinition,
    AeroReportSchemaManifest,
    load_aero_report_schema,
    load_aero_report_table,
)
from .rotor_geometry import (
    FullRotorAeroGeometry,
    RotorBladeAeroGeometry,
    build_full_rotor_aero_geometry,
    discover_rotor_blade_set_names,
)
from .runtime_participant import RotorAeroFSIParticipant, build_rotor_aero_participant
from .sharpy_backend import SharpyAeroBackend, SharpyCaseAdapter, build_sharpy_backend
from .sharpy_fsi_adapter import (
    SharpyRotorFSICaseAdapter,
    build_quasi_static_sharpy_fsi_adapter,
)
from .sharpy_numad import NuMADSharpyRotorCaseAdapter, run_sharpy_numad_rotor_case
from .surface import (
    BladeAeroSurface,
    BladeAeroSurfaceState,
    FullRotorAeroSurface,
    FullRotorAeroSurfaceState,
    build_blade_aero_surface,
    build_full_rotor_aero_surface,
)
from .types import AeroFSIRuntimeContext, RotorAerodynamicBackend, RotorAeroLoads, RotorAeroState
from .validation import (
    IEA15MW_REFERENCE_OPERATING_POINTS,
    IEA15MW_ZERO_PITCH_TSR_SWEEP,
    RotorModelComparison,
    RotorOperatingPoint,
    RotorPerformanceMetrics,
    RotorSectionalDiagnostics,
    RotorValidationObjectiveMetrics,
    build_reduced_vlm_rotor_mesh,
    compare_sharpy_to_bem_sweep,
    compare_vlm_experimental_gated_lift_to_bem_sweep,
    compare_vlm_experimental_gated_tsr_lift_to_bem_sweep,
    compare_vlm_experimental_lift_to_bem_sweep,
    compare_vlm_polar_corrected_to_bem_sweep,
    compare_vlm_to_bem_sweep,
    comparison_as_dict,
    compute_bem_metrics,
    compute_validation_objective_metrics,
    compute_vlm_polar_corrected_reduced_rotor_metrics,
    compute_vlm_reduced_rotor_metrics,
    compute_vlm_sectional_diagnostics,
)
from .vlm import QuasiSteadyVLMBackend, build_vlm_backend

__all__ = [
    "AeroFSIRuntimeContext",
    "AeroReportDefinition",
    "AeroReportSchemaManifest",
    "BladeAeroSurface",
    "BladeAeroSurfaceState",
    "DEFAULT_AERO_REPORT_SCHEMA_MANIFEST_FILENAME",
    "FullRotorAeroGeometry",
    "FullRotorAeroSurface",
    "FullRotorAeroSurfaceState",
    "IEA15MW_REFERENCE_OPERATING_POINTS",
    "IEA15MW_ZERO_PITCH_TSR_SWEEP",
    "QuasiSteadyVLMBackend",
    "NuMADSharpyRotorCaseAdapter",
    "RotorModelComparison",
    "RotorAeroFSIParticipant",
    "RotorAeroLoads",
    "RotorAerodynamicBackend",
    "RotorOperatingPoint",
    "RotorPerformanceMetrics",
    "RotorSectionalDiagnostics",
    "RotorValidationObjectiveMetrics",
    "RotorAeroState",
    "RotorBladeAeroGeometry",
    "SharpyAeroBackend",
    "SharpyCaseAdapter",
    "SharpyRotorFSICaseAdapter",
    "SUPPORTED_AERO_FSI_BACKENDS",
    "build_aero_runtime_context",
    "build_aero_participant_from_config",
    "build_blade_aero_surface",
    "build_reduced_vlm_rotor_mesh",
    "build_full_rotor_aero_geometry",
    "build_full_rotor_aero_surface",
    "build_rotor_aero_participant",
    "build_quasi_static_sharpy_fsi_adapter",
    "build_sharpy_backend",
    "build_vlm_backend",
    "compare_vlm_experimental_gated_lift_to_bem_sweep",
    "compare_vlm_experimental_gated_tsr_lift_to_bem_sweep",
    "compare_vlm_experimental_lift_to_bem_sweep",
    "compare_sharpy_to_bem_sweep",
    "compare_vlm_polar_corrected_to_bem_sweep",
    "compare_vlm_to_bem_sweep",
    "comparison_as_dict",
    "compute_bem_metrics",
    "compute_vlm_polar_corrected_reduced_rotor_metrics",
    "compute_vlm_reduced_rotor_metrics",
    "compute_vlm_sectional_diagnostics",
    "compute_validation_objective_metrics",
    "discover_rotor_blade_set_names",
    "load_aero_report_schema",
    "load_aero_report_table",
    "normalize_aero_fsi_config",
    "resolve_aero_backend",
    "run_sharpy_numad_rotor_case",
]
