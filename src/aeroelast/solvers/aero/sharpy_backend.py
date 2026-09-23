"""SHARPy-backed aerodynamic backend for rotor FSI participants.

The backend plugs SHARPy into the existing ``RotorAerodynamicBackend`` contract.
It intentionally keeps the SHARPy case/data construction behind an adapter:
SHARPy owns a rich ``PreSharpy`` data model, while aeroelast owns the preCICE
participant loop and full-rotor nodal ordering.  The adapter is the translation
boundary between those two worlds.
"""

from __future__ import annotations

import importlib
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping, Protocol

from .sharpy_fsi_adapter import build_quasi_static_sharpy_fsi_adapter
from .types import AeroFSIRuntimeContext, RotorAeroLoads, RotorAeroState

logger = logging.getLogger(__name__)


class SharpyCaseAdapter(Protocol):
    """Translation layer between aeroelast runtime state and SHARPy data."""

    def build_data(
        self,
        context: AeroFSIRuntimeContext,
        state: RotorAeroState,
        settings: Mapping[str, Any],
    ) -> Any:
        """Create the SHARPy data object used to initialise a solver."""

    def extract_loads(
        self,
        data: Any,
        context: AeroFSIRuntimeContext,
        state: RotorAeroState,
    ) -> RotorAeroLoads:
        """Map SHARPy aerodynamic forces into aeroelast full-rotor nodal loads."""


@dataclass(frozen=True)
class _SharpySolverSpec:
    model: str
    module: str
    class_name: str
    unsteady: bool


_SOLVER_SPECS: dict[str, _SharpySolverSpec] = {
    "vlm": _SharpySolverSpec(
        model="vlm",
        module="sharpy.solvers.staticuvlm",
        class_name="StaticUvlm",
        unsteady=False,
    ),
    "uvlm": _SharpySolverSpec(
        model="uvlm",
        module="sharpy.solvers.stepuvlm",
        class_name="StepUvlm",
        unsteady=True,
    ),
}


_VLM_SOLVER_SETTING_KEYS = {
    "centre_rot_g",
    "cfl1",
    "horseshoe",
    "ignore_first_x_nodes_in_force_calculation",
    "iterative_precond",
    "iterative_solver",
    "iterative_tol",
    "map_forces_on_struct",
    "n_rollup",
    "nonlifting_body_interactions",
    "num_cores",
    "only_nonlifting",
    "phantom_wing_test",
    "print_info",
    "rbm_vel_g",
    "rho",
    "rollup_aic_refresh",
    "rollup_dt",
    "rollup_tolerance",
    "velocity_field_generator",
    "velocity_field_input",
    "vortex_radius",
    "vortex_radius_wake_ind",
}

_UVLM_SOLVER_SETTING_KEYS = {
    "centre_rot",
    "centre_rot_g",
    "cfl1",
    "convection_scheme",
    "dt",
    "filter_method",
    "gamma_dot_filtering",
    "ignore_first_x_nodes_in_force_calculation",
    "interp_coords",
    "interp_method",
    "iterative_precond",
    "iterative_solver",
    "iterative_tol",
    "n_time_steps",
    "nonlifting_body_interactions",
    "num_cores",
    "only_nonlifting",
    "phantom_wing_test",
    "print_info",
    "quasi_steady",
    "rho",
    "velocity_field_generator",
    "velocity_field_input",
    "vortex_radius",
    "vortex_radius_wake_ind",
    "yaw_slerp",
}


def _normalize_model(model: str) -> str:
    normalized = str(model).strip().lower()
    if normalized not in _SOLVER_SPECS:
        supported = ", ".join(sorted(_SOLVER_SPECS))
        raise ValueError(f"Unsupported SHARPy aerodynamic model '{model}'. Supported: {supported}")
    return normalized


def _load_dotted_object(path: str) -> Any:
    module_name, separator, attr_name = path.replace(":", ".").rpartition(".")
    if not separator or not module_name or not attr_name:
        raise ValueError(
            "SHARPy case_adapter must be a dotted object path such as "
            "'my_package.my_module.MyAdapter'"
        )

    module = importlib.import_module(module_name)
    return getattr(module, attr_name)


def _build_case_adapter(adapter: Any) -> SharpyCaseAdapter | None:
    if adapter is None:
        return None
    if isinstance(adapter, str):
        adapter_obj = _load_dotted_object(adapter)
        return adapter_obj() if isinstance(adapter_obj, type) else adapter_obj
    return adapter


class SharpyAeroBackend:
    """Aerodynamic backend that delegates VLM/UVLM solves to SHARPy."""

    name = "sharpy"

    def __init__(
        self,
        *,
        model: str = "uvlm",
        solver_settings: Mapping[str, Any] | None = None,
        case_adapter: SharpyCaseAdapter | str | None = None,
        rebuild_per_step: bool = False,
    ) -> None:
        self._model = _normalize_model(model)
        self._solver_spec = _SOLVER_SPECS[self._model]
        self._solver_settings = dict(solver_settings or {})
        self._case_adapter = _build_case_adapter(case_adapter)
        self._rebuild_per_step = bool(rebuild_per_step)
        self._solver = None
        self._data = None
        self._compute_calls = 0

    @property
    def model(self) -> str:
        return self._model

    @property
    def solver_settings(self) -> dict[str, Any]:
        return dict(self._solver_settings)

    def _load_solver_class(self):
        try:
            module = importlib.import_module(self._solver_spec.module)
        except ModuleNotFoundError as exc:
            raise ModuleNotFoundError(
                "SHARPy is required for aero.backend='sharpy'. Install the optional "
                "dependency with 'pip install -e .[sharpy]' after loading the GCC module."
            ) from exc
        return getattr(module, self._solver_spec.class_name)

    def _settings_for_state(self, state: RotorAeroState) -> dict[str, Any]:
        settings = dict(self._solver_settings)
        if self._solver_spec.unsteady:
            if state.dt > 0.0:
                settings.setdefault("dt", float(state.dt))
            settings.setdefault("n_time_steps", 1)
        return settings

    def _ensure_adapter(self) -> SharpyCaseAdapter:
        if self._case_adapter is None:
            raise RuntimeError(
                "SHARPy backend requires a case_adapter to translate aeroelast rotor "
                "state into SHARPy data and map SHARPy forces back to nodal loads."
            )
        return self._case_adapter

    def _ensure_solver(
        self,
        context: AeroFSIRuntimeContext,
        state: RotorAeroState,
    ) -> None:
        if self._solver is not None:
            return

        adapter = self._ensure_adapter()
        solver_cls = self._load_solver_class()
        logger.info(
            "[AERO-SHARPY] Initializing solver model=%s rebuild_per_step=%s",
            self._model,
            self._rebuild_per_step,
        )
        logger.info("[AERO-SHARPY] Building SHARPy data model...")
        self._data = adapter.build_data(context, state, self._settings_for_state(state))
        build_initialise_settings = getattr(adapter, "build_initialise_settings", None)
        if callable(build_initialise_settings):
            init_settings = dict(build_initialise_settings(self._data, context, state) or {})
        else:
            init_settings = self._settings_for_state(state)
        logger.info(
            "[AERO-SHARPY] Creating solver class=%s with %d init settings keys",
            self._solver_spec.class_name,
            len(init_settings),
        )
        self._solver = solver_cls()
        self._solver.initialise(
            self._data,
            custom_settings=init_settings,
            restart=False,
        )
        self._data = self._solver.data
        logger.info("[AERO-SHARPY] Solver initialization completed.")

    def _update_adapter_state(
        self,
        context: AeroFSIRuntimeContext,
        state: RotorAeroState,
    ) -> None:
        adapter = self._ensure_adapter()
        update_state = getattr(adapter, "update_state", None)
        if callable(update_state):
            update_state(self._data, context, state)

    def _run_kwargs(
        self,
        context: AeroFSIRuntimeContext,
        state: RotorAeroState,
    ) -> dict[str, Any]:
        adapter = self._ensure_adapter()
        build_run_kwargs = getattr(adapter, "build_run_kwargs", None)
        if callable(build_run_kwargs):
            kwargs = dict(build_run_kwargs(self._data, context, state) or {})
        else:
            kwargs = {}

        if self._solver_spec.unsteady:
            if state.dt > 0.0:
                kwargs.setdefault("dt", float(state.dt))
            kwargs.setdefault("t", float(state.time))
            kwargs.setdefault("convect_wake", True)
        return kwargs

    def compute_loads(
        self,
        context: AeroFSIRuntimeContext,
        state: RotorAeroState,
    ) -> RotorAeroLoads:
        adapter = self._ensure_adapter()
        self._compute_calls += 1
        if self._compute_calls == 1:
            logger.info(
                "[AERO-SHARPY] First compute_loads call: t=%.6f dt=%.6f omega=%.6f",
                float(state.time),
                float(state.dt),
                float(state.omega_rad_s),
            )
        if self._rebuild_per_step:
            self._solver = None
            self._data = None
        self._ensure_solver(context, state)
        self._update_adapter_state(context, state)
        if self._compute_calls == 1:
            logger.info("[AERO-SHARPY] Running SHARPy aerodynamic solve (first call)...")
        self._data = self._solver.run(**self._run_kwargs(context, state))
        if self._compute_calls == 1:
            logger.info("[AERO-SHARPY] First SHARPy solve completed, extracting loads.")
        return adapter.extract_loads(self._data, context, state)

    def capture_state(self) -> dict[str, Any]:
        adapter = self._case_adapter
        if self._data is None:
            return {"data": None}
        capture_state = getattr(adapter, "capture_state", None)
        if callable(capture_state):
            return {"data": capture_state(self._data)}
        return {"data": None}

    def restore_state(self, checkpoint_state: Any) -> None:
        adapter = self._case_adapter
        restore_state = getattr(adapter, "restore_state", None)
        if callable(restore_state):
            data_state = (
                checkpoint_state.get("data") if isinstance(checkpoint_state, dict) else None
            )
            if self._data is None:
                if data_state is not None:
                    logger.warning(
                        "[AERO-SHARPY] Received checkpoint data before solver initialization; "
                        "state restore deferred."
                    )
                return
            restore_state(self._data, data_state)

    def finalize_time_window(self) -> None:
        adapter = self._case_adapter
        finalize = getattr(adapter, "finalize_time_window", None)
        if callable(finalize) and self._data is not None:
            finalize(self._data)

    def prime_aero_mesh_export(
        self,
        context: AeroFSIRuntimeContext,
        state: RotorAeroState,
        output_file: str | Path,
    ) -> None:
        self._ensure_solver(context, state)
        self.export_aero_mesh(output_file)

    def export_aero_mesh(self, output_file: str | Path) -> None:
        adapter = self._ensure_adapter()
        export_mesh = getattr(adapter, "export_aero_mesh", None)
        if not callable(export_mesh):
            raise RuntimeError("SHARPy case adapter does not support aerodynamic mesh export")
        if self._data is None:
            raise RuntimeError(
                "SHARPy aerodynamic mesh is not available before solver initialization"
            )
        export_mesh(self._data, output_file)


def build_sharpy_backend(runtime_context: AeroFSIRuntimeContext) -> SharpyAeroBackend:
    """Build a SHARPy backend from normalized aerodynamic config."""

    aero_cfg = runtime_context.normalized_config.get("aero", {})
    backend_cfg_raw = aero_cfg.get("backend_config", {}) if isinstance(aero_cfg, dict) else {}
    legacy_cfg_raw = runtime_context.normalized_config.get("sharpy", {})
    backend_cfg = dict(legacy_cfg_raw) if isinstance(legacy_cfg_raw, dict) else {}
    if isinstance(backend_cfg_raw, dict):
        backend_cfg.update(backend_cfg_raw)

    model = str(backend_cfg.get("model", backend_cfg.get("method", "uvlm")))
    normalized_model = _normalize_model(model)
    allowed_solver_keys = (
        _VLM_SOLVER_SETTING_KEYS if normalized_model == "vlm" else _UVLM_SOLVER_SETTING_KEYS
    )

    raw_settings = dict(backend_cfg.get("solver_settings", {}))
    settings = {key: value for key, value in raw_settings.items() if key in allowed_solver_keys}
    air_density = backend_cfg.get("air_density")
    if air_density is not None:
        settings.setdefault("rho", float(air_density))

    # Allow adapter params in either backend_config root or solver_settings for
    # backward compatibility with earlier SHARPy benchmark configs.
    adapter_cfg: dict[str, Any] = dict(backend_cfg)
    adapter_cfg.update(raw_settings)

    case_adapter = backend_cfg.get("case_adapter")
    if case_adapter is None and model.strip().lower() == "vlm":
        # V1 default: quasi-static coupled adapter for SHARPy VLM.
        case_adapter = build_quasi_static_sharpy_fsi_adapter(runtime_context, adapter_cfg)

    logger.info(
        "[AERO-SHARPY] Backend config: model=%s rebuild_per_step=%s solver_keys=%d",
        model,
        bool(backend_cfg.get("rebuild_per_step", False)),
        len(settings),
    )

    return SharpyAeroBackend(
        model=model,
        solver_settings=settings,
        case_adapter=case_adapter,
        rebuild_per_step=bool(backend_cfg.get("rebuild_per_step", False)),
    )
