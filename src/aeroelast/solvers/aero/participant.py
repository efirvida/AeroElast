"""Generic aerodynamic participant dispatch.

This layer normalizes new backend-agnostic aerodynamic participant configs
while preserving compatibility with the existing BEM-specific YAML schema.
"""

from __future__ import annotations

import copy
import logging
from typing import Any, Mapping

from aeroelast.core.mesh.model import MeshModel
from aeroelast.solvers.aero.rotor_geometry import (
    ROTOR_BLADE_SET_PREFIX,
    build_full_rotor_aero_geometry,
)
from aeroelast.solvers.aero.surface import build_full_rotor_aero_surface

from .runtime_participant import build_rotor_aero_participant
from .types import AeroFSIRuntimeContext

SUPPORTED_AERO_FSI_BACKENDS = ("bem", "vlm", "sharpy")
logger = logging.getLogger(__name__)


def _extract_rotor_config(normalized_cfg: Mapping[str, Any], backend: str) -> dict[str, Any]:
    aero_cfg = normalized_cfg.get("aero")
    rotor_cfg_raw = aero_cfg.get("rotor") if isinstance(aero_cfg, Mapping) else None
    rotor_cfg = dict(rotor_cfg_raw) if isinstance(rotor_cfg_raw, Mapping) else {}

    if backend == "bem":
        bem_cfg = normalized_cfg.get("bem")
        if isinstance(bem_cfg, Mapping):
            if "n_blades" in bem_cfg and "n_blades" not in rotor_cfg:
                rotor_cfg["n_blades"] = int(bem_cfg["n_blades"])
            if "hub_radius" in bem_cfg and "hub_radius" not in rotor_cfg:
                rotor_cfg["hub_radius"] = float(bem_cfg["hub_radius"])

    rotor_cfg.setdefault("rotation_axis", [0.0, 1.0, 0.0])
    rotor_cfg.setdefault("rotation_center", [0.0, 0.0, 0.0])
    return rotor_cfg


def resolve_aero_backend(
    cfg: Mapping[str, Any],
    default_backend: str | None = None,
) -> str:
    """Resolve the aerodynamic backend from new or legacy config layouts."""
    backend = default_backend
    if backend is None:
        aero_cfg = cfg.get("aero")
        if isinstance(aero_cfg, Mapping):
            backend = aero_cfg.get("backend")
    if backend is None and "bem" in cfg:
        backend = "bem"

    if backend is None:
        raise ValueError(
            "Unable to determine the aerodynamic backend. "
            "Provide aero.backend in the YAML or a legacy top-level bem section."
        )

    backend = str(backend).strip().lower()
    if backend not in SUPPORTED_AERO_FSI_BACKENDS:
        raise NotImplementedError(
            f"Aerodynamic backend '{backend}' is not implemented yet. "
            f"Supported backends: {', '.join(SUPPORTED_AERO_FSI_BACKENDS)}"
        )
    return backend


def normalize_aero_fsi_config(
    cfg: Mapping[str, Any],
    default_backend: str | None = None,
) -> dict[str, Any]:
    """Normalize a generic aerodynamic participant config into legacy backend input."""
    normalized = copy.deepcopy(dict(cfg))
    backend = resolve_aero_backend(normalized, default_backend=default_backend)

    aero_cfg_raw = normalized.get("aero")
    aero_cfg = dict(aero_cfg_raw) if isinstance(aero_cfg_raw, Mapping) else {}
    aero_cfg.setdefault("backend", backend)
    normalized["aero"] = aero_cfg

    blade_file = aero_cfg.get("blade_file")
    if "blade_file" not in normalized and blade_file is not None:
        normalized["blade_file"] = blade_file

    output_cfg = aero_cfg.get("output")
    if "output" not in normalized and isinstance(output_cfg, Mapping):
        normalized["output"] = dict(output_cfg)

    rotor_cfg_raw = aero_cfg.get("rotor")
    rotor_cfg = dict(rotor_cfg_raw) if isinstance(rotor_cfg_raw, Mapping) else {}
    mesh_params = normalized.get("mesh", {}).get("generator", {}).get("params")
    if isinstance(mesh_params, dict):
        if "n_blades" not in mesh_params and "n_blades" in rotor_cfg:
            mesh_params["n_blades"] = rotor_cfg["n_blades"]
        if "hub_radius" not in mesh_params and "hub_radius" in rotor_cfg:
            mesh_params["hub_radius"] = rotor_cfg["hub_radius"]

    if backend != "bem":
        legacy_backend_raw = normalized.get(backend)
        legacy_backend_cfg = (
            dict(legacy_backend_raw) if isinstance(legacy_backend_raw, Mapping) else {}
        )
        backend_cfg_raw = aero_cfg.get("backend_config", aero_cfg.get(backend, {}))
        backend_cfg = dict(backend_cfg_raw) if isinstance(backend_cfg_raw, Mapping) else {}
        merged_backend_cfg: dict[str, Any] = {}
        merged_backend_cfg.update(legacy_backend_cfg)
        merged_backend_cfg.update(backend_cfg)
        normalized[backend] = merged_backend_cfg
        return normalized

    if backend == "bem":
        legacy_bem_raw = normalized.get("bem")
        legacy_bem_cfg = dict(legacy_bem_raw) if isinstance(legacy_bem_raw, Mapping) else {}
        backend_cfg_raw = aero_cfg.get("backend_config", aero_cfg.get("bem", {}))
        backend_cfg = dict(backend_cfg_raw) if isinstance(backend_cfg_raw, Mapping) else {}

        merged_bem_cfg: dict[str, Any] = {}
        for key in ("n_blades", "hub_radius"):
            if key in rotor_cfg:
                merged_bem_cfg[key] = rotor_cfg[key]
        merged_bem_cfg.update(legacy_bem_cfg)
        merged_bem_cfg.update(backend_cfg)
        normalized["bem"] = merged_bem_cfg

    return normalized


def build_aero_runtime_context(
    mesh: MeshModel,
    cfg: Mapping[str, Any],
    default_backend: str | None = None,
) -> AeroFSIRuntimeContext:
    """Build normalized aerodynamic runtime context for a participant."""
    normalized_cfg = normalize_aero_fsi_config(cfg, default_backend=default_backend)
    backend = resolve_aero_backend(normalized_cfg, default_backend=default_backend)
    rotor_cfg = _extract_rotor_config(normalized_cfg, backend)

    geometry = None
    surface = None
    if any(name.startswith(ROTOR_BLADE_SET_PREFIX) for name in mesh.element_sets_names):
        geometry = build_full_rotor_aero_geometry(
            mesh,
            rotation_axis=rotor_cfg.get("rotation_axis", [0.0, 1.0, 0.0]),
            rotation_center=rotor_cfg.get("rotation_center", [0.0, 0.0, 0.0]),
            hub_radius=rotor_cfg.get("hub_radius"),
        )
        if backend == "sharpy":
            logger.info(
                "[AERO-SHARPY] Skipping generic mesh-derived aerodynamic surface build"
            )
        else:
            surface = build_full_rotor_aero_surface(geometry)

    return AeroFSIRuntimeContext(
        backend=backend,
        participant_name=str(normalized_cfg.get("participant", "Fluid")),
        coupling_mesh_name=str(normalized_cfg.get("coupling_mesh", "Fluid-Mesh")),
        blade_file=normalized_cfg.get("blade_file"),
        rotor=rotor_cfg,
        normalized_config=normalized_cfg,
        geometry=geometry,
        surface=surface,
    )


def build_aero_participant_from_config(
    mesh: MeshModel,
    cfg: Mapping[str, Any],
    config_file: str = "precice-config.xml",
    viz_mesh: MeshModel | None = None,
    default_backend: str | None = None,
    element_properties: dict | None = None,
):
    """Build the configured aerodynamic participant.

    The new generic config is normalized to the legacy backend-specific input
    expected by the current implementation.

    ``element_properties`` is the deck's element-set property map for the
    coupling mesh. It is consumed by the ``bem`` backend only, where it enables
    the wall-flow moment realisation of the legacy ``BEMFSIParticipant``; the
    ``vlm`` and ``sharpy`` backends ignore it (issue #16).
    """
    runtime_context = build_aero_runtime_context(
        mesh,
        cfg,
        default_backend=default_backend,
    )
    normalized_cfg = runtime_context.normalized_config
    backend = runtime_context.backend

    logger.info("[AERO-FSI] Backend dispatch selected: %s", backend)
    if backend in {"vlm", "sharpy"} and runtime_context.geometry is None:
        logger.warning(
            "[AERO-FSI] Backend '%s' expects full-rotor geometry, but none was built.",
            backend,
        )
    elif runtime_context.geometry is not None:
        logger.info(
            "[AERO-FSI] Runtime geometry ready: blades=%d nodes=%d",
            runtime_context.geometry.n_blades,
            runtime_context.geometry.rotor_mesh.node_count,
        )

    if backend == "bem":
        from aeroelast.solvers.bem.fsi_participant import build_from_config as build_bem_from_config

        logger.info("[AERO-FSI] Building legacy BEM participant...")
        participant = build_bem_from_config(
            mesh,
            normalized_cfg,
            config_file=config_file,
            viz_mesh=viz_mesh,
            element_properties=element_properties,
        )
        # Attached dynamically: the legacy BEM participant declares the attribute
        # (``aero_runtime_context``) for the checker, but only the dispatcher sets it.
        participant.aero_runtime_context = runtime_context  # type: ignore[attr-defined]
        return participant

    if backend == "vlm":
        from .vlm import build_vlm_backend

        logger.info("[AERO-FSI] Building VLM backend...")
        vlm_backend = build_vlm_backend(runtime_context)
        return build_rotor_aero_participant(
            runtime_context,
            vlm_backend,
            config_file=config_file,
        )

    if backend == "sharpy":
        from .sharpy_backend import build_sharpy_backend

        logger.info("[AERO-SHARPY] Building SHARPy backend...")
        sharpy_backend = build_sharpy_backend(runtime_context)
        return build_rotor_aero_participant(
            runtime_context,
            sharpy_backend,
            config_file=config_file,
        )

    raise NotImplementedError(
        f"Aerodynamic backend '{backend}' is not implemented yet. "
        f"Supported backends: {', '.join(SUPPORTED_AERO_FSI_BACKENDS)}"
    )
