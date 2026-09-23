"""Shared types for future aerodynamic backends.

The first implementation slice needs stable contracts before introducing a
specific VLM backend. These dataclasses define the state and load envelopes
that a full-rotor aerodynamic participant will exchange internally.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any, Protocol

import numpy as np

if TYPE_CHECKING:
    from .rotor_geometry import FullRotorAeroGeometry
    from .surface import FullRotorAeroSurface


def _coerce_optional_nodal_vectors(
    field: np.ndarray | None,
    *,
    name: str,
) -> np.ndarray | None:
    if field is None:
        return None

    arr = np.asarray(field, dtype=float)
    if arr.ndim != 2 or arr.shape[1] != 3:
        raise ValueError(f"{name} must have shape (n_nodes, 3)")
    return arr


@dataclass
class RotorAeroState:
    """Operating state passed to an aerodynamic backend.

    Parameters
    ----------
    time : float
        Physical time in seconds.
    dt : float
        Coupling time step in seconds.
    omega_rad_s : float
        Rotor angular velocity in rad/s.
    azimuth_deg : float
        Reference rotor azimuth angle in degrees.
    wind_velocity : ndarray, shape (3,)
        Inflow velocity vector in the global frame.
    nodal_displacements : ndarray, shape (n_nodes, 3), optional
        Structural nodal displacements in the full-rotor ordering.
    nodal_velocities : ndarray, shape (n_nodes, 3), optional
        Structural nodal velocities in the full-rotor ordering.
    extra : dict
        Backend-specific payload for future extensions.
    """

    time: float = 0.0
    dt: float = 0.0
    omega_rad_s: float = 0.0
    azimuth_deg: float = 0.0
    wind_velocity: np.ndarray = field(default_factory=lambda: np.zeros(3, dtype=float))
    nodal_displacements: np.ndarray | None = None
    nodal_velocities: np.ndarray | None = None
    extra: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.wind_velocity = np.asarray(self.wind_velocity, dtype=float).reshape(3)
        self.nodal_displacements = _coerce_optional_nodal_vectors(
            self.nodal_displacements,
            name="nodal_displacements",
        )
        self.nodal_velocities = _coerce_optional_nodal_vectors(
            self.nodal_velocities,
            name="nodal_velocities",
        )


@dataclass
class RotorAeroLoads:
    """Rotor loads returned by an aerodynamic backend.

    Parameters
    ----------
    nodal_forces : ndarray, shape (n_nodes, 3)
        Forces already assembled in the full-rotor structural node ordering.
    integrated_force : ndarray, shape (3,), optional
        Total aerodynamic force over the rotor.
    integrated_moment : ndarray, shape (3,), optional
        Total aerodynamic moment over the rotor.
    metadata : dict
        Auxiliary fields such as bladewise loads, sectional outputs, etc.
    """

    nodal_forces: np.ndarray
    integrated_force: np.ndarray | None = None
    integrated_moment: np.ndarray | None = None
    metadata: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.nodal_forces = np.asarray(self.nodal_forces, dtype=float)
        if self.nodal_forces.ndim != 2 or self.nodal_forces.shape[1] != 3:
            raise ValueError("nodal_forces must have shape (n_nodes, 3)")

        if self.integrated_force is not None:
            self.integrated_force = np.asarray(self.integrated_force, dtype=float).reshape(3)

        if self.integrated_moment is not None:
            self.integrated_moment = np.asarray(self.integrated_moment, dtype=float).reshape(3)


class RotorAerodynamicBackend(Protocol):
    """Protocol implemented by future full-rotor aerodynamic backends."""

    name: str

    def compute_loads(
        self,
        context: "AeroFSIRuntimeContext",
        state: RotorAeroState,
    ) -> RotorAeroLoads:
        """Return aerodynamic loads for the current full-rotor state."""


@dataclass
class AeroFSIRuntimeContext:
    """Normalized runtime context shared by aerodynamic participants.

    This captures the backend selection, the normalized config, and the
    pre-built full-rotor geometry/surface views that future backends can
    consume directly.
    """

    backend: str
    participant_name: str
    coupling_mesh_name: str
    blade_file: str | None
    rotor: dict[str, Any]
    normalized_config: dict[str, Any]
    geometry: "FullRotorAeroGeometry | None" = None
    surface: "FullRotorAeroSurface | None" = None
