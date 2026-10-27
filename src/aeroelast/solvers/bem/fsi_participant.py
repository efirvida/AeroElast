r"""BEM-based preCICE fluid participant with bidirectional FSI coupling.

Runs as an independent process alongside the structural solver, acting as
the "fluid" side of the preCICE coupling (equivalent to OpenFOAM in a
traditional FSI setup but using Blade Element Momentum theory instead of
solving the Navier–Stokes equations).

Bidirectional coupling overview
-------------------------------
Within each preCICE time window the implicit coupling loop (IQN-ILS)
iterates through the following sequence until force–displacement
equilibrium is achieved:

1. **Read displacements** :math:`\mathbf{u}_i` from the structural solver via
   preCICE (mapped from Solid-Mesh to Fluid-Mesh via RBF interpolation).

2. **Update the deformed blade geometry** from the displacement field.
   Two aerodynamically relevant quantities are extracted per BEM strip *k*:

   * The deformed radial position :math:`r_k^{\mathrm{def}}`, which captures
     flapwise and edgewise bending.
   * The elastic twist increment :math:`\Delta\theta_k`, which captures
     torsional deformation around the span axis.

   Chord length is assumed invariant under linear elastic deformation
   (no in-plane stretching).

3. **Rebuild the BEM solver** (CCBlade) with the deformed radii and twist
   distribution.  Airfoil polar assignment remains unchanged — each station
   retains its reference airfoil.

4. **Rebuild the force projector** on the deformed nodal positions so that
   distributed BEM loads are projected onto the *actual* deformed surface.

5. **Evaluate BEM** at the current operating conditions and project forces
   onto mesh nodes.

6. **Write forces** :math:`\mathbf{F}_i` to preCICE for the structural solver.

Steps 1–6 repeat within each time window until preCICE's convergence
criterion is satisfied (see ``precice-config.xml``,
``<relative-convergence-limit>``).

Mathematical formulation
========================

Deformed radial position
------------------------
For each BEM strip *k* with node index set :math:`\mathcal{I}_k`, the
deformed radial coordinate is the deck station radius plus the arithmetic
mean of the strip nodes' span-wise displacement, plus the second-order
shortening of the projected span under transverse bending:

.. math::

    r_k^{\mathrm{def}}
        = r_k^{\mathrm{ref}} + \bar{u}_{s,k}
          - \frac{1}{2}\int_{s_{\mathrm{root}}}^{s_k}
            \left|\frac{d\bar{\mathbf{u}}_\perp}{ds}\right|^2 ds,
    \qquad
    \bar{\mathbf{u}}_k = \frac{1}{|\mathcal{I}_k|}\sum_{i \in \mathcal{I}_k}\mathbf{u}_i

where :math:`\bar{u}_{s,k} = \bar{\mathbf{u}}_k\cdot\hat{\mathbf{e}}_s`,
:math:`\bar{\mathbf{u}}_\perp = \bar{\mathbf{u}} - \bar{u}_s\hat{\mathbf{e}}_s`,
:math:`\hat{\mathbf{e}}_s` is the unit span direction and :math:`s` the mesh span
coordinate of the strip centroids, integrated from the clamped root
(:math:`\bar{\mathbf{u}} = 0` there) with the trapezoid rule.

**Datum.**  The station radius comes from one source only, the deck
(``blade_aero.r``, hub-referenced, the same annulus the CCBlade solver and
the rigid short-circuit use); the mesh contributes only the displacement
increment.  At :math:`\mathbf{u} = 0` the deformed path therefore returns
``_ref_r`` exactly and the rigid short-circuit is a pure optimisation.  The
earlier definition, the mean span projection of the reference node cloud
carried onto the hub datum by one offset, differed from ``_ref_r`` by up to
0.6 m before any deformation (the node-cloud centroid of the end half-cells
and of strips with non-uniform mesh spacing is not the station radius;
``radii_datum_definition_bias``, issue #33 L5).  **Every radial quantity of
this class is hub-referenced**: ``_ref_r``, the ``r_def`` returned here, and
the radii handed to ``_rebuild_bem_solver``.  ``ForceProjector`` anchors its
own station grid on the mesh datum and re-anchors them internally.

**Why the second-order term.**  Under a geometrically linear structural
solve :math:`\bar{u}_s` is first order: at the rated point it is almost
entirely the prebent axis straightening under flap bending
(:math:`-\int y_0' u_y'\,ds`).  The shortening of a near-inextensible axis
under the same bending is second order and the linear solve cannot produce
it, yet at a 17 m tip deflection it is about twice the linear term with the
opposite sign, so without it the tip is fed 1.1 m *outward* while it really
moves inward.  The integral above is the moderate-rotation completion of
the inextensibility condition :math:`(1+w')^2 + (y_0'+v')^2 = 1 + y_0'^2`;
the linear part stays with the structural solve.  Reference
(``radii_linear_kinematics_axial_term``, ``docs/validation_closures.md``):
on BeamDyn's own geometrically exact static field under the same load the
formula reproduces BeamDyn's tip axial displacement (-1.339 m) to 0.3 % of
the second-order term.  The truncation is :math:`O(\theta^4)` in the
deflection slope.

Elastic twist extraction
------------------------
The local aerodynamic twist at strip *k* on the deformed blade is:

.. math::

    \theta_k^{\mathrm{def}} = \theta_k^{\mathrm{ref}} + \Delta\theta_k

**The section rotation, not the strip cloud's chord direction.**  A BEM strip
spans several merged mesh stations (``dr ~ 2.4 m`` against ~0.63 m here), so
the principal in-plane direction of the strip's *whole* node cloud follows the
deformed arc of a bent blade instead of its chord; the elastic twist is
measured on **one physical ring** at the strip centre.  The participant
receives only the three translational DOFs through the coupling, so the
rotation is recovered from the ring's in-plane displacement field
:math:`(u_x, u_y)` at :math:`(x, y)` in an orthonormal basis of the plane
perpendicular to :math:`\hat{\mathbf{e}}_s`, fitted by least squares to the
affine map

.. math::

    u_x = a_{11} x + a_{12} y + t_x, \qquad
    u_y = a_{21} x + a_{22} y + t_y.

Its antisymmetric part is the section rotation about the span axis:

.. math::

    \Delta\theta_k = \tfrac{1}{2}\,(a_{21} - a_{12}),

invariant to the in-plane origin and with the same sense as the shell's own
rotational DOF about the span axis (positive nose-up).  Each strip's ring is
bound to the projector's physical ring (``physical_ring_labels``: span-gap
single linkage joined along in-ring element edges) that holds the node nearest
the strip centre.  A fraction of the strip width is not a ring tolerance: it
merged the tip rings of the 0.25 m mesh, 0.075-0.109 m apart, into one fit, and
a fixed span gap sheared prebent rings (skew 0.0121 m) into arcs.  A ring that cannot define a 2-D
section keeps the reference twist: no second estimator is used, so nothing can
silently disagree with this one.

**Sense.**  The BEM twist angle and the mesh's rotation about the span axis are
opposite-sensed whenever ``normal = span x tangential`` (the default frame:
CCBlade's ``alpha = phi - theta``, while a positive mesh rotation turns the
chord toward the downwind axis).  The factor is derived from the configured
directions in ``_twist_mesh_to_bem``.

Force projection on deformed geometry
--------------------------------------
After obtaining the BEM distributed loads :math:`N_p(r_k)` and
:math:`T_p(r_k)` (normal and tangential force per unit span), the loads are
laid out by **one** ``ForceProjector`` built once on the **reference**
geometry.  The node-to-strip assignment, the aerodynamic-centre arm, the
pitching-moment axis and the moment-preserving distribution all stay on the
reference: re-deriving them from the current displacement every implicit
sub-iteration adds a deformation-dependent gain that breaks the coupled
contraction (``odd/tasks/coupled-divergence.md``).  The loads themselves do
follow the deformation, through the deformed BEM.

Aerodynamic pitching moment
----------------------------
The BEM solver computes a pitching-moment coefficient :math:`C_m` at each
station (from the airfoil polars evaluated at the local angle of attack).
The pitching moment per unit span is recovered from the distributed loads:

.. math::

    M_p(r_k) = \frac{\sqrt{N_p^2 + T_p^2}}{\sqrt{C_l^2 + C_d^2}}
               \;\cdot\; c_k \;\cdot\; C_m(r_k)

The first factor is the dynamic pressure per unit span
:math:`q_c = \tfrac{1}{2}\rho W^2 c`, back-calculated from the force
resultant and the aerodynamic coefficients.  This avoids requiring the
relative wind speed :math:`W` as an explicit BEM output.

The total moment on BEM strip *k* is :math:`M_p \, \Delta r_k` and acts
about the span axis :math:`\hat{\mathbf{e}}_s`.  It is passed to the
``ForceProjector`` as the target moment vector per strip, so that the
minimum-norm force distribution creates nodal forces whose net moment
about the strip centroid matches the aerodynamic pitching moment.

Consistent deformed strip widths
---------------------------------
The projector is not rebuilt on the deformed ``BladeAero``; only its
per-strip width :math:`\Delta r_k`, used to convert force-per-span into total
strip force, follows the stations the BEM result carries.  ``project()``
recomputes it from :math:`r_k^{\mathrm{def}}` with the same strip rule that
built the reference grid, so the deformed loads are integrated over the
deformed (shortened) span rather than the reference widths.

Checkpoint strategy
-------------------
BEM is memoryless: given a displacement field and operating conditions it
produces a deterministic force field.  The checkpoint therefore stores only
the last committed nodal forces :math:`\mathbf{F}^{n}`, managed by the
shared :class:`~aeroelast.solvers.fsi.base.Adapter` (which provides
checkpoint save/restore via preCICE's implicit-iteration protocol).  On
rollback the previous forces are restored and a fresh BEM evaluation is
performed using the new displacement iterate.

Physical assumptions and limitations
-------------------------------------
* **Linear elasticity**: chord length is invariant.  For geometrically
  nonlinear structural models the chord could be updated analogously to *r*
  and *twist*, but this is not implemented.
* **Rigid airfoil profiles**: deformation changes the angle of attack via
  twist but does not alter the airfoil shape itself.
* **Quasi-steady aerodynamics**: BEM evaluates steady-state loads at each
  implicit iteration; dynamic stall and unsteady wake effects are not
  captured.
* **Single-blade formulation**: the participant models one blade; multi-blade
  effects are handled through the BEM rotor model (rotor thrust factor *B*).
"""

from __future__ import annotations

import copy
import logging
import math
from pathlib import Path
from typing import Optional

import numpy as np

from aeroelast.core.mesh.model import MeshModel
from aeroelast.models.blade.aerodynamics import (
    AeroStation,
    BladeAero,
    load_blade_aero,
)
from aeroelast.solvers.bem.engine import BEMResult, BEMSolver
from aeroelast.solvers.bem.force_projection import ForceProjector

logger = logging.getLogger(__name__)

#: Frame the structural partner exchanges at the preCICE boundary: the structural
#: solver is always the corotational Rust rotor solver, which works in inertial axes.
PARTNER_FRAMES = ("inertial",)


# ---------------------------------------------------------------------------
# Rotor-frame transforms at the preCICE boundary
# ---------------------------------------------------------------------------
#
# The Rust rotor solver (``crates/aeroelast-solvers/src/petsc/fsi/rotor_fsi.rs``)
# writes its interface displacements as ``disps_to_inertial(u, theta_target)``
# = R(theta) . u and reads the fluid forces as global, applying
# ``forces_to_rotating`` = R^T . F.  R is the right-handed Rodrigues matrix of
# ``RotorTransforms::rotation_matrix`` (``rotor_physics.rs``) and theta_target
# is the END-of-window rotor angle, starting at 0.  The helpers below restate
# that matrix so the BEM participant (which works in the rotating blade frame)
# can undo it on read and apply it on write.


def _rodrigues(axis, theta: float) -> np.ndarray:
    """Right-handed rotation ``R = I + sin(t) K + (1 - cos(t)) K^2`` about ``axis``.

    Same matrix as Rust's ``RotorTransforms::rotation_matrix``; ``axis`` is
    normalised here, as Rust does.
    """
    n = np.asarray(axis, dtype=float).reshape(3)
    norm = float(np.linalg.norm(n))
    if norm <= 1e-12:
        raise ValueError("rotation axis must be non-zero")
    n = n / norm
    k = np.array([[0.0, -n[2], n[1]], [n[2], 0.0, -n[0]], [-n[1], n[0], 0.0]])
    return np.eye(3) + np.sin(theta) * k + (1.0 - np.cos(theta)) * (k @ k)


def _to_rotating(vectors: np.ndarray, axis, theta: float) -> np.ndarray:
    """Row-wise ``R(theta)^T . v``: inertial vectors to the rotating frame."""
    return np.asarray(vectors, dtype=float) @ _rodrigues(axis, theta)


def _to_inertial(vectors: np.ndarray, axis, theta: float) -> np.ndarray:
    """Row-wise ``R(theta) . v``: rotating vectors to the inertial frame."""
    return np.asarray(vectors, dtype=float) @ _rodrigues(axis, theta).T


def _parse_partner_frame(bem_config: dict) -> str:
    """Read ``partner_frame``; only ``"inertial"`` (the default) is accepted.

    The key is kept so existing YAMLs that declare ``inertial`` still load.  A
    ``rotating`` partner no longer exists: the structural solver is always the
    corotational Rust rotor solver.
    """
    frame = str(bem_config.get("partner_frame", "inertial")).strip().lower()
    if frame == "rotating":
        raise ValueError(
            "bem.partner_frame: 'rotating' is no longer supported; the structural partner is "
            "the corotational rotor solver and exchanges displacements and forces in inertial "
            "axes. Remove the key or set partner_frame: inertial"
        )
    if frame not in PARTNER_FRAMES:
        raise ValueError(f"bem.partner_frame must be one of {PARTNER_FRAMES}, got {frame!r}")
    return frame


def _check_theta_source(omega_mesh: str | None, omega_initial: float) -> None:
    """Fail when the structure's rotor angle cannot reach the participant.

    The corotational structure publishes theta on the omega mesh (decision D4); without
    that mesh the participant has no angle.  That is fine for a parked rotor
    (``omega == 0``, theta stays 0) and impossible for a spinning one.
    """
    if omega_mesh is None and omega_initial != 0.0:
        raise ValueError(
            "bem.omega: a spinning rotor needs the omega mesh to receive the rotor angle theta "
            f"from the structure (static omega={omega_initial} rad/s, no bem.omega.mesh); "
            "declare bem.omega.mesh / data / vertex and the Theta data in the preCICE config"
        )


def _parse_rotation_axis(bem_config: dict) -> np.ndarray:
    """Read ``rotation_axis`` (default +Y, the solid's ``rotor.rotation_axis``)."""
    axis = np.asarray(bem_config.get("rotation_axis", [0.0, 1.0, 0.0]), dtype=float).reshape(3)
    norm = float(np.linalg.norm(axis))
    if norm <= 1e-12:
        raise ValueError("bem.rotation_axis must be a non-zero 3-vector")
    return axis / norm


# ---------------------------------------------------------------------------
# Main participant class
# ---------------------------------------------------------------------------


class BEMFSIParticipant:
    r"""preCICE fluid participant that computes aerodynamic forces via BEM
    with full bidirectional deformation feedback.

    On every implicit coupling iteration the participant:

    1. Reads the structural displacement field :math:`\mathbf{u}`.
    2. Extracts deformed strip radii :math:`r_k^{\mathrm{def}}` and elastic
       twist :math:`\theta_k^{\mathrm{def}}` (see module docstring for the
       mathematical formulation).
    3. Rebuilds ``BEMSolver`` and ``ForceProjector`` on the deformed geometry.
    4. Evaluates BEM and projects nodal forces.
    5. Writes forces to preCICE for the structural solver.

    Parameters
    ----------
    mesh : MeshModel
        Blade surface mesh.  Its nodes define the Fluid-Mesh vertices
        registered with preCICE.
    blade_aero : BladeAero
        Aerodynamic blade definition (airfoil polars, radial stations).
    bem_config : dict
        BEM operating conditions and solver options.  Recognised keys:

        * ``wind_speed``           – free-stream wind speed [m/s]
        * ``omega``                – rotor angular velocity [RPM]
        * ``pitch``                – collective pitch angle [deg]
        * ``azimuth``              – blade azimuth angle [deg] (default 0)
        * ``air_density``          – ρ [kg/m³] (default 1.225)
        * ``dynamic_viscosity``    – μ [Pa·s]  (default 1.81e-5)
        * ``hub_height``           – hub centre height above ground [m]
        * ``shear_exp``            – wind shear exponent (default 0.2)
        * ``precone``              – precone angle [deg] (default 0)
        * ``tilt``                 – shaft tilt angle [deg] (default 0)
        * ``yaw``                  – nacelle yaw misalignment angle [deg] (default 0)
        * ``span_direction``       – unit vector along blade span (default Z)
        * ``normal_direction``     – sense reference for the section normal
          the BEM *Np* rides (default Y, the fluid/downwind direction); it
          must be (nearly) parallel to that normal, it is not the axis
        * ``tangential_direction`` – optional; *Tp* follows the measured leading
          edge and a value on the opposite side raises ``ValueError`` (default:
          unset, so the sense is derived from the mesh)
        * ``partner_frame``        – optional, only ``"inertial"`` (the default) is
          accepted; ``"rotating"`` raises ``ValueError``.  The structural partner is
          the corotational Rust rotor solver: displacements are un-rotated by
          R(theta)^T on read and forces rotated by R(theta) on write, theta the
          end-of-window rotor angle
        * ``rotation_axis``        – rotor axis for that transform (default
          ``[0, 1, 0]``); must match the solid's ``rotor.rotation_axis``
    participant : str
        preCICE participant name (must match the XML config).
    config_file : str or Path
        Path to ``precice-config.xml``.
    coupling_mesh : str
        Name of the mesh provided by this participant.
    displacement_data : str
        Name of the displacement data field to read.
    force_data : str
        Name of the force data field to write.
    output_folder : str or Path
        Directory for CSV diagnostics.
    log_interval : int
        How many time windows between diagnostic log lines.
    """

    #: Set by the aerodynamic dispatcher (``build_aero_participant_from_config``)
    #: when this participant is built through the generic ``aeroelast`` entry point;
    #: the dedicated ``aeroelast-bem-fsi`` CLI leaves it unset. Declared for the
    #: type checker - it is deliberately not created here.
    aero_runtime_context: object | None

    def __init__(
        self,
        mesh: MeshModel,
        blade_aero: BladeAero,
        bem_config: dict,
        participant: str = "Fluid",
        config_file: str | Path = "precice-config.xml",
        coupling_mesh: str = "Fluid-Mesh",
        displacement_data: str = "Displacement",
        force_data: str = "Force",
        output_folder: str | Path = "bem_fsi_results",
        log_interval: int = 10,
        viz_mesh: MeshModel | None = None,
        omega_mesh: str | None = None,
        omega_data: str = "AngularVelocity",
        omega_vertex: list[float] | None = None,
        element_properties: dict | None = None,
        theta_data: str = "Theta",
    ) -> None:
        self._mesh = mesh
        self._blade_aero = blade_aero
        self._bem_cfg = bem_config
        self._participant_name = participant
        self._config_file = str(config_file)
        self._coupling_mesh = coupling_mesh
        self._displacement_data = displacement_data
        self._force_data = force_data
        self._output_folder = Path(output_folder)
        self._log_interval = log_interval

        # omega preCICE read config (optional)
        self._omega_mesh: str | None = omega_mesh
        self._omega_data: str = omega_data
        # Rotor angle theta [rad] published by the structure on the omega mesh.
        self._theta_data: str = theta_data
        self._omega_vertex: np.ndarray = (
            np.asarray(omega_vertex, dtype=float) if omega_vertex is not None else np.zeros(3)
        )
        # Live omega [rad/s] — updated each time window from preCICE read or YAML fallback.
        # Converted to RPM inside _compute_forces before passing to CCBlade.
        self._current_omega: float = float(bem_config.get("omega", 0.0))

        # Reference node coordinates (n_nodes, 3) — immutable after init
        self._ref_coords: np.ndarray = mesh.coords_array.copy()
        self._n_nodes: int = self._ref_coords.shape[0]

        # -- Direction vectors (normalised once) ----------------------------
        self._span_dir = np.asarray(
            bem_config.get("span_direction", [0.0, 0.0, 1.0]),
            dtype=float,
        )
        self._span_dir /= np.linalg.norm(self._span_dir)

        self._normal_dir = np.asarray(
            bem_config.get("normal_direction", [0.0, 1.0, 0.0]),
            dtype=float,
        )
        self._normal_dir /= np.linalg.norm(self._normal_dir)

        self._configured_tangential = bem_config.get("tangential_direction")
        self._tangential_dir = np.asarray(
            self._configured_tangential
            if self._configured_tangential is not None
            else [1.0, 0.0, 0.0],
            dtype=float,
        )
        self._tangential_dir /= np.linalg.norm(self._tangential_dir)

        # -- BEM solver keyword arguments (constant across iterations) ------
        self._bem_solver_kwargs: dict = {
            "rho": bem_config.get("air_density", 1.225),
            "mu": bem_config.get("dynamic_viscosity", 1.81206e-5),
            "precone": bem_config.get("precone", 0.0),
            "tilt": bem_config.get("tilt", 0.0),
            "yaw": bem_config.get("yaw", 0.0),
            "hub_height": bem_config.get("hub_height", 150.0),
            "shear_exp": bem_config.get("shear_exp", 0.2),
        }

        # -- Reference BEM solver and projector -----------------------------
        self._bem_solver = BEMSolver(blade_aero, **self._bem_solver_kwargs)

        # ``element_properties`` carries the deck's section-property map to the
        # wall-flow moment realisation.  The CLI builds it from the same generator
        # that produced the element-bearing coupling mesh, so every wall's shear
        # stiffness ``S = G*t`` resolves and ``_ring_section_is_realisable`` accepts
        # the precomputed ring.  ``None`` (a direct caller, or a property-less mesh)
        # keeps the minimum-norm ``_distribute`` fallback.
        ref_projector = ForceProjector(
            mesh,
            blade_aero,
            span_direction=self._span_dir,
            normal_direction=self._normal_dir,
            tangential_direction=self._configured_tangential,
            element_properties=element_properties,
        )
        self._projector = ref_projector

        # -- Mesh span rotation -> BEM twist sense --------------------------
        # A positive rotation of the section about +span carries its LEADING EDGE
        # toward ``+span x leading_edge``.  The deck puts the leading edge at +x and
        # the load frame puts the downwind (normal) axis at +y, so that motion is
        # downwind: the section pitches NOSE-UP (toward stall), the angle of attack
        # RISES, and CCBlade's ``alpha = phi - theta`` therefore requires ``theta`` to
        # FALL.  The mesh rotation and the BEM twist angle are opposite-sensed and the
        # factor is -1.  ``_resolve_twist_sense`` reads the leading edge off the
        # projector, which measures it per strip against the deck's own aerofoil tables.
        self._twist_mesh_to_bem: float = self._resolve_twist_sense()

        # -- Strip-to-node mapping from the reference projector -------------
        # These index arrays define which mesh nodes belong to each BEM strip.
        # The mapping is established on the reference geometry and is used by
        # the deformed-geometry methods (radial position averaging and the
        # strip's centre-ring selection).  The projector keeps this same
        # reference assignment for the load distribution.
        self._strip_node_indices: list[np.ndarray] = [
            strip.node_indices.copy() for strip in ref_projector._strips
        ]

        # -- Reference per-strip quantities (all radial ones HUB-referenced) --
        self._ref_r: np.ndarray = blade_aero.r.copy()  # (n_strips,)
        self._ref_twist: np.ndarray = blade_aero.twist.copy()  # (n_strips,) rad
        # Mesh span coordinate (mesh datum).  Radii fed to the BEM never use it:
        # they are _ref_r plus a displacement increment (_compute_deformed_geometry).
        self._mesh_span: np.ndarray = self._ref_coords @ self._span_dir
        # Abscissa of the second-order shortening integral: the clamped root and
        # each strip's centroid on the mesh span coordinate (NaN for an empty strip).
        self._root_span: float = float(np.min(self._mesh_span))
        self._strip_span: np.ndarray = np.array(
            [self._mesh_span[idx].mean() if len(idx) else np.nan for idx in self._strip_node_indices]
        )
        # Mesh-datum strip centres and widths, from the reference projector.
        self._strip_centres: np.ndarray = np.array([s.r_center for s in ref_projector._strips])
        self._strip_widths: np.ndarray = np.array([s.dr for s in ref_projector._strips])
        # The physical ring at each strip's centre: the section the elastic
        # twist is measured on (a strip spans ~4 mesh stations, so its whole
        # node cloud is not one section).
        self._strip_ring_indices: list[np.ndarray] = self._find_strip_rings(ref_projector._strip_ring_groups)

        # -- Working copy of the mesh for projector reconstruction ----------
        # A single deep copy is made at init time; subsequent iterations
        # mutate its coords_array via the setter (O(n_nodes) per call)
        # instead of creating a new deep copy every iteration.
        self._working_mesh: MeshModel = copy.deepcopy(mesh)

        # -- Full mesh with elements for surface VTU output -----------------
        # When a coupling_node_set filter is applied in the CLI, the coupling
        # mesh has no elements. viz_mesh keeps the full mesh for ParaView output.
        self._viz_mesh: MeshModel = viz_mesh if viz_mesh is not None else mesh

        # -- Tip node index (max span coordinate) --------------------------
        # disp[-1] is NOT guaranteed to be the tip — node ordering in the
        # coupling mesh depends on the mesh generator + coupling_node_set filter.
        span_coords = self._ref_coords @ self._span_dir  # (n_nodes,)
        self._tip_node_idx: int = int(np.argmax(span_coords))

        # -- Counters and live azimuth --------------------------------------
        self._last_forces = np.zeros((self._n_nodes, 3), dtype=float)
        self._previous_iteration_forces: Optional[np.ndarray] = None
        self._window_count = 0
        self._iteration_count = 0
        # CCBlade azimuth [deg] at structure angle theta = 0 (static YAML value).
        self._azimuth_offset: float = float(bem_config.get("azimuth", 0.0))

        # -- Partner frame at the preCICE boundary (T7a) ---------------------
        # The partner (Rust rotor_fsi.rs, always corotational) exchanges global
        # vectors, so displacements are un-rotated on read and forces rotated on
        # write; the transform is unconditional.
        self._partner_frame: str = _parse_partner_frame(bem_config)
        self._rotation_axis: np.ndarray = _parse_rotation_axis(bem_config)
        # Rotor angle theta [rad] read from the structure (Rust ``theta_target``, the
        # end-of-window angle it rotates forces and displacements with).  This is the only
        # source of the angle: the participant never integrates omega itself.  It is also
        # the CCBlade azimuth (see ``_azimuth_deg``).
        self._partner_theta: float = 0.0

    # -----------------------------------------------------------------------
    # Rotor-frame boundary
    # -----------------------------------------------------------------------

    def _read_rotor_state(self, adapter) -> None:
        """Read omega and theta from the omega mesh (same time level, end of window)."""
        if self._omega_mesh is None:
            return
        omega_arr = adapter.read_data(self._omega_mesh, self._omega_data)
        if omega_arr is not None and np.asarray(omega_arr).size > 0:
            previous_omega = self._current_omega
            new_omega = float(np.asarray(omega_arr).ravel()[0])
            changed = abs(new_omega - previous_omega) > 1e-6
            logger.info(
                "[BEM-FSI] ω received from preCICE: %.6f%s rad/s (window %d, iteration %d)",
                new_omega,
                f" (was {previous_omega:.6f})" if changed else " unchanged",
                self._window_count,
                self._iteration_count,
            )
            self._current_omega = new_omega
        theta_arr = adapter.read_data(self._omega_mesh, self._theta_data)
        if theta_arr is not None and np.asarray(theta_arr).size > 0:
            self._partner_theta = float(np.asarray(theta_arr).ravel()[0])

    def _require_partner_theta(self, adapter) -> None:
        """Fail at setup when the preCICE configuration does not carry the theta data."""
        try:
            adapter.precice.get_data_dimensions(self._omega_mesh, self._theta_data)
        except Exception as exc:  # preCICE raises RuntimeError for unknown data
            raise RuntimeError(
                f"the structural partner publishes the rotor angle as data '{self._theta_data}' "
                f"on mesh '{self._omega_mesh}' (the BEM takes its frame rotation and CCBlade "
                "azimuth from it), but the preCICE configuration does not declare it: add "
                f'<data:scalar name="{self._theta_data}" />, use it in both omega meshes, '
                "write it from the structure, read it here, and exchange it from the structure "
                f"(initialize=\"true\"). preCICE said: {exc}"
            ) from exc

    def _azimuth_deg(self) -> float:
        """CCBlade azimuth [deg] for the current angle: ``bem.azimuth + degrees(theta)``.

        Both CCBlade's azimuth and the structure's theta are positive rotations of the rotor
        about its shaft (see ``tests/validation/rotor/test_bem_theta_from_structure.py``).
        Same time level as the frame rotation, the end of the window: the loads computed
        here are applied to the blade at the angle the structure rotates them with.
        """
        return (self._azimuth_offset + math.degrees(self._partner_theta)) % 360.0

    def _displacements_from_partner(self, displacements: np.ndarray) -> np.ndarray:
        """Bring the read displacements into the rotating (blade) frame."""
        if len(displacements) == 0:
            return displacements
        return _to_rotating(displacements, self._rotation_axis, self._partner_theta)

    def _forces_to_partner(self, forces: np.ndarray) -> np.ndarray:
        """Express the rotating-frame forces in the partner's frame for writing."""
        if len(forces) == 0:
            return forces
        return _to_inertial(forces, self._rotation_axis, self._partner_theta)

    # -----------------------------------------------------------------------
    # Public API
    # -----------------------------------------------------------------------

    def run(self) -> None:
        """Execute the preCICE coupling loop (blocking until finalized)."""
        # Lazy import — the Adapter (and its precice dependency) are only
        # loaded when actually running, avoiding GLIBCXX errors on login nodes.
        from aeroelast.solvers.fsi.base import Adapter

        self._output_folder.mkdir(parents=True, exist_ok=True)

        _coupling_meshes = {self._coupling_mesh: self._ref_coords}
        if self._omega_mesh is not None:
            _coupling_meshes[self._omega_mesh] = self._omega_vertex.reshape(1, 3)
            logger.info(
                "[BEM-FSI] omega preCICE read ENABLED: mesh=%s  data=%s  vertex=%s",
                self._omega_mesh,
                self._omega_data,
                self._omega_vertex.tolist(),
            )
        else:
            logger.info(
                "[BEM-FSI] omega preCICE read DISABLED "
                "(no omega_mesh configured — using static yaml omega=%.6f rad/s)",
                self._current_omega,
            )

        adapter = Adapter(
            participant=self._participant_name,
            config_file=self._config_file,
            coupling_meshes=_coupling_meshes,
        )
        if self._omega_mesh is not None:
            self._require_partner_theta(adapter)
        adapter.initialize()

        logger.info(
            "[BEM-FSI] Initialized.  Participant=%s  mesh=%s  nodes=%d  strips=%d",
            self._participant_name,
            self._coupling_mesh,
            self._n_nodes,
            len(self._strip_node_indices),
        )

        bem_result: Optional[BEMResult] = None

        while adapter.is_coupling_ongoing:
            # preCICE 3.x requires explicit checkpoint write/read guards.
            if adapter.requires_writing_checkpoint:
                adapter.store_checkpoint((self._last_forces.copy(),))

            if adapter.requires_reading_checkpoint:
                states = adapter.retrieve_checkpoint()
                self._last_forces = states[0]
                self._iteration_count += 1

            dt = adapter.dt

            # Read displacements from the structural solver (partner frame;
            # brought to the rotating frame below, once omega is current).
            displacements_read = adapter.read_data(self._coupling_mesh, self._displacement_data)

            # Live rotor speed and rotor angle from the structure.
            self._read_rotor_state(adapter)

            # From here on ``displacements`` is ROTATING-frame (blade frame):
            # the BEM geometry, disp_max, the log lines, the timestep output and
            # the bem_report.csv tip columns all use it.  With an inertial
            # partner it is R(theta_end)^T . u_read, theta_end read just above.
            displacements = self._displacements_from_partner(displacements_read)
            if len(displacements) == 0:
                disp_max = 0.0
            else:
                disp_max = float(np.max(np.linalg.norm(displacements, axis=1)))

            # BEM on (possibly deformed) geometry
            forces, bem_result = self._compute_forces(displacements)
            force_norm = float(np.linalg.norm(forces))
            if self._previous_iteration_forces is None:
                force_delta = 0.0
            else:
                force_delta = float(np.linalg.norm(forces - self._previous_iteration_forces))
            logger.info(
                "[BEM-FSI] coupling iter: window=%d iteration=%d | geometry=%s "
                "| max_disp=%.6e m | ||F||=%.6e N | ΔF_iter=%.6e N",
                self._window_count,
                self._iteration_count,
                "reference" if disp_max < 1e-12 else "deformed",
                disp_max,
                force_norm,
                force_delta,
            )
            self._previous_iteration_forces = forces.copy()

            # ``forces`` stays rotating-frame (checkpoint, _last_forces,
            # _previous_iteration_forces, output); only the written copy is
            # expressed in the partner's frame.
            adapter.write_data(
                self._coupling_mesh, self._force_data, self._forces_to_partner(forces)
            )

            adapter.advance(dt)

            if adapter.is_time_window_complete:
                self._last_forces = forces
                self._previous_iteration_forces = None
                self._window_count += 1
                self._iteration_count = 0
                current_time = self._window_count * dt

                self._write_timestep_output(
                    forces,
                    displacements,
                    bem_result,
                    self._window_count,
                    current_time,
                )

                if self._window_count % self._log_interval == 0:
                    self._log_window(
                        bem_result,
                        displacements,
                        current_time,
                    )

        adapter.finalize()
        logger.info("[BEM-FSI] Finalized after %d time windows.", self._window_count)

    # -----------------------------------------------------------------------
    # Deformed geometry extraction
    # -----------------------------------------------------------------------

    def _compute_deformed_geometry(
        self,
        displacements: np.ndarray,
    ) -> tuple[np.ndarray, np.ndarray]:
        r"""Extract deformed BEM strip radii and twist from nodal displacements.

        For the full mathematical derivation see the module-level docstring.

        **Algorithm summary** (per strip *k*):

        1. Deformed coordinates :math:`\mathbf{x}_i = \mathbf{X}_i + \mathbf{u}_i`.
        2. Radial position: the deck station radius ``_ref_r[k]`` plus the mean span-wise
           displacement of the strip nodes (see the module docstring's "Datum"), minus
           the second-order shortening of the span (``_span_shortening``).
        3. Elastic twist: the in-plane section rotation of the strip's centre physical
           ring (``_section_rotation``), converted into the BEM's twist sense.
        4. Deformed twist:
           :math:`\theta_k^{\mathrm{def}} = \theta_k^{\mathrm{ref}} + \Delta\theta_k`.

        Parameters
        ----------
        displacements : ndarray, shape (n_nodes, 3)
            Nodal displacement vectors from the structural solver.

        Returns
        -------
        r_def : ndarray, shape (n_strips,)
            Deformed radial position of each BEM strip on the **hub datum** (the same
            datum as ``_ref_r``).
        twist_def : ndarray, shape (n_strips,)
            Total aerodynamic twist at each strip [rad], reference plus the elastic
            torsion increment, in the BEM's twist sense.
        """
        n_strips = len(self._strip_node_indices)
        r_def = np.empty(n_strips)
        twist_def = np.empty(n_strips)

        s = self._span_dir
        u_mean = np.zeros((n_strips, 3))
        for k in range(n_strips):
            idx = self._strip_node_indices[k]
            if len(idx) == 0:
                r_def[k] = self._ref_r[k]
                twist_def[k] = self._ref_twist[k]
                continue

            # Deformed radial position, first order: the deck station radius plus the
            # strip's mean span-wise displacement, so u = 0 returns _ref_r exactly (see
            # the module docstring's "Datum" paragraph).
            u_mean[k] = displacements[idx].mean(axis=0)
            r_def[k] = self._ref_r[k] + float(u_mean[k] @ s)

            # Elastic twist: the rotation of the strip's centre physical ring
            # about the span axis, taken into the BEM's twist sense (see
            # _twist_mesh_to_bem).  A degenerate ring has no measurable section
            # rotation, so that strip keeps its reference twist rather than
            # falling back to a second, disagreeing estimator.
            omega = self._section_rotation(displacements, self._strip_ring_indices[k])
            twist_def[k] = self._ref_twist[k] + self._twist_mesh_to_bem * (
                0.0 if omega is None else omega
            )

        r_def -= self._span_shortening(u_mean)
        return r_def, twist_def

    def _span_shortening(self, u_mean: np.ndarray) -> np.ndarray:
        r"""Second-order shortening :math:`\frac12\int_{root}^{s_k}|\bar{\mathbf{u}}_\perp'|^2 ds`.

        The linear structural solve carries the first-order span-wise displacement only;
        this is the term it cannot produce (module docstring, "Why the second-order
        term").  The integral runs over the strip centroids sorted along the mesh span
        coordinate, starting at the clamped root where the displacement and its slope
        vanish.  Empty strips take no part and get zero.

        Parameters
        ----------
        u_mean : ndarray, shape (n_strips, 3)
            Mean nodal displacement of each strip.

        Returns
        -------
        ndarray, shape (n_strips,)
            Non-negative shortening [m] to subtract from the first-order radius.
        """
        out = np.zeros(len(u_mean))
        valid = np.flatnonzero(np.isfinite(self._strip_span))
        if len(valid) == 0:
            return out
        order = valid[np.argsort(self._strip_span[valid])]
        s = self._span_dir
        u_perp = u_mean[order] - np.outer(u_mean[order] @ s, s)
        x = np.concatenate([[self._root_span], self._strip_span[order]])
        slope = np.gradient(np.vstack([np.zeros(3), u_perp]), x, axis=0)
        g = np.sum(slope**2, axis=1)
        out[order] = 0.5 * np.cumsum(0.5 * (g[1:] + g[:-1]) * np.diff(x))
        return out

    def _resolve_twist_sense(self) -> float:
        r"""``-1`` when a positive section rotation about :math:`+\hat{e}_s` is nose-up.

        The sense is a property of the **section frame** - where the leading edge is
        relative to the downwind axis - and not of the aerodynamic load senses.  It used
        to be read from ``normal . (span x tangential)``, which conflates the two:
        ``tangential_direction`` configures the sense of the tangential force
        :math:`T_p` (``F = N_p \hat{n} + T_p \hat{c}``), and every production case YAML
        used to set it to ``-x``.  That was the reversed sense (``z x (-x) = -y`` opposes
        ``Omega = +y``; :math:`T_p` must point toward the leading edge, see
        ``ForceProjector._resolve_tangential_sense``) and it made
        ``span x tangential = -normal``, so the old expression returned ``-1`` in
        production while the static fixtures - which never pass ``tangential_direction`` -
        returned ``+1``.  A rigid +2 deg rotation measured ``+10.63%`` of
        thrust on the production frame against ``-13.05%`` on the fixtures' frame, and
        the inverted feedback is what ran the coupled rotor away at ``t ~ 1.9 s``
        (``odd/tasks/coupled-divergence.md``).

        :class:`~aeroelast.solvers.bem.force_projection.ForceProjector` already measures
        the missing piece: ``_section_ends`` identifies the leading edge on each strip's
        station ring and ``_strip_moment_axis_sign`` records whether the strip's chord
        axis runs leading-to-trailing (``+1``) or trailing-to-leading (``-1``), arbitrated
        against the deck's own WindIO aerofoils.  The leading-edge direction is therefore
        ``-_strip_moment_axis_sign * _strip_chord_dirs``, and a positive rotation about
        ``+span`` is nose-up (``-1``) when that direction has a positive component along
        ``normal x span``, nose-down (``+1``) otherwise.  The strips are weighted by their
        width, as the projector weights its own frame resolution.

        Falls back to the configured frame's expression when there is no measurable
        section outline (no chord axis): a mesh without sections has no leading edge, and
        the caller that holds one must say which way ``T_p`` points.
        """
        chords = self._projector._strip_chord_dirs
        signs = self._projector._strip_moment_axis_sign
        weights = np.array([strip.dr for strip in self._projector._strips], dtype=float)
        n = min(len(chords), len(signs), len(weights))
        leading_edge = np.zeros(3)
        for k in range(n):
            leading_edge += weights[k] * -float(signs[k]) * np.asarray(chords[k], dtype=float)
        if n == 0 or float(np.linalg.norm(leading_edge)) <= 1e-12:
            return (
                -1.0
                if float(np.dot(self._normal_dir, np.cross(self._span_dir, self._tangential_dir)))
                > 0.0
                else 1.0
            )
        sense = float(leading_edge @ np.cross(self._normal_dir, self._span_dir))
        return -1.0 if sense > 0.0 else 1.0

    def _find_strip_rings(self, ring_groups: list[list[np.ndarray]]) -> list[np.ndarray]:
        """Assign every BEM strip the physical mesh ring nearest its centre.

        ``ring_groups`` are the projector's own per-strip physical rings (indices into
        the strip's nodes, :func:`physical_ring_labels`), so the twist is fitted on the
        same ring the wall flow is realised on and no second tolerance exists.  The ring
        is the group holding the strip node nearest the strip centre.
        """
        rings: list[np.ndarray] = []
        for k, idx in enumerate(self._strip_node_indices):
            if len(idx) == 0:
                rings.append(idx)
                continue
            anchor = int(np.argmin(np.abs(self._mesh_span[idx] - self._strip_centres[k])))
            group = next(g for g in ring_groups[k] if anchor in g)
            rings.append(idx[np.sort(group)])
        return rings

    def _section_rotation(
        self,
        displacements: np.ndarray,
        ring_idx: np.ndarray,
    ) -> float | None:
        r"""Rotation of one physical ring about the span axis [rad].

        The ring's in-plane field is fitted by least squares to ``u = A x + t`` in an
        orthonormal in-plane basis ``(e1, e2)`` with ``e2 = span x e1``; the antisymmetric
        part ``0.5 (a21 - a12)`` is the section rotation, invariant to the in-plane origin
        and with the same sense as the shell's rotational DOF about the span axis.  Returns
        ``None`` when the ring cannot define a 2-D section (fewer than three nodes, or a
        rank-deficient design); that strip then keeps its reference twist, so no second
        estimator can silently disagree with this one.
        """
        if len(ring_idx) < 3:
            return None

        seed = np.array([1.0, 0.0, 0.0])
        if abs(seed @ self._span_dir) > 0.9:  # the seed must not be parallel to the span
            seed = np.array([0.0, 1.0, 0.0])
        e1 = seed - (seed @ self._span_dir) * self._span_dir
        e1 = e1 / np.linalg.norm(e1)
        e2 = np.cross(self._span_dir, e1)

        pts = self._ref_coords[ring_idx]
        x, y = pts @ e1, pts @ e2
        disp = displacements[ring_idx]
        ux, uy = disp @ e1, disp @ e2
        design = np.column_stack([x, y, np.ones_like(x)])
        if np.linalg.matrix_rank(design) < 3:
            return None
        a12 = float(np.linalg.lstsq(design, ux, rcond=None)[0][1])
        a21 = float(np.linalg.lstsq(design, uy, rcond=None)[0][0])
        return 0.5 * (a21 - a12)

    # -----------------------------------------------------------------------
    # BEM and projector reconstruction on deformed geometry
    # -----------------------------------------------------------------------

    def _rebuild_bem_solver(
        self,
        r_def: np.ndarray,
        twist_def: np.ndarray,
    ) -> tuple[BEMSolver, BladeAero]:
        r"""Construct a new ``BEMSolver`` with deformed strip radii and twist.

        Creates a new set of ``AeroStation`` objects with:

        * :math:`r_k \leftarrow r_k^{\mathrm{def}}`
        * :math:`\theta_k \leftarrow \theta_k^{\mathrm{def}}`
        * ``chord``, ``pitch_axis``, ``airfoil`` — unchanged from reference.

        A new ``BladeAero`` is assembled from these stations and passed to
        ``BEMSolver``, which internally rebuilds the ``CCBlade`` rotor object
        with the updated geometry.

        Parameters
        ----------
        r_def : ndarray, shape (n_strips,)
            Deformed radial positions [m].
        twist_def : ndarray, shape (n_strips,)
            Deformed aerodynamic twist [rad].

        Returns
        -------
        solver : BEMSolver
            Fresh solver instance on the deformed geometry.
        deformed_aero : BladeAero
            The deformed blade definition, for callers that need the deformed
            station radii.  The force projection itself deliberately stays on the
            reference geometry (see ``_compute_forces``).
        """
        ba = self._blade_aero

        # The root station may land a few micrometres inside the hub (the second-order
        # span shortening exceeds its first-order stretch).  CCBlade's hub-loss factor
        # exp(-B/2 (r - Rhub) / (Rhub sin phi)) then turns positive-exponent, the residual
        # is NaN at the lower phi bracket and the strip's load is silently zeroed
        # ("NaNs at 0/50").  A station never sits inside the hub, so it is held at Rhub,
        # which is where the rigid path puts the root and carries a finite load.
        new_stations = [
            AeroStation(
                span_fraction=st.span_fraction,
                r=max(float(r_def[k]), float(ba.hub_radius)),
                chord=st.chord,
                twist=float(twist_def[k]),
                pitch_axis=st.pitch_axis,
                airfoil=st.airfoil,
            )
            for k, st in enumerate(ba.stations)
        ]

        # Move Rtip with the tip station, so the tip gap Rtip - r_tip keeps its
        # undeformed value.  AeroDyn fixes its tip-loss constant at initialisation
        # from the undeformed geometry (BEMT.f90, B (zTip - zLocal) / (2 zLocal)),
        # so a deflected blade keeps its undeformed tip distance there.  A fixed
        # Rtip would instead widen (tip moving inward) or close (outward) the gap
        # of the tip strip, where the Prandtl factor is most sensitive to it.  The
        # deck's tip station sits at Rtip (gap 0, so F = 0 there on the rigid path
        # too); the max() only keeps every station at or inside Rtip, because a
        # station past it gives Prandtl factortip < 0 -> exp(+) > 1 -> acos(NaN).
        # It takes no buffer: a buffer would open a tip gap the rigid path does not have.
        tip = int(np.argmax(self._ref_r))
        rotor_radius = max(
            ba.rotor_radius + float(r_def[tip] - self._ref_r[tip]),
            float(np.max(r_def)),
        )

        deformed_aero = BladeAero(
            airfoils=ba.airfoils,
            stations=new_stations,
            blade_length=ba.blade_length,
            hub_radius=ba.hub_radius,
            rotor_radius=rotor_radius,
            n_blades=ba.n_blades,
        )

        return BEMSolver(deformed_aero, **self._bem_solver_kwargs), deformed_aero

    def _compute_forces(
        self,
        displacements: np.ndarray,
    ) -> tuple[np.ndarray, BEMResult]:
        r"""Evaluate BEM on the (possibly deformed) blade and project forces.

        **Pipeline:**

        1. If :math:`\max_i \|\mathbf{u}_i\| < 10^{-12}` (effectively zero
           displacement), skip geometry reconstruction and evaluate on the
           reference configuration.  This avoids unnecessary SVD and CCBlade
           rebuilds at the start of the simulation or on converged windows.

        2. Otherwise:

           a. ``_compute_deformed_geometry`` → :math:`(r^{\mathrm{def}},
              \theta^{\mathrm{def}})`.
           b. ``_rebuild_bem_solver`` → new ``BEMSolver`` with deformed
              geometry.
           c. ``BEMSolver.compute`` → distributed loads :math:`N_p, T_p`.
           d. ``ForceProjector.project`` → nodal force array, on the
              **reference** projector (see the note in the body).

        Parameters
        ----------
        displacements : ndarray, shape (n_nodes, 3)
            Nodal displacements from the structural solver.

        Returns
        -------
        forces : ndarray, shape (n_nodes, 3)
            Aerodynamic nodal forces for the current deformed BEM state, laid out
            with the reference projection geometry.
        bem_result : BEMResult
            Full BEM output (Np, Tp, alpha, Cl, Cd, integrated loads).
        """
        v_inf = float(self._bem_cfg.get("wind_speed", 45.0))
        # _current_omega is in rad/s (received from Solid via preCICE or YAML).
        # CCBlade expects RPM.
        import math as _math  # noqa: PLC0415

        omega = self._current_omega * 60.0 / (2.0 * _math.pi)
        pitch = float(self._bem_cfg.get("pitch", 0.0))
        # Azimuth from the structure's rotor angle (end of window), not integrated here.
        azimuth = self._azimuth_deg()

        disp_max = float(np.max(np.linalg.norm(displacements, axis=1)))

        if disp_max < 1e-12:
            # Zero displacement — use pre-built reference solver and projector
            bem_result = self._bem_solver.compute(v_inf, omega, pitch, azimuth=azimuth)
            forces = self._projector.project(bem_result)
            return forces, bem_result

        # -- Deformed geometry pipeline ------------------------------------
        # The BEM geometry follows the deformation; the *projection* geometry does
        # not.  Np/Tp are evaluated on the deformed annulus and twist, but the strip
        # grid, the aerodynamic-centre arm and the pitching-moment axis stay on the
        # reference mesh.  They are properties of the blade's reference geometry,
        # and re-deriving them from the deformed surface on every sub-iteration adds
        # a deformation-dependent gain to the interface that the implicit coupling
        # cannot absorb: with the rebuild HEAD converges 39 of 500 windows, and
        # keeping the reference projector converges 500 of 500.  Pinning its
        # geometry outputs one at a time does not restore contraction - no single
        # output is the lever, their sum is (odd/tasks/coupled-divergence.md, #19).
        r_def, twist_def = self._compute_deformed_geometry(displacements)
        bem_solver, _ = self._rebuild_bem_solver(r_def, twist_def)

        # _strip_node_indices is intentionally kept as the *reference* assignment and
        # is NOT updated from the deformed geometry.  Updating it would move the
        # strip's centre ring between iterations, so the same physical section would
        # be measured against a different node set from one implicit iteration to the
        # next.  The force application now follows the same rule, which is why it
        # uses the reference projector as well.

        bem_result = bem_solver.compute(v_inf, omega, pitch, azimuth=azimuth)
        forces = self._projector.project(bem_result)
        return forces, bem_result

    # -----------------------------------------------------------------------
    # Diagnostics
    # -----------------------------------------------------------------------

    def _write_timestep_output(
        self,
        forces: np.ndarray,
        displacements: np.ndarray,
        bem_result: BEMResult,
        step: int,
        time: float,
    ) -> None:
        """Write all output for one converged time window.

        Structure mirrors the solid solver::

            {output_folder}/
              results.pvd               ← PVD time collection for fields.vtu
              bem_report.csv            ← one row per window (global loads)
              {time}/
                fields.vtu              ← surface mesh: Force + Displacement
                bem_sections.vtu        ← strip centroids: BEM spanwise data
                bem_sectional.csv       ← sectional loads table
        """
        time_str = f"{time:.6g}"
        ts_dir = self._output_folder / time_str
        ts_dir.mkdir(parents=True, exist_ok=True)

        self._write_fields_vtu(forces, displacements, ts_dir, time_str, time)
        self._write_sections_vtu(bem_result, ts_dir)
        self._write_sectional_csv(bem_result, ts_dir)
        self._append_global_csv(bem_result, forces, displacements, step, time, time_str)

    def _write_fields_vtu(
        self,
        forces: np.ndarray,
        displacements: np.ndarray,
        ts_dir: "Path",
        time_str: str,
        time: float,
    ) -> None:
        """Write surface VTU with Force and Displacement nodal fields."""
        try:
            import meshio
        except ImportError:
            return

        viz = self._viz_mesh
        coords = viz.coords_array  # (N_viz, 3)

        # Build node-id → index map for the viz mesh
        node_id_to_viz_idx = {node.id: i for i, node in enumerate(viz.nodes)}

        # Map coupling-mesh forces/displacements onto the viz mesh
        # (coupling nodes are a subset of viz nodes — same Node objects)
        n_viz = len(viz.nodes)
        forces_viz = np.zeros((n_viz, 3))
        disp_viz = np.zeros((n_viz, 3))
        for i, node in enumerate(self._working_mesh.nodes):
            j = node_id_to_viz_idx.get(node.id)
            if j is not None:
                forces_viz[j] = forces[i]
                disp_viz[j] = displacements[i]

        # Build cell connectivity grouped by type
        cell_blocks: dict[str, list] = {}
        for elem in viz.elements:
            conn = [node_id_to_viz_idx[n.id] for n in elem.nodes]
            n_nds = len(conn)
            if n_nds == 3:
                ctype = "triangle"
            elif n_nds == 4:
                ctype = "quad"
            elif n_nds == 6:
                ctype = "triangle6"
            elif n_nds == 8:
                ctype = "quad8"
            else:
                continue
            cell_blocks.setdefault(ctype, []).append(conn)

        meshio_cells = [(k, np.array(v)) for k, v in cell_blocks.items()]
        if not meshio_cells:
            # Fallback to point cloud if no elements
            meshio_cells = [("vertex", np.arange(n_viz).reshape(-1, 1))]

        m = meshio.Mesh(
            points=coords,
            cells=meshio_cells,  # pyright: ignore[reportArgumentType]
            point_data={"Force": forces_viz, "Displacement": disp_viz},
        )
        vtu_path = ts_dir / "fields.vtu"
        m.write(str(vtu_path))

        # Update PVD with relative path
        rel_path = f"{ts_dir.name}/fields.vtu"
        self._update_pvd(self._output_folder, "results.pvd", rel_path, time)

    def _write_sections_vtu(self, bem_result: BEMResult, ts_dir: "Path") -> None:
        """Write strip centroids VTU with spanwise BEM aerodynamic data."""
        try:
            import meshio
        except ImportError:
            return

        strips = self._projector._strips
        n_strips = len(strips)

        centroids = np.array([s.centroid for s in strips])
        dr = np.array([s.dr for s in strips])

        # Integrated force vectors at each strip centroid
        F_section = np.zeros((n_strips, 3))
        for k, strip in enumerate(strips):
            F_n = float(bem_result.Np[k]) * strip.dr
            F_t = float(bem_result.Tp[k]) * strip.dr
            F_section[k] = F_n * self._normal_dir + F_t * self._tangential_dir

        point_data: dict[str, np.ndarray] = {
            "r_m": bem_result.r.astype(float),
            "Np_N_m": bem_result.Np.astype(float),
            "Tp_N_m": bem_result.Tp.astype(float),
            "alpha_deg": bem_result.alpha.astype(float),
            "cl": bem_result.cl.astype(float),
            "cd": bem_result.cd.astype(float),
            "a": bem_result.a.astype(float),
            "ap": bem_result.ap.astype(float),
            "dr_m": dr.astype(float),
            "Force_section": F_section,
        }
        # Optional fields — present when CCBlade version supports them
        if bem_result.cn is not None:
            point_data["cn"] = bem_result.cn.astype(float)
        if bem_result.ct is not None:
            point_data["ct"] = bem_result.ct.astype(float)
        if bem_result.W is not None:
            point_data["W_m_s"] = bem_result.W.astype(float)
        if bem_result.Re is not None:
            point_data["Re"] = bem_result.Re.astype(float)
        if bem_result.cm is not None:
            point_data["cm"] = bem_result.cm.astype(float)
        if bem_result.Mp is not None:
            point_data["Mp_Nm_m"] = bem_result.Mp.astype(float)
        if bem_result.chord is not None:
            point_data["chord_m"] = bem_result.chord.astype(float)
        if bem_result.twist_deg is not None:
            point_data["twist_deg"] = bem_result.twist_deg.astype(float)

        cells = [("vertex", np.arange(n_strips).reshape(-1, 1))]
        # meshio's stubs want ArrayLike and ndarray is one, but they are invariant, so
        # pyright reports the arguments.  The suppression has to sit on the reported
        # lines; on the closing paren it does nothing (measured at HEAD).
        m = meshio.Mesh(
            points=centroids,
            cells=cells,  # pyright: ignore[reportArgumentType]
            point_data=point_data,  # pyright: ignore[reportArgumentType]
        )
        m.write(str(ts_dir / "bem_sections.vtu"))

    def _write_sectional_csv(self, bem_result: BEMResult, ts_dir: "Path") -> None:
        """Write sectional loads CSV for this timestep."""
        strips = self._projector._strips
        dr = np.array([s.dr for s in strips])

        names = [
            "r[m]",
            "dr[m]",
            "chord[m]",
            "twist[deg]",
            "Np[N/m]",
            "Tp[N/m]",
            "alpha[deg]",
            "cl",
            "cd",
        ]
        cols = [
            bem_result.r,
            dr,
            bem_result.chord if bem_result.chord is not None else np.zeros_like(bem_result.r),
            bem_result.twist_deg
            if bem_result.twist_deg is not None
            else np.zeros_like(bem_result.r),
            bem_result.Np,
            bem_result.Tp,
            bem_result.alpha,
            bem_result.cl,
            bem_result.cd,
        ]
        for attr, col_name in [
            ("cn", "cn"),
            ("ct", "ct"),
            ("a", "a"),
            ("ap", "ap"),
            ("W", "W[m/s]"),
            ("Re", "Re"),
            ("cm", "cm"),
            ("Mp", "Mp[N.m/m]"),
        ]:
            val = getattr(bem_result, attr, None)
            if val is not None:
                names.append(col_name)
                cols.append(val)
        data = np.column_stack(cols)
        np.savetxt(
            ts_dir / "bem_sectional.csv",
            data,
            delimiter=",",
            header=",".join(names),
            comments="",
        )

    def _append_global_csv(
        self,
        bem_result: BEMResult,
        forces: np.ndarray,
        displacements: np.ndarray,
        step: int,
        time: float,
        time_str: str = "",
    ) -> None:
        """Append one row to the accumulated global loads CSV (like structural_report.csv)."""
        csv_path = self._output_folder / "bem_report.csv"
        write_header = not csv_path.exists()

        disp_mag = np.linalg.norm(displacements, axis=1)
        force_mag = np.linalg.norm(forces, axis=1)

        # Tip columns are ROTATING-frame (blade frame) displacements: with
        # the read data are always un-rotated first (T7a), so they no longer
        # match the raw preCICE (inertial) vectors.
        tip_vec = displacements[self._tip_node_idx] if len(displacements) else np.zeros(3)
        tip_disp_x, tip_disp_y, tip_disp_z = float(tip_vec[0]), float(tip_vec[1]), float(tip_vec[2])
        tip_disp_mag = float(np.linalg.norm(tip_vec))

        max_disp = float(np.max(disp_mag)) if len(disp_mag) else 0.0
        max_force = float(np.max(force_mag)) if len(force_mag) else 0.0

        CP = bem_result.CP if bem_result.CP is not None else float("nan")
        CT = bem_result.CT if bem_result.CT is not None else float("nan")
        CQ = bem_result.CQ if bem_result.CQ is not None else float("nan")
        Mb = bem_result.Mb if bem_result.Mb is not None else float("nan")

        with open(csv_path, "a") as f:
            if write_header:
                f.write(
                    "Time [s],Step,"
                    "Thrust [N],Torque [N.m],Power [W],"
                    "CT,CQ,CP,"
                    "Mb [N.m],"
                    "Max Disp [m],"
                    "Tip Disp X [m],Tip Disp Y [m],Tip Disp Z [m],Tip Disp Mag [m],"
                    "Max Nodal Force [N]\n"
                )
            f.write(
                f"{time:.6f},{step},"
                f"{bem_result.thrust:.6e},{bem_result.torque:.6e},{bem_result.power:.6e},"
                f"{CT:.6f},{CQ:.6f},{CP:.6f},"
                f"{Mb:.6e},"
                f"{max_disp:.6e},"
                f"{tip_disp_x:.6e},{tip_disp_y:.6e},{tip_disp_z:.6e},{tip_disp_mag:.6e},"
                f"{max_force:.6e}\n"
            )

    def _update_pvd(self, folder: "Path", pvd_name: str, filename: str, time: float) -> None:
        """Append a timestep entry to the PVD collection file."""
        pvd_path = folder / pvd_name

        header = (
            '<?xml version="1.0"?>\n'
            '<VTKFile type="Collection" version="0.1" byte_order="LittleEndian"'
            ' compressor="vtkZLibDataCompressor">\n'
            "  <Collection>\n"
        )
        footer = "  </Collection>\n</VTKFile>"
        entry = f'    <DataSet timestep="{time}" group="" part="0" file="{filename}"/>\n'

        if not pvd_path.exists():
            with open(pvd_path, "w") as f:
                f.write(header + entry + footer)
        else:
            with open(pvd_path, "r") as f:
                lines = f.readlines()
            valid_lines = [
                line for line in lines if "</Collection>" not in line and "</VTKFile>" not in line
            ]
            with open(pvd_path, "w") as f:
                f.writelines(valid_lines)
                f.write(entry)
                f.write(footer)

    def _log_window(
        self,
        bem_result: BEMResult,
        displacements: np.ndarray,
        current_time: float,
    ) -> None:
        """Log a diagnostic summary for this time window."""
        if len(displacements) == 0:
            tip_vec = np.zeros(3)
            max_disp = 0.0
        else:
            disp_mag = np.linalg.norm(displacements, axis=1)
            tip_vec = displacements[self._tip_node_idx]
            max_disp = float(np.max(disp_mag))
        tip_x, tip_y, tip_z = float(tip_vec[0]), float(tip_vec[1]), float(tip_vec[2])
        tip_mag = float(np.linalg.norm(tip_vec))
        logger.info(
            "[BEM-FSI] t=%.4f s | window=%d | iters=%d | T=%.2f kN "
            "| Q=%.2f kNm | P=%.2f kW "
            "| tip=(%.4f, %.4f, %.4f) mag=%.4f m | max_disp=%.4f m",
            current_time,
            self._window_count,
            self._iteration_count,
            bem_result.thrust * 1e-3,
            bem_result.torque * 1e-3,
            bem_result.power * 1e-3,
            tip_x,
            tip_y,
            tip_z,
            tip_mag,
            max_disp,
        )


# ---------------------------------------------------------------------------
# Factory helper used by the CLI
# ---------------------------------------------------------------------------


def build_from_config(
    mesh: MeshModel,
    cfg: dict,
    config_file: str | Path = "precice-config.xml",
    viz_mesh: MeshModel | None = None,
    element_properties: dict | None = None,
) -> "BEMFSIParticipant":
    """Construct a :class:`BEMFSIParticipant` from a YAML config dict.

    The dict is expected to match the schema documented in
    ``docs/cli-reference.md`` for ``aeroelast-bem-fsi``.

    Parameters
    ----------
    mesh : MeshModel
        Pre-built blade surface mesh (caller's responsibility).
    cfg : dict
        Parsed YAML content.
    config_file : str or Path
        Override for the preCICE XML path (takes precedence over
        ``cfg["config_file"]`` when explicitly provided).
    viz_mesh : MeshModel or None
        Full mesh kept for surface VTU output (optional).
    element_properties : dict or None
        The deck's section-property map for the coupling mesh, forwarded to the
        participant's ``ForceProjector``; ``None`` keeps the minimum-norm
        ``_distribute`` fallback.
    """
    bem_cfg = cfg.get("bem", {})
    output_cfg = cfg.get("output", {})

    blade_file = cfg.get("blade_file")
    if blade_file is None:
        raise ValueError("'blade_file' is required in the BEM-FSI configuration.")

    # omega can be a scalar (constant) or a dict with initial + preCICE config:
    #   omega: 0.7917                        # constant
    #   omega:
    #     initial: 0.7917                    # used until first preCICE window
    #     mesh:   "GlobalFluidMesh"
    #     data:   "AngularVelocity"
    #     vertex: [0.0, 0.0, 0.0]
    omega_cfg = bem_cfg.get("omega", 0.0)
    if isinstance(omega_cfg, dict):
        omega_initial = float(omega_cfg.get("initial", 0.0))
        omega_mesh: str | None = omega_cfg.get("mesh")
        omega_data: str = omega_cfg.get("data", "AngularVelocity")
        omega_vertex: list[float] = omega_cfg.get("vertex", [0.0, 0.0, 0.0])
        theta_data: str = omega_cfg.get("theta_data", "Theta")
        # Normalise bem_cfg so the constructor can read omega as a float
        bem_cfg = {**bem_cfg, "omega": omega_initial}
    else:
        omega_initial = float(omega_cfg)
        omega_mesh = None
        omega_data = "AngularVelocity"
        omega_vertex = [0.0, 0.0, 0.0]
        theta_data = "Theta"
    _check_theta_source(omega_mesh, omega_initial)

    blade_aero = load_blade_aero(
        blade_file,
        default_re=float(bem_cfg.get("default_re", 1e7)),
        neuralfoil_model=bem_cfg.get("neuralfoil_model", "large"),
        hub_radius=float(bem_cfg.get("hub_radius", 0.0)),
        n_blades=int(bem_cfg.get("n_blades", 3)),
        viterna_ar=float(bem_cfg.get("viterna_ar", 17.0)),
        viterna_confidence_threshold=float(bem_cfg.get("viterna_confidence_threshold", 0.5)),
    )

    return BEMFSIParticipant(
        mesh=mesh,
        blade_aero=blade_aero,
        bem_config=bem_cfg,
        participant=cfg.get("participant", "Fluid"),
        config_file=cfg.get("config_file", str(config_file)),
        coupling_mesh=cfg.get("coupling_mesh", "Fluid-Mesh"),
        displacement_data=cfg.get("displacement_data", "Displacement"),
        force_data=cfg.get("force_data", "Force"),
        output_folder=output_cfg.get("folder", "bem_fsi_results"),
        log_interval=int(output_cfg.get("log_interval", 10)),
        viz_mesh=viz_mesh,
        omega_mesh=omega_mesh,
        omega_data=omega_data,
        omega_vertex=omega_vertex,
        element_properties=element_properties,
        theta_data=theta_data,
    )
