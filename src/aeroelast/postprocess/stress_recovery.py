"""
Stress and Strain Recovery Module for Shell Elements.

This module recovers nodal stress and strain fields from a displacement
solution vector and computes the derived engineering quantities (Von Mises,
principal stresses, maximum shear) that are standard outputs in commercial
FEA post-processors such as ANSYS, Abaqus and Nastran.

Shell elements (MITC3, MITC4) are supported — Reissner–Mindlin
plate/shell theory, plane-stress assumption (σ_zz = 0).  Three components
are recovered: σ_xx, σ_yy, σ_xy.  The through-thickness variation is
captured by evaluating at TOP (+h/2), MIDDLE (0) and BOTTOM (−h/2)
surfaces.

Shell Theory
------------
For Reissner–Mindlin shells, the in-plane stress at through-thickness
coordinate *z* measured from the mid-surface is:

    σ(z) = σ_m + z · κ_b

where:

* σ_m = (C / h) · (B_m · u_e)  — **membrane stress** (constant through
  thickness).  B_m is the membrane strain–displacement matrix and C is the
  plane-stress constitutive matrix (E, ν).  Note that the element's
  ``Cm()`` method returns the *integrated* membrane stiffness D = C · h;
  we divide by h to recover the actual material stress–strain matrix.
* κ_b = B_κ · u_e — **curvature** from the bending strain–displacement
  matrix.  The bending stress contribution is (C / h) · z · κ_b, varying
  linearly through the thickness.
* h — total shell thickness.

The stress is evaluated at the element's **parametric node coordinates**
(natural coordinates of the element nodes), then contributions from
adjacent elements are averaged at shared nodes (SPR-style nodal smoothing).

Derived Quantities
------------------
Von Mises equivalent stress (3-D general form, reduces to the plane-stress
formula when σ_zz = τ_yz = τ_zx = 0):

    σ_vm = √( ½[(σ_xx−σ_yy)² + (σ_yy−σ_zz)² + (σ_zz−σ_xx)²]
              + 3[τ_xy² + τ_yz² + τ_zx²] )

Principal stresses:

* **Shells:** Mohr's circle  —  σ₁,₂ = (σ_xx+σ_yy)/2
  ± √[(σ_xx−σ_yy)²/4 + τ_xy²] ;  τ_max = (σ₁ − σ₂) / 2.

References
----------
- Bathe, K.J. (2014). *Finite Element Procedures*, 2nd Edition.
  Chapters 5 (isoparametric elements) and 6 (shell elements).
- Cook, R.D., Malkus, D.S., Plesha, M.E., Witt, R.J. (2002).
  *Concepts and Applications of Finite Element Analysis*, 4th Edition.
- Hinton, E. & Campbell, J.S. (1974). "Local and Global Smoothing of
  Discontinuous Finite Element Functions Using a Least Squares Method",
  Int. J. Num. Meth. Eng. 8, 461–480.  (Gauss-to-node extrapolation.)
- Zienkiewicz, O.C. & Zhu, J.Z. (1992). "The Superconvergent Patch
  Recovery and a posteriori error estimates", Int. J. Num. Meth. Eng.
  33, 1331–1364.  (Theoretical basis for nodal stress smoothing.)
"""

from dataclasses import dataclass
from enum import Enum
from typing import TYPE_CHECKING, Dict, Tuple

import numpy as np

if TYPE_CHECKING:
    from aeroelast.core.assembler import MeshAssembler


class StressLocation(Enum):
    """Through-thickness coordinate for shell stress evaluation.

    In Reissner–Mindlin shell theory the stress varies linearly through
    the thickness *h*.  The through-thickness coordinate *z* is measured
    from the mid-surface:

    * ``TOP``    — z = +h/2 (outer fibre, tension under positive bending).
    * ``MIDDLE`` — z = 0    (mid-surface, membrane stress only).
    * ``BOTTOM`` — z = −h/2 (inner fibre, compression under positive bending).
    """

    TOP = "top"
    MIDDLE = "middle"
    BOTTOM = "bottom"


class StressType(Enum):
    """Stress decomposition mode for shell elements.

    Shell stresses can be decomposed into two physically distinct
    contributions:

    * ``MEMBRANE`` — σ_m = (C/h) · B_m · u_e.  Constant through the
      thickness; arises from in-plane stretching/compression.
    * ``BENDING``  — σ_b(z) = (C/h) · z · B_κ · u_e.  Varies linearly;
      arises from curvature (plate bending).
    * ``TOTAL``    — σ = σ_m + σ_b(z).  Combined stress at the
      requested ``StressLocation``.
    """

    MEMBRANE = "membrane"
    BENDING = "bending"
    TOTAL = "total"


@dataclass
class StressResult:
    """
    Container for stress computation results.

    Plane-stress shell results: only ``sigma_xx``, ``sigma_yy`` and
    ``sigma_xy`` are non-zero.

    Attributes
    ----------
    sigma_xx, sigma_yy, sigma_xy : np.ndarray
        In-plane stress components.
    von_mises : np.ndarray
        Von Mises equivalent stress (general 3-D formula).
    sigma_1, sigma_2 : np.ndarray
        Maximum / minimum principal stresses (2-D Mohr circle).
    tau_max : np.ndarray
        Maximum shear stress.
    principal_angle : np.ndarray
        Angle (rad) of first principal direction from x-axis.
    """

    sigma_xx: np.ndarray
    sigma_yy: np.ndarray
    sigma_xy: np.ndarray
    von_mises: np.ndarray
    sigma_1: np.ndarray
    sigma_2: np.ndarray
    tau_max: np.ndarray
    principal_angle: np.ndarray

    def to_dict(self) -> Dict[str, np.ndarray]:
        """Convert to dictionary for VTK output."""
        d: Dict[str, np.ndarray] = {
            "sigma_xx": self.sigma_xx,
            "sigma_yy": self.sigma_yy,
            "sigma_xy": self.sigma_xy,
            "von_mises": self.von_mises,
            "sigma_1": self.sigma_1,
            "sigma_2": self.sigma_2,
            "tau_max": self.tau_max,
            "principal_angle": np.degrees(self.principal_angle),
        }
        return d


@dataclass
class StrainResult:
    """
    Container for strain computation results.

    Attributes
    ----------
    epsilon_xx, epsilon_yy : np.ndarray
        Normal strains.
    gamma_xy : np.ndarray
        In-plane engineering shear strain (= 2·ε_xy).
    epsilon_1, epsilon_2 : np.ndarray
        Principal strains (max / min).
    gamma_max : np.ndarray
        Maximum shear strain.
    """

    epsilon_xx: np.ndarray
    epsilon_yy: np.ndarray
    gamma_xy: np.ndarray
    epsilon_1: np.ndarray
    epsilon_2: np.ndarray
    gamma_max: np.ndarray

    def to_dict(self) -> Dict[str, np.ndarray]:
        """Convert to dictionary for VTK output."""
        d: Dict[str, np.ndarray] = {
            "epsilon_xx": self.epsilon_xx,
            "epsilon_yy": self.epsilon_yy,
            "gamma_xy": self.gamma_xy,
            "epsilon_1": self.epsilon_1,
            "epsilon_2": self.epsilon_2,
            "gamma_max": self.gamma_max,
        }
        return d


class StressRecovery:
    """
    Stress and strain recovery engine for shell finite elements.

    This class takes a converged displacement solution and recovers the
    complete Cauchy stress tensor at nodes or element centres, together
    with the standard derived engineering quantities (Von Mises, principal
    stresses, maximum shear) used for structural assessment.

    Workflow
    --------
    1. **Instantiation** — store a reference to the ``MeshAssembler``
       (which owns the element map, connectivity and constitutive data)
       and a copy of the displacement vector.
    2. **Stress evaluation** — for each element the displacement DOFs are
       gathered and the stress is computed:

       * *Shell elements (MITC3, MITC4)*:  Membrane and bending
         strain–displacement matrices ``B_m``, ``B_κ`` produce a plane-
         stress triplet [σ_xx, σ_yy, τ_xy] at a chosen parametric
         point (r, s) and through-thickness location z.

    3. **Nodal smoothing** — stresses from adjacent elements are averaged
       at shared nodes.

    Supported Element Topologies
    ----------------------------
    ============  =====  =======  ===================================
    Element       Nodes  GP count Quadrature rule
    ============  =====  =======  ===================================
    MITC3           3       3     Hammer (triangle)
    MITC4           4       4     2×2 Gauss–Legendre
    ============  =====  =======  ===================================

    Parameters
    ----------
    domain : MeshAssembler
        Mesh assembler that owns the element map, node coordinates,
        constitutive properties and DOF connectivity.
    u : np.ndarray or PETSc.Vec
        Full (unreduced) displacement solution vector.  It contains
        6 DOFs per node [u, v, w, θ_x, θ_y, θ_z].

    Attributes
    ----------
    n_nodes : int
        Total number of mesh nodes.
    n_elements : int
        Total number of elements in the mesh.
    dofs_per_node : int
        Number of DOFs per node (6 for shells).
    """

    def __init__(self, domain: "MeshAssembler", u):
        self.domain = domain
        if hasattr(u, "array"):
            self.u = u.array.copy()
        else:
            self.u = np.asarray(u).copy()
        self.dofs_per_node = domain.dofs_per_node
        self.n_nodes = len(domain.coords_array)
        _rust = getattr(domain, "_rust", None)
        if _rust is None:
            raise RuntimeError("StressRecovery requires a live Rust assembler (domain._rust).")
        self.n_elements = _rust.n_elems
        self._node_id_to_index = domain.node_id_to_index

    # ------------------------------------------------------------------
    # Element-level stress  (centroid / single point)
    # ------------------------------------------------------------------
    def compute_element_stresses(
        self,
        location: StressLocation = StressLocation.MIDDLE,
        stress_type: StressType = StressType.TOTAL,
        gauss_point: Tuple[float, float] = (0.0, 0.0),
    ) -> StressResult:
        """Compute a single representative stress per element (centroid value).

        This method provides one stress state per element, suitable for
        contour plots at element centres ("element results" in commercial
        codes).

        Stress is evaluated at the parametric point *gauss_point*
        (default: element centre) and through-thickness location *z*
        (from *location*).

        Parameters
        ----------
        location : StressLocation, default MIDDLE
            Through-thickness position (shells only).
        stress_type : StressType, default TOTAL
            Which shell stress contribution to include.
        gauss_point : tuple of float, default (0.0, 0.0)
            Parametric coordinates at which to evaluate shell stress.

        Returns
        -------
        StressResult
            One stress state per element.  Only the in-plane components
            (σ_xx, σ_yy, τ_xy) are populated.
        """
        r0, s0 = gauss_point

        # ------------------------------------------------------------------
        # Fast path: delegate to the Rust assembler when available and the
        # caller requests the canonical centroid point (0, 0).  The Rust
        # kernel evaluates exactly at the element centroid so it matches the
        # default gauss_point.  Non-default gauss_point values fall through
        _rust = getattr(self.domain, "_rust", None)
        if _rust is not None and r0 == 0.0 and s0 == 0.0:
            _Z_FACTOR = {
                StressLocation.BOTTOM: -0.5,
                StressLocation.MIDDLE: 0.0,
                StressLocation.TOP: +0.5,
            }
            _STRESS_TYPE = {
                StressType.MEMBRANE: 0,
                StressType.BENDING: 1,
                StressType.TOTAL: 2,
            }
            z_factor = _Z_FACTOR.get(location, 0.0)
            stype_int = _STRESS_TYPE.get(stress_type, 2)
            sigma_all, _ = _rust.compute_stress_field(self.u, z_factor, stype_int)
            return self._build_stress_result(sigma_all)

        raise NotImplementedError(
            "compute_element_stresses: non-centroid gauss_point evaluation requires "
            "r0=0.0, s0=0.0 and a live Rust assembler. "
            "The Python per-element fallback has been removed."
        )

    # ------------------------------------------------------------------
    # Nodal-averaged stress (the primary output)
    # ------------------------------------------------------------------
    def compute_nodal_stresses(
        self,
        location: StressLocation = StressLocation.MIDDLE,
        stress_type: StressType = StressType.TOTAL,
        smoothing: str = "average",
    ) -> StressResult:
        """Compute smoothed (continuous) nodal stress fields.

        This is the **primary output** method.  It produces a nodal
        stress field that can be passed directly to a VTU writer for
        visualisation.  The smoothing mimics commercial FEA practice
        (e.g. ANSYS “Nodeal Solution”, Abaqus “field output at nodes”).

        Procedure
        ---------
        * **Shell elements** — the stress is evaluated at each element's
          **parametric node coordinates** (the natural-coordinate location
          of each node).  Each element contributes its own stress value
          to the global node; contributions from all surrounding elements
          are then averaged (or area-weighted).

        Parameters
        ----------
        location : StressLocation, default MIDDLE
            Through-thickness position for shell stress evaluation.
        stress_type : StressType, default TOTAL
            Membrane, bending or total (shells only).
        smoothing : str, default ``"average"``
            Smoothing mode.  ``"average"`` gives equal weight to every
            contributing element; ``"area_weighted"`` weights by element
            area (shells).

        Returns
        -------
        StressResult
            Nodal stress field with shape ``(n_nodes,)`` for each
            component.
        """
        n_nodes = self.n_nodes

        # ------------------------------------------------------------------
        # Fast path: use Rust for element stress computation; averaging loop
        # still runs in Python but avoids the costly per-element Python
        # stress evaluation (the hot path in the original code).
        # ------------------------------------------------------------------
        _rust = getattr(self.domain, "_rust", None)
        if _rust is not None:
            _Z_FACTOR_NS = {
                StressLocation.BOTTOM: -0.5,
                StressLocation.MIDDLE: 0.0,
                StressLocation.TOP: +0.5,
            }
            _STRESS_TYPE_NS = {
                StressType.MEMBRANE: 0,
                StressType.BENDING: 1,
                StressType.TOTAL: 2,
            }
            z_factor = _Z_FACTOR_NS.get(location, 0.0)
            stype_int = _STRESS_TYPE_NS.get(stress_type, 2)
            sigma_elem, _ = _rust.compute_stress_field(self.u, z_factor, stype_int)
            # sigma_elem: (n_elems, 6) — one stress state per element

            if smoothing == "area_weighted":
                raise NotImplementedError(
                    "compute_nodal_stresses: area_weighted smoothing requires "
                    "Python element objects which have been removed. "
                    "Use smoothing='average' instead."
                )

            sigma_sum = np.zeros((n_nodes, 6))
            weight_sum = np.zeros(n_nodes)

            for elem_idx, mesh_elem in enumerate(self.domain.elements):
                node_ids = mesh_elem.node_ids
                for node_id in node_ids:
                    gn = self._node_id_to_index[node_id]
                    sigma_sum[gn, :] += sigma_elem[elem_idx, :]
                    weight_sum[gn] += 1.0

            mask = weight_sum > 0
            sigma_avg = np.zeros((n_nodes, 6))
            sigma_avg[mask] = sigma_sum[mask] / weight_sum[mask, np.newaxis]
            return self._build_stress_result(sigma_avg)

        raise NotImplementedError(
            "compute_nodal_stresses: requires a live Rust assembler (domain._rust). "
            "The Python per-element fallback has been removed."
        )

    # ------------------------------------------------------------------
    # Multi-layer shell convenience (TOP / MIDDLE / BOTTOM at once)
    # ------------------------------------------------------------------
    def compute_nodal_stresses_all_layers(
        self,
        stress_type: StressType = StressType.TOTAL,
        smoothing: str = "average",
    ) -> Dict[str, StressResult]:
        """Compute nodal stresses at TOP, MIDDLE and BOTTOM shell surfaces.

        For structures with shell elements, stress varies linearly
        through the thickness.  This convenience method evaluates
        ``compute_nodal_stresses`` three times — once at each
        through-thickness location — and returns the results in a dict.

        The three layers correspond to:

        * ``"TOP"`` — z = +h/2 (outer fibre)
        * ``"MID"`` — z = 0   (mid-surface, membrane only)
        * ``"BOT"`` — z = −h/2 (inner fibre)

        Parameters
        ----------
        stress_type : StressType, default TOTAL
            Membrane, bending or combined.
        smoothing : str, default ``"average"``
            Smoothing mode (``"average"`` or ``"area_weighted"``).

        Returns
        -------
        dict of {str: StressResult}
            Keys are ``"TOP"``, ``"MID"``, ``"BOT"``.
        """
        return {
            "TOP": self.compute_nodal_stresses(StressLocation.TOP, stress_type, smoothing),
            "MID": self.compute_nodal_stresses(StressLocation.MIDDLE, stress_type, smoothing),
            "BOT": self.compute_nodal_stresses(StressLocation.BOTTOM, stress_type, smoothing),
        }

    def compute_nodal_stresses_all_layers_dict(
        self,
        stress_type: StressType = StressType.TOTAL,
        smoothing: str = "average",
    ) -> Dict[str, np.ndarray]:
        """Return a flat dictionary of nodal stress arrays for VTU export.

        Calls ``compute_nodal_stresses_all_layers`` internally and
        flattens the per-layer ``StressResult`` objects into a single
        ``dict[str, np.ndarray]`` with prefixed keys that can be passed
        directly as ``extra_fields`` to the checkpoint writer.

        Key naming convention::

            {LAYER}_{component}

        Examples: ``TOP_von_mises``, ``MID_sigma_xx``, ``BOT_tau_max``, etc.

        Parameters
        ----------
        stress_type : StressType, default TOTAL
            Membrane, bending or combined.
        smoothing : str, default ``"average"``
            Smoothing mode.

        Returns
        -------
        dict of {str: np.ndarray}
            Flat dictionary suitable for ``CheckpointManager.write(
            extra_fields=...)``.
        """
        layers = self.compute_nodal_stresses_all_layers(stress_type, smoothing)
        out: Dict[str, np.ndarray] = {}
        for prefix, result in layers.items():
            for key, arr in result.to_dict().items():
                out[f"{prefix}_{key}"] = arr
        return out

    def compute_nodal_strains_all_layers_dict(
        self,
        smoothing: str = "average",
    ) -> Dict[str, np.ndarray]:
        """Return a flat dictionary of nodal strain arrays for VTU export.

        Evaluates ``compute_nodal_strains`` at TOP, MIDDLE and BOTTOM
        through-thickness locations (shells) and flattens the resulting
        ``StrainResult`` objects into a single ``dict[str, np.ndarray]``
        with layer-prefixed keys.

        Key naming convention::

            {LAYER}_{component}

        Examples: ``TOP_epsilon_xx``, ``MID_gamma_xy``, ``BOT_epsilon_1``.

        Parameters
        ----------
        smoothing : str, default ``"average"``
            ``"average"`` or ``"area_weighted"``.

        Returns
        -------
        dict of {str: np.ndarray}
            Flat dictionary suitable for ``CheckpointManager.write(
            extra_fields=...)``.
        """
        out: Dict[str, np.ndarray] = {}
        for prefix, loc in (
            ("TOP", StressLocation.TOP),
            ("MID", StressLocation.MIDDLE),
            ("BOT", StressLocation.BOTTOM),
        ):
            result = self.compute_nodal_strains(location=loc, smoothing=smoothing)
            for key, arr in result.to_dict().items():
                out[f"{prefix}_{key}"] = arr
        return out

    # ------------------------------------------------------------------
    # Nodal strains
    # ------------------------------------------------------------------
    def compute_element_strains(
        self,
        location: StressLocation = StressLocation.MIDDLE,
        gauss_point: Tuple[float, float] = (0.0, 0.0),
    ) -> StrainResult:
        """Compute a single representative strain per element (centroid).

        Mirrors ``compute_element_stresses`` but returns the engineering
        strain tensor instead of the Cauchy stress.

        Evaluates ε_m + z·κ at *gauss_point* / *location*.

        Parameters
        ----------
        location : StressLocation, default MIDDLE
            Through-thickness position (shells).
        gauss_point : tuple of float, default (0.0, 0.0)
            Parametric coordinates for shell evaluation.

        Returns
        -------
        StrainResult
            Element-centroid strains.  Normal strains are dimensionless;
            shear strains are *engineering* values (γ = 2ε).
        """
        r0, s0 = gauss_point

        # ------------------------------------------------------------------
        # Fast path: delegate to the Rust assembler when available and the
        # caller requests the canonical centroid point (0, 0).
        # ------------------------------------------------------------------
        _rust = getattr(self.domain, "_rust", None)
        if _rust is not None and r0 == 0.0 and s0 == 0.0:
            _Z_FACTOR = {
                StressLocation.BOTTOM: -0.5,
                StressLocation.MIDDLE: 0.0,
                StressLocation.TOP: +0.5,
            }
            z_factor = _Z_FACTOR.get(location, 0.0)
            _, eps_all = _rust.compute_stress_field(self.u, z_factor, 2)
            return self._build_strain_result(eps_all)

        raise NotImplementedError(
            "compute_element_strains: non-centroid gauss_point evaluation requires "
            "r0=0.0, s0=0.0 and a live Rust assembler. "
            "The Python per-element fallback has been removed."
        )

    def compute_nodal_strains(
        self,
        location: StressLocation = StressLocation.MIDDLE,
        smoothing: str = "average",
    ) -> StrainResult:
        """Compute smoothed (continuous) nodal strain fields.

        Applies the same procedure as ``compute_nodal_stresses`` but for
        the strain tensor:

        * **Shell** — ε(z) = ε_m + z·κ evaluated at each element's
          parametric node coordinates, then averaged at shared nodes.

        Parameters
        ----------
        location : StressLocation, default MIDDLE
            Through-thickness position (shells).
        smoothing : str, default ``"average"``
            ``"average"`` or ``"area_weighted"``.

        Returns
        -------
        StrainResult
            Nodal strain field with shape ``(n_nodes,)`` for each
            component.
        """
        n_nodes = self.n_nodes

        _rust = getattr(self.domain, "_rust", None)
        if _rust is not None:
            _Z_FACTOR_NE = {
                StressLocation.BOTTOM: -0.5,
                StressLocation.MIDDLE: 0.0,
                StressLocation.TOP: +0.5,
            }
            z_factor = _Z_FACTOR_NE.get(location, 0.0)
            _, eps_elem = _rust.compute_stress_field(self.u, z_factor, 2)
            # eps_elem: (n_elems, 6) — one strain state per element

            if smoothing == "area_weighted":
                raise NotImplementedError(
                    "compute_nodal_strains: area_weighted smoothing requires "
                    "Python element objects which have been removed. "
                    "Use smoothing='average' instead."
                )

            eps_sum = np.zeros((n_nodes, 6))
            weight_sum = np.zeros(n_nodes)

            for elem_idx, mesh_elem in enumerate(self.domain.elements):
                node_ids = mesh_elem.node_ids
                for node_id in node_ids:
                    gn = self._node_id_to_index[node_id]
                    eps_sum[gn, :] += eps_elem[elem_idx, :]
                    weight_sum[gn] += 1.0

            mask = weight_sum > 0
            eps_avg = np.zeros((n_nodes, 6))
            eps_avg[mask] = eps_sum[mask] / weight_sum[mask, np.newaxis]
            return self._build_strain_result(eps_avg)

        raise NotImplementedError(
            "compute_nodal_strains: requires a live Rust assembler (domain._rust). "
            "The Python per-element fallback has been removed."
        )

    # ------------------------------------------------------------------
    # Derived quantities helpers
    # ------------------------------------------------------------------
    @staticmethod
    def _von_mises_3d(s: np.ndarray) -> np.ndarray:
        """Compute the Von Mises equivalent stress from a Voigt array.

        The Von Mises (or Huber–Mises–Hencky) yield criterion is the
        most widely used scalar measure of the stress state.  It
        represents the distortional strain energy per unit volume and
        is used to predict the onset of yielding in ductile metals.

        The general 3-D formula is:

            σ_vm = √( ½[(σ_xx−σ_yy)² + (σ_yy−σ_zz)² + (σ_zz−σ_xx)²]
                       + 3[τ_xy² + τ_yz² + τ_zx²] )

        When σ_zz = τ_yz = τ_zx = 0 (plane stress) this reduces to
        the familiar shell formula:

            σ_vm = √(σ_xx² + σ_yy² − σ_xx·σ_yy + 3·τ_xy²)

        Parameters
        ----------
        s : np.ndarray, shape (N, 6)
            Voigt stress array with columns [σ_xx, σ_yy, σ_zz, τ_xy,
            τ_yz, τ_zx].

        Returns
        -------
        np.ndarray, shape (N,)
            Von Mises equivalent stress (scalar ≥ 0 at each point).
        """
        sxx, syy, szz = s[:, 0], s[:, 1], s[:, 2]
        txy, tyz, tzx = s[:, 3], s[:, 4], s[:, 5]
        return np.sqrt(
            0.5 * ((sxx - syy) ** 2 + (syy - szz) ** 2 + (szz - sxx) ** 2)
            + 3.0 * (txy**2 + tyz**2 + tzx**2)
        )

    @staticmethod
    def _principal_2d(sxx, syy, sxy):
        """Compute 2-D principal stresses via Mohr’s circle.

        For a plane-stress state the Cauchy stress tensor has the form:

            | σ_xx  τ_xy |
            | τ_xy  σ_yy |

        The principal stresses are the eigenvalues of this 2×2 tensor:

            σ₁,₂ = (σ_xx+σ_yy)/2  ±  R

        where R = √[(σ_xx−σ_yy)²/4 + τ_xy²] is the Mohr’s circle
        radius.  The maximum in-plane shear stress is τ_max = R, and
        the principal angle (angle from the x-axis to the first
        principal direction) is:

            θ_p = ½ arctan(2τ_xy / (σ_xx − σ_yy))

        Parameters
        ----------
        sxx, syy, sxy : np.ndarray
            In-plane stress components (broadcastable).

        Returns
        -------
        sigma_1 : np.ndarray
            Maximum principal stress.
        sigma_2 : np.ndarray
            Minimum principal stress.
        tau_max : np.ndarray
            Maximum in-plane shear stress (Mohr’s circle radius).
        theta_p : np.ndarray
            Principal angle in radians.
        """
        avg = (sxx + syy) / 2
        diff = (sxx - syy) / 2
        R = np.sqrt(diff**2 + sxy**2)
        return avg + R, avg - R, R, 0.5 * np.arctan2(2 * sxy, sxx - syy)

    # ------------------------------------------------------------------
    def _build_stress_result(self, sigma: np.ndarray) -> StressResult:
        """Assemble a ``StressResult`` from a raw Voigt stress array.

        Given an ``(N, 6)`` array of Voigt stresses (nodal or elemental),
        this factory method computes all derived quantities (Von Mises,
        principal stresses, maximum shear, principal angle) and packages
        them into a ``StressResult`` dataclass.  Principal stresses come
        from 2-D Mohr’s circle for σ₁, σ₂, θ_p, so only the in-plane
        components are populated.

        Parameters
        ----------
        sigma : np.ndarray, shape (N, 6)
            Voigt stress array [σ_xx, σ_yy, σ_zz, τ_xy, τ_yz, τ_zx].

        Returns
        -------
        StressResult
        """
        vm = self._von_mises_3d(sigma)
        s1, s2, tau, angle = self._principal_2d(sigma[:, 0], sigma[:, 1], sigma[:, 3])

        return StressResult(
            sigma_xx=sigma[:, 0],
            sigma_yy=sigma[:, 1],
            sigma_xy=sigma[:, 3],
            von_mises=vm,
            sigma_1=s1,
            sigma_2=s2,
            tau_max=tau,
            principal_angle=angle,
        )

    def _build_strain_result(self, eps: np.ndarray) -> StrainResult:
        """Assemble a ``StrainResult`` from a raw Voigt strain array.

        Computes principal strains and maximum shear strain via in-plane
        Mohr’s circle on the (ε_xx, ε_yy, γ_xy/2) components.

        Parameters
        ----------
        eps : np.ndarray, shape (N, 6)
            Voigt strain array [ε_xx, ε_yy, ε_zz, γ_xy, γ_yz, γ_zx].

        Returns
        -------
        StrainResult
        """
        exx, eyy, gxy = eps[:, 0], eps[:, 1], eps[:, 3]

        avg = (exx + eyy) / 2
        diff = (exx - eyy) / 2
        R = np.sqrt(diff**2 + (gxy / 2) ** 2)

        e1 = avg + R
        e2 = avg - R
        gmax = 2 * R

        return StrainResult(
            epsilon_xx=exx,
            epsilon_yy=eyy,
            gamma_xy=gxy,
            epsilon_1=e1,
            epsilon_2=e2,
            gamma_max=gmax,
        )
