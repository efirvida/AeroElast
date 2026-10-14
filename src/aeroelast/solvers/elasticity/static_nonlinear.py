"""Nonlinear static solver — PETSc SNES Newton-Raphson via Rust assembler."""

import logging
from time import perf_counter
from typing import List, Optional

import numpy as np

from aeroelast.core.mesh import MeshModel
from aeroelast.solvers.solver import Solver

_log = logging.getLogger(__name__)


class StaticNonlinearSolver(Solver):
    """
    Nonlinear static solver — total-Lagrangian Newton (default) or legacy updated-Lagrangian.

    The default ``scheme: newton`` solves ``f_int(u) - lam F = 0`` on the Green-Lagrange internal
    force of the Rust assembler with its consistent tangent (``assemble_kt``), load-stepped with
    step halving on failure.  It is objective for rotations about any axis (rigid rotations give
    ``|f_int| ~ 1e-8`` against a 2e7 membrane stiffness up to 2 rad; ``kt`` matches a finite
    difference of ``f_int`` to 1e-10 at multi-axis states) and agrees with CalculiX S8R NLGEOM to
    0.2-0.8 % on a pre-twisted beam under combined bending and torsion.

    Validated envelope: nodal rotations up to about ``pi``.  The rotation DOFs are a total
    rotation vector with a quadratic director, so beyond ``pi`` the corrector stops converging
    and the solve raises; it never returns a wrong answer silently.

    ``scheme: ul`` keeps the previous updated-Lagrangian loop (``K_e(x_ref) du = dlam F``,
    re-referencing after each step).  It has no equilibrium iteration and no initial-stress
    stiffness, so refining the steps does **not** converge to the dead-load equilibrium; see
    *Schemes*.  Keep it only to reproduce old results.

    YAML solver parameters (all optional)
    --------------------------------------
    scheme : str
        ``"newton"`` (default) or ``"ul"``; see *Schemes* below.
    continuation_steps : int
        Number of initial load increments (default: 8; use 16-32 for blades).
    continuation_max_steps : int
        Maximum number of step attempts allowed, halved steps included (default: 64).
    newton_rtol, newton_etol, newton_max_it
        ``newton`` only: residual tolerance relative to the applied load (default 1e-5), energy
        tolerance (default 1e-10), iterations per load step (default 25).
    diagnostics : bool
        Emit per-step diagnostics (default: False).

    ``atol``, ``rtol``, ``stol``, ``max_it`` and ``continuation`` belonged to the historical SNES
    path and are not read any more.

    Schemes
    -------
    ``"newton"`` (default) iterates ``f_int(u) - lam F = 0`` with the consistent tangent.  On
    the b = 0.1 shell strip it agrees with CalculiX S8R NLGEOM to 0.1-0.3 % (alpha = 4/7, 1, 2),
    on a narrow strip it reproduces the exact elastica to <0.4 %, and on the pre-twisted beam
    (Simo geometry, tip load 10) it is 0.2-0.5 % from CalculiX where ``"ul"`` is 23 % off.
    The tangent-consistency doubts recorded here on 2026-09-18 no longer hold: ``f_int`` is
    objective and ``kt`` consistent with it for multi-axis states.

    ``"ul"`` is a tangent flow: it never checks equilibrium and omits the initial-stress
    stiffness.  Against CalculiX on the b = 0.1 strip it reads +2.0 % too soft at 16 deg of tip
    rotation (alpha = 4/7), +4.4 % at 26 deg and +9.6 % at 45 deg; on the pre-twisted beam
    +12 % (load 5) and +23 % (load 10), independent of the mesh.  It is exact only for a pure
    end moment (no shear or axial force to rotate), and then only up to about 1 rad.

    Notes
    -----
    The assembly and the linear solve run entirely in Rust/PETSc.  The historical SNES path
    (``nonlinear_static_solve_coo``) is kept in the Rust crate but is not used here.
    """

    _DEFAULT_ATOL: float = 1e-10
    _DEFAULT_RTOL: float = 1e-8
    _DEFAULT_STOL: float = 1e-8
    _DEFAULT_MAX_IT: int = 100
    # 100 instead of 50: the large-displacement cantilever benchmark needs one
    # continuation substep that exceeds 50 Newton iterations, so the old default
    # diverged on a case the suite exercises.  Measured: with 50 the solve fails
    # (reason=-5) after burning the continuation budget; with 100 it converges
    # to the elastica case in 67 total iterations.  Larger loads (>~1.5x) still
    # need finer continuation and remain outside the solver's validated envelope.
    _DEFAULT_CONTINUATION: bool = True
    _DEFAULT_CONTINUATION_STEPS: int = 8
    _DEFAULT_CONTINUATION_MAX_STEPS: int = 64
    _SCHEMES: tuple = ("ul", "newton")
    _DEFAULT_SCHEME: str = "newton"
    # Newton scheme.  The residual norm cannot reach the SNES-style 1e-8: the shell stiffness
    # spans ~12 orders of magnitude (membrane over bending), so the direct solve leaves a
    # force residual of ~1e-6 relative.  The energy test |du.R| <= etol * |du_1.R_0| is
    # scale-free and reaches machine level, so it is the primary criterion.
    _DEFAULT_NEWTON_RTOL: float = 1e-5
    _DEFAULT_NEWTON_ETOL: float = 1e-10
    _DEFAULT_NEWTON_MAX_IT: int = 25
    _DEFAULT_DIAGNOSTICS: bool = False
    _DEFAULT_DIAGNOSTICS_EVERY: int = 1

    def __init__(self, mesh: MeshModel, fem_model_properties: dict):
        super().__init__(mesh, fem_model_properties)
        self.u: Optional[np.ndarray] = None
        self._iterations: int = 0
        self._residual_norm: float = 0.0
        self._converged_reason: int = 0

    # ------------------------------------------------------------------
    # Public interface
    # ------------------------------------------------------------------

    def solve(self) -> np.ndarray:
        """
        Run the nonlinear static solve with the configured ``scheme``.

        Returns
        -------
        np.ndarray
            Total displacement vector (n_dofs,) measured from the original
            reference configuration.

        Raises
        ------
        RuntimeError
            If the Rust assembler is not available or the load limit is not
            reached within ``continuation_max_steps``.

        Notes
        -----
        With ``scheme: ul`` each step solves the *linearized* problem on the current reference
        geometry (``K_e(x_ref) du = dlam F``, ``x_ref <- x_ref + du``); see the class docstring
        for why that does not converge to the dead-load equilibrium.
        """
        import _aeroelast  # noqa: PLC0415 — optional Rust extension

        if self.domain._rust is None:
            raise RuntimeError(
                "StaticNonlinearSolver requires the Rust assembler. "
                "Make sure _aeroelast is built and the mesh is assembled."
            )

        f_ext = self._build_f_ext()
        dirichlet_dofs = self._collect_dirichlet_dofs()
        params = self.solver_params if isinstance(self.solver_params, dict) else {}

        n_dof = self.domain.dofs_count
        free_dofs = np.array(
            sorted(set(range(n_dof)) - set(dirichlet_dofs)), dtype=np.int64
        )

        scheme = str(params.get("scheme", self._DEFAULT_SCHEME)).lower()
        if scheme not in self._SCHEMES:
            raise ValueError(
                f"Unknown StaticNonlinearSolver scheme {scheme!r}; expected one of {self._SCHEMES}"
            )

        n_steps = max(1, int(params.get("continuation_steps", self._DEFAULT_CONTINUATION_STEPS)))
        max_steps = max(
            n_steps,
            int(params.get("continuation_max_steps", self._DEFAULT_CONTINUATION_MAX_STEPS)),
        )
        diagnostics = bool(params.get("diagnostics", self._DEFAULT_DIAGNOSTICS))

        if scheme == "newton":
            return self._solve_newton(f_ext, free_dofs, n_dof, params, n_steps, max_steps, diagnostics)

        u_tot = np.zeros(n_dof, dtype=np.float64)
        lam = 0.0
        dlam = 1.0 / float(n_steps)
        step = 0

        t0 = perf_counter()
        while lam < 1.0 - 1e-14:
            if step >= max_steps:
                _log.error(
                    "[UL] max steps reached (step=%d, max=%d) at lambda=%.6f",
                    step,
                    max_steps,
                    lam,
                )
                break
            step += 1
            lam_try = min(1.0, lam + dlam)
            f_inc = f_ext * (lam_try - lam)

            # Linearized step on the current reference geometry:
            #   K_e(x_ref) · du = dλ·F_ext
            k_rows, k_cols, k_vals = self.domain._rust.assemble_k()
            du = np.asarray(
                _aeroelast.linear_static_solve_coo(
                    k_rows.astype(np.int64),
                    k_cols.astype(np.int64),
                    k_vals.astype(np.float64),
                    f_inc,
                    n_dof,
                    free_dofs,
                ),
                dtype=np.float64,
            )

            u_tot += du
            # Advance the reference configuration (rebuilds per-element
            # geometry from x_ref + du).
            self.domain._rust.update_reference(np.ascontiguousarray(du))
            lam = lam_try

            if diagnostics:
                _log.info(
                    "[UL step %3d] lambda=%.5f |du|=%.3e |u|=%.3e (%.2fs)",
                    step,
                    lam,
                    float(np.linalg.norm(du)),
                    float(np.linalg.norm(u_tot)),
                    perf_counter() - t0,
                )

        converged = lam >= 1.0 - 1e-14
        self.u = u_tot
        self._iterations = step
        self._residual_norm = 0.0
        self._converged_reason = 2 if converged else -1

        if not converged:
            raise RuntimeError(
                "StaticNonlinearSolver (UL) did not reach full load: "
                f"lambda={lam:.6f} after {step} steps (max_steps={max_steps})"
            )
        return self.u

    def _solve_newton(
        self,
        f_ext: np.ndarray,
        free_dofs: np.ndarray,
        n_dof: int,
        params: dict,
        n_steps: int,
        max_steps: int,
        diagnostics: bool,
    ) -> np.ndarray:
        """Total-Lagrangian Newton on ``f_int(u) - lam F`` with load stepping.

        A step whose corrector fails is retried with half the load increment; the attempts
        count against ``max_steps``.
        """
        import _aeroelast  # noqa: PLC0415 — optional Rust extension

        rust = self.domain._rust
        rtol = float(params.get("newton_rtol", self._DEFAULT_NEWTON_RTOL))
        etol = float(params.get("newton_etol", self._DEFAULT_NEWTON_ETOL))
        max_it = int(params.get("newton_max_it", self._DEFAULT_NEWTON_MAX_IT))
        f_free_norm = float(np.linalg.norm(f_ext[free_dofs]))

        def corrector(u0: np.ndarray, lam_target: float) -> tuple[np.ndarray, bool, int, float]:
            u = u0.copy()
            f_target = lam_target * f_ext
            ref_energy = None
            res = float("inf")
            for it in range(max_it):
                r = np.asarray(rust.assemble_fint(np.ascontiguousarray(u), True)) - f_target
                res = float(np.linalg.norm(r[free_dofs]))
                if res <= rtol * lam_target * f_free_norm:
                    return u, True, it, res
                rows, cols, vals = rust.assemble_kt(np.ascontiguousarray(u))
                du = np.asarray(
                    _aeroelast.linear_static_solve_coo(
                        rows.astype(np.int64),
                        cols.astype(np.int64),
                        vals.astype(np.float64),
                        -r,
                        n_dof,
                        free_dofs,
                    ),
                    dtype=np.float64,
                )
                if not np.all(np.isfinite(du)):
                    return u0, False, it, res
                energy = abs(float(du[free_dofs] @ r[free_dofs]))
                if ref_energy is None:
                    ref_energy = energy
                u += du
                if ref_energy > 0.0 and energy <= etol * ref_energy:
                    return u, True, it + 1, res
            return u0, False, max_it, res

        u_tot = np.zeros(n_dof, dtype=np.float64)
        lam = 0.0
        dlam = 1.0 / float(n_steps)
        attempts = 0
        total_it = 0
        res = 0.0
        t0 = perf_counter()
        while lam < 1.0 - 1e-14:
            if attempts >= max_steps:
                _log.error(
                    "[Newton] max steps reached (attempts=%d, max=%d) at lambda=%.6f",
                    attempts,
                    max_steps,
                    lam,
                )
                break
            attempts += 1
            lam_try = min(1.0, lam + dlam)
            u_new, ok, n_it, res = corrector(u_tot, lam_try)
            total_it += n_it
            if not ok:
                dlam *= 0.5
                _log.warning(
                    "[Newton] corrector failed at lambda=%.6f (|R|=%.3e); halving dlam to %.3e",
                    lam_try,
                    res,
                    dlam,
                )
                continue
            u_tot = u_new
            lam = lam_try
            if diagnostics:
                _log.info(
                    "[Newton] lambda=%.5f its=%d |R|=%.3e (%.2fs)",
                    lam,
                    n_it,
                    res,
                    perf_counter() - t0,
                )

        converged = lam >= 1.0 - 1e-14
        self.u = u_tot
        self._iterations = total_it
        self._residual_norm = res
        self._converged_reason = 2 if converged else -1
        if not converged:
            raise RuntimeError(
                "StaticNonlinearSolver (newton) did not reach full load: "
                f"lambda={lam:.6f} after {attempts} attempts (max_steps={max_steps}). "
                "The rotation DOFs are a total rotation vector, valid up to about pi of nodal "
                "rotation; beyond that the corrector cannot converge. Lower the load, or use "
                "scheme='ul' (biased, see the class docstring)."
            )
        return self.u

    def print_solver_info(self) -> None:
        """Print nonlinear solver statistics to stdout."""
        print("\n--- Nonlinear Solver (SNES / Newton-Raphson) ---")
        print(f"  Converged reason : {self._converged_reason}")
        print(f"  Iterations       : {self._iterations}")
        print(f"  Final |R|        : {self._residual_norm:.3e}")

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _build_f_ext(self) -> np.ndarray:
        """Assemble the external force vector from body + nodal loads."""
        n = self.domain.dofs_count
        f = np.zeros(n, dtype=np.float64)

        # Distributed/body loads
        for force in self.body_forces:
            fe = self.domain.assemble_load_vector(force)
            f += np.asarray(fe.getArray(), dtype=np.float64)

        # Concentrated nodal loads
        for load in self.nodal_loads:
            dofs = np.asarray(load.dofs, dtype=np.int64)
            vals = np.asarray(load.force, dtype=np.float64)
            if len(dofs) != len(vals):
                raise ValueError(
                    "Nodal load dof/value length mismatch in nonlinear solver: "
                    f"{len(dofs)} != {len(vals)}"
                )
            f[dofs] += vals

        return f

    def _collect_dirichlet_dofs(self) -> np.ndarray:
        """Collect all constrained DOF indices from Dirichlet BCs."""
        dofs: List[int] = []
        for bc in self.dirichlet_conditions:
            dofs.extend(bc.dofs)
        return np.array(dofs, dtype=np.int64)
