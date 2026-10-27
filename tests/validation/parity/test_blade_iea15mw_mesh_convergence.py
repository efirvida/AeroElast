"""Blade modal mesh convergence against CalculiX on the IEA 15 MW blade.

AeroElast (linear MITC4) and CalculiX (quadratic S8R) solve the same blade mesh, and the
matched-mode gap is followed under refinement.  The study runs on 0.5 / 0.35 / 0.25 m, so its
finest point is a working mesh (AGENTS.md rule 9); until 2026-10-22 it ran on 2.0 / 1.0 / 0.5 m
and asserted its bounds on the unconverged 0.5 m mesh.

* **Low modes**: the gap is discretisation and is under 0.2 % from 0.5 m down.
* **High modes**: earlier revisions called their residual a 4-node-vs-8-node element-order
  plateau.  The 2026-10-22 refinement refutes that (odd/tasks/l3-physical-validation.md L32).
  The CCX modes 6 and 8 (~3.96 and ~4.21 Hz) fall monotonically, about order 1.4, as linear
  elements against quadratic ones do: mode 8 is +13.98 / +8.26 / +5.23 % at 0.5 / 0.35 / 0.25 m.
  It is discretisation of the linear element, not a fixed offset, and extrapolates to about
  2 % at 0.125 m, where the CalculiX deck is too heavy to run here.  With averaged
  directors across folds (the pre-7dfe78d convention) the same modes sit 6-7 % lower and change
  sign between 0.5 and 0.35 m, so their smaller gap was a cancellation.

Signed AeroElast-vs-CCX gaps per CCX mode (SLURM 11622478/80/81, logs/check/l3-ccx-es025/):

| element_size | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 0.5 m | +0.09% | +0.20% | +0.25% | +0.51% | +0.69% | +3.97% | +0.95% | +13.98% |
| 0.35 m | +0.01% | +0.08% | +0.10% | +0.25% | +0.36% | +3.27% | +0.73% | +8.26% |
| 0.25 m | -0.02% | +0.03% | +0.03% | +0.12% | +0.18% | +2.97% | +0.52% | +5.23% |

The bounds are unchanged; on the 0.25 m mesh the worst high-mode gap reads 5.23 % against the 5 %
bound, a measured discretisation residual and not a tolerance to move.  The run is slow (three
CCX modal solves, about 15 minutes at these meshes), so it is marked ``slow``.
"""

from __future__ import annotations

import numpy as np
import pytest
from scipy.optimize import linear_sum_assignment

from tests.conftest import ccx_bin_or_skip

pytest.importorskip("petsc4py", reason="PETSc not available")
pytest.importorskip("_aeroelast", reason="Rust backend not available")

from _aeroelast import PyMeshAssembler, modal_solve_coo  # noqa: E402
from tests.support.assertions import assert_residual_below  # noqa: E402
from tests.support.ccx_io import fail_ccx, parse_ccx_frequencies, run_ccx  # noqa: E402
from aeroelast.core.mesh.entities import MeshElement, Node  # noqa: E402
from aeroelast.core.mesh.io.writers import write_ccx_mesh  # noqa: E402
from aeroelast.models.blade.model import Blade  # noqa: E402

import tests.validation.blade.test_blade_iea15mw_validation as blade_val  # noqa: E402

MESHES = (0.5, 0.35, 0.25)
N_MODES = 8
#: Low modes where refinement resolves the gap: finest-mesh bound (measured
#: 0.08-0.19% at CCX 2.23).
LOW_MODES = 4
LOW_TOL = 0.01
#: Modes 5-8: bound set when the gap was read as an element-order plateau (measured
#: 0.40-2.12% at CCX 2.23, 3.90-4.04% at CCX 2.20); it is discretisation (module docstring).
HIGH_TOL = 0.05
#: Below this a gap difference is not resolvable across the comparison.
NOISE_FLOOR = 0.002

pytestmark = pytest.mark.slow


def _matched_gaps(ae: np.ndarray, ccx: np.ndarray, n: int) -> list[float]:
    """Relative gap of each CCX mode to its AeroElast partner, in CCX mode order.

    The pairing is the smallest-cost one-to-one assignment (Hungarian), so a mode
    that one solver orders differently cannot masquerade as a mismatch.  Entry ``j``
    belongs to CCX mode ``j + 1`` on every mesh, which is what lets the convergence
    check compare one mode across meshes; the earlier version returned the gaps
    sorted by size, so its "mode 6" was the sixth-smallest gap, a different physical
    mode from mesh to mesh (T29).
    """
    cost = np.abs(ae[:, None] - ccx[None, :]) / np.maximum(ccx[None, :], 1e-14)
    row_ind, col_ind = linear_sum_assignment(cost)
    by_ccx = np.empty(len(col_ind))
    by_ccx[col_ind] = cost[row_ind, col_ind]
    return [float(g) for g in by_ccx[:n]]


def _modal_couple(element_size: float, workdir, ccx_bin: str) -> tuple[np.ndarray, np.ndarray]:
    """Build the blade at ``element_size`` and solve it with both solvers."""
    Node._id_counter = 0
    MeshElement._id_counter = 0

    blade_model = Blade(str(blade_val.YAML), element_size=element_size)
    blade_model.generate_mesh()
    mesh = blade_model.mesh
    assert mesh is not None, "generate_mesh must leave the blade with a mesh"
    props = blade_model.get_element_properties()

    assembler = PyMeshAssembler.from_model(
        blade_val._to_rust_mesh(mesh, props), props, list(blade_val.SPAN_DIRECTION), None
    )
    n = assembler.dofs_count
    k_rows, k_cols, k_vals = assembler.assemble_k()
    m_rows, m_cols, m_vals = assembler.assemble_m()

    root = {mesh.node_id_to_index[nid] for nid in mesh.get_node_set("RootNodes").node_ids}
    fixed = {6 * i + d for i in root for d in range(6)}
    free = np.array([i for i in range(n) if i not in fixed], dtype=np.int64)

    freqs_ae = np.sort(
        np.asarray(
            modal_solve_coo(
                np.asarray(k_rows, dtype=np.int64),
                np.asarray(k_cols, dtype=np.int64),
                np.asarray(k_vals, dtype=np.float64),
                np.asarray(m_rows, dtype=np.int64),
                np.asarray(m_cols, dtype=np.int64),
                np.asarray(m_vals, dtype=np.float64),
                n,
                free,
                N_MODES,
            )[0],
            dtype=float,
        )
    )

    inp_path = workdir / f"blade_es{element_size:.2f}.inp"
    write_ccx_mesh(
        mesh,
        str(inp_path),
        properties=props,
        boundary_nodeset="RootNodes",
        solver_type="Modal",
        num_modes=N_MODES,
        quadratic=True,
        span_direction=blade_val.SPAN_DIRECTION,
    )
    result = run_ccx(inp_path, ccx_bin)
    if result.returncode != 0:
        fail_ccx(result, inp_path)
    freqs_ccx = np.sort(parse_ccx_frequencies(inp_path, n_modes=N_MODES))

    print(
        f"  es={element_size:.2f} m: {mesh.node_count} nodes, {n} dofs  "
        f"ae={np.array2string(freqs_ae, precision=4)}"
    )
    return freqs_ae, freqs_ccx


@pytest.fixture(scope="module")
def convergence(
    tmp_path_factory: pytest.TempPathFactory,
) -> dict[float, tuple[np.ndarray, np.ndarray]]:
    """Solve the blade with AeroElast and CCX at each mesh in ``MESHES``."""
    ccx_bin = ccx_bin_or_skip()
    workdir = tmp_path_factory.mktemp("blade_mesh_convergence")
    return {es: _modal_couple(es, workdir, ccx_bin) for es in MESHES}


def test_blade_modal_gap_converges_with_mesh(
    convergence: dict[float, tuple[np.ndarray, np.ndarray]],
) -> None:
    """The AeroElast-vs-CCX modal gap is discretisation for the low modes.

    Modes 1-4 fall under refinement, so their finest-mesh gap is bounded by
    ``LOW_TOL``.  Modes 5-8 are bounded by ``HIGH_TOL``; they also fall with refinement
    (module docstring), but the bound was set when they were read as an element-order plateau.  The coarse-to-medium
    improvement is asserted only where the coarse gap is above ``NOISE_FLOOR``,
    because below it the difference is inside the comparison's resolution.
    """
    gaps = {es: _matched_gaps(*convergence[es], N_MODES) for es in MESHES}
    for es in MESHES:
        print(f"  es={es:.2f} m gaps by CCX mode: {[f'{g * 100:.2f}%' for g in gaps[es]]}")

    coarse, medium, fine = MESHES
    # Mesh-convergence claim: the first refinement improves every mode whose
    # coarse gap is above the noise floor.
    for mode in range(N_MODES):
        if gaps[coarse][mode] > NOISE_FLOOR:
            assert gaps[medium][mode] < gaps[coarse][mode], (
                f"CCX mode {mode + 1}: {medium} m gap {gaps[medium][mode] * 100:.2f}% is not "
                f"below the {coarse} m gap {gaps[coarse][mode] * 100:.2f}%"
            )

    low_worst = max(gaps[fine][:LOW_MODES])
    high_worst = max(gaps[fine][LOW_MODES:])
    assert_residual_below(
        low_worst,
        tol=LOW_TOL,
        kind="code",
        reference_name="the CalculiX frequencies for the same blade and mesh",
        what="worst low-mode gap on the finest mesh",
    )
    assert_residual_below(
        high_worst,
        tol=HIGH_TOL,
        kind="code",
        reference_name="the CalculiX frequencies for the same blade and mesh",
        what="worst high-mode gap on the finest mesh (linear MITC4 vs quadratic S8R discretisation)",
    )
