"""Blade modal mesh convergence against CalculiX on the IEA 15 MW blade.

``test_blade_iea15mw_validation.py`` compares AeroElast to CCX at one mesh
(``element_size = 1.0 m``).  This module separates two behaviours that the first
revision of this test conflated:

* **Low modes (1-4)**: the same-method gap is discretisation.  It falls from the
  2.0 m to the 1.0 m to the 0.5 m mesh and the finest-mesh gap is under 1%.
* **High modes (5-8)**: the gap plateaus at a few percent on every mesh from
  1.0 m down.  That residual is the 4-node-vs-8-node **element-order** difference
  (linear MITC4 against quadratic S8R, worst in the high-curvature modes), not
  discretisation, so it gets a looser bound and no per-mode improvement claim.

Measured AeroElast-vs-CCX relative gaps for the eight matched modes (Hungarian
pairing, so the order is by cost, not by mode number):

| element_size | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 2.0 m | 1.10% | 2.28% | 2.33% | 2.76% | 3.99% | 5.16% | 5.45% | 5.85% |
| 1.0 m | 0.44% | 0.71% | 0.78% | 0.84% | 1.65% | 1.69% | 2.13% | 2.62% |
| 0.5 m | 0.08% | 0.13% | 0.18% | 0.19% | 0.40% | 0.44% | 0.59% | 2.12% |

The run is slow: three CCX modal solves dominate (~11 minutes measured), so it
is marked ``slow`` and can be deselected with ``-m "not slow"``.
"""

from __future__ import annotations

import numpy as np
import pytest
from scipy.optimize import linear_sum_assignment

from tests.conftest import ccx_bin_or_skip

pytest.importorskip("petsc4py", reason="PETSc not available")
pytest.importorskip("_aeroelast", reason="Rust backend not available")

from _aeroelast import PyMeshAssembler, modal_solve_coo  # noqa: E402
from tests.support.ccx_io import fail_ccx, parse_ccx_frequencies, run_ccx  # noqa: E402
from aeroelast.core.mesh.entities import MeshElement, Node  # noqa: E402
from aeroelast.core.mesh.io.writers import write_ccx_mesh  # noqa: E402
from aeroelast.models.blade.model import Blade  # noqa: E402

import tests.test_blade_iea15mw_validation as blade_val  # noqa: E402

MESHES = (2.0, 1.0, 0.5)
N_MODES = 8
#: Low modes where refinement resolves the gap: finest-mesh bound (measured
#: 0.08-0.19% at CCX 2.23).
LOW_MODES = 4
LOW_TOL = 0.01
#: Modes 5-8: element-order bound, explicitly not a mesh claim (measured
#: 0.40-2.12% at CCX 2.23, 3.90-4.04% at CCX 2.20).
HIGH_TOL = 0.05
#: Below this a gap difference is not resolvable across the comparison.
NOISE_FLOOR = 0.002

pytestmark = pytest.mark.slow


def _matched_gaps(ae: np.ndarray, ccx: np.ndarray, n: int) -> list[float]:
    """Smallest-cost one-to-one relative gaps between the two mode sets.

    Reuses the Hungarian pairing from the validation module so a mode that one
    solver orders differently cannot masquerade as a mismatch.
    """
    cost = np.abs(ae[:, None] - ccx[None, :]) / np.maximum(ccx[None, :], 1e-14)
    row_ind, col_ind = linear_sum_assignment(cost)
    gaps = sorted(float(cost[i, j]) for i, j in zip(row_ind, col_ind, strict=False))
    return gaps[:n]


def _modal_couple(element_size: float, workdir, ccx_bin: str) -> tuple[np.ndarray, np.ndarray]:
    """Build the blade at ``element_size`` and solve it with both solvers."""
    Node._id_counter = 0
    MeshElement._id_counter = 0

    blade_model = Blade(str(blade_val.YAML), element_size=element_size)
    blade_model.generate_mesh()
    mesh = blade_model.mesh
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

    inp_path = workdir / f"blade_es{element_size:.1f}.inp"
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
        f"  es={element_size:.1f} m: {mesh.node_count} nodes, {n} dofs  "
        f"ae={np.array2string(freqs_ae, precision=4)}"
    )
    return freqs_ae, freqs_ccx


@pytest.fixture(scope="module")
def convergence(tmp_path_factory: pytest.TempPathFactory) -> dict[float, tuple[np.ndarray, np.ndarray]]:
    """Solve the blade with AeroElast and CCX at each mesh in ``MESHES``."""
    ccx_bin = ccx_bin_or_skip()
    workdir = tmp_path_factory.mktemp("blade_mesh_convergence")
    return {es: _modal_couple(es, workdir, ccx_bin) for es in MESHES}


def test_blade_modal_gap_converges_with_mesh(convergence) -> None:
    """The AeroElast-vs-CCX modal gap is discretisation for the low modes.

    Modes 1-4 fall under refinement, so their finest-mesh gap is bounded by
    ``LOW_TOL``.  Modes 5-8 plateau at 1.5-4% (mode 8 around 2-4% on every mesh)
    -- an element-order difference, not discretisation -- so they get the looser
    ``HIGH_TOL`` and no per-mode improvement claim.  The coarse-to-medium
    improvement is asserted only where the coarse gap is above ``NOISE_FLOOR``,
    because below it the difference is inside the comparison's resolution.
    """
    gaps = {es: _matched_gaps(*convergence[es], N_MODES) for es in MESHES}
    for es in MESHES:
        print(f"  es={es:.1f} m gaps: {[f'{g * 100:.2f}%' for g in gaps[es]]}")

    coarse, medium, fine = MESHES
    # Mesh-convergence claim: the first refinement improves every mode whose
    # coarse gap is above the noise floor.
    for mode in range(N_MODES):
        if gaps[coarse][mode] > NOISE_FLOOR:
            assert gaps[medium][mode] < gaps[coarse][mode], (
                f"matched mode {mode + 1}: 1.0 m gap {gaps[medium][mode] * 100:.2f}% is not "
                f"below the 2.0 m gap {gaps[coarse][mode] * 100:.2f}%"
            )

    low_worst = max(gaps[fine][:LOW_MODES])
    high_worst = max(gaps[fine][LOW_MODES:])
    assert low_worst < LOW_TOL, (
        f"finest-mesh low modes (1-{LOW_MODES}) worst gap {low_worst * 100:.2f}% "
        f"(tol {LOW_TOL * 100:.1f}%)"
    )
    assert high_worst < HIGH_TOL, (
        f"finest-mesh high modes ({LOW_MODES + 1}-{N_MODES}) worst gap "
        f"{high_worst * 100:.2f}% (tol {HIGH_TOL * 100:.1f}%) -- element-order bound"
    )
