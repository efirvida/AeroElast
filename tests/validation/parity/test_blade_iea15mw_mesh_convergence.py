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

Mode pairing (T29, T12, 2026-10-27).  Until T12 the modes were paired by frequency (Hungarian on
``|f_ae - f_ccx| / f_ccx``).  CCX modes 7 and 8 are within 1 % of each other and swap order
between 0.5 and 0.35 m: at 0.5 m CCX 7 (4.179 Hz) is the second edgewise mode, at 0.35 / 0.25 m
the edgewise mode is CCX 8 (4.183 / 4.179 Hz) and CCX 7 is a y-translation mode with the rotation
about the span axis.  The frequency pairing kept CCX 7 with the edgewise AeroElast mode, with MAC
0.32 / 0.15 at 0.35 / 0.25 m, so T29's "mode 7: 0.89 / 1.01 / 1.01 %" paired different shapes
(SLURM 11625834, ``logs/check/l3/g5/diag-11625834.log``).  By shape: the edgewise mode is
+0.89 / +0.39 / +0.25 % (converges) and the y / span-rotation mode +8.64 / +4.61 % at 0.35 / 0.25 m
with MAC 0.63 / 0.83; at 0.5 m its partner is not identifiable (best MAC 0.625).  The test now
pairs AeroElast and CCX modes by MAC of the corner-node translations and follows one physical
mode across meshes by a span-station signature (``tests/support/modal_pairing.py``).  The gaps are
indexed by the finest-mesh CCX mode order; the pairing MAC is printed per mode.  A finest-mesh
mode with pairing MAC below ``MAC_IDENTIFIABLE`` is reported as not identifiable.

CCX version: the CalculiX resolved on this cluster (``ccx_bin_or_skip``) is **2.20**, not 2.23.
The 2.23 numbers quoted in the constants below came from another binary; the frequency-paired
table that used to stand here (SLURM 11622478/80/81) is superseded by the shape-paired values
above and by the test's own printed lines.

The bounds are unchanged; they are measured discretisation residuals and not tolerances to move.
The run is slow (three CCX modal solves, about 15 minutes at these meshes), so it is marked
``slow``.
"""

from __future__ import annotations

import numpy as np
import pytest

from tests.conftest import ccx_bin_or_skip

pytest.importorskip("petsc4py", reason="PETSc not available")
pytest.importorskip("_aeroelast", reason="Rust backend not available")

from _aeroelast import PyMeshAssembler, modal_solve_coo  # noqa: E402
from tests.support.assertions import assert_residual_below  # noqa: E402
from tests.support.ccx_io import fail_ccx, parse_ccx_frequencies, run_ccx  # noqa: E402
from tests.support.modal_pairing import (  # noqa: E402
    ModalSet,
    check_corner_nodes_first,
    followed_gaps,
    parse_frd_mode_translations,
)
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
#: A pair of modes is taken as the same shape when its MAC is at least this: the usual lower
#: limit of the standard modal-correlation consistency criterion (MAC >= 0.9 consistent,
#: <= 0.05 distinct; Allemang, "The Modal Assurance Criterion - twenty years of use and abuse",
#: Sound and Vibration 37(8), 2003).
MAC_IDENTIFIABLE = 0.8

pytestmark = pytest.mark.slow


def _modal_couple(element_size: float, workdir, ccx_bin: str) -> ModalSet:
    """Build the blade at ``element_size`` and solve it with both solvers.

    Besides the frequencies it returns the translations of the corner nodes of both solvers
    (the eigenvectors are on the free DOFs in sorted order), which the pairing needs.
    """
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

    freqs, flat = modal_solve_coo(
        np.asarray(k_rows, dtype=np.int64),
        np.asarray(k_cols, dtype=np.int64),
        np.asarray(k_vals, dtype=np.float64),
        np.asarray(m_rows, dtype=np.int64),
        np.asarray(m_cols, dtype=np.int64),
        np.asarray(m_vals, dtype=np.float64),
        n,
        free,
        N_MODES,
    )
    freqs_ae = np.asarray(freqs, dtype=float)
    vecs = np.asarray(flat, dtype=float).reshape(len(freqs_ae), len(free))
    full = np.zeros((len(freqs_ae), n))
    full[:, np.sort(free)] = vecs
    coords = np.array([[nd.x, nd.y, nd.z] for nd in mesh.nodes])
    n_nodes = len(coords)
    order_ae = np.argsort(freqs_ae)
    freqs_ae = freqs_ae[order_ae]
    ae_t = full[order_ae].reshape(len(freqs_ae), n_nodes, 6)[:, :, :3]

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
    # The corner nodes keep their index in the quadratic deck: the FRD and the AeroElast
    # translations then share one node numbering.
    check_corner_nodes_first(inp_path.with_suffix(".msh"), coords)
    freqs_ccx = np.asarray(parse_ccx_frequencies(inp_path, n_modes=N_MODES), dtype=float)
    shapes_ccx = parse_frd_mode_translations(inp_path.with_suffix(".frd"), n_nodes, N_MODES)
    order_ccx = np.argsort(freqs_ccx)
    freqs_ccx, shapes_ccx = freqs_ccx[order_ccx], shapes_ccx[order_ccx]

    print(
        f"  es={element_size:.2f} m: {mesh.node_count} nodes, {n} dofs  "
        f"ae={np.array2string(freqs_ae, precision=4)}"
    )
    span = int(np.argmax(np.abs(blade_val.SPAN_DIRECTION)))
    return ModalSet(freqs_ae, freqs_ccx, ae_t, shapes_ccx, coords[:, span])


@pytest.fixture(scope="module")
def convergence(
    tmp_path_factory: pytest.TempPathFactory,
) -> dict[float, ModalSet]:
    """Solve the blade with AeroElast and CCX at each mesh in ``MESHES``."""
    ccx_bin = ccx_bin_or_skip()
    workdir = tmp_path_factory.mktemp("blade_mesh_convergence")
    return {es: _modal_couple(es, workdir, ccx_bin) for es in MESHES}


def test_blade_modal_gap_converges_with_mesh(convergence: dict[float, ModalSet]) -> None:
    """The AeroElast-vs-CCX modal gap is discretisation for the low modes.

    Modes are paired by shape (MAC) and followed across meshes by signature, so entry ``j`` of a
    mesh is the finest-mesh CCX mode ``j + 1`` (``tests/support/modal_pairing.py``).  Modes 1-4
    fall under refinement, so their finest-mesh gap is bounded by ``LOW_TOL``.  Modes 5-8 are
    bounded by ``HIGH_TOL``; they also fall with refinement (module docstring), but the bound was
    set when they were read as an element-order plateau.  The coarse-to-medium improvement is
    asserted only where the coarse gap is above ``NOISE_FLOOR``, because below it the difference
    is inside the comparison's resolution.  A finest-mesh mode whose pairing MAC is below
    ``MAC_IDENTIFIABLE`` has no partner to compare with: it is printed and listed in the failure
    message, and its gap is not taken as a residual.
    """
    fine_key = MESHES[-1]
    followed = followed_gaps(convergence, fine_key)
    gaps = {es: [float(g) for g in followed[es][0]] for es in MESHES}
    macs = {es: [float(m) for m in followed[es][1]] for es in MESHES}
    for es in MESHES:
        print(f"  es={es:.2f} m gaps by finest CCX mode: {[f'{g * 100:.2f}%' for g in gaps[es]]}")
        print(f"  es={es:.2f} m pairing MAC by finest CCX mode: {[f'{m:.3f}' for m in macs[es]]}")

    coarse, medium, fine = MESHES
    not_identifiable = [mode for mode in range(N_MODES) if macs[fine][mode] < MAC_IDENTIFIABLE]
    identifiable = [mode for mode in range(N_MODES) if mode not in not_identifiable]
    for mode in not_identifiable:
        print(
            f"  finest CCX mode {mode + 1}: not identifiable (pairing MAC "
            f"{macs[fine][mode]:.3f} < {MAC_IDENTIFIABLE}), gap {gaps[fine][mode] * 100:.2f}% not used"
        )
    # Mesh-convergence claim: the first refinement improves every mode whose
    # coarse gap is above the noise floor.  The violations are asserted last so
    # that both residual lines below always print: the store reads them from the
    # log, and a failing monotonicity check used to hide them (T12).
    not_improving = [
        f"finest CCX mode {mode + 1}: {medium} m gap {gaps[medium][mode] * 100:.2f}% is not "
        f"below the {coarse} m gap {gaps[coarse][mode] * 100:.2f}%"
        for mode in identifiable
        if gaps[coarse][mode] > NOISE_FLOOR and not gaps[medium][mode] < gaps[coarse][mode]
    ]

    low_worst = max(gaps[fine][mode] for mode in identifiable if mode < LOW_MODES)
    high_worst = max(gaps[fine][mode] for mode in identifiable if mode >= LOW_MODES)
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
    problems = not_improving + [
        f"finest CCX mode {mode + 1}: not identifiable, pairing MAC {macs[fine][mode]:.3f} "
        f"< {MAC_IDENTIFIABLE} on the {fine} m mesh"
        for mode in not_identifiable
    ]
    assert not problems, "; ".join(problems)
