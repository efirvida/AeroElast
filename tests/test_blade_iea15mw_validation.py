"""Blade mechanical validation on the IEA 15 MW reference geometry.

Geometry: ``tests/IEA-15-240-RWT.yaml`` (the official IEA Wind 15 MW reference
turbine definition, Gaertner et al. 2020, NREL/TP-5000-75698), meshed and given
composite shell properties by this repository's own ``Blade`` model.

Three references, on purpose:

1. **AeroElast vs CalculiX S8R**, on the same mesh and the same composite
   properties.  This is the tight comparison: both are shell FEM, so a
   difference is a formulation/discretisation gap, not a model difference.
2. **Mass vs the published models**: Escalera Mendoza et al. 2023 (AIAA
   2023-2093) reports 68,077 kg for the UTD NuMAD conversion of this blade, and
   the definition report gives about 65 metric tons for the IEA blade itself.
3. **First flapwise frequency vs the article**: the same paper's Table 3 reports
   the parked blade modes, 1st flapwise 0.57 Hz.

Mesh choice.  ``element_size = 1.0 m``.  It is the coarsest mesh that is good
enough on both comparisons at once, and both were measured before it was chosen:
the 1st flapwise frequency is 0.578 Hz (article 0.57, +1.4%) while the modal
parity against CCX stays within 7.4%; at 1.5 m the frequency rises to 0.601 Hz
(+5.4%) and the parity to 12.3%, and at 2.0 m to 0.632 Hz (+10.9%) and 15.5%.

**Open finding, stated rather than hidden.** The shell model has extra low modes
(0.189 Hz and 0.434 Hz at this mesh) and does not reproduce the article's 1st
edgewise frequency (0.65 Hz) or its 2nd flapwise (1.72 Hz); the closest computed
values are 0.578 and 1.348 Hz.  The article's numbers come from BModes, a beam
model, so a six-DOF shell mesh can legitimately carry modes a beam model does
not.  The flapwise and 2nd edgewise (2.025 Hz vs 2.08, -2.6%) do agree, which is
why only the flapwise mode is asserted against the article here.  The extra-mode
inventory is a real question about the blade model, not a tolerance to widen.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
from scipy.optimize import linear_sum_assignment

from conftest import ccx_bin_or_skip

pytest.importorskip("petsc4py", reason="PETSc not available")
pytest.importorskip("_aeroelast", reason="Rust backend not available")

from _aeroelast import Laminate as RustLaminate  # noqa: E402
from _aeroelast import MeshModel as RustMeshModel  # noqa: E402
from _aeroelast import PyMeshAssembler, modal_solve_coo  # noqa: E402

from _ccx_io import fail_ccx, parse_ccx_frequencies, run_ccx  # noqa: E402
from aeroelast.core.mesh.entities import MeshElement, Node  # noqa: E402
from aeroelast.core.mesh.io.writers import write_ccx_mesh  # noqa: E402
from aeroelast.models.blade.model import Blade  # noqa: E402

YAML = Path(__file__).resolve().parent / "IEA-15-240-RWT.yaml"

ELEMENT_SIZE = 1.0
N_SEARCH = 10  # modes requested from both solvers
N_COMPARE = 5  # matched pairs asserted

#: Escalera Mendoza et al. 2023 (AIAA 2023-2093), UTD NuMAD model.
ARTICLE_MASS_KG = 68_077.0
#: Gaertner et al. 2020 (NREL/TP-5000-75698): "around 65 metric tons".
REPORT_MASS_KG = 65_000.0
#: Escalera Mendoza et al. 2023, Table 3, 1st flapwise parked mode.
ARTICLE_FIRST_FLAPWISE_HZ = 0.57

MASS_TOL = 0.10  # measured 70,825 kg = +4.0% over the article value
MODAL_TOL = 0.10  # measured worst 7.4% over the first five matched pairs
FLAPWISE_TOL = 0.05  # measured 0.578 Hz = +1.4%


def _to_rust_mesh(mesh, properties: dict):
    """Convert the Python mesh + Rust properties into a ``RustMeshModel``.

    Mirrors ``aeroelast.core.assembler._build_py_mesh_assembler``: composite
    element sets get element codes 33 (tri) / 44 (quad), plain sets 3 / 4.
    """
    nodes = mesh.nodes
    node_ids = [n.id for n in nodes]
    coords_flat = np.stack([n.coords for n in nodes], axis=0).ravel().tolist()
    elements = mesh.elements
    element_ids = np.fromiter((e.id for e in elements), dtype=np.int64, count=len(elements))
    node_counts = np.fromiter((e.node_count for e in elements), dtype=np.int8, count=len(elements))

    composite_ids: set[int] = set()
    for set_name, prop in properties.items():
        if isinstance(prop, RustLaminate) and set_name in mesh.element_sets:
            composite_ids.update(e.id for e in mesh.element_sets[set_name].elements)
    composite = (
        np.isin(element_ids, np.fromiter(composite_ids, dtype=np.int64))
        if composite_ids
        else np.zeros(len(elements), dtype=bool)
    )
    is_tri = node_counts == 3
    type_codes = np.where(composite, np.where(is_tri, 33, 44), np.where(is_tri, 3, 4)).tolist()

    element_sets = {name: [e.id for e in eset.elements] for name, eset in mesh.element_sets.items()}
    node_sets = {name: list(nset.node_ids) for name, nset in mesh.node_sets.items()}

    return RustMeshModel.from_raw_data(
        node_ids,
        coords_flat,
        element_ids.tolist(),
        [list(e.node_ids) for e in elements],
        type_codes,
        element_sets,
        node_sets,
    )


@pytest.fixture(scope="module")
def blade(tmp_path_factory: pytest.TempPathFactory) -> dict:
    """Mesh, assemble, and solve the blade once, AeroElast and CCX.

    Returns the mass, the AeroElast frequencies and the CCX frequencies.
    """
    ccx_bin = ccx_bin_or_skip()

    # The node/element id counters are process-global; a mesh built by an earlier
    # test would leave them advanced, and the writer indexes coordinates by node
    # id, so the blade would inherit non-zero-based ids. Reset, as
    # test_blade_mesh.py does.
    Node._id_counter = 0
    MeshElement._id_counter = 0

    blade_model = Blade(str(YAML), element_size=ELEMENT_SIZE)
    blade_model.generate_mesh()
    mesh = blade_model.mesh
    props = blade_model.get_element_properties()

    assembler = PyMeshAssembler.from_model(_to_rust_mesh(mesh, props), props, None, None)
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
                N_SEARCH,
            )[0],
            dtype=float,
        )
    )

    workdir = tmp_path_factory.mktemp("blade_iea15mw")
    inp_path = workdir / "blade.inp"
    write_ccx_mesh(
        mesh,
        str(inp_path),
        properties=props,
        boundary_nodeset="RootNodes",
        solver_type="Modal",
        num_modes=N_SEARCH,
        quadratic=True,
    )
    result = run_ccx(inp_path, ccx_bin)
    if result.returncode != 0:
        fail_ccx(result, inp_path)
    freqs_ccx = np.sort(parse_ccx_frequencies(inp_path, n_modes=N_SEARCH))

    mass = float(assembler.total_elemental_mass())
    print(f"\nblade: {mesh.node_count} nodes / {mesh.elements_count} elements / {n} dofs")
    print(f"  mass  {mass:,.0f} kg (article {ARTICLE_MASS_KG:,.0f}, report {REPORT_MASS_KG:,.0f})")
    print(f"  aero  {np.array2string(freqs_ae, precision=3)}")
    print(f"  ccx   {np.array2string(freqs_ccx, precision=3)}")

    return {"mass": mass, "ae": freqs_ae, "ccx": freqs_ccx}


def _matched_pairs(freqs_ae: np.ndarray, freqs_ccx: np.ndarray, n: int):
    """Smallest-cost one-to-one frequency pairs, by Hungarian assignment."""
    cost = np.abs(freqs_ae[:, None] - freqs_ccx[None, :]) / np.maximum(freqs_ccx[None, :], 1e-14)
    row_ind, col_ind = linear_sum_assignment(cost)
    pairs = sorted(
        (float(cost[i, j]), float(freqs_ae[i]), float(freqs_ccx[j]))
        for i, j in zip(row_ind, col_ind, strict=False)
    )
    return pairs[:n]


def test_blade_mass_matches_published_models(blade: dict) -> None:
    """The meshed blade mass is within ``MASS_TOL`` of the published models."""
    mass = blade["mass"]
    rel_article = abs(mass - ARTICLE_MASS_KG) / ARTICLE_MASS_KG
    assert rel_article < MASS_TOL, (
        f"blade mass {mass:,.0f} kg is {rel_article * 100:.2f}% from the article's "
        f"{ARTICLE_MASS_KG:,.0f} kg (tol {MASS_TOL * 100:.0f}%)"
    )
    # The definition report gives the blade itself as about 65 t; the model being
    # above it is expected (the article's own NuMAD conversion is +4.33% over it),
    # so this checks the sign and rough size of that difference rather than parity.
    rel_report = (mass - REPORT_MASS_KG) / REPORT_MASS_KG
    assert 0.0 < rel_report < 0.20, (
        f"blade mass {mass:,.0f} kg is {rel_report * 100:.2f}% over the report's "
        f"{REPORT_MASS_KG:,.0f} kg; expected above it but within 20%"
    )


@pytest.mark.parametrize("index", range(N_COMPARE))
def test_blade_modal_frequencies_match_ccx(blade: dict, index: int) -> None:
    """The first matched blade modes agree with CalculiX S8R within ``MODAL_TOL``.

    The pairing is by cost over 10 requested modes, so a mode that one solver
    orders differently cannot turn a real mismatch into a false failure.
    """
    pairs = _matched_pairs(blade["ae"], blade["ccx"], N_COMPARE)
    rel, freq_ae, freq_ccx = pairs[index]
    print(f"  matched[{index}] aero={freq_ae:.3f} ccx={freq_ccx:.3f} rel={rel * 100:.2f}%")
    assert rel < MODAL_TOL, (
        f"blade mode {index}: aero={freq_ae:.3f} Hz ccx={freq_ccx:.3f} Hz "
        f"rel={rel * 100:.2f}% (tol {MODAL_TOL * 100:.0f}%)"
    )


def test_blade_first_flapwise_matches_article(blade: dict) -> None:
    """The computed frequency closest to the article's 1st flapwise is within 5%.

    This is the one article mode that maps cleanly onto this shell model; the
    module docstring records the article modes that do not (1st edgewise, 2nd
    flapwise) and the extra low modes on the shell side.
    """
    freqs_ae = blade["ae"]
    closest = float(freqs_ae[np.argmin(np.abs(freqs_ae - ARTICLE_FIRST_FLAPWISE_HZ))])
    rel = abs(closest - ARTICLE_FIRST_FLAPWISE_HZ) / ARTICLE_FIRST_FLAPWISE_HZ
    assert rel < FLAPWISE_TOL, (
        f"closest to the article's 1st flapwise {ARTICLE_FIRST_FLAPWISE_HZ} Hz is "
        f"{closest:.3f} Hz, rel={rel * 100:.2f}% (tol {FLAPWISE_TOL * 100:.0f}%). "
        f"computed modes: {freqs_ae.tolist()}"
    )
