"""Structural gate for the parked V50 CCX parity deck.

Gate for the campaign documented in ``odd/tasks/ccx-fsi-parity-parked-v50.md``:
the generated CalculiX deck must be a faithful, runnable model of the parked V50
solid *before* any coupled hour is spent on SLURM.  It is deliberately cheap and
runs on a coarse gate mesh; the production 0.25 m deck (``--element-size 0.25``)
is exercised by the coupled campaign itself.

The tool under test (``tools/ccx_fsi_parity_parked_v50.py``) already existed when
this module was written, so these are gate assertions rather than TDD of new
behavior: they pin the artifacts, the adapter cards and the damping mapping the
campaign depends on.
"""

from __future__ import annotations

import importlib.util
import re
import sys
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Any

# The analyzer in this environment resolves against an interpreter outside the
# workspace (the project venv lives in $SCRATCH), so third-party imports it
# cannot see are flagged although they resolve at runtime.
import numpy as np  # type: ignore[import-not-found]
import pytest  # type: ignore[import-not-found]
import yaml

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "tests"))

from tests.support.ccx_io import fail_ccx, parse_ccx_frequencies, run_ccx  # noqa: E402
from conftest import ccx_bin_or_skip  # noqa: E402


def _load_parity_tool() -> Any:
    """Load ``tools/ccx_fsi_parity_parked_v50.py`` (a script, not a package)."""
    path = REPO_ROOT / "tools" / "ccx_fsi_parity_parked_v50.py"
    spec = importlib.util.spec_from_file_location("ccx_fsi_parity_parked_v50", path)
    assert spec is not None and spec.loader is not None, f"cannot load {path}"
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


parity: Any = _load_parity_tool()

#: Coarse gate mesh: the gates catch export/assembly bugs, not discretisation
#: error.  Measured CCX-vs-AeroElast gaps at this mesh are asserted below with a
#: margin justified by the measurement (see MODAL_GAP_BOUND).
GATE_ELEMENT_SIZE = 2.0
N_MODES = 3

#: Measured at GATE_ELEMENT_SIZE (see the test's printed gaps).  The bound is
#: loose on purpose: at 2.0 m the two discretisations genuinely differ.  If this
#: ever trips, the deck export changed, not the physics.
MODAL_GAP_BOUND = 0.20


@pytest.fixture(scope="module")
def gate(tmp_path_factory: pytest.TempPathFactory) -> dict:
    """Generate every parity artifact once on the gate mesh."""
    out_dir = tmp_path_factory.mktemp("ccx_parity_gate")
    return parity.generate_parity_artifacts(out_dir, element_size=GATE_ELEMENT_SIZE)


def test_artifacts_generated_and_parse(gate: dict) -> None:
    paths = gate["paths"]
    for key, path in paths.items():
        assert Path(path).is_file(), f"{key} was not written: {path}"

    config = yaml.safe_load(Path(paths["config_yml"]).read_text())
    interfaces = config["participants"]["Solid"]["interfaces"]
    assert interfaces == [
        {
            "nodes-mesh": "Solid-Mesh",
            "patch": gate["interface_nodeset"],
            "read-data": ["Forces"],
            "write-data": ["Displacements"],
        }
    ]
    assert config["precice-config-file"] == "precice-config.xml"

    xml_text = Path(paths["precice_config_xml"]).read_text()
    # preCICE accepts bare prefixes (``data:``, ``m2n:``, ``coupling-scheme:`` …)
    # and the case files use them without namespace declarations; strict
    # ElementTree needs them removed, so parse a prefix-stripped copy while
    # asserting on the original text.
    stripped = re.sub(r"(</?)[A-Za-z][\w.-]*:", r"\1", xml_text)
    tree = ET.fromstring(stripped)
    data_names = {node.get("name") for node in tree.iter("vector")}
    assert data_names == {"Forces", "Displacements"}
    assert "Velocity" not in xml_text
    exchanges = {node.get("data") for node in tree.iter("exchange")}
    assert exchanges == {"Forces", "Displacements"}

    watch = next(iter(tree.iter("watch-point")))
    coordinate_text = watch.get("coordinate")
    assert coordinate_text is not None, "watch-point has no coordinate attribute"
    coordinate = [float(value) for value in coordinate_text.split(";")]
    assert coordinate == pytest.approx(gate["tip_coordinate"])


def test_deck_contains_adapter_cards(gate: dict) -> None:
    text = Path(gate["paths"]["solid_ccx_inp"]).read_text()
    assert "*STEP, NLGEOM" in text
    assert "*DYNAMIC" in text
    assert f"{gate['dt_s']:.6E}, {gate['t_end_s']:.6E}" in text

    # The adapter requires every interface node loaded in each spatial direction
    # at step start; the values must be the written zeros.
    nset = f"N{gate['interface_nodeset'].upper()}"
    for dof in (1, 2, 3):
        assert f"{nset}, {dof}, 0.0" in text

    alpha, beta = gate["rayleigh_damping_alpha_beta"]
    assert f"*DAMPING, ALPHA={alpha:.6E}, BETA={beta:.6E}" in text


def test_damping_reproduces_target_zeta(gate: dict) -> None:
    """The (ALPHA, BETA) mapping must give the case's 3 % in both damping modes."""
    zeta = 0.03
    frequencies = [
        gate["natural_freq_i_hz"],
        gate["natural_freq_j_hz"],
    ]
    alpha, beta = gate["rayleigh_damping_alpha_beta"]
    for frequency in frequencies:
        omega = 2.0 * np.pi * frequency
        assert alpha / (2.0 * omega) + beta * omega / 2.0 == pytest.approx(zeta, abs=1e-9)

    # ALPHA is the mass-proportional and BETA the stiffness-proportional term.
    assert alpha == pytest.approx(gate["rayleigh_eta_m_mass"])
    assert beta == pytest.approx(gate["rayleigh_eta_k_stiffness"])


def _aeroelast_frequencies(mesh, properties, span_direction, fixed_root_nodes: int) -> np.ndarray:
    """Sorted AeroElast natural frequencies on the same mesh as the deck."""
    from _aeroelast import PyMeshAssembler, modal_solve_coo

    assembler = PyMeshAssembler.from_model(
        parity._to_rust_mesh(mesh, properties), properties, list(span_direction), None
    )
    n_dofs = assembler.dofs_count
    k_rows, k_cols, k_vals = assembler.assemble_k()
    m_rows, m_cols, m_vals = assembler.assemble_m()

    root_ids = parity.require_node_set(mesh, parity.ROOT_NODESET).node_ids
    root_indices = {mesh.node_id_to_index[node_id] for node_id in root_ids}
    fixed = {6 * index + dof for index in root_indices for dof in range(6)}
    free = np.array([i for i in range(n_dofs) if i not in fixed], dtype=np.int64)

    freqs = modal_solve_coo(
        np.asarray(k_rows, dtype=np.int64),
        np.asarray(k_cols, dtype=np.int64),
        np.asarray(k_vals, dtype=np.float64),
        np.asarray(m_rows, dtype=np.int64),
        np.asarray(m_cols, dtype=np.int64),
        np.asarray(m_vals, dtype=np.float64),
        n_dofs,
        free,
        N_MODES,
    )[0]
    assert fixed_root_nodes > 0
    return np.sort(np.asarray(freqs, dtype=float))


def test_modal_gate_ccx_vs_aeroelast(gate: dict, tmp_path: Path) -> None:
    """Same mesh, same clamped root: CCX modal must track AeroElast closely."""
    ccx_bin = ccx_bin_or_skip()
    from aeroelast.core.mesh.io.writers import write_ccx_mesh

    config = parity.load_case_config()
    _, mesh, properties = parity.build_blade_mesh(config, element_size=GATE_ELEMENT_SIZE)

    inp_path = tmp_path / "gate_modal.inp"
    write_ccx_mesh(
        mesh,
        str(inp_path),
        properties=properties,
        boundary_nodeset=parity.ROOT_NODESET,
        solver_type="Modal",
        num_modes=N_MODES,
        # Quadratic S8R: the validated blade parity configuration (AeroElast MITC4
        # vs CCX S8R); S4R falls back to a smeared laminate and is 150 % off.
        quadratic=True,
        span_direction=config["span_direction"],
    )
    result = run_ccx(inp_path, ccx_bin)
    if result.returncode != 0:
        fail_ccx(result, inp_path)
    freqs_ccx = np.sort(parse_ccx_frequencies(inp_path, n_modes=N_MODES))

    root_nodes = len(parity.require_node_set(mesh, parity.ROOT_NODESET).node_ids)
    freqs_ae = _aeroelast_frequencies(mesh, properties, config["span_direction"], root_nodes)

    gaps = np.abs(freqs_ae[:N_MODES] - freqs_ccx[:N_MODES]) / freqs_ccx[:N_MODES]
    print(f"\ngate modal (element_size={GATE_ELEMENT_SIZE} m)")
    print(f"  aeroelast = {np.array2string(freqs_ae[:N_MODES], precision=4)}")
    print(f"  ccx       = {np.array2string(freqs_ccx[:N_MODES], precision=4)}")
    print(f"  gaps      = {np.array2string(100 * gaps, precision=2)} %")

    assert float(np.max(gaps)) < MODAL_GAP_BOUND


def test_dynamic_deck_runs_in_ccx(tmp_path: Path) -> None:
    """The adapter deck must parse and integrate in real CalculiX."""
    ccx_bin = ccx_bin_or_skip()
    gate = parity.generate_parity_artifacts(
        tmp_path, element_size=GATE_ELEMENT_SIZE, dt=1e-2, t_end=2e-2
    )
    inp_path = Path(gate["paths"]["solid_ccx_inp"])
    result = run_ccx(inp_path, ccx_bin)
    if result.returncode != 0:
        fail_ccx(result, inp_path)

    frd = inp_path.with_suffix(".frd")
    assert frd.is_file(), f"no FRD written next to {inp_path}"
    assert "DISP" in frd.read_text(errors="replace")
