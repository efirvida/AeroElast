"""Generate the CalculiX/preCICE parity artifacts for the parked V50 FSI case.

The parked V50 case (``tests/IEA15MW/parked_v50``) is run twice against the same
BEM fluid participant: once with the aeroelast solid and once with
CalculiX+preCICE.  This tool materialises the CalculiX-side inputs so the two
runs are apples-to-apples:

* ``solid_ccx.inp`` — DynamicFSI deck for ``ccx_preCICE`` (composite layup,
  clamped root, zero ``*CLOAD`` on the FSI interface, Rayleigh damping).
* ``config.yml`` — CalculiX-preCICE adapter interface descriptor.
* ``precice-config.xml`` — parity XML variant: data renamed to the CCX adapter
  keywords, ``Velocity`` removed, ``max-time`` set to ``t_end`` and a
  ``Flap-Tip`` watch-point on the interface tip node.
* ``solid_v50_parity.yaml`` / ``fluid_v50_parity.yaml`` — copies of the case
  YAMLs with only the parity-specific coupling/output changes.

Nothing executes on import; the module exposes pure functions plus a ``main``
entry point.  The smooth path is :func:`generate_parity_artifacts`.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path
from typing import Any

import yaml

#: Repository root, resolved from this file (``<repo>/tools/<this file>``).
REPO_ROOT = Path(__file__).resolve().parents[1]
#: Default case directory shipped with the repository.
CASE_DIR = REPO_ROOT / "tests" / "IEA15MW" / "parked_v50"

DEFAULT_ELEMENT_SIZE = 0.25
DEFAULT_DT = 1.0e-2
DEFAULT_T_END = 20.0

#: Node-set names the case relies on.  ``BladeMesh`` derives ``allOuterShellNods``
#: from the ``allOuterShellEls`` element set (see ``BladeMesh`` node-set creation)
#: and ``RootNodes`` from the section with minimum span.
INTERFACE_NODESET = "allOuterShellNods"
ROOT_NODESET = "RootNodes"

#: preCICE/CCX-adapter keyword pairs (singular case YAML -> plural adapter data).
_FORCE_DATA = "Forces"
_DISPLACEMENT_DATA = "Displacements"


# ---------------------------------------------------------------------------
# Case configuration
# ---------------------------------------------------------------------------


def _resolve_case_path(case_dir: Path, value: str) -> Path:
    """Resolve a case-relative YAML path against the case directory."""
    candidate = Path(value)
    if candidate.is_absolute():
        return candidate
    return (case_dir / candidate).resolve()


def load_case_config(case_dir: Path | str = CASE_DIR) -> dict[str, Any]:
    """Read the parked V50 case YAMLs and resolve the paths this tool needs.

    Mesh-generator paths (``yaml_file``, ``airfoil_dir``) are resolved relative
    to the case directory, exactly as the runtime does.  Returns a flat dict so
    callers do not need to know the case file layout.
    """
    case_dir = Path(case_dir).resolve()
    solid_path = case_dir / "solid_v50.yaml"
    fluid_path = case_dir / "fluid_v50.yaml"

    with open(solid_path) as handle:
        solid = yaml.safe_load(handle)
    with open(fluid_path) as handle:
        fluid = yaml.safe_load(handle)

    generator = solid["mesh"]["generator"]["params"]
    damping = solid["solver"]["damping"]
    boundaries = list(solid["coupling"]["boundaries"])

    return {
        "case_dir": case_dir,
        "solid_yaml_path": solid_path,
        "fluid_yaml_path": fluid_path,
        "precice_xml_path": case_dir / "precice-config.xml",
        "yaml_file": _resolve_case_path(case_dir, generator["yaml_file"]),
        "airfoil_dir": _resolve_case_path(case_dir, generator["airfoil_dir"]),
        "n_samples": int(generator.get("n_samples", 300)),
        "span_grading": str(generator.get("span_grading", "chord")),
        "element_size": float(generator.get("element_size", DEFAULT_ELEMENT_SIZE)),
        "span_direction": tuple(float(v) for v in solid["elements"]["span_direction"]),
        "damping": {
            "zeta": float(damping["zeta"]),
            "zeta_1": float(damping["zeta_1"]) if damping.get("zeta_1") is not None else None,
            "zeta_2": float(damping["zeta_2"]) if damping.get("zeta_2") is not None else None,
            "mode_i": int(damping.get("mode_i", 1)),
            "mode_j": int(damping.get("mode_j", 2)),
            "num_modes": int(damping.get("num_modes", 10)),
        },
        "interface_nodeset": boundaries[0] if boundaries else INTERFACE_NODESET,
        "participant": str(solid["coupling"].get("participant", "Solid")),
        "coupling_mesh": str(solid["coupling"].get("coupling_mesh", "Solid-Mesh")),
        "fluid": fluid,
    }


# ---------------------------------------------------------------------------
# Mesh + properties (same path as tests/test_blade_iea15mw_mesh_convergence.py)
# ---------------------------------------------------------------------------


def build_blade_mesh(
    config: dict[str, Any], element_size: float | None = None
) -> tuple[Any, Any, dict]:
    """Build the blade mesh and its element properties.

    Mirrors ``tests/test_blade_iea15mw_mesh_convergence.py``: reset the global
    entity-id counters, construct :class:`~aeroelast.models.blade.model.Blade`
    from the resolved blade YAML, generate the mesh and read the per-element-set
    properties (Rust-native laminates).

    Returns ``(blade_model, mesh, properties)``.
    """
    from aeroelast.core.mesh.entities import MeshElement, Node
    from aeroelast.models.blade.model import Blade

    Node._id_counter = 0
    MeshElement._id_counter = 0

    size = float(config["element_size"] if element_size is None else element_size)
    blade = Blade(str(config["yaml_file"]), element_size=size, n_samples=int(config["n_samples"]))
    blade.generate_mesh()
    return blade, blade.mesh, blade.get_element_properties()


def _to_rust_mesh(mesh: Any, properties: dict) -> Any:
    """Convert a Python :class:`MeshModel` to the Rust ``MeshModel``.

    Composite element sets (Rust laminates) get element codes 33/44, plain sets
    3/4.  Identical to the conversion used by the blade validation/parity tests.
    """
    import numpy as np
    from _aeroelast import Laminate as RustLaminate
    from _aeroelast import MeshModel as RustMeshModel

    nodes = mesh.nodes
    node_ids = [node.id for node in nodes]
    coords_flat = np.stack([node.coords for node in nodes], axis=0).ravel().tolist()
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


def require_node_set(mesh: Any, name: str) -> Any:
    """Return the named node set, raising a clear error listing what exists."""
    if name not in mesh.node_sets:
        available = sorted(mesh.node_sets.keys())
        raise ValueError(f"Mesh has no node set {name!r}. Available node sets: {available}")
    return mesh.get_node_set(name)


def interface_tip_node(mesh: Any, interface_nodeset: str) -> Any:
    """Return the interface node with the maximum span (Z) coordinate."""
    nodes = mesh.nodes
    index_of = mesh.node_id_to_index
    node_set = require_node_set(mesh, interface_nodeset)
    interface_indices = [index_of[node_id] for node_id in node_set.node_ids]
    return max((nodes[i] for i in interface_indices), key=lambda node: node.z)


# ---------------------------------------------------------------------------
# Rayleigh damping via the solver's Rust fast path
# ---------------------------------------------------------------------------


def compute_rayleigh(
    mesh: Any,
    properties: dict,
    span_direction: tuple[float, float, float],
    damping: dict[str, Any],
    rayleigh_override: tuple[float, float] | None = None,
) -> dict[str, Any]:
    """Compute Rayleigh coefficients exactly as the solver does.

    Uses ``PyMeshAssembler`` -> ``assemble_k``/``assemble_m`` -> COO, eliminates
    the clamped ``RootNodes`` DOFs (fixed = ``6*node + dof``) and calls
    ``_aeroelast.compute_rayleigh_auto``.  Also solves a short modal problem so
    the summary can report the natural frequencies of the two damping modes.

    Returns a dict with the stiffness/mass coefficients (``eta_k``/``eta_m``),
    the ``rayleigh_damping`` tuple in CCX ``(alpha, beta)`` order and the two
    target natural frequencies in Hz.
    """
    if rayleigh_override is not None:
        # CCX order: alpha is the mass-proportional and beta the stiffness-
        # proportional coefficient.  Timing runs do not need the modal solve,
        # which is the dominant cost at the production mesh size.
        alpha, beta = float(rayleigh_override[0]), float(rayleigh_override[1])
        return {
            "eta_k": beta,
            "eta_m": alpha,
            "source": "override",
            "rayleigh_damping": (alpha, beta),
            "mode_i": int(damping["mode_i"]),
            "mode_j": int(damping["mode_j"]),
            "freq_i_hz": None,
            "freq_j_hz": None,
            "num_modes": 0,
        }

    import numpy as np
    from _aeroelast import PyMeshAssembler, compute_rayleigh_auto, modal_solve_coo

    assembler = PyMeshAssembler.from_model(
        _to_rust_mesh(mesh, properties), properties, list(span_direction), None
    )
    n_dofs = assembler.dofs_count
    k_rows, k_cols, k_vals = assembler.assemble_k()
    m_rows, m_cols, m_vals = assembler.assemble_m()

    # ``compute_rayleigh_auto`` binds the COO indices as int32 (the solver's
    # PETSc-to-COO path yields int32); ``modal_solve_coo`` uses int64.
    k_rows_i32 = np.asarray(k_rows, dtype=np.int32)
    k_cols_i32 = np.asarray(k_cols, dtype=np.int32)
    k_vals_f64 = np.asarray(k_vals, dtype=np.float64)
    m_rows_i32 = np.asarray(m_rows, dtype=np.int32)
    m_cols_i32 = np.asarray(m_cols, dtype=np.int32)
    m_vals_f64 = np.asarray(m_vals, dtype=np.float64)

    root_ids = require_node_set(mesh, ROOT_NODESET).node_ids
    root_indices = {mesh.node_id_to_index[node_id] for node_id in root_ids}
    fixed = {6 * index + dof for index in root_indices for dof in range(6)}
    free = np.array([i for i in range(n_dofs) if i not in fixed], dtype=np.int64)
    free_i32 = free.astype(np.int32)

    zeta = float(damping["zeta"])
    zeta_i = float(damping["zeta_1"]) if damping.get("zeta_1") is not None else zeta
    zeta_j = float(damping["zeta_2"]) if damping.get("zeta_2") is not None else zeta
    mode_i = int(damping["mode_i"])
    mode_j = int(damping["mode_j"])
    num_modes = int(damping["num_modes"])

    eta_k: float | None
    eta_m: float | None
    rayleigh_source = "rust"
    try:
        eta_k, eta_m = compute_rayleigh_auto(
            k_rows_i32,
            k_cols_i32,
            k_vals_f64,
            m_rows_i32,
            m_cols_i32,
            m_vals_f64,
            free_i32,
            num_modes,
            mode_i,
            mode_j,
            zeta_i,
            zeta_j,
        )
    except Exception as exc:  # noqa: BLE001 — PETSc/SLEPc failures are opaque
        # The Rust fast path runs the same GHEP through the PETSc/SLEPc
        # eigenvalue solve, which fails on some hosts (seen: PETSc error 95 in
        # context 'EPSSolve').  ``modal_solve_coo`` below solves the same
        # problem and is the path the CCX parity tests use, so recompute the
        # coefficients from its frequencies with the solver's own closed-form
        # formula (linear_dynamic.py:524-534).
        print(
            f"  [warn] compute_rayleigh_auto failed ({exc}); "
            "falling back to modal_solve_coo + closed form"
        )
        eta_k = eta_m = None
        rayleigh_source = "closed-form"

    freqs = np.sort(
        np.asarray(
            modal_solve_coo(
                k_rows_i32.astype(np.int64),
                k_cols_i32.astype(np.int64),
                k_vals_f64,
                m_rows_i32.astype(np.int64),
                m_cols_i32.astype(np.int64),
                m_vals_f64,
                n_dofs,
                free,
                num_modes,
            )[0],
            dtype=float,
        )
    )

    if eta_k is None or eta_m is None:
        omega = 2.0 * np.pi * np.asarray(freqs, dtype=float)
        needed = max(mode_i, mode_j)
        if len(omega) < needed:
            raise RuntimeError(
                f"Rayleigh: modal solve returned {len(omega)} modes but mode {needed} is required."
            )
        omega_i_f = float(omega[mode_i - 1])
        omega_j_f = float(omega[mode_j - 1])
        denom = omega_i_f**2 - omega_j_f**2
        if abs(denom) < 1e-12:
            raise ValueError(
                f"Rayleigh: modes {mode_i} and {mode_j} share the same natural "
                f"frequency (omega ~= {omega_i_f:.4e} rad/s)."
            )
        eta_k = float(2.0 * (zeta_i * omega_i_f - zeta_j * omega_j_f) / denom)
        eta_m = float(
            2.0 * omega_i_f * omega_j_f * (zeta_j * omega_i_f - zeta_i * omega_j_f) / denom
        )

    def mode_frequency(index: int) -> float | None:
        return float(freqs[index - 1]) if 1 <= index <= len(freqs) else None

    return {
        "eta_k": float(eta_k),
        "eta_m": float(eta_m),
        "source": rayleigh_source,
        # CCX *DAMPING uses ALPHA=mass, BETA=stiffness.
        "rayleigh_damping": (float(eta_m), float(eta_k)),
        "mode_i": mode_i,
        "mode_j": mode_j,
        "freq_i_hz": mode_frequency(mode_i),
        "freq_j_hz": mode_frequency(mode_j),
        "num_modes": int(len(freqs)),
    }


# ---------------------------------------------------------------------------
# Deck + adapter config
# ---------------------------------------------------------------------------


def write_ccx_deck(
    mesh: Any,
    properties: dict,
    out_path: Path,
    *,
    interface_nodeset: str,
    span_direction: tuple[float, float, float],
    rayleigh_damping: tuple[float, float],
    dt: float,
    t_end: float,
) -> None:
    """Write the ``DynamicFSI`` CalculiX deck for the preCICE adapter."""
    from aeroelast.core.mesh.io.writers import write_ccx_mesh

    require_node_set(mesh, interface_nodeset)
    write_ccx_mesh(
        mesh,
        str(out_path),
        properties=properties,
        boundary_nodeset=ROOT_NODESET,
        solver_type="DynamicFSI",
        fsi_interface_nodeset=interface_nodeset,
        nlgeom=True,
        rayleigh_damping=rayleigh_damping,
        output_frequency=1,
        # S8R + a real *SHELL SECTION, COMPOSITE.  The linear (S4R) path falls
        # back to the smeared equivalent laminate material and loses the bending
        # stiffness of the layup: the modal gate measured 150 % gaps against
        # AeroElast's MITC4Composite, against <2 % for the quadratic element —
        # the same choice the validated blade parity tests make.
        quadratic=True,
        span_direction=span_direction,
        dt=dt,
        t_end=t_end,
    )


def write_adapter_config(
    out_path: Path,
    *,
    nodes_mesh: str,
    patch: str,
    precice_config_file: str = "precice-config.xml",
) -> None:
    """Write the CalculiX-preCICE adapter ``config.yml``."""
    text = (
        "# CalculiX-preCICE adapter interface (parity parked V50)\n"
        "participants:\n"
        "  Solid:\n"
        "    interfaces:\n"
        f'      - nodes-mesh: "{nodes_mesh}"\n'
        f'        patch: "{patch}"\n'
        f"        read-data: [{_FORCE_DATA}]\n"
        f"        write-data: [{_DISPLACEMENT_DATA}]\n"
        f"precice-config-file: {precice_config_file}\n"
    )
    out_path.write_text(text)


# ---------------------------------------------------------------------------
# preCICE XML parity variant
# ---------------------------------------------------------------------------


def derive_precice_xml(
    src_path: Path,
    out_path: Path,
    *,
    t_end: float,
    tip_coordinate: tuple[float, float, float],
    interface_nodeset: str = INTERFACE_NODESET,
    coupling_mesh: str = "Solid-Mesh",
) -> None:
    """Derive the parity preCICE XML from the case XML.

    The source XML uses prefixed tags (``data:vector``) without namespace
    declarations, so it is not parseable by :mod:`xml.etree`; this performs the
    known, deterministic text transforms instead.  Changes:

    * ``Force``/``Displacement`` renamed to the CCX adapter keywords
      ``Forces``/``Displacements`` (data declarations, ``use-data``,
      ``read-data``/``write-data`` and ``exchange``).
    * ``Velocity`` declaration, comments, ``use-data`` and exchange removed.
    * ``max-time`` set to ``t_end``.
    * A ``Flap-Tip`` watch-point added to participant ``Solid`` on the interface
      node with maximum Z.

    The coupling scheme, IQN-ILS acceleration, convergence limits,
    ``max-iterations`` and sockets are left untouched; ``exchange-directory`` is
    rewritten to ``.`` so the parity run directory is self-contained (both
    participants are launched from it).
    """
    text = src_path.read_text()

    text = text.replace('name="Displacement"', f'name="{_DISPLACEMENT_DATA}"')
    text = text.replace('data="Displacement"', f'data="{_DISPLACEMENT_DATA}"')
    text = text.replace('name="Force"', f'name="{_FORCE_DATA}"')
    text = text.replace('data="Force"', f'data="{_FORCE_DATA}"')

    # Drop the Velocity comment block, then every Velocity element line.
    text = re.sub(r"\n\s*<!--\s*Velocity:.*?-->\n", "\n", text, flags=re.DOTALL)
    text = re.sub(
        r"(?m)^\s*<(?:data:vector|use-data|read-data|write-data|exchange)\b"
        r'[^>]*(?:name|data)="Velocity"[^>]*/>\s*\n',
        "",
        text,
    )

    # Coupling window end time.
    text = re.sub(
        r'(<max-time value=")[^"]*(")',
        lambda match: f"{match.group(1)}{t_end}{match.group(2)}",
        text,
        count=1,
    )

    # Self-contained run directory: both participants are launched from the run
    # directory, so "." resolves to the same path for both.
    text = re.sub(
        r'(<m2n:sockets[^>]*exchange-directory=")[^"]*(")',
        lambda match: f"{match.group(1)}.{match.group(2)}",
        text,
        count=1,
    )

    # Interface tip watch-point, added at the top of participant Solid.
    coordinate = ";".join(repr(float(component)) for component in tip_coordinate)
    watch_point = (
        f'<watch-point mesh="{coupling_mesh}" name="Flap-Tip" coordinate="{coordinate}" />'
    )
    text = re.sub(
        r'(?m)^(\s*)<participant name="Solid">\n',
        lambda match: f"{match.group(0)}{match.group(1)}  {watch_point}\n",
        text,
        count=1,
    )

    if "Velocity" in text and 'name="Velocity"' in text:
        raise AssertionError(f"Velocity data still present in derived XML: {src_path}")

    out_path.write_text(text)


# ---------------------------------------------------------------------------
# Parity YAML variants
# ---------------------------------------------------------------------------

_GENERATOR_PARAM_KEYS = (
    "yaml_file:",
    "airfoil_dir:",
    "element_size:",
    "n_samples:",
    "span_grading:",
)

#: Geometry paths that a run directory outside the repository must resolve.
_GEOMETRY_PATH_KEYS = ("yaml_file", "airfoil_dir", "blade_file")

_GEOMETRY_PATH_RE = re.compile(
    r"(?m)^(\s*)(" + "|".join(_GEOMETRY_PATH_KEYS) + r'):\s*"([^"]+)"'
)


def _absolutise_geometry_paths(text: str, src_path: Path) -> str:
    """Rewrite the geometry paths to absolute (run dirs live outside the repo)."""

    def replace(match: re.Match[str]) -> str:
        indent, key, value = match.group(1), match.group(2), match.group(3)
        return f'{indent}{key}: "{_resolve_case_path(src_path.parent, value)}"'

    return _GEOMETRY_PATH_RE.sub(replace, text)


def _assert_geometry_paths_absolute(src_text: str, out_text: str, case_dir: Path) -> None:
    """Fail loudly if a geometry path did not end up absolute in the output."""
    for match in _GEOMETRY_PATH_RE.finditer(src_text):
        resolved = str(_resolve_case_path(case_dir, match.group(3)))
        if resolved not in out_text:
            raise AssertionError(f"Geometry path was not rewritten to {resolved!r}")


def _assert_generator_params_unchanged(src_text: str, out_text: str) -> None:
    """Fail loudly if a mesh-generator param line did not survive byte-identically."""
    for line in src_text.splitlines():
        if any(key in line for key in _GENERATOR_PARAM_KEYS) and line not in out_text:
            raise AssertionError(f"Mesh generator param line changed: {line!r}")


def write_solid_parity_yaml(
    src_path: Path, out_path: Path, out_dir: Path, *, abs_paths: bool = False
) -> None:
    """Write the solid parity YAML: plural data names, no ramp, out-dir output.

    ``abs_paths`` rewrites the geometry paths to absolute so the YAML works from
    a run directory outside the repository (``_GENERATOR_PARAM_KEYS``-style
    byte-identity is then asserted on the resolved values instead).
    """
    text = src_path.read_text()
    text = re.sub(r'(?m)^(\s*write_data:\s*)"Displacement"', r'\1"Displacements"', text)
    text = re.sub(r'(?m)^(\s*read_data:\s*)"Force"', r'\1"Forces"', text)
    text = re.sub(r'(?m)^(\s*config_file:\s*)"[^"]*"', r'\1"precice-config.xml"', text)
    # CCX cannot ramp the solid-side load, so both parity runs run unramped.
    text = re.sub(
        r"(?m)^\s*force_ramp_time:\s*[0-9.eE+-]+.*$",
        "  force_ramp_time: 0.0   # CCX cannot ramp solid-side; both parity runs unramped",
        text,
    )
    folder = str((out_dir / "results").resolve())
    text = re.sub(
        r'(?m)^(\s*folder:\s*)"[^"]*"',
        lambda match: f'{match.group(1)}"{folder}"',
        text,
    )
    if abs_paths:
        text = _absolutise_geometry_paths(text, src_path)
        _assert_geometry_paths_absolute(src_path.read_text(), text, src_path.parent)
    else:
        _assert_generator_params_unchanged(src_path.read_text(), text)
    out_path.write_text(text)


def write_fluid_parity_yaml(
    src_path: Path, out_path: Path, out_dir: Path, *, abs_paths: bool = False
) -> None:
    """Write the fluid parity YAML: plural data names, no velocity, out-dir output."""
    text = src_path.read_text()
    if "force_data:" not in text:
        text = re.sub(
            r"(?m)^(velocity_data:.*)$",
            'force_data:          "Forces"\n'
            'displacement_data:   "Displacements"\n'
            r"\1",
            text,
        )
    text = re.sub(r"(?m)^(velocity_data:\s*)\S+", r"\1null", text)
    text = re.sub(r'(?m)^(\s*config_file:\s*)"[^"]*"', r'\1"precice-config.xml"', text)
    folder = str((out_dir / "bem_fsi_results").resolve())
    text = re.sub(
        r'(?m)^(\s*folder:\s*)"[^"]*"',
        lambda match: f'{match.group(1)}"{folder}"',
        text,
    )
    if abs_paths:
        text = _absolutise_geometry_paths(text, src_path)
        _assert_geometry_paths_absolute(src_path.read_text(), text, src_path.parent)
    else:
        _assert_generator_params_unchanged(src_path.read_text(), text)
    out_path.write_text(text)


# ---------------------------------------------------------------------------
# Orchestration
# ---------------------------------------------------------------------------


def generate_parity_artifacts(
    out_dir: Path | str,
    *,
    element_size: float = DEFAULT_ELEMENT_SIZE,
    dt: float = DEFAULT_DT,
    t_end: float = DEFAULT_T_END,
    case_dir: Path | str = CASE_DIR,
    rayleigh_override: tuple[float, float] | None = None,
    abs_paths: bool = False,
) -> dict[str, Any]:
    """Generate every parity artifact and return a summary dict.

    ``rayleigh_override`` is a CCX-order ``(alpha, beta)`` pair that skips the
    modal solve, and ``abs_paths`` makes the parity YAMLs usable from a run
    directory outside the repository.
    """
    out_dir = Path(out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    config = load_case_config(case_dir)
    blade, mesh, properties = build_blade_mesh(config, element_size=element_size)

    interface_nodeset = config["interface_nodeset"]
    require_node_set(mesh, interface_nodeset)
    tip = interface_tip_node(mesh, interface_nodeset)
    tip_coordinate = (float(tip.x), float(tip.y), float(tip.z))

    rayleigh = compute_rayleigh(
        mesh,
        properties,
        config["span_direction"],
        config["damping"],
        rayleigh_override,
    )

    solid_inp = out_dir / "solid_ccx.inp"
    write_ccx_deck(
        mesh,
        properties,
        solid_inp,
        interface_nodeset=interface_nodeset,
        span_direction=config["span_direction"],
        rayleigh_damping=rayleigh["rayleigh_damping"],
        dt=dt,
        t_end=t_end,
    )

    config_yml = out_dir / "config.yml"
    write_adapter_config(
        config_yml,
        nodes_mesh=config["coupling_mesh"],
        patch=interface_nodeset,
    )

    precice_xml = out_dir / "precice-config.xml"
    derive_precice_xml(
        config["precice_xml_path"],
        precice_xml,
        t_end=t_end,
        tip_coordinate=tip_coordinate,
        interface_nodeset=interface_nodeset,
        coupling_mesh=config["coupling_mesh"],
    )

    solid_parity_yaml = out_dir / "solid_v50_parity.yaml"
    fluid_parity_yaml = out_dir / "fluid_v50_parity.yaml"
    write_solid_parity_yaml(
        config["solid_yaml_path"], solid_parity_yaml, out_dir, abs_paths=abs_paths
    )
    write_fluid_parity_yaml(
        config["fluid_yaml_path"], fluid_parity_yaml, out_dir, abs_paths=abs_paths
    )

    return {
        "out_dir": str(out_dir),
        "element_size_m": float(element_size),
        "dt_s": float(dt),
        "t_end_s": float(t_end),
        "mesh_nodes": int(mesh.node_count),
        "mesh_elements": int(mesh.elements_count),
        "node_sets": sorted(mesh.node_sets.keys()),
        "interface_nodeset": interface_nodeset,
        "tip_node_id": int(tip.id),
        "tip_coordinate": list(tip_coordinate),
        "rayleigh_eta_k_stiffness": rayleigh["eta_k"],
        "rayleigh_eta_m_mass": rayleigh["eta_m"],
        "rayleigh_damping_alpha_beta": list(rayleigh["rayleigh_damping"]),
        "damping_mode_i": rayleigh["mode_i"],
        "damping_mode_j": rayleigh["mode_j"],
        "natural_freq_i_hz": rayleigh["freq_i_hz"],
        "natural_freq_j_hz": rayleigh["freq_j_hz"],
        "modal_frequencies_hz": rayleigh["num_modes"],
        "case_dir": str(config["case_dir"]),
        "paths": {
            "solid_ccx_inp": str(solid_inp),
            "config_yml": str(config_yml),
            "precice_config_xml": str(precice_xml),
            "solid_v50_parity_yaml": str(solid_parity_yaml),
            "fluid_v50_parity_yaml": str(fluid_parity_yaml),
        },
    }


def main(argv: list[str] | None = None) -> int:
    """CLI entry point."""
    parser = argparse.ArgumentParser(
        prog="ccx_fsi_parity_parked_v50",
        description=("Generate the CalculiX/preCICE parity artifacts for the parked V50 FSI case."),
    )
    parser.add_argument("--out-dir", required=True, type=Path, help="Output directory")
    parser.add_argument(
        "--element-size",
        type=float,
        default=DEFAULT_ELEMENT_SIZE,
        help=f"Blade mesh target element size in metres (default {DEFAULT_ELEMENT_SIZE})",
    )
    parser.add_argument(
        "--dt",
        type=float,
        default=DEFAULT_DT,
        help=f"Dynamic step size in seconds (default {DEFAULT_DT:g})",
    )
    parser.add_argument(
        "--t-end",
        type=float,
        default=DEFAULT_T_END,
        help=f"Simulation end time in seconds (default {DEFAULT_T_END:g})",
    )
    parser.add_argument(
        "--alpha",
        type=float,
        default=None,
        help="CCX mass-proportional Rayleigh coefficient; skips the modal solve",
    )
    parser.add_argument(
        "--beta",
        type=float,
        default=None,
        help="CCX stiffness-proportional Rayleigh coefficient; skips the modal solve",
    )
    parser.add_argument(
        "--abs-paths",
        action="store_true",
        help="Rewrite the parity YAML geometry paths to absolute (run dirs outside the repo)",
    )
    args = parser.parse_args(argv)

    if (args.alpha is None) != (args.beta is None):
        parser.error("--alpha and --beta must be given together")
    override = None if args.alpha is None else (float(args.alpha), float(args.beta))

    summary = generate_parity_artifacts(
        args.out_dir,
        element_size=args.element_size,
        dt=args.dt,
        t_end=args.t_end,
        rayleigh_override=override,
        abs_paths=args.abs_paths,
    )
    print(json.dumps(summary, indent=2, default=str))
    return 0


if __name__ == "__main__":
    sys.exit(main())
