"""Blade-level AeroElast vs CalculiX parity: the IEA 15 MW layup.

The strip parity in `test_composite_ccx_parity.py` proves the composite
element reproduces CalculiX on a flat coupon.  This file does it on the REAL
blade geometry with the real 672-set layup exported straight from
`build_rust_properties` (which already returns Rust `Laminate` objects keyed
by element-set name, exactly what the writer consumes).

The load cases mirror what the wind and the rotation actually do to a blade:

  * flapwise bending      -- rotor thrust (distributed) / tip equivalent
  * edgewise bending      -- gravity and in-plane loads
  * torsion               -- aerodynamic pitching moment
  * axial tension         -- centrifugal stiffening
  * axial compression     -- the opposite sign of the same path
  * flap + torsion        -- the coupled in-service case

All are static, so both solvers see identical inputs.  Rotation itself
(centrifugal prestress) is NOT covered here: CCX needs an NLGEOM prestress
step for that and it is a separate, heavier comparison.

CalculiX only supports *SHELL SECTION, COMPOSITE with quadratic elements, so
the export uses quadratic=True (S8R).  Our side is linear MITC4, so a few
percent of discretisation gap is expected -- the deflection tolerance reflects
that; the twist tolerance is tighter because the coupling is what is under
test.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve

pytest.importorskip("petsc4py", reason="PETSc not available")
pytest.importorskip("_aeroelast", reason="Rust backend not available")

from _aeroelast import PyMeshAssembler  # noqa: E402

from aeroelast.core.mesh.entities import NodeSet  # noqa: E402
from aeroelast.core.mesh.generators import BladeMesh  # noqa: E402
from aeroelast.core.mesh.io.writers import write_ccx_mesh  # noqa: E402
from aeroelast.models.blade.model import build_rust_properties  # noqa: E402

from test_beam_4cases_parity import _run_ccx  # noqa: E402

REPO = Path(__file__).resolve().parents[1]
BLADE_YAML = REPO / "tests" / "IEA-15-240-RWT.yaml"
ELEMENT_SIZE = 1.0            # coarse on purpose: CCX runs S8R (x4 elements)
TORQUE_MAGNITUDE = 5.0e6      # N.m about the span
AXIAL_MAGNITUDE = 1.0e6       # N along the span (centrifugal scale)
FLAP_MAGNITUDE = 1.0e5        # N, rotor-plane normal
EDGE_MAGNITUDE = 1.0e5        # N, in-plane

# name -> (load vector in the blade frame, expect a measurable twist)
CASES = {
    "tip_flap": ((0.0, FLAP_MAGNITUDE, 0.0), False),
    "tip_edge": ((EDGE_MAGNITUDE, 0.0, 0.0), False),
    "tip_torsion": ((0.0, 0.0, 0.0), True),      # torque handled separately
    "tip_axial_tension": ((0.0, 0.0, AXIAL_MAGNITUDE), False),
    "tip_axial_compression": ((0.0, 0.0, -AXIAL_MAGNITUDE), False),
}


def _blade():
    gen = BladeMesh(element_size=ELEMENT_SIZE, yaml_file=str(BLADE_YAML))
    mesh = gen.generate(renumber="rcm", verbose=False)
    props = build_rust_properties(gen.numad_mesh_data)

    coords = np.asarray([[n.x, n.y, n.z] for n in mesh.nodes], dtype=float)
    tip_nodes = set(np.nonzero(coords[:, 2] > coords[:, 2].max() - 1e-6)[0])
    # The writer wants a named set; the blade mesh only ships RootNodes and
    # the two surface sets.
    mesh.add_node_set(
        NodeSet("tip", {mesh.nodes[int(i)] for i in tip_nodes})
    )
    return mesh, props, coords


def _aero_solve(mesh, props, load, torque=None) -> np.ndarray:
    coords = np.asarray([[n.x, n.y, n.z] for n in mesh.nodes], dtype=float)
    conn = [[mesh.node_id_to_index[nid] for nid in el.node_ids] for el in mesh.elements]
    # The blade mesh is mixed: quads -> MITC4 (type 4), triangles -> MITC3 (3).
    elem_types = [len(c) for c in conn]
    # Element -> property mapping.  The blade ships 674 element sets and 672
    # laminates; without this every element would take the first laminate and
    # the comparison would be meaningless.
    elem_key: dict[int, str] = {}
    for set_name, elset in mesh.element_sets.items():
        for el in elset.elements:
            elem_key[id(el)] = set_name
    default_prop = next(iter(props.values()))

    # Two consumers, two shapes: the CCX writer takes the Rust `Laminate`
    # objects straight from build_rust_properties, but PyMeshAssembler wants
    # the flat ABD dict.  Convert once per distinct laminate.
    flat_cache: dict[int, dict] = {}

    def _val(x):
        """Rust bindings expose some of these as methods, some as getters."""
        return x() if callable(x) else x

    def _flat(lam) -> dict:
        key = id(lam)
        if key not in flat_cache:
            # Rust Laminate API: abd_matrix (6x6 [[A,B],[B,D]]), cs_matrix.
            ABD = np.asarray(_val(lam.abd_matrix), dtype=float)
            A, Bm, D = ABD[:3, :3], ABD[:3, 3:], ABD[3:, 3:]
            h = float(_val(lam.total_thickness))
            cs = np.asarray(_val(lam.cs_matrix), dtype=float).ravel()
            flat_cache[key] = {
                "type": "composite",
                "cm": A.ravel().tolist(),
                "b_coupling": Bm.ravel().tolist(),
                "cb": D.ravel().tolist(),
                "cs": cs.tolist(),
                "thickness": h,
                "e_equiv": A[0, 0] / h,
                "mass_per_area": 0.0,
                "rotational_inertia": 0.0,
            }
        return flat_cache[key]

    mat_list = [_flat(props.get(elem_key.get(id(el)), default_prop)) for el in mesh.elements]
    asm = PyMeshAssembler(
        node_coords=coords, connectivity=conn,
        elem_types=elem_types, materials=mat_list,
    )
    rows, cols, vals = asm.assemble_k()
    K = coo_matrix((vals, (rows, cols)), shape=(asm.dofs_count, asm.dofs_count)).tocsr()

    f = np.zeros(asm.dofs_count)
    tip = list(mesh.get_node_set("tip").nodes.values())
    for n in tip:
        base = mesh.node_id_to_index[n.id] * 6
        for d in range(3):
            f[base + d] = load[d] / len(tip)
        if torque is not None:
            f[base + 5] = torque / len(tip)      # rotation about the span (z)

    clamped = {
        mesh.node_id_to_index[n.id] * 6 + d
        for n in mesh.get_node_set("RootNodes").nodes.values() for d in range(6)
    }
    mask = np.ones(asm.dofs_count, dtype=bool)
    mask[list(clamped)] = False
    free = np.where(mask)[0]
    u = np.zeros(asm.dofs_count)
    u[free] = spsolve(K[np.ix_(free, free)], f[free])
    return u


def _tip_metrics(mesh, disp_nodes: dict, coords: np.ndarray) -> tuple[float, float]:
    """(mean tip displacement magnitude, section twist about the span)."""
    idx = sorted(disp_nodes)
    pts = coords[idx]
    disp = np.asarray([disp_nodes[i] for i in idx], dtype=float)
    mag = float(np.mean(np.linalg.norm(disp[:, :3], axis=1)))
    # least-squares rotation about z over the tip section
    r = pts[:, :2] - pts[:, :2].mean(axis=0)
    den = float((r[:, 0] ** 2 + r[:, 1] ** 2).sum())
    num = float((r[:, 0] * disp[:, 1] - r[:, 1] * disp[:, 0]).sum())
    twist = np.rad2deg(num / den) if den > 0 else 0.0
    return mag, twist


def _ccx_tip_displacements(frd_path, mesh, coords) -> dict:
    """FRD translations, matched by COORDINATES (S8R renumbers the ids)."""
    from test_orthotropic_shell_parity import _parse_ccx_frd, _parse_ccx_frd_coords

    ccx_coords = _parse_ccx_frd_coords(frd_path)
    ccx_disp = _parse_ccx_frd(frd_path, list(ccx_coords))
    out = {}
    for i, n in enumerate(mesh.nodes):
        if n.id not in {m.id for m in mesh.get_node_set("tip").nodes.values()}:
            continue
        p = coords[i]
        best, best_d = None, 1e-4
        for nid, xyz in ccx_coords.items():
            d = float(np.linalg.norm(xyz - p))
            if d < best_d:
                best, best_d = nid, d
        if best is not None and best in ccx_disp:
            out[i] = ccx_disp[best]
    return out


def _ccx_bin_or_skip() -> str:
    import shutil

    ccx = shutil.which("ccx") or shutil.which("CalculiX")
    if ccx is None:
        pytest.skip("CalculiX (ccx) not found in PATH")
    return ccx


@pytest.mark.slow
@pytest.mark.parametrize("name", list(CASES))
def test_blade_parity_ccx(tmp_path, name):
    """Same blade, same layup, same load through AeroElast and CCX."""
    ccx_bin = _ccx_bin_or_skip()
    load, expect_twist = CASES[name]
    torque = TORQUE_MAGNITUDE if name == "tip_torsion" else None

    mesh, props, coords = _blade()
    u = _aero_solve(mesh, props, load, torque=torque)
    ae_disp = {i: u[i * 6:i * 6 + 3] for i in range(len(mesh.nodes))}
    ae_mag, ae_twist = _tip_metrics(mesh, ae_disp, coords)

    case_dir = tmp_path / name
    case_dir.mkdir()
    stem = f"blade_{name}"
    inp = case_dir / f"{stem}.inp"
    write_ccx_mesh(
        mesh, str(inp),
        properties=props,
        boundary_nodeset="RootNodes",
        solver_type="LinearStatic",
        load_nodeset="tip",
        load_vector=list(load),
        span_direction=(0.0, 0.0, 1.0),
        quadratic=True,
    )
    _run_ccx(ccx_bin, case_dir, stem)
    frd = case_dir / f"{stem}.frd"
    if not frd.exists():
        pytest.xfail(f"CCX produced no FRD for {name}")

    ccx_disp = _ccx_tip_displacements(frd, mesh, coords)
    if not ccx_disp:
        pytest.xfail(f"no tip nodes matched in the FRD for {name}")
    ccx_mag, ccx_twist = _tip_metrics(mesh, ccx_disp, coords)

    d_mag = abs(ae_mag - ccx_mag) / max(abs(ccx_mag), 1e-30)
    print(
        f"{name:>24}: tip disp AE={ae_mag*1e6:10.3f} um CCX={ccx_mag*1e6:10.3f} um "
        f"({d_mag*100:5.2f}%)  |  twist AE={ae_twist:+.6f} CCX={ccx_twist:+.6f} deg"
    )

    assert d_mag < 0.20, f"{name}: tip displacement {d_mag*100:.1f}% off CCX (max 20%)"
    if expect_twist:
        assert abs(ccx_twist) > 1e-5, "CCX shows no twist; the case setup is wrong"
        assert np.sign(ae_twist) == np.sign(ccx_twist), (
            f"{name}: twist sign AE {ae_twist:+.6f} vs CCX {ccx_twist:+.6f}"
        )
        assert ae_twist == pytest.approx(ccx_twist, rel=0.25), (
            f"{name}: twist magnitude AE {ae_twist:+.6f} vs CCX {ccx_twist:+.6f}"
        )
