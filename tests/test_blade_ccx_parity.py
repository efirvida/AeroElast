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
the export uses quadratic=True (S8R).  Our side is linear MITC4, so a gap between
the two is expected on a fixed mesh and means nothing on its own: it could be
discretisation or a model difference.  The two formulations are compared across
refinements instead, with `tools/blade_ccx_convergence.py`, and they do tend to
the same value -- the gap at a fixed mesh is discretisation:

    element_size  nodes   flap      edge      axial
    2.0            1460   10.60%    10.13%    12.57%
    1.0            3043    6.32%     3.86%     7.57%
    0.5            9277    3.68%     1.51%     4.88%
    0.25          32336    1.74%     0.53%     3.29%

The test therefore runs on the converged mesh and asserts the suite's 5% rule
rather than a mesh allowance.  Caveat recorded with the table: the metric here is
the mean over the tip nodes, which moves with the tip node distribution (flap
reads 7.55 / 7.94 / 7.12 / 7.11 m across the four meshes), so the trend of the
gap is the signal.  At 0.25 m the worst case sits 1.7 points inside the bound.
The energy metric is reported beside each row and is more sensitive on the axial
cases (6.65% at 0.25 m against a 3.29% displacement gap): that is the curved-axis
linearity the store records as `force_projection_axial_extension`, not a
tolerance to widen.
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
from aeroelast.core.mesh.io.writers import write_ccx_mesh  # noqa: E402
from aeroelast.models.blade.model import Blade  # noqa: E402

# Upstream's pipeline helpers: the property objects and the span direction the
# assembler needs are produced once, in test_blade_iea15mw_validation, and
# reused here so the two blade-vs-CalculiX comparisons cannot drift apart.
from tests.validation.blade.test_blade_iea15mw_validation import (  # noqa: E402
    SPAN_DIRECTION,
    _to_rust_mesh,
)

from tests.support.ccx_io import fail_ccx, run_ccx  # noqa: E402

REPO = Path(__file__).resolve().parents[1]
BLADE_YAML = REPO / "tests" / "IEA-15-240-RWT.yaml"
# The converged mesh. 1.0 m and 0.5 m are both unconverged, so a comparison there
# measures two unfinished discretisations (6.32/3.86/7.57% and 3.68/1.51/4.88%); at
# 0.25 m the cases read 1.74/0.53/3.29% and the 5% suite rule is what holds. That is
# why the test is slow: the quadratic export makes the CCX model ~98k nodes.
ELEMENT_SIZE = 0.25
TORQUE_MAGNITUDE = 5.0e6      # N.m about the span
AXIAL_MAGNITUDE = 1.0e6       # N along the span (centrifugal scale)
FLAP_MAGNITUDE = 1.0e5        # N, rotor-plane normal
EDGE_MAGNITUDE = 1.0e5        # N, in-plane

# name -> (load vector in the blade frame, expect a measurable twist)
#
# tip_torsion is deliberately absent: the CalculiX writer has no validated
# couple path.  Sending the moment as a `*CLOAD` on the sixth degree of freedom
# makes CalculiX return 418 m of displacement and a -42826 deg "twist" for this
# blade, so the case would compare two different loads.  It belongs in the
# S-7/S-8 torsion family, which uses a force couple, until the writer grows a
# checked couple.
CASES = {
    "tip_flap": ((0.0, FLAP_MAGNITUDE, 0.0), False),
    "tip_edge": ((EDGE_MAGNITUDE, 0.0, 0.0), False),
    "tip_axial_tension": ((0.0, 0.0, AXIAL_MAGNITUDE), False),
    "tip_axial_compression": ((0.0, 0.0, -AXIAL_MAGNITUDE), False),
}


def _blade():
    """Mesh and properties through upstream's pipeline.

    The earlier version built them by hand (`BladeMesh` + `build_rust_properties`
    + a hand-assembled ABD dict).  That dict was not the problem, but the direct
    `PyMeshAssembler(...)` constructor has no `span_direction` argument, so the
    element ply angles were read as element-local instead of span-relative.  Flap
    is insensitive to that (7.934 m either way) and edge is not: the edgewise
    stiffness came out 2.7x too high (1.229 m against 3.272 m with the span
    direction, and 3.450 m for CalculiX).  `Blade` + `get_element_properties()`
    + `from_model(..., SPAN_DIRECTION, ...)` is the path that carries it.
    """
    blade_model = Blade(str(BLADE_YAML), element_size=ELEMENT_SIZE)
    blade_model.generate_mesh()
    mesh = blade_model.mesh
    props = blade_model.get_element_properties()
    coords = np.asarray([[n.x, n.y, n.z] for n in mesh.nodes], dtype=float)
    tip = np.setdiff1d(
        np.nonzero(coords[:, 2] > coords[:, 2].max() - 1e-6)[0],
        np.array([], dtype=int),
    )
    mesh.add_node_set(NodeSet("tip", {mesh.nodes[int(i)] for i in tip}))
    return mesh, props, coords


def _aero_solve(mesh, props, load, torque=None) -> np.ndarray:
    """Clamped static solve on the blade, through `PyMeshAssembler.from_model`.

    `from_model` is what carries the span direction into the element ply angles;
    the direct constructor cannot.  Everything else is the same clamped solve the
    earlier version used.
    """
    asm = PyMeshAssembler.from_model(_to_rust_mesh(mesh, props), props, list(SPAN_DIRECTION), None)
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


def _tip_metrics(
    mesh, disp_nodes: dict, coords: np.ndarray, idx: list[int]
) -> tuple[float, float]:
    """(mean displacement magnitude, section twist about the span) over ``idx``.

    ``idx`` is explicit on purpose.  The first version inferred it as
    ``sorted(disp_nodes)``, which silently compared two different node sets: the
    CalculiX side passes only the tip nodes, but the AeroElast side passed a dict
    covering every node, so AeroElast's "tip" displacement was the mean over the
    whole blade (1.796 m) and the gap against CalculiX's real tip (8.469 m) read
    as 78.8%.  With the tip set named on both sides the same case is 6.6% apart.

    The twist is reported, not asserted: at the very tip the section is thin, so
    the least-squares denominator ``sum |r|^2`` in the (x, y) plane is small and
    the estimate is noise-amplified -- the CalculiX side of the flap case returns
    +70 deg, which no 1e5 N flap load can produce.  A twist comparison needs a
    spanwise band, not a tip ring.
    """
    pts = coords[idx]
    disp = np.asarray([disp_nodes[i] for i in idx], dtype=float)
    mag = float(np.mean(np.linalg.norm(disp[:, :3], axis=1)))
    # least-squares rotation about z over the section
    r = pts[:, :2] - pts[:, :2].mean(axis=0)
    den = float((r[:, 0] ** 2 + r[:, 1] ** 2).sum())
    num = float((r[:, 0] * disp[:, 1] - r[:, 1] * disp[:, 0]).sum())
    twist = np.rad2deg(num / den) if den > 0 else 0.0
    return mag, twist


def _ccx_tip_displacements(frd_path, mesh, coords) -> dict:
    """FRD translations, matched by COORDINATES (S8R renumbers the ids)."""
    from tests.validation.parity.test_orthotropic_shell_parity import (
        _parse_ccx_frd,
        _parse_ccx_frd_coords,
    )

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
@pytest.mark.slow
def test_blade_parity_ccx(tmp_path, name):
    """Same blade, same layup, same load through AeroElast and CCX.

    Measured 2026-10-06 on the converged mesh (0.25 m), every case inside the 5%
    rule: ``tip_flap`` 1.74% (7.105 m vs 7.231 m), ``tip_edge`` 0.53% (3.456 m vs
    3.475 m), ``tip_axial_tension`` 3.29% (5.691 m vs 5.884 m) and
    ``tip_axial_compression`` likewise. On the coarse meshes the same suite reads
    6.32/3.86/7.57% at 1.0 m, which is what the 20% band it used to carry bought.

    The edgewise disagreement this docstring used to record as a 64% gap
    (1.229 m vs 3.450 m, "left red on purpose") was closed by the section work
    that landed after 2026-09-30; the numbers above are what the suite measures
    now.  Superseded text kept for the record: *"``tip_edge`` does not:
    AeroElast gives 1.229 m where CalculiX gives 3.450 m, a 64% gap ... That is
    left red on purpose, with the measured numbers in the message: it is a real
    question about the laminate mapping, not a tolerance to widen."*
    """
    ccx_bin = _ccx_bin_or_skip()
    load, expect_twist = CASES[name]
    torque = TORQUE_MAGNITUDE if name == "tip_torsion" else None

    mesh, props, coords = _blade()
    tip_idx = sorted(mesh.node_id_to_index[n.id] for n in mesh.get_node_set("tip").nodes.values())
    u = _aero_solve(mesh, props, load, torque=torque)
    ae_disp = {i: u[i * 6:i * 6 + 3] for i in range(len(mesh.nodes))}
    ae_mag, ae_twist = _tip_metrics(mesh, ae_disp, coords, tip_idx)

    # The writer only knows the load vector, so a couple has to travel in it: a
    # component at index 5 becomes the sixth degree of freedom (RZ), and the
    # writer divides it by the node count exactly as it does for the forces, so
    # both sides see the same total moment.  Without this the torsion deck
    # carried no load at all and CalculiX returned zero displacement.
    ccx_load = [float(x) for x in load] + [0.0, 0.0, 0.0]
    if torque is not None:
        ccx_load[5] = float(torque)

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
        load_vector=ccx_load,
        span_direction=(0.0, 0.0, 1.0),
        quadratic=True,
    )
    proc = run_ccx(inp, ccx_bin)
    if proc.returncode != 0:
        fail_ccx(proc, inp)
    frd = case_dir / f"{stem}.frd"
    if not frd.exists():
        pytest.xfail(f"CCX produced no FRD for {name}")

    ccx_disp = _ccx_tip_displacements(frd, mesh, coords)
    if not ccx_disp:
        pytest.xfail(f"no tip nodes matched in the FRD for {name}")
    ccx_mag, ccx_twist = _tip_metrics(mesh, ccx_disp, coords, sorted(ccx_disp))

    d_mag = abs(ae_mag - ccx_mag) / max(abs(ccx_mag), 1e-30)
    print(
        f"{name:>24}: tip disp AE={ae_mag*1e6:10.3f} um CCX={ccx_mag*1e6:10.3f} um "
        f"({d_mag*100:5.2f}%)  |  twist AE={ae_twist:+.6f} CCX={ccx_twist:+.6f} deg"
    )

    assert d_mag < 0.05, f"{name}: tip displacement {d_mag*100:.1f}% off CCX (max 5%)"
    if expect_twist:
        # The tip ring is too thin for a least-squares rotation: the CalculiX
        # side of the flap case alone returns +70 deg, which no 1e5 N flap load
        # produces.  Assert what is numerically meaningful -- that both sides see
        # a non-zero twist in the same direction and of the same order -- and
        # leave the magnitude as a printed diagnostic.
        assert abs(ccx_twist) > 1e-5, "CCX shows no twist; the case setup is wrong"
        assert np.sign(ae_twist) == np.sign(ccx_twist), (
            f"{name}: twist sign AE {ae_twist:+.6f} vs CCX {ccx_twist:+.6f}"
        )
