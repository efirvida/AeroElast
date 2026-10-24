"""The spanwise element rows come from ``element_size`` alone; the laminate follows each row.

Candidate B of the blade input audit (``odd/tasks/blade-input-audit.md``, T7e, user L12). A NuMAD
(region, bay) holds one bay-mean laminate and the mesher cuts ``ceil(bay / size)`` rows per bay.
Inserting stations to refine the laminate (``max_bay_length``) therefore also refined the mesh
(spanwise size 0.119 m everywhere, +47 % elements). The default ``BladeMesh`` now plans the native
rows, puts a station at every row boundary and meshes one row per bay: the discretisation is the
native one and every row has its own laminate.
"""

import re
from pathlib import Path

import numpy as np
import pytest
import yaml

from aeroelast.core.mesh.generators import BladeMesh
from aeroelast.models.blade.numad.mesh_gen.mesh_gen import SpanPlan

YAML = Path(__file__).parent / "IEA-15-240-RWT.yaml"
ES = 0.25
OUTER_SET = re.compile(r"^(\d{2})_(\d{2,})_((?:HP|LP)_\w+)$")
# the spar caps are cut at the yaml web arcs into pieces ``HP_SPAR_A``, ``HP_SPAR_B``, ...
SPAR_PIECE = re.compile(r"^(?:HP|LP)_SPAR(?:_[A-Z])?$")
SLIVER = 0.01  # m, a column narrower than this is a taper artefact, not a discretisation


def _build(**kwargs):
    gen = BladeMesh(str(YAML), element_size=ES, **kwargs)
    mesh = gen.generate(verbose=False)
    return mesh, gen.numad_mesh_data


@pytest.fixture(scope="module")
def native():
    return _build(laminate_per_row=False)


@pytest.fixture(scope="module")
def decoupled():
    return _build()


def _levels(z, tol=1e-6):
    z = np.sort(np.asarray(z, dtype=float))
    return np.array([g.mean() for g in np.split(z, np.nonzero(np.diff(z) > tol)[0] + 1)])


def _outer_levels(data):
    """z of the row boundaries of the outer shell."""
    nodes = np.asarray(data["nodes"])
    elements = np.asarray(data["elements"])
    labels = np.concatenate(
        [es["labels"] for es in data["sets"]["element"] if OUTER_SET.match(es["name"])]
    )
    ids = np.unique(elements[labels])
    return _levels(nodes[ids[ids >= 0], 2])


def _outer_sets(data):
    """{name: (zmin, zmax)} of the outer-shell element sets."""
    nodes = np.asarray(data["nodes"])
    elements = np.asarray(data["elements"])
    out = {}
    for es in data["sets"]["element"]:
        if OUTER_SET.match(es["name"]):
            ids = np.unique(elements[es["labels"]])
            z = nodes[ids[ids >= 0], 2]
            out[es["name"]] = (z.min(), z.max())
    return out


def _row_nodes(mesh, levels):
    """Per row level: the sorted unique node coordinates (x, y) on that z plane."""
    xyz = np.asarray(mesh.coords_array, dtype=float)
    rows = []
    for level in levels:
        pts = xyz[np.abs(xyz[:, 2] - level) < 1e-6]
        rows.append(pts)
    return rows


def _min_gap(pts):
    """Smallest distance between two distinct nodes of one row."""
    d = np.linalg.norm(pts[:, None, :2] - pts[None, :, :2], axis=-1)
    d[np.diag_indices_from(d)] = np.inf
    return d.min()


def test_counts_and_row_positions_equal_the_native_mesh(native, decoupled):
    """Same elements, same rows; the node counts differ only by sliver cap-piece columns.

    The merged T8w spar-cap split (pieces ``_A/_B/_C`` at the yaml web arcs) decides in the native
    mesh whether a piece exists per native bay and tapers its width linearly, so in the rows where
    the yaml webs coincide (root) or nearly do (tip) the native mesh carries an extra column only
    0-5 mm wide. The decoupled mesh places the exact web arcs at the inserted stations and drops
    the zero-width piece. Blade input audit T9: the decoupled mesh is the correct one, so the
    contract is equal elements, equal row levels, and a node deficit of exactly two per row
    (one sliver column) confined to rows whose native node spacing shows a sliver.
    """
    (mesh_n, data_n), (mesh_d, data_d) = native, decoupled
    assert mesh_d.elements_count == mesh_n.elements_count
    zn, zd = _levels(np.asarray(data_n["nodes"])[:, 2]), _levels(np.asarray(data_d["nodes"])[:, 2])
    assert zd.shape == zn.shape
    assert np.abs(zd - zn).max() < 1e-6

    assert mesh_d.node_count <= mesh_n.node_count
    rows_n, rows_d = _row_nodes(mesh_n, zn), _row_nodes(mesh_d, zn)
    assert sum(len(r) for r in rows_n) == mesh_n.node_count
    assert sum(len(r) for r in rows_d) == mesh_d.node_count
    diffs = np.array([len(rn) - len(rd) for rn, rd in zip(rows_n, rows_d, strict=True)])
    differing = np.nonzero(diffs)[0]
    assert differing.size < 0.05 * zn.size  # equal counts everywhere else
    explained = 0
    for i in differing:
        assert diffs[i] > 0 and diffs[i] % 2 == 0, (zn[i], diffs[i])
        # the native row has a sliver column, the decoupled row has none
        assert _min_gap(rows_n[i]) < SLIVER, (zn[i], _min_gap(rows_n[i]))
        assert _min_gap(rows_d[i]) >= 1e-3, (zn[i], _min_gap(rows_d[i]))
        explained += diffs[i]
    assert mesh_n.node_count - mesh_d.node_count == explained


def test_every_element_row_has_its_own_laminate(native, decoupled):
    _, data_n = native
    _, data_d = decoupled
    rows_native = _outer_levels(data_n)
    sets = _outer_sets(data_d)
    n_rows = rows_native.size - 1
    # each region or cap piece is split once per row, and each set is exactly one native row
    # (its extent is two consecutive row boundaries)
    for name, (lo, hi) in sets.items():
        k = np.searchsorted(rows_native, lo - 1e-6)
        assert abs(rows_native[k] - lo) < 1e-6 and abs(rows_native[k + 1] - hi) < 1e-6, name
    per_row = {}
    for name, (lo, _) in sets.items():
        per_row.setdefault(round(lo, 6), []).append(OUTER_SET.match(name).group(3))
    assert len(per_row) == n_rows
    for lo, regions in per_row.items():
        assert len(regions) == len(set(regions)), lo  # no region or piece twice in a row
        assert len(regions) >= 12, lo  # 12 regions plus any cap pieces
    # one laminate per (region-or-piece, row)
    keys = {(OUTER_SET.match(n).group(3), round(lo, 6)) for n, (lo, _) in sets.items()}
    assert len(keys) == len(sets)
    layups = [s["elementSet"] for s in data_d["sections"] if OUTER_SET.match(s["elementSet"])]
    assert len(layups) == len(set(layups)) == len(sets)


def test_tapering_layer_equals_the_yaml_mean_over_its_row(decoupled):
    """Spar cap (CarbonUD) at r/R 0.90 against the yaml thickness averaged over the same row."""
    _, data = decoupled
    blade = yaml.safe_load(YAML.read_text())["components"]["blade"]
    length = blade["outer_shape_bem"]["reference_axis"]["z"]["values"][-1]
    cap = next(
        layer
        for layer in blade["internal_structure_2d_fem"]["layers"]
        if layer["name"] == "Spar_Cap_SS"
    )
    z_target = 0.90 * length
    pieces = {
        n: zr
        for n, zr in _outer_sets(decoupled[1]).items()
        if SPAR_PIECE.match(OUTER_SET.match(n).group(3)) and zr[0] <= z_target < zr[1]
    }
    assert pieces
    for name, (lo, hi) in pieces.items():
        layup = next(s["layup"] for s in data["sections"] if s["elementSet"] == name)
        got = sum(t for material, t, _ in layup if material == "CarbonUD")
        z = np.linspace(lo, hi, 2001)
        want = np.mean(np.interp(z / length, cap["thickness"]["grid"], cap["thickness"]["values"]))
        assert hi - lo < 0.3
        assert got == pytest.approx(want, rel=0.005), name


def test_span_plan_boundaries_and_native_chordwise_counts():
    plan = SpanPlan([0.0, 3.0, 4.0], [3, 1])
    np.testing.assert_allclose(plan.row_boundaries, [1.0, 2.0])
    plan.chord_counts[(("shell", 0), 0)] = (2, 4)
    plan.chord_counts[(("shell", 1), 0)] = (4, 2)
    plan.chord_counts[(("shell", 0), 1)] = (3, 3)
    plan.bind([0.0, 1.0, 2.0, 3.0, 4.0])
    # the transition sits in the end row next to the coarser edge, as in the native quad3 region
    assert [plan.counts(("shell", 0), i) for i in range(3)] == [(2, 4), (4, 4), (4, 4)]
    assert [plan.counts(("shell", 1), i) for i in range(3)] == [(4, 4), (4, 4), (4, 2)]
    assert plan.counts(("shell", 0), 3) == (3, 3)
    assert plan.counts(("web", 0), 0) is None
    with pytest.raises(ValueError):
        plan.bind([0.0, 1.0, 3.0, 4.0])
