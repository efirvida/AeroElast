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


def test_counts_and_row_positions_equal_the_native_mesh(native, decoupled):
    (mesh_n, data_n), (mesh_d, data_d) = native, decoupled
    assert mesh_d.node_count == mesh_n.node_count
    assert mesh_d.elements_count == mesh_n.elements_count
    zn, zd = _levels(np.asarray(data_n["nodes"])[:, 2]), _levels(np.asarray(data_d["nodes"])[:, 2])
    assert zd.shape == zn.shape
    assert np.abs(zd - zn).max() < 1e-6


def test_every_element_row_has_its_own_laminate(native, decoupled):
    _, data_n = native
    _, data_d = decoupled
    rows_native = _outer_levels(data_n)
    sets = _outer_sets(data_d)
    n_rows = rows_native.size - 1
    assert len(sets) == n_rows * 12
    # each set is exactly one native row: its extent is two consecutive row boundaries
    for name, (lo, hi) in sets.items():
        k = np.searchsorted(rows_native, lo - 1e-6)
        assert abs(rows_native[k] - lo) < 1e-6 and abs(rows_native[k + 1] - hi) < 1e-6, name
    # one laminate per (region, row)
    keys = {(OUTER_SET.match(n).group(3), round(lo, 6)) for n, (lo, _) in sets.items()}
    assert len(keys) == n_rows * 12
    layups = [s["elementSet"] for s in data_d["sections"] if OUTER_SET.match(s["elementSet"])]
    assert len(set(layups)) == n_rows * 12


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
    sets = {n: zr for n, zr in _outer_sets(data).items() if n.endswith("_HP_SPAR")}
    z_target = 0.90 * length
    name, (lo, hi) = next((n, zr) for n, zr in sets.items() if zr[0] <= z_target < zr[1])
    layup = next(s["layup"] for s in data["sections"] if s["elementSet"] == name)
    got = sum(t for material, t, _ in layup if material == "CarbonUD")
    z = np.linspace(lo, hi, 2001)
    want = np.mean(np.interp(z / length, cap["thickness"]["grid"], cap["thickness"]["values"]))
    assert hi - lo < 0.3
    assert got == pytest.approx(want, rel=0.005)


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
