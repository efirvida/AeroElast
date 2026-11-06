"""Layer drop-off steps (transverse closure segments) in the structured gmsh mesher: synthetic cylinder, no blade, no PETSc."""
import numpy as np
import pytest

gmsh = pytest.importorskip("gmsh")
shapely = pytest.importorskip("shapely")
gs = pytest.importorskip("aeroelast.models.blade.numad.mesh_gen.gmsh_skin")
base = pytest.importorskip("aeroelast.models.blade.numad.mesh_gen.gmsh_skin")

H, L, PER = 0.25, 3.0, 2.0
OPTS = {"patch_rows": 1, "extra_options": {"Mesh.Smoothing": 10.0, "Mesh.RecombineMinimumQuality": 0.7},
        "optimize": ("Relocate2D",), "optimize_iter": 5, "optimize_scope": "patches"}


def cylinder(eta, nd):
    r = PER / (2 * np.pi)
    phi = -2 * np.pi * np.asarray(nd, float)
    return np.column_stack([r * np.cos(phi), r * np.sin(phi), np.asarray(eta, float) * L])


def source(mod, lines, drops=()):
    return mod.SkinSource(lines, L, lambda e: np.full_like(np.asarray(e, float), PER), list(drops), cylinder)


def line(mod, name, nd, lo=0.0, hi=1.0):
    return mod.Line(name, np.array([0.0, 1.0]), np.array([nd, nd]), lo, hi, False)


def build(mod, lines, drops, **kw):
    return mod.generate_skin_mesh(source(mod, lines, drops), H, enforce_gate=False, structured=True, **{**OPTS, **kw})


def closure_report(res, eta, nd_a, nd_b, tol=1e-6):
    """(covered, straddling elements): chart nodes on the segment chain by mesh edges from nd_a to nd_b, none crossed."""
    from shapely.geometry import LineString, Polygon

    y, xa, xb = eta * L, nd_a * PER, nd_b * PER
    ch, el = res.chart, res.elements
    on = np.nonzero((np.abs(ch[:, 1] - y) < tol) & (ch[:, 0] > xa - tol) & (ch[:, 0] < xb + tol))[0]
    on = on[np.argsort(ch[on, 0], kind="stable")]
    edges = {tuple(sorted(e)) for e in map(tuple, gs.unique_edges(el))}
    chain = (len(on) >= 2 and abs(ch[on[0], 0] - xa) < tol and abs(ch[on[-1], 0] - xb) < tol
             and all(tuple(sorted((int(a), int(b)))) in edges for a, b in zip(on[:-1], on[1:], strict=False)))
    seg = LineString([(xa, y), (xb, y)])
    cross = 0
    for e in el:
        ids = e[e >= 0]
        yy = ch[ids, 1]
        if yy.min() < y - tol and yy.max() > y + tol and Polygon(ch[ids]).intersection(seg).length > tol:
            cross += 1
    return chain, cross


def check_mesh(res):
    assert res.metrics["jacobian"].min() > 0
    merged = np.unique(res.master)
    assert len(np.unique(np.round(res.xyz[merged], 9), axis=0)) == len(merged)  # no duplicated nodes after the seam identification
    slave = res.master != np.arange(len(res.master))
    assert slave.sum() > 0
    assert np.allclose(res.xyz[slave], res.xyz[res.master[slave]])  # the ring is closed
    area = 0.0
    for e in res.elements:
        ids = e[e >= 0]
        p = res.chart[ids]
        area += abs(0.5 * np.sum(p[:, 0] * np.roll(p[:, 1], -1) - np.roll(p[:, 0], -1) * p[:, 1]))
    assert area == pytest.approx(PER * L, rel=1e-9)  # no hole, no overlap


CASES = {
    "start_at_eta0": ([("a", 0.3), ("b", 0.7)], {}, [(0.0, "right")]),
    "end_at_eta1": ([("a", 0.3), ("b", 0.7)], {}, [(1.0, "left")]),
    "interior_start_and_end": ([("a", 0.3), ("b", 0.7), ("c", 0.1), ("d", 0.9)], {"a": (0.4, 0.7), "b": (0.4, 0.7)},
                               [(0.4, "right"), (0.7, "left")]),
    "line_between_the_ends": ([("a", 0.3), ("e", 0.5), ("b", 0.7)], {"a": (0.4, 1.0), "b": (0.4, 1.0)}, [(0.4, "right")]),
}


def make(mod, case):
    spec, span, drops = CASES[case]
    lines = [line(mod, n, nd, *span.get(n, (0.0, 1.0))) for n, nd in spec]
    return lines, [(eta, "a", "b", side) for eta, side in drops]


@pytest.mark.parametrize("case", list(CASES))
def test_drop_off_step_is_a_chain_of_mesh_edges(case):
    lines, drops = make(gs, case)
    res = build(gs, lines, drops)
    assert len(res.cleaned.closures) == len(drops)
    for eta, *_ in drops:
        chain, cross = closure_report(res, eta, 0.3, 0.7)
        assert chain, (case, eta)
        assert cross == 0, (case, eta)
    check_mesh(res)
    rep = res.info["closures"]
    assert len(rep) == len(drops) and all(r["covered"] for r in rep)


@pytest.mark.parametrize("case", ["start_at_eta0", "interior_start_and_end"])
def test_drop_off_mesh_is_deterministic(case):
    lines, drops = make(gs, case)
    a, b = build(gs, lines, drops), build(gs, lines, drops)
    assert np.array_equal(a.elements, b.elements) and np.array_equal(a.chart, b.chart)


def test_materials_differ_across_the_step():
    # the centroid of every element lies on one side of the step: the layer exists above eta = 0.4 only
    lines, drops = make(gs, "interior_start_and_end")
    res = build(gs, lines, drops[:1])
    ch, el = res.chart, res.elements
    y = 0.4 * L
    for e in el:
        yy = ch[e[e >= 0], 1]
        inside = ch[e[e >= 0], 0].mean() > 0.3 * PER and ch[e[e >= 0], 0].mean() < 0.7 * PER
        if inside:
            assert yy.max() <= y + 1e-6 or yy.min() >= y - 1e-6


@pytest.mark.parametrize("spec", [[("a", 0.3), ("b", 0.7)], [("a", 0.25), ("b", 0.5), ("c", 0.8)]])
def test_source_without_drops_is_byte_identical_to_gmsh_skin(spec):
    a = build(gs, [line(gs, n, nd) for n, nd in spec], [])
    b = build(base, [line(base, n, nd) for n, nd in spec], [])
    assert a.chart.tobytes() == b.chart.tobytes() and a.elements.tobytes() == b.elements.tobytes()
    assert a.xyz.tobytes() == b.xyz.tobytes() and a.master.tobytes() == b.master.tobytes()
    assert a.eta.tobytes() == b.eta.tobytes() and a.nd.tobytes() == b.nd.tobytes()
