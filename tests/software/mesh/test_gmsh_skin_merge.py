"""Merge criterion of the 2D cleaning (gmsh_skin): a persistent narrow band is kept as a strip, wedges still close,
near-coincident lines whose cleaned values swap order between two rows no longer raise.  Synthetic cylinder, no blade."""
import numpy as np
import pytest

gmsh = pytest.importorskip("gmsh")
shapely = pytest.importorskip("shapely")
gm = pytest.importorskip("aeroelast.models.blade.numad.mesh_gen.gmsh_skin")

H, L, PER, TAU = 0.25, 3.0, 2.0, 0.25
OPTS = {"patch_rows": 1, "extra_options": {"Mesh.Smoothing": 10.0, "Mesh.RecombineMinimumQuality": 0.7},
        "optimize": ("Relocate2D",), "optimize_iter": 5, "optimize_scope": "patches"}


def cylinder(eta, nd):
    r = PER / (2 * np.pi)
    phi = -2 * np.pi * np.asarray(nd, float)
    return np.column_stack([r * np.cos(phi), r * np.sin(phi), np.asarray(eta, float) * L])


def source(mod, lines):
    return mod.SkinSource(lines, L, lambda e: np.full_like(np.asarray(e, float), PER), [], cylinder)


def const(mod, name, nd, lo=0.0, hi=1.0, fixed=False):
    return mod.Line(name, np.array([0.0, 1.0]), np.array([nd, nd]), lo, hi, fixed)


def linear(mod, name, nd0, nd1, fixed=False):
    return mod.Line(name, np.array([0.0, 1.0]), np.array([nd0, nd1]), 0.0, 1.0, fixed)


def build(mod, lines, **kw):
    return mod.generate_skin_mesh(source(mod, lines), H, enforce_gate=False, structured=True, **{**OPTS, **kw})


def covers(res, x, tol=1e-6):
    on = np.abs(res.chart[:, 0] - x) < tol
    tot = sum(abs(res.chart[a, 1] - res.chart[b, 1]) for a, b in gm.unique_edges(res.elements) if on[a] and on[b])
    return abs(tot - L) < 1e-6


# (i) a persistent narrow band (0.6 tau h wide along the whole span: strip AR h/w = 6.7 < gate 10) is kept
def test_persistent_narrow_band_is_kept_as_a_strip():
    w = 0.6 * TAU * H
    lines = [const(gm, "a", 0.4), const(gm, "b", 0.4 + w / PER), const(gm, "c", 0.8)]
    res = build(gm, lines)
    rep = res.info["clean"]
    assert rep["n_bands_merged"] == 0 and rep["n_bands_kept"] == 1
    assert covers(res, 0.4 * PER) and covers(res, 0.4 * PER + w)  # both lines are mesh edges: the strip exists
    assert res.metrics["jacobian"].min() > 0 and res.gate["passed"], res.gate
    assert res.metrics["aspect"].max() < gm.DEFAULT_LIMITS["max_ar"]
    assert res.metrics["aspect"].max() <= 5.0 + 1e-6  # max_strip_ar = 5: rows over the strip are at most 5 w = 0.19 m long
    assert res.metrics["aspect"].max() > 3.0  # and the strip is really there (nominal 5, 4.45 measured on the chord-shortened cylinder)
    assert len(res.info["rows"]) > L / H + 1  # the rows were refined over the strip: the cost of keeping it
    old = build(gm, [const(gm, "a", 0.4), const(gm, "b", 0.4 + w / PER), const(gm, "c", 0.8)], max_strip_ar=None)
    assert old.info["clean"]["n_bands_merged"] == 1  # what the unmodified criterion does


def test_band_below_the_aspect_ratio_cap_is_still_merged_and_cap_is_a_parameter():
    w = H / 12.0  # 5 w = 0.104 m of row would be < h / 2: keeping the strip costs more than a factor 2 in rows
    lines = [const(gm, "a", 0.4), const(gm, "b", 0.4 + w / PER)]
    assert gm.clean_source(source(gm, lines), H, TAU).report["n_bands_merged"] == 1
    assert gm.clean_source(source(gm, lines), H, TAU, max_strip_ar=14.0).report["n_bands_merged"] == 0
    assert gm.clean_source(source(gm, lines), H, TAU, max_strip_ar=None).report["n_bands_merged"] == 1


# (ii) wedges still close
def test_wedge_closing_to_zero_is_still_closed_at_tau_h():
    src = source(gm, [const(gm, "a", 0.5), linear(gm, "b", 0.5, 0.5 + 0.2 / PER)])
    rep = gm.clean_source(src, H, TAU).report
    assert rep["n_bands_merged"] == 1 and rep["n_bands_kept"] == 0
    assert rep["bands"][0]["s_to_m"] == pytest.approx(TAU * H / 0.2 * L, rel=1e-6)


def test_band_that_pinches_below_the_cap_is_merged_over_its_whole_thin_stretch():
    # gap 0.2 -> 0.01 -> 0.2 (V shape): thin where < tau h, its minimum is below h / cap
    b = gm.Line("b", np.array([0.0, 0.5, 1.0]), np.array([0.5 + 0.2 / PER, 0.5 + 0.01 / PER, 0.5 + 0.2 / PER]))
    rep = gm.clean_source(source(gm, [const(gm, "a", 0.5), b]), H, TAU).report
    assert rep["n_bands_merged"] == 1 and rep["n_bands_kept"] == 0
    assert rep["bands"][0]["max_gap_m"] <= TAU * H + 1e-9


# (iii) near-coincident lines whose cleaned values swap order between two rows
def swap_lines(mod):
    """c (above) and a (0.02 m below it) merge at the start; a leaves c and its closing step at s0 coincides with the start
    of b (b copies a's source value there): the step moves a up past b between two rows."""
    s0 = 0.5
    a = mod.Line("a", np.array([0.0, s0, 1.0]), np.array([0.5 - 0.02 / PER, 0.5 - 0.0625 / PER, 0.5 - 0.3 / PER]))
    b = mod.Line("b", np.array([s0, 1.0]), np.array([0.5 - 0.0625 / PER, 0.5 - 0.3 / PER]), s0, 1.0, False)
    return [a, b, const(mod, "c", 0.5, fixed=True)]


def test_cleaned_order_swap_between_rows_does_not_raise():
    res = build(gm, swap_lines(gm))
    assert res.metrics["jacobian"].min() > 0
    area = sum(abs(0.5 * np.sum(res.chart[e[e >= 0], 0] * np.roll(res.chart[e[e >= 0], 1], -1)
                                - np.roll(res.chart[e[e >= 0], 0], -1) * res.chart[e[e >= 0], 1])) for e in res.elements)
    assert area == pytest.approx(PER * L, rel=1e-9)  # no hole, no overlap


# (iv) differential identity: sources with no narrow band and no crossing give byte-identical meshes
IDENT = {
    "parallel": lambda m: [const(m, "a", 0.25), const(m, "b", 0.5), const(m, "c", 0.75)],
    "wedge": lambda m: [const(m, "a", 0.5), linear(m, "b", 0.5, 0.5 + 0.2 / PER)],
    "thin_merge": lambda m: [const(m, "a", 0.5), const(m, "b", 0.5 + 0.01 / PER)],
    "taper": lambda m: [const(m, "a", 0.25), linear(m, "b", 0.5, 0.5 + 1.0 / PER)],
}


@pytest.mark.parametrize("case", sorted(IDENT))
def test_no_narrow_band_no_crossing_is_byte_identical_to_the_unmodified_module(case):
    new = build(gm, IDENT[case](gm), snap=0.4)
    old = build(gm, IDENT[case](gm), snap=0.4, max_strip_ar=None)
    for f in ("xyz", "elements", "master", "chart", "eta", "nd"):
        assert getattr(new, f).tobytes() == getattr(old, f).tobytes(), f
