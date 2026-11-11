"""Chord-scaled sizing of the gmsh skin mesher (odd/tasks/gmsh-chord-scaled-mesh.md, T3, S1-S3).

``chord_divisions=N`` replaces the element size: chordwise size c(eta)/N, spanwise size AR*c(eta)/N, at least
``min_band_elements`` elements across every material band.  Rules 1-4 and 6 run on a tapered synthetic ring (seconds); rule 5
(the ``element_size`` path is byte-for-byte unchanged) runs the real IEA-15 blade and is slow.
"""
import sys
from pathlib import Path

import numpy as np
import pytest

gmsh = pytest.importorskip("gmsh")
gs = pytest.importorskip("aeroelast.models.blade.numad.mesh_gen.gmsh_skin")

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "support"))
import v2_blade  # noqa: E402

L, D0, TAPER = 6.0, 2.0, 0.6  # span, root diameter, relative diameter loss at the tip


def diameter(eta):
    return D0 * (1.0 - TAPER * np.asarray(eta, float))


def ring(eta, nd):
    """Circle of diameter d(eta); clockwise nd -> outward normal for positive (eta, nd)."""
    eta = np.asarray(eta, float)
    phi = -2 * np.pi * np.asarray(nd, float)
    r = diameter(eta) / 2
    return np.column_stack([r * np.cos(phi), r * np.sin(phi), eta * L])


def source(lines=()):
    return gs.SkinSource(list(lines), L, lambda e: np.pi * diameter(e), [], ring, chord=diameter)


def straight(name, nd):
    return gs.Line(name, np.array([0.0, 1.0]), np.array([nd, nd]), 0.0, 1.0, False)


def mesh(n, ar=1.0, lines=(), **kw):
    return gs.generate_skin_mesh(source(lines), None, chord_divisions=n, span_aspect=ar, structured=True, patch_rows=1,
                                 enforce_gate=False, **kw)


def around(res, eta_lo, eta_hi):
    """Elements of the first row interval whose rows lie inside [eta_lo, eta_hi], and that interval's chart height."""
    rows = np.asarray(res.info["rows"])
    eta_rows = rows / L
    k = int(np.argmax(eta_rows >= eta_lo))
    y0, y1 = rows[k], rows[k + 1]
    assert eta_rows[k + 1] <= eta_hi
    ids = np.where(res.elements >= 0, res.elements, res.elements[:, :1])
    yc = res.chart[ids, 1].mean(1)
    return int(np.sum((yc > y0) & (yc < y1))), y1 - y0, 0.5 * (y0 + y1) / L


# --- 1. the chord helper -------------------------------------------------------------------------------------
class _Topology:
    def __init__(self, chord=None):
        shape = {} if chord is None else {"chord": {"grid": chord[0], "values": chord[1]}}
        self._definition = type("D", (), {"windio2_outer_shape": shape})()


class _OutlineLoft:
    """Loft stand-in whose station outlines are ellipses a x b (a circle when a == b); the TE sits at nd = 0 / 1."""

    def __init__(self, semi_axes):
        self.rows = np.linspace(0.0, 1.0, len(semi_axes))
        self.semi = semi_axes

    def station(self, i, nd):
        a, b = self.semi[i]
        phi = 2 * np.pi * np.asarray(nd, float)
        return np.column_stack([-a * np.cos(phi), b * np.sin(phi), np.zeros_like(phi)])


def test_chord_helper_uses_the_yaml_chord_when_the_topology_carries_it():
    c = gs.chord_profile(_Topology(([0.0, 0.5, 1.0], [5.2, 3.0, 0.5])), _OutlineLoft([(1, 1)] * 3))
    assert c(0.0) == pytest.approx(5.2) and c(1.0) == pytest.approx(0.5)
    assert 3.0 < float(c(0.25)) < 5.2


def test_chord_helper_falls_back_to_the_le_te_distance_and_is_the_diameter_at_the_circular_root():
    loft = _OutlineLoft([(1.5, 1.5), (1.0, 0.3), (0.5, 0.05)])  # circle of diameter 3 at the root, then thinner airfoil-like
    c = gs.chord_profile(_Topology(), loft)
    assert c(0.0) == pytest.approx(3.0, rel=1e-3)  # the root chord IS the diameter
    assert c(0.5) == pytest.approx(2.0, rel=1e-3)  # LE-TE distance of the outline
    assert c(1.0) == pytest.approx(1.0, rel=1e-3)


# --- 2. about 2N elements around, finer in absolute size at the tip ----------------------------------------------
def test_elements_around_a_section_do_not_depend_on_the_chord_and_the_tip_is_finer():
    n = 8
    res = mesh(n)
    root, dz_r, e_r = around(res, 0.10, 0.45)
    tip, dz_t, e_t = around(res, 0.55, 0.98)
    for count, e in ((root, e_r), (tip, e_t)):
        expected = np.pi * diameter(e) / (diameter(e) / n)  # perimeter / (c / N) = pi N for a ring
        assert count == pytest.approx(expected, rel=0.25)
    assert abs(root - tip) <= 0.25 * root  # roughly independent of the chord
    assert dz_t < dz_r  # absolute size: the tip is finer


# --- 3. spanwise spacing AR * c / N ----------------------------------------------------------------------------
def test_spanwise_spacing_follows_aspect_times_chord_over_n():
    n = 8
    for ar in (1.0, 2.0):
        rows = np.asarray(mesh(n, ar).info["rows"])
        eta = 0.5 * (rows[:-1] + rows[1:]) / L
        ratio = np.diff(rows) / (ar * diameter(eta) / n)
        assert np.median(ratio) == pytest.approx(1.0, rel=0.2), ar
        assert ratio.max() < 1.6, ar


# --- 4. the min-band floor ---------------------------------------------------------------------------------------
@pytest.mark.parametrize("min_band", [1, 2, 3])
def test_every_band_has_at_least_min_band_elements_across_at_small_n(min_band):
    n, gap = 4, 0.2  # chord 2 -> chordwise size 0.5 (merge threshold 0.125): the 0.2 m band would get one element without the floor
    lines = [straight("a", 0.5), straight("b", 0.5 + gap / (np.pi * D0))]
    res = mesh(n, lines=lines, min_band_elements=min_band, snap=0.4, min_band_fraction=None)  # None: the floor on every band
    x0, x1 = 0.5 * np.pi * D0, 0.5 * np.pi * D0 + gap  # the cylinder part of the chart: the ring is P(eta) wide, nd = x / P
    ids = np.where(res.elements >= 0, res.elements, res.elements[:, :1])
    nd = res.chart[ids, 0].mean(1) / (np.pi * diameter(res.chart[ids, 1].mean(1) / L))
    yc = res.chart[ids, 1].mean(1)
    inside = (nd > 0.5) & (nd < 0.5 + gap / (np.pi * diameter(yc / L)))
    per_row = {}
    for y in np.round(yc[inside], 6):
        per_row[y] = per_row.get(y, 0) + 1
    assert per_row and min(per_row.values()) >= min_band
    assert res.info["clean"]["n_bands_merged"] == 0  # the band is kept, not merged away
    assert res.gate["passed"], res.gate  # the floor must not buy elements the gate rejects (rows over the band are shortened)
    assert x0 < x1


# --- 6. the cache key ----------------------------------------------------------------------------------------------
def test_cache_key_differs_with_n_ar_and_min_band_and_is_stable_without_them():
    base = {"mesher": "gmsh"}
    k0 = v2_blade.cache_key(0.25, **base)
    assert v2_blade.cache_key(0.25, **base) == k0
    ks = {
        v2_blade.cache_key(0.25, chord_divisions=10, **base),
        v2_blade.cache_key(0.25, chord_divisions=20, **base),
        v2_blade.cache_key(0.25, chord_divisions=10, span_aspect=2.0, **base),
        v2_blade.cache_key(0.25, chord_divisions=10, min_band_elements=3, **base),
    }
    assert len(ks) == 4 and k0 not in ks
    assert v2_blade.cache_key(0.25, chord_divisions=10, **base) == v2_blade.cache_key(0.25, chord_divisions=10, span_aspect=1.0,
                                                                                    min_band_elements=2, **base)
    assert "gmsh_blade.py" in Path(v2_blade.__file__).read_text()  # source_digest("gmsh") hashes the blade wrapper too


# --- 7. graded strip-count transitions (T4c) -----------------------------------------------------------------------
def moving(name, nd0, nd1):
    return gs.Line(name, np.array([0.0, 1.0]), np.array([nd0, nd1]), 0.0, 1.0, False)


def test_a_strip_count_change_is_a_quad_transition_not_a_triangle_column():
    """The strip left of ``a`` widens from 4.5 to 7.5 elements: its count must change, and a change of one element
    (odd boundary parity) cannot be meshed with quads only: it left a column of triangles and 45-degree quads."""
    res = mesh(8, lines=[moving("a", 0.18, 0.30)], snap=0.4)
    steps = res.info["count_steps"]
    assert steps, "the strip must change its element count along the span"
    assert all(d % 2 == 0 for _i, _k, d in steps), steps  # an even change is the only all-quad transition
    assert int((res.elements[:, 3] < 0).sum()) == 0
    assert int((res.metrics["min_angle"] < 45.0).sum()) == 0
    assert res.gate["passed"], res.gate


def test_count_steps_never_reverse_along_the_span():
    """A monotone strip width gives a monotone count (no n -> n+1 -> n flip from rounding)."""
    res = mesh(8, lines=[moving("a", 0.18, 0.30)], snap=0.4)
    by = {}
    for i, _k, d in res.info["count_steps"]:
        by.setdefault(i, []).append(d)
    assert by and all(len({np.sign(d) for d in ds}) == 1 for ds in by.values()), by


def test_an_even_step_never_lands_outside_the_band_and_reverses_one_row_later():
    """The strip left of ``a`` narrows from 5.2 to 3.5 elements (N=8: width / (c/N) = 8 pi nd).  From n = 5 the nearest even count
    (3) is outside the band (ratio 1.27) the first time it beats staying at 5 (ratio 0.77), and below the even-step threshold the
    next rounding took the strip back to 4: 5 -> 3 -> 4, a reversal made of two transition patches in adjacent rows (the tip
    closure of IEA-15 N=20 chained them into one patch whose recombined quad had a 165.9 degree angle)."""
    res = mesh(8, lines=[moving("a", 0.207, 0.139)], snap=0.4)
    by = {}
    for i, _k, d in res.info["count_steps"]:
        by.setdefault(i, []).append(d)
    assert by, "the strip must change its element count along the span"
    assert all(len({np.sign(d) for d in ds}) == 1 for ds in by.values()), by
    assert res.gate["passed"], res.gate


# --- 8. row shortening over a band narrower than min_band * c/N (T4d) -----------------------------------------------
def _band_row_steps(res, gap_nd, min_band, size_c_ratio):
    """(row step, cap) per row interval where the band is narrower than min_band * c/N; cap = max_strip_ar * width / min_band."""
    rows = np.asarray(res.info["rows"])
    eta = rows / L
    width = gap_nd * np.pi * diameter(eta)  # the band's metric width at every row (straight lines: constant nd gap)
    out = []
    for k in range(len(rows) - 1):
        g = float(min(width[k], width[k + 1]))
        if g < min_band * size_c_ratio * float(diameter(eta[k : k + 2]).min()):
            out.append((rows[k + 1] - rows[k], 5.0 * g / min_band))
    return out


def test_row_shortening_caps_the_row_step_over_a_band_narrower_than_the_floor():
    """A 0.15 m band at N=4 is narrower than 2 * c/N (1.0 m at the root, 0.6 m at the tip): the nominal row step AR * c/N (0.5 .. 0.3 m)
    exceeds the cap 5 * 0.15 / 2 = 0.375 m near the root, so the rows there must be shorter than nominal.  Without the cap the
    two elements across the band have aspect ratio above the cap's design value (the gate rejected N=10 on IEA-15 without it)."""
    n, gap = 4, 0.15
    gap_nd = gap / (np.pi * D0)
    lines = [straight("a", 0.5), straight("b", 0.5 + gap_nd)]
    on = mesh(n, lines=lines, min_band_elements=2, snap=0.4, min_band_fraction=None)  # None: the band keeps the floor
    off = mesh(n, lines=lines, min_band_elements=2, snap=0.4, min_band_fraction=None, row_shortening=False)
    assert on.info["clean"]["n_bands_merged"] == 0 and off.info["clean"]["n_bands_merged"] == 0  # the band survives in both
    s_on = _band_row_steps(on, gap_nd, 2, 1.0 / n)
    s_off = _band_row_steps(off, gap_nd, 2, 1.0 / n)
    assert s_on and all(step <= cap * (1 + 1e-6) for step, cap in s_on), s_on
    assert any(step > cap * (1 + 1e-6) for step, cap in s_off), "the cap must bind: without it a nominal row step exceeds it"
    assert len(on.info["rows"]) > len(off.info["rows"])  # shorter rows over the band = more rows
    assert on.gate["passed"], on.gate


# --- 9. relative min-band floor (T4e) --------------------------------------------------------------------------------
def _elements_across(res, nd_lo, nd_hi):
    """Elements per row interval between the straight lines nd_lo and nd_hi (the ring is P(eta) wide, nd = x / P)."""
    ids = np.where(res.elements >= 0, res.elements, res.elements[:, :1])
    yc = res.chart[ids, 1].mean(1)
    nd = res.chart[ids, 0].mean(1) / (np.pi * diameter(yc / L))
    per_row = {}
    for y in np.round(yc[(nd > nd_lo) & (nd < nd_hi)], 6):
        per_row[y] = per_row.get(y, 0) + 1
    return per_row


def test_relative_floor_keeps_two_elements_only_in_a_band_above_the_fraction():
    """N=4: width / (c/N) = pi * N * gap_nd, constant along the span.  Band A (ratio 0.38) and band B (ratio 1.19) both get one
    element from rounding and two from the floor; with ``min_band_fraction=1.0`` only B (>= 1.0) keeps the floor."""
    n, ga, gb = 4, 0.03, 0.095
    lines = [straight("a", 0.2), straight("b", 0.2 + ga), straight("c", 0.5), straight("d", 0.5 + gb)]
    kw = {"lines": lines, "min_band_elements": 2, "snap": 0.4}
    full = mesh(n, min_band_fraction=None, **kw)
    rel = mesh(n, min_band_fraction=1.0, **kw)
    a_full, b_full = _elements_across(full, 0.2, 0.2 + ga), _elements_across(full, 0.5, 0.5 + gb)
    a_rel, b_rel = _elements_across(rel, 0.2, 0.2 + ga), _elements_across(rel, 0.5, 0.5 + gb)
    assert min(a_full.values()) >= 2 and min(b_full.values()) >= 2  # explicit None: the floor on both bands
    assert max(a_rel.values()) < 2  # below the threshold: the count rounding gives (one element)
    assert min(a_rel.values()) >= 1 and min(b_rel.values()) >= 2  # the band above the threshold keeps the floor
    for res in (full, rel):
        assert res.info["clean"]["n_bands_merged"] == 0 and res.gate["passed"], res.gate


def test_default_min_band_fraction_is_1_0_and_a_narrow_band_gets_fewer_than_min_band_elements():
    """S3 floor rule (L18): omitting the argument means 1.0, so a band narrower than 1.0 * c/N (ratio 0.38 at N=4) gets its
    rounded count, fewer than ``min_band_elements``; a band above it (1.19) keeps the floor; explicit None is the old floor of 2."""
    n, ga, gb = 4, 0.03, 0.095
    lines = [straight("a", 0.2), straight("b", 0.2 + ga), straight("c", 0.5), straight("d", 0.5 + gb)]
    kw = {"lines": lines, "min_band_elements": 2, "snap": 0.4}
    omitted, explicit, old = mesh(n, **kw), mesh(n, min_band_fraction=1.0, **kw), mesh(n, min_band_fraction=None, **kw)
    assert np.array_equal(omitted.elements, explicit.elements) and np.array_equal(omitted.xyz, explicit.xyz)
    a, b = _elements_across(omitted, 0.2, 0.2 + ga), _elements_across(omitted, 0.5, 0.5 + gb)
    assert max(a.values()) < 2 and min(a.values()) >= 1 and min(b.values()) >= 2
    assert min(_elements_across(old, 0.2, 0.2 + ga).values()) >= 2  # None: floor 2 whatever the width


def test_min_band_fraction_enters_every_chord_scaled_cache_key_and_none_cannot_collide():
    base = {"mesher": "gmsh", "chord_divisions": 10}
    k0 = v2_blade.cache_key(0.25, **base)
    assert v2_blade.cache_key(0.25, min_band_fraction=1.0, **base) == k0  # the default is 1.0
    keys = {k0, v2_blade.cache_key(0.25, min_band_fraction=None, **base), v2_blade.cache_key(0.25, min_band_fraction=0.5, **base)}
    assert len(keys) == 3
    # element_size meshes: the fraction is not part of the key
    es = {"mesher": "gmsh"}
    assert v2_blade.cache_key(0.25, min_band_fraction=None, **es) == v2_blade.cache_key(0.25, min_band_fraction=0.5, **es)


# --- 10. patch repair in the 3D metric (T4f) -----------------------------------------------------------------------
def _signed_area(chart, e):
    p = chart[[v for v in e if v >= 0]]
    return 0.5 * float(np.sum(p[:, 0] * np.roll(p[:, 1], -1) - np.roll(p[:, 0], -1) * p[:, 1]))


def test_a_patch_quad_with_a_near_straight_corner_is_cut_into_two_triangles():
    """The chart quad is a unit square, the same quad on the blade is a parallelogram with 164 / 16 degree corners (the sheared
    near-tip loft): above the 150 degree margin it becomes two triangles of the same orientation covering the same chart area."""
    chart = np.array([[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]])
    xyz = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.2, 0.35, 0.0], [1.2, 0.35, 0.0]])
    quad = np.array([[0, 1, 2, 3]])
    el, flag, rep = gs.repair_patch_elements(xyz, chart, quad, [True])
    assert rep["cut_quads"] == 1 and el.shape == (2, 4) and (el[:, 3] < 0).all() and flag.all()
    assert abs(sum(_signed_area(chart, e) for e in el) - 1.0) < 1e-12 and all(_signed_area(chart, e) > 0 for e in el)
    # a quad of the structured blocks (not a patch element) is never touched
    el2, _, rep2 = gs.repair_patch_elements(xyz, chart, quad, [False])
    assert rep2["cut_quads"] == 0 and np.array_equal(el2, quad)


def test_patch_triangle_pairs_are_flipped_only_when_it_raises_the_smallest_angle():
    chart = xyz2 = np.array([[0.0, 0.0], [1.0, 0.5], [2.0, 0.0], [1.0, -0.5]])
    xyz = np.column_stack([xyz2, np.zeros(4)])
    pair = np.array([[0, 2, 1, -1], [0, 3, 2, -1]])  # share the long diagonal 0-2: smallest angle 26.6 degrees
    el, _, rep = gs.repair_patch_elements(xyz, chart, pair, [True, True])
    assert rep["flipped_edges"] == 1 and len(el) == 2
    assert all(1 in e and 3 in e for e in el[:, :3]), el  # the short diagonal 1-3 now joins them
    assert all(_signed_area(chart, e) > 0 for e in el)
    same, _, rep0 = gs.repair_patch_elements(xyz, chart, pair, [True, False])  # one of the two is a block element: untouched
    assert rep0["flipped_edges"] == 0 and np.array_equal(same, pair)
    good = np.array([[0, 3, 1, -1], [3, 2, 1, -1]])  # already on the short diagonal
    assert gs.repair_patch_elements(xyz, chart, good, [True, True])[2]["flipped_edges"] == 0


# --- 5. the element_size path is unchanged -------------------------------------------------------------------------
@pytest.mark.slow
def test_element_size_path_keeps_the_iea15_counts_at_h_025():
    from aeroelast.models.blade.numad.mesh_gen.gmsh_skin import generate_mesh_gmsh
    from aeroelast.models.blade.numad.objects.blade import Blade

    blade = Blade()
    blade.read_yaml(str(v2_blade.DEFAULT_YAML))
    raw = generate_mesh_gmsh(blade.windio2_topology, 0.25, with_materials=True)
    assert (len(raw.nodes), len(raw.elements)) == (25012, 26018)
