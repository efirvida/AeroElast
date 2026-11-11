"""Shape optimizer for the thin-strip collapse (synthetic symmetric bands with a known answer)."""

from pathlib import Path

import numpy as np
import pytest

from aeroelast.models.blade.numad.mesh_gen.collapse_opt import (
    collapse_and_optimize_raw_mesh, collapse_capped, element_quality, mesh_stats, optimize_collapsed)
from aeroelast.models.blade.numad.mesh_gen.windio2_collapse import collapse_raw_mesh, collapse_thin_strips

ES = 0.25
DX = 0.24
RATIO = 0.15


def band(widths, dx=DX):
    """Three rows of quads: lower row, a band of half-widths +-w/2 around y=0, upper row (all ES high)."""
    nodes, ids = [], {}
    for i, w in enumerate(widths):
        for j, y in enumerate([-w / 2 - ES, -w / 2, w / 2, w / 2 + ES]):
            ids[i, j] = len(nodes)
            nodes.append((i * dx, y, 0.0))
    elems = [[ids[i, j], ids[i + 1, j], ids[i + 1, j + 1], ids[i, j + 1]] for i in range(len(widths) - 1) for j in range(3)]
    return np.array(nodes), np.array(elems), ids


def closing(w=0.12, n_open=4, n_thin=6):
    """A band of width w that ends abruptly: collapses to one triangle pair at the base ``w``."""
    nodes, elems, ids = band([w] * n_open + [0.001] * n_thin)
    cn, ce, src, rep = collapse_capped(nodes, elems, max_edge_len=RATIO * ES)
    return nodes, elems, cn, ce, rep


def tris(nodes, elems):
    sel = (elems < 0).any(axis=1)
    return elems[sel][:, :3], np.nonzero(sel)[0]


def min_angles(nodes, tri):
    out = []
    for t in tri:
        p = nodes[t]
        a = []
        for j in range(3):
            u, v = p[j - 1] - p[j], p[(j + 1) % 3] - p[j]
            a.append(np.degrees(np.arccos(np.clip(u @ v / np.linalg.norm(u) / np.linalg.norm(v), -1, 1))))
        out.append(min(a))
    return np.array(out)


def not_inverted(before, after, elems):
    for e in elems:
        e = e[e >= 0]
        n0 = np.cross(before[e][1] - before[e][0], before[e][2] - before[e][1])
        n1 = np.cross(after[e][1] - after[e][0], after[e][2] - after[e][1])
        if n1 @ n0 <= 0:
            return False
        for j in range(len(e)):
            cr = np.cross(after[e[j]] - after[e[j - 1]], after[e[(j + 1) % len(e)]] - after[e[j]])
            if cr @ n0 <= 0:
                return False
    return True


def test_equilateral_apex_distance_is_what_the_slide_reaches():
    nodes, elems, cn, ce, rep = closing(w=0.12)
    tri, _ = tris(cn, ce)
    assert len(tri) == 1                                         # the band's last quad
    before = min_angles(cn, tri)
    assert before.max() < 31                                     # (w, l, l) = (0.12, .24, .24): 29 deg
    new, orep = optimize_collapsed(cn, ce, max_slide=0.2)
    after = min_angles(new, tri)
    assert after.min() >= 53.0                                   # stated bound (60 = equilateral); measured 55.3
    apex = tri[0][np.argmax(np.abs(new[tri[0]][:, 0] - cn[tri[0]][:, 0]))]
    base_x = nodes[0][0] + 3 * DX                               # last open column
    assert new[apex, 0] - base_x == pytest.approx(np.sqrt(3) / 2 * 0.12, abs=0.02)   # 0.87 w from the base (0.104)
    assert not_inverted(cn, new, ce)
    assert orep["inverted"] == 0 and orep["nodes_moved"] >= 1
    # neighbours stay presentable: every quad within the default floor, AR <= 2.5
    st = mesh_stats(new, ce)
    assert st["max_ar"] < 2.5 and st["quad_max_angle_gt120"] <= 2
    # slide is along the merged line (y = 0): the apex does not leave it
    assert np.allclose(new[apex, 1], cn[apex, 1], atol=1e-6)


# Narrow bands cannot be made equilateral without breaking the neighbours: the two near sides become ~w
# while the quad next to them keeps a 0.25 side (AR ~ 0.25 / w).  Bounds are the measured values minus ~2 deg.
@pytest.mark.parametrize("w, bound", [(0.05, 30.0), (0.08, 40.0), (0.12, 53.0), (0.16, 58.0), (0.20, 59.0)])
def test_closing_triangle_improves_as_far_as_the_neighbours_allow(w, bound):
    nodes, elems, cn, ce, rep = closing(w=w)
    tri, _ = tris(cn, ce)
    before = min_angles(cn, tri).min()
    new, orep = optimize_collapsed(cn, ce, max_slide=0.2)
    after = min_angles(new, tri).min()
    assert after >= bound and after > before + 10
    assert not_inverted(cn, new, ce)


def test_fixed_nodes_never_move():
    nodes, elems, cn, ce, rep = closing()
    everything = set(np.unique(ce[ce >= 0]).tolist())
    new, orep = optimize_collapsed(cn, ce, fixed=everything, max_slide=0.2)
    assert np.array_equal(new, cn) and orep["nodes_moved"] == 0
    tri, _ = tris(cn, ce)
    apex = [int(i) for i in tri[0] if abs(cn[i, 0] - cn[tri[0], 0].min()) > 0.2][0]
    new, _ = optimize_collapsed(cn, ce, fixed={apex}, max_slide=0.2)
    assert np.array_equal(new[apex], cn[apex])


def test_slide_budget_shift_cap_and_quad_floor_are_respected():
    nodes, elems, cn, ce, rep = closing()
    tri, _ = tris(cn, ce)
    new, orep = optimize_collapsed(cn, ce, max_slide=0.05)
    assert np.linalg.norm(new - cn, axis=1).max() <= 0.05 + 1e-9
    # cumulative shift against the original nodes (map from the collapse)
    new2, orep2 = optimize_collapsed(cn, ce, nodes0=nodes, node_map=rep["node_map"], max_slide=0.2, max_shift=0.06)
    assert orep2["max_cumulative_shift_m"] <= max(0.06, rep["max_cumulative_shift_m"]) + 1e-9
    # a high floor keeps quads from degrading below min(floor, their own quality before)
    new3, _ = optimize_collapsed(cn, ce, max_slide=0.2, quad_floor=0.95)
    quad = (ce >= 0).all(axis=1)
    qb = element_quality(cn[ce[quad]]); qa = element_quality(new3[ce[quad]])
    assert (qa >= np.minimum(0.95, qb) - 1e-9).all()


def test_clean_mesh_is_untouched_and_the_optimizer_is_idempotent():
    nodes, elems, _ = band([0.25] * 6)
    new, rep = optimize_collapsed(nodes, elems)
    assert np.array_equal(new, nodes) and rep["nodes_moved"] == 0
    _, _, cn, ce, _ = closing()
    once, _ = optimize_collapsed(cn, ce, max_slide=0.2)
    twice, rep2 = optimize_collapsed(once, ce, max_slide=0.2)
    assert not_inverted(once, twice, ce) and rep2["quality_min_after"] >= rep2["quality_min_before"] - 1e-12
    assert np.linalg.norm(twice - once, axis=1).max() < 0.02          # converged: a second pass only polishes


def test_collapse_capped_equals_the_original_when_uncapped_and_bounds_drift_when_capped():
    w = [0.20, 0.16, 0.12, 0.08, 0.05, 0.03, 0.02, 0.01, 0.004, 0.001, 0.0005]
    nodes, elems, _ = band(w)
    ref = collapse_thin_strips(nodes, elems, max_edge_len=RATIO * ES)
    mine = collapse_capped(nodes, elems, max_edge_len=RATIO * ES)
    assert np.array_equal(ref[0], mine[0]) and np.array_equal(ref[1], mine[1]) and np.array_equal(ref[2], mine[2])
    assert ref[3]["removed_area_m2"] == pytest.approx(mine[3]["removed_area_m2"], rel=1e-12)
    capped = collapse_capped(nodes, elems, max_edge_len=RATIO * ES, max_shift=0.004)
    assert capped[3]["max_cumulative_shift_m"] <= 0.004 + 1e-12
    assert capped[3]["refused_by_cap"] >= 1
    assert (capped[3]["merged"] <= mine[3]["merged"])


def test_raw_mesh_wrapper_on_the_official_patch():
    from aeroelast.models.blade.numad.mesh_gen.windio2_mesh import generate_mesh
    from aeroelast.models.blade.numad.objects.blade import Blade

    yaml = Path(__file__).resolve().parents[2] / "IEA15MW" / "reference" / "upstream" / "IEA-15-240-RWT.yaml"
    blade = Blade()
    blade.read_yaml(str(yaml))
    raw = generate_mesh(blade.windio2_topology, ES, (0.85, 0.87), min_strip_m=1e-6)
    plain, prep = collapse_raw_mesh(raw, ES, RATIO)
    same, srep = collapse_and_optimize_raw_mesh(raw, ES, RATIO, max_shift=np.inf, optimize=False)
    assert np.array_equal(plain.elements, same.elements) and np.array_equal(plain.nodes, same.nodes)
    opt, orep = collapse_and_optimize_raw_mesh(raw, ES, RATIO, max_shift=0.05)
    n = len(opt.elements)
    assert len(opt.cells) == len(opt.surfaces) == len(opt.element_parameters) == len(opt.bases) == n
    assert orep["optimizer"]["inverted"] == 0
    assert orep["max_cumulative_shift_m"] <= 0.05 + 1e-9
    assert all(0 <= i < len(opt.nodes) for chain in opt.attachments.values() for i in chain)
    for k, chain in opt.attachments.items():                       # web attachment nodes never move
        for a, i in enumerate(chain):
            assert np.array_equal(np.asarray(opt.nodes, float)[i], np.asarray(raw.nodes, float)[raw.attachments[k][a]])
    assert mesh_stats(np.asarray(opt.nodes, float), opt.elements)["ar_gt10"] <= mesh_stats(np.asarray(plain.nodes, float), plain.elements)["ar_gt10"]
