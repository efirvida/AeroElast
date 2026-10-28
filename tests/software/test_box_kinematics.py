"""Box kinematics (odd/tasks/l3-closure.md T15, S15 (a), S8): the twist / curvature rigid fit on the box nodes.

The box is the caps (element sets whose name contains ``SPAR``) and the webs (names ending ``_SW``), the C3 rule
of ``tools/l3_twist_integral.py box_region``. The all-node fit reads the in-plane distortion of the skin and TE
panels as rotation (core-torsion C3); the box fit does not. ``fit_nodes="all"`` stays the default and unchanged.

* (i) ``box_node_mask`` picks exactly the nodes of the SPAR / _SW sets;
* (ii) a rigid rotation plus an in-plane distortion on the non-box nodes only: the box fit returns the exact
  rotation, the all-node fit does not (``station_rigid_rotation`` and the S-7 ``section_twists_deg``);
* (iii) a station with fewer than 6 box nodes is NaN, counted, never filled with all nodes;
* (iv) the S-1 GJ and tangent-frame EI extractors take the option on a short tube (finite on the box, NaN counted
  when the box has too few nodes per station, default unchanged).
"""

from __future__ import annotations

import sys

import numpy as np
import pytest

pytest.importorskip("_aeroelast", reason="Rust backend not available")

from tests.support.paths import REPO_ROOT  # noqa: E402

sys.path.insert(0, str(REPO_ROOT / "tools"))

import diagnose_s1_gj_tube as tube  # noqa: E402
import run_s1_sectional as s1  # noqa: E402
import run_s7_torsion as s7  # noqa: E402

from aeroelast.core.mesh import ElementSet  # noqa: E402
from aeroelast.postprocess.sectional import (  # noqa: E402
    SectionalExtractor,
    box_node_mask,
    station_rigid_rotation,
)

LENGTH = 6.0


def _tube_with_box(web_x=None, both_caps=False):
    """Rectangular tube (es 0.25) plus a ``WEB_SW`` set (the x = +-2 walls) and a ``CAP_SPAR_A`` set (the top wall).

    ``web_x`` picks which side walls are in the web set (default: both); ``both_caps`` puts the bottom wall in the
    cap set too, so the box is the whole tube.
    """
    g_layer = np.full(int(round(LENGTH / tube.ES)), tube.G_UNIFORM)
    mesh, props, _zl = tube.build_tube(g_layer, "rect", tube.T_WALL, length=LENGTH)
    mesh.coords_array = np.asarray([n.coords for n in mesh.nodes], dtype=float)
    webs, caps = set(), set()
    for e in mesh.elements:
        xs = np.array([mesh.nodes[mesh.node_id_to_index[n]].coords[0] for n in e.node_ids])
        ys = np.array([mesh.nodes[mesh.node_id_to_index[n]].coords[1] for n in e.node_ids])
        if np.allclose(xs, xs[0]) and (web_x is None or abs(xs[0] - web_x) < 1e-9):
            webs.add(e)
        elif np.allclose(ys, ys[0]) and (ys[0] > 0 or both_caps):
            caps.add(e)
    mesh.add_element_set(ElementSet("WEB_SW", webs))
    mesh.add_element_set(ElementSet("CAP_SPAR_A", caps))
    return mesh, props, webs, caps


def _nodes_of(mesh, elements):
    idx = mesh.node_id_to_index
    return sorted({idx[n] for e in elements for n in e.node_ids})


def test_i_box_mask_selects_spar_and_sw_sets_only():
    mesh, _props, webs, caps = _tube_with_box()
    mask = box_node_mask(mesh)
    assert mask.dtype == bool and mask.shape == (len(mesh.nodes),)
    expected = sorted(set(_nodes_of(mesh, webs)) | set(_nodes_of(mesh, caps)))
    assert np.nonzero(mask)[0].tolist() == expected
    # the other (station) sets do not count: a node on the bottom wall only is outside
    bottom = [i for i, n in enumerate(mesh.nodes) if n.coords[1] < -0.9 and abs(n.coords[0]) < 1.0]
    assert bottom and not mask[bottom].any()


def test_i_box_mask_empty_when_no_box_sets():
    g_layer = np.full(4, tube.G_UNIFORM)
    mesh, _props, _ = tube.build_tube(g_layer, "rect", tube.T_WALL, length=1.0)
    assert not box_node_mask(mesh).any()


def _ring(n=24, a=2.0, b=1.0):
    phi = 2.0 * np.pi * np.arange(n) / n
    return np.column_stack([a * np.cos(phi), b * np.sin(phi), np.zeros(n)])


def test_ii_box_fit_exact_all_fit_biased_by_distortion():
    xyz = _ring()
    theta = np.array([0.0, 0.0, 1.0e-3])
    u = np.cross(theta, xyz - xyz.mean(axis=0)) + np.array([1e-4, -2e-4, 5e-5])
    box = np.zeros(len(xyz), bool)
    box[::2] = True  # 12 box nodes
    # in-plane distortion on the non-box nodes only (a TE panel breathing)
    rng = np.random.default_rng(1)
    u[~box, :2] += 5.0e-4 * rng.standard_normal(((~box).sum(), 2))
    nodes = np.arange(len(xyz))
    rot_box = station_rigid_rotation(xyz, nodes, u, fit_nodes="box", box_mask=box)
    rot_all = station_rigid_rotation(xyz, nodes, u, fit_nodes="all")
    np.testing.assert_allclose(rot_box, theta, atol=1e-12)
    assert abs(rot_all[2] - theta[2]) > 1e-5
    # default is "all"
    np.testing.assert_allclose(station_rigid_rotation(xyz, nodes, u), rot_all, atol=0.0)


def test_iii_fewer_than_six_box_nodes_is_nan_never_all_nodes():
    xyz = _ring()
    u = np.cross([0.0, 0.0, 1e-3], xyz)
    box = np.zeros(len(xyz), bool)
    box[:5] = True
    rot = station_rigid_rotation(xyz, np.arange(len(xyz)), u, fit_nodes="box", box_mask=box)
    assert rot.shape == (3,) and np.isnan(rot).all()
    box[5] = True  # six nodes: fitted
    rot = station_rigid_rotation(xyz, np.arange(len(xyz)), u, fit_nodes="box", box_mask=box)
    np.testing.assert_allclose(rot, [0.0, 0.0, 1e-3], atol=1e-12)


def test_unknown_fit_nodes_is_refused():
    xyz = _ring()
    with pytest.raises(ValueError, match="fit_nodes"):
        station_rigid_rotation(xyz, np.arange(len(xyz)), np.zeros_like(xyz), fit_nodes="web")
    mesh, props, *_ = _tube_with_box()
    with pytest.raises(ValueError, match="fit_nodes"):
        s1.ea_gj_from_static_response(mesh, props, fit_nodes="web")


def test_ii_s7_section_twists_box_exact_all_biased():
    mesh, _props, webs, caps = _tube_with_box()
    coords = mesh.coords_array
    theta = 2.0e-3
    u = np.zeros(len(coords) * 6)
    uu = u.reshape(-1, 6)
    uu[:, 0] = -theta * coords[:, 1]
    uu[:, 1] = theta * coords[:, 0]
    mask = box_node_mask(mesh)
    rng = np.random.default_rng(3)
    uu[~mask, :2] += 1.0e-3 * rng.standard_normal(((~mask).sum(), 2))
    z_a, th_all = s7.section_twists_deg(mesh, u)
    z_b, th_box = s7.section_twists_deg(mesh, u, fit_nodes="box")
    assert np.allclose(z_a, z_b)
    np.testing.assert_allclose(th_box, np.rad2deg(theta), rtol=1e-9)
    assert np.abs(th_all - np.rad2deg(theta)).max() > 1e-3
    # default call is the all-node one, unchanged
    np.testing.assert_array_equal(s7.section_twists_deg(mesh, u)[1], th_all)


def test_iii_s7_section_without_enough_box_nodes_is_nan():
    """Box sets only below z = 3 m: the slices above have no box node, so NaN (kept, not dropped, not all-node)."""
    mesh, _props, webs, caps = _tube_with_box()
    idx = mesh.node_id_to_index
    low = {e for e in webs | caps if max(mesh.nodes[idx[n]].coords[2] for n in e.node_ids) <= 3.0 + 1e-9}
    mesh.element_sets.pop("WEB_SW")
    mesh.element_sets.pop("CAP_SPAR_A")
    mesh.add_element_set(ElementSet("LOW_SW", low))
    coords = mesh.coords_array
    u = np.zeros(len(coords) * 6)
    uu = u.reshape(-1, 6)
    uu[:, 0], uu[:, 1] = -1e-3 * coords[:, 1], 1e-3 * coords[:, 0]
    z_all, th_all = s7.section_twists_deg(mesh, u)
    z_box, th_box = s7.section_twists_deg(mesh, u, fit_nodes="box")
    assert np.allclose(z_all, z_box) and len(th_box) == len(th_all)
    above = z_box > 3.0 + 0.2
    below = z_box < 3.0 - 0.2
    assert above.any() and below.any()
    assert np.isnan(th_box[above]).all()
    np.testing.assert_allclose(th_box[below], th_all[below], rtol=1e-9)


def test_iv_gj_and_ei_extractors_take_fit_nodes_on_a_tube():
    mesh, props, webs, caps = _tube_with_box(both_caps=True)
    assert box_node_mask(mesh).all()  # the box is the whole tube: box and all-node fits are the same fit
    stations = np.arange(0.5, LENGTH, 1.0)
    info: dict = {}
    z0, _ea0, gj_all, _ = s1.ea_gj_from_static_response(mesh, props, stations=stations)
    z1, _ea1, gj_box, _ = s1.ea_gj_from_static_response(mesh, props, stations=stations, fit_nodes="box", info=info)
    assert np.allclose(z0, z1)
    mid = (z0 > 0.2 * z0.max()) & (z0 < 0.85 * z0.max())
    assert np.isfinite(gj_box[mid]).all()
    assert info["n_box_nan"] == 0
    np.testing.assert_allclose(gj_box[mid], gj_all[mid], rtol=1e-9)

    ext = SectionalExtractor(mesh, props)
    _z, fl_all, ed_all = ext.ei_from_static_response(props, frame="tangent", stations=stations)
    info_e: dict = {}
    _z, fl_box, ed_box = ext.ei_from_static_response(
        props, frame="tangent", stations=stations, fit_nodes="box", info=info_e)
    assert info_e["n_box_nan"] == 0
    np.testing.assert_allclose(ed_box, ed_all, rtol=1e-9)
    np.testing.assert_allclose(fl_box, fl_all, rtol=1e-9)

    # the solves can be shared between readings: passing the recorded solutions back skips them, same numbers
    rec_e: dict = {}
    ext.ei_from_static_response(props, frame="tangent", stations=stations, info=rec_e)
    _z, fl_re, ed_re = ext.ei_from_static_response(
        props, frame="tangent", stations=stations, solutions=rec_e["solutions"])
    np.testing.assert_array_equal(ed_re, ed_all)
    np.testing.assert_array_equal(fl_re, fl_all)
    rec_g: dict = {}
    s1.ea_gj_from_static_response(mesh, props, stations=stations, info=rec_g)
    _z, _e, gj_re, _ = s1.ea_gj_from_static_response(mesh, props, stations=stations, solutions=rec_g["solutions"])
    np.testing.assert_array_equal(gj_re, gj_all)
    with pytest.raises(ValueError, match="frame"):
        ext.ei_from_static_response(props, frame="global", fit_nodes="box")


def test_iv_extractors_count_nan_when_box_too_small():
    # webs only on the +x side and a cap set emptied: a station holds 2 rings x 9 nodes of the web = 18 >= 6, so
    # shrink the box to a single ring-node column instead
    mesh, props, webs, caps = _tube_with_box(web_x=2.0)
    keep = ElementSet("TINY_SW", set(list(webs)[:1]))  # one element: 4 nodes < 6 at every station
    mesh.element_sets.pop("WEB_SW")
    mesh.element_sets.pop("CAP_SPAR_A")
    mesh.add_element_set(keep)
    stations = np.arange(0.5, LENGTH, 1.0)
    info: dict = {}
    _z, _ea, gj, _ = s1.ea_gj_from_static_response(mesh, props, stations=stations, fit_nodes="box", info=info)
    assert info["n_box_nan"] == len(stations)
    assert np.isnan(gj).all()


def test_v_distributed_loads_tributary_weights_and_exact_beam_side():
    import l3_box_energy as be

    mesh, _props, _webs, _caps = _tube_with_box()
    coords = mesh.coords_array
    root = be.root_indices(mesh)
    nodes, forces = be.distributed_nodal_loads(coords, be.mesh_element_nodes(mesh), root, 0, 1.0e6)
    assert forces[:, 0].sum() == pytest.approx(1.0e6, rel=1e-12)
    assert np.abs(forces[:, 1:]).max() == 0.0
    assert not set(nodes) & set(root)
    # uniform ring: the force per ring is constant except the half-weight tip ring
    z = coords[nodes, 2]
    per_ring = {round(zz, 6): forces[z == zz, 0].sum() for zz in np.unique(z)}
    interior = [v for k, v in per_ring.items() if k < LENGTH - 1e-6]
    assert np.ptp(interior) < 1e-9 * max(interior)
    assert per_ring[round(LENGTH, 6)] == pytest.approx(0.5 * interior[0], rel=1e-9)
    # beam side of a uniform load on a prismatic deck: sum of the applied nodal forces = q L^5 / 40 EI exactly (continuum limit)
    kp = np.array([[0.0, 0.0, zz, 0.0] for zz in np.linspace(0.0, LENGTH, 5)])
    axis = be.reference_axis(kp)
    k = np.diag([1e12, 1e12, 1e12, 4e9, 1e9, 1e9])
    zl = (np.arange(2000) + 0.5) * LENGTH / 2000
    xyz = np.column_stack([np.zeros(2000), np.zeros(2000), zl])
    f = np.column_stack([np.full(2000, 100.0 * LENGTH / 2000), np.zeros(2000), np.zeros(2000)])
    w = be.beam_energy(axis, np.array([0.0, 1.0]), np.stack([k, k]), xyz, f, distributed=True)
    q = 100.0
    # twist 0: shell X = -BD y, so the load bends about the BD x axis = K44 (edge), as the S-1 X load reads EI_edge
    assert w["W_decoupled"] == pytest.approx(q**2 * LENGTH**5 / 40.0 / 4e9, rel=1e-5)
