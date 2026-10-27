"""Unit tests of the S-2 interior-station helpers (synthetic data; no mesh, no solver)."""

from __future__ import annotations

import json

import numpy as np
import pytest

from tests.support.s2_station import (
    bd_station_nodes,
    component_index,
    pick_station_node,
    ring_mean_at_z,
)


def _rings(z_values, fn):
    """Ring-profile rows [z, ux, uy, uz, n] for a field fn(z) -> (ux, uy, uz)."""
    return [[float(z), *fn(z), 100] for z in z_values]


def test_linear_field_is_reproduced_exactly_between_rings():
    field = lambda z: (0.3 * z - 1.0, -2.0 * z + 5.0, 0.01 * z)  # noqa: E731
    rings = _rings([110.0, 110.4, 110.9, 111.5, 112.0], field)
    for z in (110.05, 110.4, 111.2, 111.99):
        got = ring_mean_at_z(rings, z)
        assert np.allclose(got, field(z), rtol=0, atol=1e-12)


def test_ring_order_does_not_matter_and_outside_profile_raises():
    field = lambda z: (z, 2 * z, 3 * z)  # noqa: E731
    rings = _rings([112.0, 110.0, 111.0], field)
    assert np.allclose(ring_mean_at_z(rings, 110.5), field(110.5))
    with pytest.raises(ValueError, match="outside the ring profile"):
        ring_mean_at_z(rings, 109.9)
    with pytest.raises(ValueError, match="outside the ring profile"):
        ring_mean_at_z(rings, 112.1)


def _gll(z_tip, n=9):
    """Synthetic GLL-like output nodes, last one at the tip, spacing growing inboard."""
    zs = [z_tip - d for d in (8.0, 6.0, 4.3, 3.0, 2.16, 1.41, 0.7, 0.2, 0.0)][-n:]
    return [{"node": 53 + k, "z": z} for k, z in enumerate(zs)]


def test_pick_is_inner_of_the_two_nodes_nearest_the_tip_at_or_beyond_min_chords():
    nodes = _gll(120.0)
    chord = 0.5  # 2 chords = 1.0 m: candidates at 1.41 (node 58) and 2.16 (node 57)
    picked = pick_station_node(nodes, z_tip=120.0, chord=chord)
    assert picked["node"] == 57
    assert picked["z"] == pytest.approx(120.0 - 2.16)


def test_pick_ignores_nodes_closer_than_min_chords_and_the_tip():
    nodes = _gll(120.0)
    picked = pick_station_node(nodes, z_tip=120.0, chord=0.8)  # 2 chords = 1.6 m: 2.16 and 3.0 qualify
    assert picked["node"] == 56
    with pytest.raises(ValueError, match="fewer than two"):
        pick_station_node(nodes, z_tip=120.0, chord=5.0)


def test_pick_counts_a_node_exactly_at_min_chords():
    nodes = [{"node": 1, "z": 10.0}, {"node": 2, "z": 9.0}, {"node": 3, "z": 8.0}, {"node": 4, "z": 12.0}]
    # tip distances 0, 2, 3, 4 chords: the two nearest at or beyond 2 are z 10 and 9; the inner one is z 9
    assert pick_station_node(nodes, z_tip=12.0, chord=1.0)["node"] == 2


def test_component_index_is_ux_for_lc3_and_uy_otherwise():
    assert [component_index(c) for c in ("LC1", "LC2", "LC3", "LC4")] == [1, 1, 0, 1]


def _write_run(root, case, scale, uy, ux):
    d = root / f"{case}-p60-s{scale:g}-nodes"
    d.mkdir(parents=True)
    nodes = [{"node": 56, "z": 117.0, "linear_shell_uy_m": uy, "linear_shell_ux_m": ux}]
    (d / "result.json").write_text(json.dumps({"nodes": nodes}))


def test_bd_station_nodes_reports_the_scale_spread(tmp_path):
    _write_run(tmp_path, "LC2", 1e-4, 60.0, 1.0)
    _write_run(tmp_path, "LC2", 2e-4, 60.03, 1.0)
    nodes, spread = bd_station_nodes(tmp_path, "LC2")
    assert nodes[0]["node"] == 56
    assert spread == pytest.approx(0.03 / 60.03, rel=1e-9)


def test_bd_station_nodes_missing_reference_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        bd_station_nodes(tmp_path, "LC1")
