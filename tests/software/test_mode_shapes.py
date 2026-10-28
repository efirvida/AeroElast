"""Mode-shape comparison helpers (odd/tasks/l3-closure.md T8, S6): BeamDyn .lin parser, ring-station projection, MAC.

* (i) the MAC of a complex mode is invariant to its sign and phase (``real_shape`` rotates to the dominant real phase);
* (ii) the station projection of a synthetic rigid rotation plus translation of a ring set recovers both exactly;
* (iii) the Hungarian pairing returns the inverse permutation for a permuted copy of the same shapes;
* (iv) the BeamDyn ``.lin`` parser returns the frequencies of ``logs/check/of_bd/bp20_0/modes.json`` (skipped without it);
* (v) the committed BeamDyn p20 parked eigenvectors (``tests/IEA15MW/of_bd_modal/bp20_0/eigvec-parked.npz``, T21) carry
  the eight families of the group 23 arbiter, and the saved T8 shell modes pair with them as the T21 prediction
  (skipped without ``logs/check/l3/t8/shell-parked.npz``).
"""

from __future__ import annotations

import json
import sys

import numpy as np
import pytest

from tests.support.paths import REPO_ROOT

sys.path.insert(0, str(REPO_ROOT / "tools"))
sys.path.insert(0, str(REPO_ROOT))

import l3_mode_shapes as ms  # noqa: E402
from tests.support.modal_pairing import mac_matrix, match_reference_by_shape  # noqa: E402


def test_mac_invariant_to_sign_and_phase_of_complex_mode():
    rng = np.random.default_rng(1)
    v = rng.normal(size=40) + 1j * rng.normal(size=40)
    w = rng.normal(size=40) + 1j * rng.normal(size=40)
    base = ms.real_shape(v)
    for alpha in (0.3, 1.1, 2.9, -2.0):
        for sign in (1.0, -1.0):
            rot = ms.real_shape(sign * np.exp(1j * alpha) * v)
            assert mac_matrix(base[None], rot[None])[0, 0] == pytest.approx(1.0, abs=1e-12)
    # the real part of the raw vector is NOT phase invariant: that is why the rotation is needed
    raw = mac_matrix(v.real[None], (np.exp(1j * 1.1) * v).real[None])[0, 0]
    assert raw < 0.99
    assert mac_matrix(base[None], ms.real_shape(w)[None])[0, 0] < 0.5
    assert np.linalg.norm(base) == pytest.approx(1.0)


def _ring_set(n_rings=7, n_nodes=24, length=6.0):
    """A set of elliptical rings at equal z, node ring ids, arc weights and a box mask (half the nodes)."""
    phi = np.linspace(0, 2 * np.pi, n_nodes, endpoint=False)
    z = np.linspace(0.0, length, n_rings)
    coords, ring_id = [], []
    for k, zk in enumerate(z):
        coords.append(np.c_[1.5 * np.cos(phi) + 0.1 * zk, 0.4 * np.sin(phi), np.full(n_nodes, zk)])
        ring_id += [k] * n_nodes
    coords = np.vstack(coords)
    w = np.tile(1.0 + 0.5 * np.cos(phi) ** 2, n_rings)  # non-uniform arc weights
    box = np.tile(np.arange(n_nodes) % 2 == 0, n_rings)
    return coords, np.array(ring_id), w, box, z


def test_station_projection_recovers_rigid_rotation_and_translation():
    coords, ring_id, w, box, z = _ring_set()
    t = np.c_[0.01 * z + 0.002, -0.02 * z + 0.5, np.zeros_like(z)]  # linear in z
    th = 0.003 * z + 0.01
    u = np.zeros((len(coords), 6))
    for k in range(len(z)):
        sel = ring_id == k
        d = coords[sel] - coords[sel].mean(axis=0)
        u[sel, :3] = t[k] + np.cross(np.array([0.0, 0.0, th[k]]), d)
    # stations in between rings: linear interpolation of a linear field is exact
    st = np.array([0.4, 1.9, 3.3, 5.5])
    for fit in ("all", "box"):
        ux, uy, theta = ms.station_fields(coords, u, ring_id, w, box, st, fit=fit)
        assert theta == pytest.approx(0.003 * st + 0.01, abs=1e-12)
        # arc-weighted mean of a rigid rotation: t + theta x (weighted centroid - plain centroid), so test against it
        exp = []
        for s in st:
            k0 = min(int(np.searchsorted(z, s)) - 1, len(z) - 2)
            f = (s - z[k0]) / (z[k0 + 1] - z[k0])
            vals = []
            for k in (k0, k0 + 1):
                sel = ring_id == k
                d = coords[sel] - coords[sel].mean(axis=0)
                c = (d * w[sel, None]).sum(axis=0) / w[sel].sum()
                tk = t[k] + np.cross([0, 0, th[k]], c)
                vals.append(tk[:2])
            exp.append((1 - f) * vals[0] + f * vals[1])
        exp = np.array(exp)
        assert ux == pytest.approx(exp[:, 0], abs=1e-12)
        assert uy == pytest.approx(exp[:, 1], abs=1e-12)


def test_station_projection_box_fit_with_too_few_nodes_is_nan_not_filled():
    coords, ring_id, w, box, z = _ring_set()
    box = np.zeros_like(box)
    u = np.zeros((len(coords), 6))
    _, _, theta = ms.station_fields(coords, u, ring_id, w, box, np.array([1.0]), fit="box")
    assert np.isnan(theta).all()


def test_pairing_returns_inverse_permutation_of_permuted_copy():
    rng = np.random.default_rng(3)
    shapes = np.array([ms.real_shape(rng.normal(size=30) + 1j * rng.normal(size=30)) for _ in range(8)])
    perm = rng.permutation(8)
    shell = shapes[perm]  # shell mode i is BeamDyn mode perm[i]
    pair, mac = ms.pair_shell_to_bd(shell, shapes)
    assert (shell[pair] == shapes).all()
    assert pair.tolist() == np.argsort(perm).tolist()
    assert mac[pair, np.arange(8)] == pytest.approx(np.ones(8))


def test_gll_node_heights_match_the_p60_beamdyn_nodes():
    # BeamDyn p60 nodes 53..61 initial z (logs/check/worktree-archive/l3-pass2/check/bd_static/LC1-p60-s0.0001-nodes)
    z = ms.bd_node_z(60)
    assert z[52:] == pytest.approx([111.7189, 112.9087, 113.9517, 114.8453, 115.587, 116.1748, 116.6071, 116.8827, 117.0], abs=1e-4)


_MODES = REPO_ROOT / "logs/check/of_bd/bp20_0/modes.json"


@pytest.mark.skipif(not _MODES.exists(), reason="logs/check/of_bd/bp20_0/modes.json (BeamDyn p20 linearization record) is absent")
def test_lin_parser_reproduces_recorded_parked_frequencies():
    pytest.importorskip("openfast_toolbox", reason="openfast_toolbox (fx_mbc3) not available")
    lin = REPO_ROOT / "logs/check/of_bd/bp20_0/case.1.lin"
    if not lin.exists():
        pytest.skip(f"{lin} is absent")
    rec = [m["hz"] for m in json.loads(_MODES.read_text())["modes"]]
    out = ms.bd_modes_from_lin([lin], n_modes=8)
    assert out["freqs"] == pytest.approx(rec[:8], rel=1e-6)
    assert out["trans"].shape == (8, 20, 3)
    assert out["rot"].shape == (8, 20, 3)
    assert out["z"].shape == (20,)


_BD_FIXTURE = REPO_ROOT / "tests/IEA15MW/of_bd_modal/bp20_0/eigvec-parked.npz"
_T8_SHELL = REPO_ROOT / "logs/check/l3/t8/shell-parked.npz"


def test_committed_beamdyn_eigenvectors_carry_the_eight_arbiter_families():
    z, chord, shapes = ms.bd_station_modes(np.load(_BD_FIXTURE))
    assert ms.family_labels(*shapes, chord) == ["1F", "1E", "2F", "2E", "3F", "1T", "3E", "4F"]
    assert z.shape == (20,)
    assert z[-1] == pytest.approx(117.0)


@pytest.mark.skipif(not _T8_SHELL.exists(), reason="logs/check/l3/t8/shell-parked.npz (T8 shell eigenvectors) is absent")
def test_t8_shell_modes_pair_with_beamdyn_as_the_t21_prediction():
    # odd/tasks/l3-closure.md "T21 prediction": the first 10 shell modes (the group 23 N_SEARCH), fit over all ring nodes.
    z, chord, bd = ms.bd_station_modes(np.load(_BD_FIXTURE))
    sh = np.load(_T8_SHELL)
    ux, uy, th, bad = ms.shell_station_modes(sh["coords"], sh["modes"][:10], sh["ring_id"], sh["ring_w"], None, z)
    assert bad == 0
    got = match_reference_by_shape(ms.vectors(ux, uy, th, chord), ms.vectors(*bd, chord), ["F", "E", "F", "E", "F", "T", "E", "F"])
    assert [m.index for m in got] == [0, 1, 2, 3, 4, 5, 6, 8]
    assert [round(m.mac_family, 3) for m in got[5:]] == [0.990, 0.973, 0.924]
