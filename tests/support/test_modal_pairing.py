"""Unit tests of the shape-based mode pairing and cross-mesh tracking (synthetic data, no solver)."""

from __future__ import annotations

import numpy as np
import pytest

from tests.support.modal_pairing import (
    ModalSet,
    followed_gaps,
    mac_matrix,
    mode_signature,
    pair_modes_by_mac,
    parse_frd_mode_translations,
    track_modes,
)


def _shape(z: np.ndarray, axis: int, order: int) -> np.ndarray:
    """Cantilever-like translation shape along ``axis`` (n_nodes, 3)."""
    out = np.zeros((len(z), 3))
    s = z / z.max()
    out[:, axis] = np.sin((order - 0.5) * np.pi * s)
    return out


def _modal_set(n_nodes: int, specs: list[tuple[int, int, float]], f_ae_scale: float) -> ModalSet:
    """Modes given as (axis, order, f_ccx) in the CCX frequency order; AE is the same set, shuffled."""
    z = np.linspace(0.0, 100.0, n_nodes)
    ccx_t = np.array([_shape(z, ax, od) for ax, od, _ in specs])
    f_ccx = np.array([f for _, _, f in specs])
    # AeroElast lists the same shapes in another order (by its own frequencies).
    perm = np.arange(len(specs))[::-1]
    ae_t = ccx_t[perm] * 3.0  # scale is irrelevant to MAC
    f_ae = f_ccx[perm] * f_ae_scale
    order = np.argsort(f_ae)
    return ModalSet(f_ae[order], f_ccx, ae_t[order], ccx_t, z)


def test_mac_identical_and_orthogonal() -> None:
    a = np.array([[1.0, 2.0, -1.0], [2.0, -1.0, 0.0]])
    assert a[0] @ a[1] == 0.0
    mac = mac_matrix(a, a)
    np.testing.assert_allclose(np.diag(mac), 1.0)
    np.testing.assert_allclose(mac[0, 1], 0.0, atol=1e-15)
    np.testing.assert_allclose(mac_matrix(a, -5.0 * a)[0, 0], 1.0)


def test_pair_by_mac_ignores_frequency_order() -> None:
    # CCX modes 1..3; the AE modes are the same shapes with the last two swapped.
    z = np.linspace(0, 1, 40)
    ccx = np.array([_shape(z, 1, 1), _shape(z, 1, 2), _shape(z, 0, 2)])
    ae = ccx[[0, 2, 1]]
    pair, mac = pair_modes_by_mac(ae.reshape(3, -1), ccx.reshape(3, -1))
    assert pair.tolist() == [0, 2, 1]
    assert mac[pair[1], 1] == pytest.approx(1.0)


def test_signature_is_mesh_independent() -> None:
    sig_a = mode_signature(_modal_set(60, [(0, 2, 1.0)], 1.0).ccx_t, np.linspace(0, 100, 60))
    sig_b = mode_signature(_modal_set(333, [(0, 2, 1.0)], 1.0).ccx_t, np.linspace(0, 100, 333))
    assert sig_a.shape == sig_b.shape
    assert mac_matrix(sig_a, sig_b)[0, 0] > 0.999


def test_tracking_follows_a_frequency_order_swap() -> None:
    # Same two physical modes; on the coarse mesh y is the lower one, on the fine mesh x is.
    coarse = _modal_set(50, [(1, 2, 4.10), (0, 2, 4.18)], 1.0)
    fine = _modal_set(210, [(0, 2, 4.12), (1, 2, 4.15)], 1.0)
    track = track_modes(mode_signature(coarse.ccx_t, coarse.z), mode_signature(fine.ccx_t, fine.z))
    assert track.tolist() == [1, 0]  # coarse y -> fine mode 2, coarse x -> fine mode 1


def test_followed_gaps_use_the_same_physical_mode_and_report_mac() -> None:
    coarse = _modal_set(50, [(1, 2, 4.10), (0, 2, 4.18)], 1.02)
    fine = _modal_set(210, [(0, 2, 4.12), (1, 2, 4.15)], 1.005)
    out = followed_gaps({0.5: coarse, 0.25: fine}, fine_key=0.25)
    gaps_c, mac_c = out[0.5]
    gaps_f, mac_f = out[0.25]
    # Entry j is the finest CCX mode j+1: x first, y second, on both meshes.
    np.testing.assert_allclose(gaps_f, [0.005, 0.005])
    np.testing.assert_allclose(gaps_c, [0.02, 0.02])
    np.testing.assert_allclose(mac_c, 1.0)
    np.testing.assert_allclose(mac_f, 1.0)


def test_parse_frd_two_modes(tmp_path) -> None:
    def block(vals: dict[int, tuple[float, float, float]]) -> str:
        rows = "".join(
            f" -1{nid:10d}" + "".join(f"{v:12.5E}" for v in xyz) + "\n" for nid, xyz in vals.items()
        )
        return "  -4  DISP        4    1\n  -5  D1          1    2    1    0\n" + rows + " -3\n"

    text = (
        "    1C\n"
        + block({1: (1.0, 2.0, 3.0), 2: (4.0, 5.0, 6.0), 9: (9.0, 9.0, 9.0)})
        + block({1: (-1.0, 0.0, 0.5), 2: (0.0, 0.0, 2.0), 9: (9.0, 9.0, 9.0)})
    )
    frd = tmp_path / "t.frd"
    frd.write_text(text)
    modes = parse_frd_mode_translations(frd, n_nodes=2, n_modes=2)
    assert modes.shape == (2, 2, 3)
    np.testing.assert_allclose(modes[0], [[1, 2, 3], [4, 5, 6]])
    np.testing.assert_allclose(modes[1], [[-1, 0, 0.5], [0, 0, 2]])
