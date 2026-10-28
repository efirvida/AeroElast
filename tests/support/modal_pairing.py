"""Shape-based pairing of modal results: within one mesh and across meshes.

Pairing two solvers' modes by frequency fails when modes are close: CCX modes 7 and 8 of the IEA
15 MW blade sit within 1 % of each other and swap order between 0.5 and 0.35 m (T12, group 5).
Everything here pairs by shape.

* :func:`mac_matrix` / :func:`pair_modes_by_mac`: Modal Assurance Criterion between two sets of
  shapes on the same nodes, and the one-to-one assignment that maximises it (Hungarian on
  ``1 - MAC``).
* :func:`mode_signature` / :func:`track_modes`: a mesh-independent signature (mean x / y
  translation in fixed span stations) and the assignment of the modes of one mesh to those of
  another, which is what follows one physical mode across a refinement.
* :func:`followed_gaps`: the frequency gap of every pair, indexed by the finest-mesh mode order.
* :func:`parse_frd_mode_translations`, :func:`check_corner_nodes_first`: CalculiX output and the
  node-numbering fact the MAC relies on (corner nodes keep their index in the quadratic deck).
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.optimize import linear_sum_assignment

#: Span stations of the cross-mesh signature (equal bins over the blade length).
N_STATIONS = 25


def parse_frd_mode_translations(frd: Path, n_nodes: int, n_modes: int) -> np.ndarray:
    """Return (n_modes, n_nodes, 3) translations of nodes 1..n_nodes from every DISP block."""
    modes: list[np.ndarray] = []
    current: np.ndarray | None = None
    with open(frd, encoding="utf-8", errors="replace") as handle:
        for line in handle:
            stripped = line.strip()
            if stripped.startswith("-4") and "DISP" in stripped.upper():
                current = np.zeros((n_nodes, 3))
                continue
            if current is None:
                continue
            if stripped.startswith("-3"):
                modes.append(current)
                current = None
                if len(modes) == n_modes:
                    break
                continue
            if not stripped.startswith("-1"):
                continue
            node_id = int(line[3:13])
            if 1 <= node_id <= n_nodes:
                current[node_id - 1] = [float(line[13:25]), float(line[25:37]), float(line[37:49])]
    return np.array(modes)


def check_corner_nodes_first(msh: Path, coords: np.ndarray) -> float:
    """Check that the first ``len(coords)`` nodes of a CCX mesh file are the corner nodes.

    The quadratic deck stacks the corner nodes first, so the corner translations of the FRD and of
    the AeroElast eigenvectors share one index.  Returns the largest coordinate error; raises
    ``RuntimeError`` above 1e-4 of the largest coordinate.
    """
    n_nodes = len(coords)
    with open(msh, encoding="utf-8") as handle:
        lines = handle.readlines()
    start = next(i for i, ln in enumerate(lines) if ln.upper().startswith("*NODE,"))
    deck = np.array(
        [[float(v) for v in ln.split(",")[1:4]] for ln in lines[start + 1 : start + 1 + n_nodes]]
    )
    err = float(np.max(np.abs(deck - coords)))
    if err > 1e-4 * float(np.max(np.abs(coords))):
        raise RuntimeError(f"corner nodes are not first in the CCX deck (max coord err {err})")
    return err


def mac_matrix(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """MAC between the rows of ``a`` (n_a, m) and of ``b`` (n_b, m)."""
    num = np.abs(a @ b.T) ** 2
    den = np.outer(np.einsum("ij,ij->i", a, a), np.einsum("ij,ij->i", b, b))
    return num / np.maximum(den, 1e-300)


def pair_modes_by_mac(a: np.ndarray, b: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """One-to-one MAC pairing of the rows of ``a`` and ``b`` (Hungarian on ``1 - MAC``).

    Returns ``(pair, mac)`` with ``pair[j]`` the row of ``a`` assigned to row ``j`` of ``b`` and
    ``mac`` the full (n_a, n_b) matrix.  Needs ``n_a >= n_b``.
    """
    mac = mac_matrix(a, b)
    rows, cols = linear_sum_assignment(1.0 - mac)
    pair = np.empty(b.shape[0], dtype=int)
    pair[cols] = rows
    return pair, mac


def mode_signature(shapes: np.ndarray, z: np.ndarray, n_bins: int = N_STATIONS) -> np.ndarray:
    """Mesh-independent signature of each mode: binned mean x / y translation.

    ``shapes`` is (n_modes, n_nodes, 3), ``z`` the span coordinate of the nodes.  The nodes are
    binned into ``n_bins`` equal stations over ``[z.min(), z.max()]``; the signature holds the mean
    x and the mean y translation of each station (empty stations count zero), normalised to unit
    length.  Two meshes of one blade give comparable signatures whatever their node count.
    """
    edges = np.linspace(float(z.min()), float(z.max()), n_bins + 1)
    idx = np.clip(np.digitize(z, edges[1:-1]), 0, n_bins - 1)
    counts = np.bincount(idx, minlength=n_bins)
    sig = np.zeros((shapes.shape[0], 2 * n_bins))
    for m in range(shapes.shape[0]):
        for k, comp in enumerate((0, 1)):
            sums = np.bincount(idx, weights=shapes[m, :, comp], minlength=n_bins)
            sig[m, k * n_bins : (k + 1) * n_bins] = np.where(
                counts > 0, sums / np.maximum(counts, 1), 0.0
            )
    norm = np.linalg.norm(sig, axis=1, keepdims=True)
    return sig / np.maximum(norm, 1e-300)


def track_modes(sig_from: np.ndarray, sig_to: np.ndarray) -> np.ndarray:
    """Assign each mode of one mesh to a mode of another by signature MAC (Hungarian).

    Returns ``track`` with ``track[k]`` the row of ``sig_to`` of mode ``k`` of ``sig_from``.
    """
    mac = mac_matrix(sig_from, sig_to)
    rows, cols = linear_sum_assignment(1.0 - mac)
    track = np.empty(sig_from.shape[0], dtype=int)
    track[rows] = cols
    return track


@dataclass(frozen=True)
class ModalSet:
    """One mesh solved by both solvers; the shapes are translations at the corner nodes."""

    freqs_ae: np.ndarray  # (n,) Hz, ascending
    freqs_ccx: np.ndarray  # (n,) Hz, ascending
    ae_t: np.ndarray  # (n, n_nodes, 3)
    ccx_t: np.ndarray  # (n, n_nodes, 3)
    z: np.ndarray  # (n_nodes,) span coordinate


def followed_gaps(
    sets: dict[float, ModalSet], fine_key: float
) -> dict[float, tuple[np.ndarray, np.ndarray]]:
    """Gaps and pairing MAC per mesh, indexed by the finest-mesh CCX mode order.

    On every mesh the AeroElast and CCX modes are paired by MAC.  The CCX modes of a coarser mesh
    are mapped to those of ``fine_key`` by signature (:func:`track_modes`), so entry ``j`` is the
    finest-mesh CCX mode ``j + 1`` as it appears on that mesh.  The gap is
    ``|f_ae - f_ccx| / f_ccx`` of the pair, the MAC is the pairing MAC of the pair.
    """
    fine = sets[fine_key]
    sig_fine = mode_signature(fine.ccx_t, fine.z)
    out: dict[float, tuple[np.ndarray, np.ndarray]] = {}
    for key, s in sets.items():
        n = len(s.freqs_ccx)
        pair, mac = pair_modes_by_mac(s.ae_t.reshape(len(s.freqs_ae), -1), s.ccx_t.reshape(n, -1))
        gap = np.abs(s.freqs_ae[pair] - s.freqs_ccx) / np.maximum(s.freqs_ccx, 1e-14)
        pair_mac = mac[pair, np.arange(n)]
        if key == fine_key:
            track = np.arange(n)
        else:
            track = track_modes(mode_signature(s.ccx_t, s.z), sig_fine)
        gaps = np.empty(n)
        macs = np.empty(n)
        gaps[track] = gap
        macs[track] = pair_mac
        out[key] = (gaps, macs)
    return out
