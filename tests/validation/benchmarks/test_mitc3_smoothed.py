"""Strain-smoothed MITC3+ shell element (Lee & Lee 2019), end to end.

Reference:
  C. Lee, P.-S. Lee, "The strain-smoothed MITC3+ shell finite element",
  *Computers and Structures* 223 (2019) 106096.  Table 6, Mesh I, normalized
  vertical displacement at point B of the Scordelis-Lo roof.

These tests exercise the Rust union assembly exposed as
``_aeroelast.assemble_smoothed_mitc3``: the element-strain smoothing with edge
neighbours (Eqs. 15-18), the relative-frame convective transform, the union
rotation, and the global scatter.  Issue #2 was that this path over-stiffened a
curved shell by two orders of magnitude; the kernel now reproduces the paper.
"""

from __future__ import annotations

import numpy as np
import pytest
from scipy.sparse import coo_matrix

pytest.importorskip("_aeroelast", reason="Rust backend not available")
from _aeroelast import PyMeshAssembler, assemble_smoothed_mitc3

from tests.validation.benchmarks.test_ko2017_performance import (
    DOF,
    MAT_SC,
    _build_cylindrical_patch,
    _find_node_by_xyz,
    _gravity_load,
    _scordelis_fixed,
    _solve,
)

from tests.support.assertions import assert_residual_below  # noqa: E402

R = 25.0
LENGTH = 25.0
ANGLE_DEG = 40.0
THICKNESS = 0.25
WREF = 3.0240e-1

# Lee & Lee 2019, Table 6, Mesh I.  The plain MITC3+ column is the one the
# Ko et al. 2017 benchmark already pins; these are the paper's smoothed cells.
SMOOTHED_TABLE_6 = {8: 1.0323, 16: 1.0075}
MITC3_TABLE_6 = {8: 0.8793, 16: 0.9618}


def _mesh(n: int):
    return _build_cylindrical_patch(
        radius=R, length=LENGTH, angle_deg=ANGLE_DEG, nx=n, ny=n, triangular=True
    )


def _nodes(mesh):
    nodes = sorted(mesh.nodes, key=lambda node: node.id)
    return nodes, {node.id: i for i, node in enumerate(nodes)}


def _material_dict():
    return {
        "type": "isotropic",
        "e": float(MAT_SC.E),
        "nu": float(MAT_SC.nu),
        "rho": float(MAT_SC.rho),
        "thickness": THICKNESS,
        "shear_correction": 5.0 / 6.0,
        "drilling_scale": 1.0,
    }


def _assemble_smoothed(mesh):
    nodes, index = _nodes(mesh)
    ndof = len(nodes) * DOF
    connectivity = [[index[node.id] for node in elem.nodes] for elem in mesh.elements]
    coords = np.asarray([node.coords[:3] for node in nodes], dtype=float)
    rows, cols, vals = assemble_smoothed_mitc3(
        coords, connectivity, float(MAT_SC.E), float(MAT_SC.nu), THICKNESS
    )
    K = coo_matrix((vals, (rows, cols)), shape=(ndof, ndof)).tocsr()
    return K, index


def _assemble_plain_mitc3(mesh):
    nodes, index = _nodes(mesh)
    ndof = len(nodes) * DOF
    connectivity = [[index[node.id] for node in elem.nodes] for elem in mesh.elements]
    coords = np.asarray([node.coords[:3] for node in nodes], dtype=float)
    rows, cols, vals = PyMeshAssembler(
        node_coords=coords,
        connectivity=connectivity,
        elem_types=[3] * len(connectivity),
        materials=[_material_dict() for _ in mesh.elements],
    ).assemble_k()
    K = coo_matrix((vals, (rows, cols)), shape=(ndof, ndof)).tocsr()
    return K, index


def _normalized_tip(mesh, index, K) -> float:
    F = _gravity_load(mesh, index, rho=360.0, g=1.0, thickness=THICKNESS, dof=0)
    fixed = _scordelis_fixed(mesh, index, length=LENGTH)
    u = _solve(K, F, fixed)
    point = _find_node_by_xyz(
        mesh,
        float(R * np.cos(np.radians(ANGLE_DEG))),
        float(R * np.sin(np.radians(ANGLE_DEG))),
        0.0,
    )
    return abs(float(u[index[point.id] * DOF + 0])) / WREF


@pytest.mark.parametrize("n, expected", sorted(SMOOTHED_TABLE_6.items()))
def test_scordelis_lo_smoothed_mitc3_matches_lee_lee_table_6(n, expected):
    """The smoothed element must reproduce the paper's own Mesh-I cells."""
    mesh = _mesh(n)
    K, index = _assemble_smoothed(mesh)

    dense = K.todense()
    assert np.all(np.isfinite(dense)), "assembled K contains NaN/inf"
    assert np.allclose(dense, dense.T, rtol=1e-10, atol=1e-10 * float(np.abs(dense).max())), (
        "assembled K must be symmetric to round-off"
    )

    norm = _normalized_tip(mesh, index, K)
    rel = abs(norm - expected) / expected
    assert_residual_below(
        rel,
        tol=0.05,
        kind="paper",
        reference_name="Lee & Lee 2019, Table 6 Mesh I, normalized vertical displacement at B",
        what=f"N={n} Scordelis-Lo normalized displacement",
    )


def test_smoothing_relieves_scordelis_lo_membrane_locking_at_n8():
    """Smoothing must move the coarse-mesh result toward the reference plateau.

    The paper's whole point is that smoothing *improves* the triangular element
    on a curved shell instead of stiffening it.  Measured against the reference
    solution (normalized 1.0), the smoothed element must sit closer than the
    plain MITC3+ at N=8, where membrane locking is largest.
    """
    mesh = _mesh(8)
    K_smoothed, index = _assemble_smoothed(mesh)
    K_plain, index_plain = _assemble_plain_mitc3(mesh)
    assert index == index_plain

    smoothed = _normalized_tip(mesh, index, K_smoothed)
    plain = _normalized_tip(mesh, index, K_plain)

    # Guard against a vacuous comparison: the plain element must show the locking
    # and must itself sit on its own published column, so the two cells are read
    # the same way.
    assert plain < 0.95, f"expected membrane locking in the plain MITC3+, got {plain:.4f}"
    rel_plain = abs(plain - MITC3_TABLE_6[8]) / MITC3_TABLE_6[8]
    assert_residual_below(
        rel_plain,
        tol=0.05,
        kind="paper",
        reference_name="Lee & Lee 2019, Table 6 Mesh I, the plain MITC3+ column",
        what="N=8 plain MITC3+ normalized displacement",
    )
    assert abs(smoothed - 1.0) < abs(plain - 1.0), (
        f"smoothed {smoothed:.4f} must improve on MITC3+ {plain:.4f} toward 1.0"
    )
