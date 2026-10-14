"""Scordelis-Lo roof, triangular MITC3+ and quad MITC4+ against Ko et al. (2017) Table 10.

The live test ``test_3_4_scordelis_lo_tables_10_to_11`` runs only the quad path. This file drives
the same cylindrical-patch mesh, load, boundary conditions and normalization with triangles
(three-node cells select the MITC3+ element), so the triangular column of the published table is
compared too. It replaces the hand-driven, unrevisioned out-of-band measurements recorded in
``docs/validation/rows/10-out-of-band.yaml``.

Reference: Table 10 (regular mesh), normalised by ``w_ref = 3.0240e-1`` (MacNeal & Harder 1985).
"""

from __future__ import annotations

import numpy as np
import pytest

pytest.importorskip("_aeroelast", reason="Rust backend not available")

from tests.support.assertions import assert_relative_error  # noqa: E402
from tests.validation.benchmarks.test_ko2017_performance import (  # noqa: E402
    DOF,
    MAT_SC,
    _build_cylindrical_patch,
    _Case,
    _find_node_by_xyz,
    _gravity_load,
    _run_case,
    _scordelis_fixed,
)

R = 25.0
L_HALF = 25.0
ANGLE_DEG = 40.0
THICKNESS = 0.25
WREF = 3.0240e-1
TOL = 0.05  # suite rule; the paper cells are reproduced far inside it (see the store row)

# (triangular, N, Table 10 cell)
CELLS = [
    (False, 16, 0.9973),
    (True, 8, 0.8577),
    (True, 16, 0.9550),
    (True, 32, 0.9851),
]


@pytest.mark.parametrize("triangular, n, table_value", CELLS)
def test_scordelis_lo_table_10(triangular, n, table_value):
    def build():
        return _build_cylindrical_patch(
            radius=R, length=L_HALF, angle_deg=ANGLE_DEG, nx=n, ny=n, triangular=triangular
        )

    def measure(mesh, m, u):
        node = _find_node_by_xyz(
            mesh,
            float(R * np.cos(np.radians(ANGLE_DEG))),
            float(R * np.sin(np.radians(ANGLE_DEG))),
            0.0,
        )
        return float(u[m[node.id] * DOF])

    name = f"scordelis_lo_reg_{'MITC3' if triangular else 'MITC4'}_N{n}"
    case = _Case(
        name=name,
        element="MITC3" if triangular else "MITC4",
        build_mesh=build,
        material=MAT_SC,
        thickness=THICKNESS,
        load_vector=lambda mesh, m: _gravity_load(
            mesh, m, rho=360.0, g=1.0, thickness=THICKNESS, dof=0
        ),
        fixed_dofs=lambda mesh, m: _scordelis_fixed(mesh, m, length=L_HALF),
        measure=measure,
        expected_normalized=table_value,
        wref=WREF,
    )
    norm = _run_case(case)
    assert_relative_error(
        norm,
        table_value,
        tol=TOL,
        reference_name="Ko et al. 2017 Table 10 (regular mesh) cell",
        kind="paper",
        what=name,
    )
