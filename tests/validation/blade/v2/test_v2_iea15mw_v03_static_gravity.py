"""V-03 on the windIO v2 blade: integrated mass and parked gravity response (v2 counterpart of group 41).

Same references, metrics and tolerances as ``tests/test_iea15mw_v03_static_gravity.py``:

* mass: assembled total vs the IEA tabular declared 67,921 kg, ``rel_err < 0.05`` (the L3 ceiling);
* guards (not store rows in v1): tip edgewise deflection in [1.2, 4.0] m, root reaction vs the assembled
  weight within ``_REACTION_TOL = 0.02``, downward sign.

Only the mesh differs: ``build_v2_blade(0.25, skin_reference="mid")`` instead of
``BladeMesh(...).generate(renumber="rcm")`` + ``build_rust_properties``.
See docs/validation/_pending/66-v2_blade_total_mass/NOTES.md.
"""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pytest

pytest.importorskip("petsc4py", reason="PETSc not available")

from tests.support.assertions import assert_relative_error  # noqa: E402
from tests.support.v2_blade import build_v2_blade  # noqa: E402

_BLADE_MASS_KG = 67_921.0  # IEA-15-240-RWT_tabular.xlsx, Overview
_G = 9.81
_TIP_EDGE_DEFL_MIN = 1.2  # m, same order-of-magnitude bounds as v1
_TIP_EDGE_DEFL_MAX = 4.0
_REACTION_TOL = 0.02
_MASS_TOL = 0.05  # L3 ceiling, identical to v1
ELEMENT_SIZE = 0.25
SKIN_REFERENCE = "mid"


def _model(properties):
    from aeroelast.elements import ElementFamily

    return {
        "elements": {
            "element_family": ElementFamily.SHELL,
            "span_direction": (0.0, 0.0, 1.0),
            "properties": properties,
        },
    }


@pytest.mark.slow
class TestStaticGravityDeflectionV2:
    """V-03a on the v2 blade: clamped root, horizontal, gravity -Y."""

    @pytest.fixture(scope="class")
    def static_result(self):
        from scipy.sparse import coo_matrix

        from aeroelast.core.bc import BodyForce, DirichletCondition
        from aeroelast.elements import ElementFamily
        from aeroelast.solvers.elasticity.static_linear import StaticLinearSolver

        mesh, properties = build_v2_blade(ELEMENT_SIZE, skin_reference=SKIN_REFERENCE)
        cfg = {
            "solver": {},
            "elements": {
                "element_family": ElementFamily.SHELL,
                "span_direction": (0.0, 0.0, 1.0),
                "properties": properties,
            },
        }
        solver = StaticLinearSolver(mesh, cfg)
        dpn = solver.domain.dofs_per_node
        root_ids = sorted(mesh.get_node_set("RootNodes").nodes.keys())
        idx = mesh.node_id_to_index
        root_dofs = sorted([idx[n] * dpn + d for n in root_ids for d in range(dpn)])
        solver.add_dirichlet_conditions([DirichletCondition(root_dofs, 0.0)])
        solver.add_body_forces([BodyForce([0.0, -_G, 0.0])])
        u = solver.solve()

        coords = mesh.coords_array
        tip = int(np.argmax(coords[:, 2]))
        tip_disp = u[tip * dpn : tip * dpn + dpn]

        m_rows, _, m_vals = solver.domain._rust.assemble_m()
        mass = float(np.asarray(m_vals)[np.asarray(m_rows) % dpn == 0].sum())

        k_rows, k_cols, k_vals = solver.domain._rust.assemble_k()
        n_dof = solver.domain.dofs_count
        K = coo_matrix((k_vals, (k_rows, k_cols)), shape=(n_dof, n_dof)).tocsr()
        react = (K @ u - solver._build_f_ext())[root_dofs]
        return SimpleNamespace(
            tip_displacement=tip_disp,
            root_reaction_y=float(react[1::dpn].sum()),
            assembled_weight_n=mass * _G,
        )

    def test_tip_edgewise_deflection_in_range(self, static_result):
        u_edge = abs(float(static_result.tip_displacement[1]))
        assert _TIP_EDGE_DEFL_MIN <= u_edge <= _TIP_EDGE_DEFL_MAX, (
            f"Tip edgewise deflection {u_edge:.3f} m outside [{_TIP_EDGE_DEFL_MIN}, {_TIP_EDGE_DEFL_MAX}] m"
        )

    def test_root_reaction_balances_weight(self, static_result):
        r_y = abs(static_result.root_reaction_y)
        w = static_result.assembled_weight_n
        rel_err = abs(r_y - w) / w
        assert rel_err < _REACTION_TOL, (
            f"Root reaction {r_y:.0f} N vs assembled weight {w:.0f} N (err={rel_err:.2%}, tol={_REACTION_TOL:.0%})"
        )

    def test_deflection_direction_downward(self, static_result):
        u_edge = float(static_result.tip_displacement[1])
        assert u_edge < 0, f"Expected negative (downward) edgewise deflection, got {u_edge:.3f} m"


class TestBladeMassIntegrationV2:
    """V-03b on the v2 blade: assembled mass vs the declared 67,921 kg (no static solve)."""

    @pytest.fixture(scope="class")
    def assembled_mass_kg(self):
        from aeroelast.core.assembler import MeshAssembler

        mesh, properties = build_v2_blade(ELEMENT_SIZE, skin_reference=SKIN_REFERENCE)
        assembler = MeshAssembler(mesh, _model(properties))
        assembler.assemble_mass_matrix().destroy()
        m_rows, _, m_vals = assembler._rust.assemble_m()
        return float(np.asarray(m_vals)[np.asarray(m_rows) % 6 == 0].sum())

    def test_total_blade_mass_within_5pct(self, assembled_mass_kg):
        """Assembled mass within +-5 % of 67,921 kg (IEA tabular, L3 ceiling), as the v1 test."""
        assert_relative_error(
            assembled_mass_kg,
            _BLADE_MASS_KG,
            tol=_MASS_TOL,
            kind="code",
            reference_name="IEA-15-240-RWT_tabular.xlsx, Overview declared mass, 67,921 kg",
            what="V-03 v2 assembled total blade mass [kg]",
        )
