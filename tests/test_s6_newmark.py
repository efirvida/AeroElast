"""Tiny independent S-6 integration checks; no blade mesh or reference deck."""
import contextlib
import io
import sys
import unittest
from unittest.mock import patch
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from petsc4py import PETSc
from scipy.linalg import expm
from scipy.sparse import diags

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))
import run_s6_oneway_dynamic as s6
from aeroelast.core.assembler import MeshAssembler


def oscillator(integrate=s6.newmark_tip_series, stiffness=4.0, observer=None):
    """Unit mass, one free tip translation, all other shell DOFs clamped."""
    asm = SimpleNamespace(dofs_count=6, dofs_per_node=6, comm=PETSc.COMM_SELF)
    asm._coo_to_petsc = lambda r, c, v: MeshAssembler._coo_to_petsc(asm, r, c, v)
    mesh = SimpleNamespace(nodes=[SimpleNamespace(coords=(0., 0., 1.))])
    dt = .001
    forces = np.zeros((401, 6))
    forces[:, 1] = np.sin(np.linspace(0, 2 * np.pi, len(forces)))
    kwargs = {} if observer is None else {"observer": observer}
    hist = integrate(asm, mesh, diags([stiffness] * 6).tocsr(),
                     diags([1.] * 6).tocsr(), np.array([1]), forces, dt, **kwargs)
    # Independent exact continuous solution for piecewise-linear force on each
    # interval, with initial force = mean (consistent with the static IC).
    damping = s6.ETA_K * stiffness + s6.ETA_M
    generator = np.array([[0., 1., 0., 0.], [-stiffness, -damping, 1., 0.],
                          [0., 0., 0., 1.], [0., 0., 0., 0.]])
    advance = expm(dt * generator)
    mean = forces[:, 1].mean()
    state = np.array([mean / stiffness, 0., mean, 0.])
    exact = []
    for f in forces[:, 1]:
        state[3] = (f - state[2]) / dt
        state = advance @ state
        exact.append(state[0])
    return hist, np.asarray(exact)


class TestS6Newmark(unittest.TestCase):
    def test_component_capture_preserves_assembly(self):
        diagonal = np.arange(1., 7.)
        class Rust:
            def assemble_k(self):
                return np.arange(6), np.arange(6), diagonal

            def assemble_m(self):
                return np.arange(6), np.arange(6), np.ones(6)

        asm = SimpleNamespace(dofs_count=6, dofs_per_node=6, _rust=Rust())
        kg = PETSc.Mat().createAIJ(size=(6, 6), csr=(
            np.arange(7, dtype=PETSc.IntType), np.arange(6, dtype=PETSc.IntType),
            np.full(6, 2.)))
        asm.assemble_geometric_stiffness = lambda **kw: kg
        mesh = SimpleNamespace(node_count=1, node_id_to_index={},
                               get_node_set=lambda name: SimpleNamespace(nodes={}))
        baseline = s6.assemble_keff(asm, mesh, 2.)
        components = {}
        captured = s6.assemble_keff(asm, mesh, 2., components=components)
        for a, b in zip(baseline[:2], captured[:2], strict=True):
            np.testing.assert_array_equal(a.toarray(), b.toarray())
        np.testing.assert_array_equal(baseline[2], captured[2])
        np.testing.assert_array_equal(captured[0].diagonal(), diagonal + 2 - [4, 0, 4, 0, 0, 0])
        np.testing.assert_array_equal(components["Kg"].diagonal(), 2.)
        kg.destroy()

    def test_diagnostic_on_tiny_system(self):
        import diagnose_s6_divergence as diagnostic
        asm = SimpleNamespace(dofs_count=12, dofs_per_node=6, comm=PETSc.COMM_SELF)
        asm._coo_to_petsc = lambda r, c, v: MeshAssembler._coo_to_petsc(asm, r, c, v)
        mesh = SimpleNamespace(node_count=2, nodes=[
            SimpleNamespace(coords=(0., 0., 0.)), SimpleNamespace(coords=(0., 0., 1.))])
        k = diags(np.arange(1., 13.)).tocsr()
        m = diags(np.ones(12)).tocsr()

        def assemble(*args, components):
            components.update(K=k, Kg=0 * k, Ksp=0 * k)
            return k, m, np.arange(6, 12)

        output = io.StringIO()
        with patch.object(sys, "argv", ["diagnose_s6_divergence", "--eigenvalues", "2"]), \
                patch.object(diagnostic, "build_assembly", return_value=(asm, mesh)) as build, \
                patch.object(diagnostic, "assemble_keff", side_effect=assemble), \
                contextlib.redirect_stdout(output):
            diagnostic.main()
        build.assert_called_once_with(.25)
        text = output.getvalue()
        self.assertIn("extension:", text)
        self.assertIn("Krot eigenvalue=7.000000e+00", text)
        self.assertIn("Krot eigenvalue=8.000000e+00", text)
        self.assertIn("step=400:", text)
        self.assertIn("accepted_bound=25%", text)
        self.assertIn("no displacement feedback", text)
        self.assertIn("solve reasons=[", text)

    def test_forced_stable_oscillator(self):
        hist, exact = oscillator()
        np.testing.assert_allclose(hist[:, 1], exact, atol=2e-6, rtol=0)
        np.testing.assert_array_equal(hist[:, [0, 2]], 0.)

    def test_negative_stiffness_with_positive_step_matrix(self):
        # A positive dynamic step matrix does not imply physical stability.
        k = -100.
        c = s6.ETA_K * k + s6.ETA_M
        self.assertGreater(4 / .001**2 + 2 * c / .001 + k, 0.)
        hist, exact = oscillator(stiffness=k)
        np.testing.assert_allclose(hist[:, 1], exact, atol=2e-5, rtol=0)
        self.assertGreater(abs(hist[-1, 1]), .1)

    def test_observer_preserves_history_and_reports_balance(self):
        records = []
        hist, _ = oscillator(observer=lambda step, record: records.append((step, record)))
        baseline, _ = oscillator()
        np.testing.assert_array_equal(hist, baseline)
        self.assertEqual([step for step, _ in records], list(range(len(hist))))
        self.assertLess(max(r["relative_balance"] for _, r in records), 1e-9)
        self.assertTrue(all(r["solve_reason"] > 0 and r["static_solve_reason"] > 0
                            for _, r in records))
        np.testing.assert_allclose([r["translation_max"] for _, r in records],
                                   np.abs(hist[:, 1]))


if __name__ == "__main__":
    import _aeroelast._aeroelast as extension
    print(f"extension origin: {extension.__file__}", flush=True)
    unittest.main(verbosity=2)
