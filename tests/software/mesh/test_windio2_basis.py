"""Fiber-angle frame of the shell basis shared by the gmsh mesher (moved from test_windio2_web_mesh, T12, S9)."""
import numpy as np
import pytest


def test_explicit_physical_fiber_vectors_to_adapter_angles():
    from aeroelast.models.blade.numad.mesh_gen.windio2_basis import physical_fiber_angle
    # Rows are physical element x, y, normal. Independent rotated shell plane.
    basis = np.array([[0., 0., 1.], [0., 1., 0.], [-1., 0., 0.]])
    for angle in (0., 31., -47., 90.):
        a = np.deg2rad(angle)
        vector = np.array([0., np.sin(a), np.cos(a)])
        assert physical_fiber_angle(vector, basis) == pytest.approx(angle, abs=2e-14)
