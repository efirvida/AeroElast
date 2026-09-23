from __future__ import annotations

import numpy as np

from aeroelast.core.mesh import ElementSet, ElementType, MeshElement, MeshModel, Node, NodeSet
from aeroelast.solvers.aero import QuasiSteadyVLMBackend, RotorAeroState, build_aero_runtime_context


def _build_toy_rotor_mesh() -> MeshModel:
    mesh = MeshModel()

    blade_1_nodes = [
        Node([0.0, 0.0, 1.0]),
        Node([1.0, 0.0, 1.0]),
        Node([1.0, 0.0, 2.0]),
        Node([0.0, 0.0, 2.0]),
    ]
    blade_2_nodes = [
        Node([0.0, 0.0, -1.0]),
        Node([-1.0, 0.0, -1.0]),
        Node([-1.0, 0.0, -2.0]),
        Node([0.0, 0.0, -2.0]),
    ]

    for node in blade_1_nodes + blade_2_nodes:
        mesh.add_node(node)

    blade_1_elem = MeshElement(blade_1_nodes, ElementType.quad)
    blade_2_elem = MeshElement(blade_2_nodes, ElementType.quad)
    mesh.add_element(blade_1_elem)
    mesh.add_element(blade_2_elem)

    mesh.add_node_set(NodeSet("rotor_blade_1", set(blade_1_nodes)))
    mesh.add_node_set(NodeSet("rotor_blade_2", set(blade_2_nodes)))
    mesh.add_element_set(ElementSet("rotor_blade_1", {blade_1_elem}))
    mesh.add_element_set(ElementSet("rotor_blade_2", {blade_2_elem}))

    return mesh


def _build_runtime_context() -> tuple[MeshModel, object]:
    mesh = _build_toy_rotor_mesh()
    cfg = {
        "participant": "Fluid",
        "coupling_mesh": "Fluid-Mesh",
        "mesh": {
            "source": "generator",
            "generator": {"type": "RotorMesh", "params": {"yaml_file": "blade.yaml"}},
        },
        "aero": {
            "backend": "vlm",
            "blade_file": "blade.yaml",
            "rotor": {
                "n_blades": 2,
                "hub_radius": 1.0,
                "rotation_axis": [0.0, 1.0, 0.0],
                "rotation_center": [0.0, 0.0, 0.0],
            },
            "backend_config": {
                "wind_velocity": [10.0, -2.0, 0.0],
                "wake_length_scale": 15.0,
                "diagonal_regularization": 1.0e-6,
            },
        },
    }
    return mesh, build_aero_runtime_context(mesh, cfg)


def test_vlm_backend_returns_zero_for_tangential_flow_without_normal_component():
    mesh, runtime_context = _build_runtime_context()
    backend = QuasiSteadyVLMBackend(diagonal_regularization=1.0e-6)
    state = RotorAeroState(
        wind_velocity=np.array([10.0, 0.0, 0.0]),
        nodal_displacements=np.zeros((mesh.node_count, 3), dtype=float),
    )

    loads = backend.compute_loads(runtime_context, state)

    np.testing.assert_allclose(loads.nodal_forces, 0.0, atol=1.0e-9)
    np.testing.assert_allclose(loads.integrated_force, 0.0, atol=1.0e-9)


def test_vlm_backend_generates_normal_force_for_nonzero_incidence():
    mesh, runtime_context = _build_runtime_context()
    backend = QuasiSteadyVLMBackend(diagonal_regularization=1.0e-6)
    state = RotorAeroState(
        wind_velocity=np.array([10.0, -2.0, 0.0]),
        nodal_displacements=np.zeros((mesh.node_count, 3), dtype=float),
    )

    loads = backend.compute_loads(runtime_context, state)

    assert loads.nodal_forces.shape == (mesh.node_count, 3)
    assert abs(loads.integrated_force[1]) > 1.0e-6
    np.testing.assert_allclose(loads.integrated_force[[0, 2]], [0.0, 0.0], atol=1.0e-6)
    assert loads.metadata["circulation"].shape == (runtime_context.surface.n_panels,)


def test_vlm_backend_accounts_for_structural_velocity_in_relative_flow():
    mesh, runtime_context = _build_runtime_context()
    backend = QuasiSteadyVLMBackend(diagonal_regularization=1.0e-6)
    state_without_velocity = RotorAeroState(
        wind_velocity=np.array([10.0, -2.0, 0.0]),
        nodal_displacements=np.zeros((mesh.node_count, 3), dtype=float),
    )
    state_with_velocity = RotorAeroState(
        wind_velocity=np.array([10.0, -2.0, 0.0]),
        nodal_displacements=np.zeros((mesh.node_count, 3), dtype=float),
        nodal_velocities=np.tile([0.0, -2.0, 0.0], (mesh.node_count, 1)),
    )

    loads_without_velocity = backend.compute_loads(runtime_context, state_without_velocity)
    loads_with_velocity = backend.compute_loads(runtime_context, state_with_velocity)

    assert abs(loads_with_velocity.integrated_force[1]) < abs(
        loads_without_velocity.integrated_force[1]
    )


def test_vlm_backend_persists_wake_rows_across_time_windows():
    mesh, runtime_context = _build_runtime_context()
    backend = QuasiSteadyVLMBackend(
        diagonal_regularization=1.0e-6,
        wake_history_steps=2,
    )
    state = RotorAeroState(
        dt=0.25,
        omega_rad_s=1.5,
        wind_velocity=np.array([10.0, -2.0, 0.0]),
        nodal_displacements=np.zeros((mesh.node_count, 3), dtype=float),
    )

    loads_first = backend.compute_loads(runtime_context, state)
    backend.finalize_time_window()
    checkpoint_state = backend.capture_state()
    loads_second = backend.compute_loads(runtime_context, state)

    assert loads_first.metadata["wake_row_count"] == 1
    assert len(checkpoint_state["committed_wake_rows"]) == 1
    assert loads_second.metadata["wake_row_count"] == 2
    assert loads_second.metadata["wake_rows"].shape == (2, runtime_context.surface.n_panels, 2, 3)


def test_vlm_backend_restores_wake_state_from_checkpoint():
    mesh, runtime_context = _build_runtime_context()
    backend = QuasiSteadyVLMBackend(
        diagonal_regularization=1.0e-6,
        wake_history_steps=3,
    )
    state = RotorAeroState(
        dt=0.25,
        omega_rad_s=1.5,
        wind_velocity=np.array([10.0, -2.0, 0.0]),
        nodal_displacements=np.zeros((mesh.node_count, 3), dtype=float),
    )

    backend.compute_loads(runtime_context, state)
    backend.finalize_time_window()
    checkpoint_state = backend.capture_state()

    backend.compute_loads(runtime_context, state)
    backend.finalize_time_window()
    assert len(backend.capture_state()["committed_wake_rows"]) == 2

    backend.restore_state(checkpoint_state)
    restored_state = backend.capture_state()

    assert len(restored_state["committed_wake_rows"]) == 1
    np.testing.assert_allclose(
        restored_state["committed_wake_rows"][0].points,
        checkpoint_state["committed_wake_rows"][0].points,
    )


def test_vlm_backend_does_not_shed_wake_from_blade_motion_alone():
    mesh, runtime_context = _build_runtime_context()
    backend = QuasiSteadyVLMBackend(
        diagonal_regularization=1.0e-6,
        wake_history_steps=2,
    )
    state = RotorAeroState(
        dt=0.25,
        omega_rad_s=2.0,
        wind_velocity=np.zeros(3, dtype=float),
        nodal_displacements=np.zeros((mesh.node_count, 3), dtype=float),
    )

    loads = backend.compute_loads(runtime_context, state)

    assert loads.metadata["wake_row_count"] == 0
    assert loads.metadata["wake_rows"].size == 0


def test_vlm_backend_reports_wake_convection_velocities():
    mesh, runtime_context = _build_runtime_context()
    backend = QuasiSteadyVLMBackend(
        diagonal_regularization=1.0e-6,
        wake_history_steps=2,
    )
    state = RotorAeroState(
        dt=0.25,
        omega_rad_s=1.5,
        wind_velocity=np.array([10.0, -2.0, 0.0]),
        nodal_displacements=np.zeros((mesh.node_count, 3), dtype=float),
    )

    loads = backend.compute_loads(runtime_context, state)

    assert loads.metadata["wake_row_count"] == 1
    assert loads.metadata["wake_convection_velocities"].shape == (
        1,
        runtime_context.surface.n_panels,
        2,
        3,
    )


def test_vlm_backend_can_correct_wake_geometry_with_second_pass():
    mesh, runtime_context = _build_runtime_context()
    backend = QuasiSteadyVLMBackend(
        diagonal_regularization=1.0e-6,
        wake_history_steps=2,
        wake_corrector_iterations=1,
    )
    state = RotorAeroState(
        dt=0.25,
        omega_rad_s=1.5,
        wind_velocity=np.array([10.0, -2.0, 0.0]),
        nodal_displacements=np.zeros((mesh.node_count, 3), dtype=float),
    )

    backend.compute_loads(runtime_context, state)
    backend.finalize_time_window()
    loads = backend.compute_loads(runtime_context, state)

    assert loads.metadata["wake_corrector_iterations_used"] == 1
    assert loads.metadata["wake_predictor_rows"].shape == loads.metadata["wake_rows"].shape
    assert not np.allclose(loads.metadata["wake_predictor_rows"], loads.metadata["wake_rows"])
