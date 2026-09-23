"""Tests for generic aerodynamic participant config normalization."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import yaml

from aeroelast.cli.run_aero_fsi import (
    _TEMPLATE,
    _configure_sharpy_mesh_exports,
    _pin_relative_mesh_output_to_workdir,
)
from aeroelast.core.config import FSISimulationConfig, SolverType
from aeroelast.core.mesh import ElementSet, ElementType, MeshElement, MeshModel, Node, NodeSet
from aeroelast.solvers.aero import (
    RotorAeroFSIParticipant,
    build_aero_participant_from_config,
    build_aero_runtime_context,
    normalize_aero_fsi_config,
    resolve_aero_backend,
)
from aeroelast.solvers.aero.sharpy_fsi_adapter import SharpyRotorFSICaseAdapter


def _build_toy_rotor_mesh() -> MeshModel:
    mesh = MeshModel()

    blade_1_nodes = [
        Node([0.0, 0.0, 1.0]),
        Node([1.0, 0.0, 1.0]),
        Node([1.0, 1.0, 1.0]),
        Node([0.0, 1.0, 1.0]),
    ]
    blade_2_nodes = [
        Node([0.0, 0.0, -1.0]),
        Node([-1.0, 0.0, -1.0]),
        Node([-1.0, 1.0, -1.0]),
        Node([0.0, 1.0, -1.0]),
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


def test_resolve_aero_backend_from_legacy_bem_section():
    assert resolve_aero_backend({"bem": {"wind_speed": 10.0}}) == "bem"


def test_normalize_aero_fsi_config_merges_generic_bem_settings():
    cfg = {
        "participant": "Fluid",
        "config_file": "precice-config.xml",
        "mesh": {
            "source": "generator",
            "generator": {
                "type": "RotorMesh",
                "params": {"yaml_file": "blade.yaml", "element_size": 0.5},
            },
        },
        "aero": {
            "backend": "bem",
            "blade_file": "blade.yaml",
            "rotor": {
                "n_blades": 4,
                "hub_radius": 7.5,
                "rotation_axis": [0.0, 1.0, 0.0],
                "rotation_center": [0.0, 0.0, 0.0],
            },
            "backend_config": {
                "wind_speed": 10.0,
                "pitch": 2.0,
            },
        },
    }

    normalized = normalize_aero_fsi_config(cfg)

    assert normalized["aero"]["backend"] == "bem"
    assert normalized["blade_file"] == "blade.yaml"
    assert normalized["bem"]["n_blades"] == 4
    assert normalized["bem"]["hub_radius"] == 7.5
    assert normalized["bem"]["wind_speed"] == 10.0
    assert normalized["bem"]["pitch"] == 2.0
    assert normalized["mesh"]["generator"]["params"]["n_blades"] == 4
    assert normalized["mesh"]["generator"]["params"]["hub_radius"] == 7.5
    assert "blade_file" not in cfg
    assert "bem" not in cfg


def test_normalize_aero_fsi_config_preserves_legacy_bem_overrides():
    cfg = {
        "blade_file": "legacy_blade.yaml",
        "bem": {"wind_speed": 11.0, "hub_radius": 3.0},
        "aero": {
            "backend": "bem",
            "rotor": {"hub_radius": 7.5},
            "backend_config": {"pitch": 1.0},
        },
    }

    normalized = normalize_aero_fsi_config(cfg)

    assert normalized["blade_file"] == "legacy_blade.yaml"
    assert normalized["bem"]["hub_radius"] == 3.0
    assert normalized["bem"]["wind_speed"] == 11.0
    assert normalized["bem"]["pitch"] == 1.0


def test_normalize_aero_fsi_config_keeps_generic_vlm_backend_config():
    cfg = {
        "aero": {
            "backend": "vlm",
            "blade_file": "blade.yaml",
            "rotor": {"n_blades": 2, "hub_radius": 5.0},
            "backend_config": {
                "wind_velocity": [10.0, -1.0, 0.0],
                "wake_length_scale": 30.0,
                "wake_history_steps": 6,
                "wake_corrector_iterations": 1,
            },
        },
        "mesh": {
            "source": "generator",
            "generator": {"type": "RotorMesh", "params": {"yaml_file": "blade.yaml"}},
        },
    }

    normalized = normalize_aero_fsi_config(cfg)

    assert normalized["aero"]["backend"] == "vlm"
    assert normalized["vlm"]["wake_length_scale"] == 30.0
    assert normalized["vlm"]["wake_history_steps"] == 6
    assert normalized["vlm"]["wake_corrector_iterations"] == 1
    assert normalized["vlm"]["wind_velocity"] == [10.0, -1.0, 0.0]
    assert normalized["mesh"]["generator"]["params"]["n_blades"] == 2
    assert normalized["mesh"]["generator"]["params"]["hub_radius"] == 5.0


def test_normalize_aero_fsi_config_keeps_generic_sharpy_backend_config():
    cfg = {
        "aero": {
            "backend": "sharpy",
            "blade_file": "blade.yaml",
            "rotor": {"n_blades": 3, "hub_radius": 4.0},
            "backend_config": {
                "model": "uvlm",
                "air_density": 1.2,
                "solver_settings": {"convection_scheme": 3},
            },
        },
        "mesh": {
            "source": "generator",
            "generator": {"type": "RotorMesh", "params": {"yaml_file": "blade.yaml"}},
        },
    }

    normalized = normalize_aero_fsi_config(cfg)

    assert normalized["aero"]["backend"] == "sharpy"
    assert normalized["sharpy"]["model"] == "uvlm"
    assert normalized["sharpy"]["air_density"] == 1.2
    assert normalized["sharpy"]["solver_settings"]["convection_scheme"] == 3
    assert normalized["mesh"]["generator"]["params"]["n_blades"] == 3
    assert normalized["mesh"]["generator"]["params"]["hub_radius"] == 4.0


def test_fsisimulationconfig_infers_aero_fsi_solver_type():
    cfg = {
        "mesh": {
            "source": "generator",
            "generator": {
                "type": "BladeMesh",
                "params": {"yaml_file": "blade.yaml"},
            },
        },
        "aero": {"backend": "bem", "blade_file": "blade.yaml"},
    }

    parsed = FSISimulationConfig.from_dict(cfg)

    assert parsed.solver.type == SolverType.AERO_FSI.value


def test_run_aero_fsi_template_is_valid_yaml():
    parsed = yaml.safe_load(_TEMPLATE)

    assert parsed["aero"]["backend"] == "bem"
    assert parsed["aero"]["rotor"]["n_blades"] == 3
    assert parsed["output"]["sectional_bins"] == 8


def test_pin_relative_mesh_output_to_workdir_rewrites_relative_path():
    cfg = {"mesh": {"output_file": "fluid_mesh.vtu"}}

    _pin_relative_mesh_output_to_workdir(cfg, Path("/tmp/case"))

    assert cfg["mesh"]["output_file"] == "/tmp/case/fluid_mesh.vtu"


def test_configure_sharpy_mesh_exports_preserves_aerogrid_path_and_moves_coupling_mesh():
    cfg = {"mesh": {"output_file": "fluid_mesh.vtu"}}

    _configure_sharpy_mesh_exports(cfg, Path("/tmp/case"), "sharpy")

    assert cfg["mesh"]["output_file"] == "/tmp/case/fluid_mesh_coupling.vtu"
    assert cfg["output"]["aero_mesh_file"] == "/tmp/case/fluid_mesh.vtu"
    assert cfg["output"]["coupling_mesh_file"] == "/tmp/case/fluid_mesh_coupling.vtu"


def test_build_aero_runtime_context_builds_rotor_geometry_and_surface():
    mesh = _build_toy_rotor_mesh()
    cfg = {
        "participant": "Fluid",
        "coupling_mesh": "Fluid-Mesh",
        "mesh": {
            "source": "generator",
            "generator": {
                "type": "RotorMesh",
                "params": {"yaml_file": "blade.yaml"},
            },
        },
        "aero": {
            "backend": "bem",
            "blade_file": "blade.yaml",
            "rotor": {
                "n_blades": 2,
                "hub_radius": 3.5,
                "rotation_axis": [0.0, 1.0, 0.0],
                "rotation_center": [0.0, 0.0, 0.0],
            },
            "backend_config": {"wind_speed": 10.0},
        },
    }

    runtime_context = build_aero_runtime_context(mesh, cfg)

    assert runtime_context.backend == "bem"
    assert runtime_context.participant_name == "Fluid"
    assert runtime_context.coupling_mesh_name == "Fluid-Mesh"
    assert runtime_context.blade_file == "blade.yaml"
    assert runtime_context.geometry is not None
    assert runtime_context.surface is not None
    assert runtime_context.geometry.n_blades == 2
    np.testing.assert_allclose(runtime_context.geometry.rotation_axis, [0.0, 1.0, 0.0])
    np.testing.assert_allclose(runtime_context.geometry.rotation_center, [0.0, 0.0, 0.0])
    assert runtime_context.geometry.hub_radius == 3.5
    assert runtime_context.surface.n_blades == 2
    assert all(blade.root_reference_radius == 3.5 for blade in runtime_context.surface.blades)


def test_build_aero_runtime_context_skips_generic_surface_for_sharpy_backend():
    mesh = _build_toy_rotor_mesh()
    cfg = {
        "participant": "Fluid",
        "coupling_mesh": "Fluid-Mesh",
        "mesh": {
            "source": "generator",
            "generator": {
                "type": "RotorMesh",
                "params": {"yaml_file": "blade.yaml"},
            },
        },
        "aero": {
            "backend": "sharpy",
            "blade_file": "blade.yaml",
            "rotor": {
                "n_blades": 2,
                "hub_radius": 3.5,
                "rotation_axis": [0.0, 1.0, 0.0],
                "rotation_center": [0.0, 0.0, 0.0],
            },
            "backend_config": {"model": "vlm"},
        },
    }

    runtime_context = build_aero_runtime_context(mesh, cfg)

    assert runtime_context.backend == "sharpy"
    assert runtime_context.geometry is not None
    assert runtime_context.surface is None


def test_build_aero_runtime_context_skips_rotor_views_for_non_rotor_mesh():
    mesh = MeshModel(
        nodes=[Node([0.0, 0.0, 0.0]), Node([1.0, 0.0, 0.0]), Node([0.0, 1.0, 0.0])],
        elements=[],
    )
    cfg = {
        "bem": {"wind_speed": 10.0, "n_blades": 3, "hub_radius": 2.0},
        "blade_file": "blade.yaml",
    }

    runtime_context = build_aero_runtime_context(mesh, cfg)

    assert runtime_context.backend == "bem"
    assert runtime_context.geometry is None
    assert runtime_context.surface is None
    assert runtime_context.rotor["n_blades"] == 3
    assert runtime_context.rotor["hub_radius"] == 2.0


def test_build_aero_participant_from_config_builds_generic_vlm_participant():
    mesh = _build_toy_rotor_mesh()
    cfg = {
        "participant": "Fluid",
        "config_file": "precice-config.xml",
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
                "hub_radius": 3.5,
                "rotation_axis": [0.0, 1.0, 0.0],
                "rotation_center": [0.0, 0.0, 0.0],
            },
            "backend_config": {
                "wind_velocity": [10.0, -1.0, 0.0],
                "wake_length_scale": 20.0,
                "wake_history_steps": 4,
                "wake_corrector_iterations": 2,
            },
        },
    }

    participant = build_aero_participant_from_config(mesh, cfg)

    assert isinstance(participant, RotorAeroFSIParticipant)
    assert participant.backend.name == "vlm"
    assert participant.runtime_context.backend == "vlm"
    assert participant.backend._wake_history_steps == 4
    assert participant.backend._wake_corrector_iterations == 2


def test_build_aero_participant_from_config_builds_generic_sharpy_participant():
    mesh = _build_toy_rotor_mesh()
    cfg = {
        "participant": "Fluid",
        "config_file": "precice-config.xml",
        "coupling_mesh": "Fluid-Mesh",
        "mesh": {
            "source": "generator",
            "generator": {"type": "RotorMesh", "params": {"yaml_file": "blade.yaml"}},
        },
        "aero": {
            "backend": "sharpy",
            "blade_file": "blade.yaml",
            "rotor": {
                "n_blades": 2,
                "hub_radius": 3.5,
                "rotation_axis": [0.0, 1.0, 0.0],
                "rotation_center": [0.0, 0.0, 0.0],
            },
            "backend_config": {
                "model": "vlm",
                "air_density": 1.18,
                "solver_settings": {"horseshoe": True},
            },
        },
    }

    participant = build_aero_participant_from_config(mesh, cfg)

    assert isinstance(participant, RotorAeroFSIParticipant)
    assert participant.backend.name == "sharpy"
    assert participant.backend.model == "vlm"
    assert participant.backend.solver_settings["rho"] == 1.18
    assert "density" not in participant.backend.solver_settings
    assert participant.backend.solver_settings["horseshoe"] is True
    assert isinstance(participant.backend._case_adapter, SharpyRotorFSICaseAdapter)
    assert participant.runtime_context.surface is None
