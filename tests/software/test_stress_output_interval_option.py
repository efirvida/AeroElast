"""``solver.stress_output_interval`` reaches the rotor solver (odd/tasks/l5-one-way.md, T15 S1).

The key used to be dropped silently by ``SolverConfig`` and the ``FSIRunner._build_model_config`` whitelist,
so ``rotor.py`` always kept its default of 1 (stress recovery every preCICE window).
"""

from __future__ import annotations

import pytest

from aeroelast.core.config import FSISimulationConfig
from aeroelast.solvers.fsi.runner import FSIRunner


def _data(**solver_extra):
    solver = {"type": "LinearDynamicFSIRotor", "total_time": 1.0, "time_step": 0.01}
    solver.update(solver_extra)
    return {
        "mesh": {"source": "file", "file": {"path": "dummy.h5", "format": "hdf5"}},
        "material": {"type": "isotropic", "name": "m", "E": 1.0e9, "nu": 0.3, "rho": 1000.0},
        "elements": {"family": "SHELL", "thickness": 0.01},
        "solver": solver,
    }


def _model_config(cfg):
    runner = object.__new__(FSIRunner)
    runner.config = cfg
    runner.config_path = None
    runner._blade_properties = None
    runner._material = None
    from rich.console import Console

    runner._console = Console(quiet=True)
    return runner._build_model_config()


def test_value_reaches_solver_params():
    cfg = FSISimulationConfig.from_dict(_data(stress_output_interval=50))
    assert cfg.solver.stress_output_interval == 50
    assert _model_config(cfg)["solver"]["stress_output_interval"] == 50


def test_absent_key_keeps_rotor_default():
    cfg = FSISimulationConfig.from_dict(_data())
    assert cfg.solver.stress_output_interval is None
    assert "stress_output_interval" not in _model_config(cfg)["solver"]


@pytest.mark.parametrize("bad", [0, -3])
def test_non_positive_rejected(bad):
    with pytest.raises(ValueError, match="stress_output_interval must be a positive integer"):
        FSISimulationConfig.from_dict(_data(stress_output_interval=bad))
