"""A leftover ``excel_file`` key must raise, never be dropped (gmsh-official-mesher S7a, S9; T11 verifier finding).

T11 removed the excel_file forwarding from the runner and the BEM-FSI CLI; without a guard a config carrying both a v2
``yaml_file`` and an old ``excel_file`` would build silently from the yaml.
"""
from pathlib import Path

import pytest

from aeroelast.cli import run_bem_fsi
from aeroelast.solvers.fsi import runner as R

REPO = Path(__file__).resolve().parents[3]
V2 = str(REPO / "tests/IEA15MW/reference/upstream/IEA-15-240-RWT.yaml")
XLSX_CFG = REPO / "tests/archive/numad_xlsx_configs/IEA15MW/ab_test/solid_xlsx.yaml"


def _leftover_params(params):
    for key in ("airfoil_dir", "n_samples", "airfoil_spacing", "span_grading"):
        params.pop(key, None)
    params["yaml_file"] = V2
    params["excel_file"] = "blade.xlsx"
    return params


def test_runner_rejects_excel_file_next_to_v2_yaml():
    runner = R.FSIRunner(str(XLSX_CFG))
    _leftover_params(runner.config.mesh.generator.params)
    with pytest.raises(ValueError, match="excel_file"):
        runner._generate_mesh()


def test_bem_fsi_cli_rejects_excel_file_next_to_v2_yaml(tmp_path):
    cfg = {"mesh": {"generator": {"type": "GmshBladeMesh", "params": _leftover_params({})}}}
    with pytest.raises(ValueError, match="excel_file"):
        run_bem_fsi._build_mesh(cfg, tmp_path / "case.yaml")
