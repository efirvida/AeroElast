"""Chapter 6 DLC prepared cases — configuration tests.

Validates the YAML templates and run matrix for the Chapter 6 Design Load Cases
(DLC) study of the IEA-15 MW wind turbine:

    DLC 1.1  Normal power production — NTM turbulence
    DLC 6.1  Parked — EWM50 extreme wind (omega=0, pitch=90°)
    DLC 6.3  Parked with fault — EWM1 extreme wind (omega=0, pitch=90°)

Templates are parametric: wind speed, yaw, and seed are injected by the launcher
from ch6_run_matrix.csv at runtime.
"""

from __future__ import annotations

from pathlib import Path

import pytest
import yaml

# ── Paths ───────────────────────────────────────────────────────────────────────
CASES_DIR = Path("tests/IEA15MW/ch6_prepared_cases")
RUN_MATRIX = CASES_DIR / "ch6_run_matrix.csv"

DLCS = ["1.1", "6.1", "6.3"]
SOLVERS = ["corotational", "inertial"]

# IEA-15 MW asset paths expected inside every YAML (relative to YAML's parent dir)
# The canonical blade input since 2026-09-08 is the official WindIO yaml
# (tests/IEA-15-240-RWT.yaml); the UTD excel input is retired from new work.
EXPECTED_BLADE_YAML = "../../IEA-15-240-RWT.yaml"
EXPECTED_AIRFOIL_DIR = "../blade_definition/airfoils"
EXPECTED_PRECICE_CONFIG = "../precice-config.xml"

# Expected solver types per template
SOLID_SOLVER_TYPES = {
    "corotational": "LinearDynamicFSIRotor",
    "inertial": "LinearDynamicFSIRotorInertial",
}

# DLC 1.1 rated omega (rad/s); parked DLCs must be 0.
OMEGA_RATED = 0.7872360607268416
OMEGA_PARKED = 0.0


# ── Helpers ─────────────────────────────────────────────────────────────────────
def _yaml(filename: str) -> dict:
    """Load a YAML from CASES_DIR. Pass the filename exactly as it exists on disk."""
    return yaml.safe_load((CASES_DIR / filename).read_text(encoding="utf-8"))


def _fluid_yaml_name(dlc: str) -> str:
    return f"fluid_dlc_{dlc}.yaml"


def _solid_yaml_name(solver: str, dlc: str) -> str:
    return f"solid_{solver}_dlc_{dlc}.yaml"


def _read_matrix() -> list[dict[str, str]]:
    import csv

    with RUN_MATRIX.open("r", encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f))


# ── Directory and file existence ────────────────────────────────────────────────
def test_ch6_cases_dir_exists():
    assert CASES_DIR.is_dir(), f"Cases directory not found: {CASES_DIR}"


def test_ch6_run_matrix_exists():
    assert RUN_MATRIX.is_file(), f"Run matrix not found: {RUN_MATRIX}"


@pytest.mark.parametrize("dlc", DLCS)
def test_ch6_fluid_yaml_exists(dlc: str):
    f = CASES_DIR / _fluid_yaml_name(dlc)
    assert f.is_file(), f"Missing fluid YAML: {f}"


@pytest.mark.parametrize("solver,dlc", [(s, d) for s in SOLVERS for d in DLCS])
def test_ch6_solid_yaml_exists(solver: str, dlc: str):
    f = CASES_DIR / _solid_yaml_name(solver, dlc)
    assert f.is_file(), f"Missing solid YAML: {f}"


# ── Run matrix structure ─────────────────────────────────────────────────────────
def test_ch6_run_matrix_columns():
    rows = _read_matrix()
    assert rows, "Run matrix is empty"
    required = {"solver", "dlc", "wind_model", "wind_speed_mps", "seed", "case_id"}
    assert required.issubset(rows[0].keys()), f"Missing columns: {required - rows[0].keys()}"


def test_ch6_run_matrix_has_both_solvers():
    rows = _read_matrix()
    solvers = {r["solver"] for r in rows}
    assert "corotational" in solvers
    assert "inertial" in solvers


def test_ch6_run_matrix_has_all_dlcs():
    rows = _read_matrix()
    dlcs = {r["dlc"] for r in rows}
    for dlc in DLCS:
        assert dlc in dlcs, f"DLC {dlc} missing from run matrix"


def test_ch6_run_matrix_case_ids_are_unique():
    rows = _read_matrix()
    ids = [r["case_id"] for r in rows]
    assert len(ids) == len(set(ids)), "Duplicate case_id entries in run matrix"


def test_ch6_run_matrix_row_count():
    """192 data rows (96 corotational + 96 inertial) — matches the README."""
    rows = _read_matrix()
    assert len(rows) == 192, f"Expected 192 rows, got {len(rows)}"


# ── Fluid YAMLs ──────────────────────────────────────────────────────────────────
@pytest.mark.parametrize("dlc", DLCS)
def test_ch6_fluid_yaml_participant(dlc: str):
    cfg = _yaml(_fluid_yaml_name(dlc))
    assert cfg["participant"] == "Fluid"


@pytest.mark.parametrize("dlc", DLCS)
def test_ch6_fluid_yaml_asset_paths(dlc: str):
    cfg = _yaml(_fluid_yaml_name(dlc))
    gen = cfg["mesh"]["generator"]["params"]
    assert gen["yaml_file"] == EXPECTED_BLADE_YAML, (
        f"DLC {dlc} fluid: yaml_file should be '{EXPECTED_BLADE_YAML}'"
    )
    assert gen["airfoil_dir"] == EXPECTED_AIRFOIL_DIR, (
        f"DLC {dlc} fluid: airfoil_dir should be '{EXPECTED_AIRFOIL_DIR}'"
    )


@pytest.mark.parametrize("dlc", DLCS)
def test_ch6_fluid_yaml_precice_config(dlc: str):
    cfg = _yaml(_fluid_yaml_name(dlc))
    assert cfg["config_file"] == EXPECTED_PRECICE_CONFIG


@pytest.mark.parametrize("dlc", DLCS)
def test_ch6_fluid_yaml_omega_mesh_configured(dlc: str):
    """Fluid must receive omega via GlobalFluidMesh/AngularVelocity."""
    cfg = _yaml(_fluid_yaml_name(dlc))
    omega_cfg = cfg["bem"]["omega"]
    assert omega_cfg.get("mesh") == "GlobalFluidMesh"
    assert omega_cfg.get("data") == "AngularVelocity"


def test_ch6_fluid_dlc11_operational():
    """DLC 1.1 is power-production: pitch=0, omega>0."""
    cfg = _yaml("fluid_dlc_1.1.yaml")
    assert cfg["bem"]["pitch"] == 0.0
    assert cfg["bem"]["omega"]["initial"] > 0.0


@pytest.mark.parametrize("dlc", ["6.1", "6.3"])
def test_ch6_fluid_parked_conditions(dlc: str):
    """DLC 6.1 and 6.3 are parked: pitch=90, omega=0."""
    cfg = _yaml(_fluid_yaml_name(dlc))
    assert cfg["bem"]["pitch"] == 90.0, f"DLC {dlc}: pitch must be 90° (parked)"
    assert cfg["bem"]["omega"]["initial"] == 0.0, f"DLC {dlc}: omega must be 0 (parked)"


def test_ch6_fluid_output_folders_are_unique():
    """Each DLC fluid template must write to a distinct output folder."""
    folders = {}
    for dlc in DLCS:
        cfg = _yaml(_fluid_yaml_name(dlc))
        folders[dlc] = cfg["output"]["folder"]
    assert len(set(folders.values())) == len(DLCS), f"Duplicate fluid output folders: {folders}"


# ── Solid YAMLs ──────────────────────────────────────────────────────────────────
@pytest.mark.parametrize("solver,dlc", [(s, d) for s in SOLVERS for d in DLCS])
def test_ch6_solid_yaml_asset_paths(solver: str, dlc: str):
    cfg = _yaml(_solid_yaml_name(solver, dlc))
    gen = cfg["mesh"]["generator"]["params"]
    assert gen["yaml_file"] == EXPECTED_BLADE_YAML, (
        f"{solver}/DLC {dlc}: yaml_file should be '{EXPECTED_BLADE_YAML}'"
    )
    assert gen["airfoil_dir"] == EXPECTED_AIRFOIL_DIR, (
        f"{solver}/DLC {dlc}: airfoil_dir should be '{EXPECTED_AIRFOIL_DIR}'"
    )


@pytest.mark.parametrize("solver,dlc", [(s, d) for s in SOLVERS for d in DLCS])
def test_ch6_solid_yaml_precice_config(solver: str, dlc: str):
    cfg = _yaml(_solid_yaml_name(solver, dlc))
    assert cfg["coupling"]["config_file"] == EXPECTED_PRECICE_CONFIG


@pytest.mark.parametrize("solver,dlc", [(s, d) for s in SOLVERS for d in DLCS])
def test_ch6_solid_yaml_solver_type(solver: str, dlc: str):
    cfg = _yaml(_solid_yaml_name(solver, dlc))
    expected = SOLID_SOLVER_TYPES[solver]
    actual = cfg["solver"]["type"]
    assert actual == expected, (
        f"{solver}/DLC {dlc}: solver.type should be '{expected}', got '{actual}'"
    )


@pytest.mark.parametrize("solver", SOLVERS)
def test_ch6_solid_dlc11_operational_omega(solver: str):
    """DLC 1.1 must use rated omega."""
    cfg = _yaml(_solid_yaml_name(solver, "1.1"))
    omega = cfg["solver"]["rotor"]["omega"]
    assert abs(omega - OMEGA_RATED) < 1e-9, (
        f"{solver}/DLC 1.1: rotor.omega should be {OMEGA_RATED}, got {omega}"
    )


@pytest.mark.parametrize("solver,dlc", [(s, d) for s in SOLVERS for d in ["6.1", "6.3"]])
def test_ch6_solid_parked_omega_is_zero(solver: str, dlc: str):
    """Parked DLCs must have omega=0."""
    cfg = _yaml(_solid_yaml_name(solver, dlc))
    omega = cfg["solver"]["rotor"]["omega"]
    assert omega == OMEGA_PARKED, f"{solver}/DLC {dlc}: rotor.omega must be 0 (parked), got {omega}"


@pytest.mark.parametrize("solver,dlc", [(s, d) for s in SOLVERS for d in DLCS])
def test_ch6_solid_yaml_send_omega_to_precice(solver: str, dlc: str):
    """All templates must forward omega to the fluid via AngularVelocity."""
    cfg = _yaml(_solid_yaml_name(solver, dlc))
    assert cfg["solver"]["rotor"]["send_omega_to_precice"] is True, (
        f"{solver}/DLC {dlc}: send_omega_to_precice must be true"
    )


def test_ch6_solid_inertial_dlc11_send_velocity():
    """Inertial solver DLC 1.1 must also send velocity (needed by inertial coupling)."""
    cfg = _yaml(_solid_yaml_name("inertial", "1.1"))
    assert cfg["solver"]["rotor"]["send_velocity_to_precice"] is True


@pytest.mark.parametrize("solver", SOLVERS)
def test_ch6_solid_output_folders_are_unique_across_dlcs(solver: str):
    """Each DLC solid template must write to a distinct output folder."""
    folders = {}
    for dlc in DLCS:
        cfg = _yaml(_solid_yaml_name(solver, dlc))
        folders[dlc] = cfg["output"]["folder"]
    assert len(set(folders.values())) == len(DLCS), (
        f"{solver}: duplicate solid output folders: {folders}"
    )


def test_ch6_solid_corotational_inertial_output_folders_differ():
    """corotational and inertial must not share output folder names per DLC."""
    for dlc in DLCS:
        coro = _yaml(_solid_yaml_name("corotational", dlc))["output"]["folder"]
        inert = _yaml(_solid_yaml_name("inertial", dlc))["output"]["folder"]
        assert coro != inert, (
            f"DLC {dlc}: corotational and inertial share the same output folder '{coro}'"
        )
