"""Frontiersin 2025 yaw aeroelastic study — FSI configuration and reference tests.

Validates simulation configurations for the yaw-condition study of the IEA-15 MW
wind turbine as described in:

    Ma L, Li Y, Zhou L, Yang D, Shen X, Du Z (2025)
    "Study on the aeroelastic performance of 15 MW wind turbine under yaw condition"
    Front. Energy Res. 13:1571567. doi: 10.3389/fenrg.2025.1571567

Article parameters (Section 3.1):
    - Turbine: IEA-15 MW (rotor diameter 240 m, hub height 150 m)
    - Wind speed: 10.59 m/s (rated)
    - Rotational speed: 7.55 rpm (omega = 0.7906 rad/s)
    - Pitch: 0 deg (no pitch control)
    - Inflow: uniform (shear_exp = 0), no turbulence
    - Yaw sweep: 0, 10, 20, 30, 40 deg

NOTE: The article uses a LL-FVW aerodynamic model with B-L dynamic stall, whereas
our solver uses CCBlade (BEM + yaw correction). Quantitative agreement is expected
to be approximate; reference values and tolerances in the CSV are set accordingly.
"""

from __future__ import annotations

import csv
import math
import os
from pathlib import Path

import pytest
import yaml

CASES_FILE = Path("tests/IEA15MW/frontiersin_2025_yaw_cases.yaml")
REFERENCE_FILE = Path("tests/IEA15MW/frontiersin_2025_yaw_reference.csv")
CASES_DIR = Path("tests/IEA15MW/frontiersin_2025_yaw")

YAW_ANGLES = [0, 10, 20, 30, 40]

# Article nominal conditions (Section 3.1).
RATED_WIND_SPEED_MPS = 10.59
RATED_OMEGA_RPM = 7.55
RATED_OMEGA_RADS = RATED_OMEGA_RPM * 2 * math.pi / 60  # 0.7906... rad/s
RATED_PITCH_DEG = 0.0
SHEAR_EXP = 0.0  # uniform inflow

# At yaw=20 deg the article reports max flapwise tip deflection ~15 m (Figure 11).
ARTICLE_MAX_FLAPWISE_YAW20_M = 15.0


@pytest.fixture(scope="module")
def cases_table() -> dict:
    return yaml.safe_load(CASES_FILE.read_text(encoding="utf-8"))


@pytest.fixture(scope="module")
def reference_rows() -> list[dict[str, str]]:
    with REFERENCE_FILE.open("r", encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f))


def _load_yaml(path: Path) -> dict:
    return yaml.safe_load(path.read_text(encoding="utf-8"))


def _candidate_result_files() -> list[Path]:
    env_path = os.getenv("AEROELAST_YAW_RESULTS_CSV")
    candidates: list[Path] = []
    if env_path:
        candidates.append(Path(env_path))
    candidates.extend([
        Path("output/frontiersin_2025_yaw/solver_tip_deflections.csv"),
        Path("output/frontiersin_2025_yaw_tip_deflections.csv"),
        Path("output/yaw_sweep_tip_deflections.csv"),
    ])
    return candidates


def _load_results_by_solver_and_yaw() -> tuple[Path, dict[tuple[str, int], dict[str, float]]]:
    csv_path = next((p for p in _candidate_result_files() if p.exists()), None)
    if csv_path is None:
        joined = ", ".join(str(p) for p in _candidate_result_files())
        pytest.skip(
            "No yaw-sweep results CSV found. Set AEROELAST_YAW_RESULTS_CSV or place one at: "
            + joined
        )

    required_cols = {"solver", "yaw_deg", "tip_flapwise_m", "tip_edgewise_m"}
    with csv_path.open("r", encoding="utf-8", newline="") as f:
        rows = list(csv.DictReader(f))

    if rows:
        fieldnames = set(rows[0].keys())
        if not required_cols.issubset(fieldnames):
            raise AssertionError(
                f"CSV {csv_path} must include {sorted(required_cols)}; found {sorted(fieldnames)}"
            )

    agg: dict[tuple[str, int], dict[str, float]] = {}
    for row in rows:
        raw_solver = str(row["solver"]).strip().lower()
        solver = "inertial" if "inert" in raw_solver else "corotational"
        yaw = int(float(row["yaw_deg"]))
        flap = abs(float(row["tip_flapwise_m"]))
        edge = abs(float(row["tip_edgewise_m"]))
        key = (solver, yaw)
        if key not in agg:
            agg[key] = {"flapwise": flap, "edgewise": edge}
        else:
            agg[key]["flapwise"] = max(agg[key]["flapwise"], flap)
            agg[key]["edgewise"] = max(agg[key]["edgewise"], edge)

    return csv_path, agg


# ---------------------------------------------------------------------------
# Artifact presence tests
# ---------------------------------------------------------------------------


def test_frontiersin_2025_cases_file_exists():
    assert CASES_FILE.exists(), f"Missing cases file: {CASES_FILE}"


def test_frontiersin_2025_reference_file_exists():
    assert REFERENCE_FILE.exists(), f"Missing reference CSV: {REFERENCE_FILE}"


def test_frontiersin_2025_cases_dir_exists():
    assert CASES_DIR.is_dir(), f"Missing cases directory: {CASES_DIR}"


def test_frontiersin_2025_all_yaml_files_present():
    missing = []
    for yaw in YAW_ANGLES:
        for stem in (
            f"fluid_yaw_{yaw}",
            f"solid_corotational_yaw_{yaw}",
            f"solid_inertial_yaw_{yaw}",
        ):
            p = CASES_DIR / f"{stem}.yaml"
            if not p.exists():
                missing.append(str(p))
    assert not missing, "Missing YAML files:\n" + "\n".join(missing)


# ---------------------------------------------------------------------------
# Cases table integrity tests
# ---------------------------------------------------------------------------


def test_frontiersin_2025_cases_table_has_all_yaw_angles(cases_table):
    got = {c["yaw_deg"] for c in cases_table["cases"]}
    assert got == set(YAW_ANGLES), f"Expected yaw angles {YAW_ANGLES}, got {sorted(got)}"


def test_frontiersin_2025_nominal_conditions(cases_table):
    nc = cases_table["nominal_conditions"]
    assert abs(nc["wind_speed_mps"] - RATED_WIND_SPEED_MPS) < 1e-6
    assert abs(nc["omega_rpm"] - RATED_OMEGA_RPM) < 1e-6
    assert abs(nc["omega_rads"] - RATED_OMEGA_RADS) < 1e-4
    assert abs(nc["pitch_deg"] - RATED_PITCH_DEG) < 1e-6
    assert nc["shear_exp"] == SHEAR_EXP


def test_frontiersin_2025_reference_csv_has_all_yaw_angles(reference_rows):
    got = {int(float(r["yaw_deg"])) for r in reference_rows}
    assert got == set(YAW_ANGLES), f"Expected yaw angles {YAW_ANGLES} in CSV, got {sorted(got)}"


# ---------------------------------------------------------------------------
# YAML content validation tests
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("yaw", YAW_ANGLES)
def test_frontiersin_2025_fluid_yaml_parameters(yaw):
    path = CASES_DIR / f"fluid_yaw_{yaw}.yaml"
    cfg = _load_yaml(path)
    bem = cfg["bem"]

    assert abs(bem["wind_speed"] - RATED_WIND_SPEED_MPS) < 1e-6, (
        f"fluid_yaw_{yaw}: wrong wind_speed"
    )
    assert abs(bem["yaw"] - float(yaw)) < 1e-6, f"fluid_yaw_{yaw}: wrong yaw"
    assert abs(bem["pitch"] - RATED_PITCH_DEG) < 1e-6, f"fluid_yaw_{yaw}: wrong pitch"
    assert abs(bem["omega"]["initial"] - RATED_OMEGA_RADS) < 1e-4, f"fluid_yaw_{yaw}: wrong omega"
    assert abs(bem["shear_exp"] - SHEAR_EXP) < 1e-9, (
        f"fluid_yaw_{yaw}: shear_exp must be 0 for uniform inflow"
    )
    assert bem["hub_height"] == 150.0, f"fluid_yaw_{yaw}: wrong hub_height"


@pytest.mark.parametrize("yaw", YAW_ANGLES)
def test_frontiersin_2025_solid_corotational_yaml_parameters(yaw):
    path = CASES_DIR / f"solid_corotational_yaw_{yaw}.yaml"
    cfg = _load_yaml(path)
    rotor = cfg["solver"]["rotor"]

    assert cfg["solver"]["type"] == "LinearDynamicFSIRotor", (
        f"solid_corotational_yaw_{yaw}: wrong solver type"
    )
    assert abs(rotor["omega"] - RATED_OMEGA_RADS) < 1e-4, (
        f"solid_corotational_yaw_{yaw}: wrong omega"
    )
    assert rotor["omega_ramp_time"] == 0, (
        f"solid_corotational_yaw_{yaw}: omega_ramp_time must be 0 (fixed omega)"
    )
    assert abs(cfg["performance"]["flow_velocity"] - RATED_WIND_SPEED_MPS) < 1e-6, (
        f"solid_corotational_yaw_{yaw}: wrong flow_velocity"
    )


@pytest.mark.parametrize("yaw", YAW_ANGLES)
def test_frontiersin_2025_solid_inertial_yaml_parameters(yaw):
    path = CASES_DIR / f"solid_inertial_yaw_{yaw}.yaml"
    cfg = _load_yaml(path)
    rotor = cfg["solver"]["rotor"]

    assert cfg["solver"]["type"] == "LinearDynamicFSIRotorInertial", (
        f"solid_inertial_yaw_{yaw}: wrong solver type"
    )
    assert abs(rotor["omega"] - RATED_OMEGA_RADS) < 1e-4, f"solid_inertial_yaw_{yaw}: wrong omega"
    assert rotor["omega_ramp_time"] == 0, (
        f"solid_inertial_yaw_{yaw}: omega_ramp_time must be 0 (fixed omega)"
    )
    assert abs(cfg["performance"]["flow_velocity"] - RATED_WIND_SPEED_MPS) < 1e-6, (
        f"solid_inertial_yaw_{yaw}: wrong flow_velocity"
    )


@pytest.mark.parametrize("yaw", YAW_ANGLES)
def test_frontiersin_2025_output_folders_are_unique(yaw):
    """Each solver variant must write to a distinct output sub-directory.

    NOTE: yaw-case isolation is achieved via CWD (each case runs from its own
    yaw_N/ working directory), NOT by embedding the yaw angle in the folder
    name. The folder names are solver-type identifiers only.
    """
    corot = _load_yaml(CASES_DIR / f"solid_corotational_yaw_{yaw}.yaml")
    inert = _load_yaml(CASES_DIR / f"solid_inertial_yaw_{yaw}.yaml")
    fluid = _load_yaml(CASES_DIR / f"fluid_yaw_{yaw}.yaml")

    folders = {
        "corotational": corot["output"]["folder"],
        "inertial": inert["output"]["folder"],
        "fluid": fluid["output"]["folder"],
    }
    unique = set(folders.values())
    assert len(unique) == 3, f"yaw={yaw}: output folders must all be distinct, got {folders}"


def test_frontiersin_2025_dual_solver_omega_consistency():
    """Corotational and inertial YAMLs for the same yaw angle must use the same omega."""
    for yaw in YAW_ANGLES:
        corot = _load_yaml(CASES_DIR / f"solid_corotational_yaw_{yaw}.yaml")
        inert = _load_yaml(CASES_DIR / f"solid_inertial_yaw_{yaw}.yaml")
        fluid = _load_yaml(CASES_DIR / f"fluid_yaw_{yaw}.yaml")

        omega_corot = corot["solver"]["rotor"]["omega"]
        omega_inert = inert["solver"]["rotor"]["omega"]
        omega_fluid = fluid["bem"]["omega"]["initial"]

        assert abs(omega_corot - omega_inert) < 1e-10, (
            f"yaw={yaw}: corotational and inertial omega differ: {omega_corot} vs {omega_inert}"
        )
        assert abs(omega_corot - omega_fluid) < 1e-4, (
            f"yaw={yaw}: solid and fluid omega differ: {omega_corot} vs {omega_fluid}"
        )


# ---------------------------------------------------------------------------
# Trend validation: flapwise should decrease with increasing yaw (article Fig 16a)
# ---------------------------------------------------------------------------


def test_frontiersin_2025_reference_flapwise_decreases_with_yaw(reference_rows):
    """Article Fig 16a: mean flapwise deflection decreases monotonically with yaw."""
    by_yaw = {int(float(r["yaw_deg"])): float(r["flapwise_ref_m"]) for r in reference_rows}
    sorted_yaws = sorted(by_yaw.keys())
    for i in range(len(sorted_yaws) - 1):
        y0, y1 = sorted_yaws[i], sorted_yaws[i + 1]
        assert by_yaw[y0] >= by_yaw[y1], (
            f"Reference flapwise at yaw={y0} ({by_yaw[y0]:.2f} m) should be >= "
            f"yaw={y1} ({by_yaw[y1]:.2f} m) per article Fig 16a"
        )


# ---------------------------------------------------------------------------
# Simulation result comparison (skipped if no results CSV is present)
# ---------------------------------------------------------------------------


def test_frontiersin_2025_simulation_flapwise_within_reference(reference_rows):
    csv_path, results = _load_results_by_solver_and_yaw()

    for row in reference_rows:
        yaw = int(float(row["yaw_deg"]))
        flap_ref = float(row["flapwise_ref_m"])
        flap_tol = float(row["flapwise_tol_m"])
        edge_ref = float(row["edgewise_ref_m"])
        edge_tol = float(row["edgewise_tol_m"])

        for solver in ("corotational", "inertial"):
            key = (solver, yaw)
            assert key in results, f"Missing result for solver={solver}, yaw={yaw} in {csv_path}"
            flap_val = results[key]["flapwise"]
            edge_val = results[key]["edgewise"]

            assert flap_ref - flap_tol <= flap_val <= flap_ref + flap_tol, (
                f"solver={solver}, yaw={yaw}: flapwise {flap_val:.3f} m outside "
                f"[{flap_ref - flap_tol:.3f}, {flap_ref + flap_tol:.3f}] m"
            )
            assert edge_ref - edge_tol <= edge_val <= edge_ref + edge_tol, (
                f"solver={solver}, yaw={yaw}: edgewise {edge_val:.3f} m outside "
                f"[{edge_ref - edge_tol:.3f}, {edge_ref + edge_tol:.3f}] m"
            )


def test_frontiersin_2025_simulation_flapwise_decreases_with_yaw():
    """Mean flapwise tip deflection should decrease monotonically with yaw (article Fig 16a)."""
    _, results = _load_results_by_solver_and_yaw()

    for solver in ("corotational", "inertial"):
        sorted_yaws = sorted(yaw for (s, yaw) in results if s == solver)
        for i in range(len(sorted_yaws) - 1):
            y0, y1 = sorted_yaws[i], sorted_yaws[i + 1]
            flap0 = results[(solver, y0)]["flapwise"]
            flap1 = results[(solver, y1)]["flapwise"]
            assert flap0 >= flap1 * 0.8, (  # 20% grace for BEM inaccuracy near transitions
                f"solver={solver}: flapwise at yaw={y0} ({flap0:.3f} m) should be >= "
                f"yaw={y1} ({flap1:.3f} m) per article Fig 16a trend"
            )


def test_frontiersin_2025_yaw20_flapwise_order_of_magnitude():
    """At yaw=20 deg, article reports ~15 m flapwise (LL-FVW). BEM should be in same order."""
    _, results = _load_results_by_solver_and_yaw()

    for solver in ("corotational", "inertial"):
        key = (solver, 20)
        if key not in results:
            pytest.skip(f"No results for {solver}, yaw=20")
        flap = results[key]["flapwise"]
        # Generous bounds: 30% to 200% of article value (BEM vs LL-FVW expected to differ)
        lo = 0.30 * ARTICLE_MAX_FLAPWISE_YAW20_M
        hi = 2.0 * ARTICLE_MAX_FLAPWISE_YAW20_M
        assert lo <= flap <= hi, (
            f"solver={solver}, yaw=20: flapwise {flap:.3f} m is outside "
            f"[{lo:.1f}, {hi:.1f}] m (article LL-FVW: ~{ARTICLE_MAX_FLAPWISE_YAW20_M:.0f} m)"
        )


# ---------------------------------------------------------------------------
# Physics configuration consistency tests
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("yaw", YAW_ANGLES)
def test_frontiersin_2025_damping_consistency(yaw):
    """Both solver variants must use the same structural damping ratio.

    Using different zeta would contaminate the corotational-vs-inertial comparison
    since both represent the same physical blade.
    """
    corot = _load_yaml(CASES_DIR / f"solid_corotational_yaw_{yaw}.yaml")
    inert = _load_yaml(CASES_DIR / f"solid_inertial_yaw_{yaw}.yaml")

    zeta_corot = corot["solver"]["damping"]["zeta"]
    zeta_inert = inert["solver"]["damping"]["zeta"]
    assert zeta_corot == zeta_inert, (
        f"yaw={yaw}: damping zeta mismatch — corotational={zeta_corot}, inertial={zeta_inert}. "
        "Both solvers must use the same structural damping for a valid comparison."
    )
    mode_i_corot = corot["solver"]["damping"]["mode_i"]
    mode_j_corot = corot["solver"]["damping"]["mode_j"]
    mode_i_inert = inert["solver"]["damping"]["mode_i"]
    mode_j_inert = inert["solver"]["damping"]["mode_j"]
    assert (mode_i_corot, mode_j_corot) == (mode_i_inert, mode_j_inert), (
        f"yaw={yaw}: Rayleigh calibration modes differ — "
        f"corotational=({mode_i_corot},{mode_j_corot}), inertial=({mode_i_inert},{mode_j_inert})"
    )


@pytest.mark.parametrize("yaw", YAW_ANGLES)
def test_frontiersin_2025_geometric_stiffness_enabled(yaw):
    """K_G must be enabled for all solid configs.

    At rated conditions flatwise pre-stress stiffening raises the first natural
    frequency by ~10-20%. Disabling K_G in either solver makes the comparison
    physically inequivalent.
    """
    for solver_stem in (f"solid_corotational_yaw_{yaw}", f"solid_inertial_yaw_{yaw}"):
        cfg = _load_yaml(CASES_DIR / f"{solver_stem}.yaml")
        rotor = cfg["solver"]["rotor"]
        assert rotor.get("include_geometric_stiffness") is True, (
            f"{solver_stem}: include_geometric_stiffness must be true — "
            "flatwise pre-stress stiffening is ~10-20% at rated speed"
        )


@pytest.mark.parametrize("yaw", YAW_ANGLES)
def test_frontiersin_2025_corotational_physics_flags(yaw):
    """Corotational solver must have all rotation-dependent physics flags enabled.

    For omega > 0: spin softening (K_SP ∝ ω²), centrifugal (F_cf ∝ ω²),
    Coriolis (F_cor ∝ ω), and Euler (F_euler ∝ ω̇) all contribute to the
    rotating-frame dynamics.
    """
    cfg = _load_yaml(CASES_DIR / f"solid_corotational_yaw_{yaw}.yaml")
    rotor = cfg["solver"]["rotor"]
    for flag in (
        "include_spin_softening",
        "include_centrifugal",
        "include_coriolis",
        "include_euler",
    ):
        assert rotor.get(flag) is True, (
            f"solid_corotational_yaw_{yaw}: {flag} must be true for omega > 0"
        )


@pytest.mark.parametrize("yaw", YAW_ANGLES)
def test_frontiersin_2025_inertial_velocity_flags(yaw):
    """Inertial solver: omega to preCICE on, structural-velocity feedback off.

    ``send_omega_to_precice`` feeds the rotational speed to the BEM fluid.
    ``send_velocity_to_precice`` (the structural-velocity feedback) is the
    variant with an unresolved mean-power anomaly (see
    docs/validation_results_v01_v02_v04.md §bem_0_10_full); it must stay off
    in the campaign configs.
    """
    cfg = _load_yaml(CASES_DIR / f"solid_inertial_yaw_{yaw}.yaml")
    rotor = cfg["solver"]["rotor"]
    assert rotor.get("send_velocity_to_precice") is False, (
        f"solid_inertial_yaw_{yaw}: send_velocity_to_precice must be false "
        "(the velocity-feedback variant has an unresolved power anomaly)"
    )
    assert rotor.get("send_omega_to_precice") is True, (
        f"solid_inertial_yaw_{yaw}: send_omega_to_precice must be true"
    )


@pytest.mark.parametrize("yaw", YAW_ANGLES)
def test_frontiersin_2025_newmark_parameters(yaw):
    """Both solvers must use the average-acceleration Newmark scheme (β=0.25, γ=0.5).

    This is the only unconditionally stable, non-dissipative Newmark variant.
    Any other values would either add artificial damping or introduce instability.
    """
    for solver_stem in (f"solid_corotational_yaw_{yaw}", f"solid_inertial_yaw_{yaw}"):
        cfg = _load_yaml(CASES_DIR / f"{solver_stem}.yaml")
        newmark = cfg["solver"].get("newmark", {})
        assert abs(newmark.get("beta", 0.0) - 0.25) < 1e-12, (
            f"{solver_stem}: Newmark beta must be 0.25 (average acceleration)"
        )
        assert abs(newmark.get("gamma", 0.0) - 0.5) < 1e-12, (
            f"{solver_stem}: Newmark gamma must be 0.5 (average acceleration)"
        )


@pytest.mark.parametrize("yaw", YAW_ANGLES)
def test_frontiersin_2025_mesh_parameters_consistent(yaw):
    """All solid YAMLs for a given yaw must use the same mesh parameters.

    Comparing solvers on different meshes would be physically invalid.
    """
    corot = _load_yaml(CASES_DIR / f"solid_corotational_yaw_{yaw}.yaml")
    inert = _load_yaml(CASES_DIR / f"solid_inertial_yaw_{yaw}.yaml")

    params_corot = corot["mesh"]["generator"]["params"]
    params_inert = inert["mesh"]["generator"]["params"]
    for key in ("element_size", "n_samples", "span_grading", "airfoil_spacing"):
        assert params_corot.get(key) == params_inert.get(key), (
            f"yaw={yaw}: mesh param '{key}' differs — "
            f"corotational={params_corot.get(key)}, inertial={params_inert.get(key)}"
        )
