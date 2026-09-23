"""V-06: DLC FSI load-assessment cases and Chapter 6 comparisons.

This test suite does two things:
1. Builds explicit FSI case definitions for Chapter 6 DLCs using local
   simulation references (corotational + inertial).
2. Compares per-DLC flapwise/edgewise tip deflection outputs against Chapter 6
   reference envelopes.

Reference documents:
- tests/IEA15MW/75698.pdf (NREL/TP-5000-75698, Chapter 6)
- tests/IEA15MW/ch6_dlc_cases.yaml (Table 6-1 transcription)
- tests/IEA15MW/ch6_dlc_tip_deflection_reference.csv (Figure 6-2 envelopes)
"""

from __future__ import annotations

import csv
import os
from copy import deepcopy
from pathlib import Path

import pytest
import yaml

CH6_PDF = Path("tests/IEA15MW/validation papers/75698.pdf")
CH6_DLC_CASES_FILE = Path("tests/IEA15MW/ch6_dlc_cases.yaml")
CH6_DLC_REFERENCE_FILE = Path("tests/IEA15MW/ch6_dlc_tip_deflection_reference.csv")

# User-requested simulation references.
COROT_SOLID_REF = Path(
    "/scratch/leahk/eduardo.donestevez/simulations/bem_0_10/solid/simulation_xls.yaml"
)
INERTIAL_SOLID_REF_CANDIDATES = [
    Path("/scratch/leahk/eduardo.donestevez/simulations/bem_0_10_full/solid/simulation_xls.yaml"),
    Path("/scratch/leahk/eduardo.donestevez/simulations/bem_0_10/solid/simulation_xls.yaml"),
]
FLUID_REF = Path("/scratch/leahk/eduardo.donestevez/simulations/bem_0_10/fluid/simulation_xls.yaml")

CH6_MAX_OUT_OF_PLANE_TIP_DEFLECTION_M = 22.8


@pytest.fixture(scope="module")
def ch6_dlc_table() -> dict:
    return yaml.safe_load(CH6_DLC_CASES_FILE.read_text(encoding="utf-8"))


@pytest.fixture(scope="module")
def ch6_dlc_reference_rows() -> list[dict[str, str]]:
    with CH6_DLC_REFERENCE_FILE.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def _load_yaml(path: Path) -> dict:
    return yaml.safe_load(path.read_text(encoding="utf-8"))


def _resolve_inertial_solid_reference() -> Path:
    for candidate in INERTIAL_SOLID_REF_CANDIDATES:
        if candidate.exists():
            return candidate
    return INERTIAL_SOLID_REF_CANDIDATES[0]


def _build_solver_dlc_cases(
    base_solid_cfg: dict,
    base_fluid_cfg: dict,
    dlc_rows: list[dict],
    solver_tag: str,
    solver_type: str,
) -> list[dict]:
    cases: list[dict] = []
    for row in dlc_rows:
        solid_cfg = deepcopy(base_solid_cfg)
        fluid_cfg = deepcopy(base_fluid_cfg)

        solid_cfg.setdefault("solver", {})["type"] = solver_type
        solid_cfg.setdefault("dlc", {})
        solid_cfg["dlc"]["id"] = row["dlc"]
        solid_cfg["dlc"]["wind_condition"] = row["wind_condition"]
        solid_cfg["dlc"]["wind_speeds"] = row["wind_speeds"]
        solid_cfg["dlc"]["additional_settings"] = row["additional_settings"]
        solid_cfg["dlc"]["seeds"] = int(row["seeds"])
        solid_cfg["dlc"]["simulations"] = int(row["simulations"])

        # Keep BEM fluid setup and annotate with DLC metadata.
        fluid_cfg.setdefault("dlc", {})
        fluid_cfg["dlc"]["id"] = row["dlc"]
        fluid_cfg["dlc"]["wind_condition"] = row["wind_condition"]

        cases.append({
            "solver": solver_tag,
            "dlc": row["dlc"],
            "solid": solid_cfg,
            "fluid": fluid_cfg,
        })

    return cases


def _candidate_result_files() -> list[Path]:
    env_path = os.getenv("AEROELAST_DLC_RESULTS_CSV")
    candidates: list[Path] = []
    if env_path:
        candidates.append(Path(env_path))

    candidates.extend([
        Path("output/dlc/ch6_solver_dlc_tip_deflections.csv"),
        Path("output/dlc/solver_dlc_tip_deflections.csv"),
        Path("output/ch6_solver_dlc_tip_deflections.csv"),
    ])
    return candidates


def _canonical_solver_name(raw: str) -> str | None:
    val = raw.strip().lower()
    if "inert" in val:
        return "inertial"
    if "corot" in val or "rotor" in val:
        return "corotational"
    return None


def _load_solver_results_by_solver_and_dlc() -> tuple[
    Path, dict[tuple[str, str], dict[str, float]]
]:
    csv_path = next((p for p in _candidate_result_files() if p.exists()), None)
    if csv_path is None:
        joined = ", ".join(str(p) for p in _candidate_result_files())
        pytest.skip(
            "No solver DLC results CSV found. Set AEROELAST_DLC_RESULTS_CSV or place one at: "
            + joined
        )

    required_cols = {"solver", "dlc", "tip_flapwise_m", "tip_edgewise_m"}
    with csv_path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))

    fieldnames = set(rows[0].keys()) if rows else set()
    if not required_cols.issubset(fieldnames):
        raise AssertionError(
            f"CSV {csv_path} must include columns {sorted(required_cols)}; found {sorted(fieldnames)}"
        )

    agg: dict[tuple[str, str], dict[str, float]] = {}
    for row in rows:
        solver = _canonical_solver_name(str(row["solver"]))
        if solver is None:
            continue

        dlc = str(row["dlc"]).strip()
        flap = abs(float(row["tip_flapwise_m"]))
        edge = abs(float(row["tip_edgewise_m"]))
        key = (solver, dlc)

        if key not in agg:
            agg[key] = {"flapwise": flap, "edgewise": edge}
        else:
            agg[key]["flapwise"] = max(agg[key]["flapwise"], flap)
            agg[key]["edgewise"] = max(agg[key]["edgewise"], edge)

    return csv_path, agg


# ── artifact availability (data gaps are skipped with an explicit reason) ────
_HAS_DLC_REF = CH6_DLC_REFERENCE_FILE.exists()
_HAS_SIM_REFS = COROT_SOLID_REF.exists() and FLUID_REF.exists()
_RESULTS_FILE = next((p for p in _candidate_result_files() if p.exists()), None)
_HAS_RESULTS = _RESULTS_FILE is not None
_MISSING_DLC_REF = (
    f"Chapter 6 DLC reference values not transcribed at {CH6_DLC_REFERENCE_FILE}"
)
_MISSING_SIM_REFS = (
    "reference simulation configs missing (the simulations/ case layout was "
    "reorganised; expected e.g. " + str(COROT_SOLID_REF) + ")"
)
_MISSING_RESULTS = (
    "solver DLC tip-deflection results not produced yet (expected one of "
    + ", ".join(str(p) for p in _candidate_result_files()) + ")"
)


def test_ch6_reference_artifacts_exist():
    assert CH6_PDF.exists(), f"Missing Chapter 6 PDF at {CH6_PDF}"
    assert CH6_DLC_CASES_FILE.exists(), f"Missing DLC case table at {CH6_DLC_CASES_FILE}"
    if not _HAS_DLC_REF:
        pytest.skip(_MISSING_DLC_REF)


@pytest.mark.skipif(not _HAS_SIM_REFS, reason=_MISSING_SIM_REFS)
def test_reference_simulation_configs_exist_and_match_requested_sources():
    inertial_ref = _resolve_inertial_solid_reference()

    assert COROT_SOLID_REF.exists(), f"Missing corotational reference config: {COROT_SOLID_REF}"
    assert FLUID_REF.exists(), f"Missing fluid BEM reference config: {FLUID_REF}"
    assert inertial_ref.exists(), f"Missing inertial reference config: {inertial_ref}"


def test_ch6_table_6_1_dlc_settings(ch6_dlc_table):
    dlcs = ch6_dlc_table["dlcs"]
    got = {row["dlc"] for row in dlcs}
    expected = {"1.1", "6.1", "6.3"}
    assert got == expected

    row_11 = next(row for row in dlcs if row["dlc"] == "1.1")
    row_61 = next(row for row in dlcs if row["dlc"] == "6.1")
    assert row_11["wind_condition"] == "NTM"
    assert row_11["simulations"] == 72
    assert row_61["wind_condition"] == "EWM"
    assert row_61["additional_settings"] == "Yaw +/-8 deg"


@pytest.mark.skipif(not _HAS_SIM_REFS, reason=_MISSING_SIM_REFS)
def test_build_fsi_cases_from_reference_configs(ch6_dlc_table):
    inertial_ref = _resolve_inertial_solid_reference()

    corot_solid = _load_yaml(COROT_SOLID_REF)
    inertial_solid = _load_yaml(inertial_ref)
    fluid_cfg = _load_yaml(FLUID_REF)

    corot_cases = _build_solver_dlc_cases(
        base_solid_cfg=corot_solid,
        base_fluid_cfg=fluid_cfg,
        dlc_rows=ch6_dlc_table["dlcs"],
        solver_tag="corotational",
        solver_type="LinearDynamicFSIRotor",
    )
    inertial_cases = _build_solver_dlc_cases(
        base_solid_cfg=inertial_solid,
        base_fluid_cfg=fluid_cfg,
        dlc_rows=ch6_dlc_table["dlcs"],
        solver_tag="inertial",
        solver_type="LinearDynamicFSIRotorInertial",
    )

    assert len(corot_cases) == 3
    assert len(inertial_cases) == 3

    for case in corot_cases + inertial_cases:
        assert case["dlc"] in {"1.1", "6.1", "6.3"}
        assert "bem" in case["fluid"], "Fluid reference must use BEM for DLC FSI tests"
        assert case["fluid"].get("participant") == "Fluid"

    # Validate the requested dual-solver setup (same DLC set, different structural formulation).
    assert {c["dlc"] for c in corot_cases} == {c["dlc"] for c in inertial_cases}
    assert all(c["solid"]["solver"]["type"] == "LinearDynamicFSIRotor" for c in corot_cases)
    assert all(
        c["solid"]["solver"]["type"] == "LinearDynamicFSIRotorInertial" for c in inertial_cases
    )


@pytest.mark.skipif(
    not _HAS_DLC_REF or not _HAS_RESULTS,
    reason=_MISSING_DLC_REF + " / " + _MISSING_RESULTS,
)
def test_solver_dlc_results_against_ch6_reference_values(ch6_dlc_reference_rows):
    csv_path, results = _load_solver_results_by_solver_and_dlc()

    for solver in ("corotational", "inertial"):
        for ref in ch6_dlc_reference_rows:
            dlc = ref["dlc"]
            key = (solver, dlc)
            assert key in results, f"Missing result for solver={solver}, dlc={dlc} in {csv_path}"

            flap_ref = float(ref["flapwise_ref_m"])
            flap_tol = float(ref["flapwise_tol_m"])
            edge_ref = float(ref["edgewise_ref_m"])
            edge_tol = float(ref["edgewise_tol_m"])

            flap_val = results[key]["flapwise"]
            edge_val = results[key]["edgewise"]

            assert flap_ref - flap_tol <= flap_val <= flap_ref + flap_tol, (
                f"solver={solver}, dlc={dlc}: flapwise {flap_val:.3f} m outside "
                f"[{flap_ref - flap_tol:.3f}, {flap_ref + flap_tol:.3f}] m from Chapter 6 envelope"
            )
            assert edge_ref - edge_tol <= edge_val <= edge_ref + edge_tol, (
                f"solver={solver}, dlc={dlc}: edgewise {edge_val:.3f} m outside "
                f"[{edge_ref - edge_tol:.3f}, {edge_ref + edge_tol:.3f}] m from Chapter 6 envelope"
            )


@pytest.mark.skipif(not _HAS_RESULTS, reason=_MISSING_RESULTS)
def test_worst_out_of_plane_tip_deflection_matches_ch6_scale():
    _, results = _load_solver_results_by_solver_and_dlc()
    max_flapwise = max(v["flapwise"] for v in results.values())
    lo = 0.75 * CH6_MAX_OUT_OF_PLANE_TIP_DEFLECTION_M
    hi = 1.25 * CH6_MAX_OUT_OF_PLANE_TIP_DEFLECTION_M
    assert lo <= max_flapwise <= hi, (
        f"Worst flapwise tip deflection {max_flapwise:.3f} m outside [{lo:.3f}, {hi:.3f}] m "
        f"around Chapter 6 reference {CH6_MAX_OUT_OF_PLANE_TIP_DEFLECTION_M:.1f} m"
    )
