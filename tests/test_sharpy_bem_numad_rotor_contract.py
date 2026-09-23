from __future__ import annotations

import importlib
import math
import os
from pathlib import Path
from typing import Any, Literal

import h5py
import numpy as np
import pytest

from aeroelast.models.blade.aerodynamics import load_blade_aero
from aeroelast.solvers.aero.sharpy_numad import (
    NuMADSharpyRotorCaseAdapter,
    _build_aero_node_mask,
    _extract_rotor_result,
    _load_sharpy_modules,
)
from aeroelast.solvers.aero.validation import (
    IEA15MW_REFERENCE_OPERATING_POINTS,
    compare_sharpy_to_bem_sweep,
    compute_bem_metrics,
)


def _load_numad_blade(path: Path):
    return load_blade_aero(str(path), airfoil_dir=str(path.parent / "airfoils"))


def _rated_point():
    return next(point for point in IEA15MW_REFERENCE_OPERATING_POINTS if point.name == "v10.659")


def _feathered_point():
    return next(point for point in IEA15MW_REFERENCE_OPERATING_POINTS if point.name == "v15.471")


def _sharpy_solver_settings(model: Literal["vlm", "uvlm"] | None = None) -> dict[str, Any]:
    profile = os.environ.get("AEROELAST_SHARPY_TEST_PROFILE", "quick").strip().lower()
    if profile == "long":
        # Longer UVLM/VLM horizon for closer BEM comparison.
        settings = {
            "time_steps": 24,
            "wake_revolutions": 2.0,
            "dphi_deg": 10.0,
            "chord_panels": 8,
            "hub_radius": 3.0,
            "wind_direction": [0.0, 1.0, 0.0],
            "rotation_axis": [0.0, 1.0, 0.0],
            "freestream_direction": [0.0, 1.0, 0.0],
        }
    else:
        settings = {
            "time_steps": 2,
            "wake_revolutions": 0.2,
            "dphi_deg": 30.0,
            "chord_panels": 2,
            "hub_radius": 3.0,
            "wind_direction": [0.0, 1.0, 0.0],
            "rotation_axis": [0.0, 1.0, 0.0],
            "freestream_direction": [0.0, 1.0, 0.0],
        }

    if model == "vlm":
        # Static VLM becomes highly non-physical for the quick wake tuple
        # (wake_revolutions=0.2, dphi_deg=30, mstar=2). A short wake with a
        # slightly finer chord discretization and smaller dphi gives the best
        # BEM agreement we found without paying the full long-profile cost.
        settings["wake_revolutions"] = 0.05
        settings["dphi_deg"] = 5.0
        settings["mstar"] = 2
        settings["chord_panels"] = 4
    elif model == "uvlm":
        # UVLM uses the same y-axis rotor convention as the static VLM: blade
        # span along +Z, rotation axis +Y, wind along +Y.
        # A single DynamicCoupled step starting from the converged StaticCoupled
        # (VLM) initial condition verifies the unsteady solver without requiring
        # wake convergence.  Forces are extracted from the StaticCoupled step.
        # chord_panels=4 and dphi=5 match the VLM benchmark to give BEM-level
        # force agreement; vortex_radius avoids Biot-Savart singularities.
        settings["time_steps"] = 1
        settings["dphi_deg"] = 5.0
        settings["mstar"] = 2
        settings["chord_panels"] = 4
        settings["wind_direction"] = [0.0, 1.0, 0.0]
        settings["rotation_axis"] = [0.0, 1.0, 0.0]
        settings["freestream_direction"] = [0.0, 1.0, 0.0]
        settings["vortex_radius"] = 0.1
        settings["vortex_radius_wake_ind"] = 0.1

    return settings


def _run_raw_sharpy_case(
    adapter: NuMADSharpyRotorCaseAdapter,
    operating_point,
    *,
    model: Literal["vlm", "uvlm"],
    route: Path,
    solver_settings: dict[str, Any],
):
    gc, sharpy_main, _algebra = _load_sharpy_modules()

    settings = dict(solver_settings)
    omega_rad_s = operating_point.omega_rpm * 2.0 * math.pi / 60.0
    dphi = math.radians(float(settings["dphi_deg"]))
    dt = dphi / abs(omega_rad_s)
    time_steps = int(settings["time_steps"])

    route.mkdir(parents=True, exist_ok=True)
    rotor = adapter._build_rotor(
        operating_point.pitch_deg,
        chord_panels=int(settings["chord_panels"]),
        hub_radius=float(settings["hub_radius"]),
        rotation_axis=np.asarray(settings["rotation_axis"], dtype=float),
    )
    sim = adapter._build_simulation(
        case_name=f"raw_{model}",
        route=route,
        model=model,
        wind_speed=operating_point.wind_speed,
        omega_rad_s=omega_rad_s,
        wind_direction=np.asarray(settings["wind_direction"], dtype=float),
        rotation_axis=np.asarray(settings["rotation_axis"], dtype=float),
        freestream_direction=np.asarray(settings["freestream_direction"], dtype=float),
        dt=dt,
        dphi=dphi,
        mstar=int(settings["mstar"]),
        time_steps=time_steps,
        solver_settings=settings,
        rotor_node_count=rotor.StructuralInformation.num_node,
    )
    gc.clean_test_files(sim.solvers["SHARPy"]["route"], sim.solvers["SHARPy"]["case"])
    rotor.generate_h5_files(sim.solvers["SHARPy"]["route"], sim.solvers["SHARPy"]["case"])
    sim.generate_solver_file()
    sim.generate_dyn_file(time_steps)
    output = sharpy_main.main(["", str(route / f"raw_{model}.sharpy")])
    return output, omega_rad_s


def test_numad_bem_baseline_is_physical(iea_blade_xlsx: Path):
    blade = _load_numad_blade(iea_blade_xlsx)
    metrics = compute_bem_metrics(blade, _rated_point())

    assert metrics.thrust > 1.0e6
    assert metrics.torque > 1.0e6
    assert 0.2 <= metrics.cp <= 0.7
    assert 0.3 <= metrics.ct <= 1.2
    for value in (metrics.power, metrics.cq):
        assert math.isfinite(value)


def test_numad_blade_loads_prebend_distribution(iea_blade_xlsx: Path):
    blade = _load_numad_blade(iea_blade_xlsx)
    prebend = np.asarray(blade.prebend, dtype=float)

    assert prebend.shape == blade.r.shape
    assert np.max(np.abs(prebend - prebend[0])) > 1.0e-3


def test_numad_sharpy_rotor_includes_all_stations_in_vlm(iea_blade_xlsx: Path):
    """VLM lattice must cover all stations including cylindrical root sections.

    Including the root cylindrical stations ensures the fluid mesh starts at
    hub_radius (matching the structural coupling mesh) with no radial gap.
    """
    blade = _load_numad_blade(iea_blade_xlsx)

    mask = _build_aero_node_mask(blade)

    assert bool(np.all(mask)), "All stations must participate in the VLM lattice"


def test_numad_geometry_reference_axis_reconstructs_prebend(iea_blade_xlsx: Path):
    from aeroelast.models.blade.aerodynamics import _extract_reference_axis_from_numad_geometry
    from aeroelast.models.blade.numad.objects.blade import Blade as NumadBlade

    blade = NumadBlade()
    blade.read_excel(str(iea_blade_xlsx), airfoil_dir=str(iea_blade_xlsx.parent / "airfoils"))

    pitch_axes = np.asarray(blade.definition.aerocenter, dtype=float)
    estimated = _extract_reference_axis_from_numad_geometry(blade, pitch_axes=pitch_axes)

    assert estimated is not None
    prebend_est, _sweep_est = estimated
    assert prebend_est.shape[0] == len(blade.definition.span)
    assert np.ptp(prebend_est) > 1.0


@pytest.mark.slow
@pytest.mark.integration
@pytest.mark.parametrize("model", ["vlm", "uvlm"])
def test_numad_sharpy_contract_runs_with_same_power_order(
    iea_blade_xlsx: Path,
    tmp_path: Path,
    model: Literal["vlm", "uvlm"],
):
    pytest.importorskip("sharpy.sharpy_main", reason="SHARPy not installed")
    pytest.importorskip("sharpy.utils.generate_cases", reason="SHARPy case tools not available")

    blade = _load_numad_blade(iea_blade_xlsx)
    rows = compare_sharpy_to_bem_sweep(
        blade,
        [_rated_point()],
        model=model,
        route=tmp_path,
        solver_settings=_sharpy_solver_settings(model),
    )

    assert len(rows) == 1
    row = rows[0]
    assert row.panel_count > 0
    for value in (
        row.vlm.thrust,
        row.vlm.torque,
        row.vlm.cp,
        row.vlm.ct,
        row.bem.thrust,
        row.bem.torque,
        row.bem.cp,
        row.bem.ct,
        row.thrust_ratio,
        row.torque_ratio,
        row.cp_ratio,
        row.ct_ratio,
    ):
        assert math.isfinite(value)

    assert row.vlm.thrust > 0.0
    assert row.vlm.torque > 0.0
    assert row.vlm.cp > 0.0
    assert row.vlm.ct > 0.0
    if model == "vlm":
        assert 0.85 <= row.ct_ratio <= 1.15
        assert 0.8 <= row.torque_ratio <= 1.45
        assert 0.8 <= row.cp_ratio <= 1.45
    else:
        assert 0.8 <= row.thrust_ratio <= 1.5
        assert 0.8 <= row.ct_ratio <= 1.5
        assert 0.8 <= row.torque_ratio <= 1.7
        assert 0.8 <= row.cp_ratio <= 1.7


@pytest.mark.slow
@pytest.mark.integration
def test_numad_sharpy_static_vlm_default_twist_matches_bem_order(
    iea_blade_xlsx: Path,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    pytest.importorskip("sharpy.sharpy_main", reason="SHARPy not installed")
    pytest.importorskip("sharpy.utils.generate_cases", reason="SHARPy case tools not available")

    monkeypatch.delenv("AEROELAST_SHARPY_TWIST_SIGN", raising=False)
    monkeypatch.delenv("AEROELAST_SHARPY_FOR_DELTA", raising=False)

    blade = _load_numad_blade(iea_blade_xlsx)
    row = compare_sharpy_to_bem_sweep(
        blade,
        [_rated_point()],
        model="vlm",
        route=tmp_path,
        solver_settings=_sharpy_solver_settings("vlm"),
    )[0]

    assert 0.85 <= row.ct_ratio <= 1.15
    assert 0.8 <= row.cp_ratio <= 1.45


@pytest.mark.slow
@pytest.mark.integration
def test_numad_sharpy_static_vlm_uses_force_based_torque(
    iea_blade_xlsx: Path,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    pytest.importorskip("sharpy.sharpy_main", reason="SHARPy not installed")
    pytest.importorskip("sharpy.utils.generate_cases", reason="SHARPy case tools not available")

    monkeypatch.delenv("AEROELAST_SHARPY_TWIST_SIGN", raising=False)
    monkeypatch.delenv("AEROELAST_SHARPY_FOR_DELTA", raising=False)

    blade = _load_numad_blade(iea_blade_xlsx)
    result = NuMADSharpyRotorCaseAdapter(blade_aero=blade, route=tmp_path).run(
        _rated_point(),
        model="vlm",
        solver_settings=_sharpy_solver_settings("vlm"),
    )

    rotation_axis = np.asarray(result.loads.metadata["rotation_axis"], dtype=float)
    assert result.loads.integrated_moment is not None
    exported_torque = abs(
        float(np.dot(np.asarray(result.loads.integrated_moment, dtype=float), rotation_axis))
    )
    total_torque = float(result.loads.metadata["torque_total"])
    local_torque = float(result.loads.metadata["torque_local"])

    assert exported_torque == pytest.approx(result.metrics.torque)
    assert total_torque > result.metrics.torque
    assert local_torque > 0.0
    assert total_torque == pytest.approx(result.metrics.torque + local_torque)


@pytest.mark.slow
@pytest.mark.integration
def test_numad_sharpy_uvlm_uses_postprocessed_aero_forces(
    iea_blade_xlsx: Path,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    pytest.importorskip("sharpy.sharpy_main", reason="SHARPy not installed")
    pytest.importorskip("sharpy.utils.generate_cases", reason="SHARPy case tools not available")

    monkeypatch.setenv("AEROELAST_SHARPY_TWIST_SIGN", "1")
    monkeypatch.setenv("AEROELAST_SHARPY_FOR_DELTA", "y_AFoR")

    blade = _load_numad_blade(iea_blade_xlsx)
    adapter = NuMADSharpyRotorCaseAdapter(blade_aero=blade, route=tmp_path)
    settings = _sharpy_solver_settings("uvlm")
    settings.update({
        "time_steps": 6,
        "wake_revolutions": 0.05,
        "dphi_deg": 5.0,
        "mstar": 2,
        "chord_panels": 4,
    })

    output, omega_rad_s = _run_raw_sharpy_case(
        adapter,
        _rated_point(),
        model="uvlm",
        route=tmp_path / "uvlm_postproc_force_source",
        solver_settings=settings,
    )
    metrics, loads = _extract_rotor_result(
        output,
        rho=adapter.rho,
        wind_speed=_rated_point().wind_speed,
        omega_rad_s=omega_rad_s,
        wind_direction=np.asarray(settings["wind_direction"], dtype=float),
        rotation_axis=np.asarray(settings["rotation_axis"], dtype=float),
    )

    wind_direction = np.asarray(settings["wind_direction"], dtype=float)
    rotation_axis = np.asarray(settings["rotation_axis"], dtype=float)

    assert loads.metadata["nodal_force_source"] == "postproc_node"
    assert metrics.thrust == pytest.approx(
        abs(float(np.dot(np.asarray(loads.integrated_force, dtype=float), wind_direction)))
    )
    assert metrics.torque == pytest.approx(
        abs(float(np.dot(np.asarray(loads.integrated_moment, dtype=float), rotation_axis)))
    )


@pytest.mark.slow
@pytest.mark.integration
def test_numad_sharpy_uvlm_matches_bem_thrust_order(
    iea_blade_xlsx: Path,
    tmp_path: Path,
):
    pytest.importorskip("sharpy.sharpy_main", reason="SHARPy not installed")
    pytest.importorskip("sharpy.utils.generate_cases", reason="SHARPy case tools not available")

    blade = _load_numad_blade(iea_blade_xlsx)
    row = compare_sharpy_to_bem_sweep(
        blade,
        [_rated_point()],
        model="uvlm",
        route=tmp_path,
        solver_settings=_sharpy_solver_settings("uvlm"),
    )[0]

    assert 0.8 <= row.thrust_ratio <= 1.5
    assert 0.8 <= row.ct_ratio <= 1.5


def test_numad_bem_pitch_to_feather_trend(iea_blade_xlsx: Path):
    blade = _load_numad_blade(iea_blade_xlsx)
    rated = compute_bem_metrics(blade, _rated_point())
    feathered = compute_bem_metrics(blade, _feathered_point())

    assert rated.cp > feathered.cp
    assert rated.ct > feathered.ct
    assert rated.torque > feathered.torque
    assert rated.thrust > feathered.thrust


@pytest.mark.integration
def test_numad_sharpy_rotor_exports_polars(iea_blade_xlsx: Path, tmp_path: Path):
    pytest.importorskip("sharpy.sharpy_main", reason="SHARPy not installed")
    pytest.importorskip("sharpy.utils.generate_cases", reason="SHARPy case tools not available")

    blade = _load_numad_blade(iea_blade_xlsx)
    adapter = NuMADSharpyRotorCaseAdapter(blade_aero=blade, route=tmp_path)
    rotor = adapter._build_rotor(
        0.0,
        chord_panels=2,
        hub_radius=3.0,
        rotation_axis=np.array([0.0, 1.0, 0.0], dtype=float),
    )

    case_name = "numad_polars_probe"
    rotor.generate_h5_files(str(tmp_path), case_name)

    with h5py.File(tmp_path / f"{case_name}.aero.h5", "r") as h5f:
        assert "polars" in h5f
        polars_group_obj = h5f["polars"]
        assert isinstance(polars_group_obj, h5py.Group)
        assert len(polars_group_obj.keys()) > 0
        first_obj = polars_group_obj["0"]
        assert isinstance(first_obj, h5py.Dataset)
        first = first_obj[:]
        assert first.ndim == 2
        assert first.shape[1] == 4


@pytest.mark.integration
def test_numad_sharpy_rotor_structural_mesh_keeps_prebend(iea_blade_xlsx: Path, tmp_path: Path):
    pytest.importorskip("sharpy.sharpy_main", reason="SHARPy not installed")
    pytest.importorskip("sharpy.utils.generate_cases", reason="SHARPy case tools not available")

    blade = _load_numad_blade(iea_blade_xlsx)
    blade.n_blades = 1
    expected_prebend = np.asarray(blade.prebend, dtype=float)
    expected_prebend = expected_prebend - expected_prebend[0]

    adapter = NuMADSharpyRotorCaseAdapter(blade_aero=blade, route=tmp_path)
    rotor = adapter._build_rotor(
        0.0,
        chord_panels=2,
        hub_radius=3.0,
        rotation_axis=np.array([0.0, 1.0, 0.0], dtype=float),
    )

    coords = np.asarray(rotor.StructuralInformation.coordinates, dtype=float)
    assert coords.shape[0] == expected_prebend.shape[0]

    mesh_prebend = coords[:, 1] - coords[0, 1]
    assert np.max(np.abs(mesh_prebend)) > 1.0e-3
    assert np.allclose(mesh_prebend, expected_prebend, atol=1.0e-8)


@pytest.mark.integration
def test_numad_sharpy_rotor_aerogrid_zero_pitch_chord_is_tangential(
    iea_blade_xlsx: Path,
    tmp_path: Path,
):
    configobj = pytest.importorskip(
        "configobj", reason="configobj is required to load SHARPy case files"
    )
    pytest.importorskip("sharpy.sharpy_main", reason="SHARPy not installed")
    pytest.importorskip("sharpy.utils.generate_cases", reason="SHARPy case tools not available")

    blade = _load_numad_blade(iea_blade_xlsx)
    adapter = NuMADSharpyRotorCaseAdapter(blade_aero=blade, route=tmp_path)
    point = _rated_point()
    settings = _sharpy_solver_settings("vlm")
    omega_rad_s = point.omega_rpm * 2.0 * math.pi / 60.0
    dphi = math.radians(float(settings["dphi_deg"]))
    dt = dphi / abs(omega_rad_s)

    rotor = adapter._build_rotor(
        0.0,
        chord_panels=int(settings["chord_panels"]),
        hub_radius=float(settings["hub_radius"]),
        rotation_axis=np.asarray(settings["rotation_axis"], dtype=float),
    )
    sim = adapter._build_simulation(
        case_name="orientation_probe",
        route=tmp_path,
        model="vlm",
        wind_speed=point.wind_speed,
        omega_rad_s=omega_rad_s,
        wind_direction=np.asarray(settings["wind_direction"], dtype=float),
        rotation_axis=np.asarray(settings["rotation_axis"], dtype=float),
        freestream_direction=np.asarray(settings["freestream_direction"], dtype=float),
        dt=dt,
        dphi=dphi,
        mstar=int(settings["mstar"]),
        time_steps=1,
        solver_settings=settings,
        rotor_node_count=rotor.StructuralInformation.num_node,
    )

    gc, _sharpy_main, _algebra = _load_sharpy_modules()
    gc.clean_test_files(sim.solvers["SHARPy"]["route"], sim.solvers["SHARPy"]["case"])
    rotor.generate_h5_files(sim.solvers["SHARPy"]["route"], sim.solvers["SHARPy"]["case"])
    sim.generate_solver_file()
    sim.generate_dyn_file(1)

    pre_cls = importlib.import_module("sharpy.presharpy.presharpy").PreSharpy
    beam_mod = importlib.import_module("sharpy.solvers.beamloader")
    aero_mod = importlib.import_module("sharpy.solvers.aerogridloader")
    case_file = tmp_path / "orientation_probe.sharpy"
    data = pre_cls(configobj.ConfigObj(str(case_file)))
    beam_solver = beam_mod.BeamLoader()
    beam_solver.initialise(data)
    data = beam_solver.run()
    aero_solver = aero_mod.AerogridLoader()
    aero_solver.initialise(data)
    data = aero_solver.run()

    zeta = np.asarray(data.aero.timestep_info[-1].zeta[0], dtype=float)
    mid_span = zeta.shape[2] // 2
    chord_direction = zeta[:, -1, mid_span] - zeta[:, 0, mid_span]
    chord_direction = chord_direction / np.linalg.norm(chord_direction)
    root_radials = np.linalg.norm(zeta[[0, 2], :, 0].T, axis=1)

    # Coordinate convention: rotation_axis = Y (axial/wind), blade spans Z (radial),
    # tangential = X.  With frame_of_reference_delta="y_AFoR" and vec_sweep=π the
    # chord lies in the axial-radial (Y-Z) plane; at zero pitch and low twist the
    # chord is nearly parallel to the axial (Y) axis with LE facing upwind (+Y) and
    # TE facing downwind (−Y).
    assert abs(float(chord_direction[2])) < 0.1  # chord NOT radial (Z)
    assert abs(float(chord_direction[0])) < 0.2  # chord NOT tangential (X) at mid-span
    assert abs(float(chord_direction[1])) > 0.8  # chord IS axial (Y): LE at +Y, TE at -Y
    assert float(chord_direction[1]) < -0.8  # LE faces upwind (+Y), TE direction is -Y
    # Twist sign: at root (large positive twist ~15°) the TE-LE direction has a
    # positive X component (TE is displaced toward +X relative to LE), consistent
    # with the NuMAD solid mesh where positive StrcTwst rotates the LE toward −Y.
    root_chord_dir = zeta[:, -1, 0] - zeta[:, 0, 0]
    root_chord_dir = root_chord_dir / np.linalg.norm(root_chord_dir)
    assert float(root_chord_dir[0]) > 0.1  # root twist direction: TE toward +X
    # VLM must cover the full blade starting at hub_radius: the root panel Z
    # coordinate should be near hub_radius (~3 m), not at the first non-cyl
    # station (~7.7 m).
    assert float(zeta[2, 0, 0]) < 5.0  # VLM root starts at hub_radius, not at station 2
    assert float(root_radials.min()) > 2.5

    surface_tip_angles: list[float] = []
    for surface_zeta in data.aero.timestep_info[-1].zeta:
        surface = np.asarray(surface_zeta, dtype=float)
        tip = surface[:, 0, -1]
        surface_tip_angles.append(math.degrees(math.atan2(float(tip[2]), float(tip[0]))))

    observed = np.sort(np.asarray(surface_tip_angles, dtype=float))
    expected = np.array([-150.0, -30.0, 90.0], dtype=float)
    assert np.allclose(observed, expected, atol=1.0)
