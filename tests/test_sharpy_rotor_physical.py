from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest


def _build_rotor_case(case_name: str, route: Path, *, dt: float, wsp: float):
    import sharpy.cases.templates.template_wt as template_wt

    excel_file = _nrel5mw_excel_from_sharpy()
    if not excel_file.exists():
        pytest.skip(f"SHARPy NREL workbook not found at {excel_file}")

    chord_panels = np.array([8], dtype=int)
    op_params = {
        "rotation_velocity": 1.366190,
        "pitch_deg": 0.0,
        "wsp": wsp,
        "dt": dt,
    }
    geom_params = {
        "chord_panels": chord_panels,
        "tol_remove_points": 1e-8,
        "n_points_camber": 100,
        "m_distribution": "uniform",
    }
    excel_description = {
        "excel_file_name": str(excel_file),
        "excel_sheet_parameters": "parameters",
        "excel_sheet_structural_blade": "structural_blade",
        "excel_sheet_discretization_blade": "discretization_blade",
        "excel_sheet_aero_blade": "aero_blade",
        "excel_sheet_airfoil_info": "airfoil_info",
        "excel_sheet_airfoil_chord": "airfoil_coord",
    }
    options = {
        "camber_effect_on_twist": False,
        "user_defined_m_distribution_type": None,
        "include_polars": False,
        "separate_blades": False,
    }

    rotor, _hub_nodes = template_wt.rotor_from_excel_type03(
        op_params,
        geom_params,
        excel_description,
        options,
    )
    return rotor


def _nrel5mw_excel_from_sharpy() -> Path:
    import sharpy.utils.sharpydir as sharpydir

    return (
        Path(sharpydir.SharpyDir)
        / "docs"
        / "source"
        / "content"
        / "example_notebooks"
        / "source"
        / "type04_db_nrel5mw_oc3_v06.xlsx"
    )


@pytest.mark.slow
@pytest.mark.integration
def test_sharpy_rotor_modal_reference_frequencies(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    gc = pytest.importorskip("sharpy.utils.generate_cases", reason="SHARPy not installed")
    sharpy_main = pytest.importorskip(
        "sharpy.sharpy_main", reason="SHARPy CLI entrypoint not available"
    )
    pytest.importorskip(
        "sharpy.cases.templates.template_wt", reason="SHARPy turbine template not available"
    )

    from sharpy.utils.constants import deg2rad

    monkeypatch.chdir(tmp_path)

    route = tmp_path
    case = "rotor_modal_reference"
    rotation_velocity = 1.366190
    dt = (4.0 * deg2rad) / rotation_velocity
    wsp = 11.4
    air_density = 1.225

    rotor = _build_rotor_case(case, route, dt=dt, wsp=wsp)

    sim = gc.SimulationInformation()
    sim.set_default_values()
    sim.solvers["SHARPy"]["flow"] = ["BeamLoader", "Modal"]
    sim.solvers["SHARPy"]["case"] = case
    sim.solvers["SHARPy"]["route"] = str(route) + "/"
    sim.solvers["SHARPy"]["write_log"] = True
    sim.solvers["SHARPy"]["write_screen"] = "off"
    sim.set_variable_all_dicts("dt", dt)
    sim.set_variable_all_dicts("rho", air_density)
    sim.solvers["BeamLoader"]["unsteady"] = "on"
    sim.solvers["Modal"]["write_modes_vtk"] = False
    sim.solvers["Modal"]["save_data"] = True

    gc.clean_test_files(sim.solvers["SHARPy"]["route"], sim.solvers["SHARPy"]["case"])
    rotor.generate_h5_files(sim.solvers["SHARPy"]["route"], sim.solvers["SHARPy"]["case"])
    sim.generate_solver_file()

    solver_path = route / f"{case}.sharpy"
    sharpy_main.main(["", str(solver_path)])

    output_path = tmp_path / "output" / case / "beam_modal_analysis" / "frequencies.dat"
    assert output_path.exists(), "SHARPy modal output frequencies.dat was not generated"

    freq_data = np.atleast_2d(np.genfromtxt(output_path))

    flap_1 = np.average([0.6664, 0.6296, 0.6675, 0.6686, 0.6993, 0.7019]) * 2.0 * np.pi
    edge_1 = np.average([1.0793, 1.0740, 1.0898, 1.0877]) * 2.0 * np.pi
    flap_2 = np.average([1.9337, 1.6507, 1.9223, 1.8558, 2.0205, 1.9601]) * 2.0 * np.pi

    assert np.isfinite(freq_data[0, 0])
    assert np.isfinite(freq_data[0, 3])
    assert np.isfinite(freq_data[0, 6])
    assert abs(freq_data[0, 0] - flap_1) < 1.0
    assert abs(freq_data[0, 3] - edge_1) < 1.0
    assert abs(freq_data[0, 6] - flap_2) < 1.0


@pytest.mark.slow
@pytest.mark.integration
def test_sharpy_rotor_aero_produces_positive_cp_ct(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    gc = pytest.importorskip("sharpy.utils.generate_cases", reason="SHARPy not installed")
    sharpy_main = pytest.importorskip(
        "sharpy.sharpy_main", reason="SHARPy CLI entrypoint not available"
    )
    pytest.importorskip(
        "sharpy.cases.templates.template_wt", reason="SHARPy turbine template not available"
    )

    from sharpy.utils.constants import deg2rad

    monkeypatch.chdir(tmp_path)

    route = tmp_path
    case = "rotor_aero_reference"
    rotation_velocity = 1.366190
    wsp = 11.4
    air_density = 1.225
    dphi = 4.0 * deg2rad
    dt = dphi / rotation_velocity
    mstar = int(1 * 2.0 * np.pi / dphi)
    time_steps = 4

    rotor = _build_rotor_case(case, route, dt=dt, wsp=wsp)

    sim = gc.SimulationInformation()
    sim.set_default_values()
    sim.solvers["SHARPy"]["flow"] = [
        "BeamLoader",
        "AerogridLoader",
        "StaticCoupled",
        "DynamicCoupled",
        "AeroForcesCalculator",
    ]
    sim.solvers["SHARPy"]["case"] = case
    sim.solvers["SHARPy"]["route"] = str(route) + "/"
    sim.solvers["SHARPy"]["write_log"] = True
    sim.solvers["SHARPy"]["write_screen"] = "off"
    sim.set_variable_all_dicts("dt", dt)
    sim.set_variable_all_dicts("rho", air_density)

    sim.solvers["SteadyVelocityField"]["u_inf"] = wsp
    sim.solvers["SteadyVelocityField"]["u_inf_direction"] = np.array([0.0, 0.0, 1.0])
    sim.set_variable_all_dicts("velocity_field_input", sim.solvers["SteadyVelocityField"])

    sim.solvers["BeamLoader"]["unsteady"] = "on"

    sim.solvers["AerogridLoader"]["unsteady"] = "on"
    sim.solvers["AerogridLoader"]["mstar"] = mstar
    sim.solvers["AerogridLoader"]["freestream_dir"] = np.array([0.0, 0.0, 0.0])
    sim.solvers["AerogridLoader"]["wake_shape_generator"] = "HelicoidalWake"
    sim.solvers["AerogridLoader"]["wake_shape_generator_input"] = {
        "u_inf": wsp,
        "u_inf_direction": sim.solvers["SteadyVelocityField"]["u_inf_direction"],
        "rotation_velocity": rotation_velocity * np.array([0.0, 0.0, 1.0]),
        "dt": dt,
        "dphi1": dphi,
        "ndphi1": mstar,
        "r": 1.0,
        "dphimax": 10.0 * deg2rad,
    }

    sim.solvers["NonLinearStatic"]["max_iterations"] = 200
    sim.solvers["NonLinearStatic"]["num_load_steps"] = 1
    sim.solvers["NonLinearStatic"]["min_delta"] = 1e-5

    sim.solvers["StaticUvlm"]["horseshoe"] = False
    sim.solvers["StaticUvlm"]["num_cores"] = 8
    sim.solvers["StaticUvlm"]["n_rollup"] = 0
    sim.solvers["StaticUvlm"]["rollup_dt"] = dt
    sim.solvers["StaticUvlm"]["rollup_aic_refresh"] = 1
    sim.solvers["StaticUvlm"]["rollup_tolerance"] = 1e-8
    sim.solvers["StaticUvlm"]["velocity_field_generator"] = "SteadyVelocityField"
    sim.solvers["StaticUvlm"]["velocity_field_input"] = sim.solvers["SteadyVelocityField"]

    sim.solvers["StaticCoupled"]["structural_solver"] = "NonLinearStatic"
    sim.solvers["StaticCoupled"]["structural_solver_settings"] = sim.solvers["NonLinearStatic"]
    sim.solvers["StaticCoupled"]["aero_solver"] = "StaticUvlm"
    sim.solvers["StaticCoupled"]["aero_solver_settings"] = sim.solvers["StaticUvlm"]
    sim.solvers["StaticCoupled"]["tolerance"] = 1e-6
    sim.solvers["StaticCoupled"]["n_load_steps"] = 0
    sim.solvers["StaticCoupled"]["relaxation_factor"] = 0.0

    sim.solvers["StepUvlm"]["convection_scheme"] = 2
    sim.solvers["StepUvlm"]["num_cores"] = 8

    sim.solvers["DynamicCoupled"]["structural_solver"] = "RigidDynamicPrescribedStep"
    sim.solvers["DynamicCoupled"]["structural_solver_settings"] = sim.solvers[
        "RigidDynamicPrescribedStep"
    ]
    sim.solvers["DynamicCoupled"]["aero_solver"] = "StepUvlm"
    sim.solvers["DynamicCoupled"]["aero_solver_settings"] = sim.solvers["StepUvlm"]
    sim.solvers["DynamicCoupled"]["postprocessors"] = ["Cleanup"]
    sim.solvers["DynamicCoupled"]["postprocessors_settings"] = {"Cleanup": sim.solvers["Cleanup"]}
    sim.solvers["DynamicCoupled"]["minimum_steps"] = 0

    sim.solvers["AeroForcesCalculator"]["write_text_file"] = True
    sim.solvers["AeroForcesCalculator"]["text_file_name"] = "aeroforces.txt"
    sim.solvers["AeroForcesCalculator"]["screen_output"] = False
    sim.solvers["AeroForcesCalculator"]["coefficients"] = False
    sim.solvers["AeroForcesCalculator"]["u_inf_dir"] = [0.0, 0.0, 1.0]

    sim.define_num_steps(time_steps)
    sim.with_forced_vel = True
    sim.for_vel = np.zeros((time_steps, 6), dtype=float)
    sim.for_vel[:, 5] = rotation_velocity
    sim.for_acc = np.zeros((time_steps, 6), dtype=float)
    sim.with_dynamic_forces = True
    sim.dynamic_forces = np.zeros(
        (time_steps, rotor.StructuralInformation.num_node, 6), dtype=float
    )

    gc.clean_test_files(sim.solvers["SHARPy"]["route"], sim.solvers["SHARPy"]["case"])
    rotor.generate_h5_files(sim.solvers["SHARPy"]["route"], sim.solvers["SHARPy"]["case"])
    sim.generate_solver_file()
    sim.generate_dyn_file(time_steps)

    solver_path = route / f"{case}.sharpy"
    sharpy_main.main(["", str(solver_path)])

    force_file = tmp_path / "output" / case / "forces" / "forces_aeroforces.txt"
    moment_file = tmp_path / "output" / case / "forces" / "moments_aeroforces.txt"
    assert force_file.exists(), "SHARPy did not write aerodynamic force history"
    assert moment_file.exists(), "SHARPy did not write aerodynamic moment history"

    force_data = np.genfromtxt(force_file, delimiter=",", comments="#")
    moment_data = np.genfromtxt(moment_file, delimiter=",", comments="#")
    force_data = np.atleast_2d(force_data)
    moment_data = np.atleast_2d(moment_data)

    fz = force_data[-1, 3] + force_data[-1, 6]
    mz = moment_data[-1, 3] + moment_data[-1, 6]

    rotor_radius = float(
        np.max(np.linalg.norm(rotor.StructuralInformation.coordinates[:, :2], axis=1))
    )
    rotor_area = np.pi * rotor_radius**2

    cp = abs(mz) * rotation_velocity / (0.5 * air_density * rotor_area * wsp**3)
    ct = abs(fz) / (0.5 * air_density * rotor_area * wsp**2)

    assert np.isfinite(cp)
    assert np.isfinite(ct)
    assert cp > 1e-7
    assert ct > 1e-8
    assert cp < 0.2
    assert ct < 0.2
