"""SHARPy rotor case generation from aeroelast blade data.

This module builds a rigid three-bladed SHARPy wind-turbine rotor from the same
``BladeAero`` data used by the BEM/CCBlade path.  It is intended for validation
benchmarks where the physical rotor definition must be shared while SHARPy keeps
its own VLM/UVLM solver settings.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Literal, Mapping

import numpy as np

from aeroelast.models.blade.aerodynamics import BladeAero

from .types import AeroFSIRuntimeContext, RotorAeroLoads, RotorAeroState
from .validation import RotorOperatingPoint, RotorPerformanceMetrics

SharpyRotorModel = Literal["vlm", "uvlm"]


@dataclass(frozen=True)
class SharpyRotorRunResult:
    """Integrated rotor result returned by a SHARPy benchmark run."""

    loads: RotorAeroLoads
    metrics: RotorPerformanceMetrics
    panel_count: int
    case_path: Path


@dataclass(frozen=True)
class NuMADSharpyRotorCaseAdapter:
    """Build and run rigid SHARPy rotor cases from ``BladeAero``.

    The adapter deliberately treats the aeroelast contract as the source of the
    physical model: blade geometry, number of blades, operating point and metric
    normalisation.  SHARPy-specific settings such as wake length, time-step count
    and convection scheme remain configurable per solver mode.
    """

    blade_aero: BladeAero
    route: Path
    case_prefix: str = "numad_iea15mw_sharpy"
    chord_panels: int = 8
    hub_radius: float = 3.0
    n_points_camber: int = 100
    wake_revolutions: float = 5.0
    dphi_deg: float = 4.0
    time_steps: int = 90
    simulation_revolutions: float = 3.0
    polar_re: float = 1.0e7
    rho: float = 1.225
    num_cores: int = 1

    def build_data(
        self,
        context: AeroFSIRuntimeContext,
        state: RotorAeroState,
        settings: Mapping[str, Any],
    ) -> Any:
        """Build a lightweight data envelope for compatibility with SharpyCaseAdapter."""

        return {
            "context": context,
            "state": state,
            "settings": dict(settings),
        }

    def extract_loads(
        self,
        data: Any,
        context: AeroFSIRuntimeContext,
        state: RotorAeroState,
    ) -> RotorAeroLoads:
        """Return loads from a previously executed adapter data envelope."""

        loads = data.get("loads") if isinstance(data, dict) else None
        if not isinstance(loads, RotorAeroLoads):
            raise RuntimeError("No SHARPy loads are available; run the rotor case first")
        return loads

    def run(
        self,
        operating_point: RotorOperatingPoint,
        *,
        model: SharpyRotorModel,
        solver_settings: Mapping[str, Any] | None = None,
    ) -> SharpyRotorRunResult:
        """Run one SHARPy rotor case and return contract-normalised metrics."""

        model = _normalise_model(model)
        settings = dict(solver_settings or {})
        omega_rad_s = operating_point.omega_rpm * 2.0 * math.pi / 60.0
        dphi_deg = float(settings.pop("dphi_deg", self.dphi_deg))
        if dphi_deg <= 0.0:
            raise ValueError("dphi_deg must be positive")
        dphi = math.radians(dphi_deg)

        raw_dt = settings.pop("dt", None)
        if raw_dt is None:
            # For rotating cases, couple dt to angular resolution dphi; for parked
            # cases, fall back to a conservative default time step.
            dt = dphi / abs(omega_rad_s) if abs(omega_rad_s) > 1.0e-9 else 0.1
        else:
            dt = float(raw_dt)
        if dt <= 0.0:
            raise ValueError("dt must be positive")

        raw_time_steps = settings.pop("time_steps", None)
        if raw_time_steps is None:
            simulation_revolutions = float(
                settings.pop("simulation_revolutions", self.simulation_revolutions)
            )
            if simulation_revolutions <= 0.0:
                raise ValueError("simulation_revolutions must be positive")
            steps_per_revolution = max(1, int(math.ceil(2.0 * math.pi / dphi)))
            auto_time_steps = max(1, int(math.ceil(simulation_revolutions * steps_per_revolution)))
            time_steps = max(int(self.time_steps), auto_time_steps)
        else:
            time_steps = int(raw_time_steps)
        if time_steps < 1:
            raise ValueError("time_steps must be at least 1")

        wake_revolutions = float(settings.pop("wake_revolutions", self.wake_revolutions))
        mstar = max(1, int(settings.pop("mstar", wake_revolutions * 2.0 * math.pi / dphi)))
        chord_panels = int(settings.pop("chord_panels", self.chord_panels))
        hub_radius = float(settings.pop("hub_radius", self.hub_radius))
        n_points_camber = int(settings.pop("n_points_camber", self.n_points_camber))
        polar_re = float(settings.pop("polar_re", self.polar_re))
        default_axis = [0.0, 1.0, 0.0]
        # Both VLM and UVLM use the y-axis rotor convention: blade span along
        # +Z, rotation axis +Y, wind along +Y.  Callers can still override.
        wind_direction = _as_unit_vector(
            settings.pop("wind_direction", default_axis),
            name="wind_direction",
        )
        rotation_axis = _as_unit_vector(
            settings.pop("rotation_axis", default_axis),
            name="rotation_axis",
        )
        raw_freestream_direction = settings.pop("freestream_direction", None)
        freestream_direction = (
            wind_direction.copy()
            if raw_freestream_direction is None
            else _as_unit_vector(raw_freestream_direction, name="freestream_direction")
        )

        case_name = _case_name(self.case_prefix, model, operating_point.name)
        route = Path(settings.pop("route", self.route)).expanduser().resolve()
        route.mkdir(parents=True, exist_ok=True)

        rotor = self._build_rotor(
            operating_point.pitch_deg,
            chord_panels=chord_panels,
            hub_radius=hub_radius,
            n_points_camber=n_points_camber,
            rotation_axis=rotation_axis,
            polar_re=polar_re,
        )
        sim = self._build_simulation(
            case_name=case_name,
            route=route,
            model=model,
            wind_speed=operating_point.wind_speed,
            omega_rad_s=omega_rad_s,
            wind_direction=wind_direction,
            rotation_axis=rotation_axis,
            freestream_direction=freestream_direction,
            dt=dt,
            dphi=dphi,
            mstar=mstar,
            time_steps=time_steps,
            solver_settings=settings,
            rotor_node_count=rotor.StructuralInformation.num_node,
        )

        gc, sharpy_main, _algebra = _load_sharpy_modules()
        gc.clean_test_files(sim.solvers["SHARPy"]["route"], sim.solvers["SHARPy"]["case"])
        rotor.generate_h5_files(sim.solvers["SHARPy"]["route"], sim.solvers["SHARPy"]["case"])
        sim.generate_solver_file()
        sim.generate_dyn_file(time_steps)

        output = sharpy_main.main(["", str(route / f"{case_name}.sharpy")])
        # VLM: StaticCoupled modifies timestep_info[0] in-place; AeroForcesCalculator
        # writes postproc_node["aero_steady_forces"] there.  Extract with index=-1 (==0).
        # UVLM: StaticCoupled also modifies timestep_info[0] in-place, but DynamicCoupled
        # appends timestep_info[1].  Use index=0 to get the converged quasi-steady result
        # rather than the first unsteady dynamic step (index=1) which is still transient.
        ts_index = 0 if model == "uvlm" else -1
        metrics, loads = _extract_rotor_result(
            output,
            rho=self.rho,
            wind_speed=operating_point.wind_speed,
            omega_rad_s=omega_rad_s,
            wind_direction=wind_direction,
            rotation_axis=rotation_axis,
            timestep_index=ts_index,
        )

        panel_count = _panel_count(output)
        return SharpyRotorRunResult(
            loads=loads,
            metrics=metrics,
            panel_count=panel_count,
            case_path=route / f"{case_name}.sharpy",
        )

    def _build_rotor(
        self,
        pitch_deg: float,
        *,
        chord_panels: int | None = None,
        hub_radius: float | None = None,
        n_points_camber: int | None = None,
        rotation_axis: np.ndarray | None = None,
        polar_re: float | None = None,
    ) -> Any:
        gc, _sharpy_main, _algebra = _load_sharpy_modules()
        chord_panels = self.chord_panels if chord_panels is None else chord_panels
        hub_radius = self.hub_radius if hub_radius is None else hub_radius
        n_points_camber = self.n_points_camber if n_points_camber is None else n_points_camber
        polar_re = self.polar_re if polar_re is None else float(polar_re)
        rotation_axis = (
            np.array([0.0, 1.0, 0.0], dtype=float) if rotation_axis is None else rotation_axis
        )
        if hub_radius < 0.0:
            raise ValueError("hub_radius must be non-negative")
        if n_points_camber < 4:
            raise ValueError("n_points_camber must be at least 4")

        blade = gc.AeroelasticInformation()
        radii = np.asarray(self.blade_aero.r, dtype=float)
        if radii.size == 0:
            raise ValueError("BladeAero must define at least one radial station")
        # Keep the blade span distribution but place the first station at hub_radius.
        radii = radii - radii[0] + float(hub_radius)
        prebend = np.asarray(self.blade_aero.prebend, dtype=float)
        sweep = np.asarray(self.blade_aero.sweep, dtype=float)
        if prebend.shape != radii.shape:
            raise ValueError("BladeAero prebend distribution must match radial station count")
        if sweep.shape != radii.shape:
            raise ValueError("BladeAero sweep distribution must match radial station count")
        prebend = prebend - prebend[0]
        sweep = sweep - sweep[0]
        chords = np.asarray(self.blade_aero.chord, dtype=float)
        # Negate twist: SHARPy's rotation_matrix_around_axis(+Z, θ) rotates the
        # chord CCW when viewed from the tip, which would place the LE toward +Y
        # for positive θ.  For an upwind turbine the LE must face the relative
        # wind [+X, −Y], so the rotation must be CW (negative angle).  The NuMAD
        # solid mesh already has LE at [+X, −Y] for positive StrcTwst values, so
        # we negate here to align the VLM lattice with the solid mesh convention.
        twists = np.asarray(self.blade_aero.twist, dtype=float)
        pitch_axes = np.array(
            [station.pitch_axis for station in self.blade_aero.stations], dtype=float
        )
        if (len(radii) - 1) % 2 != 0:
            raise ValueError("SHARPy 3-node beam elements require an odd number of blade stations")

        node_coordinates = np.zeros((len(radii), 3), dtype=float)
        node_coordinates[:, 2] = radii
        node_coordinates[:, 1] = prebend
        node_coordinates[:, 0] = sweep
        blade.StructuralInformation.generate_1to1_from_vectors(
            num_node_elem=3,
            num_node=len(radii),
            num_elem=(len(radii) - 1) // 2,
            coordinates=node_coordinates,
            stiffness_db=_beam_stiffness_db((len(radii) - 1) // 2),
            mass_db=_beam_mass_db((len(radii) - 1) // 2),
            frame_of_reference_delta="y_AFoR",
            vec_node_structural_twist=twists,
            num_lumped_mass=0,
        )
        blade.StructuralInformation.boundary_conditions = np.zeros(len(radii), dtype=int)
        blade.StructuralInformation.boundary_conditions[0] = 1
        blade.StructuralInformation.boundary_conditions[-1] = -1

        aero_node_mask = _build_aero_node_mask(self.blade_aero)

        blade.AerodynamicInformation.create_aerodynamics_from_vec(
            blade.StructuralInformation,
            aero_node_mask,
            chords,
            np.zeros_like(twists),
            # vec_sweep = π with frame_of_reference_delta="y_AFoR":
            # yB = +Y → xB = +X (chord direction). vec_sweep=π rotates the chord
            # 180° around span, placing LE at +X (into V_rel = [-ωr, Vw, 0]).
            # Panel normal = (TE-LE) × span = (-X) × (+Z) = +Y = upwind → AoA > 0.
            # This matches SHARPy's template_wt.py convention.
            math.pi * np.ones_like(chords),
            np.array([int(chord_panels)], dtype=int),
            np.zeros(blade.StructuralInformation.num_elem, dtype=int),
            "uniform",
            pitch_axes,
            np.arange(len(radii), dtype=int),
            _station_camber_lines(self.blade_aero, n_points_camber),
            first_twist=False,
        )

        # Pitch rotates around Z (the blade span axis). rotate_around_origin
        # updates both coordinates and frame_of_reference_delta, so the B-frame
        # (and hence chord direction) is correctly feathered for non-zero pitch.
        blade.StructuralInformation.rotate_around_origin(
            np.array([0.0, 0.0, 1.0]),
            -math.radians(pitch_deg),
        )

        rotor = blade.copy()
        for blade_index in range(1, int(self.blade_aero.n_blades)):
            blade_copy = blade.copy()
            blade_copy.StructuralInformation.rotate_around_origin(
                rotation_axis,
                blade_index * 2.0 * math.pi / float(self.blade_aero.n_blades),
            )
            rotor.assembly(blade_copy)
        rotor.remove_duplicated_points(1.0e-8, skip=[0])
        n_airfoils = int(rotor.AerodynamicInformation.airfoils.shape[0])
        rotor.AerodynamicInformation.polars = _build_sharpy_polars(
            self.blade_aero,
            n_airfoils=n_airfoils,
            polar_re=polar_re,
        )
        rotor.StructuralInformation.body_number *= 0
        return rotor

    def _build_simulation(
        self,
        *,
        case_name: str,
        route: Path,
        model: SharpyRotorModel,
        wind_speed: float,
        omega_rad_s: float,
        wind_direction: np.ndarray,
        rotation_axis: np.ndarray,
        freestream_direction: np.ndarray,
        dt: float,
        dphi: float,
        mstar: int,
        time_steps: int,
        solver_settings: Mapping[str, Any],
        rotor_node_count: int,
    ) -> Any:
        gc, _sharpy_main, _algebra = _load_sharpy_modules()

        sim = gc.SimulationInformation()
        sim.set_default_values()
        sim.solvers["SHARPy"]["flow"] = _flow_for_model(model)
        sim.solvers["SHARPy"]["case"] = case_name
        sim.solvers["SHARPy"]["route"] = str(route) + "/"
        sim.solvers["SHARPy"]["write_log"] = bool(solver_settings.get("write_log", False))
        sim.solvers["SHARPy"]["write_screen"] = "off"
        sim.set_variable_all_dicts("dt", dt)
        sim.set_variable_all_dicts("rho", self.rho)

        sim.solvers["SteadyVelocityField"]["u_inf"] = wind_speed
        sim.solvers["SteadyVelocityField"]["u_inf_direction"] = wind_direction
        sim.set_variable_all_dicts("velocity_field_input", sim.solvers["SteadyVelocityField"])

        sim.solvers["BeamLoader"]["unsteady"] = "on"
        sim.solvers["AerogridLoader"]["unsteady"] = "on"
        sim.solvers["AerogridLoader"]["mstar"] = mstar
        sim.solvers["AerogridLoader"]["aligned_grid"] = bool(
            solver_settings.get("aligned_grid", False)
        )
        sim.solvers["AerogridLoader"]["initial_align"] = bool(
            solver_settings.get("initial_align", False)
        )
        sim.solvers["AerogridLoader"]["freestream_dir"] = freestream_direction
        sim.solvers["AerogridLoader"]["wake_shape_generator"] = "HelicoidalWake"
        sim.solvers["AerogridLoader"]["wake_shape_generator_input"] = {
            "u_inf": wind_speed,
            "u_inf_direction": sim.solvers["SteadyVelocityField"]["u_inf_direction"],
            "rotation_velocity": omega_rad_s * rotation_axis,
            "dt": dt,
            "dphi1": dphi,
            "ndphi1": mstar,
            "r": 1.0,
            "dphimax": 10.0 * math.pi / 180.0,
        }

        sim.solvers["StaticUvlm"].update({
            "horseshoe": False,
            "num_cores": self.num_cores,
            "n_rollup": int(solver_settings.get("n_rollup", 0)),
            "rollup_dt": float(solver_settings.get("rollup_dt", dt)),
            "rollup_aic_refresh": int(solver_settings.get("rollup_aic_refresh", 1)),
            "rollup_tolerance": float(solver_settings.get("rollup_tolerance", 1.0e-8)),
            "cfl1": bool(solver_settings.get("cfl1", True)),
            "velocity_field_generator": "SteadyVelocityField",
            "velocity_field_input": sim.solvers["SteadyVelocityField"],
        })
        sim.solvers["StaticCoupled"]["structural_solver"] = "RigidDynamicPrescribedStep"
        sim.solvers["StaticCoupled"]["structural_solver_settings"] = sim.solvers[
            "RigidDynamicPrescribedStep"
        ]
        sim.solvers["StaticCoupled"]["aero_solver"] = "StaticUvlm"
        sim.solvers["StaticCoupled"]["aero_solver_settings"] = sim.solvers["StaticUvlm"]
        sim.solvers["StaticCoupled"]["tolerance"] = 1.0e-8
        sim.solvers["StaticCoupled"]["n_load_steps"] = 0
        sim.solvers["StaticCoupled"]["relaxation_factor"] = 0.0

        if model == "uvlm":
            # DynamicCoupled drives the time-marching UVLM loop that follows
            # the StaticCoupled initialisation.
            step_uvlm_settings = {
                "num_cores": self.num_cores,
                "n_time_steps": time_steps,
                "dt": dt,
                "convection_scheme": int(solver_settings.get("convection_scheme", 2)),
                "cfl1": bool(solver_settings.get("cfl1", True)),
                "velocity_field_generator": "SteadyVelocityField",
                "velocity_field_input": sim.solvers["SteadyVelocityField"],
                "rho": self.rho,
            }
            if "vortex_radius" in solver_settings:
                step_uvlm_settings["vortex_radius"] = float(solver_settings["vortex_radius"])
            if "vortex_radius_wake_ind" in solver_settings:
                step_uvlm_settings["vortex_radius_wake_ind"] = float(
                    solver_settings["vortex_radius_wake_ind"]
                )
            sim.solvers["DynamicCoupled"] = {
                "structural_solver": "RigidDynamicPrescribedStep",
                "structural_solver_settings": sim.solvers["RigidDynamicPrescribedStep"],
                "aero_solver": "StepUvlm",
                "aero_solver_settings": step_uvlm_settings,
                "n_time_steps": time_steps,
                "dt": dt,
                "minimum_steps": 1,
                "relaxation_steps": 150,
                "final_relaxation_factor": 0.0,
                "fsi_tolerance": 1.0e-4,
                "fsi_substeps": 100,
                "relaxation_factor": 0.0,
            }

        # Integrated aero loads. u_inf_dir aligned with the wind keeps the
        # post-processor's internal flow-rotation a no-op, so the reported
        # total_*_inertial_forces are the G-frame rotor force/moment about the
        # hub (moment reference defaults to the origin = rotation centre).
        sim.solvers["AeroForcesCalculator"]["u_inf_dir"] = [float(c) for c in wind_direction]
        sim.solvers["AeroForcesCalculator"]["screen_output"] = False
        sim.solvers["AeroForcesCalculator"]["coefficients"] = False
        sim.solvers["AeroForcesCalculator"]["write_text_file"] = bool(
            solver_settings.get("write_text_file", False)
        )
        sim.solvers["AeroForcesCalculator"]["text_file_name"] = "forces"

        sim.define_num_steps(time_steps)
        sim.with_forced_vel = True
        sim.for_vel = np.zeros((time_steps, 6), dtype=float)
        sim.for_acc = np.zeros((time_steps, 6), dtype=float)
        sim.for_vel[:, 3:6] = omega_rad_s * rotation_axis
        sim.dynamic_forces = np.zeros((time_steps, rotor_node_count, 6), dtype=float)
        return sim


def _build_aero_node_mask(blade_aero: BladeAero) -> np.ndarray:
    """All stations participate in the VLM aerodynamic lattice.

    Cylindrical root stations are modelled as flat-plate panels (zero camber),
    which ensures the VLM lattice starts at hub_radius and spans the full blade
    length — matching the structural coupling mesh used in FSI.

    The earlier behaviour of masking out cylindrical sections was reverted
    because it created a radial gap between the VLM root (~7.7 m for IEA-15MW)
    and the structural mesh root (hub_radius = 3.0 m), making mesh alignment
    in Paraview impossible and breaking the radial interpolation in the FSI
    adapter.
    """
    return np.ones(len(blade_aero.stations), dtype=bool)


def run_sharpy_numad_rotor_case(
    blade_aero: BladeAero,
    operating_point: RotorOperatingPoint,
    *,
    route: Path,
    model: SharpyRotorModel,
    rho: float = 1.225,
    solver_settings: Mapping[str, Any] | None = None,
) -> SharpyRotorRunResult:
    """Convenience wrapper for one SHARPy NuMAD rotor benchmark case."""

    adapter = NuMADSharpyRotorCaseAdapter(blade_aero=blade_aero, route=route, rho=rho)
    return adapter.run(operating_point, model=model, solver_settings=solver_settings)


def _load_sharpy_modules() -> tuple[Any, Any, Any]:
    try:
        import sharpy.sharpy_main as sharpy_main
        import sharpy.utils.algebra as algebra
        import sharpy.utils.cout_utils as cout
        import sharpy.utils.generate_cases as gc
    except ModuleNotFoundError as exc:
        raise ModuleNotFoundError(
            "SHARPy is required for NuMAD rotor benchmarks. Install aeroelast[sharpy] "
            "or make the local SHARPy source importable."
        ) from exc
    cout.cout_wrap.print_file = False
    return gc, sharpy_main, algebra


def _normalise_model(model: str) -> SharpyRotorModel:
    normalised = str(model).strip().lower()
    if normalised not in {"vlm", "uvlm"}:
        raise ValueError("model must be 'vlm' or 'uvlm'")
    return normalised  # type: ignore[return-value]


def _flow_for_model(model: SharpyRotorModel) -> list[str]:
    flow = ["BeamLoader", "AerogridLoader", "StaticCoupled"]
    if model == "uvlm":
        flow.append("DynamicCoupled")
    # AeroForcesCalculator integrates the surface loads with a frame-consistent
    # moment arm (forces and positions summed in the A frame, then rotated to G)
    # and adds the unsteady contribution, rather than re-deriving the integration
    # by hand. It must sit in the top-level flow because it always reads its own
    # global settings section (it ignores the custom_settings a postprocessor
    # would receive). Cleanup is therefore disabled below so this offline pass
    # never iterates over nulled-out timesteps.
    flow.append("AeroForcesCalculator")
    return flow


def _case_name(prefix: str, model: SharpyRotorModel, point_name: str) -> str:
    safe_point = "".join(ch if ch.isalnum() else "_" for ch in point_name)
    return f"{prefix}_{model}_{safe_point}"


def _beam_stiffness_db(num_elem: int) -> np.ndarray:
    stiffness = np.zeros((num_elem, 6, 6), dtype=float)
    base = np.diag([1.0e10, 1.0e9, 1.0e9, 1.0e9, 1.0e12, 1.0e12])
    stiffness[:] = base
    return stiffness


def _beam_mass_db(num_elem: int) -> np.ndarray:
    mass = np.zeros((num_elem, 6, 6), dtype=float)
    base = np.diag([1.0e3, 1.0e3, 1.0e3, 1.0e5, 1.0e5, 1.0e5])
    mass[:] = base
    return mass


def _station_camber_lines(blade_aero: BladeAero, n_points: int) -> np.ndarray:
    gc, _sharpy_main, _algebra = _load_sharpy_modules()
    camber = np.zeros((len(blade_aero.stations), n_points, 2), dtype=float)
    camber[:, :, 0] = np.linspace(0.0, 1.0, n_points)
    for index, station in enumerate(blade_aero.stations):
        coords = np.asarray(station.airfoil.coordinates, dtype=float)
        if coords.ndim == 2 and coords.shape[0] >= 4 and coords.shape[1] >= 2:
            try:
                camber[index, :, 0], camber[index, :, 1] = gc.get_airfoil_camber(
                    coords[:, 0],
                    coords[:, 1],
                    n_points,
                )
            except Exception:
                camber[index, :, 1] = 0.0
    return camber


def _build_sharpy_polars(
    blade_aero: BladeAero,
    *,
    n_airfoils: int,
    polar_re: float,
) -> np.ndarray | None:
    if n_airfoils <= 0 or not blade_aero.stations:
        return None

    station_tables: list[np.ndarray] = []
    for station in blade_aero.stations:
        if not station.airfoil.polars:
            return None

        polar = station.airfoil.get_polar(polar_re)
        alpha = np.asarray(polar.alpha, dtype=float).reshape(-1)
        cl = np.asarray(polar.cl, dtype=float).reshape(-1)
        cd = np.asarray(polar.cd, dtype=float).reshape(-1)
        cm = np.asarray(polar.cm, dtype=float).reshape(-1)

        if not (alpha.size == cl.size == cd.size == cm.size) or alpha.size < 2:
            return None

        order = np.argsort(alpha)
        alpha = alpha[order]
        cl = cl[order]
        cd = cd[order]
        cm = cm[order]

        alpha_unique, unique_idx = np.unique(alpha, return_index=True)
        if alpha_unique.size < 2:
            return None

        alpha_max = float(np.max(np.abs(alpha_unique)))
        if alpha_max <= 2.0 * math.pi + 1.0e-6:
            alpha_export = np.rad2deg(alpha_unique)
        else:
            alpha_export = alpha_unique

        station_tables.append(
            np.column_stack((
                alpha_export,
                cl[unique_idx],
                cd[unique_idx],
                cm[unique_idx],
            ))
        )

    alpha_common = station_tables[0][:, 0]
    polars = np.zeros((n_airfoils, alpha_common.size, 4), dtype=float)
    polars[:, :, 0] = alpha_common
    for iairfoil in range(n_airfoils):
        table = station_tables[iairfoil % len(station_tables)]
        polars[iairfoil, :, 1] = np.interp(alpha_common, table[:, 0], table[:, 1])
        polars[iairfoil, :, 2] = np.interp(alpha_common, table[:, 0], table[:, 2])
        polars[iairfoil, :, 3] = np.interp(alpha_common, table[:, 0], table[:, 3])

    return polars


def _extract_rotor_result(
    sharpy_output: Any,
    *,
    rho: float,
    wind_speed: float,
    omega_rad_s: float,
    wind_direction: np.ndarray,
    rotation_axis: np.ndarray,
    timestep_index: int = -1,
) -> tuple[RotorPerformanceMetrics, RotorAeroLoads]:
    _gc, _sharpy_main, algebra = _load_sharpy_modules()
    tstep = sharpy_output.structure.timestep_info[timestep_index]
    wind_direction = _as_unit_vector(wind_direction, name="wind_direction")
    rotation_axis = _as_unit_vector(rotation_axis, name="rotation_axis")
    cga = algebra.quat2rotation(tstep.quat)

    steady_nodal_b, unsteady_nodal_b, nodal_force_source = _nodal_aero_forces_b(tstep)
    nodal_b = steady_nodal_b + unsteady_nodal_b
    nodal_a = tstep.nodal_b_for_2_a_for(nodal_b, sharpy_output.structure)

    # SHARPy stores nodal aero forces in the structural B frame. Converting
    # through nodal_b_for_2_a_for keeps the dynamic UVLM path consistent with
    # AeroForcesCalculator before rotating the resultants to the inertial G frame.
    nodal_forces_afor = np.zeros((sharpy_output.structure.num_node, 3), dtype=float)
    nodal_moments_afor = np.zeros((sharpy_output.structure.num_node, 3), dtype=float)
    for node_index in range(sharpy_output.structure.num_node):
        nodal_forces_afor[node_index, :] = np.dot(cga, nodal_a[node_index, 0:3])
        nodal_moments_afor[node_index, :] = np.dot(cga, nodal_a[node_index, 3:6])

    exported_force = np.sum(nodal_forces_afor, axis=0)
    force_based_moment_a = np.sum(
        np.cross(np.asarray(tstep.pos, dtype=float), nodal_a[:, 0:3]), axis=0
    )
    exported_moment = np.dot(cga, force_based_moment_a)
    local_aero_moment = np.sum(nodal_moments_afor, axis=0)

    # Frame-consistent integrated rotor force/moment about the hub, built from
    # the per-node aero loads (B -> A -> G). SHARPy's AeroForcesCalculator
    # total_*_inertial_forces are intentionally NOT used here: in the
    # time-marched rotating (UVLM) case they come back in an inconsistent frame
    # (axial thrust shows up in-plane). The calculator still runs in the flow
    # because it populates postproc_node["aero_steady_forces"], the per-node
    # source consumed above. The "total" moment adds the force-times-arm moment
    # to the local aerodynamic couples.
    total_force = exported_force
    total_moment = exported_moment + local_aero_moment

    radial_vectors = tstep.pos - np.outer(np.dot(tstep.pos, rotation_axis), rotation_axis)
    rotor_radius = float(np.max(np.linalg.norm(radial_vectors, axis=1)))
    rotor_area = math.pi * rotor_radius**2
    dynamic_force = 0.5 * rho * wind_speed**2 * rotor_area
    dynamic_power = dynamic_force * wind_speed

    thrust_signed = float(np.dot(exported_force, wind_direction))
    torque_signed = float(np.dot(exported_moment, rotation_axis))
    torque_total_signed = float(np.dot(total_moment, rotation_axis))
    torque_local_signed = float(np.dot(local_aero_moment, rotation_axis))

    thrust = float(abs(thrust_signed))
    torque = float(abs(torque_signed))
    power = float(torque * abs(omega_rad_s))
    cp = power / dynamic_power if dynamic_power > 0.0 else math.nan
    ct = thrust / dynamic_force if dynamic_force > 0.0 else math.nan
    cq = torque / (dynamic_force * rotor_radius) if dynamic_force > 0.0 else math.nan

    loads = RotorAeroLoads(
        nodal_forces=np.asarray(nodal_forces_afor, dtype=float),
        integrated_force=exported_force,
        integrated_moment=exported_moment,
        metadata={
            "cp": cp,
            "ct": ct,
            "cq": cq,
            "thrust": thrust,
            "torque": torque,
            "power": power,
            "rotor_radius": rotor_radius,
            "thrust_signed": thrust_signed,
            "torque_signed": torque_signed,
            "integrated_force_total": total_force.tolist(),
            "integrated_moment_total": total_moment.tolist(),
            "integrated_local_moment": local_aero_moment.tolist(),
            "torque_total": float(abs(torque_total_signed)),
            "torque_total_signed": torque_total_signed,
            "torque_local": float(abs(torque_local_signed)),
            "torque_local_signed": torque_local_signed,
            "nodal_force_source": nodal_force_source,
            "wind_direction": wind_direction.tolist(),
            "rotation_axis": rotation_axis.tolist(),
        },
    )
    return RotorPerformanceMetrics(thrust, torque, power, cp, ct, cq), loads


def _nodal_aero_forces_b(tstep: Any) -> tuple[np.ndarray, np.ndarray, str]:
    postproc_node = getattr(tstep, "postproc_node", {}) or {}
    steady_post = postproc_node.get("aero_steady_forces")
    if steady_post is not None:
        steady_b = np.asarray(steady_post, dtype=float).copy()
        unsteady_b = np.zeros_like(steady_b)
        unsteady_post = postproc_node.get("aero_unsteady_forces")
        if unsteady_post is not None:
            unsteady_b += np.asarray(unsteady_post, dtype=float)
        return steady_b, unsteady_b, "postproc_node"

    steady_b = np.asarray(tstep.steady_applied_forces, dtype=float).copy()
    unsteady_b = np.zeros_like(steady_b)
    raw_unsteady = getattr(tstep, "unsteady_applied_forces", None)
    if raw_unsteady is not None:
        unsteady_b += np.asarray(raw_unsteady, dtype=float)
    return steady_b, unsteady_b, "steady_applied_forces"


def _as_unit_vector(value: Any, *, name: str) -> np.ndarray:
    vec = np.asarray(value, dtype=float).reshape(-1)
    if vec.size != 3:
        raise ValueError(f"{name} must contain exactly 3 components")
    norm = float(np.linalg.norm(vec))
    if norm <= 1.0e-12:
        raise ValueError(f"{name} must be non-zero")
    return vec / norm


def _nodes_by_blade(structure: Any, n_blades: int) -> list[np.ndarray]:
    nodes_by_blade = []
    first_node = 0
    elem_cursor = 0
    for blade_index in range(n_blades):
        mask = np.zeros((structure.num_node,), dtype=bool)
        while (
            elem_cursor < structure.num_elem and structure.beam_number[elem_cursor] <= blade_index
        ):
            elem_cursor += 1
        end_node = int(structure.connectivities[elem_cursor - 1, 1]) + 1
        mask[first_node:end_node] = True
        nodes_by_blade.append(mask)
        first_node = end_node
    return nodes_by_blade


def _node_widths(radii: np.ndarray) -> np.ndarray:
    widths = np.zeros_like(radii, dtype=float)
    if len(radii) == 1:
        widths[0] = 1.0
        return widths
    widths[0] = 0.5 * abs(radii[1] - radii[0])
    widths[-1] = 0.5 * abs(radii[-1] - radii[-2])
    for index in range(1, len(radii) - 1):
        widths[index] = 0.5 * abs(radii[index + 1] - radii[index - 1])
    return widths


def _panel_count(sharpy_output: Any) -> int:
    try:
        return int(sum(gamma.size for gamma in sharpy_output.aero.timestep_info[-1].gamma))
    except Exception:
        return 0


__all__ = [
    "NuMADSharpyRotorCaseAdapter",
    "SharpyRotorModel",
    "SharpyRotorRunResult",
    "run_sharpy_numad_rotor_case",
]
