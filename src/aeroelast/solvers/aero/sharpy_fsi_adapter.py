"""SHARPy quasi-static adapter for aeroelast rotor FSI participants.

This adapter is intentionally conservative for V1:

- VLM quasi-static solve (`StaticUvlm`)
- SHARPy case is built once and updated in-place with translated nodal motion
- aerodynamic loads are mapped back to full-rotor FE nodes
"""

from __future__ import annotations

import importlib
import logging
import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Mapping, cast

import numpy as np

from aeroelast.models.blade.aerodynamics import load_blade_aero

from .sharpy_numad import (
    NuMADSharpyRotorCaseAdapter,
    _as_unit_vector,
    _extract_rotor_result,
    _load_sharpy_modules,
    _nodes_by_blade,
)
from .types import AeroFSIRuntimeContext, RotorAeroLoads, RotorAeroState

logger = logging.getLogger(__name__)


@dataclass
class _BladeMapCache:
    fe_indices: np.ndarray
    fe_radii_ref: np.ndarray
    sharpy_indices: np.ndarray
    sharpy_radii_ref: np.ndarray
    nearest_fe_for_sharpy: np.ndarray
    nearest_sharpy_for_fe: np.ndarray


@dataclass
class _AdapterCache:
    ref_pos: np.ndarray
    blade_maps: tuple[_BladeMapCache, ...]


@dataclass
class SharpyRotorFSICaseAdapter:
    """Quasi-static SHARPy adapter for coupled rotor FSI."""

    blade_file: str
    route: Path
    case_prefix: str = "fsi_sharpy"
    chord_panels: int = 4
    hub_radius: float = 3.0
    n_points_camber: int = 100
    dphi_deg: float = 5.0
    wake_revolutions: float = 0.05
    polar_re: float = 1.0e7
    rho: float = 1.225
    num_cores: int = 1
    default_re: float = 1.0e7
    neuralfoil_model: str = "large"
    pitch_deg: float = 0.0
    rebuild_per_step: bool = False
    _cache_by_data_id: dict[int, _AdapterCache] = field(default_factory=dict, init=False)

    @classmethod
    def from_runtime_context(
        cls,
        runtime_context: AeroFSIRuntimeContext,
        backend_cfg: Mapping[str, Any],
    ) -> "SharpyRotorFSICaseAdapter":
        blade_file = runtime_context.blade_file
        if not blade_file:
            raise ValueError("SHARPy FSI adapter requires a blade_file in aero.blade_file")
        route = Path(str(backend_cfg.get("route", "output/diagnostics/sharpy_fsi")))
        return cls(
            blade_file=str(blade_file),
            route=route,
            case_prefix=str(backend_cfg.get("case_prefix", "fsi_sharpy")),
            chord_panels=int(backend_cfg.get("chord_panels", 4)),
            hub_radius=float(
                backend_cfg.get("hub_radius", runtime_context.rotor.get("hub_radius", 3.0))
            ),
            n_points_camber=int(backend_cfg.get("n_points_camber", 100)),
            dphi_deg=float(backend_cfg.get("dphi_deg", 5.0)),
            wake_revolutions=float(backend_cfg.get("wake_revolutions", 0.05)),
            polar_re=float(backend_cfg.get("polar_re", 1.0e7)),
            rho=float(backend_cfg.get("air_density", 1.225)),
            num_cores=int(backend_cfg.get("num_cores", 1)),
            default_re=float(backend_cfg.get("default_re", 1.0e7)),
            neuralfoil_model=str(backend_cfg.get("neuralfoil_model", "large")),
            pitch_deg=float(backend_cfg.get("pitch_deg", backend_cfg.get("pitch", 0.0))),
            rebuild_per_step=bool(backend_cfg.get("rebuild_per_step", False)),
        )

    def build_data(
        self,
        context: AeroFSIRuntimeContext,
        state: RotorAeroState,
        settings: Mapping[str, Any],
    ) -> Any:
        logger.info("[AERO-SHARPY] Adapter build_data: loading blade aero from %s", self.blade_file)
        blade_aero = load_blade_aero(
            self.blade_file,
            default_re=self.default_re,
            neuralfoil_model=self.neuralfoil_model,
            hub_radius=float(context.rotor.get("hub_radius", self.hub_radius)),
            n_blades=int(context.rotor.get("n_blades", 3)),
        )
        logger.info("[AERO-SHARPY] Adapter build_data: blade aero loaded")

        numad = NuMADSharpyRotorCaseAdapter(
            blade_aero=blade_aero,
            route=self.route,
            case_prefix=self.case_prefix,
            chord_panels=int(settings.get("chord_panels", self.chord_panels)),
            hub_radius=float(settings.get("hub_radius", self.hub_radius)),
            n_points_camber=int(settings.get("n_points_camber", self.n_points_camber)),
            wake_revolutions=float(settings.get("wake_revolutions", self.wake_revolutions)),
            dphi_deg=float(settings.get("dphi_deg", self.dphi_deg)),
            time_steps=1,
            simulation_revolutions=1.0,
            polar_re=float(settings.get("polar_re", self.polar_re)),
            rho=float(settings.get("rho", self.rho)),
            num_cores=int(settings.get("num_cores", self.num_cores)),
        )

        omega = float(state.omega_rad_s)
        wind = np.asarray(state.wind_velocity, dtype=float).reshape(3)
        wind_speed = float(np.linalg.norm(wind))
        wind_dir = (
            np.array([0.0, 0.0, 1.0], dtype=float) if wind_speed < 1.0e-12 else wind / wind_speed
        )
        rotation_axis = _as_unit_vector(
            context.rotor.get("rotation_axis", [0.0, 1.0, 0.0]),
            name="rotation_axis",
        )
        dphi = math.radians(float(settings.get("dphi_deg", self.dphi_deg)))
        dt = float(state.dt) if state.dt > 0.0 else (dphi / max(abs(omega), 1.0e-9))
        mstar = max(1, int(settings.get("mstar", self.wake_revolutions * 2.0 * math.pi / dphi)))

        case_name = f"{self.case_prefix}_vlm_fsi"
        route = self.route.expanduser().resolve()
        route.mkdir(parents=True, exist_ok=True)
        logger.info(
            "[AERO-SHARPY] Adapter build_data: preparing SHARPy case '%s' in %s",
            case_name,
            route,
        )

        rotor = numad._build_rotor(
            pitch_deg=float(settings.get("pitch_deg", settings.get("pitch", self.pitch_deg))),
            chord_panels=int(settings.get("chord_panels", self.chord_panels)),
            hub_radius=float(settings.get("hub_radius", self.hub_radius)),
            n_points_camber=int(settings.get("n_points_camber", self.n_points_camber)),
            rotation_axis=rotation_axis,
            polar_re=float(settings.get("polar_re", self.polar_re)),
        )
        sim = numad._build_simulation(
            case_name=case_name,
            route=route,
            model="vlm",
            wind_speed=wind_speed,
            omega_rad_s=omega,
            wind_direction=wind_dir,
            rotation_axis=rotation_axis,
            freestream_direction=wind_dir,
            dt=dt,
            dphi=dphi,
            mstar=mstar,
            time_steps=1,
            solver_settings=settings,
            rotor_node_count=rotor.StructuralInformation.num_node,
        )
        logger.info("[AERO-SHARPY] Adapter build_data: rotor and simulation assembled")

        gc, _sharpy_main, _algebra = _load_sharpy_modules()
        gc.clean_test_files(sim.solvers["SHARPy"]["route"], sim.solvers["SHARPy"]["case"])
        rotor.generate_h5_files(sim.solvers["SHARPy"]["route"], sim.solvers["SHARPy"]["case"])
        sim.generate_solver_file()
        sim.generate_dyn_file(1)
        logger.info("[AERO-SHARPY] Adapter build_data: SHARPy input files generated")

        data = self._build_presharpy_data(route / f"{case_name}.sharpy")
        logger.info("[AERO-SHARPY] Adapter build_data: PreSharpy data loaded")
        self._prime_sharpy_data(data)
        logger.info("[AERO-SHARPY] Adapter build_data: Beam/Aerogrid primed")
        cache = self._build_cache(data, context)
        self._cache_by_data_id[id(data)] = cache
        logger.info(
            "[AERO-SHARPY] Adapter build_data: cache ready (blades=%d)",
            len(cache.blade_maps),
        )
        return data

    def _build_presharpy_data(self, case_file: Path) -> Any:
        try:
            import configobj
        except ModuleNotFoundError as exc:  # pragma: no cover - runtime dependency
            raise ModuleNotFoundError(
                "configobj is required to load SHARPy .sharpy settings files"
            ) from exc

        presharpy_module = importlib.import_module("sharpy.presharpy.presharpy")
        pre_cls = presharpy_module.PreSharpy
        settings = configobj.ConfigObj(str(case_file))
        return pre_cls(settings)

    def _prime_sharpy_data(self, data: Any) -> None:
        # StaticUvlm expects BeamLoader and AerogridLoader outputs to be available.
        beam_mod = importlib.import_module("sharpy.solvers.beamloader")
        aero_mod = importlib.import_module("sharpy.solvers.aerogridloader")
        beam_solver = beam_mod.BeamLoader()
        aero_solver = aero_mod.AerogridLoader()

        beam_solver.initialise(data)
        data = beam_solver.run()
        aero_solver.initialise(data)
        aero_solver.run()

    def export_aero_mesh(self, data: Any, output_file: str | Path) -> None:
        """Write the SHARPy aerodynamic lattice, not the 3D coupling mesh, to VTU."""
        try:
            import meshio
        except ModuleNotFoundError as exc:  # pragma: no cover - optional dependency at runtime
            raise ModuleNotFoundError(
                "meshio is required to export the SHARPy aerodynamic lattice"
            ) from exc

        aero = getattr(data, "aero", None)
        timestep_info = getattr(aero, "timestep_info", None)
        if aero is None or not timestep_info:
            raise RuntimeError("SHARPy aerodynamic timestep_info is not available for export")

        tstep = timestep_info[-1]
        zeta_blocks = getattr(tstep, "zeta", None)
        if zeta_blocks is None or len(zeta_blocks) == 0:
            raise RuntimeError("SHARPy aerodynamic lattice is missing zeta coordinates")

        points: list[np.ndarray] = []
        triangles: list[list[int]] = []
        surface_ids: list[int] = []

        for surface_index, surface_zeta in enumerate(zeta_blocks, start=1):
            zeta = np.asarray(surface_zeta, dtype=float)
            if zeta.ndim != 3 or zeta.shape[0] != 3:
                raise RuntimeError(
                    "Expected SHARPy zeta arrays with shape (3, chord_nodes, span_nodes)"
                )

            _, chord_nodes, span_nodes = zeta.shape
            point_offset = len(points)
            for chord_index in range(chord_nodes):
                for span_index in range(span_nodes):
                    points.append(zeta[:, chord_index, span_index].copy())

            for chord_index in range(chord_nodes - 1):
                row_start = point_offset + chord_index * span_nodes
                next_row_start = point_offset + (chord_index + 1) * span_nodes
                for span_index in range(span_nodes - 1):
                    p00 = zeta[:, chord_index, span_index]
                    p10 = zeta[:, chord_index + 1, span_index]
                    p11 = zeta[:, chord_index + 1, span_index + 1]
                    p01 = zeta[:, chord_index, span_index + 1]

                    a = row_start + span_index
                    b = next_row_start + span_index
                    c = next_row_start + span_index + 1
                    d = row_start + span_index + 1

                    # Some SHARPy lattice panels are slightly warped. Exporting
                    # triangles with the better diagonal avoids bow-tie quads
                    # that some viewers fail to render.
                    n1 = np.cross(p10 - p00, p11 - p00)
                    n2 = np.cross(p11 - p00, p01 - p00)
                    n1_norm = float(np.linalg.norm(n1))
                    n2_norm = float(np.linalg.norm(n2))
                    diag_02_score = -1.0
                    if n1_norm > 1.0e-12 and n2_norm > 1.0e-12:
                        diag_02_score = float(np.dot(n1 / n1_norm, n2 / n2_norm))

                    m1 = np.cross(p10 - p00, p01 - p00)
                    m2 = np.cross(p01 - p10, p11 - p10)
                    m1_norm = float(np.linalg.norm(m1))
                    m2_norm = float(np.linalg.norm(m2))
                    diag_13_score = -1.0
                    if m1_norm > 1.0e-12 and m2_norm > 1.0e-12:
                        diag_13_score = float(np.dot(m1 / m1_norm, m2 / m2_norm))

                    if diag_13_score > diag_02_score:
                        triangles.append([a, b, d])
                        triangles.append([b, c, d])
                    else:
                        triangles.append([a, b, c])
                        triangles.append([a, c, d])
                    surface_ids.extend([surface_index, surface_index])

        if not points:
            raise RuntimeError("SHARPy aerodynamic lattice export produced no points")

        output_path = Path(output_file)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        point_array = np.asarray(points, dtype=float)

        if triangles:
            cells = [("triangle", np.asarray(triangles, dtype=int))]
            cell_data = {"surface_id": [np.asarray(surface_ids, dtype=int)]}
        else:
            cells = [("vertex", np.arange(point_array.shape[0], dtype=int).reshape(-1, 1))]
            cell_data = None

        mesh = meshio.Mesh(
            points=point_array,
            cells=cast(Any, cells),
            cell_data=cast(Any, cell_data),
        )
        # ParaView has been more reliable with ASCII XML here than with the
        # default compressed binary VTU emitted by meshio.
        mesh.write(output_path, binary=False)

        preview_path: Path | None = None
        if output_path.suffix.lower() != ".vtk":
            preview_path = output_path.with_suffix(".vtk")
            mesh.write(preview_path, binary=False)

        logger.info(
            "[AERO-SHARPY] Wrote aerodynamic lattice mesh to %s (surfaces=%d triangles=%d nodes=%d)",
            output_path,
            len(zeta_blocks),
            len(triangles),
            point_array.shape[0],
        )
        if preview_path is not None:
            logger.info("[AERO-SHARPY] Wrote ParaView preview mesh to %s", preview_path)

    def _build_cache(self, data: Any, context: AeroFSIRuntimeContext) -> _AdapterCache:
        if context.geometry is None:
            raise ValueError("SHARPy FSI adapter requires runtime_context.geometry")

        tstep = data.structure.timestep_info[-1]
        ref_pos = np.asarray(tstep.pos, dtype=float).copy()
        n_blades = context.geometry.n_blades
        sharpy_nodes_by_blade = _nodes_by_blade(data.structure, n_blades)
        rotation_axis = context.geometry.rotation_axis.reshape(3)
        rotation_center = context.geometry.rotation_center.reshape(3)

        blade_maps: list[_BladeMapCache] = []
        rotor_coords = context.geometry.rotor_mesh.coords_array
        radial_ref = self._radial_distance(rotor_coords, rotation_axis, rotation_center)
        sharpy_radial = self._radial_distance(ref_pos, rotation_axis, rotation_center)

        for blade, sharpy_indices in zip(context.geometry.blades, sharpy_nodes_by_blade):
            fe_indices = np.asarray(blade.parent_node_indices, dtype=int)
            fe_radii = radial_ref[fe_indices]
            sharpy_indices = np.asarray(sharpy_indices, dtype=int)
            sharpy_radii = sharpy_radial[sharpy_indices]

            nearest_fe_for_sharpy = self._nearest_source_indices(sharpy_radii, fe_radii)
            nearest_sharpy_for_fe = self._nearest_source_indices(fe_radii, sharpy_radii)

            blade_maps.append(
                _BladeMapCache(
                    fe_indices=fe_indices,
                    fe_radii_ref=fe_radii,
                    sharpy_indices=sharpy_indices,
                    sharpy_radii_ref=sharpy_radii,
                    nearest_fe_for_sharpy=nearest_fe_for_sharpy,
                    nearest_sharpy_for_fe=nearest_sharpy_for_fe,
                )
            )

        return _AdapterCache(ref_pos=ref_pos, blade_maps=tuple(blade_maps))

    @staticmethod
    def _nearest_source_indices(target_r: np.ndarray, source_r: np.ndarray) -> np.ndarray:
        if source_r.size == 0:
            return np.zeros(target_r.shape, dtype=int)
        if target_r.size == 0:
            return np.zeros((0,), dtype=int)
        d = np.abs(target_r[:, None] - source_r[None, :])
        return np.argmin(d, axis=1)

    @staticmethod
    def _radial_distance(coords: np.ndarray, axis: np.ndarray, center: np.ndarray) -> np.ndarray:
        centered = coords - center.reshape(1, 3)
        axial = np.outer(np.dot(centered, axis), axis)
        radial = centered - axial
        return np.linalg.norm(radial, axis=1)

    def update_state(
        self,
        data: Any,
        context: AeroFSIRuntimeContext,
        state: RotorAeroState,
    ) -> None:
        if state.nodal_displacements is None:
            return

        cache = self._cache_by_data_id.get(id(data))
        if cache is None:
            cache = self._build_cache(data, context)
            self._cache_by_data_id[id(data)] = cache

        tstep = data.structure.timestep_info[-1]
        pos = np.asarray(tstep.pos, dtype=float).copy()
        displacements = np.asarray(state.nodal_displacements, dtype=float)

        for blade_map in cache.blade_maps:
            fe_disp = displacements[blade_map.fe_indices]
            if fe_disp.size == 0:
                continue
            mapped = fe_disp[blade_map.nearest_fe_for_sharpy]
            pos[blade_map.sharpy_indices] = cache.ref_pos[blade_map.sharpy_indices] + mapped

        tstep.pos[:] = pos
        if hasattr(data.structure, "pos"):
            data.structure.pos[:] = pos

    @staticmethod
    def _ensure_forces_mapped(data: Any) -> None:
        """Map VLM panel forces to structural nodes when StaticUvlm didn't do it.

        ``StaticUvlm.run()`` fills ``data.aero.timestep_info[ts].forces`` via
        ``uvlmlib.vlm_solver`` but never calls ``aero2struct_force_mapping``,
        so ``data.structure.timestep_info[ts].steady_applied_forces`` stays
        zero.  We call the mapping here before reading back those forces.
        """
        try:
            mapping_mod = importlib.import_module("sharpy.aero.utils.mapping")
            algebra_mod = importlib.import_module("sharpy.utils.algebra")
        except ImportError:
            logger.warning(
                "[AERO-SHARPY] Cannot import sharpy mapping modules — forces may be zero"
            )
            return

        ts = data.ts
        aero_tstep = data.aero.timestep_info[ts]
        str_tstep = data.structure.timestep_info[ts]

        if not aero_tstep.forces:
            logger.warning("[AERO-SHARPY] No panel forces in aero timestep — skipping mapping")
            return

        total_aero = sum(float(np.sum(np.abs(f))) for f in aero_tstep.forces)
        if total_aero < 1.0e-30:
            logger.warning(
                "[AERO-SHARPY] VLM panel forces are all zero (total=%.3e) — no mapping needed",
                total_aero,
            )
            return

        cag = algebra_mod.quat2rotation(str_tstep.quat).T
        struct_forces = mapping_mod.aero2struct_force_mapping(
            aero_tstep.forces,
            data.aero.struct2aero_mapping,
            aero_tstep.zeta,
            str_tstep.pos,
            str_tstep.psi,
            None,  # master (unused by SHARPy mapping)
            data.structure.connectivities,
            cag=cag,
        )
        str_tstep.steady_applied_forces[:] = struct_forces
        logger.debug(
            "[AERO-SHARPY] Mapped panel forces → structural nodes: "
            "force_norm=%.3e moment_norm=%.3e",
            float(np.linalg.norm(struct_forces[:, 0:3])),
            float(np.linalg.norm(struct_forces[:, 3:6])),
        )

    def extract_loads(
        self,
        data: Any,
        context: AeroFSIRuntimeContext,
        state: RotorAeroState,
    ) -> RotorAeroLoads:
        # StaticUvlm.run() does NOT call aero2struct_force_mapping — fill the
        # structural timestep from panel forces before _extract_rotor_result reads it.
        self._ensure_forces_mapped(data)
        wind = np.asarray(state.wind_velocity, dtype=float).reshape(3)
        wind_speed = float(np.linalg.norm(wind))
        wind_dir = (
            np.array([0.0, 0.0, 1.0], dtype=float) if wind_speed < 1.0e-12 else wind / wind_speed
        )
        rotation_axis = _as_unit_vector(
            context.rotor.get("rotation_axis", [0.0, 1.0, 0.0]),
            name="rotation_axis",
        )

        _metrics, sharpy_loads = _extract_rotor_result(
            data,
            rho=self.rho,
            wind_speed=max(wind_speed, 1.0e-9),
            omega_rad_s=float(state.omega_rad_s),
            wind_direction=wind_dir,
            rotation_axis=rotation_axis,
        )

        if context.geometry is None:
            return sharpy_loads

        n_fe = context.geometry.rotor_mesh.node_count
        if sharpy_loads.nodal_forces.shape[0] == n_fe:
            return sharpy_loads

        cache = self._cache_by_data_id.get(id(data))
        if cache is None:
            cache = self._build_cache(data, context)
            self._cache_by_data_id[id(data)] = cache

        fe_forces = np.zeros((n_fe, 3), dtype=float)
        for blade_map in cache.blade_maps:
            if blade_map.sharpy_indices.size == 0 or blade_map.fe_indices.size == 0:
                continue
            blade_sharpy = sharpy_loads.nodal_forces[blade_map.sharpy_indices]
            # Normalize by the number of FE nodes mapped to each SHARPy beam node so
            # that the TOTAL force applied to the FE mesh equals the SHARPy total force.
            # Without this, every FE coupling node would carry the full lumped force of
            # its nearest structural node, overcounting by N_fe/N_sharpy (~150x).
            counts = np.bincount(
                blade_map.nearest_sharpy_for_fe,
                minlength=len(blade_map.sharpy_indices),
            )
            weights = np.maximum(counts[blade_map.nearest_sharpy_for_fe], 1)
            mapped = blade_sharpy[blade_map.nearest_sharpy_for_fe] / weights[:, np.newaxis]
            fe_forces[blade_map.fe_indices] = mapped

        center = context.geometry.rotation_center.reshape(1, 3)
        arms = context.geometry.rotor_mesh.coords_array - center
        integrated_force = np.sum(fe_forces, axis=0)
        integrated_moment = np.cross(arms, fe_forces).sum(axis=0)

        metadata = dict(sharpy_loads.metadata)
        metadata["force_mapping"] = "nearest_radial_by_blade_normalized"
        metadata["source_node_count"] = int(sharpy_loads.nodal_forces.shape[0])
        metadata["target_node_count"] = int(n_fe)

        return RotorAeroLoads(
            nodal_forces=fe_forces,
            integrated_force=integrated_force,
            integrated_moment=integrated_moment,
            metadata=metadata,
        )

    def build_run_kwargs(
        self,
        data: Any,
        context: AeroFSIRuntimeContext,
        state: RotorAeroState,
    ) -> dict[str, Any]:
        return {}

    def build_initialise_settings(
        self,
        data: Any,
        context: AeroFSIRuntimeContext,
        state: RotorAeroState,
    ) -> dict[str, Any]:
        init_settings: dict[str, Any] = {}
        settings = getattr(data, "settings", None)
        if isinstance(settings, Mapping):
            static_cfg = settings.get("StaticUvlm")
            if isinstance(static_cfg, Mapping):
                logger.info(
                    "[AERO-SHARPY] Using StaticUvlm section for solver initialization settings"
                )
                init_settings = dict(static_cfg)
            else:
                static_coupled = settings.get("StaticCoupled")
                if isinstance(static_coupled, Mapping):
                    aero_solver_cfg = static_coupled.get("aero_solver_settings")
                    if isinstance(aero_solver_cfg, Mapping):
                        logger.info(
                            "[AERO-SHARPY] Using StaticCoupled.aero_solver_settings for initialization"
                        )
                        init_settings = dict(aero_solver_cfg)
                else:
                    logger.warning(
                        "[AERO-SHARPY] No SHARPy init settings found in generated case data"
                    )

        # Override rbm_vel_g with the physical rotational velocity so the VLM
        # boundary condition includes the tangential (omega × r) component.
        # Without this, the solver sees only axial wind and produces incorrect AoA.
        omega = float(state.omega_rad_s)
        rotation_axis = _as_unit_vector(
            context.rotor.get("rotation_axis", [0.0, 1.0, 0.0]),
            name="rotation_axis",
        )
        # Negate omega: SHARPy adds rbm_vel_g to U_inf in K-J force (V_eff = U_inf + v_rbm).
        # The blade moves at v_blade = +ω×r so the apparent wind is U_inf - v_blade.
        # To encode apparent wind = U_inf + v_rbm correctly: v_rbm = -ω×r → negate omega.
        rbm_rot = (-omega * rotation_axis).tolist()
        init_settings["rbm_vel_g"] = [0.0, 0.0, 0.0] + rbm_rot
        logger.info(
            "[AERO-SHARPY] Setting rbm_vel_g=[0,0,0,%.4f,%.4f,%.4f] for omega=%.4f rad/s "
            "(omega negated: SHARPy adds rbm to U_inf, so v_rbm=-omega*r gives correct apparent wind)",
            *rbm_rot,
            omega,
        )
        return init_settings

    def capture_state(self, data: Any) -> dict[str, Any]:
        if data is None or not hasattr(data, "structure"):
            return {}
        tstep = data.structure.timestep_info[-1]
        return {
            "pos": np.asarray(tstep.pos, dtype=float).copy(),
        }

    def restore_state(self, data: Any, checkpoint_state: Mapping[str, Any] | None) -> None:
        if data is None or not hasattr(data, "structure"):
            return
        if not checkpoint_state:
            return
        pos = checkpoint_state.get("pos")
        if pos is None:
            return
        tstep = data.structure.timestep_info[-1]
        tstep.pos[:] = np.asarray(pos, dtype=float)
        if hasattr(data.structure, "pos"):
            data.structure.pos[:] = np.asarray(pos, dtype=float)


def build_quasi_static_sharpy_fsi_adapter(
    runtime_context: AeroFSIRuntimeContext,
    backend_cfg: Mapping[str, Any],
) -> SharpyRotorFSICaseAdapter:
    return SharpyRotorFSICaseAdapter.from_runtime_context(runtime_context, backend_cfg)
