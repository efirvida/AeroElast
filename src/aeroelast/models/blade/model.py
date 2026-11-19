from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Tuple, Union
import warnings

import numpy as np

from aeroelast.core.material import Material, OrthotropicMaterial
from aeroelast.core.mesh import MeshModel
from aeroelast.core.mesh.generators import BladeMesh, RotorMesh
from aeroelast.elements import ElementFamily


class Blade:
    """
    Wind turbine blade model with mesh generation and visualization capabilities.

    This class wraps the BladeMesh generator and provides additional
    functionality for blade visualization and geometry inspection.

    include_webs : bool
        When *False* the shear webs are omitted and only the blade outer shell
        (outer-mold-line surface) is meshed.  Use this for CFD wall surfaces.

    Parameters
    ----------
    blade_yaml : str
        Path to the windIO 2.x blade YAML (tests/IEA15MW/reference/upstream/IEA-15-240-RWT.yaml)
    element_size : float, optional
        Uniform target element size [m]. Default None: the reference mesh (chord_divisions 20, span_aspect 1).
    n_samples : int, optional
        NuMAD-only; any value other than None makes ``generate_mesh`` raise ValueError (BladeMesh builds with gmsh).
    include_webs : bool, optional
        Set *False* to mesh only the outer shell (default: True)
    skin_reference : {"mid", "oml"}, optional
        Surface the shell is built on.  ``"mid"`` (default) is the mid-surface every structural validation anchor
        (S-1, S-2, S-4, S-7) is measured on; ``"oml"`` keeps the nodes on the outer mould line (CFD wall surfaces).
    """

    def __init__(
        self,
        blade_yaml: str,
        element_size: Optional[float] = None,
        n_samples: Optional[int] = None,
        include_webs: bool = True,
        skin_reference: str = "mid",
    ) -> None:
        if skin_reference not in ("oml", "mid"):
            raise ValueError(f"skin_reference must be 'oml' or 'mid', got {skin_reference!r}")
        self.skin_reference = skin_reference
        self.yaml_file = blade_yaml
        self.element_size = element_size
        self.n_samples = n_samples
        self.include_webs = include_webs

        self._mesh_generator: Optional[BladeMesh] = None
        self.mesh: Optional[MeshModel] = None
        self._numad_mesh: Dict = {}
        self._numad_blade = None

    def generate_mesh(self, renumber: Optional[str] = None):
        """
        Generate the blade mesh from the YAML definition.

        Parameters
        ----------
        renumber : str, optional
            Renumbering algorithm to apply after mesh generation.
            - None (default): No renumbering
            - "simple": Direct index assignment
            - "rcm": Reverse Cuthill-McKee for bandwidth reduction
        """
        options = {} if self.n_samples is None else {"n_samples": self.n_samples}  # NuMAD-only: BladeMesh rejects it
        self._mesh_generator = BladeMesh(
            yaml_file=self.yaml_file,
            element_size=self.element_size,
            include_webs=self.include_webs,
            skin_reference=self.skin_reference,
            **options,
        )
        self.mesh = self._mesh_generator.generate(renumber=renumber, verbose=True)
        self._numad_mesh = self._mesh_generator.numad_mesh_data
        self._numad_blade = self._mesh_generator.numad_blade

    def view(self) -> None:
        """Visualize the blade mesh."""
        if self.mesh is None:
            raise RuntimeError("Mesh has not been generated yet. Call generate_mesh() first.")
        self.mesh.view()

    @property
    def element_thickness_map(self) -> Dict[int, float]:
        """Compute total shell thickness per element from section layup data.

        .. deprecated::
            Use :meth:`get_element_properties` instead, which returns full
            material definitions (including laminate stacking) per element set.

        Returns
        -------
        dict[int, float]
            Mapping of element ID to total laminate thickness.
        """
        warnings.warn(
            "element_thickness_map is deprecated. Use get_element_properties() instead.",
            DeprecationWarning,
            stacklevel=2,
        )
        if self.mesh is None:
            raise RuntimeError("Mesh has not been generated yet. Call generate_mesh() first.")

        return {
            elem.id: elem.thickness for elem in self.mesh.elements if elem.thickness is not None
        }

    def get_element_properties(self) -> dict:
        """Build shell property definitions per element set from section layup data.

        Returns Rust-native types (``_aeroelast.Laminate`` or isotropic
        dict) directly — no intermediate Python ``Laminate``/``Ply`` objects.

        Returns
        -------
        dict[str, _aeroelast.Laminate | dict]
            Mapping of element-set name to the corresponding property.
            Suitable for passing directly to ``PyMeshAssembler.from_model()``.

        Raises
        ------
        RuntimeError
            If the mesh has not been generated yet.
        KeyError
            If a layer references a material name not present in the
            numad material database.
        """
        if self.mesh is None or not self._numad_mesh:
            raise RuntimeError("Mesh has not been generated yet. Call generate_mesh() first.")

        return build_rust_properties(self._numad_mesh)

    def show_plots(self) -> None:
        """Display blade geometry plots."""
        if self._numad_blade is None:
            raise RuntimeError("Mesh has not been generated yet. Call generate_mesh() first.")
        import matplotlib.pyplot as plt

        from aeroelast.core.viewer import BladeGeometryVisualizer

        visualizer = BladeGeometryVisualizer(self._numad_blade)
        visualizer.plot_airfoil_type_distribution()
        visualizer.plot_chord_distribution()
        visualizer.plot_pitch_axis_position()
        visualizer.plot_thickness_to_chord_ratio()
        visualizer.plot_twist_distribution()
        plt.show()

    def write_mesh(self, filename: str, **kwargs) -> None:
        """Write the blade mesh to a file."""
        if self.mesh is None:
            raise RuntimeError("Mesh has not been generated yet. Call generate_mesh() first.")
        self.mesh.write_mesh(filename, **kwargs)


class Rotor:
    """
    Wind turbine rotor model with mesh generation capabilities.

    This class wraps the RotorMesh generator and provides additional
    functionality for rotor visualization.

    Parameters
    ----------
    blade_yaml : str
        Path to the windIO 2.x blade YAML (tests/IEA15MW/reference/upstream/IEA-15-240-RWT.yaml)
    n_blades : int
        Number of blades in the rotor
    hub_radius : float, optional
        Radial distance from rotation axis to blade root
    element_size : float, optional
        Uniform target element size [m]. Default None: the reference mesh (chord_divisions 20, span_aspect 1).
    n_samples : int, optional
        NuMAD-only; any value other than None makes ``generate_mesh`` raise ValueError (BladeMesh builds with gmsh).
    """

    def __init__(
        self,
        blade_yaml: str,
        n_blades: int,
        hub_radius: Optional[float] = None,
        element_size: Optional[float] = None,
        n_samples: Optional[int] = None,
        include_webs: bool = True,
    ) -> None:
        self.blade_yaml = blade_yaml
        self.n_blades = n_blades
        self.hub_radius = hub_radius
        self.element_size = element_size
        self.n_samples = n_samples
        self.include_webs = include_webs

        self._mesh_generator: Optional[RotorMesh] = None
        self.mesh: Optional[MeshModel] = None

    def generate_mesh(self, renumber: Optional[str] = None):
        """
        Generate the rotor mesh.

        Parameters
        ----------
        renumber : str, optional
            Renumbering algorithm to apply after mesh generation.
        """
        self._mesh_generator = RotorMesh(
            yaml_file=self.blade_yaml,
            n_blades=self.n_blades,
            hub_radius=self.hub_radius,
            element_size=self.element_size,
            include_webs=self.include_webs,
            # NuMAD-only: RotorMesh/BladeMesh reject n_samples (S7b); forwarded only when given.
            **({} if self.n_samples is None else {"n_samples": self.n_samples}),
        )
        self.mesh = self._mesh_generator.generate(renumber=renumber, verbose=True)

    def view(self) -> None:
        """Visualize the rotor mesh."""
        if self.mesh is None:
            raise RuntimeError("Mesh has not been generated yet. Call generate_mesh() first.")
        self.mesh.view()

    def write_mesh(self, filename: str, **kwargs) -> None:
        """Write the rotor mesh to a file."""
        if self.mesh is None:
            raise RuntimeError("Mesh has not been generated yet. Call generate_mesh() first.")
        self.mesh.write_mesh(filename, **kwargs)


def _windio_to_element_order(E, G, nu):
    """Reorder windIO out-of-plane constants to the element's order.

    windIO and NuMAD's material dict list ``G = [G12, G13, G23]`` and ``nu = [nu12, nu13, nu23]``
    (1 = fibre). ``OrthotropicMaterial`` (Python and Rust) takes ``(G12, G23, G31)`` and
    ``(nu12, nu23, nu31)``, with ``nu31 = nu13 * E3 / E1`` by reciprocity.
    """
    g12, g13, g23 = G
    nu12, nu13, nu23 = nu
    return (g12, g23, g13), (nu12, nu23, nu13 * E[2] / E[0])


def material_factory(
    material_data: Dict[str, Union[Dict[str, float], Dict[str, Tuple[float, float, float]]]],
) -> Union[Material, OrthotropicMaterial]:
    """
    Factory.s function to create a material object (either Isotropic or Orthotropic)
    based on .sthe provided material data.

    Parameters
    ----------
    material_data : dict
        A dictionary containing material properties. It should include the keys:
        - 'name' : str
            The name of the material.
        - 'density' : float
            The density of the material (kg/m^3).
        - 'elastic' : dict
            A dictionary containing the elastic properties:
            - 'E' : float or list of floats
                Young's Modulus (single value for isotropic, list for orthotropic).
            - 'nu' : float or list of floats
                Poisson's ratio (single value for isotropic, list for orthotropic).
            - 'G' : list of floats
                Shear Modulus (required for orthotropic materials only).

    Returns
    -------
    Union[IsotropicMaterial, OrthotropicMaterial]
        - An instance of `IsotropicMaterial` if the elastic properties are provided as single values or
        if the orthotropic material has equal values for E, G, and nu.
        - An instance of `OrthotropicMaterial` if the elastic properties are provided as lists and they are not identical.
    """
    material_name = material_data["name"]
    density = material_data["density"]
    elastic_props = material_data["elastic"]

    # Check if elastic properties are provided as lists for E, G, and nu (orthotropic material)
    if isinstance(elastic_props["E"], Iterable):
        # Check if all components are identical (i.e., isotropic material)
        if (
            all(val == elastic_props["E"][0] for val in elastic_props["E"])
            and all(val == elastic_props["nu"][0] for val in elastic_props["nu"])
            and all(val == elastic_props["G"][0] for val in elastic_props["G"])
        ):
            # It's an isotropic material, even if it was defined as orthotropic
            E = elastic_props["E"][0]
            nu = elastic_props["nu"][0]
            return Material(name=material_name, E=E, nu=nu, rho=density)
        else:
            # It's a true orthotropic material
            E = tuple(elastic_props["E"])
            G, nu = _windio_to_element_order(E, elastic_props["G"], elastic_props["nu"])
            return OrthotropicMaterial(name=material_name, E=E, G=G, nu=nu, rho=density)
    else:
        # It's an isotropic material
        E = elastic_props["E"]
        nu = elastic_props["nu"]
        return Material(name=material_name, E=E, nu=nu, rho=density)


def build_rust_properties(numad_data: dict) -> dict:
    """Build a Rust-native properties map from NuMAD mesh data.

    Converts section layup data directly to ``_aeroelast`` types
    (``Laminate`` / isotropic dict) — no Python ``Laminate``/``Ply``/
    ``CompositeShellProperty`` objects are created.

    Parameters
    ----------
    numad_data : dict
        NuMAD mesh data with ``"materials"`` and ``"sections"`` keys.

    Returns
    -------
    dict[str, _aeroelast.Laminate | dict]
        Suitable for passing directly to ``PyMeshAssembler.from_model()``.

    Notes
    -----
    Reference surface: the blade's skin nodes sit on the outer mould line (the wet
    surface FSI couples on) and each laminate is centred on them (``z`` from ``-h/2``
    to ``+h/2``), so half of every skin laminate lies outside the aerofoil.  On the
    IEA 15 MW blade that makes the shell's global GJ about 15 % and its static flap
    stiffness about 18 % higher than the deck.  The options (mid-surface nodes, or
    plies stacked inward from the OML, which reopens ADR-6) and their measured cost
    are recorded in ``odd/tasks/shell-reference-surface.md``.
    """
    if numad_data.get("structural_ready") is False:
        raise ValueError(
            "windIO 2.0 geometry-only mesh cannot build centered properties; physical FE offsets are not supported"
        )
    from _aeroelast import Laminate as RustLaminate, Ply as RustPly, OrthotropicMaterial as RustMat  # noqa: PLC0415

    # ── material lookup: name → raw data dict ─────────────────────────────
    mat_raw: Dict[str, dict] = {m["name"]: m for m in numad_data["materials"]}

    properties: dict = {}

    for section in numad_data["sections"]:
        set_name = section["elementSet"]
        layup = section["layup"]  # [[mat_name, thickness, angle], ...]

        rust_plies = []
        for mat_name, thickness, angle in layup:
            m = mat_raw.get(mat_name)
            if m is None:
                raise KeyError(
                    f"Material '{mat_name}' referenced in section '{set_name}' "
                    "not found in numad material database."
                )
            elastic = m["elastic"]
            rho = float(m["density"])

            E_raw = elastic["E"]
            nu_raw = elastic["nu"]
            G_raw = elastic.get("G")

            if isinstance(E_raw, Iterable) and not isinstance(E_raw, (str, float, int)):
                E = tuple(float(v) for v in E_raw)
                nu = tuple(float(v) for v in nu_raw)
                G = (
                    tuple(float(v) for v in G_raw)
                    if G_raw is not None
                    else (
                        E[0] / (2.0 * (1.0 + nu[0])),
                        E[1] / (2.0 * (1.0 + nu[1])),
                        E[2] / (2.0 * (1.0 + nu[2])),
                    )
                )
            else:
                E_v = float(E_raw)
                nu_v = float(nu_raw) if not isinstance(nu_raw, Iterable) else float(list(nu_raw)[0])
                G_v = E_v / (2.0 * (1.0 + nu_v))
                E, nu, G = (E_v, E_v, E_v), (nu_v, nu_v, nu_v), (G_v, G_v, G_v)

            G, nu = _windio_to_element_order(E, G, nu)
            rust_mat = RustMat(E[0], E[1], E[2], G[0], G[1], G[2], nu[0], nu[1], nu[2], rho)
            rust_plies.append(RustPly(rust_mat, float(thickness), float(angle)))

        # Single isotropic ply at 0° → isotropic dict (lighter path in Rust)
        if len(rust_plies) == 1 and rust_plies[0].angle == 0.0:
            m = mat_raw[layup[0][0]]
            elastic = m["elastic"]
            E_raw = elastic["E"]
            nu_raw = elastic["nu"]
            is_iso = not isinstance(E_raw, Iterable) or isinstance(E_raw, (str, float, int))
            if not is_iso:
                E_list = list(E_raw)
                nu_list = list(nu_raw)
                G_list = list(elastic.get("G", []))
                is_iso = (
                    all(v == E_list[0] for v in E_list)
                    and all(v == nu_list[0] for v in nu_list)
                    and (not G_list or all(v == G_list[0] for v in G_list))
                )
            if is_iso:
                E_v = (
                    float(E_list[0])
                    if isinstance(E_raw, Iterable) and not isinstance(E_raw, (str, float, int))
                    else float(E_raw)
                )
                nu_v = (
                    float(nu_list[0])
                    if isinstance(nu_raw, Iterable) and not isinstance(nu_raw, (str, float, int))
                    else float(nu_raw)
                )
                properties[set_name] = {
                    "type": "isotropic",
                    "e": E_v,
                    "nu": nu_v,
                    "rho": float(m["density"]),
                    "thickness": float(layup[0][1]),
                    "shear_correction": 5.0 / 6.0,
                }
                continue

        scf = 5.0 / 6.0
        properties[set_name] = RustLaminate(rust_plies, scf)

    return properties


@dataclass
class LayerDef:
    material: str
    thickens: float
    angle: float


@dataclass
class BladeSection:
    set_name: str
    x_dir: np.ndarray
    xy_dir: np.ndarray
    layers: List[LayerDef]


class BladeModel:
    def __init__(self, blade_yaml: str, element_size: Optional[float] = None) -> None:
        self.blade_mesh_obj = Blade(blade_yaml, element_size)
        print(self.get_element_material(1))

    def view(self):
        self.blade_mesh_obj.view()

    def write_mesh(self, filename, **kwargs):
        self.blade_mesh_obj.mesh.write_mesh(filename, **kwargs)

    @property
    def materials(self) -> Dict[str, Union[Material, OrthotropicMaterial]]:
        if self._materials_db:
            return self._materials_db
        else:
            for material in self._numad_mesh["materials"]:
                material = material_factory(material)
                self._materials_db[material.name] = material
        return self._materials_db

    @property
    def sections(self) -> Dict[str, BladeSection]:
        if self._sections_db:
            return self._sections_db
        else:
            for section in self._numad_mesh["sections"]:
                layers = [LayerDef(*l) for l in section["layup"]]
                section = BladeSection(
                    set_name=section["elementSet"],
                    x_dir=section["xDir"],
                    xy_dir=section["xyDir"],
                    layers=layers,
                )
                self._sections_db[section.set_name] = section
            return self._sections_db

    def get_element_material(self, element_id: int) -> BladeSection:
        element = self.mesh.get_element_by_id(element_id)
        element_sets = self.mesh.get_element_associated_set(element)

        # use the first asociate element set
        valid_set = element_sets[0]
        return self.sections[valid_set.name]

    def to_model_dict(self):
        model = {}
        for element in self.mesh.elements:
            section = self.get_element_material(element_id=element.id)
            model[element.id] = {
                "element_family": ElementFamily.SHELL,
                "layers": section.layers,
                "x_dir": section.x_dir,
                "xy_dir": section.xy_dir,
            }

            return model
