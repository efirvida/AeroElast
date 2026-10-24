# for type hints
import numpy as np
from numpy import ndarray

from aeroelast.models.blade.numad.io.excel_to_blade import excel_to_blade
from aeroelast.models.blade.numad.io.yaml_to_blade import yaml_to_blade
from aeroelast.models.blade.numad.objects.bom import BillOfMaterials
from aeroelast.models.blade.numad.objects.definition import Definition
from aeroelast.models.blade.numad.objects.geometry import Geometry
from aeroelast.models.blade.numad.objects.keypoints import KeyPoints
from aeroelast.models.blade.numad.objects.materialdb import MaterialDatabase
from aeroelast.models.blade.numad.objects.settings import BladeSettings
from aeroelast.models.blade.numad.objects.stackdb import StackDatabase
from aeroelast.models.blade.numad.utils.interpolation import interpolator_wrap


class Blade:
    """Blade class

    Parameters
    ----------
    filename : str, optional
        Directory and filename of blade input file to load into the
        Blade object.

    Attributes
    ----------
    name : str
        Name of the blade
    definition : Definition
        Object containing the definition of the blade.
    ispan : ndarray
        Array defining the interpolated stations from 0 to 1
    geometry : Geometry
        Object containing the interpolated geometry of the blade
    keypoints : KeyPoints
        Object containing information about keypoint locations and areas
    bill_of_materials : BillOfMaterials
    stackdb : StackDatabase
    materialdb : MaterialDatabase
    settings : BladeSettings

    Example
    -------

    blade = Blade("path/to/blade.yaml")
    """

    def __init__(self, filename: str | None = None):
        self.name: str | None = None
        self.definition: Definition = Definition()
        self.ispan: ndarray | None = None
        self.geometry: Geometry = Geometry()
        self.keypoints: KeyPoints = KeyPoints()
        self.bill_of_materials: BillOfMaterials = BillOfMaterials()
        self.stackdb: StackDatabase = StackDatabase()
        self.materialdb: MaterialDatabase = MaterialDatabase()
        self.settings: BladeSettings = BladeSettings()

        # read input file
        if filename:
            if "yaml" in filename or "yml" in filename:
                self.read_yaml(filename)
            elif "xls" in filename or "xlsx" in filename:
                self.read_excel(filename)
            else:
                raise Exception(
                    "Unknown filetype. Currently supported inputs are excel and yaml files."
                )
            self.name = filename.split(".")[0]
        else:
            self.name = "blade"
        return

    def __str__(self):
        attributes = ""
        for attr_name, attr_value in vars(self).items():
            if isinstance(attr_value, list):
                attributes += f"{attr_name}={len(attr_value)}, "
            elif isinstance(attr_value, np.ndarray):
                attributes += f"{attr_name}={attr_value.shape}, "
            else:
                attributes += f"{attr_name}={attr_value}, "
        return f"Blade with {attributes[:-2]}"

    def read_yaml(self, filename: str):
        """Populate blade attributes with yaml file data

        Parameters
        ----------
        filename: str
            name of yaml file to be read

        Returns
        -------
        self

        """
        yaml_to_blade(self, filename)
        return self

    def read_excel(self, filename: str, airfoil_dir: str | None = None):
        """Populate blade attributes with NuMAD Excel file data.

        Parameters
        ----------
        filename : str
            Path to ``.xlsx`` file in the standard NuMAD 3-sheet format
            (Geometry / Components / Materials).
        airfoil_dir : str, optional
            Directory containing airfoil coordinate ``.txt`` files.
            If *None*, looks for an ``airfoils/`` directory next to *filename*.

        Returns
        -------
        self
        """
        excel_to_blade(self, filename, airfoil_dir=airfoil_dir)
        return self

    def update_blade(self):
        """Generates geometry, keypoints, bill of materials,
        stack database, and material database based on the
        blade definition.
        """
        self.geometry.generate(self.definition)
        self.keypoints.generate(self.definition, self.geometry)
        self.bill_of_materials.generate(self.definition, self.keypoints)
        self.stackdb.generate(self.keypoints, self.bill_of_materials, self.definition)
        self.materialdb.generate(self.definition.materials, self.stackdb)
        return self

    def expand_blade_geometry_te(self, min_edge_lengths):
        """
        TODO: docstring
        """
        self.geometry.expand_blade_geometry_te(min_edge_lengths)
        self.keypoints.generate(self.definition, self.geometry)
        return

    def add_interpolated_station(self, span_location: float):
        """Adds an interpolated station to blade geometry

        Parameters
        ----------
        span_location : float
            location along span between 0 and 1.

        Returns
        -------
        int
            integer index where the new span was inserted
        """
        x0 = self.ispan
        if x0 is None or self.definition.ispan is None:
            raise RuntimeError("Blade.add_interpolated_station() requires an interpolated span")
        insertIndex = 0
        if span_location < x0[-1] and span_location > 0:
            for i_span, spanLocation in enumerate(x0[1:]):
                if span_location < spanLocation:
                    insertIndex = i_span + 1
                    break
        else:
            raise ValueError(f"A new span location with value {span_location} is not possible.")
        self.add_interpolated_stations([span_location])
        return insertIndex

    def _interpolate_band(self, x0, band):
        """Interpolate a LE/TE band width onto ``self.ispan``, keeping its leading zero run.

        ``KeyPoints.generate`` gives every station before the first nonzero one the width of that
        first station. Linear interpolation over the 0 -> width step (yaml: r/R 0.08 -> 0.10)
        would put a tiny nonzero width at the first inserted station, so the root stations
        would inherit it and lose the band (TE reinforcement 0.08 m wide instead of 0.67 m).
        """
        band = np.asarray(band, dtype=float)
        out = interpolator_wrap(x0, band, self.ispan)
        nonzero = np.nonzero(band)[0]
        if nonzero.size:
            out[self.ispan < x0[nonzero[0]]] = 0.0
        return out

    def add_interpolated_stations(self, span_locations):
        """Insert several interpolated stations and regenerate the blade once.

        Parameters
        ----------
        span_locations : sequence of float
            span positions strictly inside the current span, not already stations.
        """
        x0 = self.ispan
        if x0 is None or self.definition.ispan is None:
            raise RuntimeError("Blade.add_interpolated_stations() requires an interpolated span")
        new = np.asarray(span_locations, dtype=float)
        if new.size == 0:
            return
        if np.any(new <= 0) or np.any(new >= x0[-1]) or np.any(np.isin(new, x0)):
            raise ValueError(f"New span locations {new} are not possible.")
        self.ispan = np.sort(np.concatenate([x0, new]))
        self.definition.ispan = np.sort(np.concatenate([self.definition.ispan, new]))

        self.definition.leband = self._interpolate_band(x0, self.definition.leband)
        self.definition.teband = self._interpolate_band(x0, self.definition.teband)
        self.definition.sparcapwidth_hp = interpolator_wrap(
            x0, self.definition.sparcapwidth_hp, self.ispan
        )
        self.definition.sparcapwidth_lp = interpolator_wrap(
            x0, self.definition.sparcapwidth_lp, self.ispan
        )
        self.definition.sparcapoffset_hp = interpolator_wrap(
            x0, self.definition.sparcapoffset_hp, self.ispan
        )
        self.definition.sparcapoffset_lp = interpolator_wrap(
            x0, self.definition.sparcapoffset_lp, self.ispan
        )

        self.update_blade()
