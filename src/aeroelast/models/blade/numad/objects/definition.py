# for type hints
import numpy as np
from numpy import ndarray

from aeroelast.models.blade.numad.objects.airfoil import Airfoil
from aeroelast.models.blade.numad.objects.component import Component
from aeroelast.models.blade.numad.objects.material import Material
from aeroelast.models.blade.numad.objects.station import Station


class Definition:
    """Definition Class

    A Definition object is designed to contain the
    basic defining features of a blade, which can
    be populated manually or by reading a blade file.

    Attributes
    ----------
    components : dict
        Dictionary of components indexed by component name
    shearweb : list
        List of shearwebs. len(shearweb) == number of shearwebs.
        shearweb[i] gives an array defining the ith shearweb.
    materials : dict
        Dictionary of material objects indexed by material name.
    stacks : ndarray
    swstacks : array
    ispan : ndarray
        Spanwise locations of interpolated output
    aerocenter : ndarray
        Aerodynamic center of airfoil (used only by NuMAD->FAST)
    chord : ndarray
        Chord distribution [m]
    chordoffset : ndarray
        Chordwise offset (in addition to natural offset)
    degreestwist : ndarray
        Twist distribution [degrees]
    percentthick : ndarray
        Percent thickness of airfoil [%]
    prebend : ndarray
        Blade prebend, reference axis location along x2 [m]
    span : ndarray
        Spanwise location of distributed properties [m]
    sparcapoffset : ndarray
    sparcapwidth : ndarray
        Locations of keypoints b & c, defines distance
        between keypoints b & c [mm]. First entry is the HP spar cap.
        Second entry is the LP spar cap
    stations : list
        List of station objects
    sweep : ndarray
        Blade sweep, reference axis location along x1 [m]
    teband : ndarray
        Location of keypoint e
    leband : ndarray
        Location of keypoint a
    te_type : list
    """

    def __init__(self):
        self.components: dict[str, Component] = {}
        self.shearweb: list = []
        self.materials: dict[str, Material] = {}
        self.stations: list[Station] = []
        self.te_type: list[str] = []
        self.stacks: ndarray | None = None
        self.swstacks: ndarray | None = None
        self.ispan: ndarray | None = None
        self.aerocenter: ndarray | None = None
        self.chord: ndarray | None = None
        self.chordoffset: ndarray | None = None
        self.degreestwist: ndarray | None = None
        self.percentthick: ndarray | None = None
        self.prebend: ndarray | None = None
        self.span: ndarray | None = None
        self.sparcapoffset: ndarray | None = None
        self.sparcapwidth: ndarray | None = None
        self.sweep: ndarray | None = None
        self.teband: ndarray | None = None
        self.leband: ndarray | None = None
        self.sparcapwidth_hp: ndarray | None = None
        self.sparcapwidth_lp: ndarray | None = None
        self.sparcapoffset_hp: ndarray | None = None
        self.sparcapoffset_lp: ndarray | None = None
        self.rotor_diameter: float | None = None
        self.hub_diameter: float | None = None
        self.hub_height: float | None = None

        # init properties
        self._natural_offset: int = 1
        self._rotorspin: int = 1
        self._swtwisted: int = 0

    def __eq__(self, other):
        attrs = vars(self).keys()
        for attr in attrs:
            self_attr = getattr(self, attr)
            other_attr = getattr(other, attr)
            if isinstance(self_attr, (int, float, str, list, dict)):
                if self_attr != other_attr:
                    return False
            elif isinstance(self_attr, ndarray):
                if (self_attr != other_attr).any():
                    return False
        return True

    @property
    def natural_offset(self):
        """
        1 = offset by max thickness location,
        0 = do not offset to max thickness
        """
        return self._natural_offset

    @natural_offset.setter
    def natural_offset(self, new_natural_offset):
        if not (new_natural_offset == 0 or new_natural_offset == 1):
            raise Exception("natural_offset must be 0 or 1")
        else:
            self._natural_offset = new_natural_offset

    @property
    def rotorspin(self):
        """
        Rotor Spin,
        1 = CW rotation looking downwind,
        -1 = CCW rotation
        """
        return self._rotorspin

    @rotorspin.setter
    def rotorspin(self, new_rotorspin):
        if not (new_rotorspin == 1 or new_rotorspin == -1):
            raise Exception("rotorspin must be 1 (cw) or -1 (ccw)")
        else:
            self._rotorspin = new_rotorspin

    @property
    def swtwisted(self):
        """
        Shear Web,
        0 = planar shear webs,
        1 = shear webs twisted by blade twist
        """
        return self._swtwisted

    @swtwisted.setter
    def swtwisted(self, new_swtwisted):
        if not (new_swtwisted == 0 or new_swtwisted == 1):
            raise Exception("swtwisted must be 0 or 1")
        else:
            self._swtwisted = new_swtwisted

    def add_station(self, af: Airfoil, spanlocation: float):
        """This method adds a station

        Specifically, a station object is created
        and appended to self.stations.

        Parameters
        ----------
        af : airfoil
        spanlocation : float

        Returns
        -------
        None

        Example
        -------
        ``blade.add_station(af,spanlocation)`` where  ``af`` = airfoil filename
        or ``Airfoil`` object
        """
        new_station = Station(af)
        new_station.spanlocation = spanlocation
        self.stations.append(new_station)

        return self

    @property
    def blade_length(self) -> float | None:
        """Blade length inferred from the available span definition."""
        span_source = self.span if self.span is not None else self.ispan
        if span_source is None or len(span_source) == 0:
            return None

        span = np.asarray(span_source, dtype=float)
        return float(span.max() - span.min())

    @property
    def rotor_radius(self) -> float | None:
        """Rotor radius derived from ``rotor_diameter`` when available."""
        if self.rotor_diameter is None:
            return None
        return float(self.rotor_diameter) / 2.0

    def resolve_hub_radius(
        self, override: float | None = None, *, tolerance: float = 1e-6
    ) -> tuple[float, str]:
        """Resolve hub radius from explicit input or NuMAD/WindIO geometry.

        Resolution priority is:
        1. Explicit ``override`` provided by the caller.
        2. ``hub_diameter / 2`` when the hub definition exists.
        3. ``rotor_diameter / 2 - blade_length`` as a fallback.

        Returns
        -------
        tuple[float, str]
            The resolved hub radius and a short string describing the source.
        """
        if override is not None:
            radius = float(override)
            if radius < 0:
                raise ValueError(f"hub_radius must be non-negative, got {radius}")
            return radius, "explicit"

        direct_radius = None
        if self.hub_diameter is not None:
            direct_radius = float(self.hub_diameter) / 2.0

        derived_radius = None
        blade_length = self.blade_length
        rotor_radius = self.rotor_radius
        if rotor_radius is not None and blade_length is not None:
            derived_radius = rotor_radius - blade_length

        if direct_radius is not None:
            radius = direct_radius
            source = "hub_diameter"
        elif derived_radius is not None:
            radius = derived_radius
            source = "rotor_diameter_minus_blade_length"
        else:
            raise ValueError(
                "Unable to determine hub radius from blade definition. "
                "Provide hub_radius explicitly or ensure the YAML defines hub diameter "
                "or a consistent rotor diameter and blade span."
            )

        if radius < -tolerance:
            raise ValueError(
                f"Resolved hub radius is negative ({radius}). Check hub and rotor geometry."
            )
        if radius < 0:
            radius = 0.0

        return float(radius), source
