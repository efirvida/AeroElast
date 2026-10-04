import numpy as np
from numpy import ndarray


class Material:
    """Material class

    Attributes
    ----------
    name : str
        User selected name of the material
    type : str
        Two options: isotropic or orthotropic
    layerthickness : float
        Layer thickness [mm]
    ex : float
                Longitudinal elastic modulus [Pa]
    ey : float
                Transverse elastic modulus [Pa]
    ez : float
                Through-the-thickness elastic modulus in the
        principal material coordinates [Pa]
    gxy : float
                In-plane shear modulus [Pa]
    gyz : float
                Transverse shear modulus [Pa]
    gxz : float
                Transverse shear modulus [Pa]
    prxy : float
                In-plane Poisson ratio [ ]
    pryz : float
                Transverse Poisson ratio [ ]
    prxz : float
                Transverse Poisson ratio [ ]
    density : float
                Cured mass density [kg/m2]
    drydensity : float
                Density of fabric
    uts : float
                1 x 3 array of ultimate tensile strength design values.
        Sequence: SL , ST, Sz, 1 x 1 for isotropic.
    ucs : float
                1 x 3 array of ultimate compressive strength design values.
      Sequence: SL , ST, Sz, 1 x 1 for isotropic.
    uss : float
                1 x 3 array of ultimate shear strength design values.
        Sequence: SLT , STz, SLz, 1 x 1 for isotropic.
    xzit : float
                Lz tensile inclination parameter for Puck failure index
    xzic : float
                Lz compressive inclination parameter for Puck failure index
    yzit : float
                Tz tensile inclination parameter for Puck failure index
    yzic : float
                Tz compressive inclination parameter for Puck failure index
    g1g2 : float
                Fracture toughness ratio between GI (mode I) and GII (mode II) [ ]
    alp0 : float
                Fracture angle under pure transverse compression [degrees]
    etat : float
                Transverse friction coefficient for Larc [ ]
    etal : float
                Longitudinal friction coefficient for Larc [ ]
    m : list
                Fatigue slope exponent [ ]
    gamma_mf : list
                from DNL-GL standard, fatigue strength reduction factor
    gamma_ms : list
                from DNV-GL standard, short term strength reduction factor
    reference : str = None

    """

    def __init__(self):
        self.name: str | None = None
        self.type: str | None = None
        self.reference: str | list | None = None
        self.layerthickness: float = 0.0
        self.ex: float | np.ndarray = 0.0
        self.ey: float | np.ndarray = 0.0
        self.ez: float | np.ndarray = 0.0
        self.gxy: float | np.ndarray = 0.0
        self.gyz: float | np.ndarray = 0.0
        self.gxz: float | np.ndarray = 0.0
        self.prxy: float | np.ndarray = 0.0
        self.pryz: float | np.ndarray = 0.0
        self.prxz: float | np.ndarray = 0.0
        self.density: float = 0.0
        self.drydensity: float = 0.0
        self.uts: float | np.ndarray = 0.0
        self.ucs: float | np.ndarray = 0.0
        self.uss: float | np.ndarray = 0.0
        self.xzit: float = 0.0
        self.xzic: float = 0.0
        self.yzit: float = 0.0
        self.yzic: float = 0.0
        self.g1g2: float = 0.0
        self.alp0: float | np.ndarray | None = None
        self.etat: float | np.ndarray | None = None
        self.etal: float = 0.0
        self.m: list[float] = []
        self.gamma_mf: list[float] = []
        self.gamma_ms: list[float] = []

    def __eq__(self, other):
        attrs = vars(self).keys()

        for attr in attrs:
            self_attr = getattr(self, attr)
            other_attr = getattr(other, attr)
            try:
                if (self_attr != self_attr) & (other_attr != other_attr):
                    continue
            except ValueError:
                pass
            if isinstance(self_attr, (int, float, str, list)):
                if self_attr != other_attr:
                    return False
            elif isinstance(self_attr, ndarray):
                if (self_attr != other_attr).any():
                    return False
        return True
