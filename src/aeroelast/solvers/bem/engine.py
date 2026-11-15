"""BEM aerodynamic solver wrapping CCBlade.

Provides a thin interface between the fem-shell aerodynamic data layer
(``BladeAero``) and the NREL CCBlade BEM solver.
"""

from dataclasses import dataclass

import numpy as np

from aeroelast.models.blade.aerodynamics import AirfoilAero, BladeAero


@dataclass
class BEMResult:
    """Results from a BEM computation.

    Attributes
    ----------
    r : ndarray
        Radial positions of blade stations (m).
    Np : ndarray
        Normal force per unit length at each station (N/m).
    Tp : ndarray
        Tangential force per unit length at each station (N/m).
    alpha : ndarray
        Angle of attack at each station (deg).
    cl : ndarray
        Lift coefficient at each station.
    cd : ndarray
        Drag coefficient at each station.
    a : ndarray
        Axial induction factor at each station.
    ap : ndarray
        Tangential induction factor at each station.
    thrust : float
        Integrated rotor thrust (N).
    torque : float
        Integrated rotor torque (N*m).
    power : float
        Integrated rotor power (W).
    cm : ndarray or None
        Moment coefficient at each station.
    Mp : ndarray or None
        Pitching moment per unit span at each station (N*m/m),
        about the aerodynamic centre.  Computed from
        ``Cm * q_c * chord`` where ``q_c`` is the dynamic
        pressure per unit span back-calculated from *Np*, *Tp*,
        *Cl*, and *Cd*.
    """

    r: np.ndarray
    Np: np.ndarray
    Tp: np.ndarray
    alpha: np.ndarray
    cl: np.ndarray
    cd: np.ndarray
    a: np.ndarray
    ap: np.ndarray
    thrust: float
    torque: float
    power: float
    # -- extended fields (industrial BEM standard) -------------------------
    cm: np.ndarray | None = None
    Mp: np.ndarray | None = None
    # Normal / tangential force coefficients
    cn: np.ndarray | None = None
    ct: np.ndarray | None = None
    # Relative inflow velocity (m/s) and chord Reynolds number
    W: np.ndarray | None = None
    Re: np.ndarray | None = None
    # Reference geometry per station (from the (possibly deformed) blade_aero)
    chord: np.ndarray | None = None
    twist_deg: np.ndarray | None = None
    # Integrated non-dimensional performance coefficients
    CP: float | None = None
    CT: float | None = None
    CQ: float | None = None
    # Blade root flap bending moment (N*m)
    Mb: float | None = None


POLAR_INTERP_MODES = ("linear", "spline", "smoothed")


def _make_table_airfoil_class():
    """Return the CCAirfoil subclass that evaluates the table without smoothing.

    CCBlade's own ``CCAirfoil`` fits ``RectBivariateSpline(..., s=0.1)`` (Cl), ``s=0.001`` (Cd)
    and ``s=0.0001`` (Cm): a *smoothing* spline that does not pass through the table nodes
    (Cl off by 0.007-0.013 outboard on the IEA 15 MW blade, odd/tasks/l4-bem-rigid-closure.md L7).
    """
    from ccblade.ccblade import CCAirfoil  # pyright: ignore[reportMissingImports]
    from scipy.interpolate import RectBivariateSpline

    class _TableAirfoil(CCAirfoil):
        """CCAirfoil whose lookup returns the table exactly at the nodes.

        ``mode="linear"``: piecewise-linear in alpha and in Reynolds number.
        ``mode="spline"``: interpolating (``s=0``) cubic spline in alpha, linear in Reynolds number.
        """

        def __init__(self, alpha_deg, Re, cl, cd, cm=(), x=(), y=(), AFName="DEFAULTAF", mode="linear"):
            super().__init__(alpha_deg, Re, cl, cd, cm=cm, x=x, y=y, AFName=AFName)
            self.interp_mode = mode
            # Re-fit through the nodes (CCAirfoil.__init__ has already padded a single Re).
            re_nodes = np.array([1e1, 1e15]) if self.one_Re else np.asarray(Re, dtype=float)
            cl = np.asarray(cl, dtype=float).reshape(len(alpha_deg), -1)
            cd = np.asarray(cd, dtype=float).reshape(len(alpha_deg), -1)
            cm_arr = np.asarray(cm, dtype=float).reshape(len(alpha_deg), -1) if self.use_cm else None
            if self.one_Re:
                cl, cd = np.c_[cl, cl], np.c_[cd, cd]
                cm_arr = np.c_[cm_arr, cm_arr] if cm_arr is not None else None
            self._re_nodes = re_nodes
            self._tables = {"cl": cl, "cd": cd}
            if cm_arr is not None:
                self._tables["cm"] = cm_arr
            if mode == "spline":
                kx = min(len(self.alpha) - 1, 3)
                ky = min(len(re_nodes) - 1, 1)
                self.cl_spline = RectBivariateSpline(self.alpha, re_nodes, cl, kx=kx, ky=ky, s=0.0)
                self.cd_spline = RectBivariateSpline(self.alpha, re_nodes, cd, kx=kx, ky=ky, s=0.0)
                if cm_arr is not None:
                    self.cm_spline = RectBivariateSpline(
                        self.alpha, re_nodes, cm_arr, kx=kx, ky=ky, s=0.0
                    )

        def _linear(self, name, alpha, Re):
            """Bilinear lookup: linear in alpha, linear in Re (clamped to the table)."""
            table = self._tables[name]
            re_nodes = self._re_nodes
            if np.ndim(alpha) == 0 and np.ndim(Re) == 0:
                re = min(max(float(Re), re_nodes[0]), re_nodes[-1])
                j1 = min(max(int(np.searchsorted(re_nodes, re)), 1), len(re_nodes) - 1)
                w = (re - re_nodes[j1 - 1]) / (re_nodes[j1] - re_nodes[j1 - 1])
                v0 = np.interp(float(alpha), self.alpha, table[:, j1 - 1])
                v1 = np.interp(float(alpha), self.alpha, table[:, j1])
                return (1.0 - w) * v0 + w * v1
            alpha = np.asarray(alpha, dtype=float)
            re = np.clip(np.asarray(Re, dtype=float), re_nodes[0], re_nodes[-1])
            j1 = np.clip(np.searchsorted(re_nodes, re), 1, len(re_nodes) - 1)
            w = (re - re_nodes[j1 - 1]) / (re_nodes[j1] - re_nodes[j1 - 1])
            a, j1, w = np.broadcast_arrays(alpha, j1, w)
            out = np.empty(a.shape)
            for k in np.ndindex(a.shape):
                v0 = np.interp(a[k], self.alpha, table[:, j1[k] - 1])
                v1 = np.interp(a[k], self.alpha, table[:, j1[k]])
                out[k] = (1.0 - w[k]) * v0 + w[k] * v1
            return out if out.ndim else float(out)

        def evaluate(self, alpha, Re, return_cm=False):
            if self.interp_mode == "linear":
                cl = self._linear("cl", alpha, Re)
                cd = self._linear("cd", alpha, Re)
                if self.use_cm and return_cm:
                    return cl, cd, self._linear("cm", alpha, Re)
                return cl, cd
            return super().evaluate(alpha, Re, return_cm=return_cm)

    return _TableAirfoil


def _build_ccairfoil(airfoil: AirfoilAero, interp: str = "linear"):
    """Convert an ``AirfoilAero`` to a CCBlade airfoil.

    ``interp`` selects the table lookup: ``"linear"`` / ``"spline"`` pass through the table
    nodes exactly (see :func:`_make_table_airfoil_class`); ``"smoothed"`` is CCBlade's stock
    smoothing spline.
    """
    from ccblade.ccblade import CCAirfoil  # pyright: ignore[reportMissingImports]

    if interp not in POLAR_INTERP_MODES:
        raise ValueError(f"polar_interp must be one of {POLAR_INTERP_MODES}, got {interp!r}")

    polars_sorted = sorted(airfoil.polars, key=lambda p: p.re)
    re_list = [p.re for p in polars_sorted]
    alpha_deg = np.rad2deg(polars_sorted[0].alpha)

    if len(polars_sorted) == 1:
        p = polars_sorted[0]
        cl = p.cl.reshape(-1, 1)
        cd = p.cd.reshape(-1, 1)
        cm = p.cm.reshape(-1, 1)
    else:
        cl = np.column_stack([p.cl for p in polars_sorted])
        cd = np.column_stack([p.cd for p in polars_sorted])
        cm = np.column_stack([p.cm for p in polars_sorted])

    coords = airfoil.coordinates
    x = coords[:, 0] if coords is not None and len(coords) > 0 else []
    y = coords[:, 1] if coords is not None and len(coords) > 0 else []

    if interp == "smoothed":
        return CCAirfoil(alpha_deg, re_list, cl, cd, cm=cm, x=x, y=y, AFName=airfoil.name)
    return _make_table_airfoil_class()(
        alpha_deg, re_list, cl, cd, cm=cm, x=x, y=y, AFName=airfoil.name, mode=interp
    )


class BEMSolver:
    """BEM solver wrapping CCBlade.

    Parameters
    ----------
    blade_aero : BladeAero
        Complete blade aerodynamic definition.
    rho : float
        Air density (kg/m^3).
    mu : float
        Dynamic viscosity (kg/(m*s)).
    precone : float
        Hub precone angle (deg).
    tilt : float
        Nacelle tilt angle (deg).
    yaw : float
        Nacelle yaw misalignment angle (deg).  Positive rotates the rotor
        clockwise when viewed from above (i.e. wind comes from the left).
    hub_height : float
        Hub height (m), used for wind shear calculation.
    shear_exp : float
        Wind shear exponent for power-law profile.
    use_prebend : bool
        Pass the blade's out-of-plane prebend (``blade_aero.prebend``, metres,
        negative = upwind: the windIO ``reference_axis.x`` and AeroDyn ``BlCrvAC``
        convention) to CCBlade as ``precurve``.  Default ``True``; it has no effect
        on a blade with zero prebend.  ``False`` runs the blade straight (sensitivity).
        CCBlade's ``precurve`` is positive upwind for negative precone
        (``x_az = -r sin(precone) + precurve cos(precone)``), hence
        ``precurve = -prebend``.  The AeroDyn cant-twist term
        ``pitch*(1/cos(cant) - 1)`` is NOT modelled (not in the AeroDyn theory
        manual; residual in odd/tasks/l4-bem-rigid-closure.md L8).
    polar_interp : {"spline", "linear", "smoothed"}
        Polar table lookup.  ``"spline"`` (default: interpolating cubic, ``s=0``; AeroDyn's
        ``InterpOrd`` default is 3, cubic spline, as the deck's polar files state) and
        ``"linear"`` return the table exactly at its nodes; ``"smoothed"`` is CCBlade's
        stock smoothing spline (``s=0.1`` on Cl), kept as an option.  On the IEA 15 MW deck
        ``spline`` and ``linear`` agree with each other to 0.001 in thrust/torque error.
    """

    def __init__(
        self,
        blade_aero: BladeAero,
        rho: float = 1.225,
        mu: float = 1.81206e-5,
        precone: float = 0.0,
        tilt: float = 0.0,
        yaw: float = 0.0,
        hub_height: float = 150.0,
        shear_exp: float = 0.2,
        use_prebend: bool = True,
        polar_interp: str = "spline",
    ):
        from ccblade.ccblade import CCBlade  # pyright: ignore[reportMissingImports]

        self.blade_aero = blade_aero

        # Build CCAirfoil objects (one per unique airfoil)
        af_cache: dict[str, object] = {}
        for airfoil in blade_aero.airfoils:
            if airfoil.name not in af_cache:
                af_cache[airfoil.name] = _build_ccairfoil(airfoil, polar_interp)

        # Per-station airfoil list
        af_list = [af_cache[s.airfoil.name] for s in blade_aero.stations]

        # Extract arrays — CCBlade expects theta in degrees
        r = blade_aero.r
        chord = blade_aero.chord
        theta_deg = np.rad2deg(blade_aero.twist)

        # Prebend: CCBlade precurve is positive upwind (see class docstring).
        self.use_prebend = bool(use_prebend)
        self.polar_interp = polar_interp
        prebend = np.asarray(blade_aero.prebend, dtype=float)
        self.precurve: np.ndarray | None = None
        precurve_kw: dict = {}
        if self.use_prebend and prebend.size == r.size and np.any(np.abs(prebend) > 0.0):
            self.precurve = -prebend
            precurve_kw = {"precurve": self.precurve, "precurveTip": float(self.precurve[-1])}

        self.rotor = CCBlade(
            r=r,
            chord=chord,
            theta=theta_deg,
            af=af_list,
            Rhub=blade_aero.hub_radius,
            Rtip=blade_aero.rotor_radius,
            B=blade_aero.n_blades,
            rho=rho,
            mu=mu,
            precone=precone,
            tilt=tilt,
            yaw=yaw,
            hubHt=hub_height,
            shearExp=shear_exp,
            **precurve_kw,
        )

    def polar(self, station: int, alpha, re: float = 1.0e7):
        """Evaluate the (cl, cd, cm) the BEM uses at ``station`` for ``alpha`` (rad)."""
        af = self.rotor.af[station]
        alpha = np.asarray(alpha, dtype=float)
        if getattr(af, "use_cm", False):
            cl, cd, cm = af.evaluate(alpha, re, return_cm=True)
        else:
            cl, cd = af.evaluate(alpha, re)
            cm = np.zeros_like(np.asarray(cl, dtype=float))
        return np.asarray(cl), np.asarray(cd), np.asarray(cm)

    def compute(
        self,
        v_inf: float,
        omega: float,
        pitch: float,
        azimuth: float = 0.0,
    ) -> BEMResult:
        """Run BEM analysis at the specified operating conditions.

        Parameters
        ----------
        v_inf : float
            Hub-height wind speed (m/s).
        omega : float
            Rotor rotational speed (RPM).
        pitch : float
            Blade collective pitch (deg, positive feather).
        azimuth : float
            Blade azimuthal position (deg).

        Returns
        -------
        BEMResult
            Distributed and integrated aerodynamic loads.
        """
        # BEM theory is singular at omega=0 (no tangential velocity).
        # Substitute a tiny rotation so the solver converges while keeping
        # the physics essentially unchanged (tip-speed ratio ~ 0).
        _OMEGA_EPS = 1e-3  # RPM
        omega_bem = omega if abs(omega) > _OMEGA_EPS else _OMEGA_EPS

        # Distributed loads at the given azimuth
        loads, _ = self.rotor.distributedAeroLoads(v_inf, omega_bem, pitch, azimuth)

        # Integrated loads (azimuth-averaged)
        outputs, _ = self.rotor.evaluate([v_inf], [omega_bem], [pitch])

        # -- Aerodynamic pitching moment per unit span ----------------------
        # Try to get Cm directly from CCBlade output; fall back to polar
        # evaluation at the computed angle of attack.
        Np = loads["Np"]
        Tp = loads["Tp"]
        cl = loads["Cl"]
        cd = loads["Cd"]
        cn = loads.get("Cn")
        ct = loads.get("Ct")
        W = loads.get("W")
        Re = loads.get("Re")

        # Integrated coefficients
        Mb = float(outputs["Mb"][0]) if "Mb" in outputs else None
        try:
            out_coef, _ = self.rotor.evaluate([v_inf], [omega_bem], [pitch], coefficients=True)
            CP = float(out_coef["CP"][0])
            CT = float(out_coef["CT"][0])
            CQ = float(out_coef["CQ"][0])
        except Exception:
            CP = CT = CQ = None

        try:
            cm = loads["Cm"]
        except (KeyError, TypeError):
            alpha_rad = np.deg2rad(loads["alpha"])
            cm = np.array(
                [
                    st.airfoil.get_polar(re=1e7).evaluate(
                        np.atleast_1d(alpha_rad[k]),
                    )[2][0]
                    for k, st in enumerate(self.blade_aero.stations)
                ]
            )

        # Back-calculate dynamic pressure per unit span:
        #   q_c = sqrt(Np² + Tp²) / sqrt(Cl² + Cd²)
        # Then: Mp = q_c · chord · Cm
        cl_cd_norm = np.sqrt(cl**2 + cd**2)
        safe = cl_cd_norm > 1e-12
        q_c = np.where(
            safe,
            np.sqrt(Np**2 + Tp**2) / np.where(safe, cl_cd_norm, 1.0),
            0.0,
        )
        Mp = q_c * self.blade_aero.chord * cm

        return BEMResult(
            r=self.blade_aero.r,
            Np=Np,
            Tp=Tp,
            alpha=loads["alpha"],
            cl=cl,
            cd=cd,
            a=loads["a"],
            ap=loads["ap"],
            thrust=float(outputs["T"][0]),
            torque=float(outputs["Q"][0]),
            power=float(outputs["P"][0]),
            cm=cm,
            Mp=Mp,
            cn=cn,
            ct=ct,
            W=W,
            Re=Re,
            chord=self.blade_aero.chord.copy(),
            twist_deg=np.rad2deg(self.blade_aero.twist),
            CP=CP,
            CT=CT,
            CQ=CQ,
            Mb=Mb,
        )
