import numpy as np
from numpy import ndarray
from scipy.interpolate import PchipInterpolator

from aeroelast.models.blade.numad.objects.airfoil import (
    get_airfoil_normals,
    get_airfoil_normals_angle_change,
)
from aeroelast.models.blade.numad.utils.affinetrans import rotation, translation
from aeroelast.models.blade.numad.utils.interpolation import interpolator_wrap


class Geometry:
    """Geometry class

    Contains the interpolated geometry of a blade

    Attributes
    ----------
    c : ndarray
        nondimensional (0-1) value of x axis
    camber : ndarray
    thickness : ndarray
    ic : ndarray
    icamber : ndarray
    ithickness : array
        interpolated thickness
    cpos : ndarray
        x-axis parametrization for airfoil curve
    ichord : ndarray
        interpolated chord
    ichordoffset : ndarray
        interpolated offset
    iaerocenter : ndarray
        interpolated aerocenter
    idegreestwist : ndarray
        interpolated twist
    ipercentthick : ndarray
        interpolated percent thickness
    profiles : array
        normalized airfoil profiles
    coordinates : ndarray
        actual x,y,z geometry
    xoffset : ndarray
        natural offset
    LEindex : ndarray
    iprebend : ndarray
        interpolated prebend
    isweep : ndarray
        interpolated sweep
    """

    def __init__(self):
        self.c: ndarray | None = None
        self.camber: ndarray | None = None
        self.thickness: ndarray | None = None
        self.ic: ndarray | None = None
        self.icamber: ndarray | None = None
        self.ithickness: ndarray | None = None
        self.cpos: ndarray | None = None
        self.ichord: ndarray | None = None
        self.ichordoffset: ndarray | None = None
        self.iaerocenter: ndarray | None = None
        self.idegreestwist: ndarray | None = None
        self.ipercentthick: ndarray | None = None
        self.profiles: ndarray | None = None
        self.coordinates: ndarray | None = None
        self.xoffset: ndarray | None = None
        self.LEindex: int | None = None
        self.iprebend: ndarray | None = None
        self.isweep: ndarray | None = None
        self._windio2_frames: ndarray | None = None
        self._windio2_reference: ndarray | None = None
        self.surface_fraction: ndarray | None = None

        # init properties
        self._natural_offset: int = 1
        self._rotorspin: int = 1
        self._swtwisted: int = 0

    def __eq__(self, other):
        attrs = [a for a in dir(self) if not a.startswith("__") and not callable(getattr(self, a))]
        for attr in attrs:
            self_attr = getattr(self, attr)
            other_attr = getattr(other, attr)
            if isinstance(self_attr, (int,)):
                if self_attr != other_attr:
                    return False
            elif isinstance(self_attr, ndarray):
                if (self_attr != other_attr).any():
                    return False
        return True

    def initialize_arrays(self, num_points: int, num_stations: int, num_istations: int):
        """_summary_

        Parameters
        ----------
        size : tuple[int]
            num_points, num_stations = size
        """
        size = (num_points, num_stations)
        isize = (num_points, num_istations)
        self.c = np.zeros(size)
        self.camber = np.zeros(size)
        self.thickness = np.zeros(size)
        self.ic = np.zeros(isize)
        self.icamber = np.zeros(isize)
        self.ithickness = np.zeros(isize)
        self.idegreestwist = np.zeros(isize)
        self.ichord = np.zeros(isize)
        self.ichordoffset = np.zeros(isize)
        self.iaerocenter = np.zeros(isize)
        self.idegreestwist = np.zeros(isize)
        self.ipercentthick = np.zeros(isize)
        self.iprebend = np.zeros(isize)
        self.isweep = np.zeros(isize)

    def generate(self, definition, *, windio2_eta=None):
        """Populates geometry attributes based on a given blade definition

        Parameters
        ----------
        definition : Definition

        Returns
        -------
        self
        """
        # Only the isolated v2 topology opts in; the v1/ordinary reader path
        # keeps its original operations and array dtypes.
        self._windio2_eta = windio2_eta
        if windio2_eta is not None:
            if definition.windio2_reference_axis is None:
                raise ValueError('Extended loft is windIO 2 only')
            from .windio2_precision import Pchip
            self._v2_pchip = Pchip
        else:
            self._v2_pchip = PchipInterpolator
        scalar = np.longdouble if windio2_eta is not None else float
        self.ispan = definition.ispan
        num_istations = self.ispan.size
        stations = definition.stations
        num_stations = len(stations)
        if num_stations > 0:
            num_points = len(stations[0].airfoil.c)
        else:
            raise Exception("BladeDef must have at least one station before updating geometry.")

        # first station must be at blade root to prevent extrapolation
        assert stations[0].spanlocation == 0, "first station must be at the blade root"

        self.initialize_arrays(num_points, num_stations, num_istations)
        if (
            self.c is None
            or self.camber is None
            or self.thickness is None
            or self.ic is None
            or self.icamber is None
            or self.ithickness is None
            or self.idegreestwist is None
            or self.ichord is None
            or self.ichordoffset is None
            or self.iaerocenter is None
            or self.ipercentthick is None
            or self.iprebend is None
            or self.isweep is None
        ):
            raise RuntimeError("Geometry.generate() requires initialised arrays")

        self._natural_offset = definition.natural_offset
        self._rotorspin = definition.rotorspin
        self._swtwisted = definition.swtwisted

        # Collect parameter tables from the stations.
        spanlocation = np.array([stations[i].spanlocation for i in range(len(stations))])
        tetype = ["_"] * num_stations
        for k in range(num_stations):
            station = stations[k]
            assert len(station.airfoil.c) == num_points, (
                "Station airfoils must have same number of samples."
            )

            self.c[:, k] = station.airfoil.c
            self.camber[:, k] = station.airfoil.camber
            self.thickness[:, k] = station.airfoil.thickness
            tetype[k] = station.airfoil.te_type

        # fix numerical issue due to precision on camber calculation
        # camber should start and end at y-values of zero
        self.camber[0, :] = np.zeros((1, num_stations))
        self.camber[-1, :] = np.zeros((1, num_stations))

        # Interpolate the station parameter tables.
        # Each column corresponds to an interpolated station.

        ## ic
        self.ic = interpolator_wrap(spanlocation, self.c, self.ispan, "pchip", axis=1)

        ## cpos
        self.cpos = np.concatenate(
            (
                -self.ic[-1, :].reshape(1, -1),
                -np.flipud(self.ic),
                self.ic[1:, :],
                self.ic[-1, :].reshape(1, -1),
            ),
            axis=0,
        )

        ## icamber
        self.icamber = interpolator_wrap(spanlocation, self.camber, self.ispan, "pchip", axis=1)

        ## ithickness
        self.ithickness = interpolator_wrap(
            spanlocation, self.thickness, self.ispan, "pchip", axis=1
        )
        # Adjust the thickness profiles based on te_type of stations.
        # This is mainly for transitions to flatbacks were the
        # interpolated airfoil needs to look like a round.
        for k in range(len(self.ispan)):
            try:
                ind = np.argwhere(self.ispan[k] < spanlocation)[0][0]
                # maybe better: ind = np.flatnonzero(self.ispan[k] < spanlocation)[0]
            except Exception:
                continue
            else:
                if ind == 1:
                    continue
            if tetype[ind] == "flat" and tetype[ind - 1] == "round":
                self.ithickness[-1, k] = 0

        # Interpolate the blade parameter curves.
        ## idegreestwist
        self.idegreestwist = interpolator_wrap(
            definition.span, definition.degreestwist, self.ispan, "pchip"
        )

        ## ichord
        if definition.chord_curve is not None:
            # the input's own chord grid, linear as the input defines it
            self.ichord = interpolator_wrap(*definition.chord_curve, self.ispan)
            self.ipercentthick = interpolator_wrap(
                definition.span, definition.percentthick, self.ispan, "pchip"
            )
        else:
            self.ichord = interpolator_wrap(definition.span, definition.chord, self.ispan, "pchip")

            ## ipercentthick
            absolutethick = np.multiply(definition.percentthick, definition.chord) / 100
            iabsolutethick = interpolator_wrap(definition.span, absolutethick, self.ispan, "pchip")
            self.ipercentthick = iabsolutethick / self.ichord * 100
        # ensure that the interpolation doesn't reduce the percent
        # thickness beneath the thinnest airfoil
        self.ipercentthick[self.ipercentthick < np.amin(definition.percentthick)] = np.amin(
            definition.percentthick
        )

        ## ichordoffset
        self.ichordoffset = interpolator_wrap(
            definition.span, definition.chordoffset, self.ispan, "pchip"
        )

        ## iaerocenter
        self.iaerocenter = interpolator_wrap(
            definition.span, definition.aerocenter, self.ispan, "pchip"
        )

        ## isweep
        if len(definition.sweep) == 0:
            definition.sweep = np.zeros(definition.span.shape)
        if len(definition.prebend) == 0:
            definition.prebend = np.zeros(definition.span.shape)

        self.isweep = interpolator_wrap(definition.span, definition.sweep, self.ispan, "pchip")

        ## iprebend
        self.iprebend = interpolator_wrap(definition.span, definition.prebend, self.ispan, "pchip")

        # Only v2 uses tangent frames. Keep the native reference curves: resampling
        # them first would change their PCHIP derivatives at inserted stations.
        self._windio2_frames = None
        self._windio2_reference = None
        if definition.windio2_reference_axis is not None:
            eta = self.ispan / definition.span[-1] if windio2_eta is None else windio2_eta
            outer = definition.windio2_outer_shape
            def native(field):
                curve = outer[field]
                return self._v2_pchip(curve["grid"], curve["values"])(eta)
            self.ichord = native("chord")
            self.ichordoffset = native("section_offset_y") / self.ichord
            self.idegreestwist = native("twist")
            self.ipercentthick = 100 * native("rthick")
            self._build_windio2_frames(definition)

        # Generate the blade surface self.
        n_istations = np.asarray(self.ispan).size
        n_areas = num_points * 2 + 1
        self.profiles = np.zeros((n_areas, 2, n_istations), dtype=scalar)
        self.coordinates = np.zeros((n_areas, 3, n_istations), dtype=scalar)
        self.xoffset = np.zeros((1, n_istations))
        self.LEindex = num_points
        if definition.windio2_reference_axis is not None:
            self._build_windio2_surface(definition, num_points)

        for k in range(n_istations):
            if definition.windio2_reference_axis is None:
                self.update_airfoil_profile(k)
            mtindex = np.argmax(self.ithickness[:, k])
            self.xoffset[0, k] = self.ic[mtindex, k]
            self.update_oml_geometry(k)

        # Calculate the arc length of each curve
        self.arclength = np.zeros((n_areas, n_istations), dtype=scalar)
        self.HParcx0 = np.zeros((1, n_istations))
        self.LParcx0 = np.zeros((1, n_istations))

        le_index = self.LEindex
        for k in range(n_istations):
            xx = self.coordinates[:, 0, k]
            yy = self.coordinates[:, 1, k]
            zz = self.coordinates[:, 2, k]
            arclen = np.sqrt(np.diff(xx) ** 2 + np.diff(yy) ** 2 + np.diff(zz) ** 2)
            arclen = np.concatenate((np.array([0]), np.cumsum(arclen)), axis=0)
            self.arclength[:, k] = arclen
            LEarcsum = self.arclength[le_index, k]
            self.arclength[:, k] = self.arclength[:, k] - LEarcsum

            # find where x=0 intersects the surface
            self.HParcx0[0, k] = (
                interpolator_wrap(xx[1 : le_index + 1], arclen[1 : le_index + 1], 0) - LEarcsum
            )
            self.LParcx0[0, k] = (
                interpolator_wrap(xx[-2 : le_index - 1 : -1], arclen[-2 : le_index - 1 : -1], 0)
                - LEarcsum
            )

        return self

    def update_airfoil_profile(self, k):
        """Updates kth airfoil profile

        Parameters
        ----------
        k : _type_
            _description_

        Returns
        -------
        _type_
            _description_
        """
        if (
            self.ithickness is None
            or self.ipercentthick is None
            or self.icamber is None
            or self.ic is None
            or self.profiles is None
        ):
            raise RuntimeError("Geometry.update_airfoil_profile() requires initialised arrays")
        thickness = self.ithickness[:, k]
        percentthick = self.ipercentthick[k]
        camber = self.icamber[:, k]
        c = self.ic[:, k]
        # jcb: note that I'm using max thickness about camber
        # instead of overall thickness of airfoil. We may need to
        # change this definition.
        maxthick = np.amax(thickness)
        tratio = percentthick / (maxthick * 100)
        thick = thickness * tratio
        hp = camber - 0.5 * thick
        lp = camber + 0.5 * thick
        profile1 = np.concatenate(([c[-1]], np.flipud(c), c[1:], [c[-1]]))
        profile2 = np.concatenate(([0], np.flipud(hp), lp[1:], [0]))
        profile = np.stack((profile1, profile2), axis=1)
        self.profiles[:, :, k] = profile
        return self

    def expand_blade_geometry_te(self, min_edge_lengths: ndarray):
        """Adjusts trailing edges for solid models.

        Opens trailing edge

        Parameters
        ----------
        min_edge_lengths : ndarray

        Returns
        -------
        Self
        """
        if (
            self.coordinates is None
            or self.ichord is None
            or self.profiles is None
            or self.ithickness is None
            or self.ipercentthick is None
            or self.ic is None
            or self.xoffset is None
        ):
            raise RuntimeError("Geometry.expand_blade_geometry_te() requires initialised arrays")
        nStations = self.coordinates.shape[2]

        for i_station in range(nStations):
            first_point = self.ichord[i_station] * self.profiles[-2, :, i_station]
            second_point = self.ichord[i_station] * self.profiles[1, :, i_station]
            edge_length = np.linalg.norm(second_point - first_point)

            # print(f"station {i_station}, edgeLength: {edge_length * 1000}\n")

            maxthick = np.amax(self.ithickness[:, i_station])
            mtindex = np.argmax(self.ithickness[:, i_station])
            tratio = self.ipercentthick[i_station] / (maxthick * 100)
            airfoil_thickness = self.ithickness[:, i_station] * tratio
            onset = self.ic[mtindex, i_station]
            if edge_length < min_edge_lengths[i_station]:
                # print(
                #     f"Updating station: {i_station} TE separation from {edge_length} to {min_edge_lengths[i_station]}"
                # )
                tet = (min_edge_lengths[i_station] - edge_length) / self.ichord[i_station]
                tes = 5 / 3 * tet  # slope of TE adjustment; 5/3*tet is "natural"
                # continuous first & second derivatives at 'onset'
                # maintain second & third derivative at mc==1 (TE)
                # adjust slope at mc==1 (TE) by tes
                A = np.array([[1, 1, 1, 1], [3, 4, 5, 6], [6, 12, 20, 30], [6, 24, 60, 120]])
                d = np.array([[tet], [tes], [0], [0]])
                p = np.linalg.solve(A, d)
                # onset = self(k).maxthick;  # start of TE modification, measured from LE
                vec = (self.ic[:, i_station] - onset) / (1 - onset)
                mc = np.maximum(vec, np.zeros(vec.shape))
                temod = np.vstack([mc**3, mc**4, mc**5, mc**6]).T @ p
                temod = temod.reshape(-1)
                airfoil_thickness = airfoil_thickness + temod

                self.ithickness[:, i_station] = airfoil_thickness / tratio
                self.update_airfoil_profile(i_station)

                mtindex = np.argmax(self.ithickness[:, i_station])
                self.xoffset[0, i_station] = self.ic[mtindex, i_station]
                self.update_oml_geometry(i_station)
                # first_point=self.ichord(i_station)*self.profiles(end-1,:,i_station);
                # second_point=self.ichord(i_station)*self.profiles(2,:,i_station);
                # edgeLength2=norm(second_point-first_point);
                # fprintf('station #i, edgeLength: #f, New edgeLength=#f, percent diff: #f\n',i_station,edgeLength*1000,edgeLength2*1000,(edgeLength2-edgeLength)/edgeLength2*100)
        return self

    def _build_windio2_surface(self, definition, count):
        """Whole-open-curve PCHIP, then pointwise master-rthick PCHIP.

        Locate the loft's minimum-x landmark only AFTER defining its shape.
        Sample that same curve around the landmark to retain NuMAD's fixed LE
        index; do not align master LEs or renormalize their separate sides.
        """
        # Use the same accurately located sampling landmark in both v2 query
        # views. Ordinary Geometry output still rounds to its float64 arrays;
        # the opt-in view retains extended profiles, coordinates and arclength.
        from .windio2_precision import Pchip, minimum_location
        interpolator = Pchip
        scalar = np.longdouble
        masters = sorted(definition.windio2_airfoils, key=lambda a: a["rthick"])
        thickness = np.array([a["rthick"] for a in masters], dtype=scalar)
        native_eta = getattr(self, '_windio2_eta', None)
        eta = self.ispan/definition.span[-1] if native_eta is None else native_eta
        field = definition.windio2_outer_shape['rthick']
        targets = interpolator(field['grid'], field['values'])(eta)
        if len(masters) < 2 or np.any(np.diff(thickness) <= 0):
            raise ValueError("windIO surface requires distinct master airfoil rthick values")
        if np.any(targets < thickness[0]) or np.any(targets > thickness[-1]):
            raise ValueError("windIO rthick lies outside the master airfoil range")
        curves = []
        for airfoil in masters:
            xy = np.array([airfoil["coordinates"][c] for c in "xy"], dtype=scalar).T
            # Source convention may run suction first; NuMAD needs pressure first.
            if np.sum(xy[:-1, 0]*xy[1:, 1] - xy[1:, 0]*xy[:-1, 1]) > 0:
                xy = xy[::-1]
            arc = np.r_[0., np.cumsum(np.linalg.norm(np.diff(xy, axis=0), axis=1))]
            if np.any(np.diff(arc) <= 0):
                raise ValueError("windIO airfoils require nonzero consecutive surface segments")
            curves.append(interpolator(arc / arc[-1], xy))
        self.surface_fraction = np.empty((2 * count - 1, len(self.ispan)), dtype=scalar)
        for k, relative_thickness in enumerate(targets):
            target = relative_thickness * 100
            def loft(fraction):
                return interpolator(thickness * 100,
                                    [curve(fraction) for curve in curves], axis=0)(target)
            # Official masters have a single leading edge. This changes only
            # sampling, not F(s,r) = PCHIP_r(PCHIP_s(source coordinates)).
            def derivative(fraction):
                values = [curve(fraction) for curve in curves]
                tangents = [curve.derivative()(fraction) for curve in curves]
                return interpolator(thickness * 100, values).directional(target, tangents)[..., 0]
            location = minimum_location(derivative)
            fraction = np.r_[np.linspace(scalar(0), location, count),
                             np.linspace(location, scalar(1), count)[1:]]
            self.surface_fraction[:, k] = fraction
            self.profiles[1:-1, :, k] = loft(fraction)
            self.profiles[0, :, k] = (self.profiles[1, :, k] + self.profiles[-2, :, k]) / 2
            self.profiles[-1, :, k] = self.profiles[0, :, k]
        # Signed chord positions remain the downstream NuMAD surface parameter.
        self.cpos = self.profiles[:, 0].copy()
        self.cpos[:self.LEindex] *= -1

    def _build_windio2_frames(self, definition):
        """windIO 63dcad6 construction guide steps 1-7, in windIO axes.

        Solve (My @ Mx) ez = t on the forward-axis branch: ax = -asin(ty),
        ay = atan2(tx, tz). This is not two independent projected slope angles.
        Native eta is normalized reference-axis length, not a z fraction.
        """
        native_eta = getattr(self, '_windio2_eta', None)
        scalar = np.longdouble if native_eta is not None else float
        eta = self.ispan / definition.span[-1] if native_eta is None else native_eta
        interpolator = getattr(self, '_v2_pchip', PchipInterpolator)
        curves = [interpolator(
            definition.windio2_reference_axis[c]["grid"],
            definition.windio2_reference_axis[c]["values"],
        ) for c in "xyz"]
        reference = np.column_stack([curve(eta) for curve in curves])
        tangent = np.column_stack([curve.derivative()(eta) for curve in curves])
        tangent /= np.linalg.norm(tangent, axis=1)[:, None]
        # windIO x points to suction, y to TE. NuMAD x points away from TE,
        # y to suction: a proper +90-degree basis change, not a reflection.
        to_numad = np.array([[0., -1., 0.], [1., 0., 0.], [0., 0., 1.]])
        self._windio2_reference = reference @ to_numad.T
        self._windio2_frames = np.empty((len(eta), 3, 3), dtype=scalar)
        for k, (tx, ty, tz) in enumerate(tangent):
            ax = -np.arcsin(np.clip(ty, -1., 1.))
            ay = np.arctan2(tx, tz)
            cx, sx, cy, sy = np.cos(ax), np.sin(ax), np.cos(ay), np.sin(ay)
            mx = np.array([[1., 0., 0.], [0., cx, -sx], [0., sx, cx]])
            my = np.array([[cy, 0., sy], [0., 1., 0.], [-sy, 0., cy]])
            mxy = my @ mx
            correction = np.arctan2(mxy[0, 1], mxy[0, 0])
            angle = np.deg2rad(self.idegreestwist[k]) + correction
            c, s = np.cos(angle), np.sin(angle)
            mz = np.array([[c, -s, 0.], [s, c, 0.], [0., 0., 1.]])
            self._windio2_frames[k] = to_numad @ mxy @ mz

    def update_oml_geometry(self, k):
        """ """
        if (
            self.profiles is None
            or self.xoffset is None
            or self.ichordoffset is None
            or self.ichord is None
            or self.idegreestwist is None
            or self.isweep is None
            or self.iprebend is None
            or self.ispan is None
            or self.coordinates is None
        ):
            raise RuntimeError("Geometry.update_oml_geometry() requires initialised arrays")
        if self._windio2_frames is not None:
            # Apply chord/offset in the local windIO section before rotating,
            # then translate by the full native reference (including root z).
            local = np.column_stack((
                self.profiles[:, 1, k] * self.ichord[k],
                (self.profiles[:, 0, k] - self.ichordoffset[k]) * self.ichord[k],
                np.zeros(self.profiles.shape[0]),
            ))
            self.coordinates[:, :, k] = (
                local @ self._windio2_frames[k].T + self._windio2_reference[k]
            )
            return self

        x = self.profiles[:, 0, k]
        y = self.profiles[:, 1, k]

        # self.xoffset[0,k] = c[mtindex]
        if self._natural_offset:
            x = x - self.xoffset[0, k]
        x = x - self.ichordoffset[k]  # apply chordwise offset
        x = x * self.ichord[k] * -1 * self._rotorspin  # scale by chord
        y = y * self.ichord[k]  # scale by chord
        twist = -1 * self._rotorspin * self.idegreestwist[k]

        # prepare for hgtransform rotate & translate
        coords = np.zeros((len(x), 4))
        coords[:, 0] = np.cos(np.deg2rad(twist)) * x - np.sin(np.deg2rad(twist)) * y
        coords[:, 1] = np.sin(np.deg2rad(twist)) * x + np.cos(np.deg2rad(twist)) * y
        coords[:, 2] = np.zeros(len(x))
        coords[:, 3] = np.ones(len(x))

        # use the generating line to translate and rotate the coordinates
        # NOTE currently, rotation is not assigned from blade properties
        # and defaults to 0
        prebend_rot = 0
        sweep_rot = 0
        """
        jcb: This code, copied from NuMAD 2.0, causes each section to rotate out
        of plane so that its normal follows the generating line direction. Need
        to replace 'twistFlag' with '-1*self.rotorspin' and calculate the slopes
        based on the available data. For now, default to parallel sections.
        if isequal(blade.PresweepRef.method,'normal')
            sweep_slope = ppval(blade.PresweepRef.dpp,sta.LocationZ);
            sweep_rot = atan(sweep_slope*twistFlag);
        end
        if isequal(blade.PrecurveRef.method,'normal')
            prebend_slope = ppval(blade.PrecurveRef.dpp,sta.LocationZ);
            prebend_rot = atan(-prebend_slope);
        endc
        """
        transX = -1 * self._rotorspin * self.isweep[k]
        transY = self.iprebend[k]
        transZ = self.ispan[k]
        Ry = rotation("y", sweep_rot)
        Rx = rotation("x", prebend_rot)
        R = Ry @ Rx
        T = translation(transX, transY, transZ)
        coords = coords @ np.transpose(R) @ np.transpose(T)
        # save the transformed coordinates
        self.coordinates[:, :, k] = coords[:, 0:3]
        return self

    def get_profile_te_type(self, k: int):
        """Get trailing edge type for the profile at the kth span location

        Parameters
        ----------
        k

        Return
        ------
        tetype : str
        """
        if self.profiles is None:
            raise RuntimeError("Geometry.get_profile_te_type() requires initialised profiles")
        xy = self.profiles[:, :, k]
        unitNormals = get_airfoil_normals(xy)
        angleChange = get_airfoil_normals_angle_change(unitNormals)
        disconts = np.flatnonzero(angleChange > 30)

        if np.std(angleChange) < 2:
            tetype = "round"
        elif len(disconts) > 1:
            tetype = "flat"
        else:
            tetype = "sharp"
        return tetype
