import warnings

import numpy as np

##from pynumad.mesh_gen.shellClasses import shellRegion, elementSet, NuMesh3D, spatialGridList2D, spatialGridList3D
from aeroelast.models.blade.numad.mesh_gen.element_utils import correct_orient
from aeroelast.models.blade.numad.mesh_gen.mesh_tools import get_direction_cosines
from aeroelast.models.blade.numad.mesh_gen.surface import Surface
from aeroelast.models.blade.numad.utils.interpolation import interpolator_wrap


class SpanPlan:
    """Element rows and chordwise counts of a native mesh, to be reproduced row by row.

    A NuMAD (region, bay) holds one laminate and the mesher cuts each bay into
    ``ceil(span / size)`` rows, so a laminate cannot follow a taper finer than the bay. The plan
    records, for the native stations, the rows per bay and the chordwise element counts at the
    bay ends. ``bind`` then maps the bays of a blade carrying extra stations (one per native row
    boundary) onto it: every bay becomes one row, and ``counts`` returns the chordwise counts that
    row has in the native mesh (the transition triangles of a bay sit in its end row).
    """

    def __init__(self, z, n_span):
        self.z = np.asarray(z, dtype=float)
        self.n_span = np.asarray(n_span, dtype=int)
        self.chord_counts = {}
        self._bay = None

    @property
    def row_boundaries(self):
        """Interior z of every native row boundary (not a native station)."""
        out = []
        for b, n in enumerate(self.n_span):
            lo, hi = self.z[b], self.z[b + 1]
            out.extend(lo + (hi - lo) * np.arange(1, n) / n)
        return np.asarray(out, dtype=float)

    def bind(self, z_new):
        """Map each bay of ``z_new`` onto (native bay, row in bay); check the row counts."""
        z_new = np.asarray(z_new, dtype=float)
        mid = 0.5 * (z_new[:-1] + z_new[1:])
        native = np.clip(np.searchsorted(self.z, mid) - 1, 0, len(self.n_span) - 1)
        rows = np.bincount(native, minlength=len(self.n_span))
        if not np.array_equal(rows, self.n_span):
            raise ValueError("blade stations do not hold the native rows of the span plan")
        first = np.concatenate([[0], np.cumsum(rows)[:-1]])
        self._bay = [(int(b), int(k - first[b])) for k, b in enumerate(native)]

    def counts(self, key, i):
        """Chordwise (inboard, outboard) counts of the row of new bay ``i``, or None."""
        if self._bay is None:
            raise RuntimeError("SpanPlan.bind() has not run")
        b, k = self._bay[i]
        ends = self.chord_counts.get((key, b))
        if ends is None:
            return None
        n1, n3 = ends
        peak = max(n1, n3)
        rows = int(self.n_span[b])
        if n1 < n3 and k == 0:
            return n1, peak
        if n3 < n1 and k == rows - 1:
            return peak, n3
        return peak, peak


def _chordwise_counts(shellKp, elementSize, key, i, tip_rule, plan, record):
    """Chordwise element counts at the inboard and outboard edges of a region."""
    if plan is not None:
        planned = plan.counts(key, i)
        if planned is not None:
            return planned
    n1 = max(1, int(np.ceil(np.linalg.norm(shellKp[1, :] - shellKp[0, :]) / elementSize)))
    n3 = max(1, int(np.ceil(np.linalg.norm(shellKp[3, :] - shellKp[2, :]) / elementSize)))
    if tip_rule:
        n3 = n1
    if record is not None:
        record.chord_counts[(key, i)] = (n1, n3)
    return n1, n3


# Spar-cap segments, split where the yaml puts the shear webs: segment -> (side, key_arcs rows).
# Side 0 is HP (the pressure side), 1 is LP; the rows are the cap edges ``c`` and ``b``.
_CAP_SEGMENTS = {3: (0, 3, 4), 8: (1, 8, 9)}
_WEB_EDGE_MARGIN = 5.0e-3  # m: a web arc never closer than this to the cap edge it would collapse on
_WEB_MERGE = 4.0e-3  # m: two webs closer than this on one side share one attachment line
_QUAD3_ORDER = (
    (0, 0), (0, 3), (3, 3), (3, 0), (0, 1), (0, 2), (1, 3), (2, 3),
    (3, 2), (3, 1), (2, 0), (1, 0), (1, 1), (1, 2), (2, 2), (2, 1),
)  # (row, column) of the 16 "quad3" control points, in ShellRegion order


def _quad3_points(columns, row0):
    """The 16 control points of a "quad3" region from 4 chordwise columns, rows ``row0..row0+3``."""
    return np.array([columns[c][row0 + r] for r, c in _QUAD3_ORDER])


def _plan_cap_columns(key_arcs, web_arcs, lo_row, hi_row):
    """Arc of every cap boundary per station: cap edge, one line per web (sorted), cap edge.

    Returns ``(arcs (n_webs + 2, n_stations), position of each web group in it)``, or ``None`` when
    the webs cross inside the cap (the mesher then keeps them on the cap edges). A web outside the
    cap is held ``_WEB_EDGE_MARGIN`` inside it, and webs closer than ``_WEB_MERGE`` become one line.
    """
    lo, hi = key_arcs[lo_row], key_arcs[hi_row]
    margin = np.minimum(_WEB_EDGE_MARGIN, 0.1 * (hi - lo))
    arcs = np.clip(np.asarray(web_arcs, dtype=float), lo + margin, hi - margin)
    pos = np.argsort(np.argsort(arcs.mean(axis=1), kind="stable"), kind="stable")
    ordered = arcs[np.argsort(pos)]
    for r in range(len(ordered) - 1):
        gap = ordered[r + 1] - ordered[r]
        if np.any(gap < -_WEB_MERGE):
            return None
        close = gap < _WEB_MERGE
        mean = 0.5 * (ordered[r] + ordered[r + 1])
        ordered[r] = np.where(close, mean, ordered[r])
        ordered[r + 1] = np.where(close, mean, ordered[r + 1])
    return np.vstack((lo, ordered, hi)), pos


def _web_control_points(shellKp):
    """Fill the interior Bezier control points of a web quad (rows 4, 5, 8, 9, 12-15) in place.

    They sit at exact thirds of the straight edges; the original 0.6666/0.3333 weights summed to
    0.9999 and shifted every control point by 1e-4 |P| towards the origin.
    """
    for lo, hi, a, b in ((4, 5, 0, 1), (8, 9, 2, 3), (12, 13, 11, 6), (14, 15, 7, 10)):
        shellKp[lo, :] = (2.0 * shellKp[a, :] + shellKp[b, :]) / 3.0
        shellKp[hi, :] = (shellKp[a, :] + 2.0 * shellKp[b, :]) / 3.0


def get_shell_mesh(
    blade, elementSize, spanGrading="chord", include_webs=True, span_plan=None, plan_only=False
):
    """
    This method generates a finite element shell mesh for the blade, based on what is
    stored in blade.geometry.coordinates, blade.keypoints.key_points,
    and blade.geometry.profiles.  Output is given as a python dictionary.

    Parameters
    -----------
    blade: Blade
    forSolid: bool
    elementSize: float
    spanGrading: str
        Spanwise element size grading strategy.
        - "uniform": constant elementSize for all sections
        - "chord":   scale by sqrt(chord/chord_max) so tip gets finer elements
    include_webs: bool
        When *False* the shear webs are skipped and only the outer shell is
        meshed.  The result is the blade outer-mold-line surface, which is what
        a CFD wall boundary needs.  No ``*SW`` element sets, no
        ``allShearWebEls`` / ``allShearWebNods`` sets are produced.
    span_plan: SpanPlan, optional
        Plan of a native mesh whose stations are a subset of the blade's. Every bay of the
        blade then gets exactly one element row, and the chordwise element counts follow the
        plan, so the discretisation equals the native one while each row has its own laminate.
    plan_only: bool
        Build the ``SpanPlan`` of this blade (native row counts per bay) and return it
        without meshing.

    Returns
    -------
    meshData:
    Nodes and elements for outer shell and shear webs:
    nodes:
        - [x, y, z]
        - [x, y, z]
        ...
    elements:
        - [n1,n2,n3,n4]
        - [n1,n2,n3,n4]
        ...
    Set list and section list for the outer shell and shear webs.
    These are companion lists with the same length and order,
    so meshData['sets']['element'][i] corresponds to meshData['sections'][i]
    sets:
        element:
            - name: set1Name
              labels: [e1, e2, e3 ...]
            - name: set2Name
              labels: [e1, e2, e3 ...]
            ...
    sections:
        - type: 'shell'
          layup:
             - [materialid,thickness,angle] ## layer 1
             - [materialid,thickness,angle] ## layer 2
             ...
          elementSet: set1Name
        - type: 'shell'
          layup:
             - [materialid,thickness,angle] ## layer 1
             - [materialid,thickness,angle] ## layer 2
             ...
          elementSet: set2Name
    """
    geometry = blade.geometry
    coordinates = geometry.coordinates
    profiles = geometry.profiles
    key_points = blade.keypoints.key_points
    stacks = blade.stackdb.stacks
    swstacks = blade.stackdb.swstacks

    geomSz = coordinates.shape
    lenGeom = geomSz[0]
    numXsec = geomSz[2]
    XSCurvePts = np.array([], dtype=int)

    # The profile polyline as read, before the key points below replace its nearest vertices
    # (their arc lengths are those of the original vertices): web attachment lines are cut on it.
    profile0 = coordinates.copy()

    ## Determine the key curve points along the OML at each cross section
    for i in range(numXsec):
        keyPts = np.array([0])
        minDist = 1
        lePt = 0
        for j in range(lenGeom):
            prof = profiles[j, :, i]
            mag = np.linalg.norm(prof)
            if mag < minDist:
                minDist = mag
                lePt = j

        for j in range(5):
            kpCrd = key_points[j, :, i]
            minDist = geometry.ichord[i]
            pti = 1
            for k in range(lePt):
                ptCrd = coordinates[k, :, i]
                vec = ptCrd - kpCrd
                mag = np.linalg.norm(vec)
                if mag < minDist:
                    minDist = mag
                    pti = k
            keyPts = np.concatenate((keyPts, [pti]))
            coordinates[pti, :, i] = np.array(kpCrd)

        keyPts = np.concatenate((keyPts, [lePt]))
        for j in range(5, 10):
            kpCrd = key_points[j, :, i]
            minDist = geometry.ichord[i]
            pti = 1
            for k in range(lePt, lenGeom):
                ptCrd = coordinates[k, :, i]
                vec = ptCrd - kpCrd
                mag = np.linalg.norm(vec)
                if mag < minDist:
                    minDist = mag
                    pti = k
            keyPts = np.concatenate((keyPts, [pti]))
            coordinates[pti, :, i] = np.array(kpCrd)

        keyPts = np.concatenate((keyPts, [lenGeom - 1]))
        allPts = np.array([keyPts[0]])
        for j in range(0, len(keyPts) - 1):
            secPts = np.linspace(keyPts[j], keyPts[j + 1], 4)
            secPts = np.round(secPts).astype(int)
            allPts = np.concatenate((allPts, secPts[1:4]))

        XSCurvePts = np.vstack((XSCurvePts, allPts)) if XSCurvePts.size else allPts
    rws, cls = XSCurvePts.shape

    ## Create longitudinal splines down the blade through each of the key X-section points

    splineX = coordinates[XSCurvePts[0, :], 0, 0]
    splineY = coordinates[XSCurvePts[0, :], 1, 0]
    splineZ = coordinates[XSCurvePts[0, :], 2, 0]
    for i in range(1, rws):
        Xrow = coordinates[XSCurvePts[i, :], 0, i]
        splineX = np.vstack((splineX, Xrow.T))
        Yrow = coordinates[XSCurvePts[i, :], 1, i]
        splineY = np.vstack((splineY, Yrow.T))
        Zrow = coordinates[XSCurvePts[i, :], 2, i]
        splineZ = np.vstack((splineZ, Zrow.T))

    # Use physical Z positions as spline parameter so that the
    # intermediate interpolation points are uniformly spaced along the
    # blade span, not in station-index space.
    z_stations = splineZ[:, 0]
    z_span = z_stations[-1] - z_stations[0]
    if z_span > 0:
        spParam = (z_stations - z_stations[0]) / z_span
    else:
        spParam = np.linspace(0, 1, rws)

    # Build spParami with 3 aligned sub-intervals per section.
    # Every 3rd row corresponds to an original station.
    nSpi = rws + 2 * (rws - 1)
    spParami = np.empty(nSpi)
    for _i in range(rws - 1):
        t0, t1 = spParam[_i], spParam[_i + 1]
        spParami[3 * _i] = t0
        spParami[3 * _i + 1] = t0 + (t1 - t0) / 3
        spParami[3 * _i + 2] = t0 + 2 * (t1 - t0) / 3
    spParami[-1] = spParam[-1]
    splineXi = interpolator_wrap(spParam, splineX[:, 0], spParami, "pchip")
    splineYi = interpolator_wrap(spParam, splineY[:, 0], spParami, "pchip")
    splineZi = interpolator_wrap(spParam, splineZ[:, 0], spParami, "pchip")
    for i in range(1, cls):
        splineXi = np.vstack(
            [
                splineXi,
                interpolator_wrap(spParam, splineX[:, i], spParami, "pchip"),
            ]
        )
        splineYi = np.vstack(
            [
                splineYi,
                interpolator_wrap(spParam, splineY[:, i], spParami, "pchip"),
            ]
        )
        splineZi = np.vstack(
            [
                splineZi,
                interpolator_wrap(spParam, splineZ[:, i], spParami, "pchip"),
            ]
        )
    splineXi = splineXi.T
    splineYi = splineYi.T
    splineZi = splineZi.T

    ## Spar-cap sub-regions and web attachment lines.
    # Each web hangs on a line at the yaml's attachment arc (``keypoints.web_attach_arcs``), inside
    # the spar cap: the cap is split there into sub-regions A, B, C, ... that keep the cap laminate,
    # so the web edge shares its nodes with the skin. Without yaml arcs the webs stay on the cap
    # edges (key points ``b`` and ``c``), as NuMAD has them.
    def _col(index):
        return np.stack((splineXi[:, index], splineYi[:, index], splineZi[:, index]), axis=1)

    def _arc_column(arcs):
        pts = np.empty((numXsec, 3))
        for k in range(numXsec):
            arc_k = geometry.arclength[1:-1, k]
            for d in range(3):
                pts[k, d] = np.interp(arcs[k], arc_k, profile0[1:-1, d, k])
        return np.stack(
            [interpolator_wrap(spParam, pts[:, d], spParami, "pchip") for d in range(3)], axis=1
        )

    # The cap split does not depend on ``include_webs``: the no-webs surface mesh must keep
    # the same skin discretisation as the structural mesh, so it is planned either way.
    capSplit = {}  # segment -> (arcs, columns, intermediate columns, web group positions)
    web_attach = getattr(blade.keypoints, "web_attach_arcs", None)
    if web_attach is not None and len(web_attach) == len(swstacks):
        for seg, (side, lo_row, hi_row) in _CAP_SEGMENTS.items():
            plan = _plan_cap_columns(
                blade.keypoints.key_arcs, [a[side] for a in web_attach], lo_row, hi_row
            )
            if plan is None:
                warnings.warn(
                    "Shear-web attachment arcs cross inside the spar cap; "
                    "the webs stay on the cap edges.",
                    stacklevel=2,
                )
                capSplit = {}
                break
            arcs, pos = plan
            cols = [_col(3 * seg)] + [_arc_column(a) for a in arcs[1:-1]] + [_col(3 * seg + 3)]
            mids = [
                [_arc_column(arcs[k] + (arcs[k + 1] - arcs[k]) * t) for t in (1.0 / 3.0, 2.0 / 3.0)]
                for k in range(len(cols) - 1)
            ]
            capSplit[seg] = (arcs, cols, mids, pos)
    if capSplit:
        webHP = [capSplit[3][1][1 + capSplit[3][3][g]] for g in range(len(swstacks))]
        webLP = [capSplit[8][1][1 + capSplit[8][3][g]] for g in range(len(swstacks))]
    else:
        webHP = [_col(12), _col(9)]
        webLP = [_col(24), _col(27)]

    ## Generate the mesh using the splines as surface guides

    # ------------------------------------------------------------------
    # Pre-scan: consistent SPANWISE element counts with chord grading.
    #
    # Adjacent chordwise segments within a section share spanwise edges
    # and MUST have the same spanwise node count.  The effective spanwise
    # element size is scaled by the local chord so that sections with a
    # smaller chord (tip) get proportionally smaller elements, keeping
    # the aspect ratio more uniform across the blade.
    #
    #   elSz_span[i] = elementSize * sqrt(chord[i] / chord_max)
    #   nSpan[i]     = max over 12 segments of ceil(span_edge / elSz_span)
    #
    # sqrt is used for moderate grading — linear would be too aggressive.
    # Chordwise counts remain per-edge (triangle transitions at chord
    # changes are the correct mechanism).
    # ------------------------------------------------------------------
    nSections = rws - 1
    nSegments = 12

    # Compute chord (chordwise extent) at each section boundary.
    # Chord is measured as the sum of the 12 chordwise edges at the
    # inboard face of each section (row _stPt in spline grid).
    _chord = np.zeros(nSections)
    _span_trail = np.zeros((nSections, nSegments))
    _span_lead = np.zeros((nSections, nSegments))
    _stPt = 0
    for _i in range(nSections):
        _stSp = 0
        _c = 0.0
        for _j in range(nSegments):
            p0 = np.array([splineXi[_stPt, _stSp], splineYi[_stPt, _stSp], splineZi[_stPt, _stSp]])
            p1 = np.array(
                [splineXi[_stPt, _stSp + 3], splineYi[_stPt, _stSp + 3], splineZi[_stPt, _stSp + 3]]
            )
            p2 = np.array(
                [
                    splineXi[_stPt + 3, _stSp + 3],
                    splineYi[_stPt + 3, _stSp + 3],
                    splineZi[_stPt + 3, _stSp + 3],
                ]
            )
            p3 = np.array(
                [splineXi[_stPt + 3, _stSp], splineYi[_stPt + 3, _stSp], splineZi[_stPt + 3, _stSp]]
            )
            _c += np.linalg.norm(p1 - p0)
            _span_trail[_i, _j] = np.linalg.norm(p2 - p1)
            _span_lead[_i, _j] = np.linalg.norm(p0 - p3)
            _stSp += 3
        _chord[_i] = _c
        _stPt += 3

    # Chord-graded or uniform spanwise element size
    if spanGrading == "chord":
        _chord_max = _chord.max()
        _elSz_span = elementSize * np.sqrt(np.clip(_chord / _chord_max, 0.1, 1.0))
    else:
        _elSz_span = np.full(nSections, elementSize)

    _nEl_span = np.maximum(
        np.ceil(_span_trail / _elSz_span[:, None]).astype(int),
        np.ceil(_span_lead / _elSz_span[:, None]).astype(int),
    )
    nSpan = np.maximum(_nEl_span.max(axis=1), 1)
    record = None
    if plan_only:
        record = SpanPlan(splineZi[0::3, 0], nSpan)
    elif span_plan is not None:
        span_plan.bind(splineZi[0::3, 0])
        nSpan = np.ones(nSections, dtype=int)

    bladeSurf = Surface()
    ## Outer shell sections
    outShES = set()
    secList = []
    stPt = 0
    for i in range(rws - 1):
        # if stPt < frstXS:
        #     stSec = 0
        #     endSec = 11
        #     stSp = 0
        # else:
        #     stSec = 1
        #     endSec = 10
        #     stSp = 3
        stSec = 0
        endSec = 11
        stSp = 0
        for j in range(stSec, endSec + 1):
            shellKp = np.array(
                [
                    [splineXi[stPt, stSp], splineYi[stPt, stSp], splineZi[stPt, stSp]],
                    [
                        splineXi[stPt, stSp + 3],
                        splineYi[stPt, stSp + 3],
                        splineZi[stPt, stSp + 3],
                    ],
                    [
                        splineXi[stPt + 3, stSp + 3],
                        splineYi[stPt + 3, stSp + 3],
                        splineZi[stPt + 3, stSp + 3],
                    ],
                    [
                        splineXi[stPt + 3, stSp],
                        splineYi[stPt + 3, stSp],
                        splineZi[stPt + 3, stSp],
                    ],
                    [
                        splineXi[stPt, stSp + 1],
                        splineYi[stPt, stSp + 1],
                        splineZi[stPt, stSp + 1],
                    ],
                    [
                        splineXi[stPt, stSp + 2],
                        splineYi[stPt, stSp + 2],
                        splineZi[stPt, stSp + 2],
                    ],
                    [
                        splineXi[stPt + 1, stSp + 3],
                        splineYi[stPt + 1, stSp + 3],
                        splineZi[stPt + 1, stSp + 3],
                    ],
                    [
                        splineXi[stPt + 2, stSp + 3],
                        splineYi[stPt + 2, stSp + 3],
                        splineZi[stPt + 2, stSp + 3],
                    ],
                    [
                        splineXi[stPt + 3, stSp + 2],
                        splineYi[stPt + 3, stSp + 2],
                        splineZi[stPt + 3, stSp + 2],
                    ],
                    [
                        splineXi[stPt + 3, stSp + 1],
                        splineYi[stPt + 3, stSp + 1],
                        splineZi[stPt + 3, stSp + 1],
                    ],
                    [
                        splineXi[stPt + 2, stSp],
                        splineYi[stPt + 2, stSp],
                        splineZi[stPt + 2, stSp],
                    ],
                    [
                        splineXi[stPt + 1, stSp],
                        splineYi[stPt + 1, stSp],
                        splineZi[stPt + 1, stSp],
                    ],
                    [
                        splineXi[stPt + 1, stSp + 1],
                        splineYi[stPt + 1, stSp + 1],
                        splineZi[stPt + 1, stSp + 1],
                    ],
                    [
                        splineXi[stPt + 1, stSp + 2],
                        splineYi[stPt + 1, stSp + 2],
                        splineZi[stPt + 1, stSp + 2],
                    ],
                    [
                        splineXi[stPt + 2, stSp + 2],
                        splineYi[stPt + 2, stSp + 2],
                        splineZi[stPt + 2, stSp + 2],
                    ],
                    [
                        splineXi[stPt + 2, stSp + 1],
                        splineYi[stPt + 2, stSp + 1],
                        splineZi[stPt + 2, stSp + 1],
                    ],
                ]
            )
            # The spar cap is cut where the webs attach; every piece keeps the cap's stack.
            # Piece k keeps its letter (and its span-plan key) when a lower piece is absent.
            pieces = [(0, stacks[j, i].name, shellKp)]
            if j in capSplit:
                arcs, cols, mids, _ = capSplit[j]
                pieces = [
                    (
                        k,
                        f"{stacks[j, i].name}_{chr(ord('A') + k)}",
                        _quad3_points([cols[k], *mids[k], cols[k + 1]], stPt),
                    )
                    for k in range(len(cols) - 1)
                    if np.max(arcs[k + 1, i : i + 2] - arcs[k, i : i + 2]) > 1.0e-6
                ]
            for k, pieceName, shellKp in pieces:
                # Per-edge chordwise + consistent nSpan spanwise
                nEl1, nEl3 = _chordwise_counts(
                    shellKp, elementSize, ("shell", j, k), i, i == rws - 2, span_plan, record
                )
                nEl = np.array([nEl1, nSpan[i], nEl3, nSpan[i]])

                bladeSurf.addShellRegion(
                    "quad3",
                    shellKp,
                    nEl,
                    name=pieceName,
                    elType="quad",
                    meshMethod="structured",
                )
                outShES.add(pieceName)
                newSec = {}
                newSec["type"] = "shell"
                layup = []
                for pg in stacks[j, i].plygroups:
                    totThick = 0.001 * pg.thickness * pg.nPlies
                    ply = [pg.materialid, totThick, pg.angle]
                    layup.append(ply)
                newSec["layup"] = layup
                newSec["elementSet"] = pieceName
                newSec["xDir"] = (shellKp[3, :] - shellKp[0, :]) + (shellKp[2, :] - shellKp[1, :])
                newSec["xyDir"] = (shellKp[1, :] - shellKp[0, :]) + (shellKp[2, :] - shellKp[3, :])
                secList.append(newSec)
            stSp = stSp + 3
        stPt = stPt + 3

    ## Shear web sections
    # ``include_webs=False`` collapses the loop range so the whole block is
    # skipped: ``swES``/``secList`` stay outer-shell only and the web node sets
    # are never created downstream.
    swES = set()
    stPt = 0
    np.array([])
    np.array([])
    for i in range(rws - 1 if include_webs else 0):
        if swstacks[0][i].plygroups:
            shellKp = np.zeros((16, 3))
            shellKp[0, :] = webHP[0][stPt]
            shellKp[1, :] = webLP[0][stPt]
            shellKp[2, :] = webLP[0][stPt + 3]
            shellKp[3, :] = webHP[0][stPt + 3]
            shellKp[6, :] = webLP[0][stPt + 1]
            shellKp[7, :] = webLP[0][stPt + 2]
            shellKp[10, :] = webHP[0][stPt + 2]
            shellKp[11, :] = webHP[0][stPt + 1]
            _web_control_points(shellKp)

            nEl1, nEl3 = _chordwise_counts(
                shellKp, elementSize, ("web", 0), i, False, span_plan, record
            )
            nEl = np.array([nEl1, nSpan[i], nEl3, nSpan[i]])

            bladeSurf.addShellRegion(
                "quad3",
                shellKp,
                nEl,
                name=swstacks[0][i].name,
                elType="quad",
                meshMethod="structured",
            )
            swES.add(swstacks[0][i].name)
            newSec = {}
            newSec["type"] = "shell"
            layup = []
            for pg in swstacks[0][i].plygroups:
                totThick = 0.001 * pg.thickness * pg.nPlies
                ply = [pg.materialid, totThick, pg.angle]
                layup.append(ply)
            newSec["layup"] = layup
            newSec["elementSet"] = swstacks[0][i].name
            newSec["xDir"] = np.array([0.0, 0.0, 1.0])
            newSec["xyDir"] = (shellKp[1, :] - shellKp[0, :]) + (shellKp[2, :] - shellKp[3, :])
            secList.append(newSec)
        if swstacks[1][i].plygroups:
            shellKp = np.zeros((16, 3))
            shellKp[0, :] = webLP[1][stPt]
            shellKp[1, :] = webHP[1][stPt]
            shellKp[2, :] = webHP[1][stPt + 3]
            shellKp[3, :] = webLP[1][stPt + 3]
            shellKp[6, :] = webHP[1][stPt + 1]
            shellKp[7, :] = webHP[1][stPt + 2]
            shellKp[10, :] = webLP[1][stPt + 2]
            shellKp[11, :] = webLP[1][stPt + 1]
            _web_control_points(shellKp)

            nEl1, nEl3 = _chordwise_counts(
                shellKp, elementSize, ("web", 1), i, False, span_plan, record
            )
            nEl = np.array([nEl1, nSpan[i], nEl3, nSpan[i]])

            bladeSurf.addShellRegion(
                "quad3",
                shellKp,
                nEl,
                name=swstacks[1][i].name,
                elType="quad",
                meshMethod="structured",
            )
            swES.add(swstacks[1][i].name)
            newSec = {}
            newSec["type"] = "shell"
            layup = []
            for pg in swstacks[1][i].plygroups:
                totThick = 0.001 * pg.thickness * pg.nPlies
                ply = [pg.materialid, totThick, pg.angle]
                layup.append(ply)
            newSec["layup"] = layup
            newSec["elementSet"] = swstacks[1][i].name
            newSec["xDir"] = np.array([0.0, 0.0, 1.0])
            newSec["xyDir"] = (shellKp[1, :] - shellKp[0, :]) + (shellKp[2, :] - shellKp[3, :])
            secList.append(newSec)
        stPt = stPt + 3

    if plan_only:
        return record

    ## Generate Shell mesh

    # print("getting blade mesh")
    shellData = bladeSurf.getSurfaceMesh()
    shellData["sections"] = secList

    ## Get local direction cosine orientations for individual elements
    # print("getting element orientations")
    nodes = shellData["nodes"]
    elements = shellData["elements"]
    numEls = len(shellData["elements"])
    elOri = np.zeros((numEls, 9), dtype=float)

    esi = 0
    for sec in secList:
        dirCos = get_direction_cosines(sec["xDir"], sec["xyDir"])
        es = shellData["sets"]["element"][esi]
        for ei in es["labels"]:
            eX = []
            eY = []
            eZ = []
            for ndi in elements[ei]:
                if ndi > -1:
                    eX.append(nodes[ndi, 0])
                    eY.append(nodes[ndi, 1])
                    eZ.append(nodes[ndi, 2])
            if len(eX) == 3:
                elType = "shell3"
            else:
                elType = "shell4"
            elCrd = np.array([eX, eY, eZ])
            elDirCos = correct_orient(dirCos, elCrd, elType)
            elOri[ei, 0:3] = elDirCos[0]
            elOri[ei, 3:6] = elDirCos[1]
            elOri[ei, 6:9] = elDirCos[2]
        esi = esi + 1

    shellData["elementOrientations"] = elOri

    ## Get all outer shell and all shear web element sets

    outerLab = []
    swLab = []
    for es in shellData["sets"]["element"]:
        nm = es["name"]
        if nm in outShES:
            outerLab.extend(es["labels"])
        elif nm in swES:
            swLab.extend(es["labels"])
    outerSet = {}
    outerSet["name"] = "allOuterShellEls"
    outerSet["labels"] = outerLab
    shellData["sets"]["element"].append(outerSet)
    swSet = {}
    swSet["name"] = "allShearWebEls"
    swSet["labels"] = swLab
    shellData["sets"]["element"].append(swSet)

    ## Get root (Zmin) node set
    minZ = np.min(splineZi)
    rootLabs = []
    lab = 0
    for nd in nodes:
        if np.abs(nd[2] - minZ) < 0.25 * elementSize:
            rootLabs.append(lab)
        lab = lab + 1
    newSet = {}
    newSet["name"] = "RootNodes"
    newSet["labels"] = rootLabs
    try:
        shellData["sets"]["node"].append(newSet)
    except Exception:
        nodeSets = []
        nodeSets.append(newSet)
        shellData["sets"]["node"] = nodeSets

    matList = []
    for mn in blade.definition.materials:
        newMat = {}
        mat = blade.definition.materials[mn]
        newMat["name"] = mat.name
        newMat["density"] = mat.density
        newMat["elastic"] = {}
        newMat["elastic"]["E"] = [mat.ex, mat.ey, mat.ez]
        if mat.type == "isotropic":
            nu = mat.prxy
            newMat["elastic"]["nu"] = [nu, nu, nu]
        else:
            newMat["elastic"]["nu"] = [mat.prxy, mat.prxz, mat.pryz]
        newMat["elastic"]["G"] = [mat.gxy, mat.gxz, mat.gyz]
        matList.append(newMat)

    shellData["materials"] = matList

    shellData["splineXi"] = splineXi
    shellData["splineYi"] = splineYi
    shellData["splineZi"] = splineZi

    return shellData
