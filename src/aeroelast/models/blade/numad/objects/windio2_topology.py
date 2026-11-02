"""Native windIO boundary arrangement, independent of NuMAD's twelve regions.

Source laminate ownership is exposed without meshing. Arc curves are piecewise linear on
native eta knots (the reader's current explicit-arc convention); their domains
are not extended. Geometry queries reuse Geometry.generate, including the whole
curve loft and tangent frame, without modifying the loaded blade.
"""
from bisect import bisect_right
from copy import copy, deepcopy
from dataclasses import dataclass
from fractions import Fraction
from functools import cached_property, lru_cache
from itertools import combinations
from types import MappingProxyType

import numpy as np
from scipy.interpolate import PchipInterpolator

from .windio2_precision import Pchip, extended, extended_array, rational

from aeroelast.models.blade.numad.objects.geometry import Geometry
from aeroelast.models.blade.windio import normalize_windio
from aeroelast.models.blade.numad.objects.windio2_laminate import (
    LaminateCell, WindIO2Laminate, span_domain,
)


def _freeze(value):
    if isinstance(value, dict):
        return MappingProxyType({k: _freeze(v) for k, v in value.items()})
    if isinstance(value, list):
        return tuple(_freeze(v) for v in value)
    return value


# Source arcs that are the same physical event reach this reader through
# different binary64 evaluation paths (e.g. a linear interpolant vs a chained
# anchor) and differ by a few ulp of 1.0. A finite band that narrow is below the
# resolution of the double-precision source and of the extended loft, so its
# quadrilateral has no meaningful orientation. Distinct source events in the
# supported inputs are >= 1e-4 apart.
EVENT_IDENTIFICATION_TOL = Fraction(8, 2**53)


def _identify_ulp_events(values):
    """Give identities whose exact arcs differ by <= EVENT_IDENTIFICATION_TOL one value.

    A cluster is measured from its lowest member (its anchor), so a chain of
    near values cannot drift beyond the tolerance. Values within the tolerance
    of the interval ends 0 and 1 join them, because sections add those two
    knots on their own. Everything farther apart is left exactly as evaluated.
    """
    result, anchor = {}, None
    for identity, nd in sorted(values.items(), key=lambda item: item[1]):
        if nd <= EVENT_IDENTIFICATION_TOL:
            result[identity] = Fraction(0)
        elif nd >= 1-EVENT_IDENTIFICATION_TOL:
            result[identity] = Fraction(1)
        else:
            if anchor is None or nd-anchor > EVENT_IDENTIFICATION_TOL:
                anchor = nd
            result[identity] = anchor
    return result


@dataclass(frozen=True)
class BoundaryCurve:
    identity: tuple[str, str, str]
    surface: str  # "skin", or the source web name for web-local coordinates
    grid: tuple[float, ...]
    values: tuple[float, ...]

    def at(self, eta):
        """Return native linear arc position, or None outside the source domain."""
        if not self.grid[0] <= eta <= self.grid[-1]:
            return None
        return float(np.interp(eta, self.grid, self.values))

    @cached_property
    def _exact_data(self):
        # Exact predicates for the reader's parsed binary64 knot values, not
        # decimal reconstruction, coordinate snapping or a new interpolant.
        grid = tuple(map(Fraction, self.grid))
        values = tuple(map(Fraction, self.values))
        slopes = tuple((b-a)/(y-x) for x, y, a, b in
                       zip(grid[:-1], grid[1:], values[:-1], values[1:], strict=True))
        return grid, values, slopes

    def _exact_at(self, eta):
        grid, values, slopes = self._exact_data
        j = min(bisect_right(grid, eta)-1, len(slopes)-1)
        return values[j]+slopes[j]*(eta-grid[j])


@dataclass(frozen=True)
class BoundaryPoint:
    nd_arc: float  # source direction: midpoint -> suction -> pressure -> midpoint
    arc: float  # dimensional NuMAD arc, zero at the completed loft's LE landmark
    xyz: tuple[float, float, float]
    identities: tuple[tuple[str, str, str], ...]


@dataclass(frozen=True)
class TopologySection:
    eta: float
    points: tuple[BoundaryPoint, ...]  # sorted in SOURCE arc direction
    intervals: tuple[tuple[BoundaryPoint, BoundaryPoint], ...]
    webs: tuple[tuple[str, tuple[BoundaryPoint, BoundaryPoint]], ...]
    skin_cells: tuple[LaminateCell, ...]
    web_cells: tuple[LaminateCell, ...]


class WindIO2Topology:
    """Immutable source records and an isolated geometry snapshot for queries.

    ``curves`` retains all anchor/layer/web identities, including unreferenced
    anchors and web-local layer bounds. ``event_stations`` includes all native
    curve knots and pairwise same-surface linear intersections. Coincident runs
    contribute their endpoints, not infinitely many events. These are topology
    events, NOT a sufficient mesh-resolution prescription.
    """

    def __init__(self, source, source_sha256, definition):
        self.source = _freeze(source)
        self.source_sha256 = source_sha256
        self._definition = deepcopy(definition)
        structure = source["components"]["blade"]["structure"]
        # Reuse the existing authoritative explicit-anchor resolver. Synthetic
        # probe layers expose every anchor handle, even when no layer uses it.
        # Only the private normalization copy receives these probe records.
        probe = deepcopy(source)
        layers = probe["components"]["blade"]["structure"]["layers"]
        records = []
        for surface, anchors in [("skin", structure.get("anchors", [])), *[
            (w["name"], w.get("anchors", [])) for w in structure.get("webs", [])
        ]]:
            for anchor in anchors:
                for handle in ("start_nd_arc", "end_nd_arc", "midpoint_nd_arc"):
                    if handle in anchor:
                        records.append((surface, ("anchor", anchor["name"], handle)))
                        layers.append({"name": f"__topology_probe_{len(records)}",
                                       "start_nd_arc": deepcopy(anchor[handle])})
        resolved = normalize_windio(probe)["components"]["blade"]["structure"]
        curves = []

        def add(surface, identity, data):
            grid, values = np.asarray(data["grid"]), np.asarray(data["values"])
            if (grid.ndim != 1 or values.shape != grid.shape or len(grid) < 2
                    or not np.all(np.isfinite(grid)) or not np.all(np.isfinite(values))
                    or np.any(np.diff(grid) <= 0) or grid[0] < 0 or grid[-1] > 1
                    or np.any(values < 0) or np.any(values > 1)):
                raise ValueError(f"Invalid windIO topology curve: {identity}")
            curves.append(BoundaryCurve(identity, surface, tuple(grid), tuple(values)))

        original_count = len(structure["layers"])
        for (surface, identity), layer in zip(records, resolved["layers"][original_count:], strict=True):
            add(surface, identity, layer["start_nd_arc"])
        for kind, items in (("layer", resolved["layers"][:original_count]),
                            ("web", resolved.get("webs", []))):
            for item in items:
                for handle in ("start_nd_arc", "end_nd_arc"):
                    if handle in item:
                        add(item.get("web", "skin"), (kind, item["name"], handle), item[handle])
        for handle, value in (("start_nd_arc", 0.), ("end_nd_arc", 1.)):
            add("skin", ("perimeter", "TE", handle), {"grid": [0., 1.], "values": [value, value]})
        self.curves = tuple(curves)
        self._laminate = WindIO2Laminate(self.source)
        events = {Fraction(x) for x in (0., 1., *self._laminate.event_stations,
                                        *(x for curve in curves for x in curve.grid))}
        crossings = []

        @lru_cache(None)
        def exact_at(curve, eta):
            return curve._exact_at(Fraction(eta))
        for a, b in combinations(curves, 2):
            if a.surface != b.surface:
                continue
            lo, hi = max(a.grid[0], b.grid[0]), min(a.grid[-1], b.grid[-1])
            knots = sorted({lo, hi, *(x for x in (*a.grid, *b.grid) if lo < x < hi)})
            if hi <= lo:
                continue
            for left, right in zip(knots[:-1], knots[1:], strict=True):
                dl = exact_at(a, left)-exact_at(b, left)
                dr = exact_at(a, right)-exact_at(b, right)
                if dl * dr < 0:
                    event = Fraction(left)+(Fraction(right)-Fraction(left))*dl/(dl-dr)
                    events.add(event)
                    crossings.append((event, a.identity, b.identity, exact_at(a, event)))
        # One rounding per proven rational event makes emission independent of
        # pair ordering and of unrelated knots splitting coincident segments.
        # Distinct roots are NEVER unioned on the basis of their float images.
        self._exact_event_stations = tuple(sorted(events))
        images = {}
        for event in self._exact_event_stations:
            images.setdefault(float(event), []).append(event)
        self._event_collisions = tuple((image, tuple(roots)) for image, roots in images.items()
                                       if len(roots) > 1)
        self._crossings = tuple((float(event), a, b) for event, a, b, _ in crossings)
        self._crossing_positions = tuple((float(event), a, b, float(value))
                                         for event, a, b, value in crossings)
        self.event_stations = tuple(sorted(images))

    def with_extended_precision(self):
        """Isolated query view: exact chart identities and extended v2 loft.

        Original public binary64 queries and their fail-closed guard remain
        available; the geometry-only mesh uses this view, never float images.
        """
        result = copy(self)
        result._extended_precision = True
        result.event_stations = self._exact_event_stations
        return result

    def _positions(self, eta, surface):
        if getattr(self, '_extended_precision', False):
            eta = rational(eta)
            return _identify_ulp_events(
                {c.identity: c._exact_at(eta) for c in self.curves
                 if c.surface == surface and c.grid[0] <= eta <= c.grid[-1]})
        if any(eta == image for image, _ in self._event_collisions):
            raise ValueError('Unsupported float64 topology event precision: '
                             'distinct exact events have the same float64 image; '
                             'source events were not merged')
        values = {c.identity: c.at(eta) for c in self.curves if c.surface == surface
                  and c.at(eta) is not None}
        # Only the proven participating identities receive the exact common
        # crossing value. A rounded coordinate is not an identity predicate.
        for event, a, b, common in self._crossing_positions:
            if eta == event and a in values and b in values:
                values[a] = values[b] = common
        return values

    def skin_points(self, eta, nd_arcs):
        """Evaluate actual completed skin geometry at source arc fractions.

        No Cartesian interpolation across stations or across cell chords is
        used. The two TE limits evaluate to the exact same source midpoint.
        """
        precise = getattr(self, '_extended_precision', False)
        if not np.isfinite(extended(eta)) or not 0 <= eta <= 1:
            raise ValueError("eta must be finite and in [0, 1]")
        eta = rational(eta) if precise else float(eta)
        nd = extended_array(nd_arcs) if precise else np.asarray(nd_arcs, dtype=float)
        if not np.isfinite(nd).all() or np.any((nd < 0) | (nd > 1)):
            raise ValueError("nd_arcs must be finite and in [0, 1]")
        arc, position = self._skin_interpolator(eta)
        return position(arc[0] + (1-nd) * (arc[-1]-arc[0]))

    def _skin_interpolator(self, eta):
        definition = copy(self._definition)
        precise = getattr(self, '_extended_precision', False)
        native_eta = np.array([extended(eta)]) if precise else None
        definition.ispan = (native_eta * definition.span[-1] if precise else
                            np.array([eta * definition.span[-1]]))
        g = Geometry().generate(definition, windio2_eta=native_eta)
        arc = g.arclength[:, 0]
        distinct, indices = np.unique(arc, return_index=True)
        interpolator = Pchip if precise else PchipInterpolator
        return arc, interpolator(distinct, g.coordinates[indices, :, 0])

    def section(self, eta):
        """Query exact eta using the existing completed geometry construction.

        Coincident identities share one point; only strictly positive source-arc
        intervals are emitted. TE arc 0 and 1 remain distinct topological points
        even though their physical midpoint is the same. Web endpoints retain
        source start/end order and are not clipped into spar caps.
        """
        precise = getattr(self, '_extended_precision', False)
        if not np.isfinite(extended(eta)) or not 0 <= eta <= 1:
            raise ValueError("eta must be finite and in [0, 1]")
        eta = rational(eta) if precise else float(eta)
        # Repeated zero-length TE legs are removed only in the interpolator.
        arc, position = self._skin_interpolator(eta)
        values = self._positions(eta, "skin")
        groups = {}
        for identity, nd in values.items():
            groups.setdefault(nd, []).append(identity)
        points = []
        for nd, identities in sorted(groups.items()):
            distance = arc[0] + (1 - extended(nd)) * (arc[-1] - arc[0]) if precise else float(
                arc[0] + (1 - nd) * (arc[-1] - arc[0]))
            points.append(BoundaryPoint(nd, distance, tuple(position(distance)), tuple(identities)))
        by_id = {identity: p for p in points for identity in p.identities}
        webs, web_cells = [], []
        for web in self.source["components"]["blade"]["structure"].get("webs", ()):
            keys = [("web", web["name"], h) for h in ("start_nd_arc", "end_nd_arc")]
            lo, hi = span_domain(web)
            if all(key in by_id for key in keys) and lo <= eta <= hi:
                ends = tuple(by_id[key] for key in keys)
                webs.append((web["name"], ends))
                local = self._positions(eta, web["name"])
                knots = sorted({0., 1., *local.values()})
                bounds = tuple(zip(knots[:-1], knots[1:], strict=True))
                start, end = (np.asarray(p.xyz) for p in ends)
                xyz = tuple(tuple(tuple(start * (1-extended(q)) + end * extended(q))
                                  for q in pair) for pair in bounds)
                web_cells.extend(self._laminate.cells(eta, web["name"], local, bounds, xyz))
        intervals = tuple(zip(points[:-1], points[1:], strict=True))
        skin = self._laminate.cells(eta, "skin", values,
                                    tuple((a.nd_arc, b.nd_arc) for a, b in intervals),
                                    tuple((a.xyz, b.xyz) for a, b in intervals))
        return TopologySection(eta, tuple(points), intervals, tuple(webs), skin, tuple(web_cells))
