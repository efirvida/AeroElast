import warnings

import numpy as np

##from pynumad.mesh_gen.shellClasses import shellRegion, elementSet, NuMesh3D, spatialGridList2D, spatialGridList3D
from aeroelast.models.blade.numad.mesh_gen.element_utils import correct_orient
from aeroelast.models.blade.numad.mesh_gen.local_refine import MeshedRegion, refine_surface_patch
from aeroelast.models.blade.numad.mesh_gen.mesh_tools import get_direction_cosines
from aeroelast.models.blade.numad.mesh_gen.shell_region import ShellRegion
from aeroelast.models.blade.numad.mesh_gen.surface import Surface
from aeroelast.models.blade.numad.utils.interpolation import interpolator_wrap


def _row_counts(n1, n3, rows, k):
    """Chordwise counts ``(inboard, outboard)`` of row ``k`` of a region cut into ``rows`` rows.

    The structured region keeps the larger count ``peak`` on every row and puts the transition in
    the end row next to the coarser edge.
    """
    peak = max(n1, n3)
    if n1 < n3 and k == 0:
        return n1, peak
    if n3 < n1 and k == rows - 1:
        return peak, n3
    return peak, peak


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
        self.web_rows = None  # per native bay: the webs meshed there (set by the planning pass)
        self.cut_webs = None  # per native bay: the webs the spar cap is cut at
        self.patches = {}  # (piece key, native bay, row) -> gmsh patch topology, see ``patch``
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

    def webs(self, i):
        """``(meshed webs, cut webs)`` of the native bay that holds new bay ``i``, or None.

        Where a web starts and ends (S6) and where the cap is cut depend on the stations the
        decision is taken at, so the bound blade, which has more stations, takes the native bay's
        decision instead of deciding again.
        """
        if self._bay is None:
            raise RuntimeError("SpanPlan.bind() has not run")
        if self.web_rows is None:
            return None
        b = self._bay[i][0]
        return self.web_rows[b], self.cut_webs[b]

    def patch(self, key, i):
        """The gmsh patch the native mesh puts on new bay ``i`` for cap piece ``key``, or None.

        ``(keep, eta, elements)`` as ``refine_surface_patch`` records it, with ``eta`` expressed in
        the one-row region of bay ``i``. The gmsh result depends chaotically on millimetres of
        boundary geometry, which differs between the native region (Bezier in the bay) and the
        bound blade (spline through the inserted stations), so the bound blade replays the native
        topology on its own boundary instead of meshing again.
        """
        if self._bay is None:
            raise RuntimeError("SpanPlan.bind() has not run")
        b, k = self._bay[i]
        found = self.patches.get((key, b, k))
        if found is None:
            return None
        keep, eta, elements, count = found
        rows = int(self.n_span[b])
        v0, v1 = -1.0 + 2.0 * k / rows, -1.0 + 2.0 * (k + 1) / rows
        eta = np.array(eta, dtype=float).reshape(-1, 2)
        eta[:, 1] = 2.0 * (eta[:, 1] - v0) / (v1 - v0) - 1.0
        return keep, eta, elements, count

    def counts(self, key, i):
        """Chordwise (inboard, outboard) counts of the row of new bay ``i``, or None."""
        if self._bay is None:
            raise RuntimeError("SpanPlan.bind() has not run")
        b, k = self._bay[i]
        ends = self.chord_counts.get((key, b))
        if ends is None:
            return None
        return _row_counts(*ends, int(self.n_span[b]), k)


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


# Key lines of each side of the section, as rows of ``key_arcs`` (13 lines bound the 12 skin
# segments): the pressure side (HP, side 0) is lines 0..6 and the suction side (LP, side 1) lines
# 6..12, line 6 being the leading edge. Segment ``j`` lies between lines ``j`` and ``j + 1``.
_SIDE_LINES = ((0, 6), (6, 12))
_N_SEGMENTS = 12
_WEB_MERGE = 4.0e-3  # m: two webs closer than this on one side share one attachment line (from c157c62)
_QUAD3_ORDER = (
    (0, 0), (0, 3), (3, 3), (3, 0), (0, 1), (0, 2), (1, 3), (2, 3),
    (3, 2), (3, 1), (2, 0), (1, 0), (1, 1), (1, 2), (2, 2), (2, 1),
)  # (row, column) of the 16 "quad3" control points, in ShellRegion order


def _quad3_points(columns, row0):
    """The 16 control points of a "quad3" region from 4 chordwise columns, rows ``row0..row0+3``."""
    return np.array([columns[c][row0 + r] for r, c in _QUAD3_ORDER])


def _snapped_arcs(key_arcs, web_arcs, snap, free=()):
    """``(arcs, onto)``: the web arcs with those closer than ``snap`` to a key line of their side
    moved onto it, and the mask ``(n_webs, 2, n_stations)`` of the arcs that moved (or sit on it).
    Stations in ``free`` never snap."""
    key_arcs = np.asarray(key_arcs, dtype=float)
    out = np.array(web_arcs, dtype=float)
    onto = np.zeros(out.shape, dtype=bool)
    stations = np.arange(key_arcs.shape[1])
    for g in range(len(out)):
        for side, (first, last) in enumerate(_SIDE_LINES):
            lines = key_arcs[first : last + 1]
            line = lines[np.argmin(np.abs(lines - out[g, side]), axis=0), stations]
            close = np.abs(line - out[g, side]) < snap
            if len(free):
                close[list(free)] = False
            out[g, side] = np.where(close, line, out[g, side])
            onto[g, side] = close
    return out, onto


def _crossings(key_arcs, arcs, web_rows):
    """``(web, side, line, row, offset at the row's two stations)`` of every key line a web is on
    opposite sides of at the two stations of a span row where it exists."""
    found = []
    for g in range(len(arcs)):
        for side, (first, last) in enumerate(_SIDE_LINES):
            for k in range(first, last + 1):
                offset = arcs[g, side] - key_arcs[k]
                found.extend(
                    (g, side, k, row, offset[row : row + 2])
                    for row in sorted(web_rows[g])
                    if offset[row] * offset[row + 1] < 0.0
                )
    return found


def snap_web_arcs(key_arcs, web_arcs, snap, web_rows, free=()):
    """The arc each web attaches at on every section: the yaml arc, or a key line close to it.

    ``key_arcs`` is ``(13, n_stations)`` (``keypoints.key_arcs``) and ``web_arcs`` ``(n_webs, 2,
    n_stations)`` (``keypoints.web_attach_arcs``, pressure side then suction side, same axis).
    A web whose arc is closer than ``snap`` to a key line of its side at a station is moved onto
    that line there (never by more than ``snap``); every other arc is returned untouched. Stations
    in ``free`` are never snapped.

    The skin is meshed in segments bounded by the key lines, and a web that crosses a line inside a
    span row (``web_rows[g]``, the rows where web ``g`` exists) is in two segments at once: neither
    piece of the cap there is a region the mesher can build. After snapping, a web that is still on
    opposite sides of a line at the two stations of a row cannot be represented without moving it
    by more than ``snap``; that is an error, never a silent clip (``web_snap_bays`` finds where
    stations make it representable).
    """
    key_arcs = np.asarray(key_arcs, dtype=float)
    out, _ = _snapped_arcs(key_arcs, web_arcs, snap, free)
    for g, side, line, row, offset in _crossings(key_arcs, out, web_rows):
        raise ValueError(
            f"shear web {g} ({('HP', 'LP')[side]} side) crosses key line {line} between stations "
            f"{row} and {row + 1} ({offset[0]:+.4f} m, {offset[1]:+.4f} m from it); it is not within "
            f"the snap distance {snap:.4f} m of the line at both, so it cannot be meshed without "
            "moving it further"
        )
    return out


# A bay where a web leaves a key line is this many snap distances long (``web_snap_spans``): the
# strip that opens in it has an apex angle of at least atan(1 / factor), and the bay is long enough
# for the transition elements of the pieces next to it. The relative slack of that length follows.
_TRANSITION_BAY = 2.0
_SNAP_BAY_SLACK = 1.0e-6


def _snap_rows(blade, elementSize, min_strip_ratio):
    """``(crossing rows, transitions, snap distance, stations z, events)`` of the webs of ``blade``.

    A crossing row has a web on opposite sides of a key line at its two stations, snapped at
    neither. A transition ``(row, end)`` has a web snapped at one station of the row (``end`` 0:
    the first, 1: the last) and not at the other. ``events`` lists them for ``web_snap_spans`` as
    ``(row, end, far)``: ``far`` is the z of the station where the stretch the web is snapped over
    ends at its other side (None if it does not end within the web rows).
    """
    web_attach = getattr(blade.keypoints, "web_attach_arcs", None)
    n_webs = len(blade.stackdb.swstacks)
    z = np.asarray(blade.ispan, dtype=float)
    snap = min_strip_ratio * elementSize
    if web_attach is None or len(web_attach) != n_webs:
        return set(), set(), snap, z, []
    exists = _web_rows_in_yaml(blade)
    web_rows = [{r for r, e in enumerate(exists) if g in e} for g in range(n_webs)]
    key_arcs = np.asarray(blade.keypoints.key_arcs, dtype=float)
    arcs, onto = _snapped_arcs(key_arcs, np.asarray(web_attach, dtype=float), snap)
    crossing = {row for *_, row, _ in _crossings(key_arcs, arcs, web_rows)}
    transition, events = set(), []
    for g in range(n_webs):
        for side in range(len(_SIDE_LINES)):
            found = []
            for row in sorted(web_rows[g]):
                first, last = onto[g, side, row], onto[g, side, row + 1]
                if last and not first:
                    found.append((row, 1))
                elif first and not last:
                    found.append((row, 0))
            transition.update(found)
            for k, (row, end) in enumerate(found):
                # a stretch entered at ``end`` 1 ends at the next transition, one left at 0 began at the last
                other = found[k + 1] if end == 1 and k + 1 < len(found) else found[k - 1] if end == 0 and k else None
                far = None if other is None else float(z[other[0] + other[1]])
                events.append((row, end, far))
    return crossing, transition, snap, z, events


def _transition_bays(events, z, snap, known):
    """``(merged intervals, inner stations)`` the web transitions want as bays.

    A transition keeps one of its two stations and wants a new one ``_TRANSITION_BAY`` snap
    distances from it. It keeps the unsnapped station (the web is where the yaml puts it there)
    when that is already a station of the blade (``known``), so a station that exists next to the
    transition is used instead of cutting a short row, and the snapped station otherwise. When the
    stretch the web is snapped over is shorter than that, the new station is the middle of the
    stretch (both ends of a crossing web share it). Transitions of different webs whose bays
    overlap or are less than a snap distance apart share one bay: stations of one inside the bay
    of another would cut it into rows too short for the elements of the pieces next to them. The
    stations inside a merged bay that sit in a snapped stretch are returned too, one per snap
    distance, so that every web sees a snapped station.
    """
    reach = _TRANSITION_BAY * snap
    known = np.sort(np.asarray(known, dtype=float))
    wanted, inner = [], []
    for row, end, far in events:
        snapped = float(z[row + end])
        sign = 1.0 if end == 1 else -1.0  # from the unsnapped station toward the snapped one
        outer = snapped - sign * reach  # the unsnapped station, one bay away
        keep, new = snapped, outer
        between = known[((known - outer) * sign > 1.0e-9) & ((snapped - known) * sign > 1.0e-9)]
        if between.size:  # a station already in the bay would cut it: start the bay there
            keep = float(between[np.argmin(np.abs(between - outer))])
            new = keep + sign * reach
        if far is not None and abs(far - snapped) < reach:
            new = 0.5 * (snapped + far)
            inner.append(new)
            keep = outer
        wanted.append((min(keep, new), max(keep, new)))
    wanted.sort()
    merged = []
    for lo, hi in wanted:
        if merged and lo <= merged[-1][1] + snap:
            merged[-1][1] = max(merged[-1][1], hi)
        else:
            merged.append([lo, hi])
    # Transitions at different places of a merged bay leave it longer than one transition wants
    # (``web_snap_bays`` flags such a row, a web snapped at one end and not at the other): split
    # every gap between the stations of the bay (its ends and the stations that already exist in
    # it) that is longer than a transition bay into equal rows no longer than that.
    longest = reach * (1.0 + 0.5 * _SNAP_BAY_SLACK)
    for lo, hi in merged:
        edges = np.concatenate(([lo], known[(known > lo + 1.0e-9) & (known < hi - 1.0e-9)], [hi]))
        for a, b in zip(edges[:-1], edges[1:], strict=True):
            parts = int(np.ceil((b - a) / longest))
            inner.extend(a + (b - a) * np.arange(1, parts) / parts)
    return merged, sorted(inner)


def web_snap_bays(blade, elementSize, min_strip_ratio=0.3):
    """Span rows (bays) that need more stations before the webs can be planned (S8).

    A web within ``min_strip_ratio * elementSize`` of a key line snaps onto it (``snap_web_arcs``).
    The rows that make that impossible or ugly are

    * a row where the web is on opposite sides of a key line at its two stations and snapped at
      neither (it would cross the line), and
    * a row where the web is snapped at one station and not at the other and which is longer than
      ``_TRANSITION_BAY`` snap distances: the skin strip between the web and the line opens from
      zero to at least the snap distance over the row, a wedge as sharp as the row is long.

    Returns the row indices (empty when the stations suffice).
    """
    crossing, transition, snap, z, _ = _snap_rows(blade, elementSize, min_strip_ratio)
    reach = _TRANSITION_BAY * snap * (1.0 + _SNAP_BAY_SLACK)
    return np.array(sorted(crossing | {r for r, _ in transition if z[r + 1] - z[r] > reach}), dtype=int)


def web_crossings(blade, elementSize, min_strip_ratio=0.3):
    """Span rows where a web still crosses a key line after snapping (they cannot be meshed)."""
    crossing, *_ = _snap_rows(blade, elementSize, min_strip_ratio)
    return np.array(sorted(crossing), dtype=int)


def web_snap_spans(blade, elementSize, min_strip_ratio=0.3, known=()):
    """The stations that bound the bay of each web transition, as a sorted array of z (S8).

    Where a web is snapped at one station of a row and not at the other (``web_snap_bays`` has
    refined the row until the transition is bracketed), the bay of the transition is
    ``_TRANSITION_BAY`` snap distances long and starts at the station the web keeps (see
    ``_transition_bays``; ``known`` are the z of the stations of the real blade, which the probe
    this runs on may have more of). Overlapping transitions share a bay, and a stretch where a web
    is snapped that is too short to hold its bays gets its station in the middle, so the snap is
    seen at a station of the bay.
    """
    _, _, snap, z, events = _snap_rows(blade, elementSize, min_strip_ratio)
    merged, inner = _transition_bays(events, z, snap, np.asarray(known, dtype=float))
    return np.unique(np.concatenate((np.array(merged, dtype=float).ravel(), inner)))


def _cuts_segment(arc, lo, hi, row):
    """Whether a line at ``arc`` splits the segment ``lo..hi`` in span row ``row``: it leaves some
    width of the segment on each side at one of the two stations of the row."""
    both = slice(row, row + 2)
    return bool(
        np.max(arc[both] - lo[both]) > _PIECE_MIN_WIDTH and np.max(hi[both] - arc[both]) > _PIECE_MIN_WIDTH
    )


def plan_cap_segments(key_arcs, web_arcs, snap, web_rows, free=()):
    """Arc of every column of each skin segment a shear web cuts, per station.

    Returns ``{segment: (arcs (n + 2, n_stations), {web: column - 1})}`` for every segment some web
    cuts in a row where it exists: the arcs are the lower key line, one line per web cutting the
    segment (sorted by arc) and the upper key line. A web is cut into whichever segment holds it
    (spar cap, leading/trailing edge reinforcement or panel) at the arc the yaml gives, snapped
    onto a key line closer than ``snap`` (``snap_web_arcs``), so its arc is clipped to the segment
    where it lies in a neighbour. Two webs closer than ``_WEB_MERGE`` become one line. ``None`` when two
    webs cross inside a segment (the mesher then keeps the webs on the cap edges).
    """
    key_arcs = np.asarray(key_arcs, dtype=float)
    snapped = snap_web_arcs(key_arcs, web_arcs, snap, web_rows, free)
    plans = {}
    for seg in range(_N_SEGMENTS):
        side = 0 if seg < _SIDE_LINES[0][1] else 1
        lo, hi = key_arcs[seg], key_arcs[seg + 1]
        members = [
            g
            for g in range(len(snapped))
            if any(_cuts_segment(snapped[g, side], lo, hi, row) for row in web_rows[g])
        ]
        if not members:
            continue
        arcs = np.clip(snapped[members, side], lo, hi)
        order = np.argsort(arcs.mean(axis=1), kind="stable")
        ordered = arcs[order]
        for r in range(len(ordered) - 1):
            gap = ordered[r + 1] - ordered[r]
            if np.any(gap < -_WEB_MERGE):
                return None
            close = gap < _WEB_MERGE
            mean = 0.5 * (ordered[r] + ordered[r + 1])
            ordered[r] = np.where(close, mean, ordered[r])
            ordered[r + 1] = np.where(close, mean, ordered[r + 1])
        plans[seg] = (np.vstack((lo, ordered, hi)), {members[k]: i for i, k in enumerate(order)})
    return plans


def _web_control_points(shellKp):
    """Fill the interior Bezier control points of a web quad (rows 4, 5, 8, 9, 12-15) in place.

    They sit at exact thirds of the straight edges; the original 0.6666/0.3333 weights summed to
    0.9999 and shifted every control point by 1e-4 |P| towards the origin.
    """
    for lo, hi, a, b in ((4, 5, 0, 1), (8, 9, 2, 3), (12, 13, 11, 6), (14, 15, 7, 10)):
        shellKp[lo, :] = (2.0 * shellKp[a, :] + shellKp[b, :]) / 3.0
        shellKp[hi, :] = (shellKp[a, :] + 2.0 * shellKp[b, :]) / 3.0


_PIECE_MIN_WIDTH = 1.0e-6  # m: narrower than this a cap piece is a line, and a strip is closed


def _piece_exists(arcs, lo, hi, row):
    """Whether the cap piece between columns ``lo`` and ``hi`` has any width in span row ``row``."""
    return np.max(arcs[hi, row : row + 2] - arcs[lo, row : row + 2]) > _PIECE_MIN_WIDTH


def _closing_strips(cap_plans, webs, row, station, step, wide, stop=frozenset()):
    """The cap strips that close to a point at ``station`` and widen to ``wide`` moving by ``step``.

    ``cap_plans`` holds one ``(arcs, web positions)`` per spar cap (``plan_cap_segments``); the
    strips are the pieces of each cap cut at ``webs`` that exist in span row ``row``. A strip closes
    at ``station`` when its width there is zero; it reaches ``wide`` at the first station, moving by
    ``step``, where it is at least that wide, or at a station listed in ``stop`` (stations placed
    where the width crosses ``wide``, see ``web_end_spans``). A strip that stops widening first does
    not close and is left out. Returns ``(previous station, station reached, widths)`` per strip.
    """
    out = []
    for arcs, pos in cap_plans:
        cuts = [0, *sorted(1 + pos[g] for g in webs if g in pos), len(arcs) - 1]
        for lo, hi in zip(cuts[:-1], cuts[1:], strict=True):
            width = arcs[hi] - arcs[lo]
            if not _piece_exists(arcs, lo, hi, row) or width[station] > _PIECE_MIN_WIDTH:
                continue
            at = station
            while (
                0 <= at + step < len(width)
                and at not in stop
                and width[at] < wide
                and width[at + step] - width[at] > _PIECE_MIN_WIDTH
            ):
                at += step
            if at != station and (at in stop or width[at] >= wide):
                out.append((at - step, at, width))
    return out


def _cap_plans(blade, n_webs, snap):
    """``(arcs, web positions)`` of each skin segment a web cuts (``plan_cap_segments``), or ``[]``."""
    web_attach = getattr(blade.keypoints, "web_attach_arcs", None)
    if web_attach is None or len(web_attach) != n_webs:
        return []
    exists = _web_rows_in_yaml(blade)
    web_rows = [{r for r, e in enumerate(exists) if g in e} for g in range(n_webs)]
    plans = plan_cap_segments(blade.keypoints.key_arcs, np.asarray(web_attach, dtype=float), snap, web_rows)
    return [] if plans is None else list(plans.values())


def _web_rows_in_yaml(blade):
    """The webs the yaml gives a laminate in each span row."""
    swstacks = blade.stackdb.swstacks
    n_rows = len(blade.ispan) - 1
    return [{g for g in range(len(swstacks)) if swstacks[g][i].plygroups} for i in range(n_rows)]


def web_end_spans(blade, elementSize, min_strip_ratio=0.3):
    """Span positions where a web starts or ends (S6 of v1-mesher-web-slivers), stations or not.

    A web starts (ends) where the last cap strip that closes at the end of its yaml domain reaches
    ``min_strip_ratio * elementSize`` (``get_shell_mesh``). That place is between two stations; it
    is found by linear interpolation of the strip width between them. Adding these positions as
    stations makes the web start there whatever the station density, so a mesh planned on the yaml
    stations and one on a blade with more stations agree. Pass them back to ``get_shell_mesh`` as
    ``web_end_spans`` so the station is taken even where interpolation leaves it a hair short.
    """
    exists = _web_rows_in_yaml(blade)
    wide = min_strip_ratio * elementSize
    plans = _cap_plans(blade, len(blade.stackdb.swstacks), wide)
    if not plans:
        return np.empty(0)
    z = np.asarray(blade.ispan, dtype=float)
    found = []
    n_rows = len(exists)
    for row in range(n_rows):
        ends = []
        if exists[row] - (exists[row - 1] if row else set()):
            ends.append((row, 1))
        if exists[row] - (exists[row + 1] if row + 1 < n_rows else set()):
            ends.append((row + 1, -1))
        for station, step in ends:
            strips = _closing_strips(plans, exists[row], row, station, step, wide)
            if not strips:
                continue
            crossings = []
            for prev, at, width in strips:
                t = (wide - width[prev]) / (width[at] - width[prev])
                crossings.append(z[prev] + t * (z[at] - z[prev]))
            found.append(max(crossings) if step > 0 else min(crossings))
    return np.unique(np.asarray(found, dtype=float))


def _meshed_web_rows(exists, apex_station):
    """The webs meshed in each span row: those the yaml gives a laminate, less the rows at each end
    of a web's domain where a skin strip closes to a point.

    ``exists[row]`` is the set of webs with a laminate in the row. ``apex_station(webs, row,
    station, step)`` is the station, moving from ``station`` by ``step`` (+1 inward from a start,
    -1 inward from an end), where the strips that close to a point at ``station`` are wide enough.
    A web that starts (ends) in ``row`` starts (ends) there instead.
    """
    n_rows = len(exists)
    kept = [set(row) for row in exists]
    for row in range(n_rows):
        starting = exists[row] - (exists[row - 1] if row else set())
        if starting:
            begin = apex_station(exists[row], row, row, 1)
            for k in range(row, min(begin, n_rows)):
                kept[k] -= starting
        ending = exists[row] - (exists[row + 1] if row + 1 < n_rows else set())
        if ending:
            end = apex_station(exists[row], row, row + 1, -1)
            for k in range(max(end, 0), row + 1):
                kept[k] -= ending
    return kept


def _cut_web_sets(web_rows, yaml_rows, plans, wide):
    """The webs each skin segment is cut at in each span row.

    ``plans`` is ``{segment: (arcs, web positions)}`` (``plan_cap_segments``); the result has, for
    every row, ``{segment: set of webs}``. A segment is cut at the webs meshed in the row that
    cut it. Where a web is not meshed yet (S6 starts it later because of the strips of the spar
    cap) its segments are still cut from the first row the yaml gives it a laminate, provided every
    piece that cut leaves is at least ``wide`` across at both stations of the row: the strips of
    that segment do not close, so there is no apex to start the web after, and the cut then does
    not begin in a row too short to carry the transition of the nodes. Two rows that share a
    station must agree on it, per segment: either one cuts at every line the other does (the
    coarser row is then remeshed to the finer row's nodes, see ``_cap_rows``) or, if neither
    contains the other, both cut at the union.
    """
    out = [{} for _ in web_rows]
    for seg, (arcs, pos) in plans.items():
        sets = []
        for r, row in enumerate(web_rows):

            def cutting(webs, r=r):
                return {g for g in webs if g in pos and _cuts_segment(arcs[1 + pos[g]], arcs[0], arcs[-1], r)}

            meshed = cutting(row)
            planned = cutting(yaml_rows[r])
            if planned != meshed:
                edges = [0, *sorted(1 + pos[g] for g in planned), len(arcs) - 1]
                widths = [arcs[hi, r : r + 2] - arcs[lo, r : r + 2] for lo, hi in zip(edges[:-1], edges[1:], strict=True)]
                if all(np.all(w >= wide) or np.all(w <= _PIECE_MIN_WIDTH) for w in widths):
                    meshed = planned
            sets.append(meshed)
        changed = True
        while changed:
            changed = False
            for r in range(len(sets) - 1):
                if not (sets[r] <= sets[r + 1] or sets[r + 1] <= sets[r]):
                    sets[r] = sets[r + 1] = sets[r] | sets[r + 1]
                    changed = True
        for r, cut in enumerate(sets):
            out[r][seg] = cut
    return out


def _edge_nodes(region, count, v):
    """3D nodes of the chordwise edge ``v = -1`` (inboard) or ``+1`` of a region: count + 1, uniform in eta."""
    u = np.linspace(-1.0, 1.0, count + 1)
    return region.XYZCoord(np.column_stack((u, np.full_like(u, v))))


def _side_nodes_between(region, u, v0, v1):
    """3D nodes of the spanwise edge ``u = -1`` or ``+1`` of a region between ``v0`` and ``v1``."""
    v = np.array([v0, v1])
    return region.XYZCoord(np.column_stack((np.full_like(v, u), v)))


def _is_line_of(imposed, line):
    """Whether the nodes ``imposed`` (a finer neighbour's, with the collapsed edges of pieces that
    are zero wide there repeating a node) are exactly the nodes ``line`` already has."""
    extent = float(np.ptp(line, axis=0).max())
    keep = [0]
    for k in range(1, len(imposed)):
        if np.linalg.norm(imposed[k] - imposed[keep[-1]]) > _PIECE_MIN_WIDTH * extent:
            keep.append(k)
    unique = imposed[keep]
    return len(unique) == len(line) and bool(np.allclose(unique, line, rtol=0.0, atol=_PIECE_MIN_WIDTH * extent))


def _cap_rows(
    shellKp,
    counts,
    n_span,
    element_size,
    min_strip_ratio,
    below=None,
    above=None,
    record=None,
    replay=None,
):
    """Mesh one cap piece, remeshing with gmsh only the rows that need it; ``None`` if none does.

    The piece is the structured grid of ``n_span`` rows (``ShellRegion``) except that a row is
    replaced by a gmsh patch, its boundary nodes exactly those of the structured row, when
    - it is the first (last) row and ``below`` (``above``) imposes on its inboard (outboard) edge
      the nodes of a neighbouring bay that cuts the cap at more lines (3D, ordered like the edge), or
    - its two chordwise counts differ and the finer one spaces its nodes closer than a strip may
      be wide (``min_strip_ratio * element_size``): the transition would fan out into slivers.
    The size next to the boundary follows each boundary segment's own length (``refine_patch``),
    so a long imposed edge does not end in needles. Patching rows, not the
    bay, gives the same discretisation whether the bay is one native bay of ``n_span`` rows or
    ``n_span`` bays of one row (``SpanPlan``): ``record`` (a dict) receives the gmsh topology of
    each patched row, and ``replay`` ({row: topology}, e.g. ``SpanPlan.patch``) patches exactly
    those rows with that topology instead of deciding and meshing.
    """
    n1, n3 = counts
    region = ShellRegion(
        "quad3", shellKp, [n1, n_span, n3, n_span], elType="quad", meshMethod="structured"
    )
    patches = {}
    for k in range(n_span):
        c1, c3 = _row_counts(n1, n3, n_span, k)
        v0, v1 = -1.0 + 2.0 * k / n_span, -1.0 + 2.0 * (k + 1) / n_span
        low, high = _edge_nodes(region, c1, v0), _edge_nodes(region, c3, v1)
        width = (np.linalg.norm(low[-1] - low[0]), np.linalg.norm(high[-1] - high[0]))
        spacing = [w / c for w, c in zip(width, (c1, c3), strict=True) if w > 0.0]
        fan = c1 != c3 and min(spacing, default=np.inf) < min_strip_ratio * element_size
        imposed = []
        if replay is not None and k not in replay:
            continue
        if k == 0 and below is not None and not _is_line_of(below, low):
            low = below
            imposed.append(True)
        if k == n_span - 1 and above is not None and not _is_line_of(above, high):
            high = above
            imposed.append(True)
        if replay is None and not (imposed or fan):
            continue
        right = _side_nodes_between(region, 1.0, v0, v1)
        left = _side_nodes_between(region, -1.0, v0, v1)
        loop = np.vstack((low[:-1], right[:-1], high[::-1][:-1], left[::-1][:-1]))
        found = [] if record is not None else None
        patches[k] = refine_surface_patch(
            loop,
            region.XYZCoord,
            element_size,
            recombine=True,
            record=found,
            replay=None if replay is None else replay[k],
        )
        if found:
            record[k] = found[0]
    if not patches:
        return None
    nodes, elements = [], []
    if len(patches) < n_span:
        base = region.createShellMesh()
        row = np.arange(len(base["elements"])) // max(n1, n3)
        kept = base["elements"][~np.isin(row, list(patches))]
        used = np.unique(kept[kept >= 0])
        renumber = -np.ones(len(base["nodes"]), dtype=int)
        renumber[used] = np.arange(len(used))
        nodes.append(base["nodes"][used])
        elements.append(np.where(kept >= 0, renumber[np.where(kept >= 0, kept, 0)], -1))
    for patch_nodes, patch_elements in patches.values():
        offset = sum(len(n) for n in nodes)
        nodes.append(patch_nodes)
        elements.append(np.where(patch_elements >= 0, patch_elements + offset, -1))
    return MeshedRegion(np.vstack(nodes), np.vstack(elements))


def get_shell_mesh(
    blade,
    elementSize,
    spanGrading="chord",
    include_webs=True,
    span_plan=None,
    plan_only=False,
    min_strip_ratio=0.3,
    web_end_spans=(),
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
    min_strip_ratio: float
        A strip of the spar cap (between two webs, or between a web and the cap edge) that closes
        to a point at an end of the yaml web domain makes a sliver at its apex whatever the mesh:
        any element there is as sharp as the wedge. Each web therefore starts and ends, moving
        inward from the ends of its yaml domain, at the first station where those strips are at
        least ``min_strip_ratio * elementSize`` wide, which makes the corner about a right angle.
        A strip that does not close (it stays narrow) does not move the web; the arcs themselves
        are never moved.
    web_end_spans: sequence of float
        Stations of the blade placed where a closing strip reaches the minimum width
        (``web_end_spans()``); the web starts (ends) at such a station even if the interpolated
        width there falls a rounding error short.

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

    ## Web attachment lines and the skin segments they cut.
    # Each web hangs on a line at the yaml's attachment arc (``keypoints.web_attach_arcs``), in
    # whichever skin segment holds it (spar cap, leading/trailing edge reinforcement or panel):
    # the segment is split there into sub-regions A, B, C, ... that keep its laminate, so the web
    # edge shares its nodes with the skin. A web closer than the minimum strip width to a key line
    # is moved onto the line (``snap_web_arcs``), never by more; one that crosses a line any other
    # way raises. Without yaml arcs the webs stay on the cap edges (key points ``b`` and ``c``),
    # as NuMAD has them.
    def _col(index):
        return np.stack((splineXi[:, index], splineYi[:, index], splineZi[:, index]), axis=1)

    def _arc_column(arcs, anchors=()):
        """Longitudinal line through the points at arc ``arcs[k]`` of every section ``k``.

        ``anchors`` are ``(arcs, column)`` pairs of lines the result may coincide with: where its
        arc equals theirs the section point is theirs exactly (the profile polyline and the key
        points differ by rounding), so a web on a key line shares its nodes with the skin.
        """
        pts = np.empty((numXsec, 3))
        for k in range(numXsec):
            held = next((col[3 * k] for a, col in anchors if arcs[k] == a[k]), None)
            if held is not None:
                pts[k] = held
                continue
            arc_k = geometry.arclength[1:-1, k]
            for d in range(3):
                pts[k, d] = np.interp(arcs[k], arc_k, profile0[1:-1, d, k])
        return np.stack(
            [interpolator_wrap(spParam, pts[:, d], spParami, "pchip") for d in range(3)], axis=1
        )

    # The cap split does not depend on ``include_webs``: the no-webs surface mesh must keep
    # the same skin discretisation as the structural mesh, so it is planned either way.
    def _web_in_row(group, row):
        """Whether shear web ``group`` is meshed in span row ``row``."""
        return group in webRows[row]

    snapDistance = min_strip_ratio * elementSize
    capSplit = {}  # segment -> (arcs, columns, intermediate columns, web group positions)
    webArcs = None  # the arcs the webs attach at, snapped onto the key lines they are close to
    onCapEdges = True  # no yaml arcs (or webs that cross): the webs stay on the cap edges
    web_attach = getattr(blade.keypoints, "web_attach_arcs", None)
    yamlRows = _web_rows_in_yaml(blade)
    if web_attach is not None and len(web_attach) == len(swstacks):
        keyArcs = blade.keypoints.key_arcs
        rowsOfWeb = [{r for r, e in enumerate(yamlRows) if g in e} for g in range(len(swstacks))]
        webArcs = snap_web_arcs(keyArcs, np.asarray(web_attach, dtype=float), snapDistance, rowsOfWeb, ())
        plans = plan_cap_segments(keyArcs, np.asarray(web_attach, dtype=float), snapDistance, rowsOfWeb, ())
        if plans is None:
            warnings.warn(
                "Shear-web attachment arcs cross inside a skin segment; "
                "the webs stay on the cap edges.",
                stacklevel=2,
            )
        else:
            onCapEdges = False
            for seg, (arcs, pos) in plans.items():
                edges = [(arcs[0], _col(3 * seg)), (arcs[-1], _col(3 * seg + 3))]
                cols = [edges[0][1]] + [_arc_column(a, edges) for a in arcs[1:-1]] + [edges[1][1]]
                mids = {}  # (lo, hi) column pair -> the two intermediate columns, built when needed
                capSplit[seg] = (arcs, cols, mids, pos)

    def _apex_station(webs, row, station, step):
        """Where the cap strips that close to a point at ``station`` reach the minimum width.

        A strip (a piece of the cap cut at ``webs``) closes to a point at the end of a web's domain
        when its width is zero there and keeps widening, station after station, until it is
        ``min_strip_ratio * elementSize`` wide; the station where it gets there is where the web
        may start (end). A strip that stops widening before that (one held at a constant width,
        e.g. by the margin that keeps a web arc off the cap edge) does not close and is ignored.
        """
        plans = [(arcs, pos) for arcs, _, _, pos in capSplit.values()]
        reached = [
            at
            for _, at, _ in _closing_strips(
                plans, webs, row, station, step, min_strip_ratio * elementSize, endStations
            )
        ]
        if not reached:
            return station
        return max(reached) if step > 0 else min(reached)

    _ispan = np.asarray(blade.ispan, dtype=float)
    endStations = frozenset(
        int(k)
        for k in np.flatnonzero(
            np.isclose(_ispan[:, None], np.asarray(web_end_spans, dtype=float)[None, :], rtol=0.0, atol=1.0e-9).any(axis=1)
        )
    ) if len(web_end_spans) else frozenset()

    def _web_line(group, side, row):
        """The longitudinal line web ``group`` attaches along on ``side`` in span row ``row``.

        It is the cut line of the skin segment that holds the web there, so the web edge and the
        skin share nodes; a web that lies on a key line along the whole row has that line.
        """
        if onCapEdges:
            return (webHP, webLP)[side][group]
        for seg, (arcs, cols, _, pos) in capSplit.items():
            if (seg >= _SIDE_LINES[0][1]) == bool(side) and group in pos:
                if _cuts_segment(arcs[1 + pos[group]], arcs[0], arcs[-1], row):
                    return cols[1 + pos[group]]
        both = slice(row, row + 2)
        for line in range(_SIDE_LINES[side][0], _SIDE_LINES[side][1] + 1):
            if np.all(webArcs[group, side, both] == blade.keypoints.key_arcs[line, both]):
                return _col(3 * line)
        raise RuntimeError(f"web {group} is neither cut into a skin segment nor on a key line in row {row}")

    webRows = _meshed_web_rows(
        yamlRows,
        _apex_station if capSplit else lambda webs, row, station, step: station,
    )
    cutWebs = _cut_web_sets(webRows, yamlRows, {seg: (c[0], c[3]) for seg, c in capSplit.items()}, snapDistance)
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
        record.web_rows, record.cut_webs = webRows, cutWebs
    elif span_plan is not None:
        span_plan.bind(splineZi[0::3, 0])
        nSpan = np.ones(nSections, dtype=int)
        if span_plan.web_rows is not None:
            # the web start/end and the cap cuts were decided at the native stations
            bound = [span_plan.webs(i) for i in range(nSections)]
            webRows = [set(w) for w, _ in bound]
            cutWebs = [{seg: set(c) for seg, c in cut.items()} for _, cut in bound]

    bladeSurf = Surface()
    ## Outer shell sections
    outShES = set()
    secList = []
    rowPieces = {}  # (row, segment) -> [(lo, hi, name, control points, inboard count, outboard count)]
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
            pieces = [(0, 0, stacks[j, i].name, shellKp)]
            if j in capSplit:
                arcs, cols, mids, pos = capSplit[j]
                # The cap is cut only at the webs meshed in this row (same predicate as the web
                # loop below); where none is, it stays one piece with the unsplit name.
                cuts = sorted(1 + pos[g] for g in cutWebs[i].get(j, ()))
                edges = [0, *cuts, len(cols) - 1]
                pieces = []
                for lo, hi in zip(edges[:-1], edges[1:], strict=True):
                    if not _piece_exists(arcs, lo, hi, i):
                        continue
                    if (lo, hi) not in mids:
                        mids[lo, hi] = [
                            _arc_column(arcs[lo] + (arcs[hi] - arcs[lo]) * t)
                            for t in (1.0 / 3.0, 2.0 / 3.0)
                        ]
                    name = stacks[j, i].name + (f"_{chr(ord('A') + lo)}" if len(edges) > 2 else "")
                    pieces.append(
                        (lo, hi, name, _quad3_points([cols[lo], *mids[lo, hi], cols[hi]], stPt))
                    )
            # Per-edge chordwise counts (consistent nSpan spanwise)
            rowPieces[i, j] = [
                (
                    lo,
                    hi,
                    pieceName,
                    shellKp,
                    *_chordwise_counts(
                        shellKp, elementSize, ("shell", j, lo, hi), i, i == rws - 2, span_plan, record
                    ),
                )
                for lo, hi, pieceName, shellKp in pieces
            ]
            stSp = stSp + 3
        stPt = stPt + 3

    ## Emit the outer shell regions, row by row.
    # A piece is meshed as a structured grid unless it needs the local remesh: it sits next to a
    # row that cuts the cap at more lines (the station line must carry that row's nodes), or its
    # two chordwise counts differ and the finer one spaces its nodes closer than a strip may be
    # wide (``min_strip_ratio * elementSize``): the transition would fan out into slivers.
    def _finer_line(j, lo, hi, row, top):
        """The nodes the finer ``row`` puts on the station line between the arcs of columns lo..hi."""
        line = []
        for flo, fhi, _, kp, n1, n3 in rowPieces[row, j]:
            if flo >= lo and fhi <= hi:
                region = ShellRegion("quad3", kp, [n1, 1, n3, 1])
                line.append(_edge_nodes(region, n3 if top else n1, 1.0 if top else -1.0))
        return np.vstack(line)

    for i in range(rws - 1):
        for j in range(12):
            for lo, hi, pieceName, shellKp, nEl1, nEl3 in rowPieces[i, j]:
                nEl = np.array([nEl1, nSpan[i], nEl3, nSpan[i]])
                patch = None
                if j in capSplit:
                    below = above = None
                    if i > 0 and cutWebs[i].get(j, set()) < cutWebs[i - 1].get(j, set()):
                        below = _finer_line(j, lo, hi, i - 1, True)
                    if i < rws - 2 and cutWebs[i].get(j, set()) < cutWebs[i + 1].get(j, set()):
                        above = _finer_line(j, lo, hi, i + 1, False)
                    replay, topology = None, {}
                    if span_plan is not None:
                        found = span_plan.patch(("cap", j, lo, hi), i)
                        replay = {} if found is None else {0: found}
                    patch = _cap_rows(
                        shellKp,
                        (nEl1, nEl3),
                        nSpan[i],
                        elementSize,
                        min_strip_ratio,
                        below=below,
                        above=above,
                        record=topology,
                        replay=replay,
                    )
                    if record is not None:
                        for k, found in topology.items():
                            record.patches[("cap", j, lo, hi), i, k] = found
                if patch is None:
                    bladeSurf.addShellRegion(
                        "quad3",
                        shellKp,
                        nEl,
                        name=pieceName,
                        elType="quad",
                        meshMethod="structured",
                    )
                else:
                    bladeSurf.shellRegions.append(patch)
                    bladeSurf.regionNames.append(pieceName)
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

    ## Shear web sections
    # ``include_webs=False`` collapses the loop range so the whole block is
    # skipped: ``swES``/``secList`` stay outer-shell only and the web node sets
    # are never created downstream.
    swES = set()
    stPt = 0
    np.array([])
    np.array([])
    for i in range(rws - 1 if include_webs else 0):
        if _web_in_row(0, i):
            hp0, lp0 = _web_line(0, 0, i), _web_line(0, 1, i)
            shellKp = np.zeros((16, 3))
            shellKp[0, :] = hp0[stPt]
            shellKp[1, :] = lp0[stPt]
            shellKp[2, :] = lp0[stPt + 3]
            shellKp[3, :] = hp0[stPt + 3]
            shellKp[6, :] = lp0[stPt + 1]
            shellKp[7, :] = lp0[stPt + 2]
            shellKp[10, :] = hp0[stPt + 2]
            shellKp[11, :] = hp0[stPt + 1]
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
        if _web_in_row(1, i):
            hp1, lp1 = _web_line(1, 0, i), _web_line(1, 1, i)
            shellKp = np.zeros((16, 3))
            shellKp[0, :] = lp1[stPt]
            shellKp[1, :] = hp1[stPt]
            shellKp[2, :] = hp1[stPt + 3]
            shellKp[3, :] = lp1[stPt + 3]
            shellKp[6, :] = hp1[stPt + 1]
            shellKp[7, :] = hp1[stPt + 2]
            shellKp[10, :] = lp1[stPt + 2]
            shellKp[11, :] = lp1[stPt + 1]
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
