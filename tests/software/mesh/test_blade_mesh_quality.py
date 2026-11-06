"""Quality of the legacy (windIO v1) blade shell mesh, production size (es 0.25 m, about 60 s each).

Commit c157c62 split the spar caps at the web attachment arcs so the webs share nodes with the skin.
It split them in every span row, also where the web does not exist: there the two web arcs converge
(IEA-15 v1, eta 0.083) and the cap sub-region between them is 1-2 cm wide, which the structured
mesh turns into outer-shell triangles of min angle 3.9 deg and scaled Jacobian 0.07.

The cap is now cut only in the rows where the web is meshed, so these tests assert on the public
mesh that

* no outer-shell element is a sliver (min interior angle >= 15 deg, scaled Jacobian >= 0.3, the
  "poor" thresholds of ``tools/mesh_quality.py``), and
* no cap sub-region (``*_SPAR_A/B/C``) reaches outside the span interval where the yaml gives the
  shear webs a laminate (read here from the yaml web layer thickness, not from the mesher),
* the shear webs, and with them the cap split, start (and end) at the first station, moving inward
  from the yaml web domain, where the skin strips that close to a point there are at least
  ``min_strip_ratio * element_size`` wide (the apex where a web arc meets the cap edge or the other
  web arc is a sliver whatever the remesh; the station is derived here from the yaml arcs),
* the outer shell is conforming: no edge used by a single element away from the root and tip
  (cutting the cap only inside the web domain changes the node pattern across the station where
  the web starts, and the coordinate merge cannot join what does not coincide), and
* no shear web is displaced from the arc the yaml gives it by more than the snap distance
  (``min_strip_ratio * element_size``) at any span station: a web lies where the yaml puts it, in
  whichever skin piece holds it, and only a web closer than that to a key line moves, onto it
  (S8, ``snap_web_arcs``); a web that crosses a key line between two stations without being close
  enough at both to snap raises instead of being clipped,
* no web has a gap inside its span extent, and
* the local refinement code holds no length literal (the sizes come from ``element_size``, the
  measured strip width and a dimensionless growth ratio).

That the webs themselves stay where the yaml puts them is ``tests/test_blade_web_position.py``.
"""

from __future__ import annotations

import inspect
import re
import sys
from pathlib import Path
from typing import Any, NamedTuple

import numpy as np
import pytest
import yaml

from aeroelast.core.mesh.generators import BladeMesh
from aeroelast.models.blade.numad.mesh_gen.mesh_gen import (
    _snapped_arcs,
    get_shell_mesh,
    plan_cap_segments,
    snap_web_arcs,
)
from tests.support.paths import DATA_DIR

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "tools"))
import audit_blade_input as ab  # noqa: E402

MIN_ANGLE_DEG = 15.0
MIN_SCALED_JACOBIAN = 0.3
ELEMENT_SIZE = 0.25
# S6b: the web end station is placed at the strip-width crossing to this relative accuracy.
CROSSING_REL = 1.0e-3
SPLIT_CAP = re.compile(r"SPAR_[ABC]$")
BLADES = {"iea15": "IEA-15-240-RWT.yaml", "iea34": "IEA-3.4-130-RWT.yaml"}


def _quality(corners):
    """Min interior angle (deg) and scaled Jacobian of one polygon, corners (k, 3), k = 3 or 4."""
    k = len(corners)
    nxt = np.roll(corners, -1, axis=0) - corners
    prv = np.roll(corners, 1, axis=0) - corners
    if k == 4:
        normal = np.cross(corners[2] - corners[0], corners[3] - corners[1])
    else:
        normal = np.cross(corners[1] - corners[0], corners[2] - corners[0])
    normal = normal / np.linalg.norm(normal)
    denom = np.linalg.norm(nxt, axis=1) * np.linalg.norm(prv, axis=1)
    cos = np.clip(np.einsum("ij,ij->i", nxt, prv) / denom, -1.0, 1.0)
    sjac = np.cross(nxt, prv) @ normal / denom
    return float(np.degrees(np.arccos(cos)).min()), float(sjac.min())


def _web_domain_z(yaml_path):
    """Span interval [z_lo, z_hi] (m) where the yaml gives the shear webs a laminate.

    The loader resamples every layer linearly onto the reference-axis grid and zero outside the
    layer's own grid, and the mesher integrates thickness over each bay, so a web exists from the
    grid point before the first positive sample to the grid point after the last one.
    """
    with open(yaml_path) as fh:
        blade = yaml.safe_load(fh)["components"]["blade"]
    axis_z = blade["outer_shape_bem"]["reference_axis"]["z"]
    grid = np.asarray(axis_z["grid"], dtype=float)
    positive = np.zeros(grid.size, dtype=bool)
    for layer in blade["internal_structure_2d_fem"]["layers"]:
        if "web" not in layer:
            continue
        g = np.asarray(layer["thickness"]["grid"], dtype=float)
        inside = (grid >= g[0]) & (grid <= g[-1])
        positive[inside] |= np.interp(grid[inside], g, layer["thickness"]["values"]) > 0.0
    first, last = np.flatnonzero(positive)[[0, -1]]
    z = np.asarray(axis_z["values"], dtype=float)
    return float(z[max(first - 1, 0)]), float(z[min(last + 1, grid.size - 1)])


class Built(NamedTuple):
    name: str
    mesh: Any
    web_domain: tuple
    builder: BladeMesh
    path: Path


@pytest.fixture(scope="module", params=sorted(BLADES))
def blade(request):
    path = DATA_DIR / BLADES[request.param]
    builder = BladeMesh(str(path), element_size=ELEMENT_SIZE)
    mesh = builder.generate(verbose=False)
    return Built(request.param, mesh, _web_domain_z(path), builder, path)


def test_outer_shell_has_no_slivers(blade):
    name, mesh = blade.name, blade.mesh
    bad = []
    n_tri = 0
    for el in mesh.element_sets["allOuterShellEls"].elements:
        corners = np.asarray(el.node_coords, dtype=float)[: el.node_count]
        n_tri += len(corners) == 3
        angle, sjac = _quality(corners)
        if angle < MIN_ANGLE_DEG or sjac < MIN_SCALED_JACOBIAN:
            bad.append((angle, sjac, el.id, corners[:, 2].mean()))
    print(f"\n{name}: {n_tri} outer-shell triangles, {len(bad)} below the thresholds")
    assert not bad, (
        f"{name}: {len(bad)} outer-shell element(s) below {MIN_ANGLE_DEG} deg / sJac "
        f"{MIN_SCALED_JACOBIAN}; worst (angle, sJac, id, z): {sorted(bad)[:5]}"
    )


def _expected_web_stations(builder, element_size):
    """First and last station of the shear webs, and of the yaml domain, derived from the yaml arcs.

    Read from the blade the mesher was given (yaml web laminates, the attachment arcs of
    ``keypoints`` and the stations as inserted). A strip is a piece of the cap cut at the webs of
    the first (last) web row; it closes to a point when it is zero wide at the domain end and keeps
    widening, station after station, up to ``min_strip_ratio * element_size``. The web starts (ends)
    at the station where the last closing strip gets there; a strip that stops widening (held at a
    constant width) does not close, and with no closing strip the web keeps the yaml domain end.

    The mesher puts a station where that width is crossed (S6b: ``BladeMesh`` samples the bay
    densely and interpolates), so the station is found here as the first one at least
    ``(1 - CROSSING_REL) * min width`` wide, and the strip that sets it must be no wider than
    ``(1 + CROSSING_REL) * min width`` there: the web starts at the crossing, not a bay later.
    Also returns that strip width over the min width at each end (None with no closing strip).
    """
    ratio = inspect.signature(get_shell_mesh).parameters["min_strip_ratio"].default
    wide = ratio * element_size
    nb = builder._numad_blade
    z = np.asarray(nb.ispan, dtype=float)
    webs = range(len(nb.stackdb.swstacks))
    exists = [{g for g in webs if nb.stackdb.swstacks[g][r].plygroups} for r in range(len(z) - 1)]
    rows = [r for r, e in enumerate(exists) if e]
    assert rows == list(range(rows[0], rows[-1] + 1)), "the yaml web domain is not one run of rows"
    web_arcs = np.asarray(nb.keypoints.web_attach_arcs, dtype=float)
    plans = list(
        plan_cap_segments(
            nb.keypoints.key_arcs, web_arcs, wide, [{r for r in rows if g in exists[r]} for g in webs]
        ).values()
    )

    def closing_station(row, station, step):
        reach, binding = station, None
        for arcs, pos in plans:
            cuts = [0, *sorted(1 + pos[g] for g in exists[row] if g in pos), len(arcs) - 1]
            for lo, hi in zip(cuts[:-1], cuts[1:], strict=True):
                width = arcs[hi] - arcs[lo]
                if width[row : row + 2].max() <= 1.0e-6 or width[station] > 1.0e-6:
                    continue
                path = width[station::step]
                grows = np.flatnonzero(np.diff(path) <= 1.0e-6)  # first station that stops widening
                run = path[: (grows[0] + 1) if grows.size else path.size]
                at_least = (1.0 - CROSSING_REL) * wide
                if run.max() >= at_least:
                    k = int(np.argmax(run >= at_least))
                    at = station + step * k
                    if (at > reach) if step > 0 else (at < reach):
                        reach, binding = at, float(run[k] / wide)
                    elif at == reach and at != station:
                        # several strips reach the same station: the narrowest one set it
                        binding = min(binding, float(run[k] / wide))
        return reach, binding

    first, w_first = closing_station(rows[0], rows[0], 1)
    last, w_last = closing_station(rows[-1], rows[-1] + 1, -1)
    return float(z[first]), float(z[last]), float(z[rows[0]]), float(z[rows[-1] + 1]), (w_first, w_last)


def test_webs_and_cap_split_start_where_the_closing_strips_are_wide_enough(blade):
    name, mesh = blade.name, blade.mesh
    z_first, z_last, z_domain_lo, z_domain_hi, binding = _expected_web_stations(
        blade.builder, ELEMENT_SIZE
    )
    print(
        f"\n{name}: webs start z {z_first:.4f} (yaml domain {z_domain_lo:.3f}), "
        f"end z {z_last:.4f} (yaml domain {z_domain_hi:.3f}); binding strip / min width {binding}"
    )
    for w in binding:
        assert w is None or w <= 1.0 + CROSSING_REL, f"{name}: web end a station past the crossing ({w})"
    z_web = np.concatenate(
        [el.node_coords[:, 2] for el in mesh.element_sets["allShearWebEls"].elements]
    )
    assert z_web.min() == pytest.approx(z_first, abs=1.0e-6), f"{name}: web starts at {z_web.min():.4f}"
    assert z_web.max() == pytest.approx(z_last, abs=1.0e-6), f"{name}: web ends at {z_web.max():.4f}"
    caps = [n for n in mesh.element_sets if SPLIT_CAP.search(n)]
    z_cap = np.concatenate([el.node_coords[:, 2] for n in caps for el in mesh.element_sets[n].elements])
    cap_lo, cap_hi = float(z_cap.min()), float(z_cap.max())
    assert cap_lo >= z_first - 1.0e-6 and cap_hi <= z_last + 1.0e-6, f"{name}: cap split outside the web"
    # S8: the cap is cut at a web only where the web is not snapped onto a key line; where the split
    # stops short of the web, every web inside the cap there lies on a cap edge
    snap = ratio_snap()
    nb = blade.builder._numad_blade
    z = np.asarray(nb.ispan, dtype=float)
    key = np.asarray(nb.keypoints.key_arcs, dtype=float)
    raw = np.asarray(nb.keypoints.web_attach_arcs, dtype=float)
    _, onto = _snapped_arcs(key, raw, snap)
    for side, (lo_row, hi_row) in enumerate(((3, 4), (8, 9))):
        inside = (raw[:, side, :] > key[lo_row] + 1.0e-9) & (raw[:, side, :] < key[hi_row] - 1.0e-9)
        for k in np.flatnonzero(((z > z_first + 1.0e-6) & (z < cap_lo - 1.0e-6)) | ((z > cap_hi + 1.0e-6) & (z < z_last - 1.0e-6))):
            bad = inside[:, k] & ~onto[:, side, k]
            assert not bad.any(), f"{name}: web {np.flatnonzero(bad).tolist()} is cut-free at z {z[k]:.3f} without lying on a cap edge"


def ratio_snap():
    return inspect.signature(get_shell_mesh).parameters["min_strip_ratio"].default * ELEMENT_SIZE


# A key point replaces the nearest polyline vertex but keeps that vertex's arc, so the arc of a web
# edge node is known to one polyline segment (measured IEA-15: 0.011-0.026 m; segments are longer).
# The nodes lie on the polyline to 1e-4 m (the replaced vertex), hence ``ON_POLYLINE``.
ON_POLYLINE = 1.0e-3


def _edge_node_arcs(mesh, geometry, k, z_k):
    """Signed arcs (from the leading edge) of the web nodes lying on the section polyline at station k."""
    ring, arclen = geometry.coordinates[:, :, k], geometry.arclength[:, k]
    found = []
    for name, s in mesh.element_sets.items():
        if not (name.endswith(f"_{k}_SW") or name.endswith(f"_{k - 1}_SW")):
            continue
        for e in s.elements:
            for p in e.node_coords[np.abs(e.node_coords[:, 2] - z_k) < 1.0e-4]:
                for i0, i1 in zip(range(len(ring) - 1), range(1, len(ring)), strict=True):
                    a, c = ring[i0], ring[i1]
                    t = float(np.clip((p - a) @ (c - a) / ((c - a) @ (c - a)), 0.0, 1.0))
                    if np.linalg.norm(a + t * (c - a) - p) < ON_POLYLINE:
                        found.append(arclen[i0] + t * (arclen[i1] - arclen[i0]))
                        break
    return np.unique(np.round(found, 4))


def test_no_web_is_displaced_beyond_the_snap_distance(blade):
    """S8: the web attachment lines sit on the yaml arc, at the stations, to within the snap.

    Measured on the mesh itself: the arc, along the section polyline, of the web nodes that lie on
    it, against the arcs ``start_nd_arc`` / ``end_nd_arc`` give, read by the audit's independent
    yaml reader (``YamlBlade.web_state``) and turned into metres with the polyline's own
    perimeter between the trailing-edge corners. (``audit_station`` is not used: its suction-side
    residual carries a ~0.1 m bias at inboard stations that is not in the mesh, L34; the clip onto the cap
    edge moved IEA-3.4 webs by up to 0.86 m.)
    """
    snap = ratio_snap()
    nb = blade.builder._numad_blade
    z = np.asarray(nb.ispan, dtype=float)
    g = nb.geometry
    yaml_blade = ab.YamlBlade.from_file(str(blade.path))
    nweb = len(yaml_blade.struct["webs"])
    rows = [r for r in range(len(z) - 1) if len(yaml_blade.web_layers(0, (z[r] - z[0]) / (z[-1] - z[0])))]
    stations = [k for k in range(1, len(z) - 1) if rows and rows[0] < k < rows[-1]][::6]
    assert len(stations) > 10, f"{blade.name}: only {len(stations)} stations to measure"
    worst = (-np.inf, 0.0, 0.0)
    for k in stations:
        eta = (z[k] - z[0]) / (z[-1] - z[0])
        arc0, perimeter = g.arclength[1, k], g.arclength[g.coordinates.shape[0] - 2, k] - g.arclength[1, k]
        have = _edge_node_arcs(blade.mesh, g, k, z[k])
        resolution = float(np.linalg.norm(np.diff(g.coordinates[:, :, k], axis=0), axis=1).max())
        for w in range(nweb):
            if not yaml_blade.web_layers(w, eta):
                continue
            ws = yaml_blade.web_state(w, eta)
            for want in (arc0 + (1.0 - ws["start"]) * perimeter, arc0 + (1.0 - ws["end"]) * perimeter):
                gap = float(np.min(np.abs(have - want))) if have.size else np.inf
                worst = max(worst, (gap - resolution, eta, gap))
    assert worst[0] <= snap, (
        f"{blade.name}: a web is {worst[2]:.3f} m off its yaml arc at span fraction {worst[1]:.3f}; "
        f"the snap distance is {snap:.3f} m (plus one polyline segment of measurement resolution)"
    )


def _web_groups_z(mesh):
    """{web group: sorted z intervals of its element sets}."""
    groups = {}
    for name, s in mesh.element_sets.items():
        if name.endswith("_SW"):
            z = np.concatenate([el.node_coords[:, 2] for el in s.elements])
            groups.setdefault(name.split("_")[0], []).append((float(z.min()), float(z.max())))
    return {g: sorted(iv) for g, iv in groups.items()}


def web_gaps(mesh):
    """{web group: [(z end, z next start)]} for the holes inside each web's span extent."""
    gaps = {}
    for group, intervals in _web_groups_z(mesh).items():
        reach = intervals[0][1]
        for lo, hi in intervals[1:]:
            if lo - reach > 1.0e-6:
                gaps.setdefault(group, []).append((round(reach, 3), round(lo, 3)))
            reach = max(reach, hi)
    return gaps


def test_no_web_gap_inside_the_yaml_domain_except_the_s6_ends(blade):
    """Each web is one run of rows between the S6 stations: no hole inside it, none off the domain."""
    name, mesh = blade.name, blade.mesh
    assert not web_gaps(mesh), f"{name}: holes inside a web {web_gaps(mesh)}"
    z_lo, z_hi = blade.web_domain
    for group, intervals in _web_groups_z(mesh).items():
        assert intervals[0][0] >= z_lo - 1.0e-6 and max(i[1] for i in intervals) <= z_hi + 1.0e-6, (
            f"{name}: web {group} outside the yaml web domain [{z_lo:.2f}, {z_hi:.2f}]"
        )


def _synthetic_keys(n_stations=6):
    """13 key lines 1 m apart (arc rows 0..12), the same at every station."""
    return np.repeat(np.arange(13, dtype=float)[:, None], n_stations, axis=1)


def test_a_crossing_that_cannot_snap_raises():
    """S8: a web that is further than the snap distance from the key line it crosses, at either
    station of the row, is neither clipped nor snapped: the mesher stops."""
    keys = _synthetic_keys()
    arc = np.full((1, 2, 6), 8.5)  # LP side mid segment 8 ...
    arc[0, 0, :] = [2.5, 2.5, 2.5, 3.4, 3.4, 3.4]  # ... HP side across line 3: 0.5 m before, 0.4 m after
    with pytest.raises(ValueError, match="crosses"):
        snap_web_arcs(keys, arc, snap=0.075, web_rows=[set(range(5))])


def test_a_crossing_within_the_snap_distance_lands_on_the_key_line():
    keys = _synthetic_keys()
    arc = np.full((1, 2, 6), 8.5)
    arc[0, 0, :] = [2.5, 2.5, 2.97, 3.05, 3.5, 3.5]  # crosses line 3, 3 cm and 5 cm from it
    snapped = snap_web_arcs(keys, arc, snap=0.075, web_rows=[set(range(5))])
    np.testing.assert_array_equal(snapped[0, 0], [2.5, 2.5, 3.0, 3.0, 3.5, 3.5])
    np.testing.assert_array_equal(snapped[0, 1], arc[0, 1])  # a web far from every line is untouched
    # a web that never comes within the snap distance of a line is the yaml arc, to the bit
    assert np.all(np.abs(snapped - arc)[0, 0, [0, 1, 4, 5]] == 0.0)


def test_a_web_outside_the_cap_is_cut_in_the_segment_that_holds_it():
    keys = _synthetic_keys()
    arc = np.full((1, 2, 6), 2.5)  # HP: segment 2 (between lines 2 and 3), not the cap
    arc[0, 1, :] = 8.5  # LP: segment 8
    plans = plan_cap_segments(keys, arc, snap=0.075, web_rows=[set(range(5))])
    assert set(plans) == {2, 8}
    arcs, pos = plans[2]
    np.testing.assert_array_equal(arcs[1], np.full(6, 2.5))  # the web line is the yaml arc, no margin
    assert pos == {0: 0}
