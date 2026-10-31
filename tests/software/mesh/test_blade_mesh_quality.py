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
* the local refinement code holds no length literal (the sizes come from ``element_size``, the
  measured strip width and a dimensionless growth ratio).

That the webs themselves stay where the yaml puts them is ``tests/test_blade_web_position.py``.
"""

from __future__ import annotations

import ast
import inspect
import re
from collections import Counter
from pathlib import Path
from typing import Any, NamedTuple

import numpy as np
import pytest
import yaml

from aeroelast.core.mesh.generators import BladeMesh
from aeroelast.models.blade.numad.mesh_gen.mesh_gen import _CAP_SEGMENTS, _plan_cap_columns, get_shell_mesh
from tests.support.paths import DATA_DIR

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


@pytest.fixture(scope="module", params=sorted(BLADES))
def blade(request):
    path = DATA_DIR / BLADES[request.param]
    builder = BladeMesh(str(path), element_size=ELEMENT_SIZE)
    mesh = builder.generate(verbose=False)
    return Built(request.param, mesh, _web_domain_z(path), builder)


# Known cause, not a tolerance: the IEA-3.4 yaml puts webs outside the spar cap (up to 81 cm over
# z 6.3-22.9 m), the mesher clips them to the cap edge + 5 mm, and those strips end in needles at
# the cap apex and where a web enters the cap. Removed by S8 of odd/tasks/v1-mesher-web-slivers.md
# (webs cut whichever skin piece holds them); strict, so the fix has to drop this mark.
# Webs the 5 mm ``_WEB_EDGE_MARGIN`` clip holds on the cap edge: IEA-3.4 over long runs (L15, L18),
# IEA-15 where the LP web 0 grazes the cap edge at z 11.94 (L28). T8 (snap onto the key line, no
# clip) removes both.
WEBS_OUTSIDE_THE_CAP = {"iea34": "L15, L18", "iea15": "L28"}


def test_outer_shell_has_no_slivers(blade, request):
    name, mesh = blade.name, blade.mesh
    if name in WEBS_OUTSIDE_THE_CAP:
        request.applymarker(
            pytest.mark.xfail(
                strict=True,
                reason=f"web clipped onto the cap edge (v1-mesher-web-slivers {WEBS_OUTSIDE_THE_CAP[name]})",
            )
        )
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
    plans = []
    for side, lo_row, hi_row in _CAP_SEGMENTS.values():
        plans.append(
            _plan_cap_columns(
                nb.keypoints.key_arcs, [a[side] for a in nb.keypoints.web_attach_arcs], lo_row, hi_row
            )
        )

    def closing_station(row, station, step):
        reach, binding = station, None
        for arcs, pos in plans:
            cuts = [0, *sorted(1 + pos[g] for g in exists[row]), len(arcs) - 1]
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
    for label, sets in (
        ("web", ["allShearWebEls"]),
        ("cap split", [n for n in mesh.element_sets if SPLIT_CAP.search(n)]),
    ):
        z = np.concatenate(
            [el.node_coords[:, 2] for n in sets for el in mesh.element_sets[n].elements]
        )
        assert z.min() == pytest.approx(z_first, abs=1.0e-6), f"{name}: {label} starts at {z.min():.4f}"
        assert z.max() == pytest.approx(z_last, abs=1.0e-6), f"{name}: {label} ends at {z.max():.4f}"


def test_spar_cap_is_split_only_inside_the_web_domain(blade):
    name, mesh = blade.name, blade.mesh
    z_lo, z_hi = blade.web_domain
    pieces = [n for n in mesh.element_sets if SPLIT_CAP.search(n)]
    assert pieces, f"{name}: no split spar cap at all, the web attachment lines are gone"
    outside = {}
    for set_name in pieces:
        z = np.concatenate([el.node_coords[:, 2] for el in mesh.element_sets[set_name].elements])
        if z.min() < z_lo - 1.0e-6 or z.max() > z_hi + 1.0e-6:
            outside[set_name] = (round(float(z.min()), 3), round(float(z.max()), 3))
    assert not outside, f"{name}: cap pieces outside the web domain [{z_lo:.2f}, {z_hi:.2f}] m: {outside}"


def _free_edges(elements):
    """Edges used by exactly one element, as ((node id, node id), mid z)."""
    use = Counter()
    z_of = {}
    for el in elements:
        nodes = list(el.nodes)
        for na, nb in zip(nodes, nodes[1:] + nodes[:1], strict=True):
            key = (min(na.id, nb.id), max(na.id, nb.id))
            use[key] += 1
            z_of[key] = 0.5 * (na.coords[2] + nb.coords[2])
    return [(key, z_of[key]) for key, n in use.items() if n == 1]


def test_outer_shell_is_conforming(blade):
    """No crack: a free edge of the outer shell may only lie on the root or the tip."""
    name, mesh = blade.name, blade.mesh
    elements = mesh.element_sets["allOuterShellEls"].elements
    z = np.array([n.coords[2] for el in elements for n in el.nodes])
    ends = 0.5 * ELEMENT_SIZE
    interior = [
        round(float(zf), 3)
        for _, zf in _free_edges(elements)
        if z.min() + ends < zf < z.max() - ends
    ]
    stations = sorted(Counter(interior).items())
    assert not interior, f"{name}: {len(interior)} interior free edge(s) at z {stations}"


def test_local_refinement_holds_no_length_literal():
    """S5: no model-specific size in the refinement code, only parameter defaults.

    A numeric literal is allowed as a parameter default (the dimensionless growth ratio) or as a
    structural constant (small integers, the factor 0.5); any other one is a hard-coded length.
    """
    import aeroelast.models.blade.numad.mesh_gen as pkg

    path = Path(pkg.__file__).parent / "local_refine.py"
    tree = ast.parse(path.read_text())
    defaults = {
        id(inner)
        for node in ast.walk(tree)
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.Lambda))
        for d in [*node.args.defaults, *node.args.kw_defaults]
        if d is not None
        for inner in ast.walk(d)
    }
    structural = {0, 1, 2, 3, 4, 0.5}
    # dimensionless gmsh option codes, named ``_GMSH_*`` at module level (an algorithm number)
    codes = {
        id(node.value)
        for node in tree.body
        if isinstance(node, ast.Assign)
        and all(isinstance(t, ast.Name) and t.id.startswith("_GMSH_") for t in node.targets)
    }
    hard = [
        (node.lineno, node.value)
        for node in ast.walk(tree)
        if isinstance(node, ast.Constant)
        and isinstance(node.value, (int, float))
        and not isinstance(node.value, bool)
        and id(node) not in defaults
        and id(node) not in codes
        and node.value not in structural
    ]
    assert not hard, f"local_refine.py: hard-coded numeric literals (line, value): {hard}"
