"""The shear webs attach to the skin where the yaml puts them (blade input audit, T8, candidate E).

NuMAD hung both webs on the spar-cap edges (keypoints ``b`` and ``c``) whatever the yaml web
``start_nd_arc`` / ``end_nd_arc`` said: 0.10-0.16 m off at r/R 0.30-0.90. The spar-cap region is now
split at the web attachment arcs into ``HP_SPAR_A/B/C`` (``LP_SPAR_A/B/C``) so the web edge shares
its nodes with the skin. The split is geometric only: every sub-region keeps the laminate of the cap.

The mesh is the production one (es 0.25 m, default bays, about 60 s), built once for the module.
The yaml side is read by the audit's independent reader (``tools/audit_blade_input.py``).
"""

from __future__ import annotations

import re
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))

import audit_blade_input as ab  # noqa: E402

YAML = str(Path(__file__).parent / "IEA-15-240-RWT.yaml")
STATIONS = (0.30, 0.60, 0.90)
TOL_ARC_M = 0.01
CAP = re.compile(r"^(\d\d)_(\d+)_((?:HP|LP)_SPAR)(?:_([ABC]))?$")


@pytest.fixture(scope="module")
def built():
    from aeroelast.models.blade.model import Blade

    blade = Blade(YAML, element_size=0.25)
    blade.generate_mesh()
    return blade, ab.add_set_zrange(ab.build_mesh(YAML, 0.25, blade=blade))


@pytest.fixture(scope="module")
def yaml_blade():
    return ab.YamlBlade.from_file(YAML)


@pytest.mark.parametrize("frac", STATIONS)
def test_web_attachment_nodes_sit_on_the_yaml_nd_arc(built, yaml_blade, frac):
    _, mesh = built
    checks = ab.audit_station(yaml_blade, mesh, frac)["checks"]
    got = [c for c in checks if c["field"] == "F4" and c["quantity"].startswith("web attachment")]
    assert len(got) == 4, "two webs, two sides"
    for c in got:
        assert abs(c["residual"]) <= TOL_ARC_M, (frac, c["quantity"], c["location"], c["residual"])


def _skin_node_ids(mesh):
    skin = {n for n in mesh["sets"] if not n.endswith("_SW") and not n.startswith("all")}
    return set(np.concatenate([mesh["elements"][e] for n in skin for e in mesh["sets"][n]]).tolist())


def test_web_edges_are_shared_with_the_skin(built):
    _, mesh = built
    els = mesh["elements"]
    skin_nodes = _skin_node_ids(mesh)
    for name in sorted(n for n in mesh["sets"] if n.endswith("_SW")):
        group = name.split("_")[0]
        edges = Counter()
        for e in mesh["sets"][name]:
            nd = els[e]
            for k in range(len(nd)):
                edges[frozenset((nd[k], nd[(k + 1) % len(nd)]))] += 1
        xyz = mesh["nodes"]
        free = [tuple(ed) for ed, n in edges.items() if n == 1]
        # a free edge across the chord is a bay end (the sections are not planar, so it is not
        # exactly at constant z); the ones that run mostly along z are the attachment lines
        along = [
            ed for ed in free
            if abs(xyz[ed[0], 2] - xyz[ed[1], 2]) > 0.5 * np.linalg.norm(xyz[ed[0]] - xyz[ed[1]])
        ]
        assert along, f"web {name} has no edge along the span"
        bad = [ed for ed in along if not (ed[0] in skin_nodes and ed[1] in skin_nodes)]
        assert not bad, f"web group {group}: {len(bad)} edge(s) along the span are free of the skin"


def test_cap_is_split_and_keeps_the_laminate_of_the_parent_region(built):
    blade, mesh = built
    stacks = blade._numad_blade.stackdb.stacks
    layups = {s["elementSet"]: s["layup"] for s in blade._numad_mesh["sections"]}
    n_cap = 0
    split = set()
    for name in mesh["sets"]:
        m = CAP.match(name)
        if not m:
            continue
        seg, st, label, part = int(m[1]), int(m[2]), m[3], m[4]
        parent = stacks[seg, st]
        assert parent.name.endswith(label)
        want = [[pg.materialid, 0.001 * pg.thickness * pg.nPlies, pg.angle] for pg in parent.plygroups]
        got = layups[name]
        assert [g[0] for g in got] == [w[0] for w in want], name
        np.testing.assert_allclose([g[1] for g in got], [w[1] for w in want], rtol=1e-12, atol=0, err_msg=name)
        n_cap += 1
        if part:
            split.add((seg, st, part))
    assert n_cap > 0
    for (seg, st) in {(a, b) for a, b, _ in split}:
        assert {p for a, b, p in split if (a, b) == (seg, st)} >= {"A", "C"}, (seg, st)
    assert any(p == "B" for _, _, p in split), "the middle sub-region exists where the two webs differ"
