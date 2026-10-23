"""The outermost yaml skin layer sits on the outward face of every shell element.

Candidate F2 of the blade input audit (``odd/tasks/blade-input-audit.md``, T12). windIO lists
``layers`` from the outside in (``UV_protection``/gelcoat first, ``Shell_skin_inner`` last) and NuMAD
adds plies in that order. The Rust laminate puts ply 0 at ``z = -h/2``, i.e. on the ``-normal`` side,
and the skin element normals point outward, so the gelcoat used to sit on the inner face. The stack
is now emitted inner-to-outer so that the last ply (the outermost yaml layer) lies on ``+normal``.
"""

import sys
from pathlib import Path

import numpy as np
import pytest
from matplotlib.path import Path as MplPath

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))
import audit_blade_input as audit  # noqa: E402

YAML = Path(__file__).parent / "IEA-15-240-RWT.yaml"
STATIONS = (0.30, 0.60, 0.90, 0.97)


@pytest.fixture(scope="module")
def mesh():
    m = audit.build_mesh(str(YAML), 0.25)
    return m


def _outward(segs, pts):
    poly = MplPath(np.vstack([pts, pts[:1]]))
    out = []
    for s in segs:
        mid = 0.5 * (np.asarray(s["pts"][0]) + np.asarray(s["pts"][1]))
        n = s["normal"][:2] / (np.linalg.norm(s["normal"][:2]) + 1e-300)
        out.append(not poly.contains_point(mid + 1e-3 * n))
    return np.array(out)


@pytest.mark.parametrize("frac", STATIONS)
def test_outermost_layer_is_on_the_outward_face(mesh, frac):
    L = mesh["nodes"][:, 2].max()
    names = [n for n in mesh["sets"] if not n.endswith("_SW") and not n.startswith("all")]
    segs = audit._slice(mesh, frac * L, names)
    pts, _ = audit._chain([s for s in segs if not s["set"].endswith("TE_FLAT")])
    outward = _outward(segs, pts)
    good = checked = 0
    for s, is_out in zip(segs, outward, strict=True):
        plies = mesh["stacks"][s["set"]]
        gel = [p for p in plies if p["mat"] == "Gelcoat"]
        if len(plies) < 2 or not gel:
            continue
        checked += 1
        # gelcoat is the outermost layer: it must be the ply farthest along the outward direction
        z_out = max(p["z_top"] for p in plies) if is_out else min(p["z_bottom"] for p in plies)
        edge = gel[0]["z_top"] if is_out else gel[0]["z_bottom"]
        good += abs(edge - z_out) < 1e-12
    assert checked >= 15
    assert good / checked >= 0.99, f"{good}/{checked} ring elements carry the gelcoat on the outer face"


def test_innermost_yaml_layer_is_on_the_inner_face(mesh):
    """``Shell_skin_inner`` is the last yaml layer: it must be the first (bottom) ply where present."""
    first = [pl[0]["mat"] for name, pl in mesh["stacks"].items() if not name.endswith("_SW") and len(pl) > 1]
    last = [pl[-1]["mat"] for name, pl in mesh["stacks"].items() if not name.endswith("_SW") and len(pl) > 1]
    assert last.count("Gelcoat") / len(last) > 0.95
    assert first.count("Gelcoat") == 0
