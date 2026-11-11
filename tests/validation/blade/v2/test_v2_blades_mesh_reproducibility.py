"""The gmsh mesher on every windIO v2 blade of the repository: quality gate and material fidelity (R1, R2 of
odd/tasks/gmsh-store-campaign.md, L6), tolerances fixed before the first number of these tests.

* R1 gate: AR <= 10, minimum corner angle >= 5 deg, no non-positive Jacobian, no coincident or unused nodes, no element
  with a repeated node (the mesher's own gate limits, evaluated here on the cached arrays, skin and webs).
* R2 material fidelity vs the YAML laminate: assigned mass of every material carrying >= 1 % of the blade mass within
  0.5 %, total skin mass within 0.5 %, at most 1 % of the skin elements with > 5 % areal-mass deviation.

Sizes are the certified ones of gmsh-validation-15-22 L22 (IEA-15 0.25, IEA-22 0.177); IEA-3.4 uses 0.1 m, the best
available, which is NOT certified (one outlier at 0.07 of unknown cause), and is marked so in its test id.
Marks: slow (mesh generation, minutes; the raw arrays are cached under $SCRATCH/tmp/v2_blade_cache).
"""

from __future__ import annotations

import numpy as np
import pytest

pytest.importorskip("gmsh", reason="gmsh not available")

from tests.support.paths import REPO_ROOT  # noqa: E402
from tests.support.v2_blade import load_raw  # noqa: E402
from tests.support.v2_fidelity import skin_material_balance  # noqa: E402

BLADES = [
    pytest.param("IEA-15", "tests/IEA15MW/reference/upstream/IEA-15-240-RWT.yaml", 0.25, id="IEA-15-h0.25"),
    pytest.param("IEA-22", "tests/IEA22MW/reference/upstream/IEA-22-280-RWT.yaml", 0.177, id="IEA-22-h0.177"),
    pytest.param("IEA-3.4", "tests/IEA34/reference/upstream/IEA-3.4-130-RWT.yaml", 0.1, id="IEA-3.4-h0.1-not-certified"),
]
ES_RANGE = (0.05, 0.25)
MASS_TOL = 0.005
FRAC_GT5_MAX = 0.01


@pytest.fixture(scope="module")
def raws():
    return {}


def _raw(raws, name, yaml_rel, h):
    if name not in raws:
        raws[name] = load_raw(h, yaml_file=REPO_ROOT / yaml_rel, mesher="gmsh", es_range=ES_RANGE)
    return raws[name]


@pytest.mark.slow
@pytest.mark.parametrize(("name", "yaml_rel", "h"), BLADES)
def test_mesh_passes_the_quality_gate(raws, name, yaml_rel, h):
    from aeroelast.models.blade.numad.mesh_gen.gmsh_skin import DEFAULT_LIMITS, element_metrics

    raw = _raw(raws, name, yaml_rel, h)
    xyz, el = np.asarray(raw["nodes"], float), np.asarray(raw["elements"])
    m = element_metrics(xyz, el)
    assert m["aspect"].max() <= DEFAULT_LIMITS["max_ar"], f"{name}: AR max {m['aspect'].max():.2f}"
    assert m["min_angle"].min() >= DEFAULT_LIMITS["min_angle"], f"{name}: min angle {m['min_angle'].min():.2f} deg"
    used = np.unique(el[el >= 0])
    assert len(used) == len(xyz), f"{name}: {len(xyz) - len(used)} unused nodes"
    assert all(len(set(r[r >= 0])) == len(r[r >= 0]) for r in el), f"{name}: element with a repeated node"
    rounded = np.round(xyz[used] / 1e-6).astype(np.int64)
    _, inv, cnt = np.unique(rounded, axis=0, return_inverse=True, return_counts=True)
    extra = int((cnt - 1).sum())
    if extra:
        surf = np.asarray(raw["surfaces"])
        owners = {}
        for k, e in enumerate(el):
            for i in e[e >= 0]:
                owners.setdefault(int(i), set()).add(str(surf[k]))
        pairs = sorted({tuple(sorted(set().union(*[owners[int(used[j])] for j in np.nonzero(np.ravel(inv) == g)[0]])))
                        for g in np.nonzero(cnt > 1)[0]})
        z_over = xyz[used][np.isin(np.ravel(inv), np.nonzero(cnt > 1)[0]), 2]
        raise AssertionError(f"{name}: {extra} extra nodes at coincident positions between surfaces {pairs}, "
                             f"z {z_over.min():.3f}-{z_over.max():.3f} m (webs sharing the same attachment line)")


@pytest.mark.slow
@pytest.mark.parametrize(("name", "yaml_rel", "h"), BLADES)
def test_material_fidelity_against_the_yaml(raws, name, yaml_rel, h):
    raw = _raw(raws, name, yaml_rel, h)
    bal = skin_material_balance(raw, REPO_ROOT / yaml_rel)
    total = bal["mass_src_kg"]
    print(f"{name} h={h}: skin mass src {total:.1f} asg {bal['mass_asg_kg']:.1f} kg "
          f"({(bal['mass_asg_kg'] / total - 1) * 100:+.3f} %), elements >5 % areal-mass deviation {bal['frac_gt5'] * 100:.2f} %")
    assert abs(bal["mass_asg_kg"] / total - 1) <= MASS_TOL
    for mat, (src, asg) in bal["per_material"].items():
        if src >= 0.01 * total:
            print(f"   {mat:24s} src {src:10.1f} asg {asg:10.1f} {(asg / src - 1) * 100:+.3f} %")
            assert abs(asg / src - 1) <= MASS_TOL, f"{name}: {mat} {(asg / src - 1) * 100:+.3f} %"
    assert bal["frac_gt5"] <= FRAC_GT5_MAX
