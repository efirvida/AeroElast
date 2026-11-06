import glob
import os

import pytest

from aeroelast.core.mesh.entities import MeshElement, Node
from aeroelast.models.blade.model import Blade
from tests.support.paths import DATA_DIR  # noqa: E402

# Reference turbine YAMLs. The historical
# `examples/reference_turbines/yamls` directory no longer exists in this tree, and
# scanning only it left `yaml_files` EMPTY -- pytest then reported "got empty
# parameter set" and this test never ran at all. The repository's reference turbine
# definition sits next to this file (it is the one tests/test_bem_engine.py uses).
_HERE = os.path.dirname(os.path.abspath(__file__))
_SEARCH_DIRS = [
    str(DATA_DIR),
    str(DATA_DIR / "reference_turbines" / "yamls"),
]
yaml_files = sorted(
    f
    for d in _SEARCH_DIRS
    if os.path.isdir(d)
    for f in glob.glob(os.path.join(d, "*.yaml"))
)

# Fail loudly rather than vanish: an empty parametrisation is a silently disabled
# test, which is worse than a red one.
if not yaml_files:
    raise RuntimeError("no reference turbine YAML found in " + ", ".join(_SEARCH_DIRS))

# Extract just the filenames for test IDs
yaml_file_ids = [os.path.basename(f) for f in yaml_files]


# Measured at 3185da6 (logs/check/blade-mesh-pins-20261027.log). The single 9277-node pin
# went stale with the blade-input-audit mesh fixes (2026-10-23..25, deliberate: TE band to
# the yaml width, cap arcs) and then matched no yaml; IEA-3.4 arrived with 188c997.
# Re-pinned for odd/tasks/v1-mesher-web-slivers.md T4c (S3, S6b: a station where each web starts and
# ends, the same for the native and the default span-plan path; logs/check/t4c-pins.log). Old -> new
# (nodes, elements, allOuterShellNods, allShearWebNods): IEA-15 and IEA-15_rigid 10637, 11248, 9446,
# 2413 -> 10612, 11245, 9429, 2403; IEA-3.4 4453, 4795, 4081, 1002 -> 3898, 4019, 3653, 428. The
# IEA-3.4 web loss is not a baseline: at element_size 0.5 the S6 walk drops web rows over z 20-59
# where its webs are clipped onto the cap edge (L25); T8 reworks that and re-pins.
PINNED_COUNTS = {
    "IEA-15-240-RWT.yaml": (
        10612, 11245, {"RootNodes": 37, "allOuterShellNods": 9429, "allShearWebNods": 2403}),
    "IEA-15-240-RWT_rigid.yaml": (
        10612, 11245, {"RootNodes": 37, "allOuterShellNods": 9429, "allShearWebNods": 2403}),
    "IEA-3.4-130-RWT.yaml": (
        3898, 4019, {"RootNodes": 20, "allOuterShellNods": 3653, "allShearWebNods": 428}),
}


@pytest.fixture(autouse=True)
def reset_counters():
    """Fixture to reset ID counters before each test"""
    Node._id_counter = 0
    MeshElement._id_counter = 0


@pytest.mark.parametrize("blade_file", yaml_files, ids=yaml_file_ids)
def test_blade_mesh_generation(blade_file):
    """Test blade mesh generation for all reference turbine YAML files"""
    blade = Blade(blade_file, element_size=0.5)
    blade.generate_mesh()

    # Explicit, so an unset mesh fails HERE with a clear message instead of an
    # AttributeError three lines down, and so a type checker can see it.
    mesh = blade.mesh
    assert mesh is not None, f"{blade_file}: generate_mesh() left blade.mesh unset"

    # Exact counts per blade at element_size=0.5, so a silent mesh change fails here instead
    # of passing on `> 0`. A yaml without an entry fails: pin it when it is added.
    name = os.path.basename(blade_file)
    assert name in PINNED_COUNTS, f"{name}: no pinned mesh counts; add them to PINNED_COUNTS"
    nodes, elements, sets = PINNED_COUNTS[name]
    assert mesh.node_count == nodes, f"{blade_file}: {mesh.node_count} nodes, expected {nodes}"
    assert mesh.elements_count == elements, (
        f"{blade_file}: {mesh.elements_count} elements, expected {elements}"
    )
    for name, expected in sets.items():
        assert name in mesh.node_sets, f"{blade_file} missing '{name}' node set"
        assert len(mesh.node_sets[name].node_ids) == expected, (
            f"{blade_file}: '{name}' has {len(mesh.node_sets[name].node_ids)} nodes, expected {expected}"
        )
