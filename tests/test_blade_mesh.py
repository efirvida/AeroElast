import glob
import os

import pytest

from aeroelast.core.mesh.entities import MeshElement, Node
from aeroelast.models.blade.model import Blade

# Reference turbine YAMLs. The historical
# `examples/reference_turbines/yamls` directory no longer exists in this tree, and
# scanning only it left `yaml_files` EMPTY -- pytest then reported "got empty
# parameter set" and this test never ran at all. The repository's reference turbine
# definition sits next to this file (it is the one tests/test_bem_engine.py uses).
_HERE = os.path.dirname(os.path.abspath(__file__))
_SEARCH_DIRS = [
    _HERE,
    os.path.join(_HERE, "..", "examples", "reference_turbines", "yamls"),
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

    # Exact counts for the reference blade at element_size=0.5, so a silent mesh
    # change fails here instead of passing on `> 0`.
    assert mesh.node_count == 9277, f"{blade_file}: {mesh.node_count} nodes, expected 9277"
    assert mesh.elements_count == 9867, (
        f"{blade_file}: {mesh.elements_count} elements, expected 9867"
    )
    for name, expected in (
        ("RootNodes", 37),
        ("allOuterShellNods", 8204),
        ("allShearWebNods", 2277),
    ):
        assert name in mesh.node_sets, f"{blade_file} missing '{name}' node set"
        assert len(mesh.node_sets[name].node_ids) == expected, (
            f"{blade_file}: '{name}' has {len(mesh.node_sets[name].node_ids)} nodes, expected {expected}"
        )
