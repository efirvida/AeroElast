"""Shared pytest fixtures for the fem-shell test suite."""

import ctypes
import os

import pytest

# Directories under tests/ that are NOT pytest test packages (papers,
# external repo clones, data collections).  Without this, pytest recurses
# into them and errors during collection.
#
# NOTE: test_vol_mesh.py is stale against the removed cap_mesh module (same
# category as test_blade_mesh.py) — exclude it on the CLI:
#   --ignore=tests/test_vol_mesh.py
collect_ignore = [
    "IEA15MW/validation papers",
    "IEA15MW/75698.pdf.txt",
]


def _prepend_ld_library_path(path: str) -> None:
    current = os.environ.get("LD_LIBRARY_PATH", "")
    entries = [entry for entry in current.split(":") if entry]
    if path not in entries:
        os.environ["LD_LIBRARY_PATH"] = f"{path}:{current}" if current else path


def _ensure_shared_lib(lib_name: str, lib_dir: str) -> None:
    try:
        ctypes.CDLL(lib_name, mode=ctypes.RTLD_GLOBAL)
        return
    except OSError:
        _prepend_ld_library_path(lib_dir)

    try:
        ctypes.CDLL(os.path.join(lib_dir, lib_name), mode=ctypes.RTLD_GLOBAL)
    except OSError:
        pass


# ---------------------------------------------------------------------------
# GCC runtime — required by preCICE / _aeroelast on this cluster
# ---------------------------------------------------------------------------
_GCC14_LIB_PATH = "/petrobr/app_sequana/gcc/14.2.0/lib64"
_ensure_shared_lib("libstdc++.so.6", _GCC14_LIB_PATH)

# ---------------------------------------------------------------------------
# GLU library — required by gmsh (equivalent to `module load glu`)
# ---------------------------------------------------------------------------
_GLU_LIB_PATH = "/scratch/app/glu/9.0.2_gnu/lib"
_ensure_shared_lib("libGLU.so.1", _GLU_LIB_PATH)

# ---------------------------------------------------------------------------
# Blade YAML fixture
# ---------------------------------------------------------------------------
# Candidate paths for the IEA-15-240-RWT turbine definition YAML.
# The file is large so it lives outside the repo in the simulations tree,
# but a copy may also exist inside examples/.
_BLADE_YAML_CANDIDATES = [
    # Co-located with tests (highest priority)
    os.path.join(
        os.path.dirname(__file__),
        "IEA-15-240-RWT.yaml",
    ),
    # Inside repo (preferred when present)
    os.path.join(
        os.path.dirname(__file__),
        "..",
        "examples",
        "reference_turbines",
        "yamls",
        "IEA-15-240-RWT.yaml",
    ),
    # Sibling simulations directory (dev environment)
    os.path.join(
        os.path.dirname(__file__),
        "..",
        "..",
        "simulations",
        "blade",
        "solid",
        "IEA-15-240-RWT.yaml",
    ),
    os.path.join(
        os.path.dirname(__file__),
        "..",
        "..",
        "simulations",
        "IEA-15-240-RWT",
        "fsi",
        "NoTower",
        "solid",
        "IEA-15-240-RWT.yaml",
    ),
]


def _find_blade_yaml():
    for p in _BLADE_YAML_CANDIDATES:
        resolved = os.path.abspath(p)
        if os.path.isfile(resolved):
            return resolved
    return None


@pytest.fixture(scope="session")
def iea_blade_yaml():
    """Return the path to IEA-15-240-RWT.yaml or skip the test."""
    path = _find_blade_yaml()
    if path is None:
        pytest.skip(
            "IEA-15-240-RWT.yaml not found. "
            "Place a copy at examples/reference_turbines/yamls/ or alongside "
            "the simulations/blade/solid/ tree."
        )
    return path


@pytest.fixture(scope="session")
def iea_numad_blade(iea_blade_yaml):
    """Return a fully initialised numadBlade for the IEA-15-240-RWT."""
    import numpy as np

    from aeroelast.models.blade.numad import Blade as numadBlade

    # Same setup sequence used by BladeMesh.generate()
    blade = numadBlade()
    blade.read_yaml(iea_blade_yaml)
    for stat in blade.definition.stations:
        stat.airfoil.resample(n_samples=150)
    blade.update_blade()
    n_stations = blade.geometry.coordinates.shape[2]
    blade.expand_blade_geometry_te(0.001 * np.ones(n_stations))
    return blade


@pytest.fixture(scope="session")
def iea_blade_xlsx():
    """Return the Path to NuMAD_utd_iea15mw.xlsx or skip the test."""
    from pathlib import Path

    path = Path(__file__).parent / "NuMAD_utd_iea15mw.xlsx"
    if not path.is_file():
        pytest.skip(
            "NuMAD_utd_iea15mw.xlsx not found. "
            "Place a copy in the tests/ directory."
        )
    return path
