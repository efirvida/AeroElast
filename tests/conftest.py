"""Shared pytest fixtures for the fem-shell test suite."""

import ctypes
import os
import shutil

import pytest


# ---------------------------------------------------------------------------
# CalculiX (ccx) executable resolution
# ---------------------------------------------------------------------------
# Legacy HPC path, kept only as a last-resort fallback.  ``CCX_BIN`` and
# ``PATH`` take precedence, so a local install is found first.
_LEGACY_CCX_BIN = "/scratch/leahk/eduardo.donestevez/venv/bin/ccx"


def ccx_bin_or_skip() -> str:
    """Return a usable CalculiX (ccx) executable, or skip the calling test.

    Resolution order:

    1. ``CCX_BIN`` environment variable, when it points at an existing file.
    2. ``ccx`` on ``PATH``.
    3. ``CalculiX`` on ``PATH``.
    4. Legacy HPC fallback ``/scratch/leahk/eduardo.donestevez/venv/bin/ccx``.
    """
    candidates = []
    env_bin = os.environ.get("CCX_BIN")
    if env_bin:
        candidates.append(env_bin)
    candidates.append(shutil.which("ccx"))
    candidates.append(shutil.which("CalculiX"))
    candidates.append(_LEGACY_CCX_BIN)

    for candidate in candidates:
        if candidate and os.path.isfile(candidate):
            return candidate

    pytest.skip("CalculiX (ccx) not found; set CCX_BIN or install calculix")


# ---------------------------------------------------------------------------
# Large-deflection beam reference (shell-vs-beam nonlinear validation)
# ---------------------------------------------------------------------------


def elastica_cantilever_tip_deflection(P, L, b, h, E) -> float:
    """Geometrically nonlinear tip deflection of an end-loaded cantilever.

    Integrates the exact elastica ODE for a cantilever bent by a vertical tip
    load ``P``:

        EI * theta''(s) = -P * cos(theta(s)),   theta(0) = 0,  theta'(L) = 0

    where ``s`` is the arc length from the clamp, ``theta`` is the tangent
    angle from the undeflected axis and ``EI`` is the bending stiffness.  The
    tip deflection is ``integral_0^L sin(theta) ds``.  This is the Bisshopp &
    Drucker (1945) large-deflection beam solution, independent of the shell
    element under test.
    """
    import numpy as np
    from scipy.integrate import solve_bvp

    inertia = b * h**3 / 12.0
    ei = E * inertia

    def ode(s, y):
        return np.vstack([y[1], -(P / ei) * np.cos(y[0])])

    def bc(ya, yb):
        return np.array([ya[0], yb[1]])

    s = np.linspace(0.0, float(L), 4001)
    y0 = np.zeros((2, s.size))
    # Linear-beam initial guess: theta(s) = P (2 L s - s^2) / (2 EI).
    y0[0] = (P / (2.0 * ei)) * (2.0 * L * s - s**2)
    y0[1] = (P / ei) * (L - s)

    sol = solve_bvp(ode, bc, s, y0, tol=1e-10, max_nodes=20000)
    if sol.status != 0:
        raise RuntimeError(f"elastica BVP did not converge: {sol.message}")

    theta = sol.sol(s)[0]
    return float(np.trapezoid(np.sin(theta), s))


# ---------------------------------------------------------------------------
# GLU library — required by gmsh (equivalent to `module load glu`)
# ---------------------------------------------------------------------------
_GLU_LIB_PATH = "/scratch/app/glu/9.0.2_gnu/lib"

try:
    ctypes.CDLL("libGLU.so.1")
except OSError:
    # Not on LD_LIBRARY_PATH yet; prepend the HPC module path and retry.
    os.environ["LD_LIBRARY_PATH"] = _GLU_LIB_PATH + ":" + os.environ.get("LD_LIBRARY_PATH", "")
    try:
        ctypes.CDLL(os.path.join(_GLU_LIB_PATH, "libGLU.so.1"))
    except OSError:
        pass  # gmsh tests will fail with a clear error if still missing

# ---------------------------------------------------------------------------
# Blade YAML fixture
# ---------------------------------------------------------------------------
# Candidate paths for the IEA-15-240-RWT turbine definition YAML.
# The file is large so it lives outside the repo in the simulations tree,
# but a copy may also exist inside examples/.
_BLADE_YAML_CANDIDATES = [
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
