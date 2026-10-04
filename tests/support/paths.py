"""Path anchors for the suite.

A test file used to resolve its data with `Path(__file__).parent...`, which ties every path to
the depth the file happens to sit at: moving a file into a subdirectory silently repointed
`../examples` and `parents[1]` at the wrong place, and two collectors failed with a
`FileNotFoundError` that named a path nobody had written. These anchors are found by walking up
to the repository marker, so a move cannot break them again.
"""

from pathlib import Path


def _find_repo_root() -> Path:
    for parent in Path(__file__).resolve().parents:
        if (parent / "pyproject.toml").exists():
            return parent
    raise RuntimeError("no pyproject.toml above tests/support/paths.py")


REPO_ROOT = _find_repo_root()
TESTS_DIR = REPO_ROOT / "tests"
# Suite data has not moved: the reference turbine definition and the OpenFAST `reference/`
# tree still live at the tests root, not beside each test file.
DATA_DIR = TESTS_DIR
SOURCES_DIR = REPO_ROOT / ".sources"
EXAMPLES_DIR = REPO_ROOT / "examples"
OUTPUT_DIR = REPO_ROOT / "output"
