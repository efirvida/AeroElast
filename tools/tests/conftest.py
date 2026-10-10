"""Shared fixtures for the tool's own tests.

The validation store under `docs/validation` is validated data, not a fixture. A test that
needs the shipped store reads it, or works on a copy (`real_store_copy`); it never writes it.
`store_is_never_written` holds every test in this directory to that rule (#32).
"""

from __future__ import annotations

import shutil
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
REAL_STORE = REPO_ROOT / "docs" / "validation"


def snapshot(root: Path) -> dict[str, bytes]:
    """Every file under `root`, keyed by its relative path, with its bytes."""
    return {
        path.relative_to(root).as_posix(): path.read_bytes()
        for path in sorted(root.rglob("*"))
        if path.is_file() and "__pycache__" not in path.parts
    }


def changes(before: dict[str, bytes], after: dict[str, bytes]) -> list[str]:
    """The paths that were added, removed or rewritten with different bytes."""
    lines = [f"added    {path}" for path in sorted(after.keys() - before.keys())]
    lines += [f"removed  {path}" for path in sorted(before.keys() - after.keys())]
    lines += [
        f"modified {path}"
        for path in sorted(before.keys() & after.keys())
        if before[path] != after[path]
    ]
    return lines


def restore(root: Path, before: dict[str, bytes], after: dict[str, bytes]) -> None:
    """Put `root` back to `before`, so a violating test does not leave damage in the tree."""
    for path in after.keys() - before.keys():
        (root / path).unlink()
    for path, data in before.items():
        if after.get(path) != data:
            (root / path).parent.mkdir(parents=True, exist_ok=True)
            (root / path).write_bytes(data)


@pytest.fixture(autouse=True)
def store_is_never_written() -> object:
    """Fail any test that changes the shipped store, and undo the change."""
    before = snapshot(REAL_STORE)
    yield
    after = snapshot(REAL_STORE)
    changed = changes(before, after)
    if changed:
        restore(REAL_STORE, before, after)
        pytest.fail(
            "the test wrote the real store (restored); run against a copy instead:\n  "
            + "\n  ".join(changed),
            pytrace=False,
        )


@pytest.fixture
def real_store_copy(tmp_path: Path) -> Path:
    """A private copy of the shipped store, for commands that write (`coherence`, `--write`)."""
    store = tmp_path / "store"
    shutil.copytree(REAL_STORE, store)
    return store
