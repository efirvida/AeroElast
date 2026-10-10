"""The guard in `conftest.py` that keeps the tool suite off the shipped store (#32).

A green suite cannot show the guard works: a test that writes identical bytes, or none, trips
nothing. These tests drive the guard's own pieces against a scratch tree instead.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path

_SPEC = importlib.util.spec_from_file_location("tools_tests_conftest", Path(__file__).with_name("conftest.py"))
assert _SPEC is not None and _SPEC.loader is not None
guard = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(guard)


def make_tree(root: Path) -> Path:
    (root / "rows").mkdir(parents=True)
    (root / "rows" / "1-a.yaml").write_text("group: '1'\n", encoding="utf-8")
    (root / "rows" / "2-b.yaml").write_text("group: '2'\n", encoding="utf-8")
    (root / "groups.yaml").write_text("groups: []\n", encoding="utf-8")
    return root


def test_an_untouched_tree_has_no_changes(tmp_path: Path) -> None:
    root = make_tree(tmp_path / "store")
    before = guard.snapshot(root)
    (root / "rows" / "1-a.yaml").write_text("group: '1'\n", encoding="utf-8")  # same bytes
    assert guard.changes(before, guard.snapshot(root)) == []


def test_added_removed_and_modified_files_are_named_and_restored(tmp_path: Path) -> None:
    root = make_tree(tmp_path / "store")
    before = guard.snapshot(root)

    (root / "rows" / "1-a.yaml").write_text("group: '1'\nrows: []\n", encoding="utf-8")
    (root / "rows" / "2-b.yaml").unlink()
    (root / "rows" / "3-c.yaml").write_text("group: '3'\n", encoding="utf-8")
    after = guard.snapshot(root)

    assert guard.changes(before, after) == [
        "added    rows/3-c.yaml",
        "removed  rows/2-b.yaml",
        "modified rows/1-a.yaml",
    ]
    guard.restore(root, before, after)
    assert guard.snapshot(root) == before


def test_the_guard_watches_the_shipped_store() -> None:
    assert guard.REAL_STORE == Path(__file__).resolve().parents[2] / "docs" / "validation"
    assert (guard.REAL_STORE / "rows").is_dir()
