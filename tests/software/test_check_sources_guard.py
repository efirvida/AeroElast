"""The gate must not read a skipped reference test as a passing one.

The full-suite registry recorded at a9ca498 (SLURM 11623882) ran in a worktree without the
git-ignored ``.sources/`` directory (the BeamDyn and ElastoDyn decks). Every test that needs a
deck skipped, ``--record`` listed none of them as red, and the next main-tree run reported three
long-standing failures (the section frame-tear guards, the pure-torsion GJ guard) as NEW-FAIL.
Two rules close that: a full run or any ``--record`` stops when ``.sources`` is missing, and the
gate keeps a known-skip registry next to the known-red one, so a test that ran before and now
skips is reported (NEW-SKIP) instead of disappearing.

``CHECK_SOURCES_DIR``, ``CHECK_REGISTRY_DIR`` and ``CHECK_SELECTION`` point the gate at a test
directory, a scratch registry and a two-test file, so these tests never touch the real
registries or run the real selection.
"""

import os
import shutil
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
CHECK_SH = ROOT / "scripts" / "check.sh"
MISSING = "/nonexistent/.sources"

pytestmark = pytest.mark.skipif(shutil.which("git") is None, reason="check.sh needs git or a snapshot")

ONE_SKIP = """\
import pytest


def test_runs():
    pass


@pytest.mark.skip(reason="no deck")
def test_needs_a_deck():
    pass
"""

TWO_SKIPS = ONE_SKIP.replace("def test_runs", '@pytest.mark.skip(reason="no deck")\ndef test_runs')
NO_SKIP = ONE_SKIP.replace('@pytest.mark.skip(reason="no deck")\n', "")


def _run(*args, **env):
    full_env = {k: v for k, v in os.environ.items() if not k.startswith("CHECK_")}
    full_env.update(env)
    return subprocess.run(["bash", str(CHECK_SH), *args], cwd=ROOT, env=full_env,
                          capture_output=True, text=True, timeout=600)


@pytest.fixture()
def sandbox(tmp_path):
    registry = tmp_path / "registry"
    registry.mkdir()
    sources = tmp_path / "sources"
    sources.mkdir()
    selection = tmp_path / "test_selection_probe.py"
    return registry, sources, selection


def _skip_entries(path: Path) -> list[str]:
    return [ln for ln in path.read_text().splitlines() if ln and not ln.startswith("#")]


def test_full_without_sources_stops_before_pytest():
    res = _run("full", CHECK_SOURCES_DIR=MISSING)
    assert res.returncode == 2, res.stdout + res.stderr
    assert MISSING in res.stderr and "skip" in res.stderr
    assert "pytest:" not in res.stdout


def test_record_without_sources_leaves_the_registries_unchanged(sandbox):
    registry, _, selection = sandbox
    selection.write_text(ONE_SKIP)
    red = registry / "known-red-quick.txt"
    red.write_text("# recorded elsewhere\ntests/x.py::test_y\n")
    before = red.read_bytes()
    res = _run("quick", "--record", CHECK_SOURCES_DIR=MISSING, CHECK_REGISTRY_DIR=str(registry),
               CHECK_SELECTION=str(selection))
    assert res.returncode == 2, res.stdout + res.stderr
    assert MISSING in res.stderr
    assert "pytest:" not in res.stdout and "recorded" not in res.stdout
    assert red.read_bytes() == before
    assert not (registry / "known-skip-quick.txt").exists()


def test_record_lists_the_skipped_tests_and_a_new_skip_fails(sandbox):
    registry, sources, selection = sandbox
    selection.write_text(ONE_SKIP)
    rec = _run("quick", "--record", CHECK_SOURCES_DIR=str(sources), CHECK_REGISTRY_DIR=str(registry),
               CHECK_SELECTION=str(selection))
    assert rec.returncode == 0, rec.stdout + rec.stderr
    assert "recorded 1 known-skip entry(ies)" in rec.stdout
    skips = registry / "known-skip-quick.txt"
    entries = _skip_entries(skips)
    assert len(entries) == 1 and entries[0].endswith("::test_needs_a_deck"), entries
    recorded = skips.read_bytes()

    # quick without --record still runs without .sources; the test that ran now skips
    selection.write_text(TWO_SKIPS)
    res = _run("quick", CHECK_SOURCES_DIR=MISSING, CHECK_REGISTRY_DIR=str(registry),
               CHECK_SELECTION=str(selection))
    assert res.returncode == 1, res.stdout + res.stderr
    new_skips = [ln for ln in res.stdout.splitlines() if ln.strip().startswith("NEW-SKIP")]
    assert len(new_skips) == 1 and new_skips[0].endswith("::test_runs"), res.stdout
    assert "1 test(s) skipped outside" in res.stdout
    assert skips.read_bytes() == recorded


def test_a_recorded_skip_that_runs_again_is_reported_not_failed(sandbox):
    registry, sources, selection = sandbox
    selection.write_text(ONE_SKIP)
    assert _run("quick", "--record", CHECK_SOURCES_DIR=str(sources), CHECK_REGISTRY_DIR=str(registry),
                CHECK_SELECTION=str(selection)).returncode == 0
    selection.write_text(NO_SKIP)
    res = _run("quick", CHECK_SOURCES_DIR=str(sources), CHECK_REGISTRY_DIR=str(registry),
               CHECK_SELECTION=str(selection))
    assert res.returncode == 0, res.stdout + res.stderr
    now_runs = [ln for ln in res.stdout.splitlines() if ln.strip().startswith("now-runs")]
    assert len(now_runs) == 1 and "::test_needs_a_deck" in now_runs[0], res.stdout
    assert "NEW-SKIP" not in res.stdout


def test_without_a_skip_registry_skips_do_not_fail(sandbox):
    registry, sources, selection = sandbox
    selection.write_text(TWO_SKIPS)
    res = _run("quick", CHECK_SOURCES_DIR=str(sources), CHECK_REGISTRY_DIR=str(registry),
               CHECK_SELECTION=str(selection))
    assert res.returncode == 0, res.stdout + res.stderr
    assert "NEW-SKIP" not in res.stdout
    assert "check[quick]: OK" in res.stdout


# ---- --only: re-measure a few tests without the full suite ------------------------------------

BOTH_FAIL = """\
def test_passes():
    assert False


def test_fails():
    assert False
"""
ONE_HEALED = BOTH_FAIL.replace("def test_passes():\n    assert False", "def test_passes():\n    pass")
KEPT = "tests/elsewhere.py::test_kept"


def _red_entries(registry: Path) -> list[str]:
    return _skip_entries(registry / "known-red-quick.txt")


def _record_both_red_plus_one_outside(registry, sources, selection):
    """A registry holding the probe's two failures and one entry the probe does not run."""
    selection.write_text(BOTH_FAIL)
    rec = _run("quick", "--record", CHECK_SOURCES_DIR=str(sources), CHECK_REGISTRY_DIR=str(registry),
               CHECK_SELECTION=str(selection))
    assert rec.returncode == 0, rec.stdout + rec.stderr
    red = registry / "known-red-quick.txt"
    red.write_text(red.read_text() + KEPT + "\n")
    entries = _red_entries(registry)
    assert len(entries) == 3, entries
    return red


def test_only_record_rewrites_its_scope_and_keeps_the_rest(sandbox):
    registry, sources, selection = sandbox
    red = _record_both_red_plus_one_outside(registry, sources, selection)
    header = [ln for ln in red.read_text().splitlines() if ln.startswith("#")]
    selection.write_text(ONE_HEALED)
    res = _run("quick", "--only", str(selection), "--record", CHECK_SOURCES_DIR=str(sources),
               CHECK_REGISTRY_DIR=str(registry))
    assert res.returncode == 0, res.stdout + res.stderr
    entries = _red_entries(registry)
    assert len(entries) == 2 and KEPT in entries, entries
    assert any(e.endswith("::test_fails") for e in entries), entries
    assert not any(e.endswith("::test_passes") for e in entries), entries
    lines = red.read_text().splitlines()
    assert lines[: len(header)] == header
    assert any(ln.startswith("# Partially re-derived") and str(selection) in ln for ln in lines), lines
    assert "recorded 1 known-red entry(ies) in scope (2 test(s)), kept 1 outside" in res.stdout
    # an empty scope (no skips here) is still written, and no scratch file is left behind
    skip_lines = (registry / "known-skip-quick.txt").read_text().splitlines()
    assert any(ln.startswith("# Partially re-derived") for ln in skip_lines), skip_lines
    assert not list(registry.glob("*.tmp")), list(registry.glob("*.tmp"))


def test_only_reports_healed_entries_inside_its_scope_alone(sandbox):
    registry, sources, selection = sandbox
    red = _record_both_red_plus_one_outside(registry, sources, selection)
    before = red.read_bytes()
    selection.write_text(ONE_HEALED)
    res = _run("quick", "--only", str(selection), CHECK_SOURCES_DIR=str(sources),
               CHECK_REGISTRY_DIR=str(registry))
    assert res.returncode == 0, res.stdout + res.stderr
    healed = [ln for ln in res.stdout.splitlines() if ln.strip().startswith("now-passes")]
    assert len(healed) == 1 and "::test_passes" in healed[0], res.stdout
    assert "--only" in healed[0], healed[0]
    assert KEPT not in res.stdout
    assert red.read_bytes() == before


def test_only_without_a_selector_is_an_argument_error():
    res = _run("quick", "--only")
    assert res.returncode == 2
    assert "--only needs a test file or node id" in res.stderr
