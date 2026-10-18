"""The repository gate on a node without ``git`` (SLURM compute nodes on SDumont).

Job 11617933 ran ``scripts/check.sh full --record`` on a compute node: the mesh-size gate died
with ``FileNotFoundError: 'git'`` and was reported as a coarse-mesh failure, and the known-red
registry was stamped "revision unknown". The gate now reads the git facts from a snapshot taken
on the login node (``scripts/check_sbatch.sh``), and without git or a snapshot it stops with an
environment error instead of a wrong verdict. ``CHECK_GIT`` names the git executable so a test
can simulate its absence.
"""

import os
import shutil
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
MESH_SIZE = ROOT / "scripts" / "check_mesh_size.py"
CHECK_SH = ROOT / "scripts" / "check.sh"
NO_GIT = "/nonexistent/git"


def _run(cmd, **env):
    full_env = {k: v for k, v in os.environ.items() if k not in ("CHECK_GIT", "CHECK_GIT_SNAPSHOT")}
    full_env.update(env)
    return subprocess.run(cmd, cwd=ROOT, env=full_env, capture_output=True, text=True, timeout=120)


def _snapshot(tmp_path: Path) -> Path:
    snap = tmp_path / "snap"
    snap.mkdir()
    (snap / "revision").write_text("abc1234\n")
    tracked = subprocess.check_output(
        ["git", "ls-files", "--cached", "--others", "--exclude-standard", "tests", "tools", "src", "scripts"],
        cwd=ROOT, text=True)
    (snap / "tracked").write_text(tracked)
    (snap / "changed_py").write_text("")
    return snap


def test_mesh_size_without_git_or_snapshot_is_an_environment_error():
    res = _run(["python3", str(MESH_SIZE)], CHECK_GIT=NO_GIT)
    assert res.returncode == 2, res.stdout + res.stderr
    assert "git" in res.stderr and "CHECK_GIT_SNAPSHOT" in res.stderr
    assert "Traceback" not in res.stderr
    assert "mesh-size: OK" not in res.stdout and "NEW coarse mesh" not in res.stdout


@pytest.mark.skipif(shutil.which("git") is None, reason="the reference run needs git (compute node)")
def test_mesh_size_from_snapshot_matches_the_git_run(tmp_path):
    with_git = _run(["python3", str(MESH_SIZE)])
    from_snapshot = _run(["python3", str(MESH_SIZE)], CHECK_GIT=NO_GIT, CHECK_GIT_SNAPSHOT=str(_snapshot(tmp_path)))
    assert (from_snapshot.returncode, from_snapshot.stdout, from_snapshot.stderr) == (
        with_git.returncode, with_git.stdout, with_git.stderr)


def test_check_sh_without_git_or_snapshot_stops_before_pytest():
    res = _run(["bash", str(CHECK_SH), "quick"], CHECK_GIT=NO_GIT)
    assert res.returncode == 2, res.stdout + res.stderr
    assert "git" in res.stderr and "scripts/check_sbatch.sh" in res.stderr
    assert "pytest" not in res.stdout


def test_check_sh_rejects_a_snapshot_without_its_files(tmp_path):
    res = _run(["bash", str(CHECK_SH), "quick"], CHECK_GIT=NO_GIT, CHECK_GIT_SNAPSHOT=str(tmp_path))
    assert res.returncode == 2, res.stdout + res.stderr
    assert "snapshot" in res.stderr
