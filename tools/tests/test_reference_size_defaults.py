"""S9 / AGENTS.md rule 9: the blade convergence tools default to element sizes below 0.3 m."""

import re
from pathlib import Path

import pytest

TOOLS = Path(__file__).resolve().parents[1]
CASES = {
    "diagnose_s4_mode_labels.py": r'"--sizes".*default=\[([^\]]+)\]',
    "run_sx_convergence.py": r'"--sizes".*default=\[([^\]]+)\]',
    "diagnose_g5_mode7.py": r'"--meshes".*default=\[([^\]]+)\]',
    "run_s8b_semichord.py": r'"--element-size".*default=([0-9.]+)',
}


@pytest.mark.parametrize(("name", "pattern"), CASES.items())
def test_default_sizes_below_0p3(name, pattern):
    m = re.search(pattern, (TOOLS / name).read_text())
    assert m, f"{name}: size default not found"
    sizes = [float(x) for x in m.group(1).split(",")]
    assert sizes and all(0.05 <= s < 0.3 for s in sizes), (name, sizes)
