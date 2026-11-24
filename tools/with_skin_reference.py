"""Run any script or pytest with the BladeMesh / GmshBladeMesh default skin_reference forced to oml or mid.

For before/after anchor measurements (odd/tasks/blade-midsurface-reference.md, T2): the same
command, once per reference surface, with no test or tool edited. Only meshes built WITHOUT an
explicit skin_reference take the forced default.

    scripts/aeroenv.sh python tools/with_skin_reference.py mid -m pytest tests/x.py -q -rA
    scripts/aeroenv.sh python tools/with_skin_reference.py oml tools/s4_modal_convergence.py 0.25
"""

from __future__ import annotations

import runpy
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(REPO_ROOT / "tools"))


def main() -> int:
    mode, rest = sys.argv[1], sys.argv[2:]
    if mode not in ("oml", "mid") or not rest:
        print(__doc__)
        return 2
    from aeroelast.core.mesh import generators
    from aeroelast.core.mesh.gmsh_blade import GmshBladeMesh

    # both blade constructors take a skin_reference default: BladeMesh (oml, builds with gmsh) and GmshBladeMesh (mid)
    for cls in (generators.BladeMesh, GmshBladeMesh):
        original = cls.__init__

        def patched(self, *args, _orig=original, **kwargs):
            kwargs.setdefault("skin_reference", mode)
            _orig(self, *args, **kwargs)

        cls.__init__ = patched
    print(f"with_skin_reference: BladeMesh / GmshBladeMesh default skin_reference = {mode}", flush=True)
    if rest[0] == "-m":
        sys.argv = rest[1:]
        runpy.run_module(rest[1], run_name="__main__", alter_sys=True)
    else:
        sys.argv = rest
        runpy.run_path(rest[0], run_name="__main__")
    return 0


if __name__ == "__main__":
    sys.exit(main())
