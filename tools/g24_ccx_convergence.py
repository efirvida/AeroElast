"""G24: paired convergence of the three CalculiX comparisons (axial, isotropic-equivalent, in-plane bending).

Predeclared criterion (written before the first run): each test is re-run on the plate mesh
refined 4x10 (the test), 8x20, 16x40, 32x80 by patching the mesh builder; the printed free-centre
displacements of both codes are read. The comparison closes when (a) both codes change < 1 % from
16x40 to 32x80, (b) the AeroElast-vs-CCX gap at 32x80 is <= the test bound (1.5 % axial and
isotropic-equivalent, 2.5 % bending), and (c) AeroElast is within 2 % of the printed closed form
at 32x80. Mesh is not an asserted quantity; a failing assertion at a refined mesh is caught and
reported, not hidden. Run through scripts/aeroenv.sh with PYTHONPATH=$SCRATCH/validenv.
"""
import contextlib
import io
import re
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, "tests/validation/parity")
sys.path.insert(0, ".")
import test_composite_beam_parity as m  # noqa: E402

MESHES = ((4, 10), (8, 20), (16, 40), (32, 80))
orig = m._build_composite_plate_mesh
TESTS = (
    ("axial", m.test_composite_axial_tension, 0.015),
    ("iso_equiv", m.test_composite_isotropic_equiv, 0.015),
    ("bending", m.test_composite_bending, 0.025),
)
NUM = r"([-+0-9.eE]+)"
PAT = re.compile(rf"(?:AeroElast|Isotropic equiv) \w+: {NUM} um, CCX: {NUM} um, (?:CLT bar|beam|[\w ]+): {NUM} um")

root = Path("logs/check/l2-audit")
for name, fn, bound in TESTS:
    res = {}
    for nx, ny in MESHES:
        m._build_composite_plate_mesh = lambda nx=nx, ny=ny: orig(nx, ny)
        buf = io.StringIO()
        status = "ok"
        with tempfile.TemporaryDirectory(dir=root) as td, contextlib.redirect_stdout(buf):
            try:
                fn(Path(td))
            except AssertionError as e:
                status = "assert: " + str(e).splitlines()[0][:80]
        mt = PAT.search(buf.getvalue())
        if not mt:
            print(f"{name} {nx}x{ny}: could not parse\n{buf.getvalue()[-400:]}")
            continue
        res[(nx, ny)] = tuple(float(v) for v in mt.groups())
        a, c, an = res[(nx, ny)]
        print(f"{name:9} {nx:>2}x{ny:<3} aero {a:9.3f} ccx {c:9.3f} closed {an:9.3f} um | "
              f"gap {100*abs(a-c)/c:5.2f}% aero-vs-closed {100*(a/an-1):+6.2f}% ccx-vs-closed {100*(c/an-1):+6.2f}% [{status}]",
              flush=True)
    (a1, c1, _), (a2, c2, an2) = res[MESHES[-2]], res[MESHES[-1]]
    da, dc = 100 * (a2 / a1 - 1), 100 * (c2 / c1 - 1)
    gap = abs(a2 - c2) / c2
    ok = abs(da) < 1 and abs(dc) < 1 and gap <= bound and abs(a2 / an2 - 1) <= 0.02
    print(f"  => {name}: 16x40->32x80 aero {da:+.3f}% ccx {dc:+.3f}%  gap {100*gap:.2f}% (bound {100*bound}%)  "
          f"{'CLOSED' if ok else 'OPEN'} under the predeclared criterion", flush=True)
