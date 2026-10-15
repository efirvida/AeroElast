"""G4 (orthotropic transverse) and G8 (isotropic transverse): paired convergence of the CalculiX comparison.

Predeclared criterion (written before the first run). Each test is re-run on the plate mesh refined
G4 4x10 (the test), 8x20, 16x40, 32x80 and G8 2x10 (the test), 4x20, 8x40, 16x80 by patching the mesh
builder; the two displacements the test hands to assert_relative_error (AeroElast, CalculiX) are
captured, not asserted. Both are an out-of-plane POINT load at the free-edge centre, so the loaded
node is the metric the row asserts. The comparison is CLOSED when (a) both codes change < 1 % from
the second-to-last to the last mesh and (b) the gap at the finest mesh is <= 5 % (the test bound).
If (a) fails the metric is not converged (point-load singularity, as measured for G17 axial) and the
row stays code-parity only. Closed forms are reported for context, not gated: Euler-Bernoulli
P L^3/(3 E I) and EB + Timoshenko shear P L/(kappa G A), kappa = 5/6. G4 reports EB for both
laminate flexural moduli (D11 and D22) because the fibre axis relative to the beam axis is not
asserted by the test. Run through scripts/aeroenv.sh with PYTHONPATH=$SCRATCH/validenv.
"""
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, "tests/validation/parity")
sys.path.insert(0, ".")
import test_isotropic_shell_parity as iso  # noqa: E402
import test_orthotropic_shell_parity as ortho  # noqa: E402

root = Path("logs/check/l2-audit")
captured = []


def grab(actual, reference, **kw):
    captured.append((float(actual), float(reference)))


def run(label, call, meshes, patch, closed):
    res = {}
    for nx, nz in meshes:
        patch(nx, nz)
        captured.clear()
        with tempfile.TemporaryDirectory(dir=root) as td:
            call(Path(td))
        if not captured:
            print(f"{label} {nx}x{nz}: nothing captured (skipped?)", flush=True)
            continue
        a, c = captured[-1]
        res[(nx, nz)] = (a, c)
        print(f"{label} {nx:>2}x{nz:<3} aero {a*1e6:11.3f} ccx {c*1e6:11.3f} um | "
              f"gap {100*abs(a-c)/abs(c):5.2f}%  aero/EB {a/closed['eb']:.4f}  aero/EB+T {a/closed['ebt']:.4f}",
              flush=True)
    k = list(res)
    (a1, c1), (a2, c2) = res[k[-2]], res[k[-1]]
    da, dc, gap = 100 * (a2 / a1 - 1), 100 * (c2 / c1 - 1), abs(a2 - c2) / abs(c2)
    ok = abs(da) < 1 and abs(dc) < 1 and gap <= 0.05
    print(f"  => {label}: {k[-2]}->{k[-1]} aero {da:+.3f}% ccx {dc:+.3f}% gap {100*gap:.2f}% "
          f"{'CLOSED' if ok else 'OPEN'} under the predeclared criterion", flush=True)


# G8: isotropic steel plate
ortho_build, iso_build = ortho._build_plate, iso.build_shell_mesh
iso.assert_relative_error = grab
E, nu, B, L, t, P = iso.MAT_STEEL.E, iso.MAT_STEEL.nu, iso.B, iso.L, iso.t, iso.FORCE
I, G, A = B * t**3 / 12, E / (2 * (1 + nu)), B * t
eb = P * L**3 / (3 * E * I)
closed = {"eb": eb, "ebt": eb + P * L / (5 / 6 * G * A)}
print(f"G8 closed form: EB {eb*1e6:.3f} um, EB+Timoshenko {closed['ebt']*1e6:.3f} um")
run("G8", lambda p: iso.TestIsotropicShellParity().test_transverse_tip_displacement(p),
    ((2, 10), (4, 20), (8, 40), (16, 80)),
    lambda nx, nz: setattr(iso, "build_shell_mesh", lambda **kw: iso_build(nx=nx, nz=nz)), closed)

# G4: single 0-degree ply laminate
import inspect  # noqa: E402
_src = inspect.getsource(ortho_build)


def _patched(nx, nz):
    ns = dict(vars(ortho))
    exec(_src.replace("nx, nz = 4, 10", f"nx, nz = {nx}, {nz}"), ns)
    return ns["_build_plate"]()


ortho.assert_relative_error = grab
lam = ortho.create_laminate_from_angles(ortho._ORTHO_MAT, ortho.thickness, [0])
D = lam.D
B4, L4, P4 = ortho.B, ortho.L, 100.0
for name, d in (("D11", D[0][0]), ("D22", D[1][1])):
    # D is the integrated bending stiffness (N m); beam EI = D * B
    eb4 = P4 * L4**3 / (3 * d * B4)
    print(f"G4 closed form EB with {name}: {eb4*1e6:.3f} um")
closed4 = {"eb": P4 * L4**3 / (3 * D[0][0] * B4)}
closed4["ebt"] = closed4["eb"]
run("G4", ortho.test_orthotropic_axial, ((4, 10), (8, 20), (16, 40), (32, 80)),
    lambda nx, nz: setattr(ortho, "_build_plate", lambda: _patched(nx, nz)), closed4)
