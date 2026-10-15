"""G4 in-plane bending (test_orthotropic_bending): mesh sweep against the closed forms.

Predeclared criterion (written before the first run). The test's [0/90/90/0] cantilever (L/B = 10,
uniform tip shear flow, mean tip Ux) is rerun on 4x10 (the test), 8x20, 16x40, 32x80 by patching the
mesh builder; the (computed, reference) pair handed to assert_relative_error is captured, not
asserted. The 3.648 % gap of the test is ATTRIBUTED when (a) the AeroElast tip changes < 1 % from
16x40 to 32x80 (converged) and (b) at the finest mesh it lies within 1 % of Euler-Bernoulli plus a
Timoshenko shear term P L / (kappa A66 B), kappa = 5/6 (reported, not gated, as kappa for an
orthotropic deep beam is approximate). If (a) fails the gap stays unattributed (bilinear membrane
locking). Run through scripts/aeroenv.sh with PYTHONPATH=$SCRATCH/validenv.
"""
import inspect
import sys

sys.path.insert(0, "tests/validation/parity")
sys.path.insert(0, ".")
import test_orthotropic_shell_parity as ortho  # noqa: E402

captured = []
ortho.assert_relative_error = lambda a, r, **kw: captured.append((float(a), float(r)))
src = inspect.getsource(ortho._build_plate)


def patched(nx, nz):
    ns = dict(vars(ortho))
    exec(src.replace("nx, nz = 4, 10", f"nx, nz = {nx}, {nz}"), ns)
    return ns["_build_plate"]()


lam = ortho.create_laminate_from_angles(ortho._ORTHO_MAT, ortho.thickness / 4, [0, 90, 90, 0])
A11, A66 = lam.A[0, 0], lam.A[2, 2]
B, L, F = ortho.B, ortho.L, 100.0
eb = F * L**3 / (3 * A11 * B**3 / 12)
ebt = eb + F * L / (5 / 6 * A66 * B)
print(f"closed forms: EB {eb*1e6:.3f} um, EB+Timoshenko {ebt*1e6:.3f} um (shear +{100*(ebt/eb-1):.2f} %)")
res = {}
for nx, nz in ((4, 10), (8, 20), (16, 40), (32, 80)):
    ortho._build_plate = lambda nx=nx, nz=nz: patched(nx, nz)
    captured.clear()
    ortho.test_orthotropic_bending()
    a, _ = captured[-1]
    res[(nx, nz)] = a
    print(f"{nx:>2}x{nz:<3} aero {a*1e6:11.3f} um | vs EB {100*(a/eb-1):+7.3f} %  vs EB+T {100*(a/ebt-1):+7.3f} %", flush=True)
k = list(res)
d = 100 * (res[k[-1]] / res[k[-2]] - 1)
g = 100 * (res[k[-1]] / ebt - 1)
ok = abs(d) < 1 and abs(g) <= 1
print(f"  => {k[-2]}->{k[-1]} change {d:+.3f} %, finest vs EB+T {g:+.3f} %: {'ATTRIBUTED' if ok else 'UNATTRIBUTED'} under the predeclared criterion")
