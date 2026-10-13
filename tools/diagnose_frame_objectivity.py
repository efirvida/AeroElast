"""Frame-objectivity of the shell tangent under a rigid rotation about z.

Compares K_t(u_rigid) from assemble_kt (TL) and assemble_kt_corotational against
T K0 T^T, with u_rigid carrying (a) translations only (what
tests/test_corotational_large_rotation_validation.py does) and (b) translations plus
theta_z = angle. Also reports the internal force of the rigid motion (must vanish).
Usage: diagnose_frame_objectivity.py <tests-dir>   (reuses the test module's builders)
"""
import importlib.util
import math
import sys

import numpy as np

spec = importlib.util.spec_from_file_location("lr", sys.argv[1] + "/test_corotational_large_rotation_validation.py")
lr = importlib.util.module_from_spec(spec)
spec.loader.exec_module(lr)

nodes, asm = lr._build_assembler()
n = nodes.shape[0]
nd = asm.dofs_count
dense = lambda t: lr._coo_to_dense(*t, nd)
k0 = dense(asm.assemble_kt(np.zeros(nd)))
fint_scale = np.linalg.norm(k0 @ np.ones(nd) * 1e-3)
print(f"dofs {nd}, nodes {n}")
for deg in (1.0, 5.0, 30.0):
    th = math.radians(deg)
    for label, rot in (("trans-only", 0.0), ("trans+theta_z", th)):
        u, r = lr._rigid_rotation_displacement(nodes, th)
        for i in range(n):
            u[6 * i + 5] = rot
        t6 = np.kron(np.eye(2 * n), r)       # rotate translation and rotation triples
        ref = t6 @ k0 @ t6.T
        ktl = dense(asm.assemble_kt(u))
        kc = dense(asm.assemble_kt_corotational(u))
        f = np.asarray(asm.assemble_fint(u, True))
        e = lambda a: np.linalg.norm(a - ref) / np.linalg.norm(ref)
        print(f"{deg:5.1f} deg {label:14s} err_coro {e(kc):.3e}  err_TL {e(ktl):.3e}  |fint|/|K0 u_scale| {np.linalg.norm(f)/fint_scale:.3e}")
