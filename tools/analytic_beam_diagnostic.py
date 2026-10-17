"""Instrument the unchanged S-4 path on a small closed shell box.

Run on sequana_cpu with OMP/OPENBLAS/RAYON thread limits. A 120-second
faulthandler timer locates native stalls; a Python trace identifies each S-4
stage. PETSC_OPTIONS='-eps_monitor -eps_converged_reason' adds native progress.
Arguments: span subdivisions, wall subdivisions, mode count, Omega [rad/s].
"""
import faulthandler
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from tests.validation.benchmarks.analytic_beam_helpers import (  # noqa: E402
    LENGTH, RHO, THICKNESS, assembled_mass, build_box,
)
from tools.run_s4_rotating_modal import assemble_rotating  # noqa: E402

faulthandler.enable()
faulthandler.dump_traceback_later(120, repeat=True)
span, wall, count = map(int, sys.argv[1:4]) if len(sys.argv) > 1 else (40, 2, 12)
omega = float(sys.argv[4]) if len(sys.argv) > 4 else 0.0
case = build_box(span, wall, triangles=True)
asm, mesh, indices, _ = case
fixed = (indices[0, :, None]*6 + np.arange(6)).ravel()
print(f"mesh={span}x{wall}; nnode={mesh.node_count}; ndof={asm.dofs_count}; "
      f"root={sorted(mesh.get_node_set('RootNodes').nodes)}; fixed={fixed.tolist()}; "
      f"Omega={omega}; mass={assembled_mass(case)}", flush=True)
start = time.monotonic()
print(f"centrifugal density argument={np.unique(asm._rho_per_elem)}; "
      f"analytic areal density={RHO*THICKNESS} kg/m2", flush=True)
if omega:
    force = np.asarray(asm._rust.centrifugal_load(
        omega, [0., 1., 0.], [0., 0., 0.], asm._rho_per_elem)).reshape(-1, 6)
    analytic = assembled_mass(case)[2] * omega**2 * LENGTH/2
    print(f"centrifugal sum Fz={force[:, 2].sum()}; analytic={analytic}; "
          f"ratio={force[:, 2].sum()/analytic}", flush=True)

if len(sys.argv) > 5 and sys.argv[5].startswith("preflight"):
    from petsc4py import PETSc
    from scipy.linalg import eigh
    from scipy.sparse import coo_matrix
    from scipy.sparse.linalg import splu
    from slepc4py import SLEPc

    free = np.setdiff1d(np.arange(asm.dofs_count), fixed)
    matrices = []
    for name, triplet in [("K", asm._rust.assemble_k()), ("M", asm._rust.assemble_m())]:
        rows, cols, vals = triplet
        matrix = coo_matrix((vals, (rows, cols)), shape=(asm.dofs_count,)*2).tocsr()[free][:, free]
        print(f"{name}: shape={matrix.shape}; diag_min={matrix.diagonal().min()}; "
              f"asymmetry={abs(matrix-matrix.T).max()}", flush=True)
        lu = splu(matrix.tocsc())
        print(f"{name}: SuperLU complete; min_pivot={abs(lu.U.diagonal()).min()}", flush=True)
        matrices.append(matrix)
    k, m = matrices
    print("dense generalized eigenproblem starting", flush=True)
    lam = eigh(k.toarray(), m.toarray(), subset_by_index=[0, 11], eigvals_only=True)
    print(f"dense frequencies={np.sqrt(lam)/(2*np.pi)}", flush=True)
    petsc_matrices = [PETSc.Mat().createAIJ(size=a.shape, csr=(a.indptr, a.indices, a.data),
                                          comm=PETSc.COMM_SELF) for a in matrices]
    eps = SLEPc.EPS().create(PETSc.COMM_SELF)
    eps.setOperators(*petsc_matrices)
    eps.setProblemType(SLEPc.EPS.ProblemType.GHEP)
    eps.setWhichEigenpairs(SLEPc.EPS.Which.TARGET_MAGNITUDE)
    eps.setTarget(0)
    eps.setDimensions(22)
    eps.setTolerances(1e-10, 1000)
    st = eps.getST()
    st.setType("sinvert")
    st.getKSP().setType("preonly")
    st.getKSP().getPC().setType("lu")
    eps.setMonitor(lambda _, iteration, nconv, eigenvalues, errors:
                   print(f"EPS iteration={iteration}; nconv={nconv}", flush=True))
    print("Python SLEPc starting", flush=True)
    eps.solve()
    print(f"Python SLEPc done: reason={eps.getConvergedReason()}; nconv={eps.getConverged()}", flush=True)
    print(f"preflight seconds={time.monotonic()-start:.3f}", flush=True)
    if sys.argv[5] == "preflight":
        faulthandler.cancel_dump_traceback_later()
        sys.exit(0)


def trace(frame, event, arg):
    if frame.f_code.co_name == "assemble_rotating" and event == "line":
        print(f"S4:{frame.f_lineno} elapsed={time.monotonic()-start:.3f}s", flush=True)
    return trace


if len(sys.argv) > 5 and sys.argv[5] == "ordering-rcm":
    from petsc4py import PETSc

    PETSc.Options()["st_pc_factor_mat_ordering_type"] = "rcm"
    print("Only change: native LU graph ordering=rcm (instead of default ND)", flush=True)

if len(sys.argv) > 5 and sys.argv[5] == "native-trace":
    import ctypes
    from petsc4py import PETSc

    lib = ctypes.CDLL("libpetsc.so.3.25")
    libc = ctypes.CDLL(None)
    stdout = ctypes.c_void_p.in_dll(libc, "stdout")
    lib.PetscLogTraceBegin.argtypes = [ctypes.c_void_p]
    lib.PetscLogTraceBegin.restype = ctypes.c_int
    print(f"PETSc trace start status={lib.PetscLogTraceBegin(stdout)}", flush=True)

sys.settrace(trace)
frequencies, modes, free = assemble_rotating(asm, mesh, omega, count)
sys.settrace(None)
print(f"frequencies={frequencies}; nfree={len(free)}; modes={modes.shape}", flush=True)
np.testing.assert_array_equal(free, np.setdiff1d(np.arange(asm.dofs_count), fixed))
print(f"finished in {time.monotonic()-start:.3f}s", flush=True)
faulthandler.cancel_dump_traceback_later()
