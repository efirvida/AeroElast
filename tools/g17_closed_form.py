"""G17: closed forms and paired convergence for the beam-shell 4-cases parity (AeroElast MITC4 vs CalculiX S4).

Predeclared criteria (written before the first run; PRE-3 analytic-first rule):

Statics, per case (tension, compression, bending_fx, transverse_fy), meshes (nx, nz) =
(4,20) [the test], (8,40), (16,80), both codes, same load/material/BC:
  C1 (closed-form gate, metric = free-face mean, the quantity beam theory predicts): both codes
     within 2 % of the closed form at 4x20 (the test bound) AND converged: the value changes
     < 1 % from 8x40 to 16x80.
  C2 (metric the CalculiX comparison actually asserts = displacement of the loaded node): the
     row closes by paired convergence only if the value changes < 1 % from 8x40 to 16x80 in BOTH
     codes. Otherwise the metric is not converged (point load) and the code-vs-code agreement is
     reported as a same-mesh coincidence; closure then rests on C1.
Modal, per mesh, both codes, 12 modes searched, 5 best globally matched pairs (as the test):
  C3 the five matched modes change < 1 % from 8x40 to 16x80 in both codes, and the matched gap
     at 16x80 is <= 3 %.
  C4 every matched mode that lies within 10 % of an Euler-Bernoulli / bar / thin-plate-torsion
     closed form (out-of-plane bending 1-3, in-plane bending 1-2, axial 1, torsion 1) is within
     5 % of it at 16x80, in both codes.
Run on a compute node through scripts/aeroenv.sh with PYTHONPATH=$SCRATCH/validenv.
"""
import sys
import tempfile
from pathlib import Path

import numpy as np
from scipy.optimize import linear_sum_assignment

sys.path.insert(0, "tests/validation/parity")
sys.path.insert(0, ".")
import test_beam_shell_4cases_parity as m  # noqa: E402

MESHES = ((4, 20), (8, 40), (16, 80))
_orig_build = m._build_shell_mesh
WORK = Path("logs/check/l2-audit")
WORK.mkdir(parents=True, exist_ok=True)


def use_mesh(nx, nz):
    m._build_shell_mesh = lambda nz=nz, nx=nx: _orig_build(nz=nz, nx=nx)


def ccx_face_and_node(case, mesh, workdir):
    center = next(iter(mesh.get_node_set("free_center").nodes.values()))
    dof = {"x": 0, "y": 1, "z": 2}[case.analytical_dof]
    stem = f"shell_{case.name}"
    d = workdir / case.name
    d.mkdir(exist_ok=True)
    m.write_ccx_mesh(
        mesh,
        str(d / f"{stem}.inp"),
        properties={
            "plate": {
                "type": "isotropic",
                "e": m.MAT_STEEL.E,
                "nu": m.MAT_STEEL.nu,
                "rho": m.MAT_STEEL.rho,
                "thickness": m.t,
            }
        },
        boundary_nodeset="clamped",
        solver_type="LinearStatic",
        load_nodeset="free_center",
        load_vector=list(case.load_vector),
    )
    m._run_ccx(m._ccx_bin(), d, stem)
    frd = d / f"{stem}.frd"
    node = m._frd_disp_at_xyz(frd, float(center.x), m.L)[dof]
    face = [
        m._frd_disp_at_xyz(frd, float(n.x), m.L)[dof]
        for n in mesh.get_node_set("free_face").nodes.values()
    ]
    return float(node), float(np.mean(face))


def statics(tmp):
    cls = m.TestBeamShell4CasesParity()
    for case in m.LOAD_CASES:
        anal = m._analytical_solution(case)
        vals = {}
        for nx, nz in MESHES:
            use_mesh(nx, nz)
            mesh, ae_node, ae_mean = cls._solve_static_case(case)
            cc_node, cc_mean = ccx_face_and_node(case, mesh, tmp)
            vals[(nx, nz)] = (ae_node, ae_mean, cc_node, cc_mean)
        print(f"== {case.name}  closed form {anal:.6E}")
        for k, (an, am, cn, cm) in vals.items():
            print(
                f"  {k[0]:>2}x{k[1]:<3} face-mean err: aero {100*(am/anal-1):+7.3f}% "
                f"ccx {100*(cm/anal-1):+7.3f}% | loaded node vs closed form: aero "
                f"{100*(an/anal-1):+7.3f}% ccx {100*(cn/anal-1):+7.3f}% | node AE-vs-CCX "
                f"{100*abs(an-cn)/abs(cn):.2f}%"
            )
        (a1, f1, c1, g1), (a2, f2, c2, g2) = vals[MESHES[1]], vals[MESHES[2]]
        ch = lambda x, y: 100 * (y / x - 1)
        c1_ok = all(
            abs(v / anal - 1) <= 0.02 for v in (vals[MESHES[0]][1], vals[MESHES[0]][3])
        ) and abs(ch(f1, f2)) < 1 and abs(ch(g1, g2)) < 1
        c2_ok = abs(ch(a1, a2)) < 1 and abs(ch(c1, c2)) < 1
        print(
            f"  8x40->16x80 change: face-mean aero {ch(f1,f2):+.3f}% ccx {ch(g1,g2):+.3f}% | "
            f"loaded node aero {ch(a1,a2):+.3f}% ccx {ch(c1,c2):+.3f}%"
        )
        print(f"  => C1 {'PASS' if c1_ok else 'FAIL'}   C2 {'converged' if c2_ok else 'NOT converged'}", flush=True)


def anchors():
    E, rho, nu, B, t, L = m.MAT_STEEL.E, m.MAT_STEEL.rho, m.MAT_STEEL.nu, m.B, m.t, m.L
    G = E / (2 * (1 + nu))
    A = B * t
    bl = (1.8751040687, 4.6940911330, 7.8547574382)
    fb = lambda I, b: b**2 / (2 * np.pi * L**2) * np.sqrt(E * I / (rho * A))
    out = {f"out{i+1}": fb(m.I_x, b) for i, b in enumerate(bl)}
    out.update({f"in{i+1}": fb(m.I_y, b) for i, b in enumerate(bl[:2])})
    out["axial1"] = np.sqrt(E / rho) / (4 * L)
    J = B * t**3 / 3
    Ip = rho * (B * t**3 + t * B**3) / 12
    out["torsion1"] = np.sqrt(G * J / Ip) / (4 * L)
    return out


def modal(tmp):
    from _aeroelast import PyMeshAssembler, modal_solve_coo

    def aero(mesh):
        xyz = np.asarray([[n.x, n.y, n.z] for n in mesh.nodes], float)
        conn = [[mesh.node_id_to_index[i] for i in e.node_ids] for e in mesh.elements]
        mats = [
            {"type": "isotropic", "e": m.MAT_STEEL.E, "nu": m.MAT_STEEL.nu,
             "rho": m.MAT_STEEL.rho, "thickness": m.t}
        ] * len(conn)
        asm = PyMeshAssembler(node_coords=xyz, connectivity=conn, elem_types=[4] * len(conn), materials=mats)
        kr, kc, kv = asm.assemble_k()
        mr, mc, mv = asm.assemble_m()
        n = asm.dofs_count
        fixed = set(m._clamped_dofs(mesh))
        free = np.array([i for i in range(n) if i not in fixed], dtype=np.int64)
        f, _ = modal_solve_coo(kr.astype(np.int64), kc.astype(np.int64), kv.astype(float),
                               mr.astype(np.int64), mc.astype(np.int64), mv.astype(float), n, free, 12)
        return np.sort(np.asarray(f, float))

    def ccx(mesh, name):
        d = tmp / f"modal_{name}"
        d.mkdir(exist_ok=True)
        m.write_ccx_mesh(
            mesh, str(d / "shell_modal.inp"),
            properties={"plate": {"type": "isotropic", "e": m.MAT_STEEL.E, "nu": m.MAT_STEEL.nu,
                                  "rho": m.MAT_STEEL.rho, "thickness": m.t}},
            boundary_nodeset="clamped", solver_type="Modal", num_modes=12,
        )
        m._run_ccx(m._ccx_bin(), d, "shell_modal")
        return np.sort(m._parse_ccx_frequencies(d / "shell_modal.dat", n_modes=12))

    an = anchors()
    print("== modal closed forms (Hz): " + "  ".join(f"{k} {v:.3f}" for k, v in an.items()))
    runs = {}
    for nx, nz in MESHES:
        mesh = _orig_build(nz=nz, nx=nx)
        runs[(nx, nz)] = (aero(mesh), ccx(mesh, f"{nx}x{nz}"))
        a, c = runs[(nx, nz)]
        print(f"  {nx:>2}x{nz:<3} aero " + " ".join(f"{v:8.3f}" for v in a))
        print(f"  {'':6} ccx  " + " ".join(f"{v:8.3f}" for v in c))
    af, cf = runs[MESHES[-1]]
    ap, cp = runs[MESHES[-2]]
    cost = np.abs(af[:, None] - cf[None, :]) / np.maximum(np.abs(cf[None, :]), 1e-14)
    ri, ci = linear_sum_assignment(cost)
    pairs = sorted(zip(ri, ci), key=lambda p: cost[p[0], p[1]])[:5]
    ok = True
    for i, j in pairs:
        da, dc, gap = 100 * (af[i] / ap[i] - 1), 100 * (cf[j] / cp[j] - 1), 100 * cost[i, j]
        name, fa = min(an.items(), key=lambda kv: abs(kv[1] - cf[j]))
        near = abs(fa - cf[j]) / fa <= 0.10
        ea, ec = 100 * (af[i] / fa - 1), 100 * (cf[j] / fa - 1)
        conv = abs(da) < 1 and abs(dc) < 1
        anc = (not near) or (abs(ea) <= 5 and abs(ec) <= 5)
        ok &= conv and gap <= 3 and anc
        tag = f"{name} aero {ea:+.2f}% ccx {ec:+.2f}%" if near else "no closed form within 10%"
        print(
            f"  matched aero#{i} {af[i]:8.3f} ccx#{j} {cf[j]:8.3f} gap {gap:5.2f}% "
            f"8x40->16x80 aero {da:+.3f}% ccx {dc:+.3f}% {'conv' if conv else 'NOT conv'} | {tag}"
        )
    print(f"  => modal comparison {'CLOSED' if ok else 'OPEN'} under C3/C4", flush=True)


with tempfile.TemporaryDirectory(dir=WORK) as td:
    statics(Path(td))
    modal(Path(td))
