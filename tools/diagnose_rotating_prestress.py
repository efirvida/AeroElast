"""Where do the rotating-only local tip modes get their softness? (odd/tasks/rotating-local-tip-modes.md)

At 7.56 rpm the IEA-15 blade (es 0.25) shows local modes (beam fraction 0.12-0.44) from 1.53 Hz,
while parked no local mode exists below 3.8 Hz. K_SP shifts omega^2 by at most Omega^2 = 0.63
rad^2/s^2, so only the centrifugal geometric stiffness K_G can pull a local mode that far down.
This tool rebuilds the exact K_G path of ``assemble_geometric_stiffness`` and reports:

  * prestress: per element the membrane resultant N = h * sigma_m (N/m) in the global frame, its
    spanwise component N_zz and principal values; per span bin the FE section force
    sum(N_zz A)/dz against the exact discrete reference T(z) = sum of the outboard centrifugal
    nodal loads along +z (T2, S2);
  * Rayleigh split of every rated mode stored by ``diagnose_mode_beam_projection.py solve``:
    omega^2 = (phi.K.phi + phi.K_G.phi + phi.K_SP.phi) / phi.M.phi, and the elements that carry
    the most negative K_G energy, with their span position, element set, thickness and N
    (T1, S1).

  python tools/diagnose_rotating_prestress.py solve --modes logs/check/modecls/modes-oml-025.npz --out X.npz
  python tools/diagnose_rotating_prestress.py analyse X.npz [--top 12]
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "tools"))

from run_s4_rotating_modal import AXIS, beam_projection, nodal_translations, omega_from_rpm  # noqa: E402

N_BINS = 200


def element_frame(p):
    """Rust local frame: MITC4 e3 = normalize(n012 + n023), e1 = edge 0-1 orthogonalised; MITC3 e1 = edge 0-1."""
    if len(p) == 4:
        e3 = np.cross(p[1] - p[0], p[2] - p[0]) + np.cross(p[2] - p[0], p[3] - p[0])
    else:
        e3 = np.cross(p[1] - p[0], p[2] - p[0])
    e3 /= np.linalg.norm(e3)
    e1 = p[1] - p[0]
    e1 = e1 - e1.dot(e3) * e3
    e1 /= np.linalg.norm(e1)
    return e1, np.cross(e3, e1), e3


def solve(args):
    import _aeroelast  # noqa: F401  (fail early on a stale extension)
    from petsc4py import PETSc

    from aeroelast.core.assembler import MeshAssembler
    from aeroelast.core.mesh.generators import BladeMesh
    from aeroelast.elements import ElementFamily
    from aeroelast.models.blade.model import build_rust_properties

    ref_modes = np.load(args.modes, allow_pickle=False)
    ref, es = str(ref_modes["ref"]), float(ref_modes["es"])
    i_rpm = int(np.flatnonzero(ref_modes["rpm"] > 0)[0])
    rpm = float(ref_modes["rpm"][i_rpm])
    omega = omega_from_rpm(rpm)

    gen = BladeMesh(element_size=es, yaml_file=str(REPO / "tests/IEA-15-240-RWT.yaml"), skin_reference=ref)
    mesh = gen.generate(renumber="rcm")
    props = build_rust_properties(gen.numad_mesh_data)
    model = {"elements": {"element_family": ElementFamily.SHELL, "span_direction": (0.0, 0.0, 1.0), "properties": props}}
    asm = MeshAssembler(mesh, model)
    xyz = np.asarray([n.coords for n in mesh.nodes])
    if xyz.shape != ref_modes["xyz"].shape or not np.allclose(xyz, ref_modes["xyz"]):
        raise SystemExit(f"mesh differs from {args.modes}: regenerate the modes first")

    n_dof, dpn = asm.dofs_count, asm.dofs_per_node
    nid2idx = mesh.node_id_to_index
    root = {nid2idx[n] for n in mesh.get_node_set("RootNodes").nodes}
    free = np.array([i * dpn + d for i in range(mesh.node_count) if i not in root for d in range(dpn)], dtype=np.int64)
    if not np.array_equal(free, ref_modes[f"free_{i_rpm}"]):
        raise SystemExit("free-DOF set differs from the stored modes")

    # 1. the pre-stress solve of assemble_geometric_stiffness, step for step
    f_cf = asm.assemble_centrifugal_load(omega, AXIS, np.zeros(3))
    kr, kc, kv = (np.asarray(a) for a in asm._rust.assemble_k())
    K = coo_matrix((kv, (kr, kc)), shape=(n_dof, n_dof)).tocsr()
    Kf = K[free][:, free].tocsr()
    Kf.sort_indices()
    K_red = PETSc.Mat().createAIJ(size=Kf.shape, csr=(Kf.indptr.astype(PETSc.IntType),
                                                      Kf.indices.astype(PETSc.IntType), Kf.data))
    K_red.assemble()
    ksp = PETSc.KSP().create()
    ksp.setOperators(K_red)
    ksp.setType("preonly")
    ksp.getPC().setType("lu")
    ksp.getPC().setFactorSolverType("mumps")
    u_red, f_red = K_red.createVecs()
    f_red.setArray(f_cf[free])
    ksp.solve(f_red, u_red)
    u = np.zeros(n_dof)
    u[free] = u_red.getArray()
    print(f"prestress solve: |u|max={np.abs(u).max():.3e} m, |K u - f|/|f| ="
          f" {np.linalg.norm((K @ u - f_cf)[free]) / np.linalg.norm(f_cf[free]):.2e}", flush=True)

    # 2. membrane stress as assemble_geometric_k_from_disp recovers it (centroid, local frame, N/h)
    sig, _ = asm._rust.compute_stress_field(np.ascontiguousarray(u), 0.0, 0)
    sig = np.asarray(sig)[:, [0, 1, 3]]

    elems = asm.elements
    n_el = len(elems)
    set_names = sorted(mesh.element_sets.keys())
    set_of = {}
    for k, name in enumerate(set_names):
        for e in mesh.element_sets[name].elements:
            set_of.setdefault(e.id, k)  # first (alphabetical) set wins; reported for orientation only
    cen, area, h, nzz, nmin, nmax, eset, n_nodes = (np.zeros((n_el, 3)), *(np.zeros(n_el) for _ in range(5)),
                                                     np.full(n_el, -1), np.zeros(n_el, dtype=np.int64))
    elem_dofs = []
    for i, el in enumerate(elems):
        p = np.asarray([n.coords for n in el.nodes], dtype=float)
        n_nodes[i] = len(p)
        cen[i] = p.mean(axis=0)
        area[i] = 0.5 * np.linalg.norm(np.cross(p[2] - p[0], p[3] - p[1])) if len(p) == 4 else \
            0.5 * np.linalg.norm(np.cross(p[1] - p[0], p[2] - p[0]))
        h[i] = float(getattr(el, "thickness", np.nan) or np.nan)
        e1, e2, _ = element_frame(p)
        sxx, syy, sxy = sig[i]
        s3 = sxx * np.outer(e1, e1) + syy * np.outer(e2, e2) + sxy * (np.outer(e1, e2) + np.outer(e2, e1))
        nzz[i] = h[i] * s3[2, 2]
        ev = np.linalg.eigvalsh(np.array([[sxx, sxy], [sxy, syy]])) * h[i]
        nmin[i], nmax[i] = ev
        eset[i] = set_of.get(el.id, -1)
        elem_dofs.append([nid2idx[n.id] * dpn + d for n in el.nodes for d in range(dpn)])

    # 3. section force vs the exact discrete reference
    z = xyz[:, 2]
    z0, span = z.min(), z.max() - z.min()
    edges = np.linspace(z0, z0 + span, N_BINS + 1)
    fz_node = f_cf[np.arange(mesh.node_count) * dpn + 2]
    order = np.argsort(z)
    cum = np.cumsum(fz_node[order][::-1])[::-1]  # outboard sum from each node
    t_ref = np.interp(0.5 * (edges[:-1] + edges[1:]), z[order], cum)
    b = np.clip(np.digitize(cen[:, 2], edges) - 1, 0, N_BINS - 1)
    t_fe = np.bincount(b, weights=nzz * area, minlength=N_BINS) / np.diff(edges)

    # 4. K_G exactly as assemble_geometric_k receives it, one contiguous COO block per element
    gr, gc, gv = (np.asarray(a) for a in asm._rust.assemble_geometric_k(np.ascontiguousarray(sig)))
    blk = np.array([len(d) ** 2 for d in elem_dofs])
    if blk.sum() != len(gv):
        raise SystemExit(f"K_G COO has {len(gv)} entries, expected {blk.sum()}: block split invalid")
    start = np.concatenate([[0], np.cumsum(blk)])
    KG = coo_matrix((gv, (gr, gc)), shape=(n_dof, n_dof)).tocsr()
    mr, mc, mv = (np.asarray(a) for a in asm._rust.assemble_m())
    M = coo_matrix((mv, (mr, mc)), shape=(n_dof, n_dof)).tocsr()
    nodal_m = np.bincount(mr[mr % dpn == 0] // dpn, weights=mv[mr % dpn == 0], minlength=mesh.node_count)
    proj = np.eye(3) - np.outer(AXIS, AXIS)

    freqs, modes = ref_modes[f"freqs_{i_rpm}"], ref_modes[f"modes_{i_rpm}"]
    n_modes = len(freqs)
    e_kg = np.zeros((n_modes, n_el))
    rq = np.zeros((n_modes, 5))  # K, KG, KSP, M, beam fraction
    for k in range(n_modes):
        phi = np.zeros(n_dof)
        phi[free] = modes[k]
        prod = gv * phi[gr] * phi[gc]
        e_kg[k] = np.add.reduceat(prod, start[:-1])
        t = phi.reshape(-1, dpn)[:, :3]
        e_sp = -omega**2 * np.sum(nodal_m * np.einsum("ij,jk,ik->i", t, proj, t))
        frac, _, _ = beam_projection(xyz, nodal_translations(modes[k], free, len(xyz)))
        rq[k] = (phi @ (K @ phi), phi @ (KG @ phi), e_sp, phi @ (M @ phi), frac)

    np.savez_compressed(
        args.out, ref=ref, es=es, rpm=rpm, omega=omega, z0=z0, span=span, cen=cen, area=area, h=h, sig=sig,
        nzz=nzz, nmin=nmin, nmax=nmax, eset=eset, set_names=np.asarray(set_names), n_nodes=n_nodes,
        edges=edges, t_ref=t_ref, t_fe=t_fe, freqs=freqs, e_kg=e_kg, rq=rq)
    print(f"wrote {args.out}", flush=True)
    analyse_file(args.out, args.top)


def analyse_file(path, top):
    d = np.load(path, allow_pickle=False)
    zl = (d["cen"][:, 2] - float(d["z0"])) / float(d["span"])
    names = d["set_names"]
    print(f"\n=== prestress ref={d['ref']} es={float(d['es'])} rpm={float(d['rpm'])} ===")
    mid = 0.5 * (d["edges"][:-1] + d["edges"][1:])
    for j in list(range(0, len(mid), 20)) + list(range(180, len(mid))):
        print(f"SECT z/L={(mid[j] - float(d['z0'])) / float(d['span']):5.3f} T_ref={d['t_ref'][j]:11.4e} N"
              f" T_fe={d['t_fe'][j]:11.4e} N ratio={d['t_fe'][j] / d['t_ref'][j] if d['t_ref'][j] else np.nan:7.3f}")
    for lo, hi in ((0.0, 0.5), (0.5, 0.9), (0.9, 0.94), (0.94, 0.97), (0.97, 1.01)):
        m = (zl >= lo) & (zl < hi)
        print(f"BAND z/L [{lo:.2f},{hi:.2f}) n={int(m.sum()):5d} Nmin p01={np.percentile(d['nmin'][m], 1):11.4e}"
              f" min={d['nmin'][m].min():11.4e} | Nzz p01={np.percentile(d['nzz'][m], 1):11.4e}"
              f" median={np.median(d['nzz'][m]):11.4e} N/m | frac Nmin<0 = {np.mean(d['nmin'][m] < 0):.3f}")
    print("\n=== Rayleigh split, omega^2 = (K + KG + KSP)/M ===")
    for k in np.argsort(d["freqs"]):
        K, KG, KSP, M, frac = d["rq"][k]
        w2 = (2 * np.pi * d["freqs"][k]) ** 2
        print(f"MODE k={int(k):2d} {float(d['freqs'][k]):8.5f} Hz frac={frac:5.3f} w2={w2:9.3f} rq={(K + KG + KSP) / M:9.3f}"
              f" K={K / M:9.3f} KG={KG / M:+9.3f} KSP={KSP / M:+7.3f}")
        e = d["e_kg"][k] / M
        neg = np.argsort(e)[:top]
        share = e[neg].sum() / e[e < 0].sum() if (e < 0).any() else 0.0
        print(f"   top {top} negative K_G elements carry {share:.1%} of the negative K_G energy;"
              f" negative total {e[e < 0].sum():+.3f}, positive total {e[e > 0].sum():+.3f}")
        for i in neg:
            print(f"   el {int(i):6d} z/L={zl[i]:5.3f} E_G={e[i]:+.4f} h={d['h'][i]:.4f} m A={d['area'][i]:.4f} m2"
                  f" nn={int(d['n_nodes'][i])} Nmin={d['nmin'][i]:+.3e} Nmax={d['nmax'][i]:+.3e} Nzz={d['nzz'][i]:+.3e}"
                  f" set={names[d['eset'][i]] if d['eset'][i] >= 0 else '-'}")


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("solve")
    s.add_argument("--modes", type=Path, required=True, help="npz written by diagnose_mode_beam_projection.py solve")
    s.add_argument("--out", type=Path, required=True)
    s.add_argument("--top", type=int, default=12)
    a = sub.add_parser("analyse")
    a.add_argument("npz", type=Path)
    a.add_argument("--top", type=int, default=12)
    args = ap.parse_args()
    if args.cmd == "solve":
        solve(args)
    else:
        analyse_file(args.npz, args.top)


if __name__ == "__main__":
    main()
