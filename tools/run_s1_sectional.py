"""S-1 — sectional property comparison: AeroElast shell vs official beam targets.

Input-chain verification (2026-09-10): the WindIO yaml publishes BOTH
structural representations side by side —
``elastic_properties_mb.six_x_six`` (the beam 6x6 matrices, used verbatim
by the OpenFAST BeamDyn blade file — identical to machine precision) and
``internal_structure_2d_fem`` (the 2D layup the shell is built from).  The
observed shell-vs-beam bands (GJ ~+17%, EI 10-25%, mass +4-8%) are the
disagreement between the yaml's OWN two representations (PreComp/BECAS 6x6
vs the NuMAD 2D layup), not an extraction or solver error: the shell is
CLT-consistent with its layup (cross-checked against an independent
Bredt-Batho closed-section estimate).  The comparisons below are therefore
a consistency check between the two published representations.

Extracts EA/GJ/EI_flap/EI_edge/mass-per-length from a shell blade mesh with
the SectionalExtractor and compares against the official reference tables
(BeamDyn 6x6 section matrices + ElastoDyn distributed properties), which are
the PreComp/VABS/BECAS results published for the IEA 15 MW.

Usage (SDumont):
  module load glu gcc/14.2.0_sequana
  python tools/run_s1_sectional.py [--blade-yaml tests/IEA-15-240-RWT.yaml]
                                   [--blade-xlsx tests/NuMAD_utd_iea15mw.xlsx]
                                   [--element-size 0.5] [--n-slices 26]

Outputs:
  docs/validation_data/generated/s1_sectional_properties.csv
  prints a spanwise comparison table.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
from scipy.interpolate import UnivariateSpline  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(REPO_ROOT / "tools"))

from tools.beam_reference import load_beamdyn_blade, load_elastodyn_blade  # noqa: E402

REF_DIR = REPO_ROOT / "tests" / "IEA15MW" / "reference"
ELASTODYN_FILE = REF_DIR / "IEA-15-240-RWT_ElastoDyn_blade.dat"
BEAMDYN_FILE = REF_DIR / "IEA-15-240-RWT_BeamDyn_blade.dat"


def _interp(x_src, y_src, x_dst):
    return np.interp(x_dst, x_src, y_src)


def ea_gj_from_static_response(mesh, props, load_ax=1e6, torque=1e6,
                                numad_data=None):
    """Extract EA(r) and GJ(r) from full-mesh static tip-load solves.

    GJ: tip torque about the span (force couple on the tip section) ->
        section twist theta(r) from a best-fit rigid rotation ->
        twist rate -> GJ = T/kappa_t.

    EA (response): tip axial load -> axial strain eps(r) = d/dz of the
        tangent-projected section-mean displacement -> EA = F/eps.  NOTE:
        on the pre-bent blade this channel is contaminated by the bending
        projection and is reported for reference only.

    EA (CLT): the direct laminate section integral — EA(z) = sum over the
        station's elements of A11(layup) * width_span — immune to the
        pre-bend contamination (the same construction as the Bredt GJ).

    NOTE: the EA extraction currently shows a systematic bias outboard
    (shell ~0.1 GN vs BeamDyn ~3.9 GN at the tip) — documented as an open
    item; use the GJ result only for assertions.
    """
    from aeroelast.core.bc import DirichletCondition, NodalLoad
    from aeroelast.elements import ElementFamily
    from aeroelast.solvers.elasticity.static_linear import StaticLinearSolver

    from aeroelast.postprocess.sectional import SectionalExtractor

    model = {
        "solver": {},
        "elements": {
            "element_family": ElementFamily.SHELL,
            "span_direction": (0.0, 0.0, 1.0),
            "properties": props,
        },
    }
    ext = SectionalExtractor(mesh, props)
    membership = ext.station_membership()
    sts = sorted(membership)
    zs = np.array([membership[st]["z"] for st in sts])
    tip_nodes = membership[sts[-1]]["nodes"]
    coords = np.asarray([n.coords for n in mesh.nodes])
    dpn = mesh.dofs_per_node if hasattr(mesh, "dofs_per_node") else 6

    root = sorted(mesh.get_node_set("RootNodes").nodes.keys())
    idx = mesh.node_id_to_index
    rootd = sorted(idx[n] * dpn + d for n in root for d in range(dpn))

    def solve_case(tip_dof_vals):
        s = StaticLinearSolver(mesh, model)
        s.add_dirichlet_conditions([DirichletCondition(rootd, 0.0)])
        for t, vals in tip_dof_vals.items():
            s.add_nodal_loads([NodalLoad([int(t) * dpn + d], [v])
                               for d, v in enumerate(vals)])
        return s.solve()

    # EA: tip axial load, displacements projected onto each station's local
    # tangent (the pre-bent reference line) — removes the bending-induced
    # axial contamination that biased the raw global-u_z derivative.
    u_ax = solve_case({t: [0.0, 0.0, load_ax / len(tip_nodes)] for t in tip_nodes})
    # station tangents from the centroid line of the pre-bent geometry
    cents = np.array([coords[membership[st]["nodes"]].mean(axis=0) for st in sts])
    tangs = np.gradient(cents, zs, axis=0)
    tangs /= np.linalg.norm(tangs, axis=1)[:, None]
    u_mean = np.array([
        [u_ax[membership[st]["nodes"] * dpn + d].mean() for d in range(3)]
        for st in sts])
    u_t = (u_mean * tangs).sum(axis=1)
    sp_u = UnivariateSpline(zs, u_t, s=1e-1)
    eps_z = sp_u.derivative(1)(zs)
    ea_static = np.where(np.abs(eps_z) > 1e-12, load_ax / np.abs(eps_z), np.nan)

    # GJ: tip torque couple about the LOCAL (pre-bent) span tangent.
    # A global-z couple on the swept blade carries a bending component that
    # contaminates the measured twist rate (same artefact as EI_flap/EA);
    # applying the pure couple about the section tangent removes it.
    tip_xyz = coords[tip_nodes]
    tip_centroid = tip_xyz.mean(axis=0)
    r_tip = tip_xyz - tip_centroid
    t_tip = tangs[-1]
    r_perp = r_tip - np.outer(r_tip @ t_tip, t_tip)
    a_couple = torque / max(np.sum(r_perp**2), 1e-12)
    f_tip = a_couple * np.cross(np.broadcast_to(t_tip, r_perp.shape), r_perp)
    tip_f = {int(t): [float(f[0]), float(f[1]), float(f[2])]
             for t, f in zip(tip_nodes, f_tip)}
    u_t = solve_case(tip_f)
    theta = []
    for i, st in enumerate(sts):
        nodes = membership[st]["nodes"]
        xyz = coords[nodes]
        r = xyz - xyz.mean(axis=0)
        t_s = tangs[i]
        r_p = r - np.outer(r @ t_s, t_s)
        u_s = np.array([[u_t[n * dpn + d] for d in range(3)] for n in nodes])
        u_p = u_s - np.outer(u_s @ t_s, t_s)
        num = float(np.sum(np.cross(r_p, u_p) @ t_s))
        den = float(np.sum(r_p**2))
        theta.append(num / max(den, 1e-12))
    theta = np.asarray(theta)
    sp_t = UnivariateSpline(zs, theta, s=1e-3)
    twist_rate = np.abs(sp_t.derivative(1)(zs))
    gj_static = np.where(twist_rate > 1e-12, torque / twist_rate, np.nan)

    # EA via the direct CLT section integral (immune to the pre-bend
    # bending contamination that biases the global-response extraction):
    #   EA(z) = sum over the station's elements of A11(layup) * width_span
    # where A11 = sum Qbar11_k * h_k (the laminate's extensional stiffness)
    # and width_span = element area / station slice length.
    if numad_data is None:
        ea_clt = np.full(len(sts), np.nan)
        return zs, ea_static, gj_static, ea_clt

    mats = {m["name"]: m for m in numad_data["materials"]}
    sections = {s2["elementSet"]: s2 for s2 in numad_data["sections"]}

    def layup_a11(layup):
        a11 = 0.0
        for mat_name, h, ang in layup:
            m = mats[mat_name]
            el = m["elastic"]
            E_raw = el["E"]
            if isinstance(E_raw, (list, tuple)):
                E1, E2 = float(E_raw[0]), float(E_raw[1])
                nu12 = float(el["nu"][0])
                G12 = float(el["G"][0]) if el.get("G") else E1 / (2 * (1 + nu12))
            else:
                E1 = E2 = float(E_raw)
                nu12 = float(el["nu"])
                G12 = float(el["G"]) if el.get("G") else E1 / (2 * (1 + nu12))
            th = np.deg2rad(ang)
            s2, c2 = np.sin(th) ** 2, np.cos(th) ** 2
            s4, c4 = s2 * s2, c2 * c2
            Q11 = E1 / (1 - nu12**2)
            Q22 = E2 / (1 - nu12**2)
            Q12 = nu12 * Q22
            Qbar11 = Q11 * c4 + 2 * (Q12 + 2 * G12) * s2 * c2 + Q22 * s4
            a11 += Qbar11 * float(h)
        return a11

    # element -> set name (the DETAILED set wins; the aggregate sets like
    # allOuterShellEls come later in the iteration and must not override)
    elem_to_set = {}
    for ename, eset in mesh.element_sets.items():
        for e in eset.elements:
            elem_to_set.setdefault(e.id, ename)
    elem_set_a11 = {}
    for eid, ename in elem_to_set.items():
        if ename in sections:
            elem_set_a11[eid] = layup_a11(sections[ename]["layup"])

    nid_to_idx = mesh.node_id_to_index
    elem_areas = np.zeros(len(mesh.elements))
    for i, e in enumerate(mesh.elements):
        c = coords[[nid_to_idx[n] for n in e.node_ids]]
        elem_areas[i] = (0.5 * np.linalg.norm(np.cross(c[2] - c[0], c[3] - c[1]))
                         if len(c) == 4 else
                         0.5 * np.linalg.norm(np.cross(c[1] - c[0], c[2] - c[0])))

    node_set = {st: set(membership[st]["nodes"].tolist()) for st in sts}
    ea_clt = np.zeros(len(sts))
    for i, st in enumerate(sts):
        zprev = membership[st - 1]["z"] if st - 1 in membership else membership[st]["z"]
        dz = max(membership[st]["z"] - zprev, 1e-9)
        total = 0.0
        ns = node_set[st]
        for e in mesh.elements:
            a11 = elem_set_a11.get(e.id)
            if a11 is None or any(nid_to_idx[nid] not in ns for nid in e.node_ids):
                continue
            total += a11 * elem_areas[e.id]
        ea_clt[i] = total / dz
    return zs, ea_static, gj_static, ea_clt


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--blade-yaml", default=str(REPO_ROOT / "tests" / "IEA-15-240-RWT.yaml"))
    ap.add_argument("--blade-xlsx", default=None)
    ap.add_argument("--element-size", type=float, default=0.5)
    ap.add_argument("--n-slices", type=int, default=26)
    ap.add_argument("--airfoil-dir", default=str(REPO_ROOT / "tests" / "airfoils"))
    args = ap.parse_args()

    from aeroelast.core.mesh.generators import BladeMesh
    from aeroelast.models.blade.model import build_rust_properties
    from aeroelast.postprocess.sectional import SectionalExtractor

    ed = load_elastodyn_blade(ELASTODYN_FILE)
    bd = load_beamdyn_blade(BEAMDYN_FILE)

    generator = BladeMesh(
        yaml_file=args.blade_yaml if not args.blade_xlsx else None,
        excel_file=args.blade_xlsx,
        airfoil_dir=args.airfoil_dir,
        element_size=args.element_size,
    )
    mesh = generator.generate(renumber="rcm")
    props = build_rust_properties(generator.numad_mesh_data)
    ext = SectionalExtractor(mesh, props)
    print(f"[S1] extracting EI(r) via static tip-load curvature...")
    sec = ext.extract(n_slices=args.n_slices)
    z_static, ei_flap_static, ei_edge_static = ext.ei_from_static_response(props)

    ei_flap_bd = _interp(bd.z, bd.ei_flap, z_static)
    ei_edge_bd = _interp(bd.z, bd.ei_edge, z_static)
    ei_flap_ed = _interp(ed.z, ed.ei_flap, z_static)
    ei_edge_ed = _interp(ed.z, ed.ei_edge, z_static)
    m_bd = _interp(bd.z, bd.mass_dens, sec.z)
    m_shell_at_static = _interp(sec.z, sec.mass_dens, z_static)

    print(f"\n{'z [m]':>7s} | {'EI_flap shell/BD/ED [GNm²]':>26s} | "
          f"{'EI_edge shell/BD/ED [GNm²]':>26s} | {'m shell/BD [kg/m]':>16s}")
    rows = []
    for i, z in enumerate(z_static):
        print(f"{z:7.1f} | {ei_flap_static[i]/1e9:8.2f} {ei_flap_bd[i]/1e9:8.2f} "
              f"{ei_flap_ed[i]/1e9:8.2f} | {ei_edge_static[i]/1e9:8.2f} "
              f"{ei_edge_bd[i]/1e9:8.2f} {ei_edge_ed[i]/1e9:8.2f} | "
              f"{m_shell_at_static[i]:7.1f} {_interp(bd.z, bd.mass_dens, [z])[0]:7.1f}")
        rows.append((z, ei_flap_static[i], ei_edge_static[i],
                     ei_flap_bd[i], ei_edge_bd[i]))

    mid = (z_static > 0.2 * z_static.max()) & (z_static < 0.85 * z_static.max())
    if mid.any():
        print(f"\nmid-span mean |rel err| (0.2–0.85R): "
              f"EI_flap={np.nanmean(np.abs(ei_flap_static[mid]/ei_flap_bd[mid]-1)):.2%}  "
              f"EI_edge={np.nanmean(np.abs(ei_edge_static[mid]/ei_edge_bd[mid]-1)):.2%}")

    # ------------------------------------------------------------------
    # V-05 closure: EA + GJ spanwise vs BeamDyn
    # ------------------------------------------------------------------
    zs, ea_static, gj_static, ea_clt = ea_gj_from_static_response(
        mesh, props, numad_data=generator.numad_mesh_data)

    ea_bd = _interp(bd.z, bd.ea, zs)
    gj_bd = _interp(bd.z, bd.gj, zs)
    print(f"\n{'z [m]':>7s} | {'EA_CLT/BD [GN]':>15s} | {'GJ shell/BD [GNm²]':>18s}")
    for i in range(len(zs)):
        print(f"{zs[i]:7.1f} | {ea_clt[i]/1e9:8.2f} {ea_bd[i]/1e9:8.2f} | "
              f"{gj_static[i]/1e9:8.2f} {gj_bd[i]/1e9:8.2f}")
    mid2 = (zs > 0.2 * zs.max()) & (zs < 0.85 * zs.max())
    if mid2.any():
        ea_err = np.nanmean(np.abs(ea_clt[mid2] / ea_bd[mid2] - 1))
        gj_err = np.nanmean(np.abs(gj_static[mid2] / gj_bd[mid2] - 1))
        print(f"\nmid-span mean |rel err| (0.2–0.85R): EA_CLT={ea_err:.2%}  GJ={gj_err:.2%}")

    out_csv = REPO_ROOT / "docs" / "validation_data" / "generated" / "s1_sectional_properties.csv"
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with open(out_csv, "w") as f:
        f.write("z_m,EA_CLT_shell,EA_beamdyn,GJ_shell,GJ_beamdyn,EI_flap_shell,"
                "EI_flap_beamdyn,EI_edge_shell,EI_edge_beamdyn,m_shell,m_beamdyn\n")
        for i, z in enumerate(zs):
            f.write(f"{z:.6e},{ea_clt[i]:.6e},{ea_bd[i]:.6e},{gj_static[i]:.6e},"
                    f"{gj_bd[i]:.6e},"
                    f"{_interp(z_static, ei_flap_static, [z])[0]:.6e},"
                    f"{_interp(bd.z, bd.ei_flap, [z])[0]:.6e},"
                    f"{_interp(z_static, ei_edge_static, [z])[0]:.6e},"
                    f"{_interp(bd.z, bd.ei_edge, [z])[0]:.6e},"
                    f"{_interp(sec.z, sec.mass_dens, [z])[0]:.6e},"
                    f"{_interp(bd.z, bd.mass_dens, [z])[0]:.6e}\n")
    print(f"\n[S1] wrote {out_csv}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
