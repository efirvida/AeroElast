"""Torque/power sensitivity of the L3 structural reds (odd/tasks/l3-closure.md, T4).

The question each L3 red has to answer before it is declared a limitation: how much does the
measured stiffness error move the rotor torque and power the thesis measures?  This tool answers
it through the **production feedback channels only** -- the two quantities the two-way coupling
hands back to the BEM (``BEMFSIParticipant._compute_deformed_geometry``): the deformed strip
radius and the elastic twist.  Nothing else of the deformation reaches the aerodynamics.

One shell solve at the rated point (V = 10.59 m/s, 7.56 rpm, pitch 0) under the rigid-blade
loads of the production projector -- the same construction as
``tests/validation/rotor/test_bem_fsi_deformed_geometry.py::fsi`` -- gives the deformed
geometry.  Each red is then removed from its own channel and the rigid-geometry BEM is
re-evaluated:

* ``twist``: the elastic twist divided by the measured twist ratio (S-7, ours / beam);
* ``flap``: the radius re-extracted with the flapwise (y) displacement divided by (1 + e_flap);
* ``edge``: the radius re-extracted with the edgewise (x) displacement divided by (1 + e_edge).

It reports thrust, torque and power for the rigid blade, the deformed blade, and each
correction, as a fraction of the deformed value and of the flexible-minus-rigid delta (the
thesis quantity).  One-way and quasi-static: a sensitivity, not a coupled answer.

Run on a compute node (es 0.25 blade solve):
    scripts/aeroenv.sh python -u tools/l3_torque_sensitivity.py --twist-ratio 1.0717 \
        --flap-error 0.0526 --edge-error 0.0801 --out logs/check/l3/sensitivity.json

``--twist-definitions`` (odd/tasks/core-torsion.md, C5) adds, after the same solve, the elastic twist
by three definitions -- (a) production ``_section_rotation``, (b) chord line, (c) box rigid fit --
and the rigid-geometry BEM torque/thrust/power with each (same deformed radius channel).  Rules:

* (b) per strip ring, ``ForceProjector._section_ends`` on the ring points (chord axis of the
  projector) gives the LE node (blunt end) and the TE extreme; the TE point is the mean of the ring
  nodes within 1 % of the chord of that extreme (the two corner nodes of a flat-back TE, one node
  of a sharp TE).  The twist is the change of the angle of ``LE - TE`` about the span axis.
* (c) ring nodes belonging to element sets whose name contains ``SPAR`` or ends with ``_SW``
  (``tools/l3_twist_integral.py --c3``); ``rigid_rotation_fit`` rotation vector projected on the
  span axis.
* (b) and (c) are converted to the BEM twist sense with ``_twist_mesh_to_bem`` exactly as (a);
  a ring without a definition keeps the reference twist.  Units: rad internally, deg reported;
  BEM sense (positive = toward feather for the production frame, see ``_resolve_twist_sense``).
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

V_RATED = 10.59
RPM_RATED = 7.56
PITCH_RATED = 0.0
BEM_CONFIG = {
    "wind_speed": V_RATED,
    "omega": RPM_RATED * 2.0 * np.pi / 60.0,
    "pitch": PITCH_RATED,
    "azimuth": 0.0,
    "air_density": 1.225,
    "dynamic_viscosity": 1.81206e-5,
    "hub_height": 150.0,
    "shear_exp": 0.0,
}


def solve_rated(element_size: float):
    """Mesh, participant and the rated displacement field (n_nodes, 3)."""
    from _aeroelast import PyMeshAssembler  # type: ignore[import-not-found]
    from scipy.sparse import coo_matrix
    from scipy.sparse.linalg import spsolve

    import tests.validation.blade.test_blade_iea15mw_validation as blade_validation
    from aeroelast.core.mesh.entities import MeshElement, Node
    from aeroelast.models.blade.model import Blade
    from aeroelast.solvers.bem.fsi_participant import BEMFSIParticipant
    from tests.support.openfast_bem import build_blade_aero_from_aerodyn
    from tests.support.paths import DATA_DIR

    ad = DATA_DIR / "reference" / "iea15mw_openfast" / "case" / "IEA-15-240-RWT_AeroDyn15.dat"
    Node._id_counter = 0
    MeshElement._id_counter = 0
    model = Blade(str(DATA_DIR / "IEA-15-240-RWT.yaml"), element_size=element_size)
    model.generate_mesh()
    mesh, props = model.mesh, model.get_element_properties()
    assembler = PyMeshAssembler.from_model(
        blade_validation._to_rust_mesh(mesh, props),
        props,
        list(blade_validation.SPAN_DIRECTION),
        None,
    )
    n = assembler.dofs_count
    rows, cols, vals = assembler.assemble_k()
    K = coo_matrix((np.asarray(vals), (np.asarray(rows), np.asarray(cols))), shape=(n, n)).tocsr()
    root = {mesh.node_id_to_index[nid] for nid in mesh.get_node_set("RootNodes").node_ids}
    fixed = {6 * r + d for r in root for d in range(6)}
    free = np.array([i for i in range(n) if i not in fixed], dtype=np.int64)

    participant = BEMFSIParticipant(mesh, build_blade_aero_from_aerodyn(ad), BEM_CONFIG)
    forces, _ = participant._compute_forces(np.zeros((n // 6, 3)))
    load = np.zeros(n)
    for dof in range(3):
        load[dof::6] = forces[:, dof]
    u = np.zeros(n)
    u[free] = np.asarray(spsolve(K[np.ix_(free, free)], load[free]))
    return participant, np.column_stack([u[0::6], u[1::6], u[2::6]])


def evaluate(participant, r, twist) -> dict[str, float]:
    solver, _ = participant._rebuild_bem_solver(r, twist)
    res = solver.compute(V_RATED, RPM_RATED, PITCH_RATED)
    return {"thrust": float(res.thrust), "torque": float(res.torque), "power": float(res.power)}


def ring_rotations(participant, disp):
    """Per-strip section rotation about +span [rad] by definition (a) production, (b) chord, (c) box.

    NaN marks a strip whose ring cannot define the measure (production treats that as zero).
    """
    from aeroelast.postprocess.sectional import rigid_rotation_fit
    from aeroelast.solvers.bem.force_projection import ForceProjector

    mesh = participant._mesh
    span = participant._span_dir
    seed = np.array([1.0, 0.0, 0.0])
    if abs(seed @ span) > 0.9:
        seed = np.array([0.0, 1.0, 0.0])
    e1 = seed - (seed @ span) * span
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(span, e1)
    box = np.zeros(participant._n_nodes, bool)
    for name, eset in mesh.element_sets.items():
        if "SPAR" in name or name.endswith("_SW"):
            for el in eset.elements:
                box[[mesh.node_id_to_index[n] for n in el.node_ids]] = True
    xyz = participant._ref_coords
    n = len(participant._strip_ring_indices)
    out = {key: np.full(n, np.nan) for key in "abc"}
    n_box = np.zeros(n, int)
    for k, ring in enumerate(participant._strip_ring_indices):
        om = participant._section_rotation(disp, ring)
        out["a"][k] = np.nan if om is None else om
        if len(ring) < 3:
            continue
        pts = xyz[ring]
        chord_hat = np.asarray(participant._projector._strip_chord_dirs[k], dtype=float)
        le_i, te_i = ForceProjector._section_ends(pts, chord_hat, span)
        proj = pts @ chord_hat
        chord = proj.max() - proj.min()
        te_set = np.abs(proj - proj[te_i]) <= 0.01 * chord
        vec0 = pts[le_i] - pts[te_set].mean(axis=0)
        vec1 = vec0 + disp[ring][le_i] - disp[ring][te_set].mean(axis=0)
        out["b"][k] = np.arctan2(vec1 @ e2, vec1 @ e1) - np.arctan2(vec0 @ e2, vec0 @ e1)
        sel = ring[box[ring]]
        n_box[k] = len(sel)
        if len(sel) >= 3:
            rot = rigid_rotation_fit(xyz[sel], xyz[sel].mean(axis=0), disp[sel])
            out["c"][k] = float(rot @ span)
    return out, n_box


def twist_definitions(participant, disp, r_def, tw_def):
    """Section (C5): rigid-geometry BEM with the elastic twist of definitions (a), (b), (c)."""
    rots, n_box = ring_rotations(participant, disp)
    sense = participant._twist_mesh_to_bem
    ref = participant._ref_twist
    twists = {k: ref + sense * np.nan_to_num(v) for k, v in rots.items()}
    check = float(np.max(np.abs(twists["a"] - tw_def)))
    cases = {k: evaluate(participant, r_def, tw) for k, tw in twists.items()}
    r_ref = participant._ref_r
    eta = r_ref / r_ref.max()
    deg = {k: np.rad2deg(tw - ref) for k, tw in twists.items()}
    picks = {f"{e}": int(np.argmin(np.abs(eta - e))) for e in (0.7, 0.85, 0.95)} | {"tip": len(eta) - 1}
    base = cases["a"]
    report = {
        "twist_sense": float(sense),
        "max_abs_a_vs_production_rad": check,
        "n_strips_without": {k: int(np.isnan(v).sum()) for k, v in rots.items()},
        "min_box_ring_nodes": int(n_box[n_box > 0].min()) if (n_box > 0).any() else 0,
        "cases": cases,
        "relative_to_a": {
            k: {q: c[q] / base[q] - 1.0 for q in ("thrust", "torque", "power")} for k, c in cases.items()
        },
        "elastic_twist_deg": {
            k: {name: float(deg[k][i]) for name, i in picks.items()} for k in deg
        },
        "strips": {
            "r": [float(v) for v in r_ref],
            "twist_deg": {k: [float(v) for v in deg[k]] for k in deg},
        },
    }
    print(f"twist definitions (BEM sense, deg; a-production b-chord c-box); a vs production max {check:.2e} rad")
    print(f"{'def':4s} {'thrust':>12s} {'torque':>12s} {'power':>12s} {'T vs a':>9s} {'Q vs a':>9s} "
          f"{'P vs a':>9s} " + " ".join(f"{'tw@' + n:>9s}" for n in picks))
    for k, c in cases.items():
        rel = report["relative_to_a"][k]
        print(f"{k:4s} {c['thrust']:12.1f} {c['torque']:12.1f} {c['power']:12.1f} {rel['thrust']:+9.3%} "
              f"{rel['torque']:+9.3%} {rel['power']:+9.3%} "
              + " ".join(f"{report['elastic_twist_deg'][k][n]:+9.3f}" for n in picks))
    sel = np.linspace(0, len(r_ref) - 1, 16).astype(int)
    print("profile r[m] | b-a, c-a [deg] | a")
    for i in sel:
        print(f"  {r_ref[i]:7.2f} | {deg['b'][i] - deg['a'][i]:+8.3f} {deg['c'][i] - deg['a'][i]:+8.3f} "
              f"| {deg['a'][i]:+8.3f}")
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--element-size", type=float, default=0.25)
    parser.add_argument("--twist-ratio", type=float, required=True, help="S-7 twist, ours / beam")
    parser.add_argument("--flap-error", type=float, required=True, help="S-2 flap tip, ours/ref - 1")
    parser.add_argument("--edge-error", type=float, required=True, help="edgewise, ours/ref - 1")
    parser.add_argument("--twist-definitions", action="store_true",
                        help="C5: also compare the production, chord-line and box twist definitions")
    parser.add_argument("--save-displacements", type=Path)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    participant, disp = solve_rated(args.element_size)
    if args.save_displacements:
        np.save(args.save_displacements, disp)
    r_ref, tw_ref = participant._ref_r, participant._ref_twist
    r_def, tw_def = participant._compute_deformed_geometry(disp)
    d_twist = tw_def - tw_ref

    flap_fixed = disp.copy()
    flap_fixed[:, 1] /= 1.0 + args.flap_error
    r_flap, _ = participant._compute_deformed_geometry(flap_fixed)
    edge_fixed = disp.copy()
    edge_fixed[:, 0] /= 1.0 + args.edge_error
    r_edge, _ = participant._compute_deformed_geometry(edge_fixed)

    cases = {
        "rigid": evaluate(participant, r_ref, tw_ref),
        "deformed": evaluate(participant, r_def, tw_def),
        "twist_only": evaluate(participant, r_ref, tw_def),
        "radius_only": evaluate(participant, r_def, tw_ref),
        "twist_corrected": evaluate(participant, r_def, tw_ref + d_twist / args.twist_ratio),
        "flap_corrected": evaluate(participant, r_flap, tw_def),
        "edge_corrected": evaluate(participant, r_edge, tw_def),
    }
    rigid, deformed = cases["rigid"], cases["deformed"]
    report = {
        "inputs": {k: v for k, v in vars(args).items() if k != "twist_definitions"} | {"save_displacements": str(args.save_displacements),
                                "out": str(args.out)},
        "geometry": {
            "tip_twist_deg": float(np.rad2deg(d_twist[-1])),
            "max_abs_twist_deg": float(np.rad2deg(np.max(np.abs(d_twist)))),
            "max_dr_m": float(np.max(r_def - r_ref)),
            "min_dr_m": float(np.min(r_def - r_ref)),
            "tip_disp_m": [float(v) for v in disp[np.argmax(np.linalg.norm(disp, axis=1))]],
        },
        "cases": cases,
        "relative": {},
    }
    for name, case in cases.items():
        rel = {}
        for q in ("thrust", "torque", "power"):
            delta_flex = deformed[q] - rigid[q]
            rel[q] = {
                "vs_deformed": case[q] / deformed[q] - 1.0,
                "vs_rigid": case[q] / rigid[q] - 1.0,
                "fraction_of_flex_delta": (case[q] - deformed[q]) / delta_flex
                if abs(delta_flex) > 0 else None,
            }
        report["relative"][name] = rel
    if args.twist_definitions:
        report["twist_definitions"] = twist_definitions(participant, disp, r_def, tw_def)
    args.out.write_text(json.dumps(report, indent=2))

    print(f"tip elastic twist {report['geometry']['tip_twist_deg']:+.3f} deg, "
          f"dr {report['geometry']['min_dr_m']:+.3f}..{report['geometry']['max_dr_m']:+.3f} m")
    print(f"{'case':16s} {'T vs def':>9s} {'Q vs def':>9s} {'P vs def':>9s} "
          f"{'Q vs rig':>9s} {'dQ/dQflex':>10s}")
    for name, rel in report["relative"].items():
        frac = rel["torque"]["fraction_of_flex_delta"]
        print(f"{name:16s} {rel['thrust']['vs_deformed']:+9.3%} {rel['torque']['vs_deformed']:+9.3%} "
              f"{rel['power']['vs_deformed']:+9.3%} {rel['torque']['vs_rigid']:+9.3%} "
              f"{frac if frac is None else format(frac, '+10.3f')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
