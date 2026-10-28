#!/usr/bin/env python
"""Box rigid fit and extractor-free energy check of the S-1 / S-7 reads (odd/tasks/l3-closure.md, T15 + T16).

Blade mode (default). Builds the mesh as the S-1 test does (``BladeMesh(skin_reference="mid")``, ``renumber="rcm"``,
``build_rust_properties``, stations = ``native_bay_centres``), runs the S-1 solves (tip edge +X and tip flap +Y of
``ei_from_static_response``; pure torque of ``ea_gj_from_static_response``) once and reports

* EI_edge, EI_flap and GJ, mid-span 0.2-0.85 R, signed mean and mean |rel| vs the deck (EI: elastic-centre principal
  values as ``tests/test_iea15mw_s1_sectional.py``), rigid fit on all station nodes and on the box nodes
  (``fit_nodes``: caps SPAR* + webs *_SW) side by side, with the count of NaN stations;
* the S-7 band ratio (0.3-0.9 of the blade length, 1000 N m, ``run_torsion_case`` + ``section_twists_deg`` vs
  ``beam_twist_profile_deg``) all vs box, on the S-1 mesh (the S-7 test builds its own sampling, so the absolute value
  is not the test's; the all -> box change is on one mesh);
* W_FE = 1/2 f.u (full DOF vectors) for the three load cases and W_beam = 1/2 int r(s)^T C(s) r(s) ds with C = K^-1 of
  the deck 6x6 (interpolated linearly in eta, as BeamDyn does), r(s) = [F; M(s)] resolved in the deck frame at s, and
  the decoupled form (M_x^2/K44 + M_y^2/K55 + M_z^2/K66). W_FE / W_beam per case and per tip cut eta_max.

Deck frame (pinned from the repository, ``tools/bd_static_reference.py`` header and the S-1 tangent frame): BeamDyn axes
(x, y, z) map to the shell axes by (x, y, z) -> (-y, x, z); z is the span, the section frame is the BeamDyn global
triad rotated about z by the key-point twist (positive twist rotates chord y toward -x), then tilted to the key-point
tangent; K44 is bending about the section x axis (edgewise, EI_edge) and K55 flapwise. The twist sign is verified
against the chord axis of the mesh sections (printed). The reference axis is the key-point line of the BeamDyn main file.

The tip is cut: the deck tail collapses (IEA-3.4 GJ 4e5 -> 2e3, IEA-15 7e6 -> 7e4 over the last 5 %), and a constant
torque makes int ds / GJ tip-dominated; W_beam is printed for several eta_max, the S-7 cut (0.9504) is the reference.

Other modes:  --mode tube       closed forms (uniform thin-walled tube, tip torque and tip force);
              --mode prismatic  GJ all vs box on a prismatic extrusion vs its mesh Bredt (core-torsion C4 / P1);
              --mode selftest   the beam-energy integral on a synthetic deck.

  python tools/l3_box_energy.py --yaml tests/IEA-15-240-RWT.yaml \\
      --deck tests/IEA15MW/reference/IEA-15-240-RWT_BeamDyn_blade.dat --element-size 0.25 --out logs/check/l3/t15/iea15.json
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(REPO_ROOT / "tools"))

from l3_bredt_yaml import interp_deck, load_beamdyn_6x6  # noqa: E402

ETA_CUTS = (0.85, 0.90, 0.95, 1.0)
S7_ETA_CUT = 0.9504  # 111.2 m of the 117 m IEA-15 blade: the last valid BeamDyn station of run_s7_torsion
S7_MOMENT = 1000.0
#: BD (x, y, z) -> shell (-y, x, z)  (tools/bd_static_reference.py): shell = BD_TO_SHELL @ bd
BD_TO_SHELL = np.array([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])


# =============================================================================================
# deck reference axis, frame, energy
# =============================================================================================
def load_keypoints(path) -> np.ndarray:
    """(n, 4) key points [x, y, z, twist deg] of a BeamDyn main file (BD axes)."""
    rows, on = [], False
    for ln in open(path):
        if "kp_xr" in ln:
            on = True
            continue
        if not on:
            continue
        tok = ln.split()
        if len(tok) == 4:
            try:
                rows.append([float(t) for t in tok])
            except ValueError:
                continue
        elif rows:
            break
    if len(rows) < 3:
        raise ValueError(f"no key points in {path}")
    return np.asarray(rows)


def reference_axis(kp, n=20001) -> dict:
    """Dense reference axis of the key points: cubic spline in chord length; arc length s, eta = s / S, tangent, twist."""
    from scipy.interpolate import CubicSpline

    t_kp = np.concatenate([[0.0], np.cumsum(np.linalg.norm(np.diff(kp[:, :3], axis=0), axis=1))])
    sp = CubicSpline(t_kp, kp[:, :3], axis=0)
    t = np.linspace(0.0, t_kp[-1], n)
    pos, d = sp(t), sp(t, 1)
    speed = np.linalg.norm(d, axis=1)
    s = np.concatenate([[0.0], np.cumsum(0.5 * (speed[1:] + speed[:-1]) * np.diff(t))])
    return {"s": s, "eta": s / s[-1], "pos": pos, "tangent": d / speed[:, None], "twist_deg": np.interp(t, t_kp, kp[:, 3]),
            "length": float(s[-1])}


def deck_triads(tangent, twist_deg, sign=1.0) -> np.ndarray:
    """(n, 3, 3) columns (x_d, y_d, z_d) in BD axes: Rz(sign * twist), then tilted so that z follows the tangent."""
    th = np.deg2rad(sign * np.asarray(twist_deg))
    c, s_ = np.cos(th), np.sin(th)
    n = len(th)
    base = np.zeros((n, 3, 3))
    base[:, 0, 0], base[:, 1, 0] = c, s_  # x_d
    base[:, 0, 1], base[:, 1, 1] = -s_, c  # y_d
    base[:, 2, 2] = 1.0  # z_d
    ez = np.array([0.0, 0.0, 1.0])
    ax = np.cross(ez, tangent)
    sin_a = np.linalg.norm(ax, axis=1)
    cos_a = tangent[:, 2]
    k = np.where(sin_a[:, None] > 1e-14, ax / np.maximum(sin_a, 1e-300)[:, None], 0.0)
    kx = np.zeros((n, 3, 3))
    kx[:, 0, 1], kx[:, 0, 2] = -k[:, 2], k[:, 1]
    kx[:, 1, 0], kx[:, 1, 2] = k[:, 2], -k[:, 0]
    kx[:, 2, 0], kx[:, 2, 1] = -k[:, 1], k[:, 0]
    eye = np.eye(3)[None]
    rot = eye + sin_a[:, None, None] * kx + (1.0 - cos_a)[:, None, None] * (kx @ kx)  # Rodrigues, ez -> tangent
    return rot @ base


def distributed_nodal_loads(coords, elements, exclude, comp, total):
    """Span-distributed transverse nodal forces: ``f_i = total * w_i / sum w`` along shell axis ``comp`` (0 = +X edge, 1 = +Y flap).

    Weight ``w_i`` = the node's tributary span length: half the mean span extent (z max - z min) of its elements on the
    root side plus half that of its elements on the tip side (a node with elements on one side only gets half of that
    side's mean). Applied to every node except ``exclude`` (the root ring). ``elements`` is a list of node-index arrays.
    Returns ``(nodes, forces)`` with ``forces`` (n, 3).
    """
    n = len(coords)
    sums = np.zeros((2, n))
    cnt = np.zeros((2, n))
    for ids in elements:
        ids = np.asarray(ids)
        dz = float(coords[ids, 2].max() - coords[ids, 2].min())
        zm = float(coords[ids, 2].mean())
        for i in ids:
            side = 0 if coords[i, 2] < zm else 1  # element lies above the node -> tip side
            sums[side, i] += dz
            cnt[side, i] += 1
    w = 0.5 * (np.divide(sums[0], cnt[0], out=np.zeros(n), where=cnt[0] > 0)
               + np.divide(sums[1], cnt[1], out=np.zeros(n), where=cnt[1] > 0))
    w[np.asarray(list(exclude), dtype=int)] = 0.0
    nodes = np.nonzero(w > 0.0)[0]
    forces = np.zeros((len(nodes), 3))
    forces[:, comp] = total * w[nodes] / w[nodes].sum()
    return nodes, forces


def beam_energy(axis, etas, ks, load_xyz, load_f, sign=1.0, eta_max=1.0, n=4000, distributed=False) -> dict:
    """W_beam = 1/2 int r^T C r ds over s in [0, eta_max S] for the nodal loads (shell axes, linear: undeformed geometry).

    ``r(s) = [F(s); M(s)]`` resolved in the deck frame at s, ``C = K^-1``, K linear in eta between the deck stations.
    ``distributed=False`` (tip loads): every load acts on every section, ``F = sum f``, ``M = sum (x_i - x_ref(s)) x f_i``.
    ``distributed=True``: only the loads outboard of the section (node z >= the reference axis z at s) act on it, so
    ``F(s) = sum_{z_i>=z(s)} f_i`` and ``M(s) = sum_{z_i>=z(s)} (x_i - x_ref(s)) x f_i``: the beam carries exactly the applied
    nodal forces (the discrete form of V = int q, M = int (r' - r) x q). Returns the full 6x6 form and the decoupled
    ``int (M_x^2/K44 + M_y^2/K55 + M_z^2/K66) ds / 2``.
    """
    s_max = eta_max * axis["length"]
    s = np.linspace(0.0, s_max, n + 1)
    eta = s / axis["length"]
    pos = np.column_stack([np.interp(s, axis["s"], axis["pos"][:, d]) for d in range(3)])
    tan = np.column_stack([np.interp(s, axis["s"], axis["tangent"][:, d]) for d in range(3)])
    tan /= np.linalg.norm(tan, axis=1, keepdims=True)
    tw = np.interp(s, axis["s"], axis["twist_deg"])
    rot = deck_triads(tan, tw, sign)  # columns x_d y_d z_d in BD axes
    ref_shell = pos @ BD_TO_SHELL.T
    xyz, f = np.asarray(load_xyz, float), np.asarray(load_f, float)
    if distributed:  # suffix sums over the loads sorted by z: loads with z_i >= z_ref(s)
        order = np.argsort(xyz[:, 2])
        xyz, f = xyz[order], f[order]
        f_suf = np.concatenate([np.cumsum(f[::-1], axis=0)[::-1], np.zeros((1, 3))])
        xf_suf = np.concatenate([np.cumsum(np.cross(xyz, f)[::-1], axis=0)[::-1], np.zeros((1, 3))])
        k0 = np.searchsorted(xyz[:, 2], ref_shell[:, 2] - 1e-9, side="left")
        f_s, moment = f_suf[k0], xf_suf[k0] - np.cross(ref_shell, f_suf[k0])
    else:
        f_s = np.tile(f.sum(axis=0), (len(s), 1))
        moment = np.zeros((len(s), 3))
        for xi, fi in zip(xyz, f, strict=True):
            moment += np.cross(xi[None, :] - ref_shell, fi[None, :])
    f_bd = f_s @ BD_TO_SHELL
    m_bd = moment @ BD_TO_SHELL
    r = np.concatenate([np.einsum("nji,nj->ni", rot, f_bd), np.einsum("nji,nj->ni", rot, m_bd)], axis=1)
    k = np.array([interp_deck(etas, ks, e) for e in eta])
    x = np.linalg.solve(k, r[:, :, None])[:, :, 0]
    full = 0.5 * np.trapezoid(np.einsum("ni,ni->n", r, x), s)
    dec = 0.5 * np.trapezoid(r[:, 3] ** 2 / k[:, 3, 3] + r[:, 4] ** 2 / k[:, 4, 4] + r[:, 5] ** 2 / k[:, 5, 5], s)
    return {"eta_max": eta_max, "W_full": float(full), "W_decoupled": float(dec)}


def solve_nodal(mesh, props, nodes, forces) -> np.ndarray:
    """Static solve, root ring clamped, nodal translation forces ``forces`` (n, 3) on ``nodes``; full solution vector."""
    from aeroelast.core.bc import DirichletCondition, NodalLoad
    from aeroelast.elements import ElementFamily
    from aeroelast.solvers.elasticity.static_linear import StaticLinearSolver

    dpn = 6
    model = {"solver": {}, "elements": {"element_family": ElementFamily.SHELL, "span_direction": (0.0, 0.0, 1.0),
                                        "properties": props}}
    idx = mesh.node_id_to_index
    root = sorted(mesh.get_node_set("RootNodes").nodes.keys())
    sol = StaticLinearSolver(mesh, model)
    sol.add_dirichlet_conditions([DirichletCondition(sorted(idx[n] * dpn + d for n in root for d in range(dpn)), 0.0)])
    dofs = [int(t) * dpn + d for t in nodes for d in range(3)]
    vals = [float(v) for fv in forces for v in fv]
    sol.add_nodal_loads([NodalLoad(dofs, vals)])
    return np.asarray(sol.solve(), float)


def mesh_element_nodes(mesh):
    idx = mesh.node_id_to_index
    return [np.array([idx[n] for n in e.node_ids]) for e in mesh.elements]


def root_indices(mesh):
    idx = mesh.node_id_to_index
    return [idx[n] for n in mesh.get_node_set("RootNodes").nodes]


def pin_twist_sign(axis, mesh_chord_deg, mesh_s_frac):
    """Twist sign (+1/-1) whose predicted chord angle (mod 180 deg, shell xy) matches the mesh chord axes best.

    BD y_d rotated by +twist maps to the shell direction (-cos t, -sin t), i.e. the chord angle is +t mod 180.
    """
    tw = np.interp(mesh_s_frac * axis["length"], axis["s"], axis["twist_deg"])

    def err(sign):
        d = (np.asarray(mesh_chord_deg) - sign * tw + 90.0) % 180.0 - 90.0
        return float(np.median(np.abs(d)))

    e_p, e_m = err(+1.0), err(-1.0)
    return (1.0 if e_p <= e_m else -1.0), e_p, e_m


def mesh_chord_angles(coords, membership, sts):
    """Chord axis angle [deg, mod 180] of each station's nodes (xy principal axis) and its z."""
    ang, z = [], []
    for st in sts:
        p = coords[membership[st]["nodes"]][:, :2]
        d = p - p.mean(axis=0)
        w, v = np.linalg.eigh(d.T @ d)
        ang.append(np.degrees(np.arctan2(v[1, 1], v[0, 1])) % 180.0)
        z.append(float(coords[membership[st]["nodes"], 2].mean()))
    return np.asarray(ang), np.asarray(z)


# =============================================================================================
# selftest of the energy integral
# =============================================================================================
def selftest() -> int:
    length = 20.0
    kp = np.array([[0.0, 0.0, z, 30.0] for z in np.linspace(0.0, length, 5)])
    axis = reference_axis(kp)
    etas = np.array([0.0, 1.0])
    k44, k55, k66, kga, kea = 4.0e9, 1.0e9, 3.0e8, 2.0e8, 5.0e9
    k = np.diag([kga, kga, kea, k44, k55, k66])
    ks = np.stack([k, k])
    force = 1.0e5
    xyz, f = np.array([[0.0, 0.0, length]]), np.array([[force, 0.0, 0.0]])
    th = np.deg2rad(30.0)
    exp_b = 0.5 * force**2 * length**3 / 3.0 * (np.cos(th) ** 2 / k44 + np.sin(th) ** 2 / k55)
    exp_s = 0.5 * force**2 * length / kga
    got = beam_energy(axis, etas, ks, xyz, f)
    err = [got["W_decoupled"] / exp_b - 1.0, got["W_full"] / (exp_b + exp_s) - 1.0]
    # a tip torque
    torque = 2.0e6
    xyz2 = np.array([[1.0, 0.0, length], [-1.0, 0.0, length]])
    f2 = np.array([[0.0, torque / 2.0, 0.0], [0.0, -torque / 2.0, 0.0]])
    got2 = beam_energy(axis, etas, ks, xyz2, f2)
    err.append(got2["W_full"] / (0.5 * torque**2 * length / k66) - 1.0)
    # uniform distributed x load q on a twisted cantilever: W = q^2 L^5 / 40 (cos^2/K44 + sin^2/K55) + shear q^2 L^3 / 6 K_GA
    q = 1.0e4
    zl = (np.arange(4000) + 0.5) * length / 4000
    xyz3 = np.column_stack([np.zeros(4000), np.zeros(4000), zl])
    f3 = np.column_stack([np.full(4000, q * length / 4000), np.zeros(4000), np.zeros(4000)])
    got3 = beam_energy(axis, etas, ks, xyz3, f3, distributed=True)
    exp3 = q**2 * length**5 / 40.0 * (np.cos(th) ** 2 / k44 + np.sin(th) ** 2 / k55)
    err.append(got3["W_decoupled"] / exp3 - 1.0)
    err.append(got3["W_full"] / (exp3 + q**2 * length**3 / 6.0 / kga) - 1.0)
    print("selftest twisted cantilever (30 deg): W_dec/exact-1, W_full/exact-1, torque, distributed dec, distributed full:", [f"{e:.2e}" for e in err])
    ok = all(abs(e) < 1e-5 for e in err)
    print("selftest", "PASS" if ok else "FAIL")
    return 0 if ok else 1


# =============================================================================================
# blade mode
# =============================================================================================
def _summary(z, ours, ref, lo=0.2, hi=0.85):
    z, ours, ref = np.asarray(z), np.asarray(ours, float), np.asarray(ref, float)
    rel = ours / ref - 1.0
    mid = (z > lo * z.max()) & (z < hi * z.max())
    n_nan = int(np.isnan(rel[mid]).sum())
    return {"signed_pct": float(100 * np.nanmean(rel[mid])), "mean_abs_pct": float(100 * np.nanmean(np.abs(rel[mid]))),
            "n_mid": int(mid.sum()), "n_nan_mid": n_nan}


def _nodal_vector(n_dof, dpn, nodes, forces) -> np.ndarray:
    f = np.zeros(n_dof)
    for t, fv in zip(nodes, forces, strict=True):
        f[int(t) * dpn: int(t) * dpn + 3] += fv
    return f


def blade_mode(args) -> dict:
    import yaml as _yaml

    from aeroelast.core.mesh.generators import BladeMesh
    from aeroelast.models.blade.model import build_rust_properties
    from aeroelast.postprocess.sectional import SectionalExtractor
    from beam_reference import load_beamdyn_blade
    from run_s1_sectional import ea_gj_from_static_response, native_bay_centres
    import run_s7_torsion as s7

    ra = _yaml.safe_load(open(args.yaml))["components"]["blade"]["outer_shape_bem"]["reference_axis"]
    length = float(args.blade_length or ra["z"]["values"][-1])
    deck = Path(args.deck)
    main_file = Path(args.bd_main) if args.bd_main else deck.with_name(deck.name.replace("_blade", ""))
    bd = load_beamdyn_blade(deck, blade_length=length)
    etas, ks = load_beamdyn_6x6(deck)

    gen = BladeMesh(yaml_file=args.yaml, element_size=args.element_size, skin_reference="mid")
    mesh = gen.generate(renumber="rcm")
    props = build_rust_properties(gen.numad_mesh_data)
    coords = np.asarray(mesh.coords_array)
    stations = native_bay_centres(args.yaml)
    ext = SectionalExtractor(mesh, props)
    dpn = ext.dpn
    n_dof = coords.shape[0] * dpn
    out: dict = {"yaml": args.yaml, "deck": str(deck), "bd_main": str(main_file), "blade_length": length,
                 "element_size": args.element_size, "n_nodes": int(coords.shape[0])}
    print(f"{Path(args.yaml).name} es {args.element_size}: {coords.shape[0]} nodes, L {length}")

    # ---- S-1 reads, all vs box (the solves are shared) ----------------------------------------
    info_ei, info_gj = {}, {}
    ip = lambda z, y: np.interp(z, bd.z, y)  # noqa: E731
    z_ei, fl_a, ed_a = ext.ei_from_static_response(props, frame="tangent", stations=stations, info=info_ei)
    info_ei_b = {}
    _z, fl_b, ed_b = ext.ei_from_static_response(props, frame="tangent", stations=stations, fit_nodes="box",
                                                 info=info_ei_b, solutions=info_ei["solutions"])
    sol_gj = {"axial": np.zeros(n_dof)}  # the response-EA channel is not read here
    z_gj, _e, gj_a, _ea = ea_gj_from_static_response(mesh, props, numad_data=gen.numad_mesh_data,
                                                     stations=stations, info=info_gj, solutions=sol_gj)
    info_gj_b = {}
    _z, _e, gj_b, _ea = ea_gj_from_static_response(mesh, props, numad_data=gen.numad_mesh_data, stations=stations,
                                                   fit_nodes="box", info=info_gj_b, solutions=info_gj["solutions"])
    reads = {"EI_edge": (z_ei, ed_a, ed_b, ip(z_ei, bd.ei_edge_ec)), "EI_flap": (z_ei, fl_a, fl_b, ip(z_ei, bd.ei_flap_ec)),
             "GJ": (z_gj, gj_a, gj_b, ip(z_gj, bd.gj))}
    out["s1"] = {}
    print("S-1 mid-span 0.2-0.85 vs deck, signed % (mean |rel| %), all nodes | box nodes")
    for name, (z, a, b, ref) in reads.items():
        sa, sb = _summary(z, a, ref), _summary(z, b, ref)
        out["s1"][name] = {"all": sa, "box": sb, "z": np.asarray(z).tolist(), "ref": np.asarray(ref).tolist(),
                           "all_values": np.asarray(a).tolist(), "box_values": np.asarray(b).tolist()}
        print(f"  {name:8s} {sa['signed_pct']:+7.2f} ({sa['mean_abs_pct']:5.2f}) | {sb['signed_pct']:+7.2f} "
              f"({sb['mean_abs_pct']:5.2f})  NaN stations box: {sb['n_nan_mid']} of {sb['n_mid']} mid")
    out["n_box_nan_all_stations"] = {"ei": info_ei_b.get("n_box_nan"), "gj": info_gj_b.get("n_box_nan")}
    bands = {}
    eta_gj = np.asarray(z_gj) / length
    for lo, hi in ((0.2, 0.3), (0.3, 0.6), (0.6, 0.85)):
        m = (eta_gj > lo) & (eta_gj < hi)
        bands[f"{lo}-{hi}"] = {"all": float(100 * np.nanmean(gj_a[m] / ip(z_gj, bd.gj)[m] - 1)),
                               "box": float(100 * np.nanmean(gj_b[m] / ip(z_gj, bd.gj)[m] - 1))}
    out["gj_band_signed_pct"] = bands
    print("  GJ band signed %:", {k: (round(v["all"], 1), round(v["box"], 1)) for k, v in bands.items()})

    # ---- S-7 band ratio ------------------------------------------------------------------------
    u7 = np.asarray(s7.run_torsion_case(mesh, props, S7_MOMENT), float)
    zs7, th_a = s7.section_twists_deg(mesh, u7)
    _z7, th_b = s7.section_twists_deg(mesh, u7, fit_nodes="box")
    th_beam = s7.beam_twist_profile_deg(S7_MOMENT, zs7, length, beamdyn_file=deck, z_cut=S7_ETA_CUT * length)
    band = (zs7 >= 0.3 * length) & (zs7 <= 0.9 * length)
    ratio = {k: float(np.nanmean(th[band] / th_beam[band])) for k, th in (("all", th_a), ("box", th_b))}
    out["s7"] = {"band": [0.3, 0.9], "ratio_all": ratio["all"], "ratio_box": ratio["box"],
                 "n_nan_box_band": int(np.isnan(th_b[band]).sum()), "n_band": int(band.sum()),
                 "mesh": "the S-1 mesh (not the S-7 test's sampling)", "z": zs7.tolist(),
                 "theta_all_deg": th_a.tolist(), "theta_box_deg": th_b.tolist(), "theta_beam_deg": th_beam.tolist()}
    print(f"S-7 band 0.3-0.9 ratio shell/beam: all {ratio['all']:.4f}  box {ratio['box']:.4f}  "
          f"(NaN box slices in band: {out['s7']['n_nan_box_band']} of {int(band.sum())})")

    # ---- energy --------------------------------------------------------------------------------
    membership = ext.station_membership()
    sts = sorted(membership)
    tip_nodes = np.asarray(info_ei["tip_nodes"])
    load = float(info_ei["load"])
    tip_note = "contaminated by tip load introduction (rectangle closed form 2.45 / 1.18)"
    cases = {}  # name -> (nodes, forces, u, distributed, note)
    for name, comp in (("tip_flap_Y", 1), ("tip_edge_X", 0)):
        forces = np.zeros((len(tip_nodes), 3))
        forces[:, comp] = load / len(tip_nodes)
        cases[name] = (tip_nodes, forces, info_ei["solutions"]["flap" if comp == 1 else "edge"], False, tip_note)
    cases["tip_torque"] = (info_gj["torque_nodes"], info_gj["torque_forces"], info_gj["solutions"]["torque"], False,
                           "clean (tube torque closed form 0.9995)")
    elem_nodes = mesh_element_nodes(mesh)
    for name, comp in (("dist_flap_Y", 1), ("dist_edge_X", 0)):
        nd, fd = distributed_nodal_loads(coords, elem_nodes, root_indices(mesh), comp, load)
        cases[name] = (nd, fd, solve_nodal(mesh, props, nd, fd), True,
                       "span-distributed, weight = tributary span length, all nodes but the root ring; beam side = same nodal forces")

    kp = load_keypoints(main_file)
    axis = reference_axis(kp)
    ang, zst = mesh_chord_angles(coords, membership, sts)
    sel = (zst > 0.1 * length) & (zst < 0.7 * length)
    tw_here = np.interp(zst[sel] / length * axis["length"], axis["s"], axis["twist_deg"])
    use = tw_here > 3.0
    sign, e_p, e_m = pin_twist_sign(axis, ang[sel][use], zst[sel][use] / length)
    pos_shell = axis["pos"] @ BD_TO_SHELL.T
    cen = np.array([coords[membership[st]["nodes"]].mean(axis=0) for st in sts])
    d_axis = np.abs(np.column_stack([np.interp(zst, pos_shell[:, 2], pos_shell[:, d]) for d in (0, 1)]) - cen[:, :2])
    out["frame"] = {"twist_sign": sign, "chord_angle_median_abs_err_deg": {"+1": e_p, "-1": e_m},
                    "axis_minus_station_centroid_xy_max_m": d_axis.max(axis=0).tolist(),
                    "keypoint_arc_length": axis["length"], "rule": "BD (x,y,z) -> shell (-y,x,z); Rz(sign*twist), tilt to tangent"}
    print(f"deck frame: twist sign {sign:+.0f} (chord-angle error median {e_p:.2f} deg for +1, {e_m:.2f} deg for -1, "
          f"{int(use.sum())} stations with twist > 3 deg); axis - station centroid xy max {d_axis.max(axis=0).round(3).tolist()} m")

    out["energy"] = {}
    print(f"{'case':>12} {'W_FE [J]':>11} " + " ".join(f"{'full@' + str(c):>12}{'dec@' + str(c):>12}" for c in ETA_CUTS))
    for name, (nodes, forces, u, dist, note) in cases.items():
        f_vec = _nodal_vector(n_dof, dpn, nodes, forces)
        w_fe = 0.5 * float(f_vec @ u)
        rec = {"W_FE": w_fe, "label": note, "distributed": dist, "F_sum": np.asarray(forces).sum(axis=0).tolist(),
               "M_about_origin": np.cross(coords[np.asarray(nodes)], np.asarray(forces)).sum(axis=0).tolist(), "cuts": []}
        for c in ETA_CUTS:
            wb = beam_energy(axis, etas, ks, coords[np.asarray(nodes)], np.asarray(forces), sign, c, distributed=dist)
            wb.update(ratio_full=w_fe / wb["W_full"], ratio_decoupled=w_fe / wb["W_decoupled"])
            rec["cuts"].append(wb)
        out["energy"][name] = rec
        print(f"{name:>12} {w_fe:11.4e} " + " ".join(f"{c['ratio_full']:12.4f}{c['ratio_decoupled']:12.4f}" for c in rec["cuts"])
              + ("   [tip load: " + note + "]" if name.startswith("tip_") and name != "tip_torque" else ""))
    out["energy_note"] = ("columns are W_FE / W_beam for eta_max in " + str(ETA_CUTS) + "; full = deck 6x6 inverse, dec = "
                          "M_x^2/K44 + M_y^2/K55 + M_z^2/K66; the torque case is tip-dominated, read eta_max 0.95 "
                          "(the S-7 cut 0.9504)")
    out["deck_frame_sign_alt"] = {}
    for name, (nodes, forces, _u, dist, _note) in cases.items():
        w_fe = out["energy"][name]["W_FE"]
        wb = beam_energy(axis, etas, ks, coords[np.asarray(nodes)], np.asarray(forces), -sign, 0.95, distributed=dist)
        out["deck_frame_sign_alt"][name] = w_fe / wb["W_full"]
    print("opposite twist sign, eta_max 0.95, W_FE/W_full:", {k: round(v, 4) for k, v in out["deck_frame_sign_alt"].items()})
    return out


# =============================================================================================
# tube closed forms
# =============================================================================================
def ring_section(tube, shape, g):
    """Bending I (about y, force along x), shear energy per V^2 per unit length, E and the Bredt GJ of the polygon ring."""
    pts = np.asarray(tube.ring_points(shape))
    q = np.roll(pts, -1, axis=0)
    seg = np.linalg.norm(q - pts, axis=1)
    ratio_t = np.array([tube.wall_thickness_ratio(shape, a, b) for a, b in zip(pts, q, strict=True)])
    th = tube.T_WALL * ratio_t
    e_mod = 2.0 * (1.0 + tube.NU) * g
    xm, dx = 0.5 * (pts[:, 0] + q[:, 0]), q[:, 0] - pts[:, 0]
    cx = np.sum(th * seg * xm) / np.sum(th * seg)
    i_bend = float(np.sum(th * seg * ((xm - cx) ** 2 + dx**2 / 12.0)))
    # shear flow of a unit force along x: q(s) = (1/I) int t (x - cx) ds from the ring point on the symmetry plane
    # (y = 0, x > cx, where q = 0); energy int q^2/(G t) ds, each polygon side subdivided
    sub, energy, q_run = 200, 0.0, 0.0
    i0 = [i for i in range(len(pts)) if abs(pts[i, 1]) < 1e-9 and pts[i, 0] > cx][0]
    for i in [(i0 + k) % len(pts) for k in range(len(pts))]:
        for k in range(sub):
            xa = pts[i, 0] + dx[i] * (k + 0.5) / sub
            ds = seg[i] / sub
            q_run += th[i] * (xa - cx) * ds / i_bend
            energy += q_run**2 / (g * th[i]) * ds
    return i_bend, 0.5 * energy, e_mod, g * tube.T_WALL * tube.bredt_factor(shape)


def discrete_beam_energy(z_load, f_load, length, ei, shear_per_v2, n=200001) -> float:
    """Euler-Bernoulli + shear energy of a cantilever under point loads ``f_load`` at ``z_load`` (the applied nodal forces)."""
    s = np.linspace(0.0, length, n)
    order = np.argsort(z_load)
    z, f = np.asarray(z_load)[order], np.asarray(f_load)[order]
    k0 = np.searchsorted(z, s - 1e-12, side="left")
    f_suf = np.concatenate([np.cumsum(f[::-1])[::-1], [0.0]])
    zf_suf = np.concatenate([np.cumsum((z * f)[::-1])[::-1], [0.0]])
    v = f_suf[k0]
    m = zf_suf[k0] - s * v
    return float(np.trapezoid(0.5 * m**2 / ei + v**2 * shear_per_v2, s))


def tube_mode(args) -> dict:
    import diagnose_s1_gj_tube as tube
    from aeroelast.solvers.bem.force_projection import section_torque_wall_forces

    out = {"cases": {}}
    for length in args.tube_length:
        nz = int(round(length / tube.ES))
        for shape in ("circle", "rectw"):
            g = tube.G_UNIFORM
            mesh, props, zl = tube.build_tube(np.full(nz, g), shape, tube.T_WALL, length=length)
            mesh.coords_array = np.asarray([n.coords for n in mesh.nodes], dtype=float)
            coords = mesh.coords_array
            n_dof = coords.shape[0] * 6
            tip = np.nonzero(coords[:, 2] >= coords[:, 2].max() - 1e-6)[0]
            i_bend, shear_v2, e_mod, gj = ring_section(tube, shape, g)

            def w_fe(nodes, forces):
                u = solve_nodal(mesh, props, nodes, forces)
                return 0.5 * float(_nodal_vector(n_dof, 6, nodes, forces) @ u)

            rec = {}
            torque = 1.0e6
            nodes, f_t = section_torque_wall_forces(mesh, tip, (0.0, 0.0, 1.0), props, torque, coords=coords)
            w_t = 0.5 * torque**2 * length / gj
            rec["torque"] = {"W_FE": w_fe(nodes, f_t), "W_exact": w_t}
            rec["torque"]["ratio"] = rec["torque"]["W_FE"] / w_t
            force = 1.0e6
            f_tip = np.zeros((len(tip), 3))
            f_tip[:, 0] = force / len(tip)
            w_ex = 0.5 * force**2 * length**3 / (3.0 * e_mod * i_bend) + force**2 * length * shear_v2
            rec["tip_force"] = {"W_FE": w_fe(tip, f_tip), "W_exact": w_ex, "label": "contaminated by tip load introduction"}
            rec["tip_force"]["ratio"] = rec["tip_force"]["W_FE"] / w_ex
            # span-distributed force: tributary-span weights, root ring excluded
            nd, fd = distributed_nodal_loads(coords, mesh_element_nodes(mesh), root_indices(mesh), 0, force)
            w_d = w_fe(nd, fd)
            w_disc = discrete_beam_energy(coords[nd, 2], fd[:, 0], length, e_mod * i_bend, shear_v2)
            q = force / length
            w_text = q**2 * length**5 / (40.0 * e_mod * i_bend) + q**2 * length**3 / 6.0 * shear_v2
            rec["dist_force"] = {"W_FE": w_d, "W_exact_discrete_loads": w_disc, "ratio": w_d / w_disc,
                                 "W_textbook_qL5_40EI_plus_shear": w_text, "ratio_textbook": w_d / w_text,
                                 "W_bending_share": q**2 * length**5 / (40.0 * e_mod * i_bend) / w_text,
                                 "load_z_range": [float(coords[nd, 2].min()), float(coords[nd, 2].max())],
                                 "n_loaded": int(len(nd))}
            r = rec
            print(f"L {length:4.0f} {shape:6s} torque {r['torque']['ratio']:.5f} | tip force {r['tip_force']['ratio']:.5f} "
                  f"(contaminated) | distributed {r['dist_force']['ratio']:.5f} vs the applied nodal loads, "
                  f"{r['dist_force']['ratio_textbook']:.5f} vs q^2L^5/40EI+shear (q = F/L)")
            out["cases"][f"{shape}-L{length:g}"] = rec
    return out


# =============================================================================================
# prismatic mode
# =============================================================================================
def prismatic_mode(args) -> dict:
    import l3_bredt_mesh as bm
    from run_s1_sectional import ea_gj_from_static_response, native_bay_centres

    from aeroelast.core.mesh.generators import BladeMesh
    from aeroelast.models.blade.model import build_rust_properties

    gen, sections, mats, elements = bm.build_mesh(args.yaml, args.element_size)
    gen2 = BladeMesh(yaml_file=args.yaml, element_size=args.element_size, skin_reference="mid")
    mesh = gen2.generate(renumber="rcm")
    props = build_rust_properties(gen2.numad_mesh_data)
    rows = np.unique(np.round([zm for _, _, zm in elements], 6))
    length = float(np.asarray(mesh.coords_array)[:, 2].max())
    stations = native_bay_centres(args.yaml)
    info = {}
    z, _e, gj_a, _ = ea_gj_from_static_response(mesh, props, numad_data=gen2.numad_mesh_data, stations=stations, info=info)
    _z, _e, gj_b, _ = ea_gj_from_static_response(mesh, props, numad_data=gen2.numad_mesh_data, stations=stations,
                                                 fit_nodes="box", solutions=info["solutions"])
    z_b = 0.5 * length
    sec = bm.mesh_section(gen, sections, mats, elements, rows, z_b)
    gj_bredt = float(sec["gj"])
    mid = (z > 0.3 * length) & (z < 0.85 * length)
    out = {"yaml": args.yaml, "element_size": args.element_size, "length": length, "mesh_bredt_gj": gj_bredt,
           "z": np.asarray(z).tolist(), "gj_all": np.asarray(gj_a).tolist(), "gj_box": np.asarray(gj_b).tolist()}
    for k, g in (("all", gj_a), ("box", gj_b)):
        r = g[mid] / gj_bredt
        out[f"ratio_{k}"] = {"mean": float(np.nanmean(r)), "min": float(np.nanmin(r)), "max": float(np.nanmax(r)),
                             "n_nan": int(np.isnan(r).sum()), "n": int(mid.sum())}
        print(f"prismatic GJ/mesh Bredt ({gj_bredt:.4e}) over eta 0.3-0.85, {k:3s}: mean {out[f'ratio_{k}']['mean']:.4f} "
              f"min {out[f'ratio_{k}']['min']:.4f} max {out[f'ratio_{k}']['max']:.4f} NaN {out[f'ratio_{k}']['n_nan']}")
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--mode", choices=("blade", "tube", "prismatic", "selftest"), default="blade")
    ap.add_argument("--yaml")
    ap.add_argument("--deck", help="BeamDyn blade file (6x6 per station)")
    ap.add_argument("--bd-main", help="BeamDyn main file with the key points (default: the deck name without _blade)")
    ap.add_argument("--blade-length", type=float, help="deck eta scale (default: the yaml reference axis z at the tip)")
    ap.add_argument("--element-size", type=float, default=0.25)
    ap.add_argument("--tube-length", type=float, nargs="+", default=[24.0, 48.0])
    ap.add_argument("--out")
    args = ap.parse_args(argv)
    if args.mode == "selftest":
        return selftest()
    if not args.out:
        ap.error("--out is required")
    if args.mode in ("blade", "prismatic") and args.element_size >= 0.3:
        raise SystemExit("element_size >= 0.3 is not converged (docs/blade_mesh_size_rule.md)")
    if args.mode in ("blade", "prismatic") and not args.yaml:
        ap.error("--yaml is required")
    if args.mode == "blade" and not args.deck:
        ap.error("--deck is required")
    out = {"blade": blade_mode, "tube": tube_mode, "prismatic": prismatic_mode}[args.mode](args)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(json.dumps(out, indent=1, default=float))
    return 0


if __name__ == "__main__":
    sys.exit(main())
