#!/usr/bin/env python
"""S-1 GJ extractor on a prismatic shell tube of known GJ (blade-input-audit T39).

Runs ``run_s1_sectional.ea_gj_from_static_response`` (T40: tip torque as the end ring's Saint-Venant wall
shear flow, twist rate from the local quadratic fit; T39 measured the previous couple + spline) on an isotropic MITC4 rectangular
tube whose torsional stiffness is known per station by thin-wall Bredt, ``GJ = G t 4 A^2 / P``
(mid-line area and perimeter). The tube copies what the S-1 extraction sees on the blade: 117 m long,
root clamped, one 0.25 m station element set per layer (``BladeMesh.max_bay_length``), the
``native_bay_centres`` stations of the official yaml and the S-1 torque of 1e6 N m.

Cases (rectangle 4.0 x 2.0 m or circle of the same 12 m perimeter, mid-line):
  A  rectangle, uniform wall, constant GJ.
  B  rectangle, wall shear modulus per station chosen so that the Bredt GJ(z) equals the BeamDyn
     K66(z) the S-1 assertion compares against: the theta(z) signal of the ideal blade. The thickness
     is kept constant so the membrane/bending ratio of the wall does not change along the span.
  C  circle, uniform (the extractor's tip couple t x r_perp is then the uniform Bredt shear flow).
  D  circle, BeamDyn GJ(z) as in B.
  E  rectangle A with a 5x thicker wall and G/5 (same GJ, wall bending stiffness x25).
  F  warping-free rectangle (T40): 4 m walls t, 2 m walls t/2 so b/t_b = h/t_h and Saint-Venant torsion
     does not warp the section; A's residual under the fixed extractor is the clamped root's
     warping-restraint layer, which F isolates away.

Printed per selected station: the extractor's GJ against Bredt, the unsmoothed central difference of
the extractor's own theta (the rigid-fit part without the spline), and theta against the exact
integral of T / GJ. The summary is the S-1 metric: mean |rel err| over 0.2-0.85 of the station range.

Usage:  scripts/aeroenv.sh python tools/diagnose_s1_gj_tube.py [--case A B C D E] [--metric s1|s7]
        [--all-stations]

``--metric s7`` runs the S-7 construction instead (``run_s7_torsion``: the tip torque as the wall shear
flow and the 40-slice rigid-fit twist) and reports the shell/exact twist ratio over its 0.3-0.9 band.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))

import run_s1_sectional as s1  # noqa: E402
from aeroelast.postprocess import sectional  # noqa: E402
from beam_reference import load_beamdyn_blade  # noqa: E402

B_MID = 4.0  # mid-line width [m]
H_MID = 2.0  # mid-line height [m]
LENGTH = 117.0  # tube length = blade length [m]
ES = 0.25  # element size around and along [m]
T_WALL = 0.02  # wall thickness [m]
NU = 0.3
TORQUE = 1.0e6  # the S-1 default [N m]
G_UNIFORM = 2.5e10  # case A wall shear modulus [Pa]
CASES = {  # name: (shape, wall thickness, GJ profile)
    "A": ("rect", T_WALL, "uniform"),
    "B": ("rect", T_WALL, "beamdyn"),
    "C": ("circle", T_WALL, "uniform"),
    "D": ("circle", T_WALL, "beamdyn"),
    "E": ("rect", 5.0 * T_WALL, "uniform"),
    "F": ("rectw", T_WALL, "uniform"),  # warping-free rectangle: 4 m walls t, 2 m walls t / 2
}


def ring_points(shape="rect"):
    if shape == "circle":
        n = int(round(2.0 * (B_MID + H_MID) / ES))
        rad = 2.0 * (B_MID + H_MID) / (2.0 * np.pi)
        ang = 2.0 * np.pi * np.arange(n) / n
        return [(rad * np.cos(a), rad * np.sin(a)) for a in ang]
    nb, nh = int(round(B_MID / ES)), int(round(H_MID / ES))
    pts = []
    for j in range(nb):
        pts.append((-B_MID / 2 + j * B_MID / nb, -H_MID / 2))
    for j in range(nh):
        pts.append((B_MID / 2, -H_MID / 2 + j * H_MID / nh))
    for j in range(nb):
        pts.append((B_MID / 2 - j * B_MID / nb, H_MID / 2))
    for j in range(nh):
        pts.append((-B_MID / 2, H_MID / 2 - j * H_MID / nh))
    return pts


def wall_thickness_ratio(shape, p, q):
    """Thickness of the wall p -> q in units of the case thickness: ``rectw`` halves the 2 m walls."""
    if shape == "rectw" and abs(q[1] - p[1]) > abs(q[0] - p[0]):
        return 0.5
    return 1.0


def bredt_factor(shape):
    """4 A^2 / sum(l_i t_wall / t_i) of the mid-line polygon (GJ = G t_wall * factor)."""
    p = np.asarray(ring_points(shape))
    q = np.roll(p, -1, axis=0)
    area = 0.5 * abs(np.sum(p[:, 0] * q[:, 1] - q[:, 0] * p[:, 1]))
    length = np.linalg.norm(q - p, axis=1)
    ratio = np.array([wall_thickness_ratio(shape, a, b) for a, b in zip(p, q, strict=True)])
    return 4.0 * area**2 / np.sum(length / ratio)


def build_tube(g_of_layer, shape="rect", t_wall=T_WALL, length=LENGTH):
    """MeshModel with one station element set per layer and the per-layer ShellProperty dict."""
    from aeroelast.core.material import IsotropicMaterial
    from aeroelast.core.mesh import ElementSet, MeshElement, MeshModel, Node, NodeSet
    from aeroelast.core.mesh.entities import ElementType
    from aeroelast.core.properties import ShellProperty

    Node._id_counter = 0
    MeshElement._id_counter = 0
    mesh = MeshModel()
    ring = ring_points("rect" if shape == "rectw" else shape)
    nz = int(round(length / ES))
    zl = np.linspace(0.0, length, nz + 1)
    rings = []
    for z in zl:
        row = []
        for x, y in ring:
            n = Node([float(x), float(y), float(z)], geometric_node=False)
            mesh.add_node(n)
            row.append(n)
        rings.append(row)
    props = {}
    nr = len(ring)
    for k in range(nz):
        els = {}
        for i in range(nr):
            j = (i + 1) % nr
            e = MeshElement(nodes=[rings[k][i], rings[k][j], rings[k + 1][j], rings[k + 1][i]],
                            element_type=ElementType.quad)
            mesh.add_element(e)
            ratio = wall_thickness_ratio(shape, ring[i], ring[j])
            els.setdefault(ratio, set()).add(e)
        g = float(g_of_layer[k])
        for comp, (ratio, group) in enumerate(sorted(els.items(), reverse=True)):
            name = f"{comp:02d}_{k:03d}_tube"  # one set per wall thickness; same station field
            mesh.add_element_set(ElementSet(name, group))
            props[name] = ShellProperty(
                material=IsotropicMaterial(name=f"m{k}_{comp}", E=2.0 * (1.0 + NU) * g, nu=NU, rho=1000.0),
                thickness=t_wall * ratio)
    mesh.add_node_set(NodeSet("RootNodes", set(rings[0])))
    return mesh, props, zl


class _Recorder:
    """Wraps ``sectional._local_quadratic`` to capture the theta the extractor differentiates."""

    def __init__(self, fn):
        self.fn, self.calls = fn, []

    def __call__(self, zs, y, half, deriv):
        self.calls.append((np.array(zs), np.array(y), half))
        return self.fn(zs, y, half, deriv)


def run_case(name, gj_layer, stations):
    shape, t_wall, _ = CASES[name]
    g_layer = gj_layer / (t_wall * bredt_factor(shape))
    mesh, props, zl = build_tube(g_layer, shape, t_wall)
    rec = _Recorder(sectional._local_quadratic)
    sectional._local_quadratic = rec
    try:
        zs, _ea, gj_ext, _ = s1.ea_gj_from_static_response(mesh, props, torque=TORQUE, stations=stations)
    finally:
        sectional._local_quadratic = rec.fn
    zt, theta, _half = rec.calls[-1]
    assert np.allclose(zt, zs)
    layer = np.clip(np.searchsorted(zl, zs) - 1, 0, len(gj_layer) - 1)
    gj_ref = gj_layer[layer]
    # exact theta at z: integral of T/GJ over the layers (the tube is clamped at z = 0)
    flex = np.concatenate([[0.0], np.cumsum(TORQUE / gj_layer * np.diff(zl))])
    theta_exact = np.interp(zs, zl, flex)
    rate_cd = np.abs(np.gradient(theta, zs))
    gj_cd = TORQUE / rate_cd
    # the S-1 test selects mid on z / z.max(); the tube starts at 0 like the blade mesh
    mid_s1 = (zs > 0.2 * zs.max()) & (zs < 0.85 * zs.max())
    err = gj_ext / gj_ref - 1
    err_cd = gj_cd / gj_ref - 1
    sign = np.sign(np.median(theta[zs > 0.1 * zs.max()]))
    th_err = sign * theta / theta_exact - 1
    print(f"\n== case {name}: {len(zs)} stations, {len(zl) - 1} layers, mesh {len(mesh.nodes)} nodes")
    print(f"{'z':>7} {'GJ_bredt':>10} {'extract%':>9} {'cdiff%':>8} {'theta%':>8}")
    for i in range(len(zs)):
        print(f"{zs[i]:7.2f} {gj_ref[i]:10.3e} {100 * err[i]:+9.2f} {100 * err_cd[i]:+8.2f} {100 * th_err[i]:+8.2f}")
    out = {
        "case": name, "shape": shape, "t_wall": t_wall, "n_stations": int(len(zs)),
        "mid_mean_abs_extract": float(np.nanmean(np.abs(err[mid_s1]))),
        "mid_mean_signed_extract": float(np.nanmean(err[mid_s1])),
        "mid_mean_abs_cdiff": float(np.nanmean(np.abs(err_cd[mid_s1]))),
        "mid_mean_signed_cdiff": float(np.nanmean(err_cd[mid_s1])),
        "mid_max_abs_theta": float(np.nanmax(np.abs(th_err[mid_s1]))),
        "theta_sign": float(sign),
    }
    print(f"mid-span (0.2-0.85 z.max) mean |rel err|: extractor {out['mid_mean_abs_extract']:.2%} "
          f"(signed {out['mid_mean_signed_extract']:+.2%}), unsmoothed cdiff {out['mid_mean_abs_cdiff']:.2%} "
          f"(signed {out['mid_mean_signed_cdiff']:+.2%}), max |theta err| {out['mid_max_abs_theta']:.2%}")
    print("TUBE_JSON " + json.dumps(out))
    return out


def run_s7_case(name, gj_layer):
    """The S-7 construction (``run_s7_torsion``: tip couple F_y ~ x_off, 40 slice rigid fit) on the tube."""
    from run_s7_torsion import run_torsion_case, section_twists_deg

    shape, t_wall, _ = CASES[name]
    g_layer = gj_layer / (t_wall * bredt_factor(shape))
    mesh, props, zl = build_tube(g_layer, shape, t_wall)
    if getattr(mesh, "coords_array", None) is None:
        mesh.coords_array = np.asarray([n.coords for n in mesh.nodes], dtype=float)
    from aeroelast.solvers.bem.force_projection import section_torque_wall_forces

    z = mesh.coords_array[:, 2]
    ring, f_ring = section_torque_wall_forces(
        mesh, np.nonzero(z >= z.max() - 1e-6)[0], (0.0, 0.0, 1.0), props, TORQUE)
    r_ring = mesh.coords_array[ring, :2] - mesh.coords_array[ring, :2].mean(axis=0)
    delivered = float(np.sum(r_ring[:, 0] * f_ring[:, 1] - r_ring[:, 1] * f_ring[:, 0]))
    u = np.asarray(run_torsion_case(mesh, props, TORQUE), dtype=float)
    zs, th = section_twists_deg(mesh, u)
    flex = np.concatenate([[0.0], np.cumsum(TORQUE / gj_layer * np.diff(zl))])
    th_exact = np.rad2deg(np.interp(zs, zl, flex))
    ratio = th / th_exact
    band = (zs >= 0.3 * LENGTH) & (zs <= 0.9 * LENGTH)
    print(f"\n== S-7 construction, case {name}: delivered torque {delivered:.6e}, {len(zs)} slices")
    print(f"{'z':>7} {'theta_exact':>12} {'ratio':>8}")
    for i in range(len(zs)):
        print(f"{zs[i]:7.2f} {th_exact[i]:12.5e} {ratio[i]:8.4f}")
    out = {"case": name, "metric": "s7", "band_ratio_mean": float(np.nanmean(ratio[band])),
           "band_ratio_min": float(np.nanmin(ratio[band])), "band_ratio_max": float(np.nanmax(ratio[band]))}
    print(f"S-7 band 0.3-0.9 mean ratio {out['band_ratio_mean']:.4f} "
          f"(min {out['band_ratio_min']:.4f}, max {out['band_ratio_max']:.4f}); exact 1")
    print("TUBE_JSON " + json.dumps(out))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--case", nargs="+", choices=sorted(CASES), default=["A", "B"])
    ap.add_argument("--metric", choices=["s1", "s7"], default="s1",
                    help="s1: ea_gj_from_static_response; s7: run_s7_torsion tip couple + slice twist")
    ap.add_argument("--all-stations", action="store_true", help="stations=None (every 0.25 m layer)")
    args = ap.parse_args()
    yaml = ROOT / "tests" / "IEA-15-240-RWT.yaml"
    stations = None if args.all_stations else s1.native_bay_centres(yaml)
    nz = int(round(LENGTH / ES))
    zmid = (np.arange(nz) + 0.5) * ES
    bd = load_beamdyn_blade(ROOT / "tests" / "IEA15MW" / "reference" / "IEA-15-240-RWT_BeamDyn_blade.dat")
    gj_uniform = G_UNIFORM * T_WALL * bredt_factor("rect")
    for name in args.case:
        prof = CASES[name][2]
        gj = np.full(nz, gj_uniform) if prof == "uniform" else np.interp(zmid, bd.z, bd.gj)
        if args.metric == "s7":
            run_s7_case(name, gj)
        else:
            run_case(name, gj, stations)
    return 0


if __name__ == "__main__":
    sys.exit(main())
