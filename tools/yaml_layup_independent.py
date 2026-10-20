"""Mesh-free WindIO section integrals (S4 of blade-midsurface-reference).

Run: python -B tools/yaml_layup_independent.py [--self-check] [--output FILE]
Only NumPy, PyYAML and SciPy (plus the standard library) are imported. No mesh,
NuMAD or aeroelast code is used. Default input is the official IEA15MW YAML.

The YAML coordinates run TE -> suction -> LE -> pressure -> TE: nd_arc is
cumulative open-contour arc / total arc, NOT chord fraction or point index.
The finite TE gap is not part of that arc. Neighbouring profiles are blended
on a common arc grid, including all original vertices. Chord scales lengths;
pitch-axis translation and twist (radians) are rigid transforms and leave
these scalar integrals invariant. Web heights are straight endpoint distances
from their explicit nd_arc definitions; redundant offset/rotation data are
reported as a geometry cross-check, not added as extra webs.

EA_Q is the restrained membrane integral sum(Qbar11*t*ds), with reciprocal
nu21=nu12*E2/E1. EA_Ex uses free-transverse-stress ply Ex=1/Sbar11 instead.
Neither is a full relaxed/warping beam-section homogenization. EA_S1 uses
the S-1 diagnostic's nonreciprocal 1-nu12**2 denominator, solely to isolate
that convention from mapping differences; it is NOT physical orthotropic Q.
Thickness is integrated continuously (no integer-ply rounding). Spanwise
interpolation is PCHIP by default; --linear measures interpolation sensitivity.
Output retains every layer contribution, web height and beam reference.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import yaml
from scipy.interpolate import PchipInterpolator

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_YAML = ROOT / "tests/IEA-15-240-RWT.yaml"


def interp(field, eta, linear=False):
    """Interpolate only within the declared support; absent layers are zero."""
    grid = np.asarray(field["grid"], dtype=float)
    if eta < grid[0] or eta > grid[-1]:
        return 0.0
    values = np.asarray(field["values"], dtype=float)
    return float(np.interp(eta, grid, values) if linear else
                 PchipInterpolator(grid, values)(eta))


def moduli(material, theta):
    """Return physical Qbar11, ply Ex(theta), and legacy S-1 Qbar11, in Pa."""
    if material["orth"]:
        e1, e2 = map(float, material["E"][:2])
        nu = float(material["nu"][0])
        g = float(material["G"][0])
    else:
        e1 = e2 = float(material["E"])
        nu = float(material["nu"])
        g = e1 / (2 * (1 + nu))
    c2, s2 = np.cos(theta)**2, np.sin(theta)**2
    numerator = e1*c2*c2 + 2*(nu*e2 + 2*g*(1-nu*nu*e2/e1))*s2*c2 + e2*s2*s2
    q = numerator / (1-nu*nu*e2/e1)
    ex = 1 / (c2*c2/e1 + s2*s2/e2 + (1/g-2*nu/e1)*s2*c2)
    legacy = ((e1*c2*c2 + 2*nu*e2*s2*c2 + e2*s2*s2) /
              (1-nu*nu) + 4*g*s2*c2)
    return np.array([q, ex, legacy])


def contour(data, eta, linear=False, resolution=4001):
    """Blend profiles independently on normalized perimeter arc."""
    outer = data["components"]["blade"]["outer_shape_bem"]
    position = outer["airfoil_position"]
    grid = np.asarray(position["grid"])
    j = int(np.clip(np.searchsorted(grid, eta, side="right")-1, 0, len(grid)-2))
    weight = (eta-grid[j])/(grid[j+1]-grid[j])
    profiles = {a["name"]: a["coordinates"] for a in data["airfoils"]}
    arcs, points = [], []
    for label in position["labels"][j:j+2]:
        p = profiles[label]
        xy = np.column_stack([p["x"], p["y"]]).astype(float)
        # WindIO ordering is visible directly in the YAML: positive y first.
        if np.mean(xy[1:np.argmin(xy[:, 0]), 1]) < 0:
            xy = xy[::-1]
        arc = np.r_[0, np.cumsum(np.linalg.norm(np.diff(xy, axis=0), axis=1))]
        keep = np.r_[True, np.diff(arc) > 0]
        arcs.append(arc[keep]/arc[-1])
        points.append(xy[keep])
    u = np.unique(np.r_[np.linspace(0, 1, resolution), *arcs])
    pp = [np.column_stack([np.interp(u, a, p[:, k]) for k in (0, 1)])
          for a, p in zip(arcs, points)]
    xy = ((1-weight)*pp[0] + weight*pp[1])*interp(outer["chord"], eta, linear)
    lengths = np.r_[0, np.cumsum(np.linalg.norm(np.diff(xy, axis=0), axis=1))]
    # Reparametrize the blended shape by its own arc, as nd_arc requires.
    return lengths/lengths[-1], xy, float(lengths[-1])


def section(data, eta, linear=False, resolution=4001):
    blade = data["components"]["blade"]
    structure = blade["internal_structure_2d_fem"]
    u, xy, perimeter = contour(data, eta, linear, resolution)
    point = lambda a: np.array([np.interp(a, u, xy[:, k]) for k in (0, 1)])
    webs, web_checks = {}, {}
    outer = blade["outer_shape_bem"]
    for web in structure["webs"]:
        a, b = [point(interp(web[k], eta, linear))
                for k in ("start_nd_arc", "end_nd_arc")]
        webs[web["name"]] = float(np.linalg.norm(b-a))
        chord = interp(outer["chord"], eta, linear)
        # Check the redundant line definition in the profile frame. The
        # published rotation=-twist gives normal (cos(twist), -sin(twist)).
        rotation = interp(web["rotation"], eta, linear)
        normal = np.array([np.cos(rotation), np.sin(rotation)])
        origin = np.array([interp(outer["pitch_axis"], eta, linear)*chord, 0.])
        projected = np.array([(p-origin) @ normal for p in (a, b)])
        offset = interp(web["offset_y_pa"], eta, linear)
        web_checks[web["name"]] = dict(endpoint_offset_m=projected.tolist(),
                                      specified_offset_m=offset,
                                      max_offset_residual_m=float(np.max(abs(projected-offset))),
                                      rotation_plus_twist_rad=rotation+interp(outer["twist"], eta, linear))
    mats = {m["name"]: m for m in data["materials"]}
    layers = []
    for layer in structure["layers"]:
        h = interp(layer["thickness"], eta, linear)
        if "web" in layer:
            width = webs[layer["web"]]
        elif "start_nd_arc" not in layer or "end_nd_arc" not in layer:
            # TE bands have one fixed TE endpoint and a physical width on
            # the declared side, rather than a second normalized endpoint.
            anchor = layer.get("start_nd_arc", layer.get("end_nd_arc"))
            if anchor.get("fixed") != "TE":
                raise ValueError(f"Unsupported single-ended layer: {layer['name']}")
            width = interp(layer["width"], eta, linear)
        else:
            start, end = [interp(layer[k], eta, linear) for k in ("start_nd_arc", "end_nd_arc")]
            if not 0 <= start <= end <= 1:
                raise ValueError(f"Invalid nd_arc interval in {layer['name']}: {start}, {end}")
            width = (end-start)*perimeter
        mat = mats[layer["material"]]
        theta = interp(layer["fiber_orientation"], eta, linear)
        ea = moduli(mat, theta)*h*width
        layers.append(dict(name=layer["name"], material=layer["material"],
                           width_m=width, thickness_m=h, mass=float(mat["rho"])*h*width,
                           ea_q=ea[0], ea_ex=ea[1], ea_s1=ea[2]))
    target = blade["elastic_properties_mb"]["six_x_six"]
    # WindIO stores the upper triangle in row-major order: K33 index 11.
    beam_ea = dict(grid=target["stiff_matrix"]["grid"],
                   values=[r[11] for r in target["stiff_matrix"]["values"]])
    beam_mass = dict(grid=target["inertia_matrix"]["grid"],
                     values=[r[0] for r in target["inertia_matrix"]["values"]])
    return dict(eta=float(eta), perimeter_m=perimeter, web_heights_m=webs,
                web_checks=web_checks, layers=layers,
                **{k:sum(l[k] for l in layers) for k in ("mass", "ea_q", "ea_ex", "ea_s1")},
                beam_ea=interp(beam_ea, eta, True), beam_mass=interp(beam_mass, eta, True))


def load(path=DEFAULT_YAML):
    with open(path) as stream:
        return yaml.safe_load(stream)


def self_check():
    """Analytical checks, no fitted blade reference or mutable test stores."""
    iso = dict(orth=0, E=10e9, nu=0.3)
    for theta in (0, 0.3, np.pi/2):
        np.testing.assert_allclose(moduli(iso, theta), [10e9/0.91, 10e9, 10e9/0.91])
    orth = dict(orth=1, E=[100e9, 10e9], G=[5e9], nu=[0.3])
    np.testing.assert_allclose(moduli(orth, 0), [100e9/0.991, 100e9, 100e9/0.91])
    np.testing.assert_allclose(moduli(orth, np.pi/2), [10e9/0.991, 10e9, 10e9/0.91])
    constant = lambda v: dict(grid=[0., 1.], values=[v, v])
    outer = dict(airfoil_position=dict(grid=[0., 1.], labels=['square', 'square']),
                 chord=constant(2.), pitch_axis=constant(.5), twist=constant(0.))
    layer = dict(name='skin', material='iso', thickness=constant(.01),
                 start_nd_arc=constant(0.), end_nd_arc=constant(1.),
                 fiber_orientation=constant(.3))
    web = dict(name='web', start_nd_arc=constant(.125), end_nd_arc=constant(.625),
               rotation=constant(0.), offset_y_pa=constant(0.))
    web_layer = dict(layer, name='web_skin', web='web')
    synthetic = dict(components=dict(blade=dict(outer_shape_bem=outer,
        internal_structure_2d_fem=dict(layers=[layer, web_layer], webs=[web]),
        elastic_properties_mb=dict(six_x_six=dict(stiff_matrix=dict(grid=[0., 1.], values=[[1.]*21]*2),
                                                inertia_matrix=dict(grid=[0., 1.], values=[[1.]*21]*2))))),
        materials=[dict(iso, name='iso', rho=1000.)],
        airfoils=[dict(name='square', coordinates=dict(x=[1., 1., 0., 0., 1.], y=[0., 1., 1., 0., 0.]))])
    result = section(synthetic, .5)
    np.testing.assert_allclose([result['perimeter_m'], result['web_heights_m']['web'], result['mass']], [8., 2., 100.])
    np.testing.assert_allclose([result['ea_q'], result['ea_ex']], [1e9/.91, 1e9])
    print("SELF-CHECK PASS: isotropic rotation invariance; orthotropic 0/90 Q/Ex; square skin+web geometry/mass/EA")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--blade-yaml", type=Path, default=DEFAULT_YAML)
    parser.add_argument("--output", type=Path, help="Optional JSON output (choose an authorized log path)")
    parser.add_argument("--linear", action="store_true", help="Use linear rather than PCHIP spanwise interpolation")
    parser.add_argument("--self-check", action="store_true", help="Run analytical constitutive checks before the blade integral")
    args = parser.parse_args()
    if args.self_check:
        self_check()
    data = load(args.blade_yaml)
    rows = [section(data, eta, args.linear) for eta in np.linspace(0.2, 0.85, 12)]
    print("eta perimeter[m] web0[m] web1[m] mass[kg/m] EA_Q[GN] EA_Ex[GN] EA_S1[GN] mass/beam-1[%] Q/beam-1[%]")
    for r in rows:
        print(f"{r['eta']:.6f} {r['perimeter_m']:.5f} "
              f"{r['web_heights_m']['web0']:.5f} {r['web_heights_m']['web1']:.5f} "
              f"{r['mass']:.5f} {r['ea_q']/1e9:.6f} {r['ea_ex']/1e9:.6f} {r['ea_s1']/1e9:.6f} "
              f"{100*(r['mass']/r['beam_mass']-1):+.4f} {100*(r['ea_q']/r['beam_ea']-1):+.4f}")
    for key in ("mass", "ea_q", "ea_ex", "ea_s1"):
        ratios = np.array([r[key]/r['beam_mass' if key == 'mass' else 'beam_ea']-1 for r in rows])
        print(f"MEAN {key}: signed={100*ratios.mean():+.6f}% absolute={100*np.abs(ratios).mean():.6f}%")
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(rows, indent=2)+"\n")


if __name__ == "__main__":
    main()
