"""Mesh-free WindIO 2.x section integrals (S4 of blade-midsurface-reference; v2 port T18).

Run: python -B tools/yaml_layup_independent.py [--self-check] [--output FILE]
Only NumPy, PyYAML and SciPy (plus the standard library) are imported. No mesh,
gmsh, NuMAD or aeroelast code is used; in particular the mesher's anchor/region
resolution (gmsh_skin, windio2_basis) and normalize_windio are NOT used: the
anchors below are resolved here, independently. Default input is the official
windIO 2.x IEA-15-240-RWT.yaml (tests/IEA15MW/reference/upstream). windIO 1.x
input (outer_shape_bem / internal_structure_2d_fem) is rejected.

windIO 2 schema read here: components.blade.outer_shape.{chord, twist [deg],
airfoils[{name, spanwise_position, weight}]}, components.blade.structure.{webs,
layers, anchors, elastic_properties}. A layer/web arc endpoint is either a
{grid, values} field or {anchor: {name, handle}}; the named anchor (a structure
anchor or a web-owned anchor) holds the field (or another reference) under that
handle. TE reinforcement anchors also carry a physical `width` (+ `defines`),
kept as a cross-check against the resolved nd_arc interval, not used to define it.

The YAML coordinates run TE -> suction -> LE -> pressure -> TE: nd_arc is
cumulative open-contour arc / total arc, NOT chord fraction or point index.
The finite TE gap is not part of that arc. Neighbouring profiles are blended
on a common arc grid, including all original vertices. Chord scales lengths;
pitch-axis translation and twist (radians) are rigid transforms and leave
these scalar integrals invariant. Web heights are straight endpoint distances
from their resolved nd_arc anchors (windIO 2 has no redundant offset/rotation).

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
DEFAULT_YAML = ROOT / "tests/IEA15MW/reference/upstream/IEA-15-240-RWT.yaml"


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
    outer = data["components"]["blade"]["outer_shape"]
    stations = outer["airfoils"]
    for st in stations:
        if list(st.get("weight", [1.0])) != [1.0]:
            raise ValueError(f"Blended airfoil weights are not supported: {st['name']}")
    grid = np.array([st["spanwise_position"] for st in stations], dtype=float)
    j = int(np.clip(np.searchsorted(grid, eta, side="right")-1, 0, len(grid)-2))
    weight = (eta-grid[j])/(grid[j+1]-grid[j])
    profiles = {a["name"]: a["coordinates"] for a in data["airfoils"]}
    arcs, points = [], []
    for label in (st["name"] for st in stations[j:j+2]):
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


def anchor_table(structure):
    """Name -> anchor for the structure anchors and the web-owned anchors (names must be unique)."""
    table = {}
    owned = list(structure.get("anchors", []))
    for web in structure.get("webs", []):
        owned += web.get("anchors", [])
    for anchor in owned:
        if anchor["name"] in table:
            raise ValueError(f"Duplicate anchor name: {anchor['name']}")
        table[anchor["name"]] = anchor
    return table


def arc(table, spec, eta, linear=False, _seen=()):
    """Resolve an nd_arc endpoint: a {grid, values} field or {anchor: {name, handle}} (recursive)."""
    if "anchor" not in spec:
        return interp(spec, eta, linear)
    name, handle = spec["anchor"]["name"], spec["anchor"]["handle"]
    if (name, handle) in _seen:
        raise ValueError(f"Cyclic anchor reference at {name}.{handle}")
    if name not in table:
        raise ValueError(f"Unknown anchor {name!r}")
    if handle not in table[name]:
        raise ValueError(f"Anchor {name!r} has no handle {handle!r}")
    return arc(table, table[name][handle], eta, linear, _seen + ((name, handle),))


def layer_interval(table, layer, eta, linear=False):
    """(start, end) nd_arc of a layer at eta; inverted only where the layer is absent (h = 0) is tolerated."""
    start, end = (arc(table, layer[k], eta, linear) for k in ("start_nd_arc", "end_nd_arc"))
    if end < start:
        if interp(layer["thickness"], eta, linear) > 0:
            raise ValueError(f"Invalid nd_arc interval in {layer['name']}: {start}, {end}")
        end = start
    if not -1e-12 <= start <= end <= 1 + 1e-12:
        raise ValueError(f"Invalid nd_arc interval in {layer['name']}: {start}, {end}")
    return start, end


def web_endpoints(table, web, eta, u, xy, linear=False):
    point = lambda a: np.array([np.interp(a, u, xy[:, k]) for k in (0, 1)])
    return [point(arc(table, web[k], eta, linear)) for k in ("start_nd_arc", "end_nd_arc")]


def section(data, eta, linear=False, resolution=4001):
    blade = data["components"]["blade"]
    structure = blade["structure"]
    table = anchor_table(structure)
    u, xy, perimeter = contour(data, eta, linear, resolution)
    webs = {}
    for web in structure["webs"]:
        a, b = web_endpoints(table, web, eta, u, xy, linear)
        webs[web["name"]] = float(np.linalg.norm(b-a))
    # physical TE widths declared next to the anchors: cross-check of the resolved nd_arc
    declared = {a["name"]: (a["width"], a["width"].get("defines", [])) for a in table.values() if "width" in a}
    mats = {m["name"]: m for m in data["materials"]}
    layers = []
    for layer in structure["layers"]:
        h = interp(layer["thickness"], eta, linear)
        start, end = layer_interval(table, layer, eta, linear)
        extent = webs[layer["web"]] if "web" in layer else perimeter
        width = (end-start)*extent
        check = {}
        if layer["name"] in declared:
            check["declared_width_m"] = interp(declared[layer["name"]][0], eta, linear)
            check["width_residual_m"] = width - check["declared_width_m"]
        mat = mats[layer["material"]]
        theta = interp(layer["fiber_orientation"], eta, linear)
        ea = moduli(mat, theta)*h*width
        layers.append(dict(name=layer["name"], material=layer["material"],
                           width_m=width, thickness_m=h, mass=float(mat["rho"])*h*width,
                           ea_q=ea[0], ea_ex=ea[1], ea_s1=ea[2], nd_arc=[start, end], **check))
    props = structure["elastic_properties"]
    beam_ea = dict(grid=props["stiffness_matrix"]["grid"], values=props["stiffness_matrix"]["K33"])
    beam_mass = dict(grid=props["inertia_matrix"]["grid"], values=props["inertia_matrix"]["mass"])
    beam = {k: interp(dict(grid=props["stiffness_matrix"]["grid"], values=props["stiffness_matrix"][k]), eta, True)
            for k in ("K44", "K55", "K66")}
    return dict(eta=float(eta), perimeter_m=perimeter, web_heights_m=webs, layers=layers,
                **{k:sum(l[k] for l in layers) for k in ("mass", "ea_q", "ea_ex", "ea_s1")},
                beam_ea=interp(beam_ea, eta, True), beam_mass=interp(beam_mass, eta, True),
                beam_k44=beam["K44"], beam_k55=beam["K55"], beam_k66=beam["K66"])


def blade_length(data):
    """Arc length of the reference axis (m), straight segments between its x, z samples."""
    axis = data["components"]["blade"]["reference_axis"]
    grid = np.asarray(axis["z"]["grid"], dtype=float)
    x = np.interp(grid, axis["x"]["grid"], axis["x"]["values"])
    y = np.interp(grid, axis["y"]["grid"], axis["y"]["values"])
    z = np.asarray(axis["z"]["values"], dtype=float)
    return float(np.sum(np.linalg.norm(np.diff(np.column_stack([x, y, z]), axis=0), axis=1)))


def blade_mass(data, n=101, linear=False):
    """Blade mass (kg): layup mass per length integrated over reference-axis arc length (trapezoid, n stations),
    and the same integral of the declared elastic_properties mass. Returns (layup, declared)."""
    etas = np.linspace(0.0, 1.0, n)
    length = blade_length(data)
    m = np.array([section(data, e, linear)["mass"] for e in etas])
    props = data["components"]["blade"]["structure"]["elastic_properties"]["inertia_matrix"]
    declared = np.interp(etas, props["grid"], props["mass"])
    return float(np.trapezoid(m, etas*length)), float(np.trapezoid(declared, etas*length))


def load(path=DEFAULT_YAML):
    with open(path) as stream:
        data = yaml.safe_load(stream)
    if not str(data.get("windIO_version", "")).startswith("2"):
        raise ValueError(f"{path}: windIO 2.x required (windIO_version={data.get('windIO_version')!r}); "
                         f"use {DEFAULT_YAML.relative_to(ROOT)}")
    return data


def self_check():
    """Analytical checks, no fitted blade reference or mutable test stores."""
    iso = dict(orth=0, E=10e9, nu=0.3)
    for theta in (0, 0.3, np.pi/2):
        np.testing.assert_allclose(moduli(iso, theta), [10e9/0.91, 10e9, 10e9/0.91])
    orth = dict(orth=1, E=[100e9, 10e9], G=[5e9], nu=[0.3])
    np.testing.assert_allclose(moduli(orth, 0), [100e9/0.991, 100e9, 100e9/0.91])
    np.testing.assert_allclose(moduli(orth, np.pi/2), [10e9/0.991, 10e9, 10e9/0.91])
    constant = lambda v: dict(grid=[0., 1.], values=[v, v])
    ref = lambda name, handle: dict(anchor=dict(name=name, handle=handle))
    outer = dict(airfoils=[dict(name='square', spanwise_position=0., weight=[1.]),
                           dict(name='square', spanwise_position=1., weight=[1.])],
                 chord=constant(2.), twist=constant(0.))
    anchors = [dict(name='skin', start_nd_arc=constant(0.), end_nd_arc=constant(1.)),
               dict(name='band', start_nd_arc=ref('skin', 'start_nd_arc'), end_nd_arc=constant(.25),
                    width=dict(grid=[0., 1.], values=[2., 2.], defines=['end_nd_arc'])),
               dict(name='web', start_nd_arc=constant(.125), end_nd_arc=constant(.625))]
    layer = dict(name='skin', material='iso', thickness=constant(.01),
                 start_nd_arc=ref('skin', 'start_nd_arc'), end_nd_arc=ref('skin', 'end_nd_arc'),
                 fiber_orientation=constant(.3))
    web = dict(name='web', start_nd_arc=ref('web', 'start_nd_arc'), end_nd_arc=ref('web', 'end_nd_arc'),
               anchors=[dict(name='web_att', start_nd_arc=constant(0.), end_nd_arc=constant(1.))])
    web_layer = dict(layer, name='web_skin', web='web', start_nd_arc=ref('web_att', 'start_nd_arc'),
                     end_nd_arc=ref('web_att', 'end_nd_arc'))
    band = dict(layer, name='band', start_nd_arc=ref('band', 'start_nd_arc'), end_nd_arc=ref('band', 'end_nd_arc'),
                thickness=constant(.02))
    k = lambda v: dict(grid=[0., 1.], K33=[v, v], K44=[v, v], K55=[v, v], K66=[v, v])
    synthetic = dict(components=dict(blade=dict(outer_shape=outer,
        structure=dict(layers=[layer, web_layer, band], webs=[web], anchors=anchors,
        elastic_properties=dict(stiffness_matrix=k(1.), inertia_matrix=dict(grid=[0., 1.], mass=[1., 1.]))))),
        materials=[dict(iso, name='iso', rho=1000.)],
        airfoils=[dict(name='square', coordinates=dict(x=[1., 1., 0., 0., 1.], y=[0., 1., 1., 0., 0.]))])
    result = section(synthetic, .5)
    # skin 8 m * .01 + web 2 m * .01 + band (.25*8 = 2 m) * .02, rho 1000
    np.testing.assert_allclose([result['perimeter_m'], result['web_heights_m']['web'], result['mass']], [8., 2., 100.+40.])
    np.testing.assert_allclose(result['layers'][2]['width_residual_m'], 0., atol=1e-9)
    np.testing.assert_allclose([result['ea_q'], result['ea_ex']], [1.4e9/.91, 1.4e9])
    print("SELF-CHECK PASS: isotropic rotation invariance; orthotropic 0/90 Q/Ex; square skin+web+anchored band geometry/mass/EA")


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
    layup, declared = blade_mass(data, linear=args.linear)
    print(f"BLADE_MASS layup={layup:.1f} kg declared(elastic_properties)={declared:.1f} kg ratio-1={100*(layup/declared-1):+.3f}%")
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(rows, indent=2)+"\n")


if __name__ == "__main__":
    main()
