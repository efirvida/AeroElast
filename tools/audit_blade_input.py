#!/usr/bin/env python
"""Blade input audit: does the 0.25 m shell mesh carry what the official yaml says?

Spec and pre-registered tolerances: ``odd/tasks/blade-input-audit.md`` (S2-S4). The audit reports
and decides nothing: it writes ``report.json`` and ``stations.csv``, prints a summary, exits 0.

Two sides, kept apart on purpose:

* the **yaml side** is ``YamlBlade`` and the pure helpers around it. They need only PyYAML and
  numpy and import nothing from ``aeroelast`` (S2); spans are interpolated linearly on the yaml's
  own grids;
* the **mesh side** (``build_mesh`` onward) builds the production ``Blade`` mesh and reads what the
  solver is given: node coordinates, the element sets, and the laminates from
  ``Blade.get_element_properties()``.

Conventions (each is cited where it is used)
--------------------------------------------
``nd_arc``  arc length measured on the airfoil from the suction-side trailing edge (0), over the
            leading edge, to the pressure-side trailing edge (1), normalised by the arc length of
            that open polyline (the trailing-edge gap is not part of it). Evidence in the repo
            (the windIO ontology schema text itself is not vendored, see ``convention_note``):
            the airfoil coordinates of ``tests/IEA-15-240-RWT.yaml`` run TE (y>0, suction) -> LE
            -> TE (y<0); the spar caps tagged ``side: suction`` carry nd_arc 0.23-0.40 and the
            ``pressure`` ones 0.58-0.77; ``TE_reinforcement_SS`` has ``start_nd_arc`` only and
            ``TE_reinforcement_PS`` ``end_nd_arc`` only, i.e. 0 and 1 are the two trailing edges.
span grid   the yaml ``grid`` is the fraction of the blade length ``reference_axis.z[-1]``
            (117 m); ``r/R`` in S3 is read as that fraction, not as (r + hub radius)/rotor radius.
mesh frame  z along the span; chord line ``x`` (LE at +x, TE towards -x for zero twist); ``y``
            towards the downwind side for positive twist. Prebend (yaml ``reference_axis.x``,
            negative = upwind at the tip) is ``y``; sweep (``reference_axis.y``) is ``-x``.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import subprocess
import sys
from pathlib import Path

import numpy as np
import yaml

STATIONS = (0.30, 0.60, 0.90, 0.97)

# S4 tolerances, pre-registered
TOL_CHORD_REL = 0.005
TOL_TWIST_DEG = 0.1
TOL_POS_ABS_M = 0.01
TOL_POS_REL_C = 0.005
TOL_AIRFOIL_C = 0.005
TOL_THICK_REL = 0.01
TOL_ANGLE_DEG = 0.1
TOL_INTEGRAL_REL = 0.01
TOL_MATERIAL_REL = 1e-9

ACTIVE_T_MIN = 1.0e-12  # m; a layer thinner than this is "not there"

CONVENTION_NOTE = (
    "windIO ontology schema text is not vendored in the repo (only the instance yaml); nd_arc "
    "direction is inferred from the airfoil coordinate order and the side-tagged spar caps"
)


# =============================================================================================
# yaml side: standalone reader, no aeroelast imports
# =============================================================================================
def interp(grid, values, x) -> float:
    """Linear interpolation on the quantity's own grid, clamped at the ends."""
    return float(np.interp(x, np.asarray(grid, float), np.asarray(values, float)))


def _poly_length(pts: np.ndarray) -> float:
    return float(np.linalg.norm(np.diff(pts, axis=0), axis=1).sum())


def _chord_frame(le: np.ndarray, te: np.ndarray, mirrored: bool):
    """Unit vectors (e1 along LE->TE, e2 towards the suction side) and the chord length."""
    v = te - le
    c = float(np.linalg.norm(v))
    e1 = v / c
    e2 = np.array([e1[1], -e1[0]]) if mirrored else np.array([-e1[1], e1[0]])
    return e1, e2, c


def _to_chord_coords(pts: np.ndarray, le: np.ndarray, e1: np.ndarray, e2: np.ndarray, c: float):
    d = np.asarray(pts, float) - le
    return np.c_[d @ e1, d @ e2] / c


def _resample_by_arc(pts: np.ndarray, n: int) -> np.ndarray:
    s = np.r_[0.0, np.cumsum(np.linalg.norm(np.diff(pts, axis=0), axis=1))]
    s /= s[-1]
    t = np.linspace(0.0, 1.0, n)
    return np.c_[np.interp(t, s, pts[:, 0]), np.interp(t, s, pts[:, 1])]


class YamlBlade:
    """Reads the windIO v1 blade of ``tests/IEA-15-240-RWT.yaml``; linear on the yaml grids."""

    def __init__(self, data: dict):
        blade = data["components"]["blade"]
        self.outer = blade["outer_shape_bem"]
        self.struct = blade["internal_structure_2d_fem"]
        self.elastic = blade.get("elastic_properties_mb", {})
        self._airfoils = {a["name"]: a for a in data["airfoils"]}
        self._materials = {m["name"]: m for m in data["materials"]}
        self._layers = {layer["name"]: layer for layer in self.struct["layers"]}

    @classmethod
    def from_file(cls, path) -> YamlBlade:
        with open(path) as fh:
            return cls(yaml.safe_load(fh))

    @classmethod
    def from_dict(cls, data: dict) -> YamlBlade:
        return cls(data)

    @property
    def length(self) -> float:
        return float(self.outer["reference_axis"]["z"]["values"][-1])

    # -- outer shape ---------------------------------------------------------------------
    def outer_shape(self, frac: float) -> dict:
        o, ra = self.outer, self.outer["reference_axis"]
        at = lambda d: interp(d["grid"], d["values"], frac)  # noqa: E731
        return {
            "chord": at(o["chord"]),
            "twist_deg": math.degrees(at(o["twist"])),
            "pitch_axis": at(o["pitch_axis"]),
            "prebend": at(ra["x"]),
            "sweep": at(ra["y"]),
            "z": at(ra["z"]),
        }

    # -- airfoils ------------------------------------------------------------------------
    def _airfoil_chord_coords(self, name: str) -> np.ndarray:
        co = self._airfoils[name]["coordinates"]
        p = np.c_[np.asarray(co["x"], float), np.asarray(co["y"], float)]
        le = p[int(np.argmin(p[:, 0]))]
        te = 0.5 * (p[0] + p[-1])
        e1, e2, c = _chord_frame(le, te, mirrored=False)
        return _to_chord_coords(p, le, e1, e2, c)

    def airfoil_at(self, frac: float):
        """(polygon in chord units TE-SS -> LE -> TE-PS, defined, scheme).

        ``airfoil_position`` gives labels at grid points only. Between two *identical* labels the
        airfoil is that airfoil (defined). Between two different labels the yaml states no
        blending rule: the polygon returned is a linear blend of both, resampled on normalised
        arc length, and ``defined`` is False.
        """
        g, lab = self.outer["airfoil_position"]["grid"], self.outer["airfoil_position"]["labels"]
        i = int(np.clip(np.searchsorted(g, frac, side="right") - 1, 0, len(g) - 2))
        a, b = lab[i], lab[i + 1]
        if a == b:
            return self._airfoil_chord_coords(a), True, f"{a} (both neighbours)"
        w = (frac - g[i]) / (g[i + 1] - g[i])
        pa = _resample_by_arc(self._airfoil_chord_coords(a), 400)
        pb = _resample_by_arc(self._airfoil_chord_coords(b), 400)
        scheme = f"linear blend {a} -> {b}, weight {w:.3f}, arc-length resampled (yaml states none)"
        return (1 - w) * pa + w * pb, False, scheme

    def perimeter(self, frac: float) -> float:
        """Arc length of the open airfoil polyline in metres (chord x polygon length)."""
        poly, _, _ = self.airfoil_at(frac)
        return self.outer_shape(frac)["chord"] * _poly_length(poly)

    # -- layers ----------------------------------------------------------------------------
    def layer_state(self, name: str, frac: float) -> dict:
        layer = self._layers[name]
        at = lambda k: interp(layer[k]["grid"], layer[k]["values"], frac) if k in layer else None  # noqa: E731
        fo = at("fiber_orientation")
        return {
            "name": name,
            "material": layer["material"],
            "thickness": at("thickness"),
            "start": at("start_nd_arc"),
            "end": at("end_nd_arc"),
            "width": at("width"),
            "side": layer.get("side"),
            "angle_deg": math.degrees(fo) if fo is not None else 0.0,
            "web": layer.get("web"),
        }

    def layer_bounds(self, name: str, frac: float, perimeter: float):
        """(start, end) nd_arc. A missing bound is ``start + width/perimeter`` (suction side) or
        ``end - width/perimeter`` (pressure side), as the yaml's TE reinforcements require."""
        st = self.layer_state(name, frac)
        s, e, w = st["start"], st["end"], st["width"]
        if s is None and e is not None and w is not None:
            s = e - w / perimeter
        if e is None and s is not None and w is not None:
            e = s + w / perimeter
        return s, e

    def bay_mean_thickness(self, name: str, z0: float, z1: float) -> float:
        """Mean of the yaml thickness of a layer over the span interval [z0, z1] (metres).

        Context only: the mesh holds one laminate per spanwise bay, the yaml a point value.
        """
        fr = np.linspace(z0, z1, 201) / self.length
        layer = self._layers[name]["thickness"]
        return float(np.mean(np.interp(fr, layer["grid"], layer["values"])))

    def _skin_layers(self):
        return [name for name, layer in self._layers.items() if "web" not in layer]

    def active_layers(self, frac: float, a: float, perimeter: float) -> list[dict]:
        """Shell layers active at (frac, nd_arc a), in yaml order (outer to inner)."""
        out = []
        for name in self._skin_layers():
            st = self.layer_state(name, frac)
            s, e = self.layer_bounds(name, frac, perimeter)
            if st["thickness"] <= ACTIVE_T_MIN or s is None or e is None:
                continue
            if s <= a <= e:
                out.append({**st, "start": s, "end": e})
        return out

    def all_bounds(self, frac: float, perimeter: float) -> list[float]:
        """Every start/end nd_arc of the layers that exist at ``frac`` (for the boundary rule)."""
        bounds = []
        for name in self._skin_layers():
            if self.layer_state(name, frac)["thickness"] <= ACTIVE_T_MIN:
                continue
            bounds += [b for b in self.layer_bounds(name, frac, perimeter) if b is not None]
        return bounds

    # -- webs ------------------------------------------------------------------------------
    def web_state(self, k: int, frac: float) -> dict:
        w = self.struct["webs"][k]
        return {
            "name": w["name"],
            "start": interp(w["start_nd_arc"]["grid"], w["start_nd_arc"]["values"], frac),
            "end": interp(w["end_nd_arc"]["grid"], w["end_nd_arc"]["values"], frac),
        }

    def web_layers(self, k: int, frac: float) -> list[dict]:
        wname = self.struct["webs"][k]["name"]
        out = []
        for name, layer in self._layers.items():
            if layer.get("web") == wname:
                st = self.layer_state(name, frac)
                if st["thickness"] > ACTIVE_T_MIN:
                    out.append(st)
        return out

    # -- materials -------------------------------------------------------------------------
    def material(self, name: str) -> dict:
        m = self._materials[name]
        if m.get("orth", 0) in (0, 0.0, False):
            e, g, nu = m["E"], m["G"], m["nu"]
            return {"e1": e, "e2": e, "e3": e, "g12": g, "g13": g, "g23": g,
                    "nu12": nu, "nu13": nu, "nu23": nu, "rho": m["rho"]}
        e, g, nu = m["E"], m["G"], m["nu"]  # windIO order: xx,yy,zz / xy,xz,yz / xy,xz,yz
        return {"e1": e[0], "e2": e[1], "e3": e[2], "g12": g[0], "g13": g[1], "g23": g[2],
                "nu12": nu[0], "nu13": nu[1], "nu23": nu[2], "rho": m["rho"]}

    # -- deck context ----------------------------------------------------------------------
    def deck_m11_k11_k33(self, frac: float):
        """(M11 kg/m, K11 N, K33 N) of the yaml 6x6 deck, or None.

        ``values`` hold the 21 upper-triangle entries row by row, so index 0 is K11/M11 and
        index 11 is K33. z is the span axis, hence K33 is the axial stiffness EA; K11 is a shear
        stiffness. S4 asks for "M11/K11" and K11 is printed literally; K33 is added as the
        quantity that F5's axial integral can be compared with.
        """
        six = self.elastic.get("six_x_six")
        if not six:
            return None
        k, m = six["stiff_matrix"], six["inertia_matrix"]
        kv, mv = np.asarray(k["values"], float), np.asarray(m["values"], float)
        return (interp(m["grid"], mv[:, 0], frac), interp(k["grid"], kv[:, 0], frac),
                interp(k["grid"], kv[:, 11], frac))


# -- pure helpers --------------------------------------------------------------------------
def ring_arcs(points: np.ndarray):
    """Cumulative arc length over total (nd_arc of each vertex) and the total length."""
    seg = np.linalg.norm(np.diff(np.asarray(points, float), axis=0), axis=1)
    cum = np.r_[0.0, np.cumsum(seg)]
    return cum / cum[-1], float(cum[-1])


def axial_modulus(e1, e2, g12, nu12, theta_deg) -> float:
    """Uniaxial modulus along the span of a ply whose fibre is ``theta`` off the span axis."""
    th = math.radians(theta_deg)
    c, s = math.cos(th), math.sin(th)
    compliance = c**4 / e1 + (1.0 / g12 - 2.0 * nu12 / e1) * s * s * c * c + s**4 / e2
    return 1.0 / compliance


def section_integrals(elements) -> tuple[float, float]:
    """(mass, axial) per unit span: sum of rho*t*ds and Ex*t*ds over (plies, ds) elements."""
    mass = axial = 0.0
    for plies, ds in elements:
        for p in plies:
            mass += p["rho"] * p["t"] * ds
            axial += axial_modulus(p["e1"], p["e2"], p["g12"], p["nu12"], p["angle_deg"]) * p["t"] * ds
    return mass, axial


def normal_distance(points: np.ndarray, polyline: np.ndarray):
    """(max, mean) unsigned distance from each point to the polyline."""
    pts = np.asarray(points, float)
    a, b = polyline[:-1], polyline[1:]
    ab = b - a
    denom = np.maximum((ab**2).sum(1), 1e-300)
    d = np.empty(len(pts))
    for i, p in enumerate(pts):
        t = np.clip(((p - a) * ab).sum(1) / denom, 0.0, 1.0)
        d[i] = np.linalg.norm(a + t[:, None] * ab - p, axis=1).min()
    return float(d.max()), float(d.mean())


def classify(residual, tol, *, boundary=False, gated=True) -> str:
    if not gated:
        return "not-gated"
    if boundary:
        return "boundary"
    return "pass" if abs(residual) <= tol else "fail"


def is_boundary(a, bounds, *, perimeter, element_size) -> bool:
    """Centroid within half an element (metres) of any layer start/end."""
    half = 0.5 * element_size / perimeter
    return any(abs(a - b) <= half + 1e-12 for b in bounds)


# -- mesh-frame forward model and ring measurement ----------------------------------------
def profile_to_mesh_ring(prof, chord, pitch_axis, twist_deg=0.0, prebend=0.0, sweep=0.0):
    """The documented mesh convention applied to a unit-chord profile (used to test the reader).

    x_m = -(x - pitch_axis)*chord, y_m = y*chord, rotated by -twist, shifted by (-sweep, prebend).
    """
    p = np.asarray(prof, float)
    x = -(p[:, 0] - pitch_axis) * chord
    y = p[:, 1] * chord
    th = math.radians(-twist_deg)
    xr = math.cos(th) * x - math.sin(th) * y
    yr = math.sin(th) * x + math.cos(th) * y
    return np.c_[xr - sweep, yr + prebend]


def _ring_frame(ring: np.ndarray):
    te = 0.5 * (ring[0] + ring[-1])
    le = ring[int(np.argmax(np.linalg.norm(ring - te, axis=1)))]
    e1, e2, c = _chord_frame(le, te, mirrored=True)
    return le, te, e1, e2, c


def orient_ring(ring: np.ndarray) -> np.ndarray:
    """Return the ring starting at the suction-side trailing edge (w > 0 on the first half)."""
    le, _, e1, e2, c = _ring_frame(ring)
    w = _to_chord_coords(ring, le, e1, e2, c)[:, 1]
    k = len(ring) // 2
    return ring if w[:k].mean() >= w[k:].mean() else ring[::-1]


def measure_outer(ring: np.ndarray, pitch_axis: float, ref_point=None) -> dict:
    """Chord, twist and pitch-axis point of an oriented ring.

    ``twist_deg`` is the angle of LE->TE to the rotor plane, positive when the TE is on the
    +y (downwind) side. ``px, py`` is where the yaml pitch axis sits on the measured chord line.
    ``pa_from_ref`` projects ``ref_point`` (the yaml reference axis) on the chord line.
    """
    le, te, e1, _, c = _ring_frame(ring)
    p = le + pitch_axis * (te - le)
    out = {"chord": c, "twist_deg": math.degrees(math.atan2(e1[1], -e1[0])),
           "px": float(p[0]), "py": float(p[1]), "le": le, "te": te}
    if ref_point is not None:
        out["pa_from_ref"] = float((np.asarray(ref_point) - le) @ e1 / c)
    return out


def parse_args(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--yaml", default="tests/IEA-15-240-RWT.yaml")
    ap.add_argument("--element-size", type=float, default=0.25)
    ap.add_argument("--output-dir", default="logs/check/input-audit")
    ap.add_argument("--mesh-cache", default=None,
                    help="pickle of the extracted mesh data (development shortcut; built if absent)")
    return ap.parse_args(argv)


# =============================================================================================
# mesh side: the production path
# =============================================================================================
def build_mesh(yaml_path: str, element_size: float) -> dict:
    """Build the mesh with ``Blade`` and extract what the audit reads (plain arrays and dicts).

    Reads: node coordinates (``MeshModel.nodes``), element node ids and element sets
    (``MeshModel.elements`` / ``element_sets``), the laminate per element set from
    ``Blade.get_element_properties()`` (what is handed to the solver) and the material names from
    the NuMAD section layup (``Blade._numad_mesh['sections']``), which the Rust laminate lacks.
    """
    from aeroelast.models.blade.model import Blade  # noqa: PLC0415

    blade = Blade(yaml_path, element_size=element_size)
    blade.generate_mesh()
    m = blade.mesh
    nodes = np.array([n.coords for n in m.nodes], float)
    idx = {n.id: k for k, n in enumerate(m.nodes)}
    elements = [[idx[n.id] for n in e.nodes] for e in m.elements]
    eid_pos = {e.id: k for k, e in enumerate(m.elements)}
    sets = {name: [eid_pos[e.id] for e in s.elements] for name, s in m.element_sets.items()}
    layups = {s["elementSet"]: s["layup"] for s in blade._numad_mesh["sections"]}
    props = blade.get_element_properties()
    stacks = {}
    for name, lam in props.items():
        layup = layups.get(name)
        if layup is None:
            continue
        if isinstance(lam, dict):  # isotropic single-ply path
            e, nu = lam["e"], lam["nu"]
            plies = [{"mat": layup[0][0], "t": lam["thickness"], "angle_deg": 0.0, "z_bottom": -lam["thickness"] / 2,
                      "z_top": lam["thickness"] / 2, "rho": lam["rho"], "e1": e, "e2": e, "e3": e,
                      "g12": e / (2 * (1 + nu)), "g13": e / (2 * (1 + nu)), "g23": e / (2 * (1 + nu)),
                      "nu12": nu, "nu23": nu, "nu31": nu, "iso_dict": True}]
        else:
            plies = []
            for (mat, _t, _a), d in zip(layup, lam.plies_data(), strict=True):
                mt = d["material"]
                plies.append({"mat": mat, "t": d["thickness"], "angle_deg": d["angle"],
                              "z_bottom": d["z_bottom"], "z_top": d["z_top"], "rho": mt["rho"],
                              "e1": mt["e1"], "e2": mt["e2"], "e3": mt["e3"], "g12": mt["g12"],
                              "g13": mt["g13"], "g23": mt["g23"], "nu12": mt["nu12"],
                              "nu23": mt["nu23"], "nu31": mt["nu31"]})
        stacks[name] = plies
    return {"nodes": nodes, "elements": elements, "sets": sets, "stacks": stacks,
            "element_size": element_size}


def add_set_zrange(mesh: dict) -> dict:
    """Span extent (z min, z max) of every element set: the bay a laminate is held over."""
    nodes, els = mesh["nodes"], mesh["elements"]
    mesh["set_zrange"] = {}
    for name, members in mesh["sets"].items():
        if members:
            z = nodes[np.unique(np.concatenate([els[e] for e in members])), 2]
            mesh["set_zrange"][name] = (float(z.min()), float(z.max()))
    return mesh


def _slice(mesh: dict, z: float, names) -> list[dict]:
    """Segments where the planar cut z=const crosses the quads/tris of the given sets."""
    nodes, els = mesh["nodes"], mesh["elements"]
    zt = z + 1e-9 if np.any(np.abs(nodes[:, 2] - z) < 1e-9) else z
    segs = []
    for name in names:
        for ei in mesh["sets"][name]:
            nd = els[ei]
            zz = nodes[nd, 2]
            if zz.min() > zt or zz.max() < zt:
                continue
            pts, keys = [], []
            for k in range(len(nd)):
                i, j = nd[k], nd[(k + 1) % len(nd)]
                zi, zj = nodes[i, 2], nodes[j, 2]
                if (zi - zt) * (zj - zt) < 0:
                    a, b = (i, j) if i < j else (j, i)
                    t = (zt - nodes[a, 2]) / (nodes[b, 2] - nodes[a, 2])
                    pts.append(nodes[a, :2] + t * (nodes[b, :2] - nodes[a, :2]))
                    keys.append((a, b))
            if len(keys) != 2:
                continue
            q = nodes[nd[:3]]
            n = np.cross(q[1] - q[0], q[2] - q[0])
            segs.append({"set": name, "elem": ei, "keys": keys, "pts": pts,
                         "normal": n / (np.linalg.norm(n) + 1e-300)})
    return segs


def _chain(segs: list[dict]):
    """Order the segments into one open path. Returns (points, ordered segments) or raises."""
    adj: dict = {}
    for i, s in enumerate(segs):
        for k in s["keys"]:
            adj.setdefault(k, []).append(i)
    ends = [k for k, v in adj.items() if len(v) == 1]
    if len(ends) != 2 or any(len(v) > 2 for v in adj.values()):
        raise ValueError(f"cut is not a simple open path ({len(ends)} endpoints, {len(segs)} segments)")
    key, used, order, pts = ends[0], set(), [], []
    point = {}
    for s in segs:
        for k, p in zip(s["keys"], s["pts"], strict=True):
            point[k] = p
    pts.append(point[key])
    while True:
        nxt = [i for i in adj[key] if i not in used]
        if not nxt:
            break
        i = nxt[0]
        used.add(i)
        s = segs[i]
        other = s["keys"][1] if s["keys"][0] == key else s["keys"][0]
        order.append({**s, "a": pts[-1], "b": point[other]})
        pts.append(point[other])
        key = other
    return np.array(pts), order


def _signed_offset(ring_uw: np.ndarray, poly_uw: np.ndarray):
    """(mean, max) signed distance of the ring from the yaml polygon, outside positive."""
    from shapely.geometry import LineString, Point, Polygon  # noqa: PLC0415

    poly = Polygon(poly_uw)
    line = LineString(poly_uw)
    d = np.array([line.distance(Point(p)) * (1.0 if not poly.contains(Point(p)) else -1.0)
                  for p in ring_uw])
    return float(d.mean()), float(d[np.argmax(np.abs(d))])


# =============================================================================================
# the audit
# =============================================================================================
def _check(station, field, quantity, mesh_v, yaml_v, residual, tol, status, location=None, note=None):
    return {"station": station, "field": field, "quantity": quantity, "mesh": mesh_v,
            "yaml": yaml_v, "residual": residual, "tolerance": tol, "status": status,
            "location": location or {}, "note": note}


def _relabs(got, want):
    return abs(got - want) / abs(want) if want else abs(got - want)


def _rel(a, b):
    return (a - b) / b if b else float("inf") if a else 0.0


def audit_station(yb: YamlBlade, mesh: dict, frac: float) -> dict:
    es = mesh["element_size"]
    o = yb.outer_shape(frac)
    z = o["z"]
    outer_names = [n for n in mesh["sets"] if not n.endswith("_SW") and not n.startswith("all")]
    web_names = [n for n in mesh["sets"] if n.endswith("_SW")]
    checks: list[dict] = []
    st = f"{frac:.2f}"
    info: dict = {"r_over_R": frac, "z": z, "yaml_outer_shape": o}

    cut = _slice(mesh, z, outer_names)
    # the flat-back face closes the loop but is not on the yaml's open airfoil arc: kept apart
    flats = [s for s in cut if s["set"].endswith("TE_FLAT")]
    try:
        pts, segs = _chain([s for s in cut if not s["set"].endswith("TE_FLAT")])
    except ValueError as exc:
        info["error"] = f"outer ring: {exc}"
        return {"info": info, "checks": checks}
    flip = orient_ring(pts) is not pts
    if flip:
        pts, segs = pts[::-1], [{**s, "a": s["b"], "b": s["a"]} for s in segs[::-1]]
    ring = pts
    arcs, p_mesh = ring_arcs(ring)
    info.update(n_ring_points=len(ring), n_ring_elements=len(segs), mesh_perimeter=p_mesh,
                ring_flipped_to_start_at_suction_te=bool(flip))

    # ---- F1 outer shape ----------------------------------------------------------------------
    ref = np.array([-o["sweep"], o["prebend"]])
    m = measure_outer(ring, o["pitch_axis"], ref_point=ref)
    pos_tol = max(TOL_POS_ABS_M, TOL_POS_REL_C * o["chord"])
    pos_tol_strict = min(TOL_POS_ABS_M, TOL_POS_REL_C * o["chord"])
    for q, mv, yv, tol, rel in (
        ("chord", m["chord"], o["chord"], TOL_CHORD_REL, True),
        ("twist_deg", m["twist_deg"], o["twist_deg"], TOL_TWIST_DEG, False),
        ("pitch_axis_fraction (yaml reference point projected on mesh chord line)",
         m["pa_from_ref"], o["pitch_axis"], TOL_POS_REL_C, False),
        ("pitch_point_x (= -sweep)", m["px"], -o["sweep"], pos_tol, False),
        ("pitch_point_y (= prebend)", m["py"], o["prebend"], pos_tol, False),
    ):
        res = _rel(mv, yv) if rel else mv - yv
        c = _check(st, "F1", q, mv, yv, res, tol, classify(res, tol))
        if q.startswith("pitch_point"):
            c["note"] = f"tolerance max(1 cm, 0.5% c) = {pos_tol:.4f} m; strict min = {pos_tol_strict:.4f} m"
            c["pass_strict"] = abs(res) <= pos_tol_strict
        checks.append(c)

    # ---- F2 airfoil shape --------------------------------------------------------------------
    le, te, e1, e2, c_m = _ring_frame(ring)
    ring_uw = _to_chord_coords(ring, le, e1, e2, c_m)
    poly, defined, scheme = yb.airfoil_at(frac)
    dmax, dmean = normal_distance(ring_uw, poly)
    checks.append(_check(st, "F2", "max normal distance / chord", dmax, 0.0, dmax, TOL_AIRFOIL_C,
                         classify(dmax, TOL_AIRFOIL_C, gated=defined), note=scheme))
    info["airfoil_scheme"] = scheme
    info["airfoil_defined_by_yaml"] = defined
    info["airfoil_mean_distance_over_chord"] = dmean

    # ---- element rows ------------------------------------------------------------------------
    p_yaml = yb.perimeter(frac)
    bounds = yb.all_bounds(frac, p_yaml)
    info["yaml_perimeter"] = p_yaml
    yaml_rows, mesh_rows = [], []
    # nd_arc of each element centroid = cumulative arc up to the cut-segment midpoint
    cum = np.r_[0.0, np.cumsum([np.linalg.norm(s["b"] - s["a"]) for s in segs])]
    f3_elems = [(s, float((cum[i] + 0.5 * np.linalg.norm(s["b"] - s["a"])) / cum[-1]))
                for i, s in enumerate(segs)]

    # ---- F3 laminate stack, F6 materials ------------------------------------------------------
    mat_worst: dict = {}
    for s, a_c in f3_elems:
        stack = mesh["stacks"].get(s["set"])
        want = yb.active_layers(frac, a_c, p_yaml)
        bnd = is_boundary(a_c, bounds, perimeter=p_yaml, element_size=es)
        loc = {"set": s["set"], "elem": int(s["elem"]), "nd_arc": a_c}
        if stack is None:
            checks.append(_check(st, "F3", "stack present", None, [w["name"] for w in want], 1.0, 0.0,
                                 classify(1.0, 0.0, boundary=bnd), loc))
            continue
        got_mats = [p["mat"] for p in stack]
        want_mats = [w["material"] for w in want]
        same = got_mats == want_mats
        checks.append(_check(st, "F3", "material and order", got_mats, want_mats, 0.0 if same else 1.0, 0.0,
                             classify(0.0 if same else 1.0, 0.0, boundary=bnd), loc,
                             note=None if same else "yaml layers: " + ",".join(w["name"] for w in want)))
        tt_got, tt_want = sum(p["t"] for p in stack), sum(w["thickness"] for w in want)
        r = _rel(tt_got, tt_want)
        checks.append(_check(st, "F3", "total thickness", tt_got, tt_want, r, TOL_THICK_REL,
                             classify(r, TOL_THICK_REL, boundary=bnd), loc))
        if same:
            for k, (p, w) in enumerate(zip(stack, want, strict=True)):
                r = _rel(p["t"], w["thickness"])
                checks.append(_check(st, "F3", f"ply thickness [{k}] {w['name']}", p["t"], w["thickness"], r,
                                     TOL_THICK_REL, classify(r, TOL_THICK_REL, boundary=bnd), loc))
                ra = p["angle_deg"] - w["angle_deg"]
                checks.append(_check(st, "F3", f"fibre angle [{k}] {w['name']}", p["angle_deg"],
                                     w["angle_deg"], ra, TOL_ANGLE_DEG,
                                     classify(ra, TOL_ANGLE_DEG, boundary=bnd), loc))
        z0, z1 = mesh["set_zrange"][s["set"]]
        tt_bay = sum(yb.bay_mean_thickness(w["name"], z0, z1) for w in want)
        checks.append(_check(st, "F3", "total thickness vs yaml bay-mean (context)", tt_got, tt_bay, _rel(tt_got, tt_bay),
                             None, "not-gated", loc, note=f"bay z {z0:.3f}-{z1:.3f}"))
        for p in stack:
            ym = yb.material(p["mat"])
            for cname in ("e1", "e2", "g12", "nu12", "rho"):
                r = _relabs(p[cname], ym[cname])
                key = (p["mat"], cname)
                if key not in mat_worst or r > mat_worst[key][0]:
                    mat_worst[key] = (r, p[cname], ym[cname], s["set"])
        yaml_rows.append(([{**yb.material(w["material"]), "t": w["thickness"], "angle_deg": w["angle_deg"]}
                           for w in want], float(np.linalg.norm(s["b"] - s["a"]))))
        mesh_rows.append(([{k: p[k] for k in ("rho", "t", "e1", "e2", "g12", "nu12", "angle_deg")} for p in stack],
                          float(np.linalg.norm(s["b"] - s["a"]))))
    for (mat, cname), (r, got, want, where) in sorted(mat_worst.items()):
        checks.append(_check(st, "F6", f"{mat}.{cname}", got, want, r, TOL_MATERIAL_REL,
                             classify(r, TOL_MATERIAL_REL), {"set": where}))
    extras = _material_extras(yb, mesh, outer_names, f3_elems)
    checks += [_check(st, "F6", f"{mat}.{q} (not in S4 list)", g, w, _relabs(g, w), TOL_MATERIAL_REL,
                      "not-gated", {"set": where}, note="out-of-plane constant, reported")
               for (mat, q), (g, w, where) in extras.items()]

    # ---- F5 integrals ------------------------------------------------------------------------
    mm, ma = section_integrals(mesh_rows)
    flat_rows = [([{k: p[k] for k in ("rho", "t", "e1", "e2", "g12", "nu12", "angle_deg")} for p in mesh["stacks"][s["set"]]],
                  float(np.linalg.norm(s["pts"][1] - s["pts"][0])))
                 for s in flats if s["set"] in mesh["stacks"]]
    fm, fa = section_integrals(flat_rows)
    checks.append(_check(st, "F5", "TE flat-back face, mass [kg/m] (off the yaml arc, excluded above)", fm, 0.0, fm / mm if mm else None,
                         None, "not-gated", note=f"axial {fa:.4e} N; {len(flats)} elements; fraction of ring mass in residual"))
    ym_, ya = section_integrals(yaml_rows)
    # the yaml side integrates on the yaml perimeter: scale each row's ds by perimeter ratio
    scale = p_yaml / p_mesh
    ym_, ya = ym_ * scale, ya * scale
    for q, g, w in (("mass per unit span sum(rho t ds) [kg/m]", mm, ym_),
                    ("axial sum(Ex t ds) [N]", ma, ya)):
        r = _rel(g, w)
        checks.append(_check(st, "F5", q, g, w, r, TOL_INTEGRAL_REL, classify(r, TOL_INTEGRAL_REL)))
    deck = yb.deck_m11_k11_k33(frac)
    if deck:
        for q, got, want in (("deck M11 [kg/m] (context: deck includes webs; mesh = ring only)", mm, deck[0]),
                             ("deck K33 axial [N] (context: deck includes webs; mesh = ring only)", ma, deck[2]),
                             ("deck K11 [N] (literal S4 context; a shear term, not comparable)", None, deck[1])):
            checks.append(_check(st, "F5", q, got, want, _rel(got, want) if got is not None else None, None,
                                 "not-gated", note="mesh column = ring integral, yaml column = deck value"))
    # context: the same yaml integral with each layer's thickness averaged over the mesh bay
    ctx = _bay_mean_integrals(yb, mesh, segs, f3_elems, frac, p_yaml, p_mesh)
    for q, v in ctx.items():
        base = ym_ if q.startswith("mass") else ya
        mv = mm if q.startswith("mass") else ma
        checks.append(_check(st, "F5", q + " (yaml thickness averaged over the mesh bay; context)", mv, v, _rel(mv, v),
                             None, "not-gated", note=f"point-value yaml integral is {base:.6g}"))

    # ---- F4 webs -----------------------------------------------------------------------------
    checks += _audit_webs(yb, mesh, web_names, z, frac, ring, arcs, p_mesh, es, st, c_m)

    # ---- F7 reference surface ----------------------------------------------------------------
    mean_off, max_off = _signed_offset(ring_uw, poly)
    cen = ring.mean(axis=0)
    out_dir = np.array([np.sign(np.dot(s["normal"][:2], 0.5 * (s["a"] + s["b"]) - cen)) for s in segs])
    sample = mesh["stacks"].get(f3_elems[len(f3_elems) // 2][0]["set"], [])
    checks.append(_check(st, "F7", "ring offset from yaml OML / chord (signed, outside +)", mean_off, 0.0, mean_off,
                         None, "not-gated", note=f"largest signed offset {max_off:+.5f}; ring nodes sit on the OML if ~0"))
    zb = min((p["z_bottom"] for p in sample), default=None)
    zt_ = max((p["z_top"] for p in sample), default=None)
    checks.append(_check(st, "F7", "laminate z range about the node surface [m]", [zb, zt_], None, None, None,
                         "not-gated", note="symmetric about 0 means centred on the nodes"))
    frac_out = float((out_dir > 0).mean())
    checks.append(_check(st, "F7", "fraction of ring elements whose normal points outward", frac_out, None, None, None,
                         "not-gated", note="ply 0 (yaml UV, outermost layer) sits at z_bottom, i.e. on the -normal side"))
    return {"info": info, "checks": checks}


def _bay_mean_integrals(yb, mesh, segs, f3_elems, frac, p_yaml, p_mesh):
    """Mass and axial integrals of the yaml stack with bay-averaged layer thickness (context)."""
    rows = []
    for s, a_c in f3_elems:
        z0, z1 = mesh["set_zrange"][s["set"]]
        plies = []
        for w in yb.active_layers(frac, a_c, p_yaml):
            plies.append({**yb.material(w["material"]), "t": yb.bay_mean_thickness(w["name"], z0, z1),
                          "angle_deg": w["angle_deg"]})
        rows.append((plies, float(np.linalg.norm(s["b"] - s["a"])) * p_yaml / p_mesh))
    mass, axial = section_integrals(rows)
    return {"mass per unit span [kg/m]": mass, "axial [N]": axial}


def _material_extras(yb, mesh, names, f3_elems):
    out = {}
    for s, _ in f3_elems:
        for p in mesh["stacks"].get(s["set"], []):
            ym = yb.material(p["mat"])
            # same names on both sides (windIO: G = xy,xz,yz; nu = xy,xz,yz with x the fibre axis)
            for q, got, want in (("e3", p["e3"], ym["e3"]), ("g13", p["g13"], ym["g13"]),
                                 ("g23", p["g23"], ym["g23"]), ("nu23", p["nu23"], ym["nu23"])):
                r = _relabs(got, want)
                if (p["mat"], q) not in out or r > _relabs(*out[(p["mat"], q)][:2]):
                    out[(p["mat"], q)] = (got, want, s["set"])
    return out


def _audit_webs(yb, mesh, web_names, z, frac, ring, arcs, p_mesh, es, st, c_m):
    checks = []
    groups: dict = {}
    for n in web_names:
        groups.setdefault(n.split("_")[0], []).append(n)
    nweb = len(yb.struct["webs"])
    mesh_webs = []
    for gid, names in sorted(groups.items()):
        cut = _slice(mesh, z, names)
        if not cut:
            continue  # no element of this web at the station: counted below
        try:
            wpts, wsegs = _chain(cut)
        except ValueError as exc:
            checks.append(_check(st, "F4", f"web group {gid} chain", str(exc), None, None, None, "not-gated"))
            continue
        ends = []
        for p in (wpts[0], wpts[-1]):
            d = np.linalg.norm(ring - p, axis=1)
            k = int(np.argmin(d))
            # refine on the two adjacent ring segments
            best = (d[k], arcs[k])
            for k0, k1 in ((max(k - 1, 0), k), (k, min(k + 1, len(ring) - 1))):
                if k0 == k1:
                    continue
                a, b = ring[k0], ring[k1]
                t = float(np.clip((p - a) @ (b - a) / ((b - a) @ (b - a)), 0.0, 1.0))
                dist = float(np.linalg.norm(a + t * (b - a) - p))
                if dist < best[0]:
                    best = (dist, arcs[k0] + t * (arcs[k1] - arcs[k0]))
            ends.append(best[1])
        mesh_webs.append({"gid": gid, "arcs": sorted(ends), "set": wsegs[0]["set"]})
    yaml_active = [k for k in range(nweb) if yb.web_layers(k, frac)]
    same_count = len(mesh_webs) == len(yaml_active)
    checks.append(_check(st, "F4", "number of webs present", len(mesh_webs), len(yaml_active),
                         float(len(mesh_webs) - len(yaml_active)), 0.0,
                         classify(0.0 if same_count else 1.0, 0.0)))
    if not mesh_webs or not same_count:
        return checks
    # assign mesh webs to yaml webs by smallest total attachment error (independent of numad names)
    import itertools  # noqa: PLC0415

    best_assign, best_cost = None, float("inf")
    for perm in itertools.permutations(range(nweb), len(mesh_webs)):
        cost = sum(abs(mw["arcs"][0] - yb.web_state(k, frac)["start"]) + abs(mw["arcs"][1] - yb.web_state(k, frac)["end"])
                   for mw, k in zip(mesh_webs, perm, strict=True))
        if cost < best_cost:
            best_assign, best_cost = perm, cost
    tol_arc = 0.5 * es / p_mesh
    for mw, k in zip(mesh_webs, best_assign, strict=True):
        ws = yb.web_state(k, frac)
        loc = {"mesh_web_group": mw["gid"], "yaml_web": ws["name"]}
        for side, got, want in (("suction", mw["arcs"][0], ws["start"]), ("pressure", mw["arcs"][1], ws["end"])):
            r_m = (got - want) * p_mesh
            checks.append(_check(st, "F4", f"web attachment nd_arc ({side}) [m along arc]", got, want, r_m, 0.5 * es,
                                 classify(r_m, 0.5 * es), loc, note=f"tolerance es/2 = {0.5 * es} m = {tol_arc:.4f} nd_arc"))
        stack = mesh["stacks"].get(mw["set"]) or []
        want_l = yb.web_layers(k, frac)
        same = [p["mat"] for p in stack] == [w["material"] for w in want_l]
        checks.append(_check(st, "F4", "web material and order", [p["mat"] for p in stack],
                             [w["material"] for w in want_l], 0.0 if same else 1.0, 0.0,
                             classify(0.0 if same else 1.0, 0.0), loc))
        tg, tw = sum(p["t"] for p in stack), sum(w["thickness"] for w in want_l)
        r = _rel(tg, tw)
        checks.append(_check(st, "F4", "web total thickness", tg, tw, r, TOL_THICK_REL, classify(r, TOL_THICK_REL), loc))
        if same:
            for j, (p, w) in enumerate(zip(stack, want_l, strict=True)):
                r = _rel(p["t"], w["thickness"])
                checks.append(_check(st, "F4", f"web ply thickness [{j}] {w['name']}", p["t"], w["thickness"], r,
                                     TOL_THICK_REL, classify(r, TOL_THICK_REL), loc))
                ra = p["angle_deg"] - w["angle_deg"]
                checks.append(_check(st, "F4", f"web fibre angle [{j}] {w['name']}", p["angle_deg"], w["angle_deg"], ra,
                                     TOL_ANGLE_DEG, classify(ra, TOL_ANGLE_DEG), loc))
    return checks


# =============================================================================================
# reporting
# =============================================================================================
def summarize(checks: list[dict]) -> dict:
    """Per (station, field): status counts and the worst gated or boundary residual/tol."""
    out: dict = {}
    for c in checks:
        e = out.setdefault((c["station"], c["field"]),
                           {"pass": 0, "fail": 0, "boundary": 0, "not-gated": 0, "worst": None})
        e[c["status"]] += 1
        if c["status"] in ("pass", "fail") and c["tolerance"]:
            ratio = abs(c["residual"]) / c["tolerance"]
            if e["worst"] is None or ratio > e["worst"]["ratio"]:
                e["worst"] = {"ratio": ratio, "quantity": c["quantity"], "residual": c["residual"],
                              "tolerance": c["tolerance"], "location": c["location"]}
        elif c["status"] in ("pass", "fail"):
            if e["worst"] is None:
                e["worst"] = {"ratio": 0.0, "quantity": c["quantity"], "residual": c["residual"],
                              "tolerance": c["tolerance"], "location": c["location"]}
    return out


def _json_default(o):
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, np.bool_):
        return bool(o)
    raise TypeError(type(o))


def write_outputs(out_dir: Path, meta: dict, stations: list[dict]) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    all_checks = [c for s in stations for c in s["checks"]]
    summ = summarize(all_checks)
    report = {
        "meta": meta,
        "stations": [{**s["info"], "checks": s["checks"]} for s in stations],
        "summary": [{"station": k[0], "field": k[1], **v} for k, v in sorted(summ.items())],
    }
    (out_dir / "report.json").write_text(json.dumps(report, indent=1, default=_json_default))
    with open(out_dir / "stations.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["station", "field", "quantity", "location", "mesh", "yaml", "residual", "tolerance", "status", "note"])
        for c in all_checks:
            w.writerow([c["station"], c["field"], c["quantity"], json.dumps(c["location"], default=_json_default),
                        json.dumps(c["mesh"], default=_json_default), json.dumps(c["yaml"], default=_json_default),
                        c["residual"], c["tolerance"], c["status"], c["note"] or ""])


def print_summary(stations: list[dict]) -> None:
    all_checks = [c for s in stations for c in s["checks"]]
    summ = summarize(all_checks)
    print(f"{'r/R':5} {'F':3} {'pass':>5} {'fail':>5} {'bnd':>5} {'n/g':>5}  worst gated residual (res/tol) at")
    for (st, f), v in sorted(summ.items()):
        w = v["worst"]
        tail = ""
        if w:
            loc = w["location"].get("set") or w["location"].get("yaml_web") or ""
            tail = f"{w['quantity'][:44]} {w['residual']:+.3e} ({w['ratio']:.2f}x) {loc} a={w['location'].get('nd_arc', float('nan')):.3f}" \
                if "nd_arc" in w["location"] else f"{w['quantity'][:44]} {w['residual']:+.3e} ({w['ratio']:.2f}x) {loc}"
        print(f"{st:5} {f:3} {v['pass']:5d} {v['fail']:5d} {v['boundary']:5d} {v['not-gated']:5d}  {tail}")
    for s in stations:
        if "error" in s["info"]:
            print(f"station {s['info']['r_over_R']}: {s['info']['error']}")


def main(argv=None) -> int:
    args = parse_args(argv)
    yb = YamlBlade.from_file(args.yaml)
    mesh = None
    cache = Path(args.mesh_cache) if args.mesh_cache else None
    if cache and cache.exists():
        import pickle  # noqa: PLC0415

        mesh = pickle.loads(cache.read_bytes())
    if mesh is None:
        mesh = build_mesh(args.yaml, args.element_size)
        if cache:
            import pickle  # noqa: PLC0415

            cache.write_bytes(pickle.dumps(mesh))
    add_set_zrange(mesh)
    try:
        rev = subprocess.run(["git", "rev-parse", "--short", "HEAD"], capture_output=True, text=True).stdout.strip()
    except OSError:
        rev = ""
    meta = {"yaml": args.yaml, "element_size": args.element_size, "git_rev": rev,
            "stations": list(STATIONS), "convention_note": CONVENTION_NOTE,
            "blade_length_m": yb.length, "n_nodes": len(mesh["nodes"]), "n_elements": len(mesh["elements"])}
    stations = []
    for frac in STATIONS:
        try:
            stations.append(audit_station(yb, mesh, frac))
        except Exception as exc:  # the audit reports; it must not die on one station
            stations.append({"info": {"r_over_R": frac, "error": f"{type(exc).__name__}: {exc}"}, "checks": []})
    write_outputs(Path(args.output_dir), meta, stations)
    print_summary(stations)
    print(f"wrote {args.output_dir}/report.json and stations.csv")
    return 0


if __name__ == "__main__":
    sys.exit(main())
