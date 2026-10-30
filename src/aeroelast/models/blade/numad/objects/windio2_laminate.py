"""Source-order continuum laminates for the isolated windIO 2 topology API.

No physical ply splitting, stack-normal direction or element axes are inferred.
Native one-dimensional distributions use linear interpolation without extension.
"""
from dataclasses import dataclass
from collections.abc import Mapping

import numpy as np
from fractions import Fraction
from bisect import bisect_right
from .windio2_precision import extended, rational


@dataclass(frozen=True)
class LayerSample:
    name: str
    material_name: str
    material: Mapping  # immutable original record: E/G [Pa], rho [kg/m^3]
    thickness_m: float
    fiber_angle_deg: float  # source section-local convention, NOT a global angle
    n_plies: float | None  # only when explicitly supplied; never inferred/rounded


@dataclass(frozen=True)
class LaminateCell:
    surface: str  # skin or source web name
    nd_bounds: tuple[float, float]
    xyz: tuple[tuple[float, float, float], tuple[float, float, float]]  # metres
    layers: tuple[LayerSample, ...]  # original source list order

    @property
    def thickness_m(self):
        return sum(layer.thickness_m for layer in self.layers)


class NativeField:
    def __init__(self, data, label, nonnegative=False):
        grid, values = np.asarray(data["grid"]), np.asarray(data["values"])
        if (grid.ndim != 1 or values.shape != grid.shape or len(grid) < 2
                or not np.all(np.isfinite(grid)) or not np.all(np.isfinite(values))
                or np.any(np.diff(grid) <= 0) or grid[0] < 0 or grid[-1] > 1
                or (nonnegative and np.any(values < 0))):
            raise ValueError(f"Unsupported windIO laminate field: {label}; expected finite native 1D data")
        self.grid, self.values = tuple(grid), tuple(values)

    def at(self, eta):
        if self.grid[0] <= eta <= self.grid[-1]:
            if isinstance(eta, Fraction):
                j = min(bisect_right(self.grid, eta)-1, len(self.grid)-2)
                a, b = map(rational, self.grid[j:j+2])
                x, y = map(rational, self.values[j:j+2])
                return extended(x+(y-x)*(eta-a)/(b-a))
            return float(np.interp(eta, self.grid, self.values))
        return None


def span_domain(record):
    lo, hi = record.get("start_nd_grid", 0.), record.get("end_nd_grid", 1.)
    if not np.isfinite([lo, hi]).all() or not 0 <= lo <= hi <= 1:
        raise ValueError("Invalid windIO laminate span domain")
    return lo, hi


class WindIO2Laminate:
    def __init__(self, source):
        structure = source["components"]["blade"]["structure"]
        self.materials = {m["name"]: m for m in source["materials"]}
        self.layers = []
        events = set()
        for record in (*structure.get("webs", ()), *structure["layers"]):
            events.update(span_domain(record))
        webs = {w["name"] for w in structure.get("webs", ())}
        for layer in structure["layers"]:
            if layer["material"] not in self.materials:
                raise ValueError(f"Unknown windIO laminate material: {layer['material']}")
            if "web" in layer and layer["web"] not in webs:
                raise ValueError(f"Unknown windIO laminate web: {layer['web']}")
            if "thickness" not in layer:
                raise ValueError(f"windIO laminate requires dimensional thickness: {layer['name']}")
            fields = [NativeField(layer["thickness"], layer["name"] + ".thickness", True),
                      NativeField(layer.get("fiber_orientation", {"grid": [0., 1.],
                                                                  "values": [0., 0.]}),
                                  layer["name"] + ".fiber_orientation")]
            if "n_plies" in layer:
                fields.append(NativeField(layer["n_plies"], layer["name"] + ".n_plies", True))
            events.update(x for field in fields for x in field.grid)
            self.layers.append((layer, tuple(fields)))
        self.event_stations = tuple(sorted(events))

    def cells(self, eta, surface, positions, bounds, xyz):
        """Assign complete open cells; boundaries already include every layer edge."""
        samples = []
        for layer, fields in self.layers:
            if layer.get("web", "skin") != surface:
                continue
            lo, hi = span_domain(layer)
            if not lo <= eta <= hi:
                continue
            edges = [positions.get(("layer", layer["name"], h))
                     for h in ("start_nd_arc", "end_nd_arc")]
            if None in edges or edges[0] == edges[1]:
                continue
            a, b = edges
            if surface != "skin" and a > b:
                raise ValueError(f"Reversed windIO web-local layer bounds: {layer['name']}")
            values = [field.at(eta) for field in fields]
            if None in values or values[0] == 0:
                continue
            t, angle = values[:2]
            sample = LayerSample(layer["name"], layer["material"], self.materials[layer["material"]],
                                 t, angle, values[2] if len(values) == 3 else None)
            samples.append((a, b, sample))
        return tuple(LaminateCell(surface, (lo, hi), ends,
                                 tuple(sample for a, b, sample in samples
                                       if ((a <= lo and hi <= b) if a < b
                                           else (lo >= a or hi <= b))))
                     for (lo, hi), ends in zip(bounds, xyz, strict=True))
