"""Normalize supported windIO 2.0 blade fields for the existing NuMAD reader.

The source file is never rewritten. Explicit anchor arcs are authoritative; plane
intersection metadata is not used to recompute them. See windIO's turbine-model
construction guide, revision 63dcad638a9ac41d24f7a68172e24bde378c23e6.
"""

from copy import deepcopy

import numpy as np


def normalize_windio(data: dict) -> dict:
    """Return the legacy reader view of a v2.0 turbine, or v1 data unchanged.

    This is not a general schema converter. Weighted airfoil blends, nonzero
    normal section offsets and unresolved geometric anchors are not supported.
    They fail explicitly rather than silently changing the physical model.
    """
    version = str(data.get("windIO_version", "1.0"))
    if version.startswith("1."):
        return data
    if version != "2.0":
        raise ValueError(f"Unsupported windIO version: {version}")
    result = deepcopy(data)
    blade = result["components"]["blade"]
    outer = blade["outer_shape"]
    structure = blade["structure"]
    normal_offset = outer.get("section_offset_x", {}).get("values", [0.0])
    if np.any(np.asarray(normal_offset) != 0):
        raise ValueError("windIO section_offset_x is not supported by the blade reader")
    anchors = {a["name"]: a for a in structure.get("anchors", [])}
    for web in structure.get("webs", []):
        for anchor in web.get("anchors", []):
            if anchor["name"] in anchors:
                raise ValueError(f"Duplicate windIO anchor: {anchor['name']}")
            anchors[anchor["name"]] = anchor

    def arc(value, chain=()):
        if "grid" in value and "values" in value:
            return {"grid": deepcopy(value["grid"]), "values": deepcopy(value["values"])}
        reference = value.get("anchor", {})
        name, handle = reference.get("name"), reference.get("handle")
        if name not in anchors:
            raise ValueError(f"Unknown windIO anchor: {name}")
        key = (name, handle)
        if key in chain:
            raise ValueError(f"Cyclic windIO anchor: {name}.{handle}")
        if handle not in anchors[name]:
            raise ValueError(f"Missing explicit windIO anchor arc: {name}.{handle}")
        return arc(anchors[name][handle], (*chain, key))

    for item in [*structure.get("layers", []), *structure.get("webs", [])]:
        reference = item.get("start_nd_arc", {}).get("anchor", {})
        anchor = anchors.get(item["name"], anchors.get(reference.get("name"), {}))
        for field in ("start_nd_arc", "end_nd_arc"):
            if field in item:
                item[field] = arc(item[field])
        # Width and side remain construction inputs to NuMAD's named regions.
        if "width" in anchor:
            item["width"] = {
                k: deepcopy(anchor["width"][k]) for k in ("grid", "values")
            }
        side = anchor.get("plane_intersection", {}).get("side")
        if side in ("suction", "pressure"):
            item["side"] = side

    chord = outer["chord"]
    offset = outer["section_offset_y"]
    grid = np.union1d(chord["grid"], offset["grid"])
    chord_values = np.interp(grid, chord["grid"], chord["values"])
    if np.any(chord_values <= 0):
        raise ValueError("windIO chord must be positive for section_offset_y conversion")
    pitch_axis = np.interp(grid, offset["grid"], offset["values"]) / chord_values
    profiles = outer["airfoils"]
    if any(p.get("weight", [1.0]) != [1.0] for p in profiles):
        raise ValueError("Weighted windIO airfoil blends are not supported by the blade reader")
    blade["outer_shape_bem"] = {
        "reference_axis": deepcopy(blade["reference_axis"]),
        "chord": deepcopy(chord),
        "twist": {
            "grid": deepcopy(outer["twist"]["grid"]),
            "values": np.deg2rad(outer["twist"]["values"]).tolist(),
        },
        "pitch_axis": {"grid": grid.tolist(), "values": pitch_axis.tolist()},
        "airfoil_position": {
            "grid": [p["spanwise_position"] for p in profiles],
            "labels": [p["name"] for p in profiles],
        },
        "rthick": deepcopy(outer["rthick"]),
    }
    blade["internal_structure_2d_fem"] = structure
    for airfoil in result["airfoils"]:
        airfoil["relative_thickness"] = airfoil["rthick"]
        polars = []
        for polar in airfoil.get("polars", []):
            for table in polar.get("re_sets", []):
                entry = {"configuration": polar["configuration"], "re": table["re"]}
                for old, new in (("c_l", "cl"), ("c_d", "cd"), ("c_m", "cm")):
                    entry[old] = {
                        "grid": np.deg2rad(table[new]["grid"]).tolist(),
                        "values": deepcopy(table[new]["values"]),
                    }
                polars.append(entry)
        airfoil["polars"] = polars
    return result
