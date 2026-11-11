"""Cached windIO v2 IEA-15 blade for the v2 validation groups (odd/tasks/v2-validation-store.md, T1).

``build_v2_blade(element_size, collapse_ratio, skin_reference) -> (mesh, properties)`` is the v2 counterpart of the
v1 pair ``BladeMesh(yaml_file=..., element_size=..., skin_reference="mid").generate(renumber="rcm")`` +
``aeroelast.models.blade.model.build_rust_properties(generator.numad_mesh_data)``.

Mesh source: ``windio2_mesh.generate_mesh(topology, es, min_strip_m=1e-6, collapse_ratio=...)`` + ``map_mesh_materials`` on the windIO v2
YAML ``tests/IEA15MW/reference/upstream/IEA-15-240-RWT.yaml``. ``element_size`` must lie in [0.15, 0.25] m (AGENTS.md
rule 9; the v2 mesher rejects other values).

Reference-surface rule (OML -> mid-surface inset)
-------------------------------------------------
The v2 mesher keeps every skin node on the outer mould line (OML) and stacks the plies inward
(``MappedStack.reference == "outer_mold_line"``). With ``skin_reference="mid"`` (default) each skin node moves to the
laminate mid-surface; with ``"oml"`` nothing moves. For a node ``i`` shared by the skin elements ``k``::

    n_i = normalise( sum_k A_k * n_k )                 (n_k: stack.physical_normal, skin outward; A_k: element area)
    d_i = sum_k A_k * stack_k.midplane_offset_m / sum_k A_k      (signed; -h_k/2 for a stack from -h to 0)
    x_i <- x_i + d_i * n_i

i.e. the area-weighted mean of the real per-stack mid-plane offset (``MappedStack.midplane_offset_m``) along the
area-weighted skin normal. Web-only nodes stay; skin/web shared nodes move (as in the v1 ``apply_midsurface_inset``).
Unlike v1 the normal is NOT projected into the section plane and there is no ``*_TE_FLAT`` exception (the v2 mesh has
no such sets). Every laminate is then centred on its surface by ``build_rust_properties`` (z from -h/2 to +h/2), which
is the mid-surface reading. Mass is computed on the moved mesh by the consumers.

Properties
----------
One element set ``SEC{i}`` per distinct layup (order of first appearance, deterministic), and ``properties[SEC{i}]``
is produced by ``build_rust_properties`` itself on a numad-style dict (materials + sections), so material ordering
(windIO ``G=[G12,G13,G23]``), isotropic single-ply shortcuts and the 5/6 shear factor are exactly the v1 ones. Plies
are listed in the mapped bottom-to-top order; fibre angles are the windIO SOURCE angles (all 0 in the official blade),
referred to the solver's projected span direction.

MeshModel attributes provided (what the v1 consumers read)
----------------------------------------------------------
* node sets ``RootNodes`` (every node at the lowest eta, skin + webs), ``TipNodes`` (highest eta),
  ``allOuterShellNods`` (all skin nodes), ``allWebNods`` when webs exist;
* element sets ``SEC{i}``, ``allOuterShellEls`` (skin), ``allWebEls`` (webs), ``<web name>`` per web;
* ``mesh.coords_array``, ``mesh.node_id_to_index`` (renumbered with ``"rcm"`` like v1);
* ``mesh.numad_mesh_data``: the numad-style dict (``materials``, ``sections``) fed to ``build_rust_properties``,
  for the consumers that took ``generator.numad_mesh_data``;
* ``mesh.v2_info``: dict (cache key, element size, collapse ratio, skin reference, counts).

Not reproduced: the v1 station-named element sets ``{cc:02d}_{ss:02d}_{surface}_{component}`` (the v2 mesh has no
NuMAD station/segment structure), the ``*_SPAR*`` / ``*_SW`` / ``*_TE_FLAT`` named sets, and a per-ring node-count
guarantee. Consumers that need sections by span must bin ``coords_array[:, 2]`` or use the eta rings.

Cache
-----
Generation takes ~10 min at element size 0.25, so the RAW arrays (never live objects) are cached as an ``.npz`` under
``$SCRATCH/tmp/v2_blade_cache/`` (override: ``V2_BLADE_CACHE_DIR``), keyed by sha256(YAML) + a digest of
``src/aeroelast/models/blade/numad/mesh_gen/windio2_*.py`` + this module's ``CACHE_VERSION`` + element size, collapse
ratio, eta bounds and ``min_strip_m`` (``skin_reference`` is applied after the cache). A process-level memo sits on top. Produce the
cache with an sbatch job (never on a login node): ``python tests/support/v2_blade.py <es> [collapse_ratio]``.
"""

from __future__ import annotations

import hashlib
import json
import os
import sys
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_YAML = REPO_ROOT / "tests/IEA15MW/reference/upstream/IEA-15-240-RWT.yaml"
MESH_GEN_DIR = REPO_ROOT / "src/aeroelast/models/blade/numad/mesh_gen"
GMSH_BLADE = REPO_ROOT / "src/aeroelast/core/mesh/gmsh_blade.py"
CACHE_VERSION = 3  # 3: chord-scaled gmsh sizing (N, AR, min band in the key; gmsh_blade.py in the gmsh digest)
MESHERS = ("windio2", "gmsh")
ES_RANGE = (0.15, 0.25)
# Resolution guard of the public v2 adapter (core/mesh/windio2.py::geometry_mesh): stations closer than 1 micrometre of
# span are one. generate_mesh's own default (0) hits "Nonpositive shell Jacobian" at es 0.25 (z ~ 94 m, SLURM 11641645).
MIN_STRIP_M = 1e-6
SKIN_REFERENCES = ("mid", "oml")

_MEMO: dict = {}  # raw arrays by cache key
_BUILT: dict = {}  # (mesh, properties) by (cache key, skin_reference)


# --------------------------------------------------------------------------- cache key / location
def _sha256_file(path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def source_digest(mesher="windio2") -> str:
    """Digest of the mesher sources, names and bytes, sorted: ``windio2_*.py`` + ``collapse_opt.py`` for the windio2 mesher,
    plus ``gmsh_skin.py`` and ``core/mesh/gmsh_blade.py`` for the gmsh one (which also runs ``windio2_*`` code)."""
    h = hashlib.sha256()
    files = [*MESH_GEN_DIR.glob("windio2_*.py"), MESH_GEN_DIR / "collapse_opt.py"]
    if mesher == "gmsh":
        files.extend([MESH_GEN_DIR / "gmsh_skin.py", GMSH_BLADE])
    for p in sorted(files):
        h.update(p.name.encode())
        h.update(p.read_bytes())
    return h.hexdigest()


def cache_key(element_size, collapse_ratio=None, eta_bounds=(0.0, 1.0), yaml_file=DEFAULT_YAML,
              min_strip_m=MIN_STRIP_M, mesher="windio2", chord_divisions=None, span_aspect=1.0, min_band_elements=2,
              min_band_fraction=1.0) -> str:
    """Stable key of the raw mesh: YAML bytes + mesher sources + generator arguments.  With ``chord_divisions`` (gmsh only) the
    key carries N, the spanwise aspect, the minimum band elements and ``min_band_fraction`` (None, the no-relative-floor
    behaviour, is encoded as null so it cannot collide with a number) and the element size drops out (it is not used)."""
    args = {
        "v": CACHE_VERSION,
        "es": float(element_size),
        "collapse": None if not collapse_ratio else float(collapse_ratio),
        "eta": [float(eta_bounds[0]), float(eta_bounds[1])],
        "min_strip": float(min_strip_m),
    }
    if chord_divisions is not None:
        if mesher != "gmsh":
            raise ValueError("chord_divisions is a gmsh mesher option")
        del args["es"]
        args["chord"] = {"N": float(chord_divisions), "AR": float(span_aspect), "min_band": int(min_band_elements),
                         "min_band_fraction": None if min_band_fraction is None else float(min_band_fraction)}
    if mesher != "windio2":
        args["mesher"] = mesher  # the default key stays what it was for the windio2 mesher (no extra field)
    blob = _sha256_file(yaml_file) + (source_digest() if mesher == "windio2" else source_digest(mesher)) + json.dumps(args, sort_keys=True)
    return hashlib.sha256(blob.encode()).hexdigest()


def cache_dir() -> Path:
    explicit = os.environ.get("V2_BLADE_CACHE_DIR")
    if explicit:
        return Path(explicit)
    scratch = os.environ.get("SCRATCH")
    if not scratch:
        raise RuntimeError("set SCRATCH (or V2_BLADE_CACHE_DIR): the v2 blade cache lives outside the repository")
    return Path(scratch) / "tmp" / "v2_blade_cache"


def cache_path(key: str) -> Path:
    return cache_dir() / f"v2blade-{key[:24]}.npz"


# --------------------------------------------------------------------------- raw generation (slow) and (de)serialisation
def _element_areas(x, elements):
    areas = np.empty(len(elements))
    for k, el in enumerate(elements):
        p = x[el[el >= 0]]
        cross = np.cross(p[2] - p[0], p[-1] - p[1]) / 2 if len(p) == 4 else np.cross(p[1] - p[0], p[2] - p[0]) / 2
        areas[k] = np.linalg.norm(cross)
    return areas


def _generate_raw(element_size, collapse_ratio, eta_bounds, yaml_file, min_strip_m=MIN_STRIP_M, mesher="windio2", es_range=None,
                  chord=None) -> dict:
    """Run the mesher (slow) and return plain arrays + JSON-able tables."""
    from aeroelast.models.blade.numad.mesh_gen.windio2_mesh import generate_mesh, map_mesh_materials
    from aeroelast.models.blade.numad.objects.blade import Blade

    blade = Blade()
    blade.read_yaml(str(yaml_file))
    topology = blade.windio2_topology
    if mesher == "gmsh":
        from aeroelast.models.blade.numad.mesh_gen.gmsh_skin import generate_mesh_gmsh

        extra = {} if es_range is None else {"h_range": tuple(es_range)}  # the mesher has its own range check
        extra.update(chord or {})
        raw = generate_mesh_gmsh(topology, element_size, with_materials=True, **extra)  # quality gate REJECTS, never repairs
    else:
        raw = map_mesh_materials(
            topology, generate_mesh(topology, element_size, tuple(eta_bounds), min_strip_m=min_strip_m,
                          collapse_ratio=collapse_ratio or None)
        )
    ne = len(raw.elements)
    params = np.full((ne, 4, 2), np.nan)
    for k in range(ne):
        p = np.asarray(raw.element_parameters[k], float)
        params[k, : len(p)] = p
    surfaces = [str(s) for s in raw.surfaces]

    layups, layup_index = [], {}
    layup_id = np.empty(ne, dtype=np.int64)
    normals = np.empty((ne, 3))
    offsets = np.empty(ne)
    materials = {}
    for k, stack in enumerate(raw.material_stacks):
        layers = []
        for layer in stack.layers:
            s = layer.source
            layers.append([s.material_name, round(float(s.thickness_m), 14), float(s.fiber_angle_deg)])
            if s.material_name not in materials:
                m = s.material
                tup = lambda v: [float(t) for t in np.atleast_1d(np.asarray(v, float))]  # noqa: E731
                materials[s.material_name] = {
                    "E": tup(m["E"]) if np.ndim(m["E"]) else float(m["E"]),
                    "nu": tup(m["nu"]) if np.ndim(m["nu"]) else float(m["nu"]),
                    "G": (tup(m["G"]) if np.ndim(m["G"]) else float(m["G"])) if "G" in m else None,
                    "rho": float(np.mean(np.asarray(m["rho"], float))),
                }
        key = json.dumps(layers)
        if key not in layup_index:
            layup_index[key] = len(layups)
            layups.append(layers)
        layup_id[k] = layup_index[key]
        normals[k] = np.asarray(stack.physical_normal, float)
        offsets[k] = float(stack.midplane_offset_m)
    nodes = np.asarray(raw.nodes, float)
    tables = {"layups": layups, "materials": materials, "surfaces": sorted(set(surfaces), key=surfaces.index)}
    return {
        "nodes": nodes,
        "elements": np.asarray(raw.elements, dtype=np.int64),
        "surfaces": np.array(surfaces),
        "params": params,
        "layup_id": layup_id,
        "normals": normals,
        "offsets": offsets,
        "areas": _element_areas(nodes, np.asarray(raw.elements, dtype=np.int64)),
        "tables_json": np.frombuffer(json.dumps(tables).encode(), dtype=np.uint8),
    }


def _write_cache(path: Path, arrays: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f".{path.name}.{os.getpid()}.tmp.npz")
    np.savez(tmp, **arrays)
    os.replace(tmp, path)


def _resolve_collapse(collapse_ratio):
    """An explicit ratio wins; otherwise the environment variable ``V2_BLADE_COLLAPSE_RATIO`` (e.g. ``0.15``), so a
    whole test campaign runs on the thin-strip-collapsed mesh without editing tests; otherwise ``None`` (plain)."""
    if collapse_ratio is None and os.environ.get("V2_BLADE_COLLAPSE_RATIO"):
        return float(os.environ["V2_BLADE_COLLAPSE_RATIO"])
    return collapse_ratio


def _resolve_element_size(element_size):
    """``V2_BLADE_ELEMENT_SIZE`` (campaign convergence runs) replaces the size a test hard-codes; unset = unchanged.
    The range check of ``load_raw`` still applies, so an override cannot reach a size the rule forbids."""
    override = os.environ.get("V2_BLADE_ELEMENT_SIZE")
    return float(override) if override else element_size


def _resolve_mesher(mesher):
    """An explicit mesher wins; otherwise ``V2_BLADE_MESHER`` (``gmsh``), so a campaign runs on the gmsh mesh without
    editing tests; otherwise ``windio2``.  The gmsh mesher builds the full blade only and takes no collapse ratio."""
    mesher = mesher or os.environ.get("V2_BLADE_MESHER") or "windio2"
    if mesher not in MESHERS:
        raise ValueError(f"mesher {mesher!r} not in {MESHERS}")
    return mesher


def load_raw(element_size=0.25, collapse_ratio=None, eta_bounds=(0.0, 1.0), *, yaml_file=DEFAULT_YAML,
             generate=True, min_strip_m=MIN_STRIP_M, mesher=None, es_range=None, chord_divisions=None, span_aspect=1.0,
             min_band_elements=2, min_band_fraction=1.0) -> dict:
    """Raw mesh arrays: process memo -> disk cache -> (if ``generate``) the slow mesher. ``generate=False`` raises
    ``FileNotFoundError`` on a cold cache instead of spending ~10 min. ``collapse_ratio=None`` follows
    ``V2_BLADE_COLLAPSE_RATIO`` (see ``_resolve_collapse``), so every consumer of the raw arrays agrees with
    ``build_v2_blade`` on which mesh is meant."""
    mesher = _resolve_mesher(mesher)
    element_size = _resolve_element_size(element_size)
    collapse_ratio = _resolve_collapse(collapse_ratio)
    if mesher == "gmsh":
        if collapse_ratio or tuple(eta_bounds) != (0.0, 1.0):
            raise ValueError("the gmsh mesher builds the full blade and takes no collapse ratio")
    # ``es_range`` lets a blade whose converged size is outside the IEA-15 window (IEA-3.4: 0.1 m) say so explicitly
    rng = ES_RANGE if es_range is None else tuple(es_range)
    chord = None
    if chord_divisions is not None:  # chord-scaled sizing replaces the element size (and its rule-9 range)
        if mesher != "gmsh":
            raise ValueError("chord_divisions is a gmsh mesher option")
        chord = {"chord_divisions": chord_divisions, "span_aspect": span_aspect, "min_band_elements": min_band_elements,
                 "min_band_fraction": min_band_fraction}
    elif not rng[0] <= float(element_size) <= rng[1]:
        raise ValueError(f"element_size {element_size} outside {rng} (v2 mesher range; AGENTS.md rule 9)")
    key = cache_key(element_size, collapse_ratio, eta_bounds, yaml_file, min_strip_m, mesher, **(chord or {}))
    if key in _MEMO:
        return _MEMO[key]
    path = cache_path(key)
    if path.exists():
        with np.load(path, allow_pickle=False) as z:
            arrays = {name: z[name] for name in z.files}
    elif generate:
        arrays = _generate_raw(element_size, collapse_ratio, eta_bounds, yaml_file, min_strip_m, mesher, es_range=es_range,
                               chord=chord)
        _write_cache(path, arrays)
    else:
        raise FileNotFoundError(f"cold v2 blade cache {path}; produce it with an sbatch job")
    arrays["key"] = key
    _MEMO[key] = arrays
    return arrays


# --------------------------------------------------------------------------- MeshModel / properties from raw arrays (fast)
def _numad_data(tables: dict, layup_names: list) -> dict:
    materials = []
    for name, m in tables["materials"].items():
        elastic = {"E": m["E"], "nu": m["nu"]}
        if m["G"] is not None:
            elastic["G"] = m["G"]
        materials.append({"name": name, "density": m["rho"], "elastic": elastic})
    sections = [{"elementSet": name, "layup": layup} for name, layup in zip(layup_names, tables["layups"])]
    return {"materials": materials, "sections": sections}


def _inset(raw: dict) -> np.ndarray:
    """OML -> mid-surface node coordinates (rule in the module docstring)."""
    x = np.array(raw["nodes"], float)
    nsum, osum, wsum = np.zeros_like(x), np.zeros(len(x)), np.zeros(len(x))
    for k, el in enumerate(raw["elements"]):
        if raw["surfaces"][k] != "skin":
            continue
        ids = el[el >= 0]
        a = raw["areas"][k]
        nsum[ids] += a * raw["normals"][k]
        osum[ids] += a * raw["offsets"][k]
        wsum[ids] += a
    moved = wsum > 0
    norm = np.maximum(np.linalg.norm(nsum[moved], axis=1), 1e-300)[:, None]
    x[moved] += (osum[moved] / wsum[moved])[:, None] * nsum[moved] / norm
    return x


def mesh_from_raw(raw: dict, skin_reference: str = "mid"):
    """(MeshModel, properties) from the cached arrays; seconds."""
    from aeroelast.core.mesh.entities import ElementSet, ElementType, MeshElement, Node, NodeSet
    from aeroelast.core.mesh.model import MeshModel
    from aeroelast.models.blade.model import build_rust_properties

    if skin_reference not in SKIN_REFERENCES:
        raise ValueError(f"skin_reference must be one of {SKIN_REFERENCES}, got {skin_reference!r}")
    tables = json.loads(bytes(raw["tables_json"]).decode())
    x = _inset(raw) if skin_reference == "mid" else np.array(raw["nodes"], float)
    mesh = MeshModel()
    for p in x:
        mesh.add_node(Node(p))
    nodes = list(mesh.nodes)
    surfaces = [str(s) for s in raw["surfaces"]]
    elements = []
    for el in raw["elements"]:
        ids = [int(i) for i in el if i >= 0]
        e = MeshElement([nodes[i] for i in ids], ElementType.quad if len(ids) == 4 else ElementType.triangle)
        mesh.add_element(e)
        elements.append(e)

    layup_names = [f"SEC{i}" for i in range(len(tables["layups"]))]
    by_layup, by_surface = {}, {}
    for k, e in enumerate(elements):
        by_layup.setdefault(int(raw["layup_id"][k]), set()).add(e)
        by_surface.setdefault(surfaces[k], set()).add(e)
    for i, es in sorted(by_layup.items()):
        mesh.add_element_set(ElementSet(name=layup_names[i], elements=es))
    skin_els = by_surface.get("skin", set())
    web_els = set().union(*[v for s, v in by_surface.items() if s != "skin"]) if len(by_surface) > 1 else set()
    mesh.add_element_set(ElementSet(name="allOuterShellEls", elements=skin_els))
    if web_els:
        mesh.add_element_set(ElementSet(name="allWebEls", elements=web_els))
        for s, v in by_surface.items():
            if s != "skin":
                mesh.add_element_set(ElementSet(name=s, elements=v))

    eta = np.full(len(x), np.nan)
    for k, el in enumerate(raw["elements"]):
        for j, i in enumerate(el[el >= 0]):
            eta[i] = raw["params"][k, j, 0]
    lo, hi = np.nanmin(eta), np.nanmax(eta)
    mesh.add_node_set(NodeSet(name="RootNodes", nodes={nodes[i] for i in np.nonzero(eta == lo)[0]}))
    mesh.add_node_set(NodeSet(name="TipNodes", nodes={nodes[i] for i in np.nonzero(eta == hi)[0]}))
    mesh.add_node_set(NodeSet(name="allOuterShellNods", nodes={n for e in skin_els for n in e.nodes}))
    if web_els:
        mesh.add_node_set(NodeSet(name="allWebNods", nodes={n for e in web_els for n in e.nodes}))

    numad = _numad_data(tables, layup_names)
    properties = build_rust_properties(numad)
    mesh.renumber_mesh(algorithm="rcm")
    mesh.numad_mesh_data = numad
    mesh.v2_info = {
        "cache_key": str(raw["key"]), "n_nodes": len(x), "n_elements": len(elements),
        "n_sections": len(layup_names), "skin_reference": skin_reference,
    }
    return mesh, properties


def build_v2_blade(element_size=0.25, collapse_ratio=None, skin_reference="mid", *, eta_bounds=(0.0, 1.0),
                   yaml_file=DEFAULT_YAML, min_strip_m=MIN_STRIP_M, mesher=None, es_range=None, chord_divisions=None,
                   span_aspect=1.0, min_band_elements=2, min_band_fraction=1.0):
    """windIO v2 IEA-15 blade ``(mesh, properties)``; cached on disk and per process (see module docstring).

    ``collapse_ratio=None`` follows ``V2_BLADE_COLLAPSE_RATIO`` (see ``_resolve_collapse``); an explicit argument
    always wins, and with neither set the mesh is the plain one.  ``chord_divisions`` N (gmsh mesher) sizes the mesh by the local
    chord instead of ``element_size``: chordwise c/N, spanwise ``span_aspect`` * c/N, >= ``min_band_elements`` across a band of width >= ``min_band_fraction`` (default 1.0) * c/N; None = every band.
    """
    raw = load_raw(element_size, collapse_ratio, eta_bounds, yaml_file=yaml_file, min_strip_m=min_strip_m, mesher=mesher,
                   es_range=es_range, chord_divisions=chord_divisions, span_aspect=span_aspect,
                   min_band_elements=min_band_elements, min_band_fraction=min_band_fraction)
    memo_key = (raw["key"], skin_reference)
    if memo_key not in _BUILT:
        _BUILT[memo_key] = mesh_from_raw(raw, skin_reference)
    return _BUILT[memo_key]


if __name__ == "__main__":  # sbatch entry: build (and cache) one mesh, print a summary
    import time

    es = float(sys.argv[1]) if len(sys.argv) > 1 else 0.25
    cr = float(sys.argv[2]) if len(sys.argv) > 2 and sys.argv[2] != "none" else None
    t0 = time.time()
    mesh, props = build_v2_blade(es, cr)
    print(f"es={es} collapse={cr} nodes={mesh.node_count} elements={mesh.elements_count} sections={len(props)} "
          f"cache={cache_path(mesh.v2_info['cache_key'])} t={time.time() - t0:.0f}s")
