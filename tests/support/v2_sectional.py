"""Sectional extractor for the windIO v2 blade mesh (odd/tasks/v2-validation-store.md, T4; user authorisation L5).

The v1 ``SectionalExtractor`` (``aeroelast.postprocess.sectional``) reads stations from the NuMAD element sets
``{cc}_{ss}_{surface}_{component}``, which the v2 mesh does not have (``tests/support/v2_blade.py``).  The v2 mesh is a
stack of node rings, one per eta row, and every element spans two consecutive rings.  ``V2SectionalExtractor`` keeps
EVERY definition of the v1 class (lumped-mass partition, tangent-frame curvature EI, ``ei_from_static_response``,
``select_stations``, the slice stiffness) and replaces one method, ``station_membership``:

    station r  =  the nodes of the elements whose lower ring is r  (ring r and ring r+1),   key r = 0 .. n_rings-2

An element that does not span exactly two consecutive rings (collapsed / collapsed+optimised meshes, where merged nodes
and optimiser-slid apex nodes make it touch one or three rings) is assigned to the one bay containing its centroid z
(ring z = node mean of the ring).  Plain-mesh elements never take that branch, so the plain membership is unchanged.

which is the v1 station semantics (one bay between two section rings; neighbouring stations share a ring, the first
station's lower face is ``RootNodes``).  ``z`` of a station is its node mean, as in v1.

Ring identification
-------------------
The mesh handed to the consumers is renumbered (rcm) and the skin moved to the mid-surface, so the eta row of a node is
recovered from the raw builder arrays by EXACT coordinate match (``ring_ids_for_v2_mesh``: the mid-surface positions
are recomputed with the builder's own ``_inset``, bit-identical to the mesh).  ``ring_ids_from_z`` is the geometry-only
fallback (split sorted z at gaps wider than a tolerance) used by the synthetic unit tests.  Rings of the official blade
are slightly non-planar in z (spread <= 0.063 m against >= 0.043 m to the next ring), so the exact match is the one used
for the validation test; consecutive rings never overlap in z.

Differences from v1, stated
---------------------------
* ``extract`` computes the MASS path only by default (``stiffness=False``): the v1 test reads ``sec.z`` and
  ``sec.mass_dens`` only and the slice stiffness of ``extract`` is flagged experimental in the v1 docstring.  The other
  ``SectionProperties`` fields are NaN.  ``stiffness=True`` runs the inherited slice homogenisation.
* ``ea_gj_from_static_response_v2`` is the v1 ``tools/run_s1_sectional.ea_gj_from_static_response`` run unchanged with
  this extractor substituted for the station membership (module attribute swap for the duration of the call).
* ``bay_centres_v2`` gives the reference stations (v1: ``native_bay_centres``): centres of the 49 bays of the
  windIO v2 reference-axis z grid (2.3878 m, identical to the grid of the v1 yaml).
"""

from __future__ import annotations

import contextlib
from typing import Dict, Optional

import numpy as np

from tests.support.paths import REPO_ROOT

V2_YAML = REPO_ROOT / "tests/IEA15MW/reference/upstream/IEA-15-240-RWT.yaml"


# --------------------------------------------------------------------------- ring identification
def ring_ids_from_z(coords, tol: float) -> np.ndarray:
    """Ring index per node by splitting the sorted z at gaps wider than ``tol`` (geometry only)."""
    z = np.asarray(coords, float)[:, 2]
    order = np.argsort(z, kind="stable")
    cuts = np.nonzero(np.diff(z[order]) > tol)[0] + 1
    ring = np.empty(len(z), dtype=np.int64)
    for r, chunk in enumerate(np.split(order, cuts)):
        ring[chunk] = r
    return ring


def ring_ids_for_v2_mesh(mesh) -> np.ndarray:
    """Eta-row index per mesh node (0 = root ring) for a mesh built by ``tests.support.v2_blade.build_v2_blade``.

    Exact coordinate match against the raw builder arrays (the cache memo of the same process or the disk cache).
    Raises ``ValueError`` if a node has no unique raw counterpart.
    """
    from tests.support import v2_blade as vb

    info = mesh.v2_info
    raw = vb._MEMO.get(info["cache_key"])
    if raw is None:  # same cache key as the build: load without regenerating
        raw = vb.load_raw(info.get("element_size", 0.25), generate=False)
        if raw["key"] != info["cache_key"]:
            raise ValueError("raw v2 blade arrays for this mesh are not available (cache key mismatch)")
    x = vb._inset(raw) if info["skin_reference"] == "mid" else np.asarray(raw["nodes"], float)
    el, params = raw["elements"], raw["params"]
    eta = np.full(len(x), np.nan)
    for k, e in enumerate(el):
        for j, i in enumerate(e[e >= 0]):
            eta[i] = params[k, j, 0]
    _, raw_ring = np.unique(eta, return_inverse=True)
    lookup = {}
    for i, p in enumerate(x):
        key = p.tobytes()
        if lookup.setdefault(key, raw_ring[i]) != raw_ring[i]:  # coincident nodes on one ring are fine
            raise ValueError("coincident raw nodes on different rings: exact match is ambiguous")
    coords = np.ascontiguousarray(mesh.coords_array, dtype=float)
    try:
        return np.array([lookup[c.tobytes()] for c in coords], dtype=np.int64)
    except KeyError as exc:  # pragma: no cover - means the mesh is not the builder's
        raise ValueError("mesh node without a raw v2 counterpart") from exc


# --------------------------------------------------------------------------- the extractor
def _base():
    from aeroelast.postprocess.sectional import SectionalExtractor

    return SectionalExtractor


def make_v2_extractor_class():
    """Build the class lazily (importing ``aeroelast`` needs the bootstrapped environment)."""
    Base = _base()

    class V2SectionalExtractor(Base):
        """``SectionalExtractor`` whose stations are the ring-to-ring bays of a v2 mesh (see module docstring)."""

        def __init__(self, mesh, properties=None, *, ring_ids=None, ring_tol: Optional[float] = None, **kw):
            super().__init__(mesh, properties, **kw)
            if ring_ids is None:
                ring_ids = (ring_ids_for_v2_mesh(mesh) if ring_tol is None
                            else ring_ids_from_z(self.coords, ring_tol))
            self.ring_ids = np.asarray(ring_ids, dtype=np.int64)
            if self.ring_ids.shape != (mesh.node_count,):
                raise ValueError("ring_ids must hold one entry per mesh node")

        def station_membership(self) -> Dict[int, Dict[str, np.ndarray]]:
            idx = self.mesh.node_id_to_index
            ring = self.ring_ids
            n_rings = int(ring.max()) + 1
            zr = None  # ring z positions, built only if some element is not a plain two-ring element
            bays: Dict[int, set] = {}
            for e in self.mesh.elements:
                nodes = [idx[n] for n in e.node_ids]
                rs = ring[nodes]
                lo = int(rs.min())
                if int(rs.max()) - lo != 1:
                    # collapsed / optimised meshes: an element can touch one or three consecutive rings.  It goes to
                    # the single bay whose z interval [zr[r], zr[r+1]) holds its centroid z (partition, deterministic).
                    if zr is None:
                        zr = np.array([self.coords[ring == r, 2].mean() for r in range(n_rings)])
                    cz = float(self.coords[nodes, 2].mean())
                    lo = int(np.clip(np.searchsorted(zr, cz, side="right") - 1, 0, n_rings - 2))
                bays.setdefault(lo, set()).update(nodes)
            out: Dict[int, Dict[str, np.ndarray]] = {}
            for r in sorted(bays):
                arr = np.array(sorted(bays[r]), dtype=np.int64)
                out[r] = {"nodes": arr, "z": float(self.coords[arr, 2].mean())}
            return out

        def extract(self, n_slices=None, stations=None, stiffness: bool = False):
            """Mass per station (v1 definition); ``stiffness=True`` adds the inherited slice homogenisation."""
            from aeroelast.postprocess.sectional import SectionProperties, select_stations

            if stiffness:
                return super().extract(n_slices=n_slices, stations=stations)
            membership = self.station_membership()
            order = sorted(membership)
            keep = select_stations(np.array([membership[st]["z"] for st in order]), stations)
            zs, m, cg = [], [], []
            for st in (order[i] for i in keep):
                ms = self.section_mass(st, membership)
                zs.append(membership[st]["z"])
                m.append(ms["mass_dens"])
                cg.append(ms["cg"])
            nan = np.full(len(zs), np.nan)
            return SectionProperties(z=np.asarray(zs), ea=nan, gj=nan.copy(), ei_flap=nan.copy(), ei_edge=nan.copy(),
                                     mass_dens=np.asarray(m), cg=np.asarray(cg))

    return V2SectionalExtractor


def V2SectionalExtractor(mesh, properties=None, **kw):  # noqa: N802 - class-like factory
    """Construct the v2 extractor (see :func:`make_v2_extractor_class`)."""
    return make_v2_extractor_class()(mesh, properties, **kw)


# --------------------------------------------------------------------------- reference stations and EA/GJ
def bay_centres_v2(yaml_file=V2_YAML) -> np.ndarray:
    """Centres of the bays of the windIO v2 reference-axis z grid [m] (counterpart of ``native_bay_centres``)."""
    import yaml

    with open(yaml_file) as f:
        d = yaml.safe_load(f)
    z = np.asarray(d["components"]["blade"]["reference_axis"]["z"]["values"], float)
    return 0.5 * (z[:-1] + z[1:])


@contextlib.contextmanager
def v2_extractor_substituted(ring_ids=None):
    """Within the block, ``aeroelast.postprocess.sectional.SectionalExtractor(mesh, props)`` builds the v2 extractor.

    ``tools/run_s1_sectional.ea_gj_from_static_response`` imports the class at call time, so swapping the module
    attribute runs that v1 function unchanged on the v2 station membership.
    """
    import aeroelast.postprocess.sectional as sectional

    cls = make_v2_extractor_class()
    original = sectional.SectionalExtractor

    def factory(mesh, props=None, *a, **kw):
        return cls(mesh, props, *a, ring_ids=ring_ids, **kw)

    sectional.SectionalExtractor = factory
    try:
        yield
    finally:
        sectional.SectionalExtractor = original


def ea_gj_from_static_response_v2(mesh, props, *, numad_data=None, stations=None, ring_ids=None, **kw):
    """v1 ``ea_gj_from_static_response`` (GJ static torsion + EA CLT) on the v2 station membership."""
    import sys

    sys.path.insert(0, str(REPO_ROOT / "tools"))
    from run_s1_sectional import ea_gj_from_static_response

    with v2_extractor_substituted(ring_ids):
        return ea_gj_from_static_response(mesh, props, numad_data=numad_data, stations=stations, **kw)
