"""Static pure-torsion stiffness of the IEA-15MW shell, per station (route (c), verdict section 22.10).

Three independent routes to the blade's torsional stiffness:

* **(a) the anchor deck.** ``GKt = 8.748569e10 N.m^2`` at the root, from the official
  IEA-15-240-RWT BeamDyn blade deck via ``K66toPropsDecoupled`` (section 22.8).  Read from the
  deck, not re-derived.
* **(b) the modal ratio.** ``GJ_shell/GJ_ref = (4.000/4.290)^2 = 0.870`` against NuMAD's 1st
  torsion, whose **reference scatter is 15.1%** (section 15.1).  Quoted, never re-run.
* **(c) the new static measurement.** Apply a self-equilibrated pure torque to the real shell
  blade - ``+T`` as a wall shear flow ``q = T/(2A)`` on the tip ring and ``-T`` on the root ring
  (the realisation validated on the exact Bredt tube, ``tests/test_thin_walled_tube_torsion.py``
  / ``d4fec33``; **not** the production ``_distribute``, which section 22.9 proved is
  circle-tangential and not a shear flow), remove the six rigid modes, solve once, and read the
  local twist rate ``dphi/dz`` between consecutive physical stations from the ring section
  rotation - the antisymmetric part of the ring's affine fit (section 22.8), never ``mean(theta_z)``.

The applied system is checked with the independent torque ruler ``sum(x Fy - y Fx)`` (section
22.9's pattern) and with ``sum F = 0``; the enclosed area ``A`` is a hand-written shoelace of the
ordered contour, independent of the assembled force vector.  The assertions are only the
physics-pinned invariants: the ruler identity, self-equilibration, the removed rigid modes, and the
sign of the accumulating twist.  No bound is asserted against the deck, the modal ratio or Zhou.
"""

from __future__ import annotations

import sys
from collections import Counter
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("_aeroelast", reason="Rust backend not available")

from _aeroelast import PyMeshAssembler  # noqa: E402
from scipy.sparse import coo_matrix  # noqa: E402

from aeroelast.core.mesh.entities import MeshElement, Node  # noqa: E402
from aeroelast.models.blade.model import Blade  # noqa: E402

import tests.test_blade_iea15mw_validation as blade_validation  # noqa: E402
import tests.test_thin_walled_tube_torsion as tube  # noqa: E402
from tests.test_blade_rated_twist import STATION_GAP_TOLERANCE, _physical_stations  # noqa: E402

_OPENFAST_TOOLBOX = Path.home() / "openfast_toolbox"
if _OPENFAST_TOOLBOX.is_dir() and str(_OPENFAST_TOOLBOX) not in sys.path:
    sys.path.insert(0, str(_OPENFAST_TOOLBOX))
from openfast_toolbox.converters.beam import K66toPropsDecoupled  # noqa: E402

YAML = Path(__file__).resolve().parent / "IEA-15-240-RWT.yaml"
DECK = (Path(__file__).resolve().parent.parent / ".sources" / "openfast" / "iea15mw"
        / "IEA-15-240-RWT" / "IEA-15-240-RWT_BeamDyn_blade.dat")
TORQUE = 1.0e7  # [N.m] the applied pure couple
RULER_TOL = 1e-9  # exact by construction: the flow is assembled from the contour, the ruler from it
GJ_REF_ROOT = 8.748569e10  # [N.m^2] the deck GKt at the root (section 22.8), read not re-derived
MODAL_RATIO = (4.000 / 4.290) ** 2  # (section 15.1) recorded, the modal analysis is NOT re-run
MODAL_SCATTER = 0.151  # (section 15.1) reference scatter on the 1st torsional frequency
DECK_TIP_FLOOR = 1.0e3  # GKt below median/1e3 is the degenerate last station (section 22.8)
WINDOW = (20.0, 100.0)  # [m] the trustworthy band; the ends and the z~112 tear are reported only


def _deck_gkt() -> tuple[np.ndarray, np.ndarray]:
    """The BeamDyn 6x6 deck -> ``(span_fraction, GKt)`` via ``K66toPropsDecoupled`` (section 22.8)."""
    lines = DECK.read_text(errors="replace").splitlines()
    start = next(i for i, ln in enumerate(lines) if "DISTRIBUTED PROPERTIES" in ln)
    rows: list[list[float]] = []
    for ln in lines[start + 1:]:
        body = ln.split("!", 1)[0].strip()
        if not body:
            continue
        try:
            rows.append([float(t) for t in body.split()])
        except ValueError:
            continue
    n = len(rows) // 13  # one span + 6x6 stiffness + 6x6 mass per station
    frac = np.array([rows[i * 13][0] for i in range(n)])
    K = np.array([np.asarray(rows[i * 13 + 1:i * 13 + 7]) for i in range(n)])
    return frac, np.array([K66toPropsDecoupled(K[i], convention="BeamDyn")[5] for i in range(n)])


def _outer_contour(mesh, coords, ring) -> list[int]:
    """The ordered, consistently oriented outer boundary cycle of one physical ring.

    An in-plane wall edge shared by two elements is interior (a web seam or a split prebent ring);
    an edge that occurs once is on the boundary.  Walking that boundary in one direction makes
    ``sum d = 0`` exact, so the shear flow built on it is self-equilibrated by construction.
    """
    ids = {int(i) for i in ring}
    index = mesh.node_id_to_index
    count: Counter = Counter()
    for element in mesh.elements:
        nodes = [index[nd.id] for nd in element.nodes]
        for a, b in zip(nodes, nodes[1:] + nodes[:1], strict=True):
            if a in ids and b in ids and abs(coords[a, 2] - coords[b, 2]) < STATION_GAP_TOLERANCE:
                count[tuple(sorted((a, b)))] += 1
    boundary = [edge for edge, n in count.items() if n == 1]
    adjacency: dict[int, list[int]] = {}
    for a, b in boundary:
        adjacency.setdefault(a, []).append(b)
        adjacency.setdefault(b, []).append(a)
    assert boundary and all(len(v) == 2 for v in adjacency.values()), (
        f"the ring boundary is not a simple closed cycle: {len(boundary)} edges, degrees "
        f"{sorted(len(v) for v in adjacency.values())}"
    )
    start = min(adjacency)
    cycle, prev, cur = [start], None, start
    while True:
        nxt = adjacency[cur][0] if adjacency[cur][0] != prev else adjacency[cur][1]
        if nxt == start:
            break
        cycle.append(nxt)
        prev, cur = cur, nxt
    assert len(cycle) == len(boundary), (len(cycle), len(boundary))
    return cycle


def _enclosed_area(cycle: list[int], coords: np.ndarray) -> float:
    """Hand-written shoelace area of the ordered contour [m^2], independent of the force vector."""
    x, y = coords[cycle, 0], coords[cycle, 1]
    return 0.5 * abs(float(np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y)))


def _contour_shear_flow(coords: np.ndarray, cycle: list[int], torque: float) -> np.ndarray:
    """The Bredt wall flow ``q = T/(2A)`` along one ring's contour, half to each edge end."""
    q = torque / (2.0 * _enclosed_area(cycle, coords))
    f = np.zeros(6 * len(coords))
    for a, b in zip(cycle, cycle[1:] + cycle[:1], strict=True):
        pa, pb = coords[a], coords[b]
        ell = float(np.hypot(pb[0] - pa[0], pb[1] - pa[1]))
        tx, ty = (pb[0] - pa[0]) / ell, (pb[1] - pa[1]) / ell
        for nd in (a, b):
            f[6 * nd] += 0.5 * q * ell * tx
            f[6 * nd + 1] += 0.5 * q * ell * ty
    return f


def _realised_torque(coords: np.ndarray, f: np.ndarray, centre=None) -> float:
    """The torque ruler ``sum((x-xc) Fy - (y-yc) Fx)`` [N.m] read back from the nodal system."""
    xc, yc = (0.0, 0.0) if centre is None else centre
    return float(np.sum((coords[:, 0] - xc) * f[1::6] - (coords[:, 1] - yc) * f[0::6]))


def _ring_section_rotation(coords: np.ndarray, u: np.ndarray, ring) -> float:
    """The antisymmetric part of the ring's in-plane affine fit [rad] (section 22.8)."""
    pts = coords[ring]
    design = np.column_stack([pts[:, 0], pts[:, 1], np.ones(len(ring))])
    a12 = np.linalg.lstsq(design, u[6 * ring], rcond=None)[0][1]
    a21 = np.linalg.lstsq(design, u[6 * ring + 1], rcond=None)[0][0]
    return float(0.5 * (a21 - a12))


@pytest.fixture(scope="module")
def static_torsion():
    """The shell mesh, one self-equilibrated pure-torque solve, and the per-ring section rotation."""
    if not DECK.exists():
        pytest.skip(f"BeamDyn deck not present: {DECK}")
    Node._id_counter = 0
    MeshElement._id_counter = 0
    model = Blade(str(YAML), element_size=1.0)
    model.generate_mesh()
    mesh, props = model.mesh, model.get_element_properties()
    assembler = PyMeshAssembler.from_model(
        blade_validation._to_rust_mesh(mesh, props), props,
        list(blade_validation.SPAN_DIRECTION), None)
    dofs = assembler.dofs_count
    rows, cols, vals = assembler.assemble_k()
    K = coo_matrix((np.asarray(vals), (np.asarray(rows), np.asarray(cols))),
                   shape=(dofs, dofs)).tocsr()
    coords = np.array([[nd.x, nd.y, nd.z] for nd in mesh.nodes], dtype=float)
    stations = _physical_stations(coords)
    rings = {zz: np.where(np.abs(coords[:, 2] - zz) < STATION_GAP_TOLERANCE)[0] for zz in stations}
    tip, root = stations[-1], stations[0]
    f_tip = _contour_shear_flow(coords, _outer_contour(mesh, coords, rings[tip]), +TORQUE)
    f_root = _contour_shear_flow(coords, _outer_contour(mesh, coords, rings[root]), -TORQUE)
    f = f_tip + f_root
    modes = tube._rigid_modes(coords)
    u = tube._solve_rigid_removed(K, f, modes)
    phi = np.array([_ring_section_rotation(coords, u, rings[zz]) for zz in stations])
    distortion = np.array([tube._parallelogram(coords, u, rings[zz])[2] for zz in stations])
    return {"coords": coords, "u": u, "f_tip": f_tip, "f_root": f_root, "f": f, "modes": modes,
            "stations": stations, "rings": rings, "tip": tip, "root": root, "dofs": dofs,
            "phi": phi, "distortion": distortion, "deck": _deck_gkt()}


def test_static_pure_torsion_sets_gj_per_station(static_torsion):
    """Assert the physics invariants, print the per-station ``GJ_static`` and the three routes."""
    coords, u, f, modes = (static_torsion[k] for k in ("coords", "u", "f", "modes"))
    f_tip, f_root = static_torsion["f_tip"], static_torsion["f_root"]
    stations = np.asarray(static_torsion["stations"], dtype=float)
    phi, rings, distortion = static_torsion["phi"], static_torsion["rings"], static_torsion["distortion"]

    # --- the torque ruler, read back from the applied nodal system (section 22.9's pattern) ---
    tip_centre = coords[rings[static_torsion["tip"]]][:, :2].mean(axis=0)
    ruler = {
        "tip about the origin": _realised_torque(coords, f_tip),
        "root about the origin": _realised_torque(coords, f_root),
        "tip about the tip centroid": _realised_torque(coords, f_tip, tip_centre),
    }
    print(f"\ntorque ruler T = {TORQUE:.6e} N.m: "
          + "; ".join(f"{k} {v:.10e} ({v / TORQUE:+.2e})" for k, v in ruler.items()))
    assert abs(ruler["tip about the origin"] - TORQUE) / TORQUE < RULER_TOL
    assert abs(ruler["root about the origin"] + TORQUE) / TORQUE < RULER_TOL
    assert abs(ruler["tip about the tip centroid"] - TORQUE) / TORQUE < RULER_TOL

    # --- self-equilibration: a pure couple carries no net force, so its torque is centre-free ---
    net_force = np.abs(f.reshape(-1, 6)[:, :3].sum(axis=0)).max()
    force_scale = np.abs(f.reshape(-1, 6)[:, :3]).max()
    resultants = np.array([abs(modes[a] @ f) / np.linalg.norm(f) for a in range(6)])
    print(f"net force {net_force:.3e} N against max nodal force {force_scale:.3e} N "
          f"({net_force / force_scale:.2e} relative); rigid-mode resultants "
          f"{' '.join(f'{v:.1e}' for v in resultants)}")
    assert net_force / force_scale < RULER_TOL
    assert resultants.max() < RULER_TOL, "the load is not self-equilibrated"
    # The rigid-mode gauge is a solver tolerance, not a physics invariant (the saddle-point system
    # is ill-scaled against K ~ 1e11).  It is reported: a residual rigid rotation about z shifts
    # every ring equally and cannot change a rate, and rigid rotations about x/y are invisible to
    # the z-section fit, so the measurement is gauge-independent.
    print(f"rigid-mode gauge residual |Phi^T u| / |u| = "
          f"{np.linalg.norm(modes @ u) / np.linalg.norm(u):.2e}")

    # --- sign: a positive torque about +span accumulates a positive, monotone twist ---
    assert np.all(np.diff(phi) > 0.0), "the twist does not accumulate monotonically along +span"
    assert phi[-1] > 0.0

    # --- route (c): the per-station static stiffness ---
    frac, gkt = static_torsion["deck"]
    r_deck = frac * float(stations[-1])
    rate = np.diff(phi) / np.diff(stations)
    gj_static = TORQUE / rate
    mid = 0.5 * (stations[1:] + stations[:-1])
    gj_ref = np.interp(mid, r_deck, gkt)
    ratio = gj_static / gj_ref
    print(f"\nroute (c): static pure torsion, q = T/(2A) at the tip and -q at the root, T = "
          f"{TORQUE:.3e} N.m; GJ_ref = deck GKt interpolated to r\n")
    print(f"  {'r[m]':>7} {'n':>3} {'GJ_static[N.m^2]':>16} {'GJ_ref[N.m^2]':>14} {'ratio':>7} "
          f"{'distort':>8}")
    for k in range(len(mid)):
        print(f"  {mid[k]:>7.2f} {len(rings[stations[k]]):>3} {gj_static[k]:>16.6e} "
              f"{gj_ref[k]:>14.6e} {ratio[k]:>7.3f} {distortion[k]:>8.3f}")

    # --- the three routes side by side and the branch they land in ---
    in_window = (mid >= WINDOW[0]) & (mid <= WINDOW[1])
    med = float(np.median(ratio[in_window]))
    refs = {"(a) deck at the root": 1.0, "(b) modal (0.870)": MODAL_RATIO,
            "(c) stiffer by 1.15-1.35": 1.25}
    nearest = min(refs, key=lambda k: abs(med - refs[k]))
    print(f"\nthree routes: (a) deck GKt root {GJ_REF_ROOT:.6e} N.m^2 (GJ_static/GJ_ref = 1 by "
          f"definition); (b) modal GJ_shell/GJ_ref {MODAL_RATIO:.3f} with {MODAL_SCATTER:.1%} "
          f"reference scatter (section 15.1, quoted);")
    print(f"(c) GJ_static/GJ_ref over r = {WINDOW[0]:.0f}-{WINDOW[1]:.0f} m: median {med:.3f}, "
          f"range {ratio[in_window].min():.3f}-{ratio[in_window].max():.3f}; nearest route "
          f"{nearest}")
    if abs(med - 1.0) <= MODAL_SCATTER:
        print("  branch 1: GJ_static ~ GJ_ref -> the shell's torsional stiffness is fine and the "
              "deficit would have to be in the delivered torque (contradicts section 22.7)")
    elif med > 1.0:
        print("  branch 2: the shell is STIFFER in static torsion than the deck and than the modal "
              "route says -> the section 15.1 modal ratio is the suspect input")
    else:
        print("  branch 3: static and modal agree -> the deficit is in the blade load/twist "
              "measurement of section 22.8, not in the stiffness")
    floor = float(np.median(gkt)) / DECK_TIP_FLOOR
    print(f"  deck degeneracy: median(GKt)/1e3 = {floor:.3e}; the last deck station GKt = "
          f"{gkt[-1]:.3e} is excluded, so GJ_ref stops being a section stiffness at r > "
          f"{r_deck[-2]:.2f} m")
    print(f"  distortion (ring residual after the best rigid rotation) at the window ends: "
          f"r={WINDOW[0]:.0f} m -> {distortion[int(np.argmin(np.abs(stations - WINDOW[0])))]:.3f}, "
          f"r={WINDOW[1]:.0f} m -> {distortion[int(np.argmin(np.abs(stations - WINDOW[1])))]:.3f}, "
          f"tip -> {distortion[-1]:.3f}")
