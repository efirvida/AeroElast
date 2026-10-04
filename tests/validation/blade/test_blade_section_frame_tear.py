"""The IEA-15MW section frame on the real mesh, and whether the z~112 m tear limits it.

Verdict section 22.12, on the repo's own blade mesh.

**The tear or the frame?**  Verdict 14.3 recorded free edges at ``z ~ 11.94 m`` and ``~ 112.22 m``;
section 22.10 excluded ``r >~ 111 m`` from its trustworthy band.  The discriminator is
tolerance-free: at the tear station the ring's chordwise extents, enclosed area and declared chord
must lie **inside the interval spanned by their two immediate neighbours** (a connectivity tear
would move them off the local trend), while the fitted section rotation's non-affine residual
(``distort``, section 22.10) is printed beside them as the estimator side of the discriminator.

**The frame, pinned once.**  ``x = (pitch_axis - f) * c``: ``x = 0`` is the pitch axis, the leading
edge sits at ``+pitch_axis * c``, the trailing edge at ``-(1 - pitch_axis) * c``; ``pitch_axis`` is a
chord fraction from the leading edge and agrees exactly between ``tests/IEA-15-240-RWT.yaml`` line 24
and column 2 of the ElastoDyn blade deck.  The blunt end (the larger in-plane spread inside the outer
quarter of the chord) is the leading edge at ``+x``; the declared load sense is fluid ``+Y`` /
rotor clockwise viewed from behind.  The frame is derived from the ring outline and the turbine
definition, never from the projector's configured axes; ``tests/validation/bem/
test_force_projection_load_frame.py`` is the P5 guard whose 10 deg axis bound the inter-station chord
angle must clear.  One shell solve (the section 22.10 pure-torque realisation) supplies ``distort``.
Reported, never asserted: section 21.3's 9.35 deg "sharp corner" does not reproduce as an in-plane
ring angle (the trailing edge is blunt, reading 102.2 deg), while its 97.3 deg companion reproduces.
"""
from __future__ import annotations

from collections import Counter

import numpy as np
import pytest
import yaml as _yaml

pytest.importorskip("_aeroelast", reason="Rust backend not available")

from _aeroelast import PyMeshAssembler  # noqa: E402
from scipy.sparse import coo_matrix  # noqa: E402
from scipy.spatial import ConvexHull  # noqa: E402

from aeroelast.core.mesh.entities import MeshElement, Node  # noqa: E402
from aeroelast.models.blade.model import Blade  # noqa: E402

import tests.validation.blade.test_blade_iea15mw_validation as blade_validation  # noqa: E402
import tests.validation.blade.test_blade_section_torsion_stiffness as torsion  # noqa: E402
import tests.validation.parity.test_thin_walled_tube_torsion as tube  # noqa: E402
from tests.validation.blade.test_blade_rated_twist import (  # noqa: E402
    STATION_GAP_TOLERANCE, _physical_stations)
from tests.validation.bem.test_force_projection_load_frame import (  # noqa: E402
    _angle_deg, _orient_leading_to_trailing)
from tests.support.paths import DATA_DIR, SOURCES_DIR  # noqa: E402

YAML = DATA_DIR / "IEA-15-240-RWT.yaml"
ELASTODYN = (SOURCES_DIR / "openfast" / "iea15mw" / "IEA-15-240-RWT"
             / "IEA-15-240-RWT_ElastoDyn_blade.dat")
L_BLADE = 117.0           # reference_axis z end [m]
PITCH_TOL = 1e-12         # yaml vs the deck: the same authored numbers, to round-off
SPLIT_TOL = 0.05          # the suite rule: pitch axis within 5 % of its declared fraction
TIE_REL = 0.10            # blunt/sharp tie: <=10 % spread difference is a tie
END_SLAB = 0.25           # outer quarter of the chord, per end
AXIS_BOUND_DEG = 10.0     # the P5 per-section axis bound (section 22.3)
TEAR_STATIONS = (11.9379, 112.2245)
SPAN = np.array([0.0, 0.0, 1.0])


def _clusters(values, gap=STATION_GAP_TOLERANCE):
    """Merge a value list into clusters no wider than ``gap`` -> list of value lists."""
    groups: list[list[float]] = []
    for value in np.sort(np.asarray(values, dtype=float)):
        if groups and value - groups[-1][-1] <= gap:
            groups[-1].append(float(value))
        else:
            groups.append([float(value)])
    return groups


def _edge_maps(mesh):
    """``(free_edges, adjacency)``: element edges used once, and node index -> neighbour indices."""
    index, count, adj = mesh.node_id_to_index, Counter(), {}
    for element in mesh.elements:
        ids = [index[nd.id] for nd in element.nodes]
        for a, b in zip(ids, ids[1:] + ids[:1], strict=True):
            count[tuple(sorted((a, b)))] += 1
            adj.setdefault(a, set()).add(b)
            adj.setdefault(b, set()).add(a)
    return [edge for edge, n in count.items() if n == 1], adj


def _declared_frame(stations):
    """Turbine definition interpolated to the stations: chord, pitch axis, twist, prebend."""
    outer = _yaml.safe_load(YAML.read_text())["components"]["blade"]["outer_shape_bem"]

    def at(key):
        return np.interp(stations, np.asarray(outer[key]["grid"], dtype=float) * L_BLADE,
                         np.asarray(outer[key]["values"], dtype=float))

    axis = outer["reference_axis"]["x"]
    return {"chord": at("chord"), "pitch_axis": at("pitch_axis"), "twist": at("twist"),
            "prebend": np.interp(stations, np.asarray(axis["grid"], dtype=float) * L_BLADE,
                                 np.asarray(axis["values"], dtype=float)),
            "grid": outer["pitch_axis"]["grid"], "values": outer["pitch_axis"]["values"]}


def _deck_pitch_axis():
    """``(span_fraction, pitch_axis)`` from the official ElastoDyn blade deck (column 2)."""
    lines = ELASTODYN.read_text(errors="replace").splitlines()
    start = next(i for i, ln in enumerate(lines) if "BlFract" in ln)
    rows = []
    for ln in lines[start + 2:]:
        try:
            rows.append([float(t) for t in ln.split("!", 1)[0].split()])
        except ValueError:
            continue
    deck = np.asarray(rows)
    return deck[:, 0], deck[:, 1]


def _blunt_at_plus_x(x, y):
    """Production blunt-end rule: blunt (LE) at ``+x`` unless the ``-x`` end is >10 % thicker."""
    slab = END_SLAB * (x.max() - x.min())
    t_lo, t_hi = float(np.ptp(y[x <= x.min() + slab])), float(np.ptp(y[x >= x.max() - slab]))
    return not (t_lo > t_hi and (t_lo - t_hi) > TIE_REL * max(t_lo, t_hi))


def _station(coords, blade, k):
    """Station ``k`` about its pitch-axis origin: ``(x, y, geometry dict, chord axis)``.

    The chord axis is the ring's principal in-plane axis (SVD), oriented by the blunt-end rule.
    """
    ring = blade["rings"][blade["stations"][k]]
    origin = np.array([0.0, blade["declared"]["prebend"][k]])
    pts = coords[ring][:, :2]
    axis = np.linalg.svd(pts - pts.mean(axis=0), full_matrices=False)[2][0]
    off = pts - origin
    c_hat = axis if _blunt_at_plus_x(off @ axis, off @ np.array([-axis[1], axis[0]])) else -axis
    prog = off @ c_hat
    thick = off @ np.array([-c_hat[1], c_hat[0]])
    centre = pts.mean(axis=0)
    order = np.argsort(np.arctan2(pts[:, 1] - centre[1], pts[:, 0] - centre[0]))
    px, py = pts[order, 0], pts[order, 1]
    geom = {"n": len(ring), "xmin": float(prog.min()), "xmax": float(prog.max()),
            "thick": float(np.ptp(thick)), "hull_area": float(ConvexHull(pts).volume),
            "shoelace": 0.5 * abs(float(np.sum(px * np.roll(py, -1) - np.roll(px, -1) * py)))}
    return prog, thick, geom, c_hat


def _corner_angle(coords, ring, pos, adj):
    """In-plane interior angle [deg] at a ring node with exactly two in-ring neighbours."""
    node, members = int(ring[pos]), {int(i) for i in ring}
    neighbours = [n for n in adj.get(node, ()) if n in members]
    if len(neighbours) != 2:
        return None
    v1, v2 = coords[neighbours[0]] - coords[node], coords[neighbours[1]] - coords[node]
    return float(np.rad2deg(np.arccos(np.clip(
        np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2)), -1.0, 1.0))))


@pytest.fixture(scope="module")
def blade():
    """Real mesh + declared frame + one pure-torque solve (the section 22.10 realisation)."""
    if not YAML.exists() or not ELASTODYN.exists():
        pytest.skip("IEA-15-240-RWT definition or ElastoDyn deck not present")
    Node._id_counter = 0
    MeshElement._id_counter = 0
    model = Blade(str(YAML), element_size=1.0)
    model.generate_mesh()
    mesh, props = model.mesh, model.get_element_properties()
    coords = np.array([[nd.x, nd.y, nd.z] for nd in mesh.nodes], dtype=float)
    stations = np.asarray(_physical_stations(coords), dtype=float)
    rings = {zz: np.where(np.abs(coords[:, 2] - zz) < STATION_GAP_TOLERANCE)[0] for zz in stations}
    free, adj = _edge_maps(mesh)
    assembler = PyMeshAssembler.from_model(
        blade_validation._to_rust_mesh(mesh, props), props,
        list(blade_validation.SPAN_DIRECTION), None)
    rows, cols, vals = assembler.assemble_k()
    K = coo_matrix((np.asarray(vals), (np.asarray(rows), np.asarray(cols))),
                   shape=(assembler.dofs_count,) * 2).tocsr()
    f = torsion._contour_shear_flow(coords, torsion._outer_contour(mesh, coords, rings[stations[-1]]),
                                    +torsion.TORQUE)
    f += torsion._contour_shear_flow(coords, torsion._outer_contour(mesh, coords, rings[stations[0]]),
                                     -torsion.TORQUE)
    u = tube._solve_rigid_removed(K, f, tube._rigid_modes(coords))
    del K, f
    return {"mesh": mesh, "coords": coords, "stations": stations, "rings": rings, "free": free,
            "adj": adj, "declared": _declared_frame(stations), "u": u,
            "phi": np.array([torsion._ring_section_rotation(coords, u, rings[zz]) for zz in stations]),
            "distort": np.array([tube._parallelogram(coords, u, rings[zz])[2] for zz in stations])}


def test_section_frame_invariants_per_station(blade):
    """Pin the section frame on every ring, and the yaml/deck pitch axis that defines it.

    ``x = 0`` is the pitch axis and the extents about it must be ``[-(1 - pitch_axis) c,
    + pitch_axis c]``; the blunt (leading) end must be at ``+x``.
    """
    coords, stations, d = blade["coords"], blade["stations"], blade["declared"]
    grid, values = np.asarray(d["grid"], dtype=float), np.asarray(d["values"], dtype=float)
    deck_err = float(np.max(np.abs(np.interp(grid, *_deck_pitch_axis()) - values)))
    print(f"\npitch_axis yaml line 24 vs the ElastoDyn deck column 2: max |diff| = {deck_err:.3e} "
          f"over {len(grid)} stations")
    assert deck_err < PITCH_TOL, "the yaml pitch_axis is not the ElastoDyn deck's PitchAxis column"

    split, chord_dev, twist_dev, ratios, xmins, xmaxs, blunts = [], [], [], [], [], [], []
    for k, zz in enumerate(stations):
        pts = coords[blade["rings"][zz]][:, :2]
        ratios.append(np.linalg.svd(pts - pts.mean(axis=0), compute_uv=False).tolist())
        x, y, geom, c_hat = _station(coords, blade, k)
        extent = geom["xmax"] - geom["xmin"]
        xmins.append(geom["xmin"])
        xmaxs.append(geom["xmax"])
        split.append(max(abs(geom["xmax"] / extent - d["pitch_axis"][k]),
                         abs(-geom["xmin"] / extent - (1.0 - d["pitch_axis"][k]))))
        chord_dev.append(abs(extent / d["chord"][k] - 1.0))
        twist = np.array([np.cos(-d["twist"][k]), np.sin(-d["twist"][k])])
        twist_dev.append(float(np.rad2deg(np.arccos(np.clip(c_hat @ twist, -1.0, 1.0)))))
        blunts.append(_blunt_at_plus_x(x, y))
    split, chord_dev, twist_dev = map(np.asarray, (split, chord_dev, twist_dev))
    defined = np.asarray(ratios)[:, 0] > 1.2 * np.asarray(ratios)[:, 1]  # near-circular: no axis
    print(f"pitch-axis split: median {np.median(split):.3%}, max {split.max():.3%} at "
          f"z={stations[np.argmax(split)]:.3f}; declared-chord extent: median "
          f"{np.median(chord_dev):.3%}, max {chord_dev.max():.1%} at "
          f"z={stations[np.argmax(chord_dev)]:.3f} ({int(np.sum(chord_dev > 0.05))} stations >5 %); "
          f"ring axis vs declared twist: median {np.median(twist_dev):.2f} deg, max "
          f"{twist_dev[defined].max():.2f} deg over {int(defined.sum())} stations (near-circular "
          f"root stations without a defined axis: {[round(z, 3) for z in stations[~defined]]}); "
          f"blunt at +x on all: {all(blunts)}")

    assert np.all(np.asarray(xmins) < 0.0) and np.all(np.asarray(xmaxs) > 0.0), (
        "the pitch axis is not strictly inside the chordwise extent on some ring")
    assert split.max() <= SPLIT_TOL, (
        f"the pitch axis is not at the declared chord fraction: max {split.max():.3%} at "
        f"z={stations[np.argmax(split)]:.3f}")
    assert np.all(twist_dev[defined] <= AXIS_BOUND_DEG) and all(blunts), (
        f"the mesh chord axis departs from the declared twist by up to {twist_dev[defined].max():.2f}"
        f" deg (bound {AXIS_BOUND_DEG:.0f}), or the blunt (LE) end is at -x")


def test_blunt_sharp_corners_and_the_global_axis(blade):
    """The mesh's own axis vs the P5 bound, and section 21.3's recorded blunt/sharp kink angles."""
    coords, stations = blade["coords"], blade["stations"]

    def axis3(zz):
        pts = coords[blade["rings"][zz]]
        off = pts - pts.mean(axis=0)
        off -= np.outer(off @ SPAN, SPAN)
        axis = np.linalg.svd(off, full_matrices=False)[2][0]
        return _orient_leading_to_trailing(pts, axis / np.linalg.norm(axis), SPAN)

    chosen = [min(stations, key=lambda s: abs(s - stations[0] - f * (stations[-1] - stations[0])))
              for f in (0.03, 0.15, 0.25, 0.50, 0.75, 0.90, 1.00)]
    angles = [(zz, axis3(zz)) for zz in chosen]
    best = max((_angle_deg(ci, cj), (zi, zj)) for i, (zi, ci) in enumerate(angles)
               for zj, cj in angles[i + 1:])
    print(f"\ninter-station chord direction: max {best[0]:.3f} deg between z={best[1][0]:.3f} and "
          f"z={best[1][1]:.3f} m (bound {2 * AXIS_BOUND_DEG:.0f} deg = 2 x {AXIS_BOUND_DEG:.0f}); "
          f"root-ish c_hat={np.array2string(angles[0][1], precision=5)}, "
          f"tip c_hat={np.array2string(angles[-1][1], precision=5)}")
    assert best[0] >= 2.0 * AXIS_BOUND_DEG, (
        f"the chosen chords differ by only {best[0]:.3f} deg; no single global axis is excluded")

    k = int(np.argmin(np.abs(stations - 68.051020408)))
    x, y, _, _ = _station(coords, blade, k)
    ring = blade["rings"][stations[k]]
    slab = END_SLAB * (x.max() - x.min())
    print(f"corners at z={stations[k]:.4f} m: outer-quarter in-plane spread "
          f"{np.ptp(y[x <= x.min() + slab]):.4f} m at min-x (TE) vs "
          f"{np.ptp(y[x >= x.max() - slab]):.4f} m at max-x (LE, blunt); interior angles "
          f"{_corner_angle(coords, ring, int(np.argmax(x)), blade['adj']):.2f} deg (section 21.3's "
          f"97.3 reproduces) and {_corner_angle(coords, ring, int(np.argmin(x)), blade['adj']):.2f} "
          f"deg (not 9.35: the trailing edge is blunt)")
    assert _blunt_at_plus_x(x, y), "the blunt (leading) end is not at +x"

    tip = axis3(stations[-1])
    print(f"LE/TE sense at the tip: c_hat={np.array2string(tip, precision=5)}, "
          f"c_hat.[1,0,0]={tip[0]:+.5f} (>0: mesh +x is the LE direction), "
          f"(span x c_hat).[0,1,0]={np.cross(SPAN, tip)[1]:+.5f} (>0: section normal is +Y = the "
          f"fluid direction)")
    assert tip[0] > 0.0 and np.cross(SPAN, tip)[1] > 0.0, (
        "the tip chord axis is not toward +x (LE) with its normal along +Y (the declared sense)")


def test_the_grid_and_the_tear_do_not_move_the_section_geometry(blade):
    """Report the station grid and the free-edge tears, then apply the tear discriminator.

    A connectivity tear would move the ring's chordwise extents, enclosed area or chord off the
    local trend; a purely estimator/frame limit would leave them and move only ``distort``.  The
    bracket assertion is tolerance-free: the tear value must lie inside the interval spanned by its
    two immediate neighbours.
    """
    coords, stations = blade["coords"], blade["stations"]
    spacing = np.diff(stations)
    sizes = np.array([len(blade["rings"][zz]) for zz in stations])
    edge_z = np.array([0.5 * (coords[a, 2] + coords[b, 2]) for a, b in blade["free"]])
    clusters = _clusters(edge_z)
    interior = [g for g in clusters
                if STATION_GAP_TOLERANCE < np.mean(g) < stations[-1] - STATION_GAP_TOLERANCE]
    print(f"\nmerged physical stations: {len(stations)} "
          f"(from {len(_clusters(np.unique(np.round(coords[:, 2], 9))))} raw z buckets); span "
          f"{stations[0]:.4f}..{stations[-1]:.4f} m; spacing min {spacing.min():.4f} / median "
          f"{np.median(spacing):.4f} / max {spacing.max():.4f} m; ring nodes min {sizes.min()} / max "
          f"{sizes.max()}\nfree edges: {len(blade['free'])} -> "
          + ", ".join(f"{len(g)} at z={np.mean(g):.4f}" for g in clusters)
          + " | tear clusters: " + ", ".join(f"z={np.mean(g):.4f} n={len(g)}" for g in interior))
    for anchor in TEAR_STATIONS:
        for k in [k for k, zz in enumerate(stations) if abs(zz - anchor) < 2.0]:
            zz = stations[k]
            print(f"  z={zz:9.4f} gap={float(stations[k] - stations[k - 1]) if k else 0.0:.4f} "
                  f"ring_n={len(blade['rings'][zz])} "
                  f"buckets={len(_clusters(np.unique(np.round(coords[blade['rings'][zz], 2], 9))))} "
                  f"free={int(np.sum(np.abs(edge_z - zz) < STATION_GAP_TOLERANCE))}")
    assert len(interior) == 2 and len(interior[0]) == 8 and len(interior[1]) == 2, (
        "the recorded tear signature (8 free edges at z~11.94, 2 at z~112.22) is gone; update "
        "verdict sections 14.3 / 22.10 and this guard")
    assert np.allclose([np.mean(g) for g in interior], TEAR_STATIONS, atol=2e-3)
    assert np.all(spacing > 0.0), "physical stations are not strictly increasing"

    rows = [{"z": zz, "chord": blade["declared"]["chord"][k], "distort": blade["distort"][k],
             "geom": _station(coords, blade, k)[2]}
            for k, zz in enumerate(stations) if zz >= 108.0]
    print("outer band: raw-node geometry vs the fitted section rotation")
    print(f"  {'z[m]':>9} {'n':>3} {'extent':>8} {'xmin':>8} {'xmax':>8} {'thick':>7} "
          f"{'hull A':>8} {'shoelace':>8} {'c':>7} {'distort':>8}")
    for row in rows:
        g = row["geom"]
        print(f"  {row['z']:>9.4f} {g['n']:>3} {g['xmax'] - g['xmin']:>8.4f} {g['xmin']:>8.4f} "
              f"{g['xmax']:>8.4f} {g['thick']:>7.4f} {g['hull_area']:>8.5f} {g['shoelace']:>8.5f} "
              f"{row['chord']:>7.4f} {row['distort']:>8.3f}")
    tear = int(np.argmin(np.abs(np.asarray([r["z"] for r in rows]) - TEAR_STATIONS[1])))
    assert abs(rows[tear]["z"] - TEAR_STATIONS[1]) < STATION_GAP_TOLERANCE
    before, after = rows[tear - 1], rows[tear + 1]
    quantities = {"chordwise extent": lambda r: r["geom"]["xmax"] - r["geom"]["xmin"],
                  "xmin": lambda r: r["geom"]["xmin"], "xmax": lambda r: r["geom"]["xmax"],
                  "hull area": lambda r: r["geom"]["hull_area"],
                  "shoelace area": lambda r: r["geom"]["shoelace"],
                  "declared chord": lambda r: r["chord"]}
    print(f"\ntear station z={rows[tear]['z']:.4f} m between z={before['z']:.4f} and "
          f"z={after['z']:.4f}; neighbour bracket:")
    for name, fn in quantities.items():
        value, low, high = fn(rows[tear]), fn(before), fn(after)
        inside = min(low, high) <= value <= max(low, high)
        print(f"  {name:>16}: prev {low:.5f}  tear {value:.5f}  next {high:.5f}  -> "
              f"{'inside' if inside else 'OUTSIDE the bracket'}")
        assert inside, (f"{name} jumps at the z={rows[tear]['z']:.3f} m tear: {value:.5f} outside "
                        f"[{min(low, high):.5f}, {max(low, high):.5f}]")
    print(f"  {'fitted distort':>16}: prev {before['distort']:.3f}  tear {rows[tear]['distort']:.3f}  "
          f"next {after['distort']:.3f} (smooth: the tear does not disturb the affine fit)")
    assert np.all(np.diff(blade["phi"]) > 0.0), "the section rotation does not accumulate"
