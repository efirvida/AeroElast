"""Does the two-realisation divergence of a section moment depend on section SHAPE?

WU: one physics question, one script.  A BEM strip's span couple is realised on the shell
nodes in two ways:

* **wall flow** (physical, production) - the Bredt shear flow ``q = T/(2A)`` along the closed
  wall, through :func:`aeroelast.solvers.bem.force_projection.realise_section_load`;
* **minimum-norm** (legacy default) - the constrained minimum-norm nodal field
  :meth:`ForceProjector._distribute`, whose closed form for a pure span couple is
  ``f_j = omega x d_j`` (magnitude ~ ``|d_j|``, direction perpendicular to the radius, not
  tangent to the wall).

Hypothesis under test: on a **circular** ring the two coincide; on an **elongated or
non-circular** ring they do not (minimum-norm starves the mid-side walls and dumps the couple
near the far corners), so the divergence is a function of section **shape**, not of FSI.

What this script measures
-------------------------
**Experiment 1 - CONTROL (mandatory, first).**  The already-recorded numbers on the existing
rectangular thin-walled tube fixture are reproduced with the fixture's own helpers imported,
not reimplemented.  The record (``odd/tasks/composite-bend-twist-verdict.md`` section 22.9) is
``rate/Bredt = 1.0031`` for the wall flow and ``31.6933`` for the minimum-norm field, both on
the self-equilibrated (+T tip / -T root, rigid modes removed) solve.  A clamped-tip variant is
also printed for completeness (the record notes ``8.9x`` there).

**Experiment 2 - TREATMENT.**  The same experiment on the deck's own airfoil section: the real
outer-skin ring of ``tests/IEA-15-240-RWT.yaml`` at a mid-span station, extruded into a short
tube, clamped at the root, loaded with the same pure span torque at the tip through BOTH
realisations, rate measured from interior stations over a window away from both ends.  The
reference is Bredt ``GJ = 4 A^2 / (integral ds/A66)`` of that ring's own wall graph with each
edge's ``A66`` resolved from the deck's own laminate (:func:`_membrane_shear_stiffness`).

For both shapes and both realisations the script prints: the realised resultant (sum F, full
moment vector, span component), the interior twist rate and its ratio to Bredt, and an
explicit in-plane distortion metric (RMS of the displacement minus its best-fit rigid
translation+rotation over an interior ring, both absolute and normalised by the section
dimension).  Distortion is the quantity the hypothesis is about, so it is measured.

Diagnostic only: no assertions against tolerances, not collected by the test suite, writes
nothing.  Run with the repository bootstrap, e.g.::

    scripts/aeroenv.sh python tools/diagnose_section_moment_realization.py
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix

# A script puts its own directory on sys.path, not the repository root, so `import tests.*`
# below fails when this is run as `python tools/diagnose_section_moment_realization.py`.
_REPO_ROOT = Path(__file__).resolve().parents[1]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from _aeroelast import PyMeshAssembler  # noqa: E402
from aeroelast.core.mesh.entities import ElementType  # noqa: E402
from aeroelast.core.mesh.generators import BladeMesh  # noqa: E402
from aeroelast.models.blade.model import build_rust_properties  # noqa: E402
from aeroelast.solvers.bem.force_projection import (  # noqa: E402
    RING_GAP_FRACTION,
    ForceProjector,
    _element_property_by_id,
    _membrane_shear_stiffness,
    _RingSection,
    _Strip,
    realise_section_load,
    ring_section,
)
from aeroelast.solvers.bem.section_contour import (  # noqa: E402
    SectionCell,
    cell_adjacency,
    contour_report,
    rings_in_band,
    section_plane_axes,
)
from tests.support.paths import DATA_DIR  # noqa: E402

import tests.validation.parity.test_thin_walled_tube_moment_realization as mr  # noqa: E402
import tests.validation.parity.test_thin_walled_tube_torsion as tube  # noqa: E402

# ─────────────────────────────────────────────────────────────────────────────
# Shared helpers (the estimators come from the tube fixture; these are new here)
# ─────────────────────────────────────────────────────────────────────────────

SPAN = np.array([0.0, 0.0, 1.0])
BLADE_YAML = DATA_DIR / "IEA-15-240-RWT.yaml"
AIRFOIL_DIR = DATA_DIR / "airfoils"
STATION_TARGET = 40.0  # m along the blade, mid-span-ish (r/R ~ 0.34 of a 117 m blade)
ELEMENT_SIZE = 0.25
TORQUE = 1.0e4  # N.m, the tube fixture's own value
L_OVER_DIM = 6.0  # span / max section dimension, matching the rectangle (dim 1 m, L 6 m)
N_Z_TREATMENT = 40


def _resultants(
    f: np.ndarray, coords: np.ndarray, span: np.ndarray
) -> tuple[np.ndarray, np.ndarray, float]:
    """Net force, net moment about the origin, and the span component of the moment."""
    fv = f.reshape(-1, 6)[:, :3]
    force = fv.sum(axis=0)
    moment = np.cross(coords, fv).sum(axis=0)
    return force, moment, float(moment @ span)


def _distortion(coords: np.ndarray, u: np.ndarray, ring: list[int]) -> tuple[float, float, float]:
    """In-plane distortion of one ring: RMS of ``u`` minus its best-fit rigid translation+rotation.

    Construction: for the ring's translations ``(ux, uy)`` and in-plane coordinates ``(xc, yc)``
    measured from the ring centroid, remove the mean translation ``(tx, ty)`` and the
    least-squares rotation ``theta = sum(xc*ay - yc*ax) / sum(xc^2 + yc^2)`` (``ax, ay`` the
    translation-removed field).  Returns

    ``(rms_residual, rms_residual / section_dimension, theta)``

    with ``section_dimension = max(ptp(xc), ptp(yc))``.  The normalised value is dimensionless
    and comparable between sections; the rotation is reported so rigid twist is separated out.
    """
    idx = np.asarray(ring, dtype=int)
    xy = coords[idx]
    c = xy.mean(axis=0)
    xc = xy[:, 0] - c[0]
    yc = xy[:, 1] - c[1]
    ux = u[6 * idx]
    uy = u[6 * idx + 1]
    tx = float(ux.mean())
    ty = float(uy.mean())
    ax = ux - tx
    ay = uy - ty
    den = float(np.sum(xc**2 + yc**2))
    theta = float(np.sum(xc * ay - yc * ax) / den) if den > 0.0 else 0.0
    res = np.sqrt(np.mean((ux - (tx - theta * yc)) ** 2 + (uy - (ty + theta * xc)) ** 2))
    dim = float(max(np.ptp(xc), np.ptp(yc)))
    return float(res), float(res / dim) if dim > 0.0 else float("nan"), theta


def _embed(nodes_local: np.ndarray, level_ids: list[int], n_nodes: int) -> np.ndarray:
    """Embed an ``(n_ring, 3)`` nodal force field onto one mesh ring in the 6-DOF layout."""
    f = np.zeros(6 * n_nodes)
    for i, nd in enumerate(level_ids):
        f[6 * nd : 6 * nd + 3] = nodes_local[i]
    return f


def _slope_report(label, coords, u, rings, ref_rate, window):
    """theta_fit rate over ``window``, its ratio to Bredt, and the interior distortion."""
    z = np.asarray([coords[r[0], 2] for r in rings])
    theta = np.asarray([tube._theta_fit(coords, u, r) for r in rings])
    rate = tube._lsq_slope(z, theta, *window)
    station = int(0.75 * (len(rings) - 1))
    rms, norm, th = _distortion(coords, u, rings[station])
    print(
        f"    {label:22s} rate={rate: .6e} rad/m  rate/Bredt={rate / ref_rate: .5f}"
        f"   distortion: rms/dim={norm:.4e} rms={rms:.3e} theta={th: .3e}"
    )
    return rate, rate / ref_rate, norm


# ─────────────────────────────────────────────────────────────────────────────
# Experiment 1 - CONTROL on the rectangular thin-walled tube
# ─────────────────────────────────────────────────────────────────────────────


def experiment_1_control() -> None:
    print("=" * 96)
    print("EXPERIMENT 1 - CONTROL: rectangular thin-walled tube (fixture helpers imported)")
    print("=" * 96)
    coords, conn, rings, n_ring = tube._tube_mesh()
    _, K = tube._assemble(coords, conn, 4, tube._iso_prop())
    modes = tube._rigid_modes(coords)
    gj, ref_rate = tube._bredt_isotropic()
    window = (0.4 * tube.L, 0.9 * tube.L)
    print(
        f"  tube: B={tube.B} H={tube.H} L={tube.L} t={tube.THICKNESS}  ring nodes={n_ring}  "
        f"T={tube.TORQUE:.3e} N.m"
    )
    print(f"  Bredt: GJ={gj:.6e} N.m^2  theta'=T/GJ={ref_rate:.6e} rad/m  window={window}")

    # ---- self-equilibrated (+T tip, -T root), rigid modes removed - the record's configuration
    print("\n  [self-equilibrated +T tip / -T root, rigid modes removed - the 22.9 record]")
    cases = {
        "wall flow (hand)": tube._ring_shear_flow(coords, rings[-1], tube.TORQUE)
        + tube._ring_shear_flow(coords, rings[0], -tube.TORQUE),
        "wall flow (prod.)": mr._production_tip_moment(coords, rings[-1], +tube.TORQUE)
        + mr._production_tip_moment(coords, rings[0], -tube.TORQUE),
        "minimum-norm": _min_norm_tip(coords, rings[-1], +tube.TORQUE)
        + _min_norm_tip(coords, rings[0], -tube.TORQUE),
    }
    for label, f in cases.items():
        F, M, mz = _resultants(f, coords, SPAN)
        u = tube._solve_rigid_removed(K, f, modes)
        print(
            f"    {label:18s} sum F={np.round(F, 12).tolist()}  M={np.round(M, 6).tolist()}"
            f"  span={mz:.6e}"
        )
        _slope_report(label, coords, u, rings, ref_rate, window)

    # ---- clamped tip, for the record's own 8.9x note and to bridge to experiment 2
    print("\n  [clamped root, +T tip only - the same BC experiment 2 uses]")
    clamped = tube._clamped_dofs(rings)
    for label, f in (
        ("wall flow (prod.)", mr._production_tip_moment(coords, rings[-1], +tube.TORQUE)),
        ("minimum-norm", _min_norm_tip(coords, rings[-1], +tube.TORQUE)),
    ):
        F, M, mz = _resultants(f, coords, SPAN)
        u = tube._solve_clamped(K, f, clamped)
        print(
            f"    {label:18s} sum F={np.round(F, 12).tolist()}  M={np.round(M, 6).tolist()}"
            f"  span={mz:.6e}"
        )
        _slope_report(label, coords, u, rings, ref_rate, window)


def _min_norm_tip(coords: np.ndarray, ring: list[int], torque: float) -> np.ndarray:
    """The minimum-norm realisation of a pure span couple about one ring's centroid."""
    strip = mr._ring_strip(coords, ring)
    nodes = ForceProjector._distribute(strip, np.zeros(3), np.asarray([0.0, 0.0, torque]))
    f = np.zeros(6 * len(coords))
    for j, nd in enumerate(ring):
        f[6 * nd : 6 * nd + 3] = nodes[j]
    return f


# ─────────────────────────────────────────────────────────────────────────────
# Experiment 2 - TREATMENT on the deck's own airfoil section
# ─────────────────────────────────────────────────────────────────────────────


def _edge_owner_map(mesh, ring_global: np.ndarray, props: dict) -> dict:
    """Map every undirected ring edge (local indices) to its owning mesh element's deck property.

    Mirrors :func:`ring_section`'s own enumeration exactly (corner cycle of tri/quad families,
    only edges whose two endpoints are both in the ring), so the per-edge wall property is the
    one the production wall-flow solver would read.
    """
    tri = (ElementType.triangle, ElementType.triangle6)
    quad = (ElementType.quad, ElementType.quad8, ElementType.quad9)
    ring = np.asarray(ring_global, dtype=np.intp)
    local_of_global = {int(g): i for i, g in enumerate(ring.tolist())}
    element_property = _element_property_by_id(mesh, props)
    owner: dict[tuple[int, int], object] = {}
    for element in mesh.elements:
        et = element.element_type
        if et in tri:
            corner_count = 3
        elif et in quad:
            corner_count = 4
        else:
            continue
        indices = [mesh.node_id_to_index[nid] for nid in element.node_ids]
        corners = indices[:corner_count]
        for a, b in zip(corners, corners[1:] + corners[:1], strict=True):
            if a not in local_of_global or b not in local_of_global or a == b:
                continue
            la, lb = local_of_global[a], local_of_global[b]
            key = (la, lb) if la < lb else (lb, la)
            owner.setdefault(key, element_property.get(int(element.id)))
    return owner


def _generic_case(
    label: str,
    base: np.ndarray,
    cells: list,
    prop_for_edge: dict,
    edge_a66: dict,
    length: float,
    nz: int,
) -> None:
    """Extruded-tube measurement of both realisations against Bredt of the ring's own wall graph.

    ``base`` is the ring in section-plane coordinates (centroid at the origin), ``cells`` the
    closed face(s) of its wall graph in local node indices, ``prop_for_edge`` / ``edge_a66`` the
    per-edge isotropic substitute dict and its A66 [N/m].  The reference is the single-cell Bredt
    ``GJ = 4 A^2 / integral ds/A66`` of the same wall graph.  Both realisations are applied at the
    tip ring (clamped root and self-equilibrated variants), and the twist rate is read from
    interior stations over ``(0.4L, 0.9L)``.
    """
    n_ring = len(base)
    dz = length / nz
    adjacency = cell_adjacency(cells)
    area = float(sum(float(cell.area) for cell in cells))
    edge_length: dict[tuple[int, int], float] = {}
    wall_edges: list[tuple[int, int]] = []
    directed_edges: list[tuple[int, int]] = []
    for cell in cells:
        for a, b in cell.edges:
            directed_edges.append((int(a), int(b)))
            key = (min(int(a), int(b)), max(int(a), int(b)))
            if key not in edge_length:
                edge_length[key] = float(np.linalg.norm(base[key[0]] - base[key[1]]))
                wall_edges.append(key)
    fallback_key = next(iter(prop_for_edge))
    for key in wall_edges:
        edge_a66.setdefault(key, edge_a66[fallback_key])
        prop_for_edge.setdefault(key, prop_for_edge[fallback_key])
    harmonic = sum(edge_length[k] / edge_a66[k] for k in wall_edges)
    gj = 4.0 * area**2 / harmonic
    ref_rate = TORQUE / gj
    print(
        f"  [{label}] wall graph: A={area:.4f} m^2  perimeter={sum(edge_length.values()):.4f} m  "
        f"GJ={gj:.6e} N.m^2  theta'=T/GJ={ref_rate:.6e} rad/m"
    )

    # ---- extruded tube
    n_nodes = (nz + 1) * n_ring
    tube_coords = np.zeros((n_nodes, 3))
    levels: list[list[int]] = []
    for k in range(nz + 1):
        ids = []
        for i in range(n_ring):
            nid = len(ids) + k * n_ring
            tube_coords[nid] = base[i] + (k * dz) * SPAN
            ids.append(nid)
        levels.append(ids)
    conn: list[list[int]] = []
    mats: list[dict] = []
    for k in range(nz):
        for a, b in directed_edges:
            conn.append([levels[k][a], levels[k][b], levels[k + 1][b], levels[k + 1][a]])
            mats.append(prop_for_edge[(min(a, b), max(a, b))])
    asm = PyMeshAssembler(tube_coords, conn, [4] * len(conn), mats)
    rows, cols, vals = asm.assemble_k()
    K = coo_matrix((vals, (rows, cols)), shape=(asm.dofs_count, asm.dofs_count)).tocsr()
    clamped = [6 * nd + d for nd in levels[0] for d in range(6)]
    modes = tube._rigid_modes(tube_coords)
    window = (0.4 * length, 0.9 * length)
    print(
        f"    extruded tube: {n_nodes} nodes, {len(conn)} elements, {asm.dofs_count} dof, "
        f"L={length:.3f} dz={dz:.3f} window={window}"
    )

    strip = _Strip(
        node_indices=np.arange(n_ring, dtype=np.intp),
        r_center=0.0,
        dr=dz,
        centroid=base.mean(axis=0),
        offsets=base - base.mean(axis=0),
    )
    section = _RingSection(
        cells,
        adjacency,
        edge_length,
        {k: edge_a66[k] for k in wall_edges},
        np.arange(n_ring, dtype=np.intp),
        from_element_properties=True,
    )
    wall_nodes = realise_section_load(
        strip, np.zeros(3), np.asarray([0.0, 0.0, TORQUE]), SPAN, ring_sections=[section]
    )
    min_nodes = ForceProjector._distribute(strip, np.zeros(3), np.asarray([0.0, 0.0, TORQUE]))
    for name, nodal in (("wall flow (prod.)", wall_nodes), ("minimum-norm", min_nodes)):
        f_tip = _embed(nodal, levels[-1], n_nodes)
        force, moment, mz = _resultants(f_tip, tube_coords, SPAN)
        print(
            f"    [{name}] tip load: sum F={np.round(force, 10).tolist()}  "
            f"M={np.round(moment, 6).tolist()}  span={mz:.6e}"
        )
        print("      clamped root:")
        _slope_report(
            name, tube_coords, tube._solve_clamped(K, f_tip, clamped), levels, ref_rate, window
        )
        f_eq = f_tip + _embed(-nodal, levels[0], n_nodes)
        print("      self-equilibrated (rigid modes removed):")
        _slope_report(
            name, tube_coords, tube._solve_rigid_removed(K, f_eq, modes), levels, ref_rate, window
        )


def experiment_2_treatment() -> None:
    print()
    print("=" * 96)
    print("EXPERIMENT 2 - TREATMENT: deck airfoil outer-skin ring, extruded tube")
    print("=" * 96)
    t0 = time.time()
    generator = BladeMesh(
        yaml_file=str(BLADE_YAML),
        airfoil_dir=str(AIRFOIL_DIR),
        element_size=ELEMENT_SIZE,
        n_samples=300,
        airfoil_spacing="constant",
        span_grading="chord",
    )
    mesh = generator.generate(renumber="rcm", verbose=False)
    props = build_rust_properties(generator.numad_mesh_data)
    coords = mesh.coords_array
    print(
        f"  blade mesh: {coords.shape[0]} nodes, {len(mesh.elements)} elements, "
        f"element_size={ELEMENT_SIZE}  ({time.time() - t0:.1f} s)"
    )
    print(f"  node sets available: {sorted(mesh.node_sets.keys())[:12]}")

    set_name = "allOuterShellNods"
    if set_name not in mesh.node_sets:
        print(
            f"  OBSTACLE: node set '{set_name}' not present on the generated mesh; cannot continue."
        )
        return
    skin_ids = sorted(mesh.node_sets[set_name].nodes.keys())
    skin_idx = np.asarray([mesh.node_id_to_index[nid] for nid in skin_ids], dtype=np.intp)
    zspan = coords @ SPAN
    gap = RING_GAP_FRACTION * float(zspan.max() - zspan.min())
    groups = rings_in_band(zspan[skin_idx], gap)
    means = np.asarray([float(zspan[skin_idx[g]].mean()) for g in groups])
    gi = int(np.argmin(np.abs(means - STATION_TARGET)))
    ring_global = skin_idx[groups[gi]]
    ring_pts = coords[ring_global]
    station_z = float(means[gi])
    print(
        f"  station used: z={station_z:.3f} m (target {STATION_TARGET} m; the nearest of "
        f"{len(groups)} skin rings).  ring nodes={len(ring_global)}"
    )

    # ---- the ring's own wall graph, read from the production helper, so every wall edge is a
    #      real mesh edge and every wall carries its own laminate/A66.
    cells, _, _, _ = ring_section(mesh, ring_global, SPAN, props, coords=coords)
    print(f"  ring_section: {len(cells)} closed cell(s)")
    report = contour_report(ring_pts, SPAN)
    print(
        f"  contour_report (independent ordering check): usable={report.usable} "
        f"simple={report.simple} area={report.area:.4f} m^2 perimeter={report.perimeter:.4f} m"
    )

    # ---- per-edge deck property (for the FEM) and A66 (for the reference)
    #
    # The deck's own laminate is defined in the blade's layup frame (its ``xDir`` is the span,
    # ``xyDir`` the chord), not in the shell element's local frame.  Feeding its raw ``cm`` to a
    # hand-built tube would rotate A66 relative to the wall tangent and the reference would no
    # longer be the material the element sees.  Since the hypothesis is about section SHAPE,
    # the orientation confound is removed with an isotropic substitute per wall whose shear
    # stiffness is the deck's own resolved A66, ``G * t = A66``, with the deck's own thickness.
    # A66 (and so the whole Bredt reference) is therefore the deck's; only E/nu are a stated
    # substitute.  This is exact for the tangential shear that carries the torsion.
    owner = _edge_owner_map(mesh, ring_global, props)
    prop_dicts: dict[tuple[int, int], dict] = {}
    edge_a66: dict[tuple[int, int], float] = {}
    for key, prop in owner.items():
        if prop is None:
            continue
        value = _membrane_shear_stiffness(prop)
        if value is None or value <= 0.0:
            continue
        if isinstance(prop, dict):
            thickness = float(prop["thickness"])
        else:
            thickness = float(prop.total_thickness)
        if thickness <= 0.0:
            continue
        shear = float(value) / thickness
        poisson = 0.3
        prop_dicts[key] = {
            "type": "isotropic",
            "e": 2.0 * shear * (1.0 + poisson),
            "nu": poisson,
            "rho": 1.0,
            "thickness": thickness,
            "shear_correction": 5.0 / 6.0,
        }
        edge_a66[key] = float(value)
    if not prop_dicts:
        print("  OBSTACLE: no deck property could be resolved for the ring edges.  Stopping.")
        return
    print(
        "  FEM wall material: isotropic substitute with G*t = the deck's own A66 per edge "
        "(orientation-independent); A66 and thickness are the deck's, E/nu are substitutes."
    )

    # ---- geometry: the ring in its own section plane, centroid at the origin
    u_ax, v_ax = section_plane_axes(SPAN)
    centroid = ring_pts.mean(axis=0)
    rel = ring_pts - centroid
    pu = rel @ u_ax
    pv = rel @ v_ax
    base = np.outer(pu, u_ax) + np.outer(pv, v_ax)  # (n_ring, 3), z = 0
    max_dim = float(max(np.ptp(pu), np.ptp(pv)))
    length = L_OVER_DIM * max_dim
    print(
        f"  section frame: u={np.round(u_ax, 4).tolist()} v={np.round(v_ax, 4).tolist()}; "
        f"dimensions: chord-extent={np.ptp(pu):.3f} m thickness-extent={np.ptp(pv):.3f} m "
        f"max_dim={max_dim:.3f} m"
    )
    print(
        f"  extruded tube: L={length:.3f} m = {L_OVER_DIM:.0f} x max_dim, nz={N_Z_TREATMENT}, "
        f"dz={length / N_Z_TREATMENT:.3f} m  (matching the rectangle's L/dim = {tube.L / 1.0:.0f})"
    )

    # ---- A66 variation (for the record), then run the generic tube case
    a66_vals = np.asarray([edge_a66[k] for k in edge_a66])
    a66_vals = a66_vals[np.isfinite(a66_vals) & (a66_vals > 0.0)]
    if a66_vals.size == 0:
        print("  OBSTACLE: no positive A66 resolved on any edge.  Stopping.")
        return
    print(
        f"  A66 per edge [N/m]: min={a66_vals.min():.4e} max={a66_vals.max():.4e} "
        f"mean={a66_vals.mean():.4e}  (variation max/min = {a66_vals.max() / a66_vals.min():.2f})"
    )
    _generic_case("airfoil", base, cells, prop_dicts, edge_a66, length, N_Z_TREATMENT)


def experiment_2b_rectangle_generic() -> None:
    """The rectangle through the SAME generic tube pipeline, to trust the treatment harness.

    The control (experiment 1) uses the fixture's own hand-built constructor; if the generic
    pipeline also reproduces the rectangle's Bredt rate then the airfoil number is a shape result,
    not a pipeline artefact.
    """
    print()
    print("-" * 96)
    print("  [harness cross-check + discretisation control] the rectangle through the SAME")
    print("  generic extruded-tube pipeline, at several ring node counts")
    print("-" * 96)
    a66 = tube.G_ISO * tube.THICKNESS
    for n_seg in (4, 8, 12):
        ring = tube._midline_ring(n_seg=n_seg)
        base = np.zeros((len(ring), 3))
        base[:, :2] = np.asarray(ring)
        n_ring = len(ring)
        cell = SectionCell(
            boundary=np.arange(n_ring, dtype=np.intp),
            area=float(tube.ENCLOSED_AREA),
            edges=tuple((i, (i + 1) % n_ring) for i in range(n_ring)),
        )
        prop = _isotropic_prop(a66, tube.THICKNESS)
        edge_a66 = {
            (min(i, (i + 1) % n_ring), max(i, (i + 1) % n_ring)): a66 for i in range(n_ring)
        }
        prop_for_edge = dict.fromkeys(edge_a66, prop)
        _generic_case(
            f"rectangle n_seg={n_seg} ({n_ring} nodes)",
            base,
            [cell],
            prop_for_edge,
            edge_a66,
            tube.L,
            tube.N_Z,
        )


def _isotropic_prop(a66: float, thickness: float) -> dict:
    """Isotropic property whose membrane shear stiffness is exactly ``a66 = G * t``."""
    shear = a66 / thickness
    poisson = 0.3
    return {
        "type": "isotropic",
        "e": 2.0 * shear * (1.0 + poisson),
        "nu": poisson,
        "rho": 1.0,
        "thickness": thickness,
        "shear_correction": 5.0 / 6.0,
    }


def main() -> int:
    np.set_printoptions(precision=6, suppress=True, linewidth=140)
    experiment_1_control()
    experiment_2_treatment()
    experiment_2b_rectangle_generic()
    return 0


if __name__ == "__main__":
    sys.exit(main())
