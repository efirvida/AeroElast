"""Diagnostic: does the multi-cell realisation converge to the minimum-norm one?

Issue #11.  The multi-cell Bredt-Batho realisation is now wired into
``ForceProjector`` (see ``force_projection.py``).  On the IEA-15MW blade at
``element_size=1.0`` (14-16 nodes per physical ring, i.e. 1-2 per wall) the wired
production path measures a **nose-up** tip section rotation at the rated point,
``omega = +0.1729 deg``, while the pre-wiring minimum-norm realisation carried the
physical nose-down sense.  Refining the mesh to ``element_size=0.5`` (≈28 nodes per
ring, 3-4 per wall) flips the sign back to nose-down with the *same* wiring, so the
coarse-mesh sign is a discretisation artefact rather than a formulation defect.
This tool measures that directly.

For each requested ``element_size`` it builds the IEA-15MW blade, projects the
**same** rated BEM load through both realisations, solves the shell once per
realisation, and prints:

* the mesh size (nodes, elements) and the nodes per physical ring at the mid and
  tip bands;
* the tip section rotation ``omega`` [deg] and ``distortion/|omega|`` from the
  suite's own estimator
  (:func:`tests.validation.blade.test_blade_rated_twist._ring_kinematics`), for
  **both** the minimum-norm realisation (``ForceProjector._distribute`` for every
  strip - what ``main`` does, because every strip holds several physical rings and
  the old single-ring path fell back to it) and the wired multi-cell one
  (``ForceProjector.project``);
* the wall-clock seconds and the process peak RSS for each solve.

The multi-cell realisation is fed the blade's laminate map
(``Blade.get_element_properties()``), so each wall's Bredt-Batho flow sees its own
membrane shear stiffness ``S = G*t``: this is the ``standalone.py`` production
wiring.  The FSI participant and the validation-guard fixture do **not** hold that
map and fall back to the geometric-only ``S = 1.0``; pass ``--uniform-stiffness``
to reproduce that variant.  Both variants are measured here to converge to the
same nose-down sign once the mesh is refined.

**This is a diagnostic, not a guard, and it writes nothing.**  It asserts no
tolerance, changes no repository file, and is not collected by the test suite; the
guards live in ``tests/validation/blade/``.  The rated point comes from the same
helpers the validation module uses (``V_RATED``, ``RPM_RATED``, ``PITCH_RATED``,
``build_blade_aero_from_aerodyn``, ``BEMSolver``), so the load is the one the guards
use.

The default element-size list is ``1.0, 0.5, 0.25``.  ``0.25`` is *attempted* and
skipped gracefully when a stated limit would be exceeded: the per-case wall-clock
limit (``--max-seconds``, default 1200 s) and the peak-RSS limit (``--max-rss-gb``,
default 96 GB).  Because the direct solve's memory and time both grow at least with
the node count, the tool projects the next case from the last measured one and
skips only when even that lower-bound projection is already over a limit; a
``MemoryError`` during a case is caught and reported the same way.

Usage:
    python tools/diagnose_multicell_convergence.py
    python tools/diagnose_multicell_convergence.py --element-sizes 1.0 0.5
    python tools/diagnose_multicell_convergence.py --max-seconds 600 --max-rss-gb 64
"""

from __future__ import annotations

import argparse
import resource
import sys
import time
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from _aeroelast import PyMeshAssembler  # noqa: E402
from aeroelast.core.mesh.entities import MeshElement, Node  # noqa: E402
from aeroelast.models.blade.model import Blade  # noqa: E402
from aeroelast.solvers.bem.engine import BEMSolver  # noqa: E402
from aeroelast.solvers.bem.force_projection import ForceProjector  # noqa: E402

import tests.validation.blade.test_blade_iea15mw_validation as blade_validation  # noqa: E402
from tests.support.openfast_bem import build_blade_aero_from_aerodyn  # noqa: E402
from tests.validation.blade.test_blade_rated_twist import (  # noqa: E402
    AD_PRIMARY,
    PITCH_RATED,
    RPM_RATED,
    STATION_GAP_TOLERANCE,
    V_RATED,
    YAML,
    ZHOU_TIP_TORSION_DEG,
    _physical_stations,
    _ring_kinematics,
)

#: Effective default element sizes.  ``0.25`` is attempted last and skipped
#: gracefully when the stated time / memory limits would be exceeded.
DEFAULT_ELEMENT_SIZES = (1.0, 0.5, 0.25)


def _peak_rss_gb() -> float:
    """Process peak resident set size in GB (Linux ``ru_maxrss`` is in KiB)."""
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / (1024.0**2)


def _build_blade(element_size: float):
    """Build the blade mesh, its sparse stiffness and the physical-ring map.

    Returns everything the two realisations need: the Rust assembler, the free
    DOFs (root clamped, as in the validation fixtures), the snapshot stiffness
    ``Kff``, the node coordinates, the merged physical stations / rings, and the
    tip ring.
    """
    Node._id_counter = 0
    MeshElement._id_counter = 0
    model = Blade(str(YAML), element_size=element_size)
    model.generate_mesh()
    mesh = model.mesh
    props = model.get_element_properties()
    assembler = PyMeshAssembler.from_model(
        blade_validation._to_rust_mesh(mesh, props),
        props,
        list(blade_validation.SPAN_DIRECTION),
        None,
    )
    n = assembler.dofs_count
    rows, cols, vals = assembler.assemble_k()
    K = coo_matrix((np.asarray(vals), (np.asarray(rows), np.asarray(cols))), shape=(n, n)).tocsr()
    root = {mesh.node_id_to_index[nid] for nid in mesh.get_node_set("RootNodes").node_ids}
    free = np.array(
        [i for i in range(n) if i not in {6 * r + d for r in root for d in range(6)}],
        dtype=np.int64,
    )
    coords = np.array([[nd.x, nd.y, nd.z] for nd in mesh.nodes], dtype=float)

    # The blade is prebent, so the raw unique-z set splits each physical ring into
    # near-duplicate buckets; merge them exactly as the validation module does, so
    # the estimator runs on complete rings.
    phys_stations = _physical_stations(coords)
    phys_rings = {
        zz: np.where(np.abs(coords[:, 2] - zz) < STATION_GAP_TOLERANCE)[0] for zz in phys_stations
    }
    tip_z = phys_stations[-1]
    mid_z = phys_stations[
        int(np.argmin(np.abs(np.asarray(phys_stations, dtype=float) - 0.5 * tip_z)))
    ]
    return {
        "mesh": mesh,
        "n": n,
        "Kff": K[np.ix_(free, free)],
        "free": free,
        "coords": coords,
        "mid_count": len(phys_rings[mid_z]),
        "tip_count": len(phys_rings[tip_z]),
        "tip": phys_rings[tip_z],
        "nodes": len(mesh.nodes),
        "elements": len(mesh.elements),
        "properties": props,
    }


def _minimum_norm_forces(projector: ForceProjector, bem) -> np.ndarray:
    """The pre-wiring realisation: ``_distribute`` for every strip.

    Mirrors :meth:`ForceProjector.project`'s per-strip load construction (same
    frames, same AC-to-centroid arm, same ``Mp``) and routes each strip through the
    minimum-norm solve directly.  On this blade every strip holds 2-5 physical
    rings, so the old ``realise_section_load(..., ring_groups=...)`` call fell back
    to ``_distribute`` for every strip; this is that load, with no wall flow.
    """
    forces = np.zeros((projector._n_nodes, 3))
    for k, strip in enumerate(projector._strips):
        if len(strip.node_indices) == 0:
            continue
        f_strip = (
            float(bem.Np[k]) * strip.dr * projector._strip_normal_dirs[k]
            + float(bem.Tp[k]) * strip.dr * projector._strip_chord_dirs[k]
        )
        m_ac = (
            float(bem.Mp[k]) * strip.dr * projector._span_dir if bem.Mp is not None else np.zeros(3)
        )
        m_strip = m_ac - np.cross(projector._strip_ac_offsets[k], f_strip)
        forces[strip.node_indices] = ForceProjector._distribute(strip, f_strip, m_strip)
    return forces


def _solve_and_measure(blade, force, label: str):
    """Solve ``K u = f`` (root clamped) and measure the tip section kinematics.

    ``force`` is the ``(n_nodes, 3)`` nodal field both realisations return; it is
    scattered into the assembler's 6-DOF layout (``[u, v, w, tx, ty, tz]``) exactly
    as the validation fixtures do.  Returns ``(omega_deg, distortion_ratio,
    seconds, peak_rss_gb)``.
    """
    f = np.zeros(blade["n"])
    f[0::6] = force[:, 0]
    f[1::6] = force[:, 1]
    f[2::6] = force[:, 2]

    start = time.perf_counter()
    u = np.zeros(blade["n"])
    u[blade["free"]] = spsolve(blade["Kff"], f[blade["free"]])
    seconds = time.perf_counter() - start
    kin = _ring_kinematics(blade["coords"], u, blade["tip"])
    omega = float(np.rad2deg(kin["omega"]))
    ratio = kin["distortion"] / abs(kin["omega"]) if abs(kin["omega"]) > 1e-30 else float("inf")
    print(
        f"    {label:10s} omega={omega:+9.4f} deg   distortion/|omega|={ratio:8.4f}   "
        f"solve={seconds:8.2f} s   peakRSS={_peak_rss_gb():7.2f} GB",
        flush=True,
    )
    return omega, ratio, seconds


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--element-sizes",
        type=float,
        nargs="+",
        default=list(DEFAULT_ELEMENT_SIZES),
        help="chordwise/spanwise mesh sizes to sweep (default: 1.0 0.5 0.25)",
    )
    parser.add_argument(
        "--max-seconds",
        type=float,
        default=1200.0,
        help="per-case wall-clock limit in seconds; a projected case above it is skipped",
    )
    parser.add_argument(
        "--max-rss-gb",
        type=float,
        default=96.0,
        help="peak-RSS limit in GB; a projected case above it is skipped",
    )
    parser.add_argument(
        "--uniform-stiffness",
        action="store_true",
        help=(
            "multi-cell walls use the uniform S=1.0 geometric-only fallback "
            "(the FSI-participant / validation-guard wiring) instead of the "
            "laminate S=G*t map the standalone production path passes"
        ),
    )
    args = parser.parse_args()

    if not AD_PRIMARY.exists():
        raise SystemExit(f"AeroDyn reference not present: {AD_PRIMARY}")

    print("multi-cell vs minimum-norm rated-tip convergence (diagnostic; writes nothing)")
    print(
        f"rated point: V={V_RATED} m/s, {RPM_RATED} rpm, pitch {PITCH_RATED} deg; "
        f"Zhou tip torsion {ZHOU_TIP_TORSION_DEG:+.2f} deg"
    )
    print(f"mesh sizes: {args.element_sizes}")
    print(
        f"limits: per-case {args.max_seconds:.0f} s, peak RSS {args.max_rss_gb:.0f} GB "
        f"(the last size is attempted and skipped if a projected case exceeds them)"
    )
    print(
        "multi-cell stiffness: "
        + (
            "uniform S=1.0 (participant/guard fixture)"
            if args.uniform_stiffness
            else "laminate S=G*t (standalone production path)"
        )
    )
    print()

    blade_aero = build_blade_aero_from_aerodyn(AD_PRIMARY)
    bem = BEMSolver(blade_aero, rho=1.225, mu=1.81206e-5, hub_height=150.0, shear_exp=0.0).compute(
        V_RATED, RPM_RATED, PITCH_RATED
    )

    reference = None  # (nodes, peak_rss_gb, seconds) of the previous completed case
    size_prev: float | None = None  # element_size of the previous completed case
    summary: list[tuple[float, float, float]] = []  # (size, omega_min, omega_multi)
    for size in args.element_sizes:
        if reference is not None and size_prev is not None:
            nodes_prev, rss_prev, seconds_prev = reference
            growth = (size_prev / size) ** 2
            nodes_est = nodes_prev * growth
            # Lower-bound projections: assembly/solve scale at least with the DOF
            # count, so skipping only when the projection is already over a limit
            # never discards a case that could have fit.
            rss_est = rss_prev * growth
            seconds_est = seconds_prev * growth
            if rss_est > args.max_rss_gb or seconds_est > args.max_seconds:
                print(
                    f"[size={size:.2f}] skipped: projected from size={size_prev:.2f} "
                    f"(~{nodes_est:.0f} nodes, ~{rss_est:.1f} GB, ~{seconds_est:.0f} s) "
                    f"exceeds the stated limits ({args.max_rss_gb:.0f} GB / "
                    f"{args.max_seconds:.0f} s)"
                )
                continue

        print(f"[size={size:.2f}] building blade and assembling K ...", flush=True)
        case_start = time.perf_counter()
        try:
            blade = _build_blade(size)
        except MemoryError:
            print(f"[size={size:.2f}] skipped: MemoryError while building the mesh")
            continue
        print(
            f"[size={size:.2f}] nodes={blade['nodes']} elements={blade['elements']} | "
            f"mid ring={blade['mid_count']} nodes | tip ring={blade['tip_count']} nodes",
            flush=True,
        )

        try:
            projector = ForceProjector(
                blade["mesh"],
                blade_aero,
                span_direction=blade_validation.SPAN_DIRECTION,
                element_properties=None if args.uniform_stiffness else blade["properties"],
            )
            multi_cell = projector.project(bem)
            minimum_norm = _minimum_norm_forces(projector, bem)
            omega_min, ratio_min, sec_min = _solve_and_measure(blade, minimum_norm, "min-norm")
            omega_multi, ratio_multi, sec_multi = _solve_and_measure(
                blade, multi_cell, "multi-cell"
            )
        except MemoryError:
            print(f"[size={size:.2f}] skipped: MemoryError during projection/solve")
            continue

        geometry_build = time.perf_counter() - case_start
        print(f"[size={size:.2f}] case wall-clock {geometry_build:.2f} s", flush=True)
        summary.append((size, omega_min, omega_multi))
        reference = (blade["nodes"], _peak_rss_gb(), max(sec_min, sec_multi))
        size_prev = size

    print()
    print("summary: tip omega [deg] (min-norm / multi-cell) as the mesh refines")
    for size, omega_min, omega_multi in summary:
        same = "same sign" if np.sign(omega_min) == np.sign(omega_multi) else "OPPOSITE signs"
        gap = abs(omega_min - omega_multi)
        print(
            f"  element_size={size:.2f}: min-norm {omega_min:+9.4f}  "
            f"multi-cell {omega_multi:+9.4f}  seed gap {gap:.4f} deg  ({same})"
        )
    if len(summary) >= 2:
        first, last = summary[0], summary[-1]
        gap_first = abs(first[1] - first[2])
        gap_last = abs(last[1] - last[2])
        same_last = np.sign(last[1]) == np.sign(last[2])
        sign_note = (
            "share a sign at the finest mesh"
            if same_last
            else "still differ in sign at the finest mesh"
        )
        trend = "shrinks" if gap_last < gap_first else "does not shrink"
        print(
            f"  from size={first[0]:.2f} to size={last[0]:.2f}: the two realisations "
            f"{sign_note}, the gap {trend} from {gap_first:.4f} to {gap_last:.4f} deg; "
            f"min-norm {first[1]:+.4f} -> {last[1]:+.4f}, "
            f"multi-cell {first[2]:+.4f} -> {last[2]:+.4f}"
        )


if __name__ == "__main__":
    main()
