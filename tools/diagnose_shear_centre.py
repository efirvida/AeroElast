"""Per-station 6x6 section stiffness of the shell blade, for WU-C of issue #14.

Why this exists
---------------
Issue #14's rated twist is dominated by an eccentricity: the aerodynamic normal force acts at the
25 %-chord aerodynamic centre while the section's shear centre sits elsewhere. Whether that is a
geometry difference or a load-application artefact is WU-C, and it needs a **per-station shear
centre of the shell mesh** to compare against the BeamDyn deck's `xS`
(`K66toPropsDecoupled(K66[i], convention="BeamDyn")[8]`, metres, against the section origin).

The house extractor cannot give it: `SectionalExtractor.section_stiffness` returns four scalars
(`EA`, `GJ`, `EI_flap`, `EI_edge`) and never computes the **couplings** a shear centre is made of.
This diagnostic assembles the full 6x6 through the same Lagrange-multiplier route that method uses.

What it does and does not claim
-------------------------------
It assembles the 6x6 as the column-wise response of the section's net generalized reactions to six
**unit generalized strains** prescribed on the upper face ring (lower face clamped, interior free),
expressed in the station's own section frame. It then validates that matrix two independent ways:

1. **Diagonals against the existing scalars.** `section_stiffness` computes `EA`, `GJ`, `EI_flap`,
   `EI_edge` from the same states; this matrix's corresponding diagonal terms must reproduce them.
   That check also pins the ordering, because the four scalars are pinned by that method's own
   construction.
2. **Symmetry (Maxwell-Betti).** A stiffness matrix must be symmetric. An assembly error in the
   couplings shows up here, so this is the check that makes the couplings believable.

It deliberately does **not** call `ComputeStiffnessProps().ComputeShearCenter` yet. That function
takes `K[0:3, 0:3]` and `K[0:3, 3:6]`, solves `K1 Y = -K3` and returns `[-Y[1,2], Y[0,2]]`, and it
does not document which ordering of the six DOFs it assumes. Guessing that permutation is exactly
the failure mode this project has now hit four times (a convention left implicit), so the ordering
question is named here as open and the full matrix is printed so it can be re-mapped without
re-running. The deck comparison that settles it is WU-C's next task.

Measured status, 2026-10-07 (element_size 1.0, stations z = 5.97 / 44.16 / 77.60 / 115.21 m)
--------------------------------------------------------------------------------------------

**The diagonals are validated to machine precision.** `EA`, `GJ`, `EI_flap` and `EI_edge` from this
matrix reproduce `SectionalExtractor.section_stiffness` on the same station with a relative
difference between `0.00e+00` and `2.19e-16`. That pins the section frame, the datum, the DOF
mapping and the scaling of all six prescribed states. (An earlier revision carried a double `dz`
scaling on the rotation states and was off by exactly `dz`; the cross-check caught it, which is
what it is for.)

**The couplings are NOT trustworthy, and this scheme does not produce them.** `max|K - K^T| /
max|K|` measures `0.50`, `0.92`, `0.99` and `0.78` at those four stations, so the off-diagonal
terms fail Maxwell-Betti reciprocity by a factor of order one. Independent prescription of the six
generalized strains on a **fully clamped** lower face does not yield a section stiffness matrix:
the clamp absorbs work that is not in the section's generalized coordinates, so the mutual terms
are not the duals of the diagonal ones. **Do not compute a shear centre from this matrix.**

What follows from it: the 6x6 has to come from a formulation whose mutual terms are reciprocal by
construction (a proper influence-coefficient / Saint-Venant treatment, or the energy route), not
from six separately constrained slices. That is a design decision for WU-C's next task, recorded in
`odd/tasks/shear-centre-arbitration.md`; this diagnostic stays as the evidence that the naive route
is closed, and as the guard that the diagonals are right before anyone blames the couplings.

The deck profile decides WU-C without the 6x6 (2026-10-07)
----------------------------------------------------------

`--deck-only` prints the deck's own shear centre against the point the rated load is applied at,
both in the deck's datum (the pitch axis, the frame the anchor beam's `(x_ac - xS) * Np` term
already uses and whose sign was arbitrated on a measurement). Measured, ten stations:

| frac | chord [m] | pitch_axis | `xS` [m] | `xS/c` | `x_AC/c` | arm `/c` |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 0.000 | 5.200 | 0.505 | 0.006 | **0.0011** | 0.2545 | 0.2534 |
| 0.150 | 5.647 | 0.376 | 0.040 | **0.0071** | 0.1258 | 0.1186 |
| 0.300 | 5.367 | 0.312 | 0.077 | **0.0142** | 0.0624 | 0.0482 |
| 0.550 | 3.904 | 0.288 | 0.040 | **0.0102** | 0.0379 | 0.0276 |
| 0.850 | 2.525 | 0.323 | 0.044 | **0.0175** | 0.0731 | 0.0555 |
| 1.000 | 0.500 | 0.368 | 0.006 | **0.0115** | 0.1182 | 0.1067 |

`|xS/c|` never exceeds 0.018, so the deck's shear centre sits essentially **on its own reference
axis**, and the eccentricity the rated twist rides is the **application point** (`x_AC`, up to
0.2545 c at the root), not the shear centre. A mesh-side shear centre could differ from the deck's
by a percent or two of chord at most, which cannot move a twist that differs by a factor 2.8. So the
shear-centre hypothesis is eliminated from the deck side alone and **the 6x6 is not needed to decide
WU-C** -- it would only quantify a term that is already known to be second order.

This also reads against the docstring claim that "the mesh's measured shear centre sits at 0.477 of
the chord" (`tests/validation/blade/test_blade_rated_twist.py`): the deck puts the pitch axis between
0.288 and 0.505 of the chord and the shear centre within 0.018 c of it, and the shell's sectional
stiffness tracks the deck (S-1), so a mesh shear centre at 0.477 c is far from anything measured here
and remains unbacked until someone measures it.

Frame and datum (stated once, as the two defects this builds on demand)
----------------------------------------------------------------------
`SectionalExtractor._section_axes` measures the frame from the mesh by PCA and forces
`tangent = [0, 0, 1]` (global z), with `chord` and `flap` the transverse principal directions.
`_station_slice` refers moments to the station's **slice node centroid**. The global `(x, y)` axes
are *not* the chord/flap axes. This tool prints the three axes per station so the chord direction
(toward the leading edge or the trailing edge) is visible rather than assumed.

Diagnostic only: no assertions against tolerances, not collected by the test suite, writes nothing.

Usage:
    module load glu gcc/14.2.0_sequana
    python tools/diagnose_shear_centre.py [--element-size 1.0] [--stations 6]
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
from scipy.sparse.linalg import spsolve

# A script puts its own directory on sys.path, not the repository root, so `import tests.*` below
# fails when this is run as `python tools/diagnose_shear_centre.py`.
_REPO_ROOT = Path(__file__).resolve().parents[1]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from aeroelast.models.blade.model import Blade  # noqa: E402
from aeroelast.postprocess.sectional import (  # noqa: E402
    SectionalExtractor,
    _dofs_of_nodes,
)
from tests.support.openfast_bem import build_blade_aero_from_aerodyn  # noqa: E402
from tests.support.paths import DATA_DIR  # noqa: E402
from tests.validation.blade.test_blade_twist_anchor_beam import (  # noqa: E402
    BEAMDYN_BLADE,
    ELASTODYN_BLADE,
    _numeric_rows,
)
from tests.validation.blade.test_blade_rated_twist import AD_PRIMARY  # noqa: E402

YAML = DATA_DIR / "IEA-15-240-RWT.yaml"

#: Ordering of the six generalized strains this tool prescribes, and of the resulting generalized
#: forces that make up a column of the 6x6. Both are the same order, which is what makes the matrix
#: symmetric when the assembly is right.
ORDERING = ("axial", "shear_chord", "shear_flap", "torsion", "bend_chord", "bend_flap")


def _build_shell(element_size: float):
    model = Blade(str(YAML), element_size=element_size)
    model.generate_mesh()
    mesh = model.mesh
    if mesh is None:
        raise SystemExit("Blade.generate_mesh() produced no mesh")
    props = model.get_element_properties()
    return mesh, props


def _section_6x6(ext: SectionalExtractor, st: int, membership) -> dict:
    """Assemble one station's 6x6 by prescribed unit generalized strains.

    Each state clamps the lower face ring, prescribes a rigid motion on the upper face ring with
    translations AND rotations (unlike `section_stiffness`, which prescribes rotations only for the
    bending and torsion states), and reads the net reaction force and moment of the constrained
    slice. Columns are the generalized forces per unit generalized strain, with the same `dz`
    scaling `section_stiffness` uses so that a rotation becomes a curvature and a translation a
    shear strain.
    """
    sl = ext._station_slice(st, membership)
    lower, upper, interior = sl["lower"], sl["upper"], sl["interior"]
    if upper.size == 0 or lower.size == 0:
        raise ValueError(f"station {st}: empty face ring")
    ref = sl["centroid"]
    tangent, chord, flap = ext._section_axes(membership[st]["nodes"])
    dz = max(sl["z1"] - sl["z0"], 1e-9)

    all_nodes = np.concatenate([lower, interior, upper])
    dofs = _dofs_of_nodes(all_nodes, ext.dpn)
    K_loc = ext._submatrix(dofs).tocsc()
    n_all = all_nodes.size

    eps = 1e-6  # prescribed generalized strain (dimensionless for shear, rad for a face rotation)

    def reactions(upper_trans: np.ndarray, upper_rots: np.ndarray | None):
        pres_dofs = list(np.concatenate([lower * ext.dpn + d for d in range(6)]))
        pres_vals = [0.0] * (lower.size * 6)
        for d in range(3):
            pres_dofs.extend((upper * ext.dpn + d).tolist())
        pres_vals.extend(np.asarray(upper_trans, dtype=np.float64).ravel().tolist())
        if upper_rots is not None:
            for d in range(3, 6):
                pres_dofs.extend((upper * ext.dpn + d).tolist())
            pres_vals.extend(np.asarray(upper_rots, dtype=np.float64).ravel().tolist())
        pres_dofs = np.asarray(pres_dofs, dtype=np.int64)
        pres_vals = np.asarray(pres_vals, dtype=np.float64)
        free = np.setdiff1d(dofs, pres_dofs)

        local = {int(d): k for k, d in enumerate(dofs)}
        u_loc = np.zeros(dofs.size)
        u_loc[[local[int(d)] for d in pres_dofs]] = pres_vals
        if free.size:
            f_idx = np.array([local[int(d)] for d in free])
            p_idx = np.array([local[int(d)] for d in pres_dofs])
            u_p = np.array([u_loc[local[int(d)]] for d in pres_dofs])
            u_loc[f_idx] = np.asarray(
                spsolve(K_loc[f_idx][:, f_idx].tocsc(), -(K_loc[f_idx][:, p_idx].tocsc() @ u_p)),
                dtype=np.float64,
            ).ravel()
        r = (K_loc @ u_loc).reshape(n_all, ext.dpn)
        f_glob = r[:, :3].sum(axis=0)
        m_glob = np.cross(ext.coords[all_nodes] - ref, r[:, :3]).sum(axis=0)
        return f_glob, m_glob

    def generalized(f_glob, m_glob):
        return np.array(
            [
                float(f_glob @ tangent),
                float(f_glob @ chord),
                float(f_glob @ flap),
                float(m_glob @ tangent),
                float(m_glob @ chord),
                float(m_glob @ flap),
            ]
        )

    def rigid_rot(axis, theta):
        return theta * np.cross(axis, ext.coords[upper] - ref)

    def pure_bend(axis, theta):
        d_xy = (ext.coords[upper] - ref).copy()
        d_xy[:, 2] = 0.0
        vals = np.zeros((upper.size, 3))
        vals[:, 2] = theta * (axis[0] * d_xy[:, 1] - axis[1] * d_xy[:, 0])
        return vals

    # Six unit generalized strains. Transverse translations are `eps * dz` along the corresponding
    # section axis so that the engineering shear strain is `eps`; face rotations are `eps`, so the
    # curvature is `eps / dz` and the torsion/bending columns are scaled by `dz` below.
    states = {
        "axial": (np.outer(np.full(upper.size, eps * dz), tangent), None, eps),
        "shear_chord": (np.outer(np.full(upper.size, eps * dz), chord), None, eps),
        "shear_flap": (np.outer(np.full(upper.size, eps * dz), flap), None, eps),
        "torsion": (rigid_rot(tangent, eps), ext._projected_rots(upper, tangent, eps), eps / dz),
        "bend_chord": (pure_bend(chord, eps), ext._projected_rots(upper, chord, eps), eps / dz),
        "bend_flap": (pure_bend(flap, eps), ext._projected_rots(upper, flap, eps), eps / dz),
    }

    K = np.zeros((6, 6))
    for k, name in enumerate(ORDERING):
        trans, rots, strain = states[name]
        f_glob, m_glob = reactions(np.asarray(trans), rots)
        gen = generalized(f_glob, m_glob)
        # `strain` is already the curvature (eps / dz) for the rotation states, so the column is
        # the generalized force per unit generalized strain with no further scaling. Multiplying
        # by dz here as well double-counted the station length and inflated GJ/EI by exactly dz.
        K[:, k] = gen / strain

    # a pure bending state must produce no net force; a pure axial state no net moment about the
    # tangent. Printed, not asserted: a diagnostic reports.
    return {
        "K": K,
        "dz": dz,
        "ref": ref,
        "axes": (tangent, chord, flap),
        "chord_len": float(np.ptp(ext.coords[membership[st]["nodes"]] @ chord)),
    }


def _deck_profile(n_show: int) -> None:
    """The deck's own shear centre against the point the rated load is applied at.

    Everything here is expressed in the deck's datum, the **pitch axis**, which is the frame the
    anchor beam's `m = Mp + (x_ac - xS) * Np - ...` already uses and whose sign was arbitrated on a
    measurement (the opposite one gave a ~10x torque error). Comparing `xS` against the aerodynamic
    centre in that same frame needs no new convention decision:

        x_AC = (pitch_axis - 0.25) * chord          the AC's offset from the pitch axis
        xS                                          the shear centre's offset from the pitch axis
        x_AC - xS                                   the eccentricity arm the beam integrates

    If `xS` is a small fraction of the chord, the eccentricity that dominates the rated twist is the
    force acting at 0.25 c from the leading edge (a definition/application matter). If `xS` is a
    large fraction, the section's shear centre genuinely sits far from the aerodynamic centre
    (a geometry matter). That is WU-C's question, answered from the deck alone -- no mesh 6x6 needed.
    """
    if not BEAMDYN_BLADE.exists() or not ELASTODYN_BLADE.exists():
        print("deck not present under tests/IEA15MW/reference; deck profile skipped")
        return
    from openfast_toolbox.converters.beam import K66toPropsDecoupled

    aero = build_blade_aero_from_aerodyn(AD_PRIMARY)
    r_hub = np.asarray(aero.r, dtype=float)
    chord = np.asarray(aero.chord, dtype=float)
    hub = float(aero.hub_radius)
    blade_len = float(aero.blade_length)

    rows = _numeric_rows(BEAMDYN_BLADE, "DISTRIBUTED PROPERTIES")
    n_deck = len(rows) // 13
    frac = np.array([rows[i * 13][0] for i in range(n_deck)])
    K66 = np.array([np.asarray(rows[i * 13 + 1 : i * 13 + 7]) for i in range(n_deck)])
    deck = [K66toPropsDecoupled(K66[i], convention="BeamDyn") for i in range(n_deck)]
    x_shear = np.array([p[8] for p in deck])  # metres, offset from the pitch axis
    x_cent = np.array([p[6] for p in deck])

    ed = np.asarray(_numeric_rows(ELASTODYN_BLADE, "DISTRIBUTED BLADE PROPERTIES"))
    pitch_axis = np.interp(frac, ed[:, 0], ed[:, 1])  # chord fraction from the leading edge
    chord_d = np.interp(hub + frac * blade_len, r_hub, chord)
    x_ac = (pitch_axis - 0.25) * chord_d  # the anchor beam's own arbitrated expression

    print()
    print("deck shear centre vs the aerodynamic centre, both in the deck's datum (the pitch axis):")
    print("  (metres; x_AC is the anchor beam's own (pitch_axis - 0.25) * chord, whose sign was")
    print("   arbitrated on a measurement, so the two are directly comparable)")
    print(
        f"  {'frac':>6} {'chord[m]':>9} {'pitch_axis':>10} {'xS[m]':>9} {'xS/c':>7} "
        f"{'xC[m]':>8} {'x_AC[m]':>9} {'x_AC/c':>7} {'x_AC-xS[m]':>11} {'arm/c':>7}"
    )
    idx = np.linspace(0, n_deck - 1, n_show).astype(int)
    for i in idx:
        c = chord_d[i]
        print(
            f"  {frac[i]:6.3f} {c:9.3f} {pitch_axis[i]:10.3f} {x_shear[i]:9.4f} "
            f"{x_shear[i] / c:7.4f} {x_cent[i]:8.4f} {x_ac[i]:9.4f} {x_ac[i] / c:7.4f} "
            f"{x_ac[i] - x_shear[i]:11.4f} {(x_ac[i] - x_shear[i]) / c:7.4f}"
        )
    print()
    print("  reading: |xS/c| small means the eccentricity is the 0.25 c application point at work,")
    print("  not the section's shear centre; |xS/c| of order 0.2 means the deck puts the shear")
    print("  centre a fifth of a chord from its own reference axis, which is a geometry fact.")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--element-size", type=float, default=1.0)
    ap.add_argument("--stations", type=int, default=6)
    ap.add_argument(
        "--deck-only",
        action="store_true",
        help="print only the deck profile (no mesh 6x6, which needs a second run)",
    )
    args = ap.parse_args()

    _deck_profile(args.stations)
    if args.deck_only:
        return 0

    mesh, props = _build_shell(args.element_size)
    ext = SectionalExtractor(mesh, props)
    membership = ext.station_membership()
    stations = sorted(membership)
    print(
        f"mesh: element_size={args.element_size}, {mesh.coords_array.shape[0]} nodes, "
        f"{len(stations)} stations"
    )

    picks = [
        stations[int(round(f * (len(stations) - 1)))]
        for f in np.linspace(0.05, 0.95, args.stations)
    ]
    picks = sorted(set(picks))

    print()
    print("ordering of both the generalized strains and the generalized forces:")
    print("   " + ", ".join(f"{i}={n}" for i, n in enumerate(ORDERING)))
    print("   generalized forces are (F.tangent, F.chord, F.flap, M.tangent, M.chord, M.flap)")
    print()

    for st in picks:
        try:
            out = _section_6x6(ext, st, membership)
            sec = ext.section_stiffness(st, membership)
        except (ValueError, RuntimeError) as exc:
            print(f"station {st} (z={membership[st]['z']:.2f}): skipped ({exc})")
            continue
        K = out["K"]
        asymmetry = float(np.max(np.abs(K - K.T)) / max(np.max(np.abs(K)), 1e-30))
        print(
            f"station {st}  z={membership[st]['z']:.2f} m  dz={out['dz']:.3f} m  "
            f"chord={out['chord_len']:.3f} m  ref={np.round(out['ref'], 3)}"
        )
        t, c, f = out["axes"]
        print(f"   axes: tangent={np.round(t, 4)} chord={np.round(c, 4)} flap={np.round(f, 4)}")
        print(
            f"   max|K - K^T| / max|K| = {asymmetry:.3e}   "
            f"(Maxwell-Betti: must be ~0 for a correct assembly)"
        )
        print("   K =")
        for row in K:
            print("     " + "  ".join(f"{v: .4e}" for v in row))
        checks = (
            ("EA", K[0, 0], sec["EA"]),
            ("GJ", K[3, 3], sec["GJ"]),
            ("EI_flap", K[4, 4], sec["EI_flap"]),
            ("EI_edge", K[5, 5], sec["EI_edge"]),
        )
        print("   diagonals vs section_stiffness (the ordered cross-check):")
        for name, mine, theirs in checks:
            rel = abs(abs(mine) - theirs) / theirs if theirs else float("nan")
            print(
                f"     {name:8} this tool {abs(mine): .6e}   section_stiffness {theirs: .6e}"
                f"   rel.diff {rel:.2e}"
            )
        print()

    return 0


if __name__ == "__main__":
    sys.exit(main())
