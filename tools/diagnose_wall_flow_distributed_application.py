"""#30 T5: does the **production application mode** of a section moment deliver Bredt?

The validated fixture (group 31, ``test_thin_walled_tube_moment_realization.py``) applies the wall
shear flow as an **end couple**: ``+T`` on the tip ring and ``-T`` on the root ring, so the internal
torque is constant along the tube and the interior window is far from any application point.
Production does something else.  ``ForceProjector.project`` hands ``realise_section_load`` one
**strip** (a BEM band, ~13 physical rings on the campaign mesh) and
``_realise_multi_cell_section_load`` splits the strip's torsion over its rings with
``_ring_tributary_weights`` (trapezoidal), each ring realising its own share as its own wall flow -
so the internal torque **ramps** and a self-equilibrated couple is applied at *every* ring.  On the
campaign mesh every ring is single-cell, so the per-ring flow is the validated single-ring Bredt
``q = T_i/(2A)``.

T5 asks whether that application mode is Bredt-consistent for the internal torque it actually
carries, because the answer decides the reading of #30's coupled movement: Bredt-consistent means the
outer-span amplification T3b part 2 measured is the deck section's own physics and the activation is
justified; over-delivery means the production *application* carries a defect that the fixture with
the reference can name.

Metric, insensitive to the torque ramp::

    r(z) = theta'(z) * GJ_Bredt / M_cum(z)

``theta'(z)`` is a local least-squares slope of the fixture's own ``_theta_fit`` profile (exact for a
quadratic profile, so the ramp itself does not bias it), and ``M_cum(z)`` is the **outboard
cumulative applied moment read back from the applied field itself**
(``z_hat . sum_{outboard nodes} (x_j - c_z) x f_j``): no assumption about how the load was spread.
``r = 1`` is Bredt and the suite's 5 % rule is the band.

Modes, root clamped (the case's own boundary condition):

* **end** - the validated application, ``+T`` on the tip ring only; its ``M_cum`` is constant and
  ``r`` reproduces the record's clamped ``rate/Bredt = 0.99156``;
* **distributed wall flow** - one share per ring, ``T * _ring_tributary_weights(ring positions)``,
  each ring through the production single-ring route; this is production's structure;
* **distributed minimum-norm** - the parked default, applied as its branch applies it: **one**
  ``ForceProjector._distribute`` over the whole node set for the whole moment.

``--sweep`` refines ``n_z`` to separate a genuine application defect from a discretisation effect of
applying a discrete couple per ring: a Saint-Venant disturbance decays over ~one section dimension
(here ~1 m), so at ``dz = 0.2`` m the per-ring disturbances overlap along the whole tube.

Diagnostic only: no assertions, not collected by the suite, writes nothing.  Run with the repository
bootstrap::

    scripts/aeroenv.sh python tools/diagnose_wall_flow_distributed_application.py
    scripts/aeroenv.sh python tools/diagnose_wall_flow_distributed_application.py --sweep 30,60,120
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

_REPO_ROOT = Path(__file__).resolve().parents[1]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from aeroelast.solvers.bem.force_projection import (  # noqa: E402
    ForceProjector,
    _ring_tributary_weights,
    _Strip,
)

import tests.validation.parity.test_thin_walled_tube_moment_realization as mr  # noqa: E402
import tests.validation.parity.test_thin_walled_tube_torsion as tube  # noqa: E402

SPAN = np.array([0.0, 0.0, 1.0])
LOCAL_HALF_WIDTH = 2  # rings each side of a station for the local slope
RAMP_FLOOR = 0.25  # stations whose M_cum is below this share of T are reported separately


def _embed(nodes_local: np.ndarray, ring, n_nodes: int) -> np.ndarray:
    f = np.zeros(6 * n_nodes)
    ring = np.asarray(ring, dtype=np.intp)
    for j, node in enumerate(ring):
        f[6 * node : 6 * node + 3] = nodes_local[j]
    return f


def _ring_field(coords: np.ndarray, ring, torque: float) -> np.ndarray:
    """The validated single-ring production route, embedded in the mesh's 6-DOF layout."""
    nodes = mr._production_tip_moment(coords, ring, torque).reshape(-1, 6)
    return _embed(nodes[np.asarray(ring, dtype=np.intp), :3], ring, len(coords))


def _distributed_wall_flow(coords: np.ndarray, rings, z: np.ndarray) -> np.ndarray:
    """Production's structure: one tributary share of the total torsion per ring, its own wall flow."""
    weights = _ring_tributary_weights(z)
    f = np.zeros(6 * len(coords))
    for ring, weight in zip(rings, weights, strict=True):
        f += _ring_field(coords, ring, tube.TORQUE * float(weight))
    return f


def _distributed_min_norm(coords: np.ndarray, dr: float) -> np.ndarray:
    """The parked default's branch: one minimum-norm solve for the whole moment over all nodes."""
    centroid = coords.mean(axis=0)
    strip = _Strip(
        node_indices=np.arange(len(coords), dtype=np.intp),
        r_center=float(centroid[2]),
        dr=dr,
        centroid=centroid,
        offsets=coords - centroid,
    )
    nodes = ForceProjector._distribute(strip, np.zeros(3), np.array([0.0, 0.0, tube.TORQUE]))
    f = np.zeros(6 * len(coords))
    f[0::6] = nodes[:, 0]
    f[1::6] = nodes[:, 1]
    f[2::6] = nodes[:, 2]
    return f


def _cumulative_moment(coords, f, ring, i: int) -> float:
    """``z_hat . sum_{nodes outboard of station i} (x_j - c_i) x f_j`` - read from the field itself."""
    centre = coords[ring].mean(axis=0)
    outboard = coords[:, 2] > float(coords[ring][:, 2].mean())
    force = f.reshape(-1, 6)[:, :3]
    return float(np.cross(coords[outboard] - centre, force[outboard]).sum(axis=0) @ SPAN)


def _work(f: np.ndarray, u: np.ndarray) -> float:
    """Work done by the applied field on the clamped structure (the clamped root does no work)."""
    return 0.5 * float(f @ u)


def _bredt_work(coords, f, rings, z: np.ndarray, gj: float) -> float:
    """Bredt's strain energy for the internal-torque profile the applied field itself produces.

    ``integral M_cum(z)^2 / (2 GJ) dz`` with ``M_cum`` read back from the field.  This is a
    global, estimator-free arbiter: it does not use any local slope, so it can confirm or
    refute the local ``r(z)`` metric instead of repeating it.
    """
    m_cum = np.asarray([_cumulative_moment(coords, f, rings[i], i) for i in range(len(rings))])
    return float(np.trapezoid(m_cum**2, z)) / (2.0 * gj)


def _local_slope(z: np.ndarray, theta: np.ndarray, i: int) -> float:
    lo = max(0, i - LOCAL_HALF_WIDTH)
    hi = min(len(z), i + LOCAL_HALF_WIDTH + 1)
    return float(np.polyfit(z[lo:hi], theta[lo:hi], 1)[0])


def _mode_ratios(coords, u, rings, z, f, gj, stations) -> tuple[np.ndarray, np.ndarray]:
    theta = np.asarray([tube._theta_fit(coords, u, ring) for ring in rings])
    ratios = []
    for i in stations:
        m_cum = _cumulative_moment(coords, f, rings[i], i)
        rate = _local_slope(z, theta, i)
        ratios.append(rate * gj / m_cum)
    return np.asarray(ratios), np.asarray(
        [_cumulative_moment(coords, f, rings[i], i) for i in stations]
    )


def run_mesh(n_z: int, verbose: bool) -> dict:
    coords, conn, rings, _n_ring = tube._tube_mesh(n_z=n_z)
    dr = tube.L / n_z
    _, K = tube._assemble(coords, conn, 4, tube._iso_prop())
    gj, ref_rate = tube._bredt_isotropic()
    clamped = tube._clamped_dofs(rings)
    z = np.asarray([coords[ring][:, 2].mean() for ring in rings])
    interior = [
        i
        for i in range(len(rings))
        if 0.4 * tube.L <= z[i] <= 0.9 * tube.L
        and LOCAL_HALF_WIDTH <= i < len(rings) - LOCAL_HALF_WIDTH
    ]

    # harness control: the validated end application, self-equilibrated (record 1.00309)
    f_end_eq = _ring_field(coords, rings[-1], +tube.TORQUE) + _ring_field(
        coords, rings[0], -tube.TORQUE
    )
    u = tube._solve_rigid_removed(K, f_end_eq, tube._rigid_modes(coords))
    theta = np.asarray([tube._theta_fit(coords, u, ring) for ring in rings])
    control = tube._lsq_slope(z, theta, 0.4 * tube.L, 0.9 * tube.L) / ref_rate

    out = {"n_z": n_z, "dz": dr, "control": control, "rings": len(rings), "modes": {}}
    cases = (
        ("end", _ring_field(coords, rings[-1], +tube.TORQUE)),
        ("wall flow distributed", _distributed_wall_flow(coords, rings, z)),
        ("min-norm distributed", _distributed_min_norm(coords, dr)),
    )
    for label, f in cases:
        realized = tube._realised_torque(coords, f)
        u = tube._solve_clamped(K, f, clamped)
        ratios, m_cum = _mode_ratios(coords, u, rings, z, f, gj, interior)
        loaded = m_cum > RAMP_FLOOR * tube.TORQUE
        energy = _work(f, u) / _bredt_work(coords, f, rings, z, gj)
        out["modes"][label] = {
            "realized": realized,
            "energy_ratio": energy,
            "median": float(np.median(ratios)),
            "worst": float(np.max(np.abs(ratios - 1.0))),
            "median_loaded": float(np.median(ratios[loaded])) if loaded.any() else float("nan"),
            "worst_loaded": float(np.max(np.abs(ratios[loaded] - 1.0))) if loaded.any() else float("nan"),
        }
        if verbose:
            print(f"\n[{label}] realised torque {realized:.6e} N.m of {tube.TORQUE:.1e}")
            print(f"    {'z [m]':>7} {'M_cum [N.m]':>13} {'r':>10}")
            for i in interior[:: max(1, len(interior) // 8)]:
                m_cum_i = _cumulative_moment(coords, f, rings[i], i)
                theta_i = np.asarray([tube._theta_fit(coords, u, ring) for ring in rings])
                print(
                    f"    {z[i]:>7.3f} {m_cum_i:>13.4e} "
                    f"{_local_slope(z, theta_i, i) * gj / m_cum_i:>10.5f}"
                )
    return out


def _verdict(metrics: dict) -> str:
    return "INSIDE" if metrics["worst"] <= tube.TOL else "OUTSIDE"


def main() -> int:
    parser = argparse.ArgumentParser(
        description="#30 T5: the production distributed application mode against Bredt"
    )
    parser.add_argument("--n-z", type=int, default=tube.N_Z, help="rings along the tube")
    parser.add_argument(
        "--sweep",
        default=None,
        help="comma-separated n_z values; prints one comparison line per mesh",
    )
    args = parser.parse_args()
    np.set_printoptions(precision=6, suppress=True, linewidth=160)

    if args.sweep:
        values = [int(v) for v in args.sweep.split(",")]
        runs = [run_mesh(n_z, verbose=False) for n_z in values]
        print(
            "refinement sweep, interior window (0.4L, 0.9L).  W/Wbredt = the estimator-free energy "
            "balance against integral M_cum^2/(2GJ) (1 = Bredt); r(z) = theta'(z)*GJ/M_cum(z), "
            "whose tip blow-up is a metric artefact of dividing by M_cum -> 0"
        )
        print(f"{'mode':>23} {'n_z':>5} {'dz [m]':>8} {'W/Wbredt':>10} {'median r':>10} {'worst |r-1|':>12}")
        for run in runs:
            for label, metrics in run["modes"].items():
                print(
                    f"{label:>23} {run['n_z']:>5} {run['dz']:>8.4f} {metrics['energy_ratio']:>10.5f} "
                    f"{metrics['median']:>10.5f} {metrics['worst']:>12.4%}"
                )
        print(
            "\n(control below = the validated self-equilibrated end couple, record 1.00309; the "
            "validated bound covers the end-couple mode, production runs the distributed one)"
        )
        for run in runs:
            print(f"  control n_z={run['n_z']:>3}: window rate/Bredt = {run['control']:.5f}")
        return 0

    run = run_mesh(args.n_z, verbose=True)
    print(
        f"\ntube: {run['rings']} rings, L={tube.L} m, dz={run['dz']:.4f} m, "
        f"T={tube.TORQUE:.3e} N.m"
    )
    print(f"control (validated end couple, self-equilibrated) window rate/Bredt = {run['control']:.5f}")
    for label, metrics in run["modes"].items():
        print(
            f"[{label:22s}] r median {metrics['median']:.5f}  worst |r-1| {metrics['worst']:.4%} "
            f"({_verdict(metrics)} the {tube.TOL:.0%} rule)  |  above {RAMP_FLOOR:.0%}T: "
            f"median {metrics['median_loaded']:.5f}, worst {metrics['worst_loaded']:.4%}"
        )
        print(
            f"{'':24s}energy balance W / integral M_cum^2/(2GJ) = {metrics['energy_ratio']:.5f} "
            f"({(metrics['energy_ratio'] - 1.0) * 100:+.3f} %)"
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
