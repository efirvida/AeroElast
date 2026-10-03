"""An independent 1D beam twist from the IEA-15MW anchor's own section stiffness.

The beam's ``GJ`` and shear centre come from the official IEA-15-240-RWT BeamDyn blade deck via
``K66toPropsDecoupled(convention='BeamDyn')``, its loads from the production BEM +
``ForceProjector`` the shell is solved with, and its twist from the root-fixed cantilever form
``theta(z) = int_0^z T(s)/GJ(s) ds``, ``T(z) = int_z^R m ds``, with
``m = Mp + (x_AC - xS) Np - (y_AC - yS) Tp``.  ``xS, yS, GJ`` are the **anchor's**, not the
shell's: the two models share *loads, not section properties*.

The 6x6 reading is pinned by the ElastoDyn cross-check (both bending stiffnesses within 2 % of
``FlpStff``/``EdgStff``).  ``x_AC`` uses the measured sign ``(pitch_axis - 0.25) * chord`` of the
production AC's chordwise projection (the task's written sign inverts it, a 10x torque error).
The point-like tip interval is excluded by a ``GKt < median/1e3`` guard.  Zhou et al. 2025
Table 4 is printed beside the columns, never asserted.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("ccblade", reason="ccblade not installed (pip install -e '.[bem]')")
pytest.importorskip("_aeroelast", reason="Rust backend not available")

_OPENFAST_TOOLBOX = Path.home() / "openfast_toolbox"
if _OPENFAST_TOOLBOX.is_dir() and str(_OPENFAST_TOOLBOX) not in sys.path:
    sys.path.insert(0, str(_OPENFAST_TOOLBOX))

from openfast_toolbox.converters.beam import K66toPropsDecoupled  # noqa: E402

from _aeroelast import PyMeshAssembler  # noqa: E402
from scipy.sparse import coo_matrix  # noqa: E402
from scipy.sparse.linalg import spsolve  # noqa: E402

from aeroelast.core.mesh.entities import MeshElement, Node  # noqa: E402
from aeroelast.models.blade.model import Blade  # noqa: E402
from aeroelast.solvers.bem.fsi_participant import BEMFSIParticipant  # noqa: E402

import tests.validation.blade.test_blade_iea15mw_validation as blade_validation  # noqa: E402
from tests.support.openfast_bem import build_blade_aero_from_aerodyn  # noqa: E402
from tests.validation.blade.test_blade_rated_twist import STATION_GAP_TOLERANCE, _physical_stations  # noqa: E402

from tests.support.paths import DATA_DIR, SOURCES_DIR  # noqa: E402
_SOURCES = (SOURCES_DIR / "openfast" / "iea15mw"
            / "IEA-15-240-RWT")
BEAMDYN_BLADE = _SOURCES / "IEA-15-240-RWT_BeamDyn_blade.dat"
ELASTODYN_BLADE = _SOURCES / "IEA-15-240-RWT_ElastoDyn_blade.dat"
YAML = DATA_DIR / "IEA-15-240-RWT.yaml"
AD_PRIMARY = (DATA_DIR / "reference" / "iea15mw_openfast" / "case"
              / "IEA-15-240-RWT_AeroDyn15.dat")

V_RATED, RPM_RATED, PITCH_RATED = 10.59, 7.56, 0.0
BEM_CONFIG = {
    "wind_speed": V_RATED, "omega": RPM_RATED * 2.0 * np.pi / 60.0, "pitch": PITCH_RATED,
    "azimuth": 0.0, "air_density": 1.225, "dynamic_viscosity": 1.81206e-5,
    "hub_height": 150.0, "shear_exp": 0.0,
}
EI_BENDING_TOL = 0.02          # the ElastoDyn cross-check that fixes the 6x6 reading
GKT_MEDIAN_FACTOR = 1.0e3      # a GKt this far below the span median is degenerate
ZHOU_TIP_TORSION_DEG = -3.60    # Zhou et al. 2025 Table 4 - reported, never asserted
ZHOU_TIP_FLAP_M, ZHOU_TIP_EDGE_M = 13.86, -1.22


def _numeric_rows(path: Path, marker: str) -> list[list[float]]:
    """Every all-float line after the header line containing ``marker``."""
    lines = path.read_text(errors="replace").splitlines()
    start = next(i for i, ln in enumerate(lines) if marker in ln)
    rows: list[list[float]] = []
    for ln in lines[start + 1:]:
        body = ln.split("!", 1)[0].strip()
        if not body:
            continue
        try:
            rows.append([float(t) for t in body.split()])
        except ValueError:
            continue
    return rows


def _ring_section_rotation(coords: np.ndarray, disp: np.ndarray, ring: np.ndarray) -> float:
    """Antisymmetric part of the ring's least-squares in-plane affine fit [rad]."""
    pts = coords[ring]
    design = np.column_stack([pts[:, 0], pts[:, 1], np.ones(len(ring))])
    a12 = np.linalg.lstsq(design, disp[ring, 0], rcond=None)[0][1]
    a21 = np.linalg.lstsq(design, disp[ring, 1], rcond=None)[0][0]
    return 0.5 * (a21 - a12)


@pytest.fixture(scope="module")
def anchor():
    """BeamDyn 6x6 stiffness -> decoupled props, plus the ElastoDyn cross-check deck."""
    if not BEAMDYN_BLADE.exists() or not ELASTODYN_BLADE.exists():
        pytest.skip(f"OpenFAST deck not present under {_SOURCES}")
    rows = _numeric_rows(BEAMDYN_BLADE, "DISTRIBUTED PROPERTIES")
    n = len(rows) // 13  # 1 span + 6x6 stiffness + 6x6 mass per station
    frac = np.array([rows[i * 13][0] for i in range(n)])
    K = np.array([np.asarray(rows[i * 13 + 1:i * 13 + 7]) for i in range(n)])
    names = ["EA", "EIxp", "EIyp", "kxsGA", "kysGA", "GKt", "xC", "yC", "xS", "yS",
             "theta_p", "theta_s"]
    props = [K66toPropsDecoupled(K[i], convention="BeamDyn") for i in range(n)]
    out = {name: np.array([p[i] for p in props]) for i, name in enumerate(names)}
    ed = np.asarray(_numeric_rows(ELASTODYN_BLADE, "DISTRIBUTED BLADE PROPERTIES"))
    out.update(frac=frac, K=K, pitch_frac=ed[:, 0], pitch_axis=ed[:, 1],
               flp_stff=ed[:, 4], edg_stff=ed[:, 5])
    return out


@pytest.fixture(scope="module")
def production():
    """Real mesh + production participant + one rated shell solve under its own loads."""
    if not AD_PRIMARY.exists():
        pytest.skip(f"AeroDyn reference not present: {AD_PRIMARY}")
    Node._id_counter = 0
    MeshElement._id_counter = 0
    model = Blade(str(YAML), element_size=1.0)
    model.generate_mesh()
    mesh = model.mesh
    props = model.get_element_properties()
    assembler = PyMeshAssembler.from_model(
        blade_validation._to_rust_mesh(mesh, props), props,
        list(blade_validation.SPAN_DIRECTION), None)
    n = assembler.dofs_count
    rows, cols, vals = assembler.assemble_k()
    K = coo_matrix((np.asarray(vals), (np.asarray(rows), np.asarray(cols))),
                   shape=(n, n)).tocsr()
    root = {mesh.node_id_to_index[nid] for nid in mesh.get_node_set("RootNodes").node_ids}
    free = np.array([i for i in range(n) if i not in {6 * r + d for r in root for d in range(6)}],
                    dtype=np.int64)
    coords = np.array([[nd.x, nd.y, nd.z] for nd in mesh.nodes], dtype=float)
    blade_aero = build_blade_aero_from_aerodyn(AD_PRIMARY)
    participant = BEMFSIParticipant(mesh, blade_aero, BEM_CONFIG)
    forces_rigid, bem_rigid = participant._compute_forces(np.zeros((len(coords), 3)))
    load = np.zeros(n)
    for dof in range(3):
        load[dof::6] = forces_rigid[:, dof]
    u = np.zeros(n)
    u[free] = np.asarray(spsolve(K[np.ix_(free, free)], load[free]))
    phys_stations = _physical_stations(coords)
    phys_rings = {zz: np.where(np.abs(coords[:, 2] - zz) < STATION_GAP_TOLERANCE)[0]
                  for zz in phys_stations}
    return {"blade_aero": blade_aero, "coords": coords, "u": u,
            "displacements": np.column_stack([u[0::6], u[1::6], u[2::6]]),
            "bem_rigid": bem_rigid, "phys_stations": phys_stations, "phys_rings": phys_rings}


def _anchor_beam(anchor: dict, production: dict) -> dict:
    """The anchor-stiffness beam driven by the production loads, at the anchor stations."""
    bem, blade_aero = production["bem_rigid"], production["blade_aero"]
    r_hub = np.asarray(blade_aero.r, dtype=float)
    frac = anchor["frac"]
    r_anchor = blade_aero.hub_radius + frac * blade_aero.blade_length
    # Explicit interpolation of the production loads onto the anchor stations.
    Np = np.interp(r_anchor, r_hub, bem.Np)
    Tp = np.interp(r_anchor, r_hub, bem.Tp)
    Mp = np.interp(r_anchor, r_hub, bem.Mp)
    chord = np.interp(r_anchor, r_hub, bem.chord)
    pitch = np.interp(frac, anchor["pitch_frac"], anchor["pitch_axis"])
    x_ac = (pitch - 0.25) * chord          # AC offset from the pitch axis; y_AC = 0
    m = Mp + (x_ac - anchor["xS"]) * Np - (0.0 - anchor["yS"]) * Tp

    GKt = anchor["GKt"]
    floor = float(np.median(GKt)) / GKT_MEDIAN_FACTOR
    interval_ok = np.minimum(GKt[:-1], GKt[1:]) >= floor
    dr = np.diff(r_anchor)
    torque = np.zeros_like(m)              # T(z) = int_z^R m ds (root-fixed, tip-guarded)
    for i in range(len(m) - 2, -1, -1):
        torque[i] = torque[i + 1] + (0.5 * (m[i] + m[i + 1]) * dr[i] if interval_ok[i] else 0.0)
    twist = np.zeros_like(m)               # theta(z) = int_0^z T/GJ ds
    for i in range(1, len(m)):
        # Interval [i-1, i]: its length is dr[i-1] and its guard is interval_ok[i-1].  Using
        # dr[i] here measures a different interval than the one being admitted.
        step = (0.5 * (torque[i - 1] / GKt[i - 1] + torque[i] / GKt[i]) * dr[i - 1]
                if interval_ok[i - 1] else 0.0)
        twist[i] = twist[i - 1] + step
    i = len(m) - 2
    tip_torque = 0.5 * (m[i] + m[i + 1]) * dr[i]
    return {"frac": frac, "m": m, "torque": torque, "twist_deg": np.rad2deg(twist),
            "GKt": GKt, "interval_ok": interval_ok, "floor": floor,
            "tip_artefact_deg": float(np.rad2deg(0.5 * (tip_torque / GKt[i]) * dr[i]))}


def test_anchor_sections_cross_check_and_degenerate_tip(anchor):
    """Pin the 6x6 reading against ElastoDyn and guard the point-like tip interval.

    Asserted identities: ``K[2,2]`` is the axial term, ``K[5,5]`` the raw torsion, both
    decoupled bending stiffnesses reproduce ``FlpStff``/``EdgStff`` within 2 % (which is what
    rejects the torsion-at-index-3 reading), and only the last interval falls below the
    ``median/1e3`` degenerate-stiffness floor.
    """
    i = 0
    EA, GKt = anchor["EA"][i], anchor["GKt"][i]
    flp, edg = anchor["flp_stff"][i], anchor["edg_stff"][i]
    err_x = abs(anchor["EIxp"][i] - flp) / flp
    err_y = abs(anchor["EIyp"][i] - edg) / edg
    bad = abs(anchor["K"][i, 5, 5] - flp) / flp
    print(f"\nanchor root: EA {EA:.6e} N, GKt {GKt:.6e} N.m^2; "
          f"EIxp vs FlpStff {err_x:+.2%}; EIyp vs EdgStff {err_y:+.2%}; "
          f"torsion-at-index-3 reading misses FlpStff by {bad:.1%}")
    assert abs(EA - 4.605e10) / 4.605e10 < 0.01
    assert abs(GKt - 8.749e10) / 8.749e10 < 0.01
    assert abs(anchor["K"][i, 2, 2] - EA) / EA < 1e-12
    assert abs(anchor["K"][i, 5, 5] - GKt) / anchor["K"][i, 5, 5] < 0.01  # decoupled vs raw
    assert err_x < EI_BENDING_TOL and err_y < EI_BENDING_TOL and bad > EI_BENDING_TOL

    GKt = anchor["GKt"]
    floor = float(np.median(GKt)) / GKT_MEDIAN_FACTOR
    interval_ok = np.minimum(GKt[:-1], GKt[1:]) >= floor
    print(f"degenerate-tip guard: median/1e3 = {floor:.3e} N.m^2, GKt[-2] = {GKt[-2]:.3e}, "
          f"GKt[-1] = {GKt[-1]:.3e}; excluded intervals {np.where(~interval_ok)[0].tolist()}")
    assert GKt[-1] < floor <= GKt[-2]
    assert not interval_ok[-1] and interval_ok[:-1].all()


def test_anchor_beam_arbitrates_the_shell_section_rotation(anchor, production):
    """Print the anchor-beam / shell table and state which estimator the beam supports.

    The only assertion is an invariant: a root-fixed beam's twist magnitude is monotone
    non-decreasing toward the tip.  Every shell number is reported, never fitted to a bound.
    """
    beam = _anchor_beam(anchor, production)
    coords, u, disp = production["coords"], production["u"], production["displacements"]
    stations = production["phys_stations"]
    L = production["blade_aero"].blade_length

    rows = []
    for i in range(len(anchor["frac"])):
        r_root = float(anchor["frac"][i] * L)
        zz = stations[int(np.argmin(np.abs(np.asarray(stations) - r_root)))]
        ring = production["phys_rings"][zz]
        rows.append((r_root, zz, len(ring), beam["twist_deg"][i],
                     np.rad2deg(_ring_section_rotation(coords, disp, ring)),
                     np.rad2deg(float(np.mean(u[6 * ring + 5])))))

    print("\nanchor beam (BeamDyn K66) vs the shell's two section estimators "
          "(degrees about +span):")
    print(f"  {'r_root[m]':>9} {'r_shell[m]':>10} {'nodes':>5} {'phi_beam':>9} "
          f"{'ring rot':>9} {'mean tz':>9} {'ring/beam':>9} {'meanz/beam':>10}")
    for r_root, zz, count, phi, ring_rot, mean_tz in rows:
        denom = phi if abs(phi) > 1e-12 else float("nan")
        print(f"  {r_root:>9.2f} {zz:>10.2f} {count:>5} {phi:>9.4f} {ring_rot:>9.4f} "
              f"{mean_tz:>9.4f} {ring_rot / denom:>9.4f} {mean_tz / denom:>10.4f}")
    tip = rows[-1]
    print(f"  tip: phi_beam {tip[3]:+.4f}, ring rotation {tip[4]:+.4f} "
          f"(ratio {tip[4] / tip[3]:.4f}), mean theta_z {tip[5]:+.4f} "
          f"(ratio {tip[5] / tip[3]:.4f})")
    print(f"  excluded tip interval would have added {beam['tip_artefact_deg']:+.4f} deg")
    print(f"  Zhou 2025 Table 4 (reported, NOT asserted): torsion {ZHOU_TIP_TORSION_DEG:+.2f} deg,"
          f" flap {ZHOU_TIP_FLAP_M:+.2f} m, edge {ZHOU_TIP_EDGE_M:+.2f} m")

    assert np.all(np.diff(np.abs(beam["twist_deg"])) >= -1e-9), "twist is not monotone"
    assert beam["tip_artefact_deg"] < 0.0  # the excluded tip interval is nose-down
