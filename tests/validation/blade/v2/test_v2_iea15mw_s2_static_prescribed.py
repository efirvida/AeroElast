"""S-2 on the windIO v2 blade: static prescribed loads vs OpenFAST BeamDyn p60 (v2 counterpart of group 39).

Same cases (LC1 gravity, LC2/LC3 1 MN tip, LC4 5 kN/m distributed), same load definitions
(``run_s2_static_cases.run_shell_case``: tip ``equal``/``arc``, LC4 ``equal``/``area``, 60 slices), same BeamDyn
p60 references (tip uy/ux, interior station nodes of tests/IEA15MW/bd_static), same metrics, same tolerance
(K3 = 5 %, the L3 ceiling) and the same predeclared station rule as ``tests/test_iea15mw_s2_static_prescribed.py``
(official blade only; there is no UTD input on the v2 side).  Only the mesh differs:
``build_v2_blade(0.25, skin_reference="mid")`` instead of ``BladeMesh(...).generate(renumber="rcm")``.

What had to be adapted (v2 rings are constant-eta rows, not constant-z planes: the twisted blade's rings are
not planar, tip ring spans ~0.1 m in z):
  * tip load: ``run_shell_case`` selects ``z >= z.max() - 1e-6`` (one node on a v2 mesh); here the tip ring is
    ``TipNodes`` (highest eta) and the ``arc`` weights are the v1 ``ring_arc_weights`` of that ring.
  * ring profile (S12): rings are the eta rows (nodes sharing an eta, skin and web); the weights are the v1 ones
    (half of each element side with both corners on the ring), the ring z is the arc-weighted mean z.
  * tip chord: largest in-plane extent of the TipNodes ring.
  * the in-house EB beam comparison (printed only in v1) is not repeated; it is built from z-sliced v1 geometry.
LC1 and LC4 are asserted at the tip node and at the station (as v1); LC2 at the station only (T5i); LC3 is printed.

See docs/validation/_pending/64-v2_blade_static_prescribed/NOTES.md.
"""

from __future__ import annotations

import sys

import numpy as np
import pytest

pytest.importorskip("petsc4py", reason="PETSc not available")

from scipy.spatial import cKDTree  # noqa: E402

from tests.support.assertions import assert_residual_below  # noqa: E402
from tests.support.paths import REPO_ROOT  # noqa: E402
from tests.support.s2_station import (  # noqa: E402
    bd_station_nodes,
    component_index,
    pick_station_node,
    ring_mean_at_z,
)
from tests.support.v2_blade import build_v2_blade, load_raw  # noqa: E402

sys.path.insert(0, str(REPO_ROOT / "tools"))

from run_s2_static_cases import ring_arc_weights, run_shell_case, shell_tip_displacement  # noqa: E402

ELEMENT_SIZE = 0.25
SKIN_REFERENCE = "mid"
GRAVITY = (0.0, -9.81, 0.0)
TIP_FORCE_N = 1.0e6
DISTRIB_Q_NM = 5.0e3

_TOL = 0.05  # K3, identical to v1 (rel_tol_dominant and _STATION_TOL)
# OpenFAST v5.0.0 BeamDyn driver, p60, linear static, official deck: tip uy [m] (K9), same numbers as v1.
_BD_STATIC_UY = {"LC1_gravity": -2.231483, "LC2_tip_flap": 82.273278, "LC4_distributed": 9.062631}
_NOT_ASSERTED_VS_BD = {"LC2_tip_flap"}
_BD_STATIC_LC3_UX = 34.469441  # K11.3: printed, not asserted (v1 _LC3_ASSERTED_VS_BD = False)
_BD_STATIC_DIR = REPO_ROOT / "tests" / "IEA15MW" / "bd_static"
_BD_SCALE_SPREAD = 1.0e-3
# LC2 joins at the interior station (odd/tasks/gmsh-chord-scaled-mesh.md T5i, L63-L64): the tip-introduction effect decays
# within ~one section depth, the station sits 4.35 tip chords inboard and is load-distribution independent; measured +0.2271 %
# vs BeamDyn node 56 (es 0.25, 2026-11-13).  The LC2 tip ring stays reported only.  LC3 stays diagnostic: +10.22 % (open).
_STATION_ASSERTED = ("LC1", "LC2", "LC4")

_TIP_CASES = {
    "LC1_gravity": {"type": "gravity", "value": GRAVITY},
    "LC2_tip_flap": {"type": "tip", "value": (0.0, TIP_FORCE_N, 0.0)},
    "LC3_tip_edge": {"type": "tip", "value": (TIP_FORCE_N, 0.0, 0.0)},
    "LC4_distributed": {"type": "distributed", "value": (0.0, DISTRIB_Q_NM, 0.0), "n_slices": 60},
}
_STATION_CASES = {
    "LC1": ("LC1_gravity", {"type": "gravity", "value": GRAVITY}),
    "LC2": ("LC2_tip_flap", {"type": "tip", "value": (0.0, TIP_FORCE_N, 0.0), "distribution": "arc"}),
    "LC3": ("LC3_tip_edge", {"type": "tip", "value": (TIP_FORCE_N, 0.0, 0.0), "distribution": "arc"}),
    "LC4": ("LC4_distributed", {"type": "distributed", "value": (0.0, DISTRIB_Q_NM, 0.0),
                                "n_slices": 60, "distribution": "area"}),
}

_CACHE: dict = {}


def _mesh():
    if "mesh" not in _CACHE:
        _CACHE["mesh"] = build_v2_blade(ELEMENT_SIZE, skin_reference=SKIN_REFERENCE)
    return _CACHE["mesh"]


def _tip_idxs(mesh) -> np.ndarray:
    idx = mesh.node_id_to_index
    return np.array(sorted(idx[nid] for nid in mesh.get_node_set("TipNodes").nodes), dtype=int)


def _eta_rows(mesh) -> np.ndarray:
    """eta of every mesh node (renumbered); the raw v2 elements are matched to mesh elements by centroid."""
    raw = load_raw(ELEMENT_SIZE)
    idx = mesh.node_id_to_index
    coords = mesh.coords_array
    cen_m, nodes_m = [], []
    for e in mesh.elements:
        n = [idx[nd.id] for nd in e.nodes]
        nodes_m.append(n)
        cen_m.append(coords[n].mean(axis=0))
    from tests.support.v2_blade import _inset

    x = _inset(raw)
    cen_r = np.array([x[el[el >= 0]].mean(axis=0) for el in raw["elements"]])
    dist, match = cKDTree(np.array(cen_m)).query(cen_r)
    assert dist.max() < 1e-9 and len(np.unique(match)) == len(match), "raw/mesh element matching is not a bijection"
    eta = np.full(len(coords), np.nan)
    for k, el in enumerate(raw["elements"]):
        ids = el[el >= 0]
        nm = nodes_m[match[k]]
        for j, i in enumerate(ids):
            # corner j of the raw element -> the mesh node of that element at the same position
            c = int(nm[int(np.argmin(np.linalg.norm(coords[nm] - x[i], axis=1)))])
            eta[c] = raw["params"][k, j, 0]
    assert not np.isnan(eta).any()
    return eta


def _ring_profile_eta(mesh, z_min):
    """v1 ``ring_profile`` on eta rows: (ring z, node indices, arc weights) of every row with mean z >= z_min."""
    coords = mesh.coords_array
    eta = _eta_rows(mesh)
    levels = np.unique(eta)
    row_of = np.searchsorted(levels, eta)
    idx = mesh.node_id_to_index
    edges = set()
    for e in mesh.elements:
        nodes = [idx[n.id] for n in e.nodes]
        for a, b in zip(nodes, nodes[1:] + nodes[:1], strict=True):
            if row_of[a] == row_of[b]:
                edges.add((min(a, b), max(a, b)))
    ed = np.array(sorted(edges), dtype=int)
    half = 0.5 * np.linalg.norm(coords[ed[:, 0]] - coords[ed[:, 1]], axis=1)
    w = np.zeros(len(coords))
    np.add.at(w, ed[:, 0], half)
    np.add.at(w, ed[:, 1], half)
    rings = []
    for k in range(len(levels)):
        nodes = np.nonzero((row_of == k) & (w > 0))[0]
        if len(nodes) == 0:
            continue
        wk = w[nodes]
        zk = float((coords[nodes, 2] * wk).sum() / wk.sum())
        if zk >= z_min:
            rings.append((zk, nodes, wk))
    return rings


def _profile_means(rings, u):
    uu = np.asarray(u).reshape(-1, 6)[:, :3]
    return [[z, *((uu[n] * w[:, None]).sum(axis=0) / w.sum()).tolist(), int(len(n))] for z, n, w in rings]


_SLAB_HALF_WIDTH_M = 0.25


def _slab_mean(mesh, u, z: float, half_width: float = _SLAB_HALF_WIDTH_M) -> np.ndarray:
    """Tributary-area weighted mean translation of the nodes with ``|z - z0| < half_width``.

    A v2 ring is a constant-eta row, not a constant-z plane: it is twisted and its z spread and spacing change with the element
    size, so the interpolated ring mean drifts with h (0.25 / 0.177 / 0.15 m: LC4 station -2.82 / -2.98 / -3.48 % while this slab
    mean stays -2.65 / -2.69 / -2.56 %, logs/check/gmsh_conv/tmetric-*.json).  On a v1 mesh the rings are z planes and both agree."""
    coords, idx = mesh.coords_array, mesh.node_id_to_index
    area = np.zeros(len(coords))
    for e in mesh.elements:
        n = [idx[nd.id] for nd in e.nodes]
        p = coords[n]
        a = 0.5 * np.linalg.norm(np.cross(p[2] - p[0], p[1] - p[0])) if len(n) == 3 else 0.5 * np.linalg.norm(np.cross(p[2] - p[0], p[3] - p[1]))
        area[n] += a / len(n)
    sel = np.abs(coords[:, 2] - z) < half_width
    if not sel.any():
        raise ValueError(f"no nodes within {half_width} m of z {z:.3f}")
    uu = np.asarray(u).reshape(-1, 6)[:, :3]
    return (uu[sel] * area[sel, None]).sum(axis=0) / area[sel].sum()


def _run(mesh, properties, load) -> np.ndarray:
    """``run_shell_case``; tip loads act on the v2 tip ring (TipNodes) instead of ``z >= z.max() - 1e-6``."""
    if load["type"] != "tip":
        return np.asarray(run_shell_case(mesh, properties, load))
    from aeroelast.core.bc import DirichletCondition, NodalLoad
    from aeroelast.elements import ElementFamily
    from aeroelast.solvers.elasticity.static_linear import StaticLinearSolver

    cfg = {"solver": {}, "elements": {"element_family": ElementFamily.SHELL, "span_direction": (0.0, 0.0, 1.0),
                                      "properties": properties}}
    solver = StaticLinearSolver(mesh, cfg)
    dpn = solver.domain.dofs_per_node
    idx = mesh.node_id_to_index
    root = sorted(idx[nid] * dpn + d for nid in mesh.get_node_set("RootNodes").nodes for d in range(dpn))
    solver.add_dirichlet_conditions([DirichletCondition(root, 0.0)])
    tip = _tip_idxs(mesh)
    force = np.asarray(load["value"], dtype=np.float64)
    if load.get("distribution", "equal") == "arc":
        tip, w = ring_arc_weights(mesh, tip)
        f_nodes = np.outer(w / w.sum(), force)
    else:
        f_nodes = np.tile(force / len(tip), (len(tip), 1))
    dofs, vals = [], []
    for i, f in zip(tip, f_nodes, strict=True):
        for d in range(3):
            dofs.append(int(i) * dpn + d)
            vals.append(float(f[d]))
    solver.add_nodal_loads([NodalLoad(dofs, vals)])
    return np.asarray(solver.solve())


def _tip_results():
    if "tip" not in _CACHE:
        mesh, props = _mesh()
        out = {}
        for name, load in _TIP_CASES.items():
            u = _run(mesh, props, load)
            tip = _tip_idxs(mesh)
            out[name] = {"shell": shell_tip_displacement(mesh, u),
                         "ring": np.asarray(u).reshape(-1, 6)[tip, :3].mean(axis=0)}
        _CACHE["tip"] = out
    return _CACHE["tip"]


@pytest.mark.slow
class TestV2S2StaticPrescribed:
    @pytest.mark.parametrize("case_name", ["LC1_gravity", "LC2_tip_flap", "LC4_distributed"])
    def test_dominant_component_matches_beam(self, case_name):
        d_shell = _tip_results()[case_name]["shell"][1]
        d_bd = _BD_STATIC_UY[case_name]
        print(f"S-2 v2 {case_name} uy [m]: shell {d_shell:.6f} vs OpenFAST BeamDyn p60 {d_bd:.6f} "
              f"(signed error: {d_shell / d_bd - 1:+.4%}, bound: {_TOL:.4%})")
        if case_name in _NOT_ASSERTED_VS_BD:
            return
        rel = abs(d_shell - d_bd) / abs(d_bd)
        assert rel < _TOL, (f"[v2] {case_name} uy: shell={d_shell:.4f} BeamDyn p60={d_bd:.4f} "
                            f"(rel={rel:.2%}, tol={_TOL:.0%})")
        if case_name == "LC1_gravity":
            assert d_shell < 0, "[v2] LC1_gravity: expected downward deflection"

    def test_edgewise_lc3_tip_ux(self):
        """LC3 edgewise tip ux vs BeamDyn p60 (K11.3): printed, not asserted (as v1)."""
        res = _tip_results()["LC3_tip_edge"]
        d_shell, d_ring, d_bd = res["shell"][0], res["ring"][0], _BD_STATIC_LC3_UX
        print(f"S-2 v2 LC3_tip_edge ux [m]: shell {d_shell:.6f} vs OpenFAST BeamDyn p60 {d_bd:.6f} "
              f"(signed error: {d_shell / d_bd - 1:+.4%}, bound: {_TOL:.4%}); tip ring mean {d_ring:.6f} "
              f"(node vs ring {d_shell / d_ring - 1:+.4%})")


@pytest.mark.slow
class TestV2S2InteriorStation:
    """S12: BeamDyn p60 node displacement vs the shell z-slab mean at the node's z (the eta-row ring mean is printed, not asserted)."""

    @pytest.fixture(scope="class")
    def station_profiles(self):
        if "station" not in _CACHE:
            mesh, props = _mesh()
            z_tip = float(mesh.coords_array[:, 2].max())
            c = mesh.coords_array[_tip_idxs(mesh), :2]
            chord = float(np.max(np.linalg.norm(c[:, None, :] - c[None, :, :], axis=-1)))
            rings = _ring_profile_eta(mesh, z_tip - 8.0)
            sol = {case: _run(mesh, props, load) for case, (_n, load) in _STATION_CASES.items()}
            out = {case: _profile_means(rings, u) for case, u in sol.items()}
            _CACHE["station"] = (out, z_tip, chord, mesh, sol)
        return _CACHE["station"]

    @pytest.mark.parametrize("case", list(_STATION_CASES))
    def test_station_matches_beamdyn_node(self, station_profiles, case):
        profiles, z_tip, chord, mesh, sol = station_profiles
        nodes, spread = bd_station_nodes(_BD_STATIC_DIR, case)
        assert spread < _BD_SCALE_SPREAD, f"{case}: BeamDyn scales 1e-4 / 2e-4 differ by {spread:.3e}"
        node = pick_station_node(nodes, z_tip=z_tip, chord=chord)
        k = component_index(case)
        comp = "xy"[k]
        d_bd = node[f"linear_shell_u{comp}_m"]
        d_ring = ring_mean_at_z(profiles[case], node["z"])[k]
        d_shell = _slab_mean(mesh, sol[case], node["z"])[k]
        print(f"S-2 v2 {case} station u{comp} [m]: z-slab mean {d_shell:.6f} (asserted metric), eta-row ring mean {d_ring:.6f} "
              f"(comparison, drifts with h), BeamDyn {d_bd:.6f}")
        rel = abs(d_shell - d_bd) / abs(d_bd)
        if case not in _STATION_ASSERTED:
            print(f"S-2 v2 {case} station u{comp} [m] diagnostic (S13, not asserted: tip point load without MPC): "
                  f"shell {d_shell:.6f} vs BeamDyn p60 node {node['node']} {d_bd:.6f}, signed {d_shell / d_bd - 1:+.4%}")
            return
        assert_residual_below(
            rel,
            tol=_TOL,
            reference_name=f"OpenFAST v5.0.0 BeamDyn p60 node {node['node']} (z {node['z']:.3f} m, {(z_tip - node['z']) / chord:.2f} chords)",
            kind="code",
            what=f"S-2 v2 {case} station u{comp} [m] shell z-slab mean (|dz|<{_SLAB_HALF_WIDTH_M} m, area weighted) {d_shell:.6f} (BeamDyn {d_bd:.6f}, signed {d_shell / d_bd - 1:+.4%})",
        )
