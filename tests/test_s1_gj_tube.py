"""S-1 GJ extraction and S-7 load on short prismatic tubes of known GJ (blade-input-audit T40).

T39 (``tools/diagnose_s1_gj_tube.py``) ran the S-1 extractor on a tube whose Saint-Venant torsion
constant is exact thin-wall Bredt, ``GJ = G t 4 A^2 / P``, and found two defects with no blade in
the loop: the tip torque ``a (t x r_perp)`` is not the section's wall shear flow on a rectangle, so it
excites section distortion that the rigid-fit twist reads as rotation (A: 257 % mean error), and the
absolute smoothing spline on theta costs ~5 points on a varying-GJ signal even when theta is exact
(D). The fix applies the torque as the multi-cell Saint-Venant wall flow
(``section_torque_wall_forces``) and differentiates theta with the local quadratic fit.

Truth tests (MITC4, element size 0.25 m, 24 m long, root clamped; reference = Bredt, which the MITC4
tube benchmarks reproduce to ~1 %, ``test_thin_walled_tube_torsion``; T39 case C reached 0.01 %):

The rectangle of (a) and (c) is the warping-free one (4 m walls t = 0.02, 2 m walls t = 0.01, so
b/t_b = h/t_h): the uniform-thickness rectangle warps in Saint-Venant torsion, the clamped root restrains
it, and that end layer couples into distortion that decays over tens of metres. That is physics, not the
extractor; ``test_uniform_rectangle_root_layer_diagnostic`` prints it without asserting.

* (a) rectangle 4 x 2 m, uniform wall: S-1 GJ within 1 % of Bredt at every mid-span station;
* (b) circle of the same perimeter, GJ(z) falling ~10x: S-1 GJ mid-span mean |err| < 1 %;
* (c) the S-7 construction (``run_torsion_case`` + ``section_twists_deg``) on the rectangle:
  per-slice twist ratio to the exact ``T int dz/GJ`` within 1 % over its 0.3-0.9 band;
* (e) oblique circular tube (axis 20 deg off z in xz, rings cut by z-planes: ellipses): S-1 GJ
  within 1 % of ``G t 2 pi R^3`` at every mid-span station. Needs the end-ring-normal load and the
  6-DOF tangent-frame twist (blade-input-audit T40b);
* (d) the load helper: delivered torque exact, zero net force, geometric-only refused.
"""

from __future__ import annotations

import sys

import numpy as np
import pytest

pytest.importorskip("_aeroelast", reason="Rust backend not available")

from tests.support.paths import REPO_ROOT  # noqa: E402

sys.path.insert(0, str(REPO_ROOT / "tools"))

import diagnose_s1_gj_tube as tube  # noqa: E402
import run_s1_sectional as s1  # noqa: E402
import run_s7_torsion as s7  # noqa: E402

LENGTH = 24.0  # tube length [m]
TORQUE = 1.0e6  # the S-1 / S-7 torque [N m]
BOUND = 0.01  # Bredt / MITC4 tube benchmark agreement
STATIONS = np.arange(0.5, LENGTH, 1.0)  # one S-1 station per metre, on the 0.25 m layers


def _layers() -> int:
    return int(round(LENGTH / tube.ES))


def _gj_uniform(shape: str) -> np.ndarray:
    return np.full(_layers(), tube.G_UNIFORM * tube.T_WALL * tube.bredt_factor(shape))


def _gj_falling(shape: str) -> np.ndarray:
    zmid = (np.arange(_layers()) + 0.5) * tube.ES
    return tube.G_UNIFORM * tube.T_WALL * tube.bredt_factor(shape) * np.exp(-np.log(10.0) * zmid / LENGTH)


def _build(shape: str, gj_layer: np.ndarray):
    g_layer = gj_layer / (tube.T_WALL * tube.bredt_factor(shape))
    mesh, props, zl = tube.build_tube(g_layer, shape, tube.T_WALL, length=LENGTH)
    mesh.coords_array = np.asarray([n.coords for n in mesh.nodes], dtype=float)
    return mesh, props, zl


def _s1_error(shape: str, gj_layer: np.ndarray):
    mesh, props, zl = _build(shape, gj_layer)
    zs, _ea, gj, _ = s1.ea_gj_from_static_response(mesh, props, torque=TORQUE, stations=STATIONS)
    layer = np.clip(np.searchsorted(zl, zs) - 1, 0, len(gj_layer) - 1)
    mid = (zs > 0.2 * zs.max()) & (zs < 0.85 * zs.max())
    return zs[mid], gj[mid] / gj_layer[layer][mid] - 1.0


def test_a_rectangle_uniform_gj_within_bredt():
    zs, err = _s1_error("rectw", _gj_uniform("rectw"))
    assert np.all(np.isfinite(err)), f"GJ is NaN at z={zs[~np.isfinite(err)]}"
    assert np.max(np.abs(err)) < BOUND, f"max |err| {np.max(np.abs(err)):.4%} (per station {err})"


def test_b_circle_varying_gj_mean_error():
    zs, err = _s1_error("circle", _gj_falling("circle"))
    assert np.all(np.isfinite(err)), f"GJ is NaN at z={zs[~np.isfinite(err)]}"
    assert np.mean(np.abs(err)) < BOUND, f"mean |err| {np.mean(np.abs(err)):.4%} (per station {err})"


def test_c_s7_construction_twist_ratio_per_slice():
    gj_layer = _gj_uniform("rectw")
    mesh, props, zl = _build("rectw", gj_layer)
    u = np.asarray(s7.run_torsion_case(mesh, props, TORQUE), dtype=float)
    zs, th = s7.section_twists_deg(mesh, u)
    flex = np.concatenate([[0.0], np.cumsum(TORQUE / gj_layer * np.diff(zl))])
    ratio = th / np.rad2deg(np.interp(zs, zl, flex))
    band = (zs >= 0.3 * LENGTH) & (zs <= 0.9 * LENGTH)
    assert band.sum() >= 10
    assert np.max(np.abs(ratio[band] - 1.0)) < BOUND, f"ratio min/max {ratio[band].min():.4f}/{ratio[band].max():.4f}"


BETA = np.deg2rad(20.0)


def _build_oblique():
    """Circular tube of radius R whose axis is tilted BETA from z in the xz-plane; rings in z-planes."""
    from aeroelast.core.material import IsotropicMaterial
    from aeroelast.core.mesh import ElementSet, MeshElement, MeshModel, Node, NodeSet
    from aeroelast.core.mesh.entities import ElementType
    from aeroelast.core.properties import ShellProperty

    Node._id_counter = 0
    MeshElement._id_counter = 0
    rad = 2.0 * (tube.B_MID + tube.H_MID) / (2.0 * np.pi)
    nr = int(round(2.0 * (tube.B_MID + tube.H_MID) / tube.ES))
    phi = 2.0 * np.pi * np.arange(nr) / nr
    nz = _layers()
    mesh = MeshModel()
    rings = []
    for z in np.linspace(0.0, LENGTH, nz + 1):
        s_ax = (z + rad * np.cos(phi) * np.sin(BETA)) / np.cos(BETA)
        x = s_ax * np.sin(BETA) + rad * np.cos(phi) * np.cos(BETA)
        row = []
        for xi, yi in zip(x, rad * np.sin(phi), strict=True):
            n = Node([float(xi), float(yi), float(z)], geometric_node=False)
            mesh.add_node(n)
            row.append(n)
        rings.append(row)
    g = tube.G_UNIFORM
    props = {}
    for k in range(nz):
        els = set()
        for i in range(nr):
            j = (i + 1) % nr
            e = MeshElement(nodes=[rings[k][i], rings[k][j], rings[k + 1][j], rings[k + 1][i]],
                            element_type=ElementType.quad)
            mesh.add_element(e)
            els.add(e)
        name = f"00_{k:03d}_tube"
        mesh.add_element_set(ElementSet(name, els))
        props[name] = ShellProperty(
            material=IsotropicMaterial(name=f"m{k}", E=2.0 * (1.0 + tube.NU) * g, nu=tube.NU, rho=1000.0),
            thickness=tube.T_WALL)
    mesh.add_node_set(NodeSet("RootNodes", set(rings[0])))
    mesh.coords_array = np.asarray([n.coords for n in mesh.nodes], dtype=float)
    return mesh, props, g * tube.T_WALL * 2.0 * np.pi * rad**3


def test_e_oblique_circular_tube_gj_within_bredt():
    mesh, props, gj_ref = _build_oblique()
    zs, _ea, gj, _ = s1.ea_gj_from_static_response(mesh, props, torque=TORQUE, stations=STATIONS)
    mid = (zs > 0.2 * zs.max()) & (zs < 0.85 * zs.max())
    err = gj[mid] / gj_ref - 1.0
    assert np.all(np.isfinite(err))
    assert np.max(np.abs(err)) < BOUND, f"max |err| {np.max(np.abs(err)):.4%} (per station {np.round(err, 5)})"


def test_uniform_rectangle_root_layer_diagnostic(capsys):
    """Printed, not asserted: the uniform-thickness clamped rectangle keeps the root warping layer."""
    zs, err = _s1_error("rect", _gj_uniform("rect"))
    with capsys.disabled():
        print(f"\n[diagnostic] uniform rect S-1 GJ error at z={np.round(zs, 1).tolist()}: "
              f"{np.round(100 * err, 2).tolist()} %")
    assert np.all(np.isfinite(err))


def _tip_ring(mesh):
    z = mesh.coords_array[:, 2]
    return np.nonzero(z >= z.max() - 1e-6)[0]


def test_d_helper_delivers_exact_torque_and_no_net_force():
    from aeroelast.solvers.bem.force_projection import section_torque_wall_forces

    mesh, props, _ = _build("rect", _gj_uniform("rect"))
    ring = _tip_ring(mesh)
    nodes, forces = section_torque_wall_forces(mesh, ring, (0.0, 0.0, 1.0), props, TORQUE)
    xyz = mesh.coords_array[nodes]
    r = xyz - xyz.mean(axis=0)
    delivered = float(np.sum(r[:, 0] * forces[:, 1] - r[:, 1] * forces[:, 0]))
    assert delivered == pytest.approx(TORQUE, rel=1e-9)
    assert np.abs(forces.sum(axis=0)).max() <= 1e-9 * np.abs(forces).max()
    assert np.abs(forces[:, 2]).max() == 0.0
    # the wall flow is the uniform Bredt q = T / 2A on every wall node pair, not a tangential field
    area = tube.B_MID * tube.H_MID
    assert np.abs(forces).sum(axis=1).max() <= TORQUE / (2.0 * area) * tube.ES * (1.0 + 1e-9)


def test_d_helper_refuses_geometric_only_stiffness():
    from aeroelast.solvers.bem.force_projection import section_torque_wall_forces

    mesh, _props, _ = _build("rect", _gj_uniform("rect"))
    with pytest.raises(ValueError, match="stiffness"):
        section_torque_wall_forces(mesh, _tip_ring(mesh), (0.0, 0.0, 1.0), None, TORQUE)
