"""MITC4 K on a folded, unsupported shell must stay positive semi-definite (T5, L14/L16).

Physical reference: an unconstrained linear elastic body has a PSD stiffness whose nullspace
is the 6 rigid-body modes plus MITC4's known drilling zero modes.  Two flat strips of quads
meeting along a common edge at dihedral angle theta (between the strip normals) must obey
that for every theta.  The drilling baseline is the zero-mode count of the flat theta = 0
control; a fold may not add negative eigenvalues.

Hypothesis under test: ``mitc4_nodal_directors`` (``assembly/assembler.rs``) area-averages
the element directors of *all* MITC4 elements at a node without a crease criterion, so at a
fold the director is theta/2 off each face normal and the assembled K is no longer a sum of
PSD element matrices.  The element-alone matrices are the control (a one-element model has
only its own director).
"""

from __future__ import annotations

import numpy as np
import pytest

_aeroelast = pytest.importorskip("_aeroelast")
from _aeroelast import (  # noqa: E402
    Laminate,
    MeshModel as RustMeshModel,
    OrthotropicMaterial,
    Ply,
    PyMeshAssembler,
)

NEG_TOL = 1e-10  # negative if lambda < -NEG_TOL * lambda_max
ZERO_TOL = 1e-9  # zero mode if |lambda| < ZERO_TOL * lambda_max

_MAT = OrthotropicMaterial(40e9, 10e9, 10e9, 5e9, 4e9, 5e9, 0.3, 0.3, 0.3, 1800.0)
PROPS = {
    "iso": dict(type="isotropic", e=70e9, nu=0.3, rho=2700.0, thickness=0.01,
                shear_correction=5 / 6, drilling_scale=1.0),
    "comp": Laminate([Ply(_MAT, 0.004, 0.0), Ply(_MAT, 0.004, 45.0), Ply(_MAT, 0.004, -45.0)], 0.75),
}
SPAN = [1.0, 0.0, 0.0]


def fold_mesh(theta_deg, length=1.0, width=1.0, nx=2, ns=2):
    """Two strips of nx * ns quads sharing the edge s = 0 (the x axis); returns (xyz, conn)."""
    th = np.radians(theta_deg)
    xs = np.linspace(0.0, length, nx + 1)
    ss = np.linspace(-width, width, 2 * ns + 1)
    xyz, idx = [], {}
    for i, x in enumerate(xs):
        for j, s in enumerate(ss):
            idx[i, j] = len(xyz)
            xyz.append([x, s, 0.0] if s <= 0 else [x, s * np.cos(th), s * np.sin(th)])
    conn = [[idx[i, j], idx[i + 1, j], idx[i + 1, j + 1], idx[i, j + 1]]
            for j in range(2 * ns) for i in range(nx)]
    return np.array(xyz), conn


def assemble(xyz, conn, set_name, type_code):
    ne = len(conn)
    mesh = RustMeshModel.from_raw_data(
        list(range(len(xyz))), np.asarray(xyz, float).ravel().tolist(),
        list(range(ne)), conn, [type_code] * ne, {set_name: list(range(ne))}, {})
    asm = PyMeshAssembler.from_model(mesh, {set_name: PROPS[set_name]}, SPAN, None)
    rows, cols, vals = asm.assemble_k()
    k = np.zeros((6 * len(xyz), 6 * len(xyz)))
    np.add.at(k, (np.asarray(rows), np.asarray(cols)), np.asarray(vals))
    return k


def min_ratio_and_zeros(k):
    w = np.linalg.eigvalsh(0.5 * (k + k.T))
    lmax = float(np.abs(w).max())
    return float(w[0] / lmax), int(np.sum(np.abs(w) < ZERO_TOL * lmax))


def element_ratios(xyz, conn, set_name, type_code):
    return [min_ratio_and_zeros(assemble(xyz[c], [[0, 1, 2, 3]], set_name, type_code))[0] for c in conn]


CASES = [
    ("iso", 4, 0.0, 1.0, 1.0), ("iso", 4, 30.0, 1.0, 1.0), ("iso", 4, 60.0, 1.0, 1.0),
    ("iso", 4, 90.0, 1.0, 1.0), ("comp", 44, 90.0, 1.0, 1.0),
    # TE sliver: short edge 1.8e-3 m across the fold, long edge 0.3 m along it (aspect 167)
    ("iso", 4, 90.0, 0.3, 1.8e-3), ("comp", 44, 90.0, 0.3, 1.8e-3),
]


@pytest.mark.parametrize("pname,code,theta,length,width", CASES,
                         ids=[f"{c[0]}-{int(c[2])}deg-w{c[4]:g}" for c in CASES])
def test_folded_assembled_k_is_psd(pname, code, theta, length, width):
    xyz, conn = fold_mesh(theta, length, width)
    k = assemble(xyz, conn, pname, code)
    ratio, zeros = min_ratio_and_zeros(k)
    elem = element_ratios(xyz, conn, pname, code)
    _, base_zeros = min_ratio_and_zeros(assemble(*fold_mesh(0.0, length, width), pname, code))
    print(f"\n[fold] {pname} theta={theta:g} w={width:g}: assembled min/lmax={ratio:.3e} "
          f"zeros={zeros} (flat baseline {base_zeros}); element-alone min/lmax "
          f"min={min(elem):.3e} max={max(elem):.3e}")
    assert min(elem) >= -NEG_TOL, "control: each element alone must be PSD"
    assert zeros >= 6, "rigid-body modes must remain in the nullspace"
    assert ratio >= -NEG_TOL, f"assembled K indefinite: min/lmax = {ratio:.3e}"


# --------------------------------------------------------------------------------------
# Fold criterion (T5, L18).  Reference: Abaqus, "Defining the initial geometry of
# conventional shell elements", *Calculation of average nodal normals*: the elements at a
# node are grouped by normals within 20 deg of each other (transitively); a group is
# averaged only if all its normals are pairwise within 20 deg, otherwise every element in
# it keeps its own normal.  Without the criterion the IEA-15 TE closure (0.13 m^2 REINF
# quad + 5.4e-4 m^2 TE_FLAT sliver at 80-92 deg) hands the sliver a director 79-85 deg off
# its own normal, and the blade K becomes indefinite (SLURM 11622060).

FOLD_DEG = 20.0
THICK = Laminate([Ply(_MAT, 0.004, 0.0), Ply(_MAT, 0.004, 45.0), Ply(_MAT, 0.004, -45.0)], 0.75)


def strip_mesh(strips, length=0.3, nx=2):
    """Strips sharing the edge s = 0 along x; strip k = (phi_deg, width) runs from the edge
    along d = (0, cos phi, sin phi) with normal x × d, so two strips are |dphi| apart.  A
    negative width reverses the strip (normal -(x × d)): ``(180, -w)`` is the panel that a
    strip at phi = theta folds away from by theta, with consistently oriented normals."""
    xs = np.linspace(0.0, length, nx + 1)
    xyz = [[x, 0.0, 0.0] for x in xs]
    conn = []
    for phi, width in strips:
        d = np.array([0.0, np.cos(np.radians(phi)), np.sin(np.radians(phi))])
        base = len(xyz)
        xyz += [[x, *(abs(width) * d[1:])] for x in xs]
        quads = [[i, i + 1, base + i + 1, base + i] for i in range(nx)]
        conn += [q[::-1] for q in quads] if width < 0 else quads
    return np.array(xyz), conn


def scattered_element_alone(xyz, conn, set_name, type_code):
    k = np.zeros((6 * len(xyz), 6 * len(xyz)))
    for c in conn:
        dofs = np.array([6 * n + d for n in c for d in range(6)])
        k[np.ix_(dofs, dofs)] += assemble(xyz[c], [[0, 1, 2, 3]], set_name, type_code)
    return k


def director_effect(strips, pname, code):
    """||K_assembled - sum of element-alone K|| / ||K_assembled||: zero iff every element
    keeps its own director (a one-element model has nothing to average with)."""
    xyz, conn = strip_mesh(strips)
    k = assemble(xyz, conn, pname, code)
    return float(np.abs(k - scattered_element_alone(xyz, conn, pname, code)).max() / np.abs(k).max())


SEPARATE = [  # every element keeps its own normal
    ("fold-21deg", [(180.0, -0.45), (21.0, 0.45)]),
    ("fold-92deg", [(180.0, -0.45), (92.0, 1.8e-3)]),
    ("fold-150deg", [(180.0, -0.45), (150.0, 0.45)]),
    # chain 0/15/30 deg: one group by transitivity, but 0 vs 30 deg exceeds 20 deg
    ("chain-0-15-30", [(0.0, 0.45), (15.0, 0.45), (30.0, 0.45)]),
]
AVERAGED = [  # smooth: one averaged director (ADR-4 option A unchanged)
    ("fold-10deg", [(180.0, -0.45), (10.0, 0.45)]),
    ("fold-19deg", [(180.0, -0.45), (19.0, 0.45)]),
    ("chain-0-10-19", [(0.0, 0.45), (10.0, 0.45), (19.0, 0.45)]),
]


@pytest.mark.parametrize("pname,code", [("iso", 4), ("comp", 44)])
@pytest.mark.parametrize("name,strips", SEPARATE, ids=[s[0] for s in SEPARATE])
def test_fold_beyond_20deg_keeps_element_directors(name, strips, pname, code):
    eff = director_effect(strips, pname, code)
    print(f"\n[fold-criterion] {pname} {name}: |K - sum K_alone|/|K| = {eff:.2e}")
    assert eff < 1e-12, f"director averaged across a fold > {FOLD_DEG:g} deg: {eff:.2e}"


@pytest.mark.parametrize("pname,code", [("iso", 4), ("comp", 44)])
@pytest.mark.parametrize("name,strips", AVERAGED, ids=[s[0] for s in AVERAGED])
def test_smooth_fold_keeps_averaged_directors(name, strips, pname, code):
    eff = director_effect(strips, pname, code)
    print(f"\n[fold-criterion] {pname} {name}: |K - sum K_alone|/|K| = {eff:.2e}")
    assert eff > 1e-6, f"smooth-surface director averaging lost: {eff:.2e}"


# Thick laminate (h = 12 mm) on a 1.8 mm sliver: before the criterion these folds gave
# min/lmax down to -5e-7 (logs/check/l3-te-fold-sweep.log).
@pytest.mark.parametrize("theta", [30.0, 35.0, 40.0, 45.0, 145.0, 160.0])
def test_thick_sliver_fold_is_psd(theta, monkeypatch):
    monkeypatch.setitem(PROPS, "thick", THICK)
    xyz, conn = strip_mesh([(180.0, -0.45), (theta, 1.8e-3)])
    ratio, zeros = min_ratio_and_zeros(assemble(xyz, conn, "thick", 44))
    print(f"\n[thick-sliver] theta={theta:g}: min/lmax={ratio:.3e} zeros={zeros}")
    assert zeros >= 6
    assert ratio >= -NEG_TOL, f"assembled K indefinite: min/lmax = {ratio:.3e}"
