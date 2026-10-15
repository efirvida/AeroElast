"""Feasibility probe: can OpenSeesPy act as a second arbiter for the CalculiX-only L2 rows?

Run with an interpreter that has ``openseespy`` (it needs no aeroelast): cases with a CLOSED-FORM
answer, each stated before it is read. Predeclared pass criteria (relative error to the reference):

* ``iso_tip``   isotropic clamped strip, tip load, mean tip w vs Euler-Bernoulli: expected plate
                stiffening of 1-3 % (our CalculiX S8R gives 1.6 % at 32 x 2), accept [0.95, 1.00].
* ``iso_modal`` first bending frequency of the same strip vs the Euler-Bernoulli cantilever:
                accept |err| <= 3 %.
* ``clt_axial`` laminate coupon, edge resultant along x, tip u vs CLT ``1 / (a11 h)`` (a = inverse
                of A), plies [0], [90], quasi-isotropic: accept |err| <= 1 %.
* ``clt_bend``  [0] coupon, tip load, tip w vs Euler-Bernoulli with ``E1 h^3 / 12``: [0.95, 1.00].
* ``clt_coupling`` [0x4/90x4], axial resultant, tip w vs ``kappa_x L^2 / 2`` from the ABD inverse:
                accept |err| <= 5 %.

Plies other than 0 and 90 deg are NOT available: this build lacks ``OrthotropicRotation`` and an
anisotropic ``nDMaterial``, so ``quasi_iso`` is expected to report NotImplementedError.

Usage: python tools/ops_shell_probe.py [--elements ASDShellQ4 ShellMITC4] [--nx 32 --ny 2]
"""

from __future__ import annotations

import argparse
import math

import numpy as np
import openseespy.opensees as ops

E, NU, RHO = 2.1e11, 0.3, 7850.0
L, B, H = 1.0, 0.1, 1.0e-3
# CFRP-like ply, fixed here; the CLT reference is computed from the same numbers.
E1, E2, G12, NU12, PLY_T = 140e9, 10e9, 5e9, 0.3, 0.25e-3


def _grid(nx: int, ny: int, length: float, width: float):
    nodes = {}
    tag = 0
    for i in range(nx + 1):
        for j in range(ny + 1):
            tag += 1
            nodes[(i, j)] = tag
            ops.node(tag, length * i / nx, width * j / ny, 0.0)
    return nodes


def _mesh(element: str, section_tag: int, nx: int, ny: int, length: float, width: float):
    nodes = _grid(nx, ny, length, width)
    etag = 0
    for i in range(nx):
        for j in range(ny):
            etag += 1
            n = [nodes[(i, j)], nodes[(i + 1, j)], nodes[(i + 1, j + 1)], nodes[(i, j + 1)]]
            ops.element(element, etag, *n, section_tag)
    for j in range(ny + 1):
        ops.fix(nodes[(0, j)], 1, 1, 1, 1, 1, 1)
    return nodes


def _fresh():
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)


def _static(nodes, nx, ny, dof, total_force):
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    # consistent for a bilinear edge: 1/2 on the end nodes, 1 on the interior ones (equal spacing)
    w = np.ones(ny + 1)
    w[0] = w[-1] = 0.5
    w /= w.sum()
    for j in range(ny + 1):
        f = [0.0] * 6
        f[dof] = total_force * w[j]
        ops.load(nodes[(nx, j)], *f)
    ops.system("BandGeneral")
    ops.numberer("RCM")
    ops.constraints("Plain")
    ops.test("NormDispIncr", 1e-12, 10)
    ops.algorithm("Linear")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    assert ops.analyze(1) == 0, "static analysis failed"


def _tip_mean(nodes, nx, ny, dof):
    return float(np.mean([ops.nodeDisp(nodes[(nx, j)], dof + 1) for j in range(ny + 1)]))


def case_iso_tip(element, nx, ny):
    _fresh()
    ops.section("ElasticMembranePlateSection", 1, E, NU, H, RHO)
    nodes = _mesh(element, 1, nx, ny, L, B)
    p = 0.05 * E * (B * H**3 / 12) / L**2
    _static(nodes, nx, ny, 2, p)
    return _tip_mean(nodes, nx, ny, 2) / (p * L**3 / (3 * E * B * H**3 / 12))


def case_iso_modal(element, nx, ny):
    _fresh()
    ops.section("ElasticMembranePlateSection", 1, E, NU, H, RHO)
    _mesh(element, 1, nx, ny, L, B)
    lam = ops.eigen(6)
    f = sorted(math.sqrt(x) / (2 * math.pi) for x in lam if x > 0)
    f_eb = (1.8751**2 / (2 * math.pi)) * math.sqrt(E * (B * H**3 / 12) / (RHO * B * H * L**4))
    return f[0] / f_eb, f[:3], f_eb


def clt_abd(angles):
    """A, B, D of a laminate of equal plies (bottom to top), plain numpy CLT."""
    nu21 = NU12 * E2 / E1
    d = 1 - NU12 * nu21
    q = np.array([[E1 / d, NU12 * E2 / d, 0], [NU12 * E2 / d, E2 / d, 0], [0, 0, G12]])
    n = len(angles)
    z = np.linspace(-n * PLY_T / 2, n * PLY_T / 2, n + 1)
    a, bm, dm = np.zeros((3, 3)), np.zeros((3, 3)), np.zeros((3, 3))
    for k, ang in enumerate(angles):
        t = math.radians(ang)
        c, s = math.cos(t), math.sin(t)
        # engineering-shear transformation: Qbar = T^-1 Q R T R^-1
        t_mat = np.array([[c * c, s * s, 2 * s * c], [s * s, c * c, -2 * s * c], [-s * c, s * c, c * c - s * s]])
        r = np.diag([1.0, 1.0, 2.0])
        qb = np.linalg.inv(t_mat) @ q @ r @ t_mat @ np.linalg.inv(r)
        a += qb * (z[k + 1] - z[k])
        bm += qb * (z[k + 1] ** 2 - z[k] ** 2) / 2
        dm += qb * (z[k + 1] ** 3 - z[k] ** 3) / 3
    return a, bm, dm


def layered_section(angles):
    """Cross-ply laminate: this build has no ply-rotation material, so only 0 and 90 deg plies."""
    _fresh()
    nu23, g23 = 0.4, E2 / (2 * 1.4)
    nu21 = NU12 * E2 / E1
    # ElasticOrthotropic(Ex, Ey, Ez, nuxy, nuyz, nuzx, Gxy, Gyz, Gzx) in the global axes
    ops.nDMaterial("ElasticOrthotropic", 100, E1, E2, E2, NU12, nu23, nu21, G12, g23, G12)
    ops.nDMaterial("ElasticOrthotropic", 101, E2, E1, E2, nu21, NU12, nu23, G12, G12, g23)
    layers = []
    for k, ang in enumerate(angles):
        if ang not in (0, 90):
            raise NotImplementedError(f"ply angle {ang}: no rotation material in this OpenSees build")
        ops.nDMaterial("PlateFiber", 300 + k, 100 if ang == 0 else 101)
        layers += [300 + k, PLY_T]
    ops.section("LayeredShell", 2, len(angles), *layers)
    return 2


def case_clt_axial(element, nx, ny, angles):
    layered_section(angles)
    nodes = _mesh(element, 2, nx, ny, L, B)
    n_load = 1.0e4  # N/m resultant
    _static(nodes, nx, ny, 0, n_load * B)
    u_ref = n_load * L * np.linalg.inv(clt_abd(angles)[0])[0, 0]
    return _tip_mean(nodes, nx, ny, 0) / u_ref


def case_clt_bend(element, nx, ny):
    angles = [0] * 8
    layered_section(angles)
    nodes = _mesh(element, 2, nx, ny, L, B)
    h = len(angles) * PLY_T
    p = 0.01 * E1 * B * h**3 / 12 / L**2
    _static(nodes, nx, ny, 2, p)
    return _tip_mean(nodes, nx, ny, 2) / (p * L**3 / (3 * E1 * B * h**3 / 12))


def case_clt_coupling(element, nx, ny):
    """[0x4/90x4]: axial resultant bends the strip; tip w vs 1/2 kappa_x L^2 of the ABD inverse."""
    angles = [0] * 4 + [90] * 4
    layered_section(angles)
    nodes = _mesh(element, 2, nx, ny, L, B)
    n_load = 1.0e4
    _static(nodes, nx, ny, 0, n_load * B)
    a, bm, dm = clt_abd(angles)
    abd = np.block([[a, bm], [bm, dm]])
    kappa_x = np.linalg.inv(abd)[3, 0] * n_load
    return _tip_mean(nodes, nx, ny, 2) / (0.5 * kappa_x * L**2)


def _try(label, fn):
    try:
        print(f"  {label}: {fn()}")
    except Exception as exc:  # report the failure, never hide it
        print(f"  {label}: FAILED {type(exc).__name__}: {exc}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--elements", nargs="*", default=["ASDShellQ4", "ShellMITC4"])
    ap.add_argument("--nx", type=int, default=32)
    ap.add_argument("--ny", type=int, default=2)
    args = ap.parse_args()
    quasi = [0, 45, -45, 90, 90, -45, 45, 0]
    for el in args.elements:
        print(f"== {el}  mesh {args.nx} x {args.ny}")
        for nx, ny in ((args.nx, args.ny), (2 * args.nx, 2 * args.ny)):
            _try(f"iso_tip {nx}x{ny} w/EB", lambda: f"{case_iso_tip(el, nx, ny):.4f}")
        _try("iso_modal f1/f_EB", lambda: f"{case_iso_modal(el, args.nx, args.ny)[0]:.4f}")
        for name, ang in (("uni_0", [0] * 8), ("uni_90", [90] * 8), ("cross_sym", [0, 90, 0, 90, 90, 0, 90, 0]), ("quasi_iso", quasi)):
            _try(f"clt_axial {name} u/CLT", lambda: f"{case_clt_axial(el, args.nx, args.ny, ang):.4f}")
        _try("clt_bend uni_0 w/EB", lambda: f"{case_clt_bend(el, args.nx, args.ny):.4f}")
        _try("clt_coupling asym w/CLT", lambda: f"{case_clt_coupling(el, args.nx, args.ny):.4f}")


if __name__ == "__main__":
    main()
