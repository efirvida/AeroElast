"""3D layer-wise solid judge for the bend-twist coupon (WU5 of the issue-#9 ODD).

Answers hypothesis **H6**: is CLT/FSDT itself insufficient for this strip, i.e. does
a real 3D section distortion make the shell answer differ from 3D elasticity?

The model is a hand-authored CalculiX deck (never ``write_ccx_mesh``): one C3D20R
element per ply, each layer carrying its own ``*ORIENTATION`` ply angle, a clamped
``x = 0`` face and a loaded free tip. It reads the **same metric** as the shell -
the least-squares slope of ``w`` against ``y`` over the tip edge - but on the
mid-plane nodes.

Two traps this file pins down, because both produced a wrong answer first:

1. **A uniform force per node over a brick face is not a uniform traction.** It
   over-weights the perimeter nodes and injects a spurious through-thickness
   moment, which moved the 3D twist ~21% away from the shell - a thickness- and
   mesh-independent gap that looks exactly like a theory verdict and is not one.
   ``consistent_traction_loads`` integrates the 8-node serendipity face shape
   functions against a constant traction instead, which is the honest 3D
   counterpart of a shell's mid-surface line load.
2. **A structured C3D20R grid has orphan node positions.** Nodes at an odd index in
   both in-plane directions are the would-be centre of a 20-node brick's face; no
   element owns them, so CalculiX reports no result for them and a naive sampler
   dies on the missing key. Sampling only element-corner stations avoids them.

Result (see ``test_3d_matches_the_shell_within_2_percent``): with a consistent load
the 3D model reproduces the composite MITC4 shell to ~0.2%, so H6 is refuted at
coupon scale and the CLT reference, not the theory, was the original problem.
"""

from __future__ import annotations

import subprocess
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("_aeroelast", reason="Rust backend not available")

import test_laminate_bend_twist as coupon
from conftest import ccx_bin_or_skip

#: Wall-clock bound for one CalculiX run [s]; a hung solver must fail fast.
_CCX_TIMEOUT_S = 900

# The same orthotropic ply as the coupon suite, now with the full 3D constants the
# shell's CLT cannot see (E3, nu13, nu23, G13, G23).
_PLY_3D = dict(E1=120e9, E2=10e9, E3=10e9, NU12=0.3, NU13=0.3, NU23=0.3,
               G12=5e9, G13=5e9, G23=3e9)


def serendipity8(xi: float, eta: float) -> list[float]:
    """Shape functions of the 8-node serendipity quad, corners then mid-sides."""
    return [
        0.25 * (1 - xi) * (1 - eta) * (-xi - eta - 1),
        0.25 * (1 + xi) * (1 - eta) * (xi - eta - 1),
        0.25 * (1 + xi) * (1 + eta) * (xi + eta - 1),
        0.25 * (1 - xi) * (1 + eta) * (-xi + eta - 1),
        0.5 * (1 - xi * xi) * (1 - eta),
        0.5 * (1 + xi) * (1 - eta * eta),
        0.5 * (1 - xi * xi) * (1 + eta),
        0.5 * (1 - xi) * (1 - eta * eta),
    ]


def consistent_traction_loads(nid, nx: int, ny: int, nz: int, ply_t: float,
                              total: float, width: float) -> dict[int, float]:
    """Nodal forces for a uniform traction ``total / (width * h)`` on the ``+x`` face.

    Trap 1 of this module: a uniform force *per node* is not a uniform traction.
    """
    gauss, weights = np.polynomial.legendre.leggauss(3)
    traction = total / (width * (2 * nz * ply_t))
    dy, dz = width / ny, 2 * ply_t
    jacobian = 0.25 * dy * dz
    loads: dict[int, float] = {}
    for layer in range(nz):
        for ey in range(ny):
            j0, k0 = 2 * ey, 2 * layer
            tip_i = 2 * nx  # the loaded face is the last node station, not element index nx
            face = [nid(tip_i, j0, k0), nid(tip_i, j0 + 2, k0),
                    nid(tip_i, j0 + 2, k0 + 2), nid(tip_i, j0, k0 + 2),
                    nid(tip_i, j0 + 1, k0), nid(tip_i, j0 + 2, k0 + 1),
                    nid(tip_i, j0 + 1, k0 + 2), nid(tip_i, j0, k0 + 1)]
            for a in range(3):
                for b in range(3):
                    weight = weights[a] * weights[b] * jacobian * traction
                    for node, shape in zip(face, serendipity8(gauss[a], gauss[b]), strict=True):
                        loads[node] = loads.get(node, 0.0) + weight * shape
    return loads


def write_3d_deck(path, angles: list[float], ply_t: float, nx: int, ny: int,
                  load_mode: str = "consistent", length: float = coupon.L,
                  width: float = coupon.B, total: float = coupon.P_TIP):
    """Write the layer-wise C3D20R coupon deck; return ``(nid, NX, NY, NZ, axes)``."""
    nz = len(angles)
    NX, NY, NZ = 2 * nx + 1, 2 * ny + 1, 2 * nz + 1

    def nid(i: int, j: int, k: int) -> int:
        return k * NY * NX + j * NX + i + 1

    xs = np.linspace(0.0, length, NX)
    ys = np.linspace(-0.5 * width, 0.5 * width, NY)
    zs = np.linspace(-nz * ply_t / 2.0, nz * ply_t / 2.0, NZ)

    def nset(handle, name, ids):
        handle.write(f"*NSET, NSET={name}\n")
        for start in range(0, len(ids), 16):  # CalculiX: at most 16 entries per line
            handle.write(", ".join(str(c) for c in ids[start:start + 16]) + "\n")

    with open(path, "w", encoding="ascii") as fh:
        fh.write("*HEADING\nbend-twist coupon: 3D layer-wise C3D20R (WU5 judge)\n")
        fh.write("*NODE\n")
        for k in range(NZ):
            for j in range(NY):
                for i in range(NX):
                    fh.write(f"{nid(i,j,k)}, {xs[i]:.10e}, {ys[j]:.10e}, {zs[k]:.10e}\n")
        for layer in range(nz):
            fh.write(f"*ELEMENT, TYPE=C3D20R, ELSET=L{layer + 1}\n")
            for ey in range(ny):
                for ex in range(nx):
                    i0, j0, k0 = 2 * ex, 2 * ey, 2 * layer
                    conn = [nid(i0, j0, k0), nid(i0 + 2, j0, k0),
                            nid(i0 + 2, j0 + 2, k0), nid(i0, j0 + 2, k0),
                            nid(i0, j0, k0 + 2), nid(i0 + 2, j0, k0 + 2),
                            nid(i0 + 2, j0 + 2, k0 + 2), nid(i0, j0 + 2, k0 + 2),
                            nid(i0 + 1, j0, k0), nid(i0 + 2, j0 + 1, k0),
                            nid(i0 + 1, j0 + 2, k0), nid(i0, j0 + 1, k0),
                            nid(i0 + 1, j0, k0 + 2), nid(i0 + 2, j0 + 1, k0 + 2),
                            nid(i0 + 1, j0 + 2, k0 + 2), nid(i0, j0 + 1, k0 + 2),
                            nid(i0, j0, k0 + 1), nid(i0 + 2, j0, k0 + 1),
                            nid(i0 + 2, j0 + 2, k0 + 1), nid(i0, j0 + 2, k0 + 1)]
                    eid = layer * nx * ny + ey * nx + ex + 1
                    tokens = [str(eid)] + [str(c) for c in conn]
                    for start in range(0, len(tokens), 15):
                        chunk = tokens[start:start + 15]
                        fh.write(", ".join(chunk) + ("," if start + 15 < len(tokens) else "") + "\n")
        mid = NZ // 2
        nset(fh, "CLAMP", [nid(0, j, k) for k in range(NZ) for j in range(NY)])
        nset(fh, "MIDPLANE", [nid(i, j, mid) for i in range(NX) for j in range(NY)])
        fh.write("*MATERIAL, NAME=PLY\n*ELASTIC, TYPE=ENGINEERING CONSTANTS\n")
        fh.write(f"{_PLY_3D['E1']}, {_PLY_3D['E2']}, {_PLY_3D['E3']}, {_PLY_3D['NU12']}, "
                 f"{_PLY_3D['NU13']}, {_PLY_3D['NU23']}, {_PLY_3D['G12']}, {_PLY_3D['G13']}\n"
                 f"{_PLY_3D['G23']}, 0.0\n")
        fh.write("*DENSITY\n1500.0\n")
        for layer, angle in enumerate(angles):
            fh.write(f"*ORIENTATION, NAME=O{layer + 1}, SYSTEM=RECTANGULAR\n"
                     "1., 0., 0., 0., 1., 0.\n")
            fh.write(f"3, {float(angle):.4f}\n")
            fh.write(f"*SOLID SECTION, ELSET=L{layer + 1}, MATERIAL=PLY, "
                     f"ORIENTATION=O{layer + 1}\n")
        fh.write("*STEP\n*STATIC\n*BOUNDARY\nCLAMP, 1, 3, 0.0\n*CLOAD\n")
        if load_mode == "consistent":
            loads = consistent_traction_loads(nid, nx, ny, nz, ply_t, total, width)
        elif load_mode == "mid":
            face = [nid(NX - 1, j, mid) for j in range(NY)]
            loads = {nd: total / len(face) for nd in face}
        elif load_mode == "face":
            face = [nid(NX - 1, j, k) for k in range(NZ) for j in range(NY)]
            loads = {nd: total / len(face) for nd in face}
        else:
            raise ValueError(load_mode)
        for node, force in loads.items():
            fh.write(f"{node}, 3, {force:.10e}\n")
        fh.write("*NODE FILE, NSET=MIDPLANE\nU\n*END STEP\n")
    return nid, NX, NY, NZ, (xs, ys, zs)


def parse_all_disp(frd: Path) -> dict[int, list[float]]:
    """Every DISP record in an FRD file.

    ``tests/_ccx_io.py::parse_frd_disp`` keeps only the *last* ``-4 DISP`` block,
    which is fine for small models and silently incomplete for large ones.
    """
    disps: dict[int, list[float]] = {}
    in_disp = False
    for line in open(frd, errors="replace"):
        stripped = line.strip()
        if "-4" in line and "DISP" in line.upper():
            in_disp = True
            continue
        if not in_disp:
            continue
        if stripped.startswith("-3"):
            in_disp = False
            continue
        if not stripped.startswith("-1"):
            continue
        try:
            node = int(line[3:13])
            disps[node] = [float(line[13 + 12 * j:25 + 12 * j]) for j in range(3)]
        except ValueError:
            continue
    return disps


def solid_tip_twist(half: list[float], ply_t: float = coupon.PLY_T, nx: int = 10,
                    ny: int = 6, load_mode: str = "consistent", workdir=None) -> float:
    """Tip twist of the 3D model, read with the coupon's own LSQ metric [rad]."""
    angles = [float(a) for a in half] + [float(a) for a in reversed(half)]
    work = Path(workdir) if workdir else Path("/tmp/ccx3d_coupon")
    work.mkdir(parents=True, exist_ok=True)
    nid, NX, NY, NZ, (_, ys, zs) = write_3d_deck(work / "c.inp", angles, ply_t, nx, ny, load_mode)
    try:
        result = subprocess.run([ccx_bin_or_skip(), "c"], cwd=work,
                                capture_output=True, text=True, timeout=_CCX_TIMEOUT_S)
    except subprocess.TimeoutExpired:
        pytest.fail(f"CalculiX did not finish within {_CCX_TIMEOUT_S}s on the 3D coupon deck")
    if result.returncode != 0:
        pytest.fail(f"CCX failed (rc={result.returncode}):\n{result.stdout[-2000:]}")
    disp = parse_all_disp(work / "c.frd")
    mid = NZ // 2
    w = np.array([disp[nid(NX - 1, j, mid)][2] for j in range(NY)])
    return float(np.polyfit(np.asarray(ys), w, 1)[0])


# ─────────────────────────────────────────────────────────────────────────────
# Deck controls
# ─────────────────────────────────────────────────────────────────────────────


def test_deck_controls(tmp_path):
    """A decoupled layup must not twist, and a +/-45 pair must be antisymmetric.

    The second control is the one that actually validates the per-layer
    ``*ORIENTATION`` handling: an all-zero layup would give zero twist under any
    angle convention.
    """
    decoupled = solid_tip_twist([0.0, 0.0, 0.0, 0.0], workdir=tmp_path / "dec")
    assert abs(decoupled) < 1e-9, f"decoupled 3D control twisted: {decoupled}"

    plus = solid_tip_twist([45.0, 0.0, 0.0, 45.0], workdir=tmp_path / "plus")
    minus = solid_tip_twist([-45.0, 0.0, 0.0, -45.0], workdir=tmp_path / "minus")
    assert plus * minus < 0.0, (plus, minus)
    assert abs(plus + minus) < 0.01 * abs(plus), (plus, minus)


# ─────────────────────────────────────────────────────────────────────────────
# The load-discretisation trap, pinned so it cannot come back
# ─────────────────────────────────────────────────────────────────────────────


def test_a_uniform_nodal_force_on_the_face_is_not_a_uniform_traction(tmp_path):
    """The artifact that first looked like a theory verdict.

    A uniform force per node over the whole tip face over-weights the perimeter
    and moves the 3D twist ~21% away from the shell - the same order as a
    plausible "CLT is wrong" conclusion. A consistent traction removes it. This
    test exists so the difference is measured, not rediscovered.
    """
    lam = coupon._laminate([45.0, 0.0, 0.0, 45.0])
    coords, u, tips = coupon._run_coupon(lam, nx=16, ny=96)
    shell = coupon.twist_lsq(coords, u, tips)

    face = solid_tip_twist([45.0, 0.0, 0.0, 45.0], load_mode="face", workdir=tmp_path / "face")
    consistent = solid_tip_twist([45.0, 0.0, 0.0, 45.0], load_mode="consistent",
                                 workdir=tmp_path / "cons")
    print(f"\nshell theta(L)          = {np.rad2deg(shell):+.6f} deg")
    print(f"3D, nodal force / face  = {np.rad2deg(face):+.6f} deg  ratio {face / shell:.4f}")
    print(f"3D, consistent traction = {np.rad2deg(consistent):+.6f} deg  ratio {consistent / shell:.4f}")

    assert 0.70 < face / shell < 0.85, "the documented face-load artifact moved"
    assert abs(consistent / shell - 1.0) < 0.02


# ─────────────────────────────────────────────────────────────────────────────
# The H6 verdict
# ─────────────────────────────────────────────────────────────────────────────


def test_3d_matches_the_shell_within_2_percent(tmp_path):
    """H6 refuted: 3D layer-wise elasticity reproduces the shell/CLT twist.

    With a consistent load and a matched in-plane refinement, the 3D model and the
    composite MITC4 shell agree to ~0.2%, so CLT/FSDT is not the coupon's problem.
    """
    lam = coupon._laminate([45.0, 0.0, 0.0, 45.0])
    coords, u, tips = coupon._run_coupon(lam, nx=16, ny=96)
    shell = coupon.twist_lsq(coords, u, tips)
    three_d = solid_tip_twist([45.0, 0.0, 0.0, 45.0], nx=16, ny=10,
                              workdir=tmp_path / "fine")
    ratio = three_d / shell
    print(f"\n3D layer-wise = {np.rad2deg(three_d):+.6f} deg, shell = "
          f"{np.rad2deg(shell):+.6f} deg, ratio = {ratio:.4f}")
    assert abs(ratio - 1.0) < 0.02, f"3D/shell = {ratio:.4f}"


def test_3d_mesh_convergence(tmp_path):
    """The 3D answer must settle under in-plane refinement, toward the shell."""
    lam = coupon._laminate([45.0, 0.0, 0.0, 45.0])
    coords, u, tips = coupon._run_coupon(lam, nx=16, ny=96)
    shell = coupon.twist_lsq(coords, u, tips)
    print("\n3D in-plane convergence, consistent traction:")
    values = []
    for nx, ny in ((6, 4), (10, 6), (16, 10)):
        value = solid_tip_twist([45.0, 0.0, 0.0, 45.0], nx=nx, ny=ny,
                                workdir=tmp_path / f"c{nx}")
        values.append(value)
        print(f"  {nx:2d}x{ny:2d}: {np.rad2deg(value):+.6f} deg   ratio 3D/shell = {value / shell:.4f}")
    coarse_gap = abs(values[1] - values[0])
    fine_gap = abs(values[2] - values[1])
    assert fine_gap < coarse_gap, (values, "3D twist is not converging with refinement")
    assert abs(values[2] / shell - 1.0) < 0.02
