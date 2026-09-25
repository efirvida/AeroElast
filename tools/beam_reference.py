"""S-2 beam reference — 3D Euler-Bernoulli cantilever from IEA 15 MW props.

Builds a 1D beam model of the IEA 15 MW blade whose GEOMETRY is derived from
the AeroElast shell mesh (section centroids + measured structural twist per
spanwise slice) and whose STIFFNESS/MASS come from the official reference
files (ElastoDyn blade file + BeamDyn 6x6 section matrices).

This isolates shell-vs-beam discretization physics: same geometry, same
section orientation, reference stiffness.

Conventions (from the AeroElast mesh):
  - span direction: +Z (blade root at z=0, tip at z~117 m)
  - flap axis: section thickness direction (minor PCA spread)
  - edge axis: chord direction (major PCA spread)

Usage:
    from tools.beam_reference import BeamReference, load_elastodyn_blade
    ed = load_elastodyn_blade("tests/IEA15MW/reference/IEA-15-240-RWT_ElastoDyn_blade.dat")
    beam = BeamReference.from_mesh(mesh, ed)
    u_tip = beam.solve_tip_load((0.0, 1.0, 0.0), 1e6)
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Tuple

import numpy as np


# ---------------------------------------------------------------------------
# Reference file parsers
# ---------------------------------------------------------------------------

@dataclass
class ElastoDynBlade:
    """Distributed properties from an ElastoDyn 6-column blade file."""
    bl_fract: np.ndarray      # station fraction of blade length (-)
    z: np.ndarray             # station spanwise coordinate (m)
    twist_deg: np.ndarray     # structural twist (deg)
    mass_dens: np.ndarray     # kg/m
    ei_flap: np.ndarray       # N m^2
    ei_edge: np.ndarray       # N m^2


@dataclass
class BeamDynBlade:
    """Section matrices from a BeamDyn blade file.

    File matrix layout (verified against ElastoDyn values):
      1,2 = shear, 3 = axial EA, 4 = EI_edge, 5 = EI_flap, 6 = GJ.
    """
    z: np.ndarray
    ea: np.ndarray
    gj: np.ndarray
    ei_flap: np.ndarray
    ei_edge: np.ndarray
    mass_dens: np.ndarray


def load_elastodyn_blade(path: str | Path, blade_length: float = 117.0) -> ElastoDynBlade:
    """Parse the 6-column ElastoDyn blade file (50 stations)."""
    rows = []
    with open(path) as f:
        for line in f:
            tokens = line.split()
            if len(tokens) == 6:
                try:
                    rows.append([float(t) for t in tokens])
                except ValueError:
                    continue
    arr = np.asarray(rows, dtype=np.float64)
    if arr.shape[0] < 10:
        raise ValueError(f"Expected >=10 property rows in {path}, got {arr.shape[0]}")
    bl_fract = arr[:, 0]
    return ElastoDynBlade(
        bl_fract=bl_fract,
        z=bl_fract * blade_length,
        twist_deg=arr[:, 2],
        mass_dens=arr[:, 3],
        ei_flap=arr[:, 4],
        ei_edge=arr[:, 5],
    )


def load_beamdyn_blade(path: str | Path, blade_length: float = 117.0) -> BeamDynBlade:
    """Parse a BeamDyn blade file: 26 stations with 6x6 K and M matrices."""
    lines = [l.strip() for l in open(path)]
    i = 0
    while i < len(lines) and "DISTRIBUTED PROPERTIES" not in lines[i]:
        i += 1
    i += 1
    zs, ea, gj, eif, eie, md = [], [], [], [], [], []
    while i < len(lines):
        parts = lines[i].split()
        if len(parts) != 1:
            i += 1
            continue
        try:
            s = float(parts[0])
        except ValueError:
            i += 1
            continue
        j = i + 1
        while j < len(lines) and not lines[j]:
            j += 1
        K = np.array([[float(x) for x in lines[j + k].split()] for k in range(6)])
        j += 6
        while j < len(lines) and not lines[j]:
            j += 1
        M = np.array([[float(x) for x in lines[j + k].split()] for k in range(6)])
        zs.append(s * blade_length)
        ea.append(K[2, 2])
        eie.append(K[3, 3])
        eif.append(K[4, 4])
        gj.append(K[5, 5])
        md.append(M[0, 0])
        i = j + 6
    if not zs:
        raise ValueError(f"No stations parsed from {path}")
    return BeamDynBlade(
        z=np.asarray(zs),
        ea=np.asarray(ea),
        gj=np.asarray(gj),
        ei_flap=np.asarray(eif),
        ei_edge=np.asarray(eie),
        mass_dens=np.asarray(md),
    )


def _interp(x_src: np.ndarray, y_src: np.ndarray, x_dst: np.ndarray) -> np.ndarray:
    return np.interp(x_dst, x_src, y_src)


# ---------------------------------------------------------------------------
# Section orientation from the AeroElast mesh
# ---------------------------------------------------------------------------

def extract_section_orientation(
    coords: np.ndarray, n_slices: int = 60
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Slice the blade mesh along Z and return (z, centroid, chord_dir, flap_dir).

    Returns
    -------
    z : (n_slices,) slice midpoints
    centroids : (n_slices, 3)
    chord_dir : (n_slices, 3) unit vectors in the chord (edge) direction
    flap_dir : (n_slices, 3) unit vectors in the flap (thickness) direction
    """
    z_all = coords[:, 2]
    z_min, z_max = float(z_all.min()), float(z_all.max())
    edges = np.linspace(z_min, z_max, n_slices + 1)
    z_mid, cents, chord, flap = [], [], [], []
    for a, b in zip(edges[:-1], edges[1:]):
        mask = (z_all >= a) & (z_all < b) if b < z_max else (z_all >= a) & (z_all <= b)
        pts = coords[mask]
        if pts.shape[0] < 3:
            continue
        cent = pts.mean(axis=0)
        d = pts - cent
        cov = d.T @ d / pts.shape[0]
        evals, evecs = np.linalg.eigh(cov)
        # major axis (chord) = largest eigenvalue; minor (flap) = smallest
        ch = evecs[:, 2]
        fl = evecs[:, 0]
        # make chord/flap a right-handed frame with local tangent +Z
        tang = np.array([0.0, 0.0, 1.0])
        ch = ch - (ch @ tang) * tang
        ch /= np.linalg.norm(ch)
        fl = fl - (fl @ tang) * tang
        fl /= np.linalg.norm(fl)
        # enforce continuity of direction along the blade (avoid per-slice flips)
        if chord and np.dot(chord[-1], ch) < 0:
            ch = -ch
        if flap and np.dot(flap[-1], fl) < 0:
            fl = -fl
        z_mid.append(0.5 * (a + b))
        cents.append(cent)
        chord.append(ch)
        flap.append(fl)
    return (
        np.asarray(z_mid),
        np.asarray(cents),
        np.asarray(chord),
        np.asarray(flap),
    )


# ---------------------------------------------------------------------------
# 3D Euler-Bernoulli beam element (12 DOF)
# ---------------------------------------------------------------------------

def _beam_local_k(ea: float, eif: float, eie: float, gj: float, L: float) -> np.ndarray:
    """12x12 local stiffness. Local axes: x=tangent, y=chord(edge), z=flap.

    DOF order per node: [ux, uy, uz, rx, ry, rz].
    Flapwise bending (deflection uz, rotation about y) uses EI_flap;
    edgewise bending (deflection uy, rotation about z) uses EI_edge.
    """
    k = np.zeros((12, 12))
    # axial
    k[0, 0] = k[6, 6] = ea / L
    k[0, 6] = k[6, 0] = -ea / L
    # torsion
    k[3, 3] = k[9, 9] = gj / L
    k[3, 9] = k[9, 3] = -gj / L
    # bending in xy plane (deflection uy, rotation rz) — EI_edge
    a, b, c, d = 12 * eie / L**3, 6 * eie / L**2, 4 * eie / L, 2 * eie / L
    k[1, 1] = k[7, 7] = a
    k[1, 7] = k[7, 1] = -a
    k[1, 5] = k[5, 1] = k[1, 11] = k[11, 1] = b
    k[7, 5] = k[5, 7] = k[7, 11] = k[11, 7] = -b
    k[5, 5] = k[11, 11] = c
    k[5, 11] = k[11, 5] = d
    # bending in xz plane (deflection uz, rotation ry) — EI_flap
    a, b, c, d = 12 * eif / L**3, 6 * eif / L**2, 4 * eif / L, 2 * eif / L
    k[2, 2] = k[8, 8] = a
    k[2, 8] = k[8, 2] = -a
    k[2, 4] = k[4, 2] = k[2, 10] = k[10, 2] = -b
    k[8, 4] = k[4, 8] = k[8, 10] = k[10, 8] = b
    k[4, 4] = k[10, 10] = c
    k[4, 10] = k[10, 4] = d
    return k


class BeamReference:
    """3D Euler-Bernoulli beam model of the blade, mesh-derived geometry."""

    def __init__(
        self,
        z: np.ndarray,
        centroids: np.ndarray,
        chord_dir: np.ndarray,
        flap_dir: np.ndarray,
        ed: ElastoDynBlade,
        bd: Optional[BeamDynBlade] = None,
        g: Tuple[float, float, float] = (0.0, 0.0, -9.81),
    ):
        self.z = z
        self.centroids = centroids
        self.chord_dir = chord_dir
        self.flap_dir = flap_dir
        self.ed = ed
        self.g = np.asarray(g, dtype=np.float64)
        # Interpolate the reference data by SPAN FRACTION.  The mesh may live
        # in the rotor frame (BladeMesh translates the blade by the hub radius,
        # matching the official ElastoDyn HubRad 3.97 / TipRad 120.97) while
        # the ElastoDyn/BeamDyn blade files are blade-local.  Using absolute z
        # silently mixes the two datums and shifts every property by the hub
        # radius.  Note the fraction must be (x - x.min())/(x.max() - x.min()):
        # when the origins differ, x/x.max() is NOT the span fraction.
        def _frac(x):
            x = np.asarray(x, dtype=np.float64)
            return (x - x.min()) / (x.max() - x.min())

        z_n = _frac(z)
        ed_n = _frac(ed.z)
        self.ei_flap = _interp(ed_n, ed.ei_flap, z_n)
        self.ei_edge = _interp(ed_n, ed.ei_edge, z_n)
        self.mass_dens = _interp(ed_n, ed.mass_dens, z_n)
        if bd is not None:
            bd_n = _frac(bd.z)
            self.gj = _interp(bd_n, bd.gj, z_n)
            self.ea = _interp(bd_n, bd.ea, z_n)
        else:
            self.gj = np.full_like(z, 1e12)
            self.ea = np.full_like(z, 1e12)
        self._assemble()

    @classmethod
    def from_mesh(
        cls, mesh, ed: ElastoDynBlade, bd: Optional[BeamDynBlade] = None,
        n_slices: int = 60,
    ) -> "BeamReference":
        coords = mesh.coords_array
        z, cents, chord, flap = extract_section_orientation(coords, n_slices)
        return cls(z, cents, chord, flap, ed, bd)

    # ------------------------------------------------------------------ FEM --
    def _assemble(self) -> None:
        n_nodes = len(self.z)
        n_dof = 6 * n_nodes
        K = np.zeros((n_dof, n_dof))
        for e in range(n_nodes - 1):
            L = self.z[e + 1] - self.z[e]
            if L <= 0:
                continue
            tang = self.centroids[e + 1] - self.centroids[e]
            tang = tang / np.linalg.norm(tang)
            ch = self.chord_dir[e] + self.chord_dir[e + 1]
            fl = self.flap_dir[e] + self.flap_dir[e + 1]
            # orthonormalize against tangent
            ch = ch - (ch @ tang) * tang
            ch /= np.linalg.norm(ch)
            fl = fl - (fl @ tang) * tang - (fl @ ch) * ch
            fl /= np.linalg.norm(fl)
            Rm = np.column_stack([tang, ch, fl])  # columns = local axes in global
            ei_f = 0.5 * (self.ei_flap[e] + self.ei_flap[e + 1])
            ei_e = 0.5 * (self.ei_edge[e] + self.ei_edge[e + 1])
            gj = 0.5 * (self.gj[e] + self.gj[e + 1])
            ea = 0.5 * (self.ea[e] + self.ea[e + 1])
            k_loc = _beam_local_k(ea, ei_f, ei_e, gj, L)
            T = np.zeros((12, 12))
            T[0:3, 0:3] = Rm
            T[3:6, 3:6] = Rm
            T[6:9, 6:9] = Rm
            T[9:12, 9:12] = Rm
            k_glob = T @ k_loc @ T.T
            idx = np.array([6 * e + d for d in range(12)])
            K[np.ix_(idx, idx)] += k_glob
        self.K = K
        self.n_nodes = n_nodes

    def _solve(self, F: np.ndarray) -> np.ndarray:
        free = np.arange(6, 6 * self.n_nodes)
        K_ff = self.K[np.ix_(free, free)]
        u = np.zeros(6 * self.n_nodes)
        u[free] = np.linalg.solve(K_ff, F[free])
        return u

    def tip_dofs(self) -> slice:
        return slice(-6, None)

    # ------------------------------------------------------------- load cases --
    def solve_tip_load(self, force: Tuple[float, float, float]) -> np.ndarray:
        """Tip nodal force [Fx, Fy, Fz] in global coords."""
        F = np.zeros(6 * self.n_nodes)
        F[-6:-3] = force
        return self._solve(F)

    def solve_distributed_load(self, q: Tuple[float, float, float]) -> np.ndarray:
        """Uniform distributed load q [N/m] in global coords (consistent nodal)."""
        q = np.asarray(q, dtype=np.float64)
        F = np.zeros(6 * self.n_nodes)
        for e in range(self.n_nodes - 1):
            L = self.z[e + 1] - self.z[e]
            if L <= 0:
                continue
            w = 0.5 * L * q
            F[6 * e : 6 * e + 3] += w
            F[6 * e + 6 : 6 * e + 9] += w
        return self._solve(F)

    def solve_gravity(
        self, g: Optional[Tuple[float, float, float]] = None
    ) -> np.ndarray:
        """Self-weight consistent nodal loads from BMassDen.

        g defaults to the acceleration passed at construction.
        """
        g = np.asarray(self.g if g is None else g, dtype=np.float64)
        F = np.zeros(6 * self.n_nodes)
        for e in range(self.n_nodes - 1):
            L = self.z[e + 1] - self.z[e]
            if L <= 0:
                continue
            m = 0.5 * (self.mass_dens[e] + self.mass_dens[e + 1])
            F[6 * e : 6 * e + 3] += 0.5 * L * m * g
            F[6 * e + 6 : 6 * e + 9] += 0.5 * L * m * g
        return self._solve(F)

    @property
    def total_mass_kg(self) -> float:
        m = 0.0
        for e in range(self.n_nodes - 1):
            L = self.z[e + 1] - self.z[e]
            m += 0.5 * (self.mass_dens[e] + self.mass_dens[e + 1]) * L
        return float(m)


# ---------------------------------------------------------------------------
# Self-test: uniform isotropic cantilever, tip load, analytic solution
# ---------------------------------------------------------------------------

def _self_test() -> None:
    L, EI, EA, GJ = 10.0, 1.0e8, 1.0e9, 1.0e8
    n = 21
    z = np.linspace(0.0, L, n)
    cents = np.column_stack([np.zeros(n), np.zeros(n), z])
    chord = np.tile([1.0, 0.0, 0.0], (n, 1))
    flap = np.tile([0.0, 1.0, 0.0], (n, 1))
    ed = ElastoDynBlade(
        bl_fract=z / L, z=z, twist_deg=np.zeros(n),
        mass_dens=np.zeros(n), ei_flap=np.full(n, EI), ei_edge=np.full(n, EI),
    )
    bd = BeamDynBlade(
        z=z, ea=np.full(n, EA), gj=np.full(n, GJ),
        ei_flap=np.full(n, EI), ei_edge=np.full(n, EI), mass_dens=np.zeros(n),
    )
    beam = BeamReference(z, cents, chord, flap, ed, bd)
    u = beam.solve_tip_load((0.0, 1.0, 0.0))  # Fy = 1 N
    P = 1.0
    delta_analytic = P * L**3 / (3 * EI)
    delta_fem = float(u[-5])
    rel = abs(delta_fem - delta_analytic) / delta_analytic
    print(f"[beam self-test] tip Fy=1N: FEM={delta_fem:.6e} analytic={delta_analytic:.6e} rel={rel:.2%}")
    assert rel < 0.01, f"beam self-test failed: rel err {rel:.2%}"


if __name__ == "__main__":
    _self_test()
    print("OK")
