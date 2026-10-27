"""S-1 — sectional property extractor for shell blade meshes.

Computes spanwise-distributed sectional properties (EA, GJ, EI_flap,
EI_edge, mass per unit length, section CG) directly from the ASSEMBLED
stiffness and mass of a shell blade model.

Method (Lagrange-multiplier homogenization, BECAS-style):

  For each spanwise slice [z0, z1]:
    1. The slice stiffness K_slice is extracted from the assembled K.
    2. Six rigid section-strain states (3 translations + 3 rotations) are
       imposed on ALL slice nodes via Lagrange multipliers on the free slice.
    3. The multipliers for unit strains ARE the section forces, so the
       resulting 6x6 stiffness falls out column by column.  Interior slice
       nodes stay free and absorb section warping.
    4. The 6x6 is rotated into the section frame (tangent/chord/flap
       measured from the mesh via PCA) and the diagonal entries reported.

The section frame is measured from the mesh itself (chord = major PCA
spread direction, flap = minor), so no twist sign convention is assumed.

Note: properties are reported about the slice node centroid.  The reference
BeamDyn/ElastoDyn tables use the elastic/mass centers; offsets are cm-scale
except near the root, so centroid-based values are within a few % of the
beam targets for EI comparisons.

Usage:
    from aeroelast.postprocess.sectional import SectionalExtractor
    ext = SectionalExtractor(mesh, properties)
    props = ext.extract()               # mass per station (station membership)
    z, ei_flap, ei_edge = ext.ei_from_static_response(properties)

STATUS (2026-09-09): the MASS path (station membership via the detailed
element sets) and the EI path (static tip-load curvature method) are
validated against the BeamDyn targets (mid-span: EI_edge ~10%, EI_flap
~20-25% systematically high from global-load direction mixing, m(r) ~1-10%).
The slice-based stiffness homogenization (section_stiffness) remains
experimental: the prescribed-rigid-motion reactions carry fixed-end effects
and axis-mixing biases that the curvature method avoids.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Optional, Tuple

import numpy as np
from scipy.interpolate import UnivariateSpline
from scipy.sparse import coo_matrix, csr_matrix
from scipy.sparse.linalg import spsolve

from aeroelast.core.mesh import MeshModel


@dataclass
class SectionProperties:
    z: np.ndarray          # slice mid-span coordinate (m)
    ea: np.ndarray         # axial stiffness (N)
    gj: np.ndarray         # torsional stiffness (N m^2)
    ei_flap: np.ndarray    # flapwise bending stiffness (N m^2)
    ei_edge: np.ndarray    # edgewise bending stiffness (N m^2)
    mass_dens: np.ndarray  # mass per unit length (kg/m)
    cg: np.ndarray         # (n, 3) section centroids (m)


def _dofs_of_nodes(node_indices: np.ndarray, dpn: int) -> np.ndarray:
    return np.concatenate([node_indices * dpn + d for d in range(dpn)])


class SectionalExtractor:
    """Sectional property extraction from an assembled shell blade model."""

    def __init__(self, mesh: MeshModel, properties: Optional[dict] = None,
                 assembler=None, span_direction=(0.0, 0.0, 1.0)):
        from aeroelast.core.assembler import MeshAssembler
        from aeroelast.elements import ElementFamily

        self.mesh = mesh
        if assembler is not None:
            self._assembler = assembler
        else:
            model = {
                "elements": {
                    "element_family": ElementFamily.SHELL,
                    "span_direction": tuple(span_direction),
                    "properties": properties,
                },
            }
            self._assembler = MeshAssembler(mesh, model)
        assembler = self._assembler
        self.dpn = int(getattr(assembler, "dofs_per_node", 6))
        n_dof = assembler.dofs_count

        k_rows, k_cols, k_vals = assembler._rust.assemble_k()
        self._k_rows = np.asarray(k_rows, dtype=np.int64)
        self._k_cols = np.asarray(k_cols, dtype=np.int64)
        self._k_vals = np.asarray(k_vals, dtype=np.float64)
        self._n_dof = n_dof
        self.K: csr_matrix = coo_matrix(
            (self._k_vals, (self._k_rows, self._k_cols)),
            shape=(n_dof, n_dof),
        ).tocsr()

        m_rows, m_cols, m_vals = assembler._rust.assemble_m()
        M = coo_matrix(
            (np.asarray(m_vals), (np.asarray(m_rows), np.asarray(m_cols))),
            shape=(n_dof, n_dof),
        ).tocsr()
        # lumped nodal translational mass: row-sum of the x-dof rows
        n_nodes = mesh.node_count
        self.node_mass = np.zeros(n_nodes)
        for i in range(n_nodes):
            self.node_mass[i] = float(M.getrow(i * self.dpn).sum())

        self.coords = np.asarray(mesh.coords_array, dtype=np.float64)
        self._node_normals = self._compute_node_normals()

    def _compute_node_normals(self) -> np.ndarray:
        """Per-node shell normals (area-weighted element normals)."""
        normals = np.zeros((self.mesh.node_count, 3))
        for elem in self.mesh.elements:
            ids = elem.node_ids
            p0 = self.coords[ids[0]]
            p1 = self.coords[ids[1]]
            p2 = self.coords[ids[2]]
            n = np.cross(p1 - p0, p2 - p0)
            for i in ids:
                normals[i] += n
        lens = np.linalg.norm(normals, axis=1, keepdims=True)
        lens[lens < 1e-12] = 1.0
        return normals / lens

    def _submatrix(self, dofs: np.ndarray) -> "csr_matrix":
        """Slice the assembled K to a DOF subset via direct COO filtering.

        Avoids scipy csr fancy indexing entirely (the raw Rust COO carries
        duplicate entries whose summed structure is not safely indexable).
        """
        keep = np.isin(self._k_rows, dofs) & np.isin(self._k_cols, dofs)
        remap = {int(d): k for k, d in enumerate(dofs)}
        rr = np.array([remap[int(r)] for r in self._k_rows[keep]], dtype=np.int64)
        cc = np.array([remap[int(c)] for c in self._k_cols[keep]], dtype=np.int64)
        vv = self._k_vals[keep]
        return coo_matrix((vv, (rr, cc)), shape=(dofs.size, dofs.size)).tocsc()

    # ------------------------------------------------------------------ utils --    # ------------------------------------------------------------------ utils --
    def station_membership(self) -> Dict[int, Dict[str, np.ndarray]]:
        """Classify nodes by blade station via the detailed element sets
        named ``{station:02d}_{seg:02d}_{surface}_{component}``.

        Returns {station: {"nodes": idx array, "z": mean station z}}.
        """
        nid_to_idx = self.mesh.node_id_to_index
        stations: Dict[int, set] = {}
        for name, eset in self.mesh.element_sets.items():
            st = _station_of(name)
            if st is None:
                continue
            if st not in stations:
                stations[st] = set()
            for e in eset.elements:
                for nid in e.node_ids:
                    stations[st].add(nid_to_idx[nid])
        out: Dict[int, Dict[str, np.ndarray]] = {}
        for st, nodes in stations.items():
            arr = np.array(sorted(nodes), dtype=np.int64)
            out[st] = {"nodes": arr, "z": float(self.coords[arr, 2].mean())}
        return out

    def _station_slice(self, st: int, membership: Dict[int, Dict[str, np.ndarray]]
                       ) -> Dict[str, np.ndarray]:
        """Face/interior classification for one blade station.

        Faces are the node rings shared with the neighbouring stations
        (twisted sections have non-planar rings, so z-band classification is
        unreliable); the first station's lower face is the mesh RootNodes set.
        """
        nodes = membership[st]["nodes"]
        if st - 1 in membership:
            lower = np.intersect1d(nodes, membership[st - 1]["nodes"])
        else:
            root_set = self.mesh.get_node_set("RootNodes")
            root = np.array(sorted(root_set.nodes.keys()), dtype=np.int64)
            lower = np.intersect1d(nodes, root)
        if st + 1 in membership:
            upper = np.intersect1d(nodes, membership[st + 1]["nodes"])
        else:
            upper = np.array([], dtype=np.int64)
        interior = np.setdiff1d(nodes, np.union1d(lower, upper))
        centroid = self.coords[nodes].mean(axis=0)
        z_lo = self.coords[lower, 2].mean() if lower.size else self.coords[nodes, 2].min()
        z_hi = self.coords[upper, 2].mean() if upper.size else self.coords[nodes, 2].max()
        return {"lower": lower, "upper": upper, "interior": interior,
                "centroid": centroid, "z0": z_lo, "z1": z_hi}

    def _section_axes(self, nodes: np.ndarray):
        """Return (tangent, chord, flap) unit axes for a station's nodes."""
        pts = self.coords[nodes]
        d = pts - pts.mean(axis=0)
        _, evecs = np.linalg.eigh(d.T @ d / max(pts.shape[0], 1))
        chord = evecs[:, 2]
        flap = evecs[:, 0]
        tangent = np.array([0.0, 0.0, 1.0])
        chord -= (chord @ tangent) * tangent
        chord /= np.linalg.norm(chord)
        flap -= (flap @ tangent) * tangent
        flap -= (flap @ chord) * chord
        flap /= np.linalg.norm(flap)
        return tangent, chord, flap

    def _projected_rots(self, nodes: np.ndarray, axis: np.ndarray,
                        theta: float) -> np.ndarray:
        """Section rotation theta·axis projected onto each node's shell plane
        (removes the drilling component)."""
        n = self._node_normals[nodes]
        rot = theta * np.asarray(axis, dtype=np.float64)
        rot = rot[None, :] - (np.sum(rot * n, axis=1)[:, None]) * n
        return rot

    # ------------------------------------------------------------- extraction --
    def ei_from_static_response(
        self, properties: Optional[dict] = None, load: float = 1e6,
        frame: str = "global",
        stations: Optional[np.ndarray] = None,
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Extract EI_flap(r) and EI_edge(r) from full-mesh static tip-load
        solves (classical curvature method: EI = M(r)/kappa(r)).

        Runs two linear static solves (tip load in +Y for flapwise, +X for
        edgewise, distributed over the tip-section nodes) and derives the
        sectional bending stiffness from the smoothed curvature of the
        section-mean transverse deflection.  Returns (z, EI_flap, EI_edge).

        Note: the load direction is global, so the extracted EI in each
        direction mixes flap/edge through the blade twist exactly like the
        BeamDyn reference convention (which is also about global-pitch axes).

        ``frame`` selects the frame in which moment and curvature are read:

        ``"global"`` (default, the historical path, unchanged): curvature is
        the second derivative of the section-mean global X/Y deflection.  On
        a PREBENT, torsionally soft blade the global-Y moment of the edge
        case has a torsion component about the tilted beam axis, so the
        reading is 1/EI_g = cos^2(phi)/EI_loc + sin^2(phi)/GJ, which biases
        EI_edge low (about 3.6-6 pp on the IEA-15 blade, odd/tasks/
        l3-physical-validation.md L27/L28).

        ``"tangent"``: the BeamDyn 6x6 frame.  Per station, the beam-axis
        tangent t is the slope of the station-centroid line, the edge axis
        is global X projected normal to t and the flap axis is t x edge.
        The rigid rotation vector theta(z) of each station (least-squares
        fit of u = u0 + theta x (r - c) to its nodal displacements) is
        differentiated with a local quadratic fit over +-2 station spacings,
        kappa = d theta/ds = (d theta/dz) t_z.  The moment of the tip force
        about the station centroid, M = (c_tip - c) x F, is resolved with
        kappa in the same local axes and only the bending components are
        used: EI_edge = |M.flap| / |kappa.flap| (X load) and EI_flap =
        |M.edge| / |kappa.edge| (Y load).  Torsion enters neither M nor
        kappa component, so it cannot soften the bending reading; bend-twist
        coupling of the section, if any, is not removed.

        ``stations``: reference span coordinates [m]. When given, the curvature is fitted and
        returned only at the mesh stations nearest to them (``select_stations``), so the result
        does not depend on how many stations the mesher inserted (``BladeMesh.max_bay_length``;
        the spline smoothing and the local-fit window both scale with the station count). The
        static solves still run on the full mesh. ``None`` keeps every station.
        """
        if frame not in ("global", "tangent"):
            raise ValueError(f"frame must be 'global' or 'tangent', got {frame!r}")
        from aeroelast.core.bc import DirichletCondition, NodalLoad
        from aeroelast.elements import ElementFamily
        from aeroelast.solvers.elasticity.static_linear import StaticLinearSolver

        model = {
            "solver": {},
            "elements": {
                "element_family": ElementFamily.SHELL,
                "span_direction": (0.0, 0.0, 1.0),
                "properties": properties,
            },
        }
        membership = self.station_membership()
        sts = sorted(membership)
        zs = np.array([membership[st]["z"] for st in sts])
        L = zs.max()
        tip_nodes = membership[sts[-1]]["nodes"]

        root = sorted(self.mesh.get_node_set("RootNodes").nodes.keys())
        idx = self.mesh.node_id_to_index
        dpn = self.dpn
        rootd = sorted(idx[n] * dpn + d for n in root for d in range(dpn))

        profiles = {}
        solutions = {}
        for name, comp in (("flap", 1), ("edge", 0)):
            s = StaticLinearSolver(self.mesh, model)
            s.add_dirichlet_conditions([DirichletCondition(rootd, 0.0)])
            s.add_nodal_loads([NodalLoad([int(t) * dpn + comp],
                                         [load / len(tip_nodes)])
                               for t in tip_nodes])
            u = s.solve()
            u_mean = np.array([u[membership[st]["nodes"] * dpn + comp].mean()
                               for st in sts])
            profiles[name] = u_mean
            solutions[name] = u

        if frame == "tangent":
            return self._ei_tangent_frame(zs, [membership[st]["nodes"] for st in sts],
                                          solutions, load, stations)

        # Note: the global load directions bend the twisted blade with a
        # small flap/edge mixing inboard (twist ~15 deg at the root, ~0-2
        # deg outboard).  Outboard — where the reference EI values are
        # smallest and dominate the comparison — the +X/+Y loads are nearly
        # pure edgewise/flapwise, so the raw global profiles are used.
        ei = {}
        sel = select_stations(zs, stations)
        zs = zs[sel]
        for name in ("flap", "edge"):
            sp = UnivariateSpline(zs, profiles[name][sel], s=1e-2)
            kappa = np.abs(sp.derivative(2)(zs))
            M = load * (L - zs)
            ei[name] = np.where(kappa > 1e-12, M / np.maximum(kappa, 1e-12), np.nan)
        return zs, ei["flap"], ei["edge"]

    def _ei_tangent_frame(self, zs, station_nodes, solutions, load, stations=None):
        """Tangent-frame EI from the two tip-load displacement fields (see
        ``ei_from_static_response``, ``frame="tangent"``)."""
        dpn = self.dpn
        sel = select_stations(zs, stations)
        zs = np.asarray(zs)[sel]
        station_nodes = [station_nodes[i] for i in sel]
        n_st = len(zs)
        cen = np.array([self.coords[n].mean(axis=0) for n in station_nodes])
        half = 2.0 * float(np.median(np.diff(zs)))

        def local_fit(z, y, deriv):
            return _local_quadratic(zs, y, half, deriv)

        slope_xy = np.column_stack([local_fit(zs, cen[:, d], True) for d in (0, 1)])
        t = np.column_stack([slope_xy, np.ones(n_st)])
        t /= np.linalg.norm(t, axis=1, keepdims=True)
        x_hat = np.array([1.0, 0.0, 0.0])
        edge = x_hat[None, :] - (t @ x_hat)[:, None] * t
        edge /= np.linalg.norm(edge, axis=1, keepdims=True)
        flap = np.cross(t, edge)

        tip_c = cen[-1]
        ei = {}
        for name, comp, bend_axis in (("flap", 1, edge), ("edge", 0, flap)):
            u = solutions[name].reshape(-1, dpn)[:, :3]
            theta = np.zeros((n_st, 3))
            for i, nodes in enumerate(station_nodes):
                theta[i] = rigid_rotation_fit(self.coords[nodes], cen[i], u[nodes])
            kappa = np.column_stack([local_fit(zs, theta[:, k], True)
                                     for k in range(3)]) * t[:, 2:3]
            force = np.zeros(3)
            force[comp] = load
            moment = np.cross(tip_c[None, :] - cen, force[None, :])
            m_b = np.einsum("ij,ij->i", moment, bend_axis)
            k_b = np.abs(np.einsum("ij,ij->i", kappa, bend_axis))
            ei[name] = np.where(k_b > 1e-12, np.abs(m_b) / np.maximum(k_b, 1e-12), np.nan)
        return zs, ei["flap"], ei["edge"]

    def section_stiffness(self, st: int, membership=None) -> Dict[str, float]:
        """Return {EA, GJ, EI_flap, EI_edge} for one blade station.

        Method: the LOWER face ring (shared with the previous station) is
        clamped; a rigid section motion is prescribed on the UPPER face ring
        (relative extension / twist / curvature, rotations projected onto
        each node's shell plane); interior nodes are free.  The section
        resultants are the net reactions of the constrained slice.

        Section frame: x = tangent (span), y = chord (edge), z = flap
        (thickness).  Bending about y is flapwise (EI_flap), bending about z
        is edgewise (EI_edge).  Curvature = rotation / station length.
        """
        if membership is None:
            membership = self.station_membership()
        sl = self._station_slice(st, membership)
        lower = sl["lower"]
        upper = sl["upper"]
        if upper.size == 0 or lower.size == 0:
            raise ValueError(f"station {st}: empty face ring")
        interior = sl["interior"]
        ref = sl["centroid"]
        tangent, chord, flap = self._section_axes(membership[st]["nodes"])
        dz = max(sl["z1"] - sl["z0"], 1e-9)

        all_nodes = np.concatenate([lower, interior, upper])
        dofs = _dofs_of_nodes(all_nodes, self.dpn)
        K_loc = self._submatrix(dofs).tocsc()
        n_all = all_nodes.size

        eps = 1e-6      # prescribed axial strain
        kappa = 1e-6    # prescribed face relative rotation (rad)

        def run_state(upper_trans: np.ndarray,
                      upper_rots: Optional[np.ndarray]) -> tuple:
            pres_dofs = list(np.concatenate([lower * self.dpn + d for d in range(6)]))
            pres_vals = [0.0] * (lower.size * 6)
            for d in range(3):
                pres_dofs.extend((upper * self.dpn + d).tolist())
            pres_vals.extend(np.asarray(upper_trans, dtype=np.float64).ravel().tolist())
            if upper_rots is not None:
                for d in range(3, 6):
                    pres_dofs.extend((upper * self.dpn + d).tolist())
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
                u_loc[f_idx] = spsolve(
                    K_loc[f_idx][:, f_idx].tocsc(),
                    -(K_loc[f_idx][:, p_idx].tocsc() @ u_p))

            r = (K_loc @ u_loc).reshape(n_all, self.dpn)
            # net reaction of the constrained slice = the section force
            F = r[:, :3].sum(axis=0)
            M = np.cross(self.coords[all_nodes] - ref, r[:, :3]).sum(axis=0)
            return F, M

        def rigid_rot(axis: np.ndarray, pts: np.ndarray, r0: np.ndarray,
                      theta: float) -> np.ndarray:
            return theta * np.cross(axis, pts - r0)

        def pure_bend(axis: np.ndarray, pts: np.ndarray, r0: np.ndarray,
                      theta: float) -> np.ndarray:
            d_xy = (pts - r0).copy()
            d_xy[:, 2] = 0.0
            vals = np.zeros((pts.shape[0], 3))
            vals[:, 2] = theta * (axis[0] * d_xy[:, 1] - axis[1] * d_xy[:, 0])
            return vals

        upper_pts = self.coords[upper]

        # axial: relative extension eps over the station
        u_ax = np.column_stack([np.zeros((upper.size, 2)),
                                np.full(upper.size, eps * dz)])
        F, _ = run_state(u_ax, None)
        ea = abs(float(F[2])) / eps

        # torsion: relative twist kappa about the tangent
        u_t = rigid_rot(tangent, upper_pts, ref, kappa)
        rots_t = self._projected_rots(upper, tangent, kappa)
        _, M_t = run_state(u_t, rots_t)
        gj = abs(float(np.dot(M_t, tangent))) * dz / kappa

        # flap curvature about the chord axis
        u_b = pure_bend(chord, upper_pts, ref, kappa)
        rots_b = self._projected_rots(upper, chord, kappa)
        _, M_c = run_state(u_b, rots_b)
        ei_flap = abs(float(np.dot(M_c, chord))) * dz / kappa

        # edge curvature about the flap axis
        u_b2 = pure_bend(flap, upper_pts, ref, kappa)
        rots_b2 = self._projected_rots(upper, flap, kappa)
        _, M_f = run_state(u_b2, rots_b2)
        ei_edge = abs(float(np.dot(M_f, flap))) * dz / kappa

        return {"EA": ea, "GJ": gj, "EI_flap": ei_flap, "EI_edge": ei_edge}

    def section_mass(self, st: int, membership=None) -> Dict[str, float]:
        """Return {mass_dens, cg} for one blade station."""
        if membership is None:
            membership = self.station_membership()
        sl = self._station_slice(st, membership)
        m = 0.5 * self.node_mass[sl["lower"]].sum() + \
            0.5 * self.node_mass[sl["upper"]].sum() + \
            self.node_mass[sl["interior"]].sum()
        dz = max(sl["z1"] - sl["z0"], 1e-9)
        return {"mass_dens": float(m / dz), "cg": sl["centroid"]}

    def extract(self, n_slices: Optional[int] = None,
                stations: Optional[np.ndarray] = None) -> SectionProperties:
        """Extract sectional properties per blade station (via the detailed
        element-set membership — robust for twisted sections).

        Nothing here is fitted across stations (each value is a per-bay quantity), so the result
        does not depend on the station count; ``stations`` (reference span coordinates, see
        ``select_stations``) only restricts the evaluation to the nearest mesh stations so the
        sampled z set matches the other S-1 channels. ``None`` keeps every station.
        """
        membership = self.station_membership()
        zs, ea, gj, eif, eie, m, cg = [], [], [], [], [], [], []
        order = sorted(membership)
        keep = select_stations(np.array([membership[st]["z"] for st in order]), stations)
        for st in (order[i] for i in keep):
            try:
                sec = self.section_stiffness(st, membership)
                ms = self.section_mass(st, membership)
            except (ValueError, RuntimeError):
                continue
            zs.append(membership[st]["z"])
            ea.append(sec["EA"])
            gj.append(sec["GJ"])
            eif.append(sec["EI_flap"])
            eie.append(sec["EI_edge"])
            m.append(ms["mass_dens"])
            cg.append(ms["cg"])
        return SectionProperties(
            z=np.asarray(zs), ea=np.asarray(ea), gj=np.asarray(gj),
            ei_flap=np.asarray(eif), ei_edge=np.asarray(eie),
            mass_dens=np.asarray(m), cg=np.asarray(cg),
        )


def select_stations(zs, ref_z=None):
    """Indices (ascending, unique) of the stations of ``zs`` nearest to each ``ref_z``.

    The S-1 curvature, twist-rate and strain fits depend on the station count (smoothing bound,
    window of twice the median spacing). Fitting on the stations nearest to a fixed reference grid
    (the yaml ``reference_axis`` grid) makes them independent of the stations the mesher inserted.
    ``ref_z=None`` selects every station; when the stations already are the reference grid the
    selection is the identity.
    """
    zs = np.asarray(zs, dtype=float)
    if ref_z is None:
        return np.arange(zs.size)
    ref = np.asarray(ref_z, dtype=float)
    return np.unique(np.abs(zs[:, None] - ref[None, :]).argmin(axis=0))


def rigid_rotation_fit(xyz, centre, u):
    """Rotation vector of the least-squares rigid motion of a node set.

    Fits ``u_k = u0 + theta x d_k`` (``d_k = xyz_k - centre``) to the nodal translations ``u`` (n, 3) by
    ``lstsq`` on the 6 unknowns ``(u0, theta)`` and returns ``theta``.
    """
    d = np.asarray(xyz) - centre
    n = d.shape[0]
    # u_k = u0 + theta x d_k  ->  (I | -[d_k]x) (u0, theta) = u_k
    A = np.zeros((3 * n, 6))
    A[:, :3] = np.tile(np.eye(3), (n, 1))
    A[0::3, 4], A[0::3, 5] = d[:, 2], -d[:, 1]
    A[1::3, 3], A[1::3, 5] = -d[:, 2], d[:, 0]
    A[2::3, 3], A[2::3, 4] = d[:, 1], -d[:, 0]
    return np.linalg.lstsq(A, np.asarray(u).ravel(), rcond=None)[0][3:]


def _local_quadratic(zs, y, half, deriv):
    """Value or slope of a local quadratic fit about every station.

    The window is ``|z - z_i| <= half``, widened to the three nearest stations when it holds fewer:
    non-uniform stations (``BladeMesh.max_bay_length``) can leave a station alone in the window.
    """
    zs = np.asarray(zs, dtype=float)
    out = np.empty(zs.size)
    for i in range(zs.size):
        d = np.abs(zs - zs[i])
        m = d <= half + 1e-9
        if m.sum() < 3:
            m = d <= np.sort(d)[min(2, zs.size - 1)] + 1e-12
        c = np.polyfit(zs[m] - zs[i], y[m], min(2, int(m.sum()) - 1))
        out[i] = c[-2] if deriv else c[-1]
    return out


def _station_of(name):
    """Station index of a detailed NuMAD set ``{component:02d}_{station:02d}_{...}``, else None.

    ``:02d`` grows past two digits beyond station 99 (``BladeMesh.max_bay_length``), so the field is
    read up to the next underscore, not as ``name[3:5]``.
    """
    parts = name.split("_", 2)
    if len(parts) < 3 or len(parts[0]) != 2 or not parts[0].isdigit():
        return None
    if len(parts[1]) < 2 or not parts[1].isdigit():
        return None
    return int(parts[1])
