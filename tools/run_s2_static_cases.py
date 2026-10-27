"""S-2 — static prescribed-load validation: AeroElast shell vs 1D beam reference.

Runs four static load cases on the IEA 15 MW blade (clamped root):
  LC1  gravity                 (0, -9.81, 0) m/s²
  LC2  tip load, flap-ish      (0, +1, 0) × 1 MN at tip
  LC3  tip load, edge-ish      (+1, 0, 0) × 1 MN at tip
  LC4  distributed load        (0, +1, 0) × 5 kN/m uniform over the span

Both models share the SAME geometry source (the AeroElast BladeMesh): the
beam reference derives its node positions and section orientation from the
shell mesh, while its stiffness/mass come from the official ElastoDyn +
BeamDyn distributed properties.  This isolates shell-vs-beam discretization.

Output:
  docs/validation_data/generated/s2_static_prescribed.csv  (case summary)
  prints a comparison table with relative errors.

Usage (SDumont):
  module load glu gcc/14.2.0_sequana
  python tools/run_s2_static_cases.py [--element-size 0.5] [--n-slices 60]
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(REPO_ROOT / "tools"))

from tools.beam_reference import (  # noqa: E402
    BeamReference,
    load_beamdyn_blade,
    load_elastodyn_blade,
)

REF_DIR = REPO_ROOT / "tests" / "IEA15MW" / "reference"
ELASTODYN_FILE = REF_DIR / "IEA-15-240-RWT_ElastoDyn_blade.dat"
BEAMDYN_FILE = REF_DIR / "IEA-15-240-RWT_BeamDyn_blade.dat"
NEG_Z_THRESH = 0.0  # blade spans z in [0, 117]

GRAVITY = (0.0, -9.81, 0.0)
TIP_FORCE_N = 1.0e6
DISTRIB_Q_NM = 5.0e3


def run_shell_case(mesh, properties, load) -> np.ndarray:
    """Run one StaticLinearSolver case. load: dict with keys type/params."""
    from aeroelast.core.bc import BodyForce, DirichletCondition, NodalLoad
    from aeroelast.elements import ElementFamily
    from aeroelast.models.blade.model import build_rust_properties
    from aeroelast.solvers.elasticity.static_linear import StaticLinearSolver

    cfg = {
        "solver": {},
        "elements": {
            "element_family": ElementFamily.SHELL,
            "span_direction": (0.0, 0.0, 1.0),
            "properties": properties,
        },
    }
    solver = StaticLinearSolver(mesh, cfg)

    dpn = solver.domain.dofs_per_node
    root_node_ids = sorted(mesh.get_node_set("RootNodes").nodes.keys())
    node_id_to_idx = mesh.node_id_to_index
    root_dofs = sorted(
        node_id_to_idx[nid] * dpn + d for nid in root_node_ids for d in range(dpn)
    )
    solver.add_dirichlet_conditions([DirichletCondition(root_dofs, 0.0)])

    if load["type"] == "gravity":
        solver.add_body_forces([BodyForce(list(load["value"]))])
    elif load["type"] == "tip":
        coords = mesh.coords_array
        z = coords[:, 2]
        tip_mask = z >= z.max() - 1e-6
        tip_idxs = np.nonzero(tip_mask)[0]
        force = np.asarray(load["value"], dtype=np.float64)
        if load.get("distribution", "equal") == "arc":
            # work-equivalent uniform line traction on the free tip edge: a line load stays finite on a
            # Reissner-Mindlin shell, equal nodal forces are point loads (l3-closure T10c)
            tip_idxs, w = ring_arc_weights(mesh, tip_idxs)
            f_nodes = np.outer(w / w.sum(), force)
        else:
            f_nodes = np.tile(force / len(tip_idxs), (len(tip_idxs), 1))
        dofs, vals = [], []
        for i, f_node in zip(tip_idxs, f_nodes, strict=True):
            for d in range(3):
                dofs.append(int(i) * dpn + d)
                vals.append(float(f_node[d]))
        solver.add_nodal_loads([NodalLoad(dofs, vals)])
    elif load["type"] == "distributed":
        coords = mesh.coords_array
        z = coords[:, 2]
        q = np.asarray(load["value"], dtype=np.float64)
        n_slices = load["n_slices"]
        edges = np.linspace(z.min(), z.max(), n_slices + 1)
        area = nodal_areas(mesh) if load.get("distribution", "equal") == "area" else None
        for a, bnd in zip(edges[:-1], edges[1:]):
            mask = (z >= a) & (z < bnd) if bnd < z.max() else (z >= a) & (z <= bnd)
            idxs = np.nonzero(mask)[0]
            if len(idxs) == 0:
                continue
            # "area": a uniform surface traction per slice (tributary-area nodal forces), not point loads
            share = area[idxs] / area[idxs].sum() if area is not None else np.full(len(idxs), 1.0 / len(idxs))
            dofs, vals = [], []
            for i, c in zip(idxs, share, strict=True):
                for d in range(3):
                    dofs.append(int(i) * dpn + d)
                    vals.append(float(q[d] * (bnd - a) * c))
            solver.add_nodal_loads([NodalLoad(dofs, vals)])
    else:
        raise ValueError(f"unknown load type {load['type']}")

    return solver.solve()


def shell_tip_displacement(mesh, u: np.ndarray) -> np.ndarray:
    dpn = mesh.dofs_per_node if hasattr(mesh, "dofs_per_node") else 6
    coords = mesh.coords_array
    tip_idx = int(np.argmax(coords[:, 2]))
    return np.asarray(u[tip_idx * dpn : tip_idx * dpn + 3], dtype=np.float64)


def _element_node_indices(mesh) -> list[list[int]]:
    idx = mesh.node_id_to_index
    return [[idx[n.id] for n in e.nodes] for e in mesh.elements]


def ring_arc_weights(mesh, ring: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Arc-length weight of each ring node: half the length of the element edges joining it to ring neighbours.

    The edges are the element sides with both corner nodes on the ring, so the weights sum to the ring perimeter
    whatever the node spacing; nodes with no ring edge are dropped (they carry no arc length).
    """
    coords = mesh.coords_array
    on_ring = {int(i) for i in ring}
    w = {}
    seen = set()
    for nodes in _element_node_indices(mesh):
        for a, b in zip(nodes, nodes[1:] + nodes[:1], strict=True):
            if a in on_ring and b in on_ring and (key := (min(a, b), max(a, b))) not in seen:
                seen.add(key)
                half = 0.5 * float(np.linalg.norm(coords[a] - coords[b]))
                w[a] = w.get(a, 0.0) + half
                w[b] = w.get(b, 0.0) + half
    idxs = np.array(sorted(w), dtype=int)
    return idxs, np.array([w[i] for i in idxs])


def nodal_areas(mesh) -> np.ndarray:
    """Tributary shell area per node: each element's area split equally among its corner nodes."""
    coords = mesh.coords_array
    area = np.zeros(len(coords))
    for nodes in _element_node_indices(mesh):
        x = coords[nodes]
        a = 0.5 * np.linalg.norm(np.cross(x[1] - x[0], x[2] - x[0]))
        if len(nodes) == 4:
            a += 0.5 * np.linalg.norm(np.cross(x[2] - x[0], x[3] - x[0]))
        area[nodes] += a / len(nodes)
    return area


def ring_mean_at(mesh, u: np.ndarray, z_target: float, weighting="node") -> tuple[float, np.ndarray, int]:
    """Mean translation of the node ring on the mesh z-plane nearest ``z_target``: (z, mean u, node count).

    ``weighting="node"`` is the arithmetic node mean (biased toward dense parts of the ring, and stepping where the
    ring node count changes: T10b); ``"arc"`` weights each node by its arc length (``ring_arc_weights``).
    """
    z = mesh.coords_array[:, 2]
    z_s = float(z[np.argmin(np.abs(z - z_target))])
    ring = np.nonzero(np.abs(z - z_s) < 1e-6)[0]
    uu = np.asarray(u).reshape(-1, 6)
    if weighting == "arc":
        ring, w = ring_arc_weights(mesh, ring)
        return z_s, (uu[ring, :3] * w[:, None]).sum(axis=0) / w.sum(), len(ring)
    return z_s, uu[ring, :3].mean(axis=0), len(ring)


def saint_venant_tip(mesh, u: np.ndarray, scale: float, depths=(1.0, 1.5, 2.0, 2.5, 3.0), degree=3,
                     weighting="node") -> np.ndarray:
    """Tip translation extrapolated from ring means ``depths * scale`` inboard of the tip (l3-closure T10a).

    A tip load split equally over the tip ring bends flat walls out of their plane; that local term decays within
    about one section depth (Saint-Venant), so the rings inboard carry the beam translation. A polynomial in z
    through them, evaluated at the tip, estimates the beam-equivalent tip translation. ``scale`` is the section
    size that sets the decay length (the largest tip-ring extent).
    """
    z_tip = float(mesh.coords_array[:, 2].max())
    pts = [ring_mean_at(mesh, u, z_tip - d * scale, weighting) for d in depths]
    zs = np.array([p[0] for p in pts])
    us = np.array([p[1] for p in pts])
    return np.array([np.polyval(np.polyfit(zs - z_tip, us[:, k], degree), 0.0) for k in range(3)])


def tip_ring_scale(mesh) -> float:
    """Largest in-plane extent of the tip ring [m] (the chord at a blade tip)."""
    c = mesh.coords_array
    ring = c[c[:, 2] >= c[:, 2].max() - 1e-6, :2]
    return float(np.max(np.linalg.norm(ring[:, None, :] - ring[None, :, :], axis=-1)))


def beam_tip_displacement(beam: BeamReference, load) -> np.ndarray:
    if load["type"] == "gravity":
        u = beam.solve_gravity(load["value"])
    elif load["type"] == "tip":
        u = beam.solve_tip_load(load["value"])
    elif load["type"] == "distributed":
        u = beam.solve_distributed_load(load["value"])
    else:
        raise ValueError(f"unknown load type {load['type']}")
    return np.asarray(u[-6:-3], dtype=np.float64)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--element-size", type=float, default=0.5)
    ap.add_argument("--out-csv", default=None, help="default: docs/validation_data/generated/s2_static_prescribed.csv")
    ap.add_argument("--ei-source", choices=("elastodyn", "beamdyn"), default="elastodyn")
    ap.add_argument("--n-slices", type=int, default=60)
    ap.add_argument("--blade-xlsx", default=str(REPO_ROOT / "tests" / "NuMAD_utd_iea15mw.xlsx"))
    ap.add_argument("--blade-yaml", default=None,
                    help="WindIO blade yaml (overrides --blade-xlsx)")
    ap.add_argument("--airfoil-dir", default=str(REPO_ROOT / "tests" / "airfoils"))
    args = ap.parse_args()

    from aeroelast.core.mesh.generators import BladeMesh
    from aeroelast.models.blade.model import build_rust_properties

    print("[S2] generating blade mesh (element_size=%.2f)..." % args.element_size)
    generator = BladeMesh(
        excel_file=args.blade_xlsx if not args.blade_yaml else None,
        yaml_file=args.blade_yaml,
        airfoil_dir=args.airfoil_dir,
        element_size=args.element_size,
    )
    mesh = generator.generate(renumber="rcm")
    properties = build_rust_properties(generator.numad_mesh_data)
    print(f"[S2] mesh: {mesh.coords_array.shape[0]} nodes")

    ed = load_elastodyn_blade(ELASTODYN_FILE)
    bd = load_beamdyn_blade(BEAMDYN_FILE)
    beam = BeamReference.from_mesh(mesh, ed, bd, n_slices=args.n_slices, ei_source=args.ei_source)
    print(f"[S2] beam reference: {beam.n_nodes} slices, "
          f"integrated mass={beam.total_mass_kg:.0f} kg")

    cases = [
        ("LC1_gravity", {"type": "gravity", "value": GRAVITY}),
        ("LC2_tip_flap", {"type": "tip", "value": (0.0, TIP_FORCE_N, 0.0)}),
        ("LC3_tip_edge", {"type": "tip", "value": (TIP_FORCE_N, 0.0, 0.0)}),
        ("LC4_distributed", {"type": "distributed", "value": (0.0, DISTRIB_Q_NM, 0.0),
                             "n_slices": args.n_slices}),
    ]

    rows = []
    print(f"\n{'case':18s} {'component':10s} {'shell [m]':>12s} {'beam [m]':>12s} {'rel err':>10s}")
    print("-" * 66)
    for name, load in cases:
        u_shell = run_shell_case(mesh, properties, load)
        d_shell = shell_tip_displacement(mesh, u_shell)
        d_beam = beam_tip_displacement(beam, load)
        for comp, i in (("ux", 0), ("uy", 1), ("uz", 2)):
            if abs(d_beam[i]) < 1e-9 and abs(d_shell[i]) < 1e-9:
                rel = np.nan
            else:
                denom = abs(d_beam[i]) if abs(d_beam[i]) > 1e-9 else abs(d_shell[i])
                rel = (d_shell[i] - d_beam[i]) / denom
            rel_s = "n/a" if (rel is None or not np.isfinite(rel)) else f"{rel:+.2%}"
            print(f"{name:18s} {comp:10s} {d_shell[i]:12.4e} {d_beam[i]:12.4e} {rel_s:>10s}")
            rows.append({"case": name, "component": comp,
                         "shell_m": d_shell[i], "beam_m": d_beam[i],
                         "rel_err": rel})
        print()

    out_csv = (Path(args.out_csv) if args.out_csv else
               REPO_ROOT / "docs" / "validation_data" / "generated" / "s2_static_prescribed.csv")
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with open(out_csv, "w") as f:
        f.write("case,component,shell_m,beam_m,rel_err\n")
        for r in rows:
            f.write(f"{r['case']},{r['component']},{r['shell_m']:.8e},"
                    f"{r['beam_m']:.8e},{r['rel_err']}\n")
    print(f"[S2] wrote {out_csv}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
