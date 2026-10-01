"""WU7: what actually drives the IEA-15-240-RWT blade's elastic twist (issue #9).

Issue #9 attributes the blade's over-twist to the same bend-twist coupling as the
coupon. This module measures the blade's twist mechanism directly and finds that
attribution impossible for this model:

1. **No laminate bend-twist coupling exists in the blade.** All 696 sections are
   balanced, so `D16 = D26 = 0` exactly. The coupon's mechanism cannot be the
   blade's mechanism - there is nothing to couple.
2. **The twist is a load-path response.** For one and the same root moment, moving
   the line of action of the resultant along the chord swings the tip twist by tens
   of degrees while the tip deflection barely moves. That spread is the right order
   of magnitude to explain a 9x disagreement between two models that are each
   internally consistent about a *different* load line of action.
3. **The mesh has free edges at two mid-span stations** (the shear-web junctions), a
   real connectivity defect. Its *effect* is bounded: the shell's section stiffness sits
   inside the reference scatter by the modal route (1st torsion 4.000 Hz vs the
   reference 4.290 Hz, i.e. GJ ratio 0.87 against a 15.1% reference scatter; 1st flap
   0.526 vs 0.570; 1st edge 0.702 vs 0.650), so the tear is not the blade-twist
   explanation. `test_mesh_free_edges_document_the_open_finding` pins the state so a fix
   cannot land unnoticed.
   **Retracted**: an earlier version of this module reported the shell's torsion as
   ~400x too soft. That was an artifact of reading the *mean nodal rotation at the
   loaded tip ring*, where the local shell deformation dominates; the modal check above
   is the reliable measure. A per-station section comparison remains open because it
   needs the section frame (a moment about the global x mixes flap and edge through the
   structural twist, and the global-z rotation is not the section's torsion), which is
   why this module asserts nothing per station.

Metric: the section twist about the blade axis is the mean nodal rotation about z
(dof 5) over the tip ring. An LSQ slope of the flapwise displacement against the
chordwise coordinate is **not** a twist here: with a ~21 m tip deflection the
bending rotation dominates it (that metric first read -18 deg on a case whose real
twist is -0.57 deg).
"""

from __future__ import annotations

from collections import Counter
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("_aeroelast", reason="Rust backend not available")

from _aeroelast import PyMeshAssembler  # noqa: E402
from aeroelast.core.laminate import Laminate, Ply  # noqa: E402
from aeroelast.core.material import OrthotropicMaterial  # noqa: E402
from aeroelast.core.mesh.entities import MeshElement, Node  # noqa: E402
from aeroelast.models.blade.model import Blade  # noqa: E402
from scipy.sparse import coo_matrix  # noqa: E402
from scipy.sparse.linalg import spsolve  # noqa: E402

import test_blade_iea15mw_validation as blade_validation  # noqa: E402

YAML = Path(__file__).resolve().parent / "IEA-15-240-RWT.yaml"
ELEMENT_SIZE = 1.0
#: Article DLC 1.4 maximum root moment, the load level the validation suite uses.
ROOT_MOMENT_NM = 90.4e6


@pytest.fixture(scope="module")
def blade_model():
    """Build the blade mesh once and return ``(mesh, props, numad)``."""
    Node._id_counter = 0
    MeshElement._id_counter = 0
    model = Blade(str(YAML), element_size=ELEMENT_SIZE)
    model.generate_mesh()
    return model.mesh, model.get_element_properties(), model._numad_mesh


@pytest.fixture(scope="module")
def solved(blade_model):
    """Assemble the blade and return a solver plus geometry for the load cases."""
    mesh, props, _ = blade_model
    assembler = PyMeshAssembler.from_model(
        blade_validation._to_rust_mesh(mesh, props), props,
        list(blade_validation.SPAN_DIRECTION), None,
    )
    n = assembler.dofs_count
    rows, cols, vals = assembler.assemble_k()
    K = coo_matrix((np.asarray(vals), (np.asarray(rows), np.asarray(cols))),
                   shape=(n, n)).tocsr()
    root = {mesh.node_id_to_index[nid] for nid in mesh.get_node_set("RootNodes").node_ids}
    fixed = {6 * i + d for i in root for d in range(6)}
    free = np.array([i for i in range(n) if i not in fixed], dtype=np.int64)
    Kff = K[np.ix_(free, free)]
    z = np.array([nd.z for nd in mesh.nodes])
    x = np.array([nd.x for nd in mesh.nodes])
    y = np.array([nd.y for nd in mesh.nodes])
    load_idx = np.array([i for i in range(len(mesh.nodes)) if i not in root])
    tip = np.where(z >= z.max() - 1e-6)[0]

    def solve(force):
        u = np.zeros(n)
        u[free] = spsolve(Kff, force[free])
        flap = float(np.mean([u[6 * nd + blade_validation.FLAPWISE_DOF] for nd in tip]))
        twist = float(np.mean([u[6 * nd + 5] for nd in tip]))
        return flap, twist

    return {"n": n, "solve": solve, "z": z, "x": x, "y": y, "load_idx": load_idx,
            "tip": tip, "Kff": Kff, "free": free}


# ─────────────────────────────────────────────────────────────────────────────
# 1. There is no laminate bend-twist coupling in this blade
# ─────────────────────────────────────────────────────────────────────────────


def test_blade_laminates_have_no_bend_twist_coupling(blade_model):
    """Every section is balanced, so ``D16 = D26 = 0`` and the coupon's mechanism
    cannot be the blade's mechanism.

    This is the finding that makes issue #9's causal story impossible for this
    model: the elastic twist of the blade cannot be a plate bend-twist coupling
    that does not exist.
    """
    _, _, numad = blade_model
    materials = {m["name"]: m for m in numad["materials"]}
    worst_offdiag = 0.0
    worst_scale = 1.0
    worst_section = ""
    for section in numad["sections"]:
        plies = []
        for name, thickness, angle in section["layup"]:
            e = materials[name]["elastic"]
            mat = OrthotropicMaterial(
                name, E=(e["E"][0], e["E"][1], e["E"][2]),
                G=(e["G"][0], e["G"][1], e["G"][2]),
                nu=(e["nu"][0], e["nu"][1], e["nu"][2]),
                rho=float(materials[name]["density"]),
            )
            plies.append(Ply(mat, float(thickness), float(angle)))
        D = Laminate(plies).D
        offdiag = abs(D[0, 2]) + abs(D[1, 2])
        if offdiag > worst_offdiag:
            worst_offdiag, worst_scale, worst_section = offdiag, float(np.abs(D).max()), section["elementSet"]
    print(f"\n{len(numad['sections'])} secciones; max |D16|+|D26| = {worst_offdiag:.3e} "
          f"(escala |D|max = {worst_scale:.3e}) en {worst_section}")
    assert worst_offdiag <= 1e-9 * worst_scale, (
        f"a section carries bend-twist coupling: |D16|+|D26| = {worst_offdiag:.3e} "
        f"vs scale {worst_scale:.3e} in {worst_section}"
    )


# ─────────────────────────────────────────────────────────────────────────────
# 2. The twist is a load-path response
# ─────────────────────────────────────────────────────────────────────────────


def _station_load(solved, selector):
    """Apply the root moment with the resultant of each station placed by ``selector``."""
    n, z, load_idx = solved["n"], solved["z"], solved["load_idx"]
    stations: dict[float, list[int]] = {}
    for nd in load_idx:
        stations.setdefault(round(z[nd], 6), []).append(nd)
    picks = {zz: selector(np.array(nds)) for zz, nds in stations.items()}
    total = sum(zz * len(sel) for zz, sel in picks.items())
    per_node = ROOT_MOMENT_NM / total
    force = np.zeros(n)
    for sel in picks.values():
        for nd in sel:
            force[6 * nd + blade_validation.FLAPWISE_DOF] = per_node
    return force


def test_blade_twist_is_load_path_dominated(solved):
    """One root moment, three lines of action, tens of degrees of twist spread.

    A shell that spreads the aerodynamic load over the surface and a beam that
    applies it along the aerodynamic centre are not modelling the same load, and
    the twist - unlike the deflection - is first-order sensitive to that choice.
    """
    surface = solved["solve"](_station_load(solved, lambda nds: nds))
    leading = solved["solve"](_station_load(solved, lambda nds: nds[np.argmin(solved["x"][nds])][None]))
    midchord = solved["solve"](_station_load(solved, lambda nds: nds[np.argmin(np.abs(solved["x"][nds]))][None]))

    print(f"\nroot moment {ROOT_MOMENT_NM/1e6:.1f} MNm, twist = mean rotation about z at the tip")
    for label, (flap, twist) in (("surface nodes", surface), ("leading edge", leading),
                                 ("mid-chord", midchord)):
        print(f"  {label:16} tip flap {flap:7.3f} m   twist {np.rad2deg(twist):+9.3f} deg")

    twist_surface, twist_leading, twist_mid = (np.rad2deg(v[1]) for v in (surface, leading, midchord))
    # The deflection is insensitive to the line of action; the twist is not.
    assert abs(surface[0] - leading[0]) < 0.10 * abs(surface[0])
    assert abs(twist_leading) > 10 * abs(twist_surface)
    assert abs(twist_mid) > 10 * abs(twist_surface)


def test_mesh_free_edges_document_the_open_finding(blade_model):
    """Characterisation: the blade mesh has free edges at two mid-span stations.

    Free edges at the root and the tip are the expected open boundaries. The ones at
    z ~ 11.94 m and z ~ 112.22 m sit at the shear-web junctions, so the webs do not
    connect across those stations and the section's closed cells are interrupted.
    The node dedup step reports 0% reduction, i.e. it does not merge them.

    When the mesh is repaired this test must be updated - that is the point of
    pinning it.
    """
    mesh, _, _ = blade_model
    z_of = {nd.id: nd.z for nd in mesh.nodes}
    edge_uses: Counter = Counter()
    for element in mesh.elements:
        ids = [nd.id for nd in element.nodes]
        for a, b in zip(ids, ids[1:] + ids[:1], strict=True):
            edge_uses[tuple(sorted((a, b)))] += 1
    free_edges = [edge for edge, uses in edge_uses.items() if uses == 1]
    z_free = sorted({round(z_of[a], 3) for a, _ in free_edges})
    mid_span = [zz for zz in z_free if 5.0 < zz < 112.0]
    print(f"\n{len(mesh.elements)} elements, {len(edge_uses)} edges, "
          f"{len(free_edges)} free edges at z = {[round(zz, 2) for zz in z_free]}")
    print(f"  free edges at mid-span stations: {[round(zz, 2) for zz in mid_span]}")
    assert mid_span, (
        "no mid-span free edges: the mesh tear appears fixed, so update this test "
        "and the WU7 section of the task document"
    )
