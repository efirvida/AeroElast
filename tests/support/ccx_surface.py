"""CalculiX S8R decks for a structured shell patch on a parametric surface.

The patch is ``nu x nv`` eight-node shell elements over ``point(u, v)`` with ``u, v`` in [0, 1].
It exists to give benchmark tests an arbiter that is not our element: CalculiX S8R is a different
code and a different (quadratic) formulation, so agreement is evidence and not a restatement.
"""

from __future__ import annotations

from collections.abc import Callable, Iterable
from pathlib import Path

import numpy as np

Point = Callable[[float, float], tuple[float, float, float]]
Boundary = Callable[[float, float, float], Iterable[int]]  # xyz -> constrained dofs, 1-based


def node_id(i: int, j: int, nu: int, *, periodic: bool = False) -> int:
    """Id of the grid point ``(i, j)``, ``i`` in ``0..2nu``, ``j`` in ``0..2nv``.

    With ``periodic`` the column ``2nu`` is the column ``0`` (a closed loop in ``u``).
    """
    if periodic:
        return j * (2 * nu) + (i % (2 * nu)) + 1
    return j * (2 * nu + 1) + i + 1


def write_s8r_surface_deck(
    inp_path: Path,
    point: Point,
    nu: int,
    nv: int,
    *,
    young: float,
    poisson: float,
    thickness: float,
    boundary: Boundary,
    density: float = 0.0,
    gravity: tuple[float, float, float] | None = None,
    nodal_loads: dict[int, tuple[int, float]] | None = None,
    print_nodes: Iterable[int] = (),
    periodic_u: bool = False,
    pressure: tuple[Callable[[int, int], bool], float] | None = None,
    nlgeom: int | None = None,
) -> dict[int, np.ndarray]:
    """Write the deck and return the node coordinates by id.

    ``gravity`` is the acceleration vector applied with ``*DLOAD ... GRAV`` (needs ``density``).
    ``nodal_loads`` maps a node id to ``(dof, value)``. ``print_nodes`` are written to the ``.dat``
    with ``*NODE PRINT`` (read them with ``read_dat_displacements``): the ``.frd`` of a shell model
    holds the expanded 3D nodes, not the ids of this deck. ``periodic_u`` closes the patch into a
    loop in ``u``. ``pressure`` is ``(select, value)``: ``select(ei, ej)`` picks the elements that
    get ``*DLOAD ... P`` (positive pressure acts against the element normal ``dr/du x dr/dv``).
    ``nlgeom`` is the number of load increments of a geometrically nonlinear step
    (``*STEP, NLGEOM``); ``CLOAD`` keeps its global direction, i.e. a dead load.
    """
    coords: dict[int, np.ndarray] = {}
    for j in range(2 * nv + 1):
        for i in range(2 * nu if periodic_u else 2 * nu + 1):
            coords[node_id(i, j, nu, periodic=periodic_u)] = np.asarray(
                point(i / (2 * nu), j / (2 * nv)), dtype=float
            )

    lines = ["*NODE"]
    for nid, xyz in coords.items():
        lines.append(f"{nid}, {xyz[0]:.12e}, {xyz[1]:.12e}, {xyz[2]:.12e}")

    lines.append("*ELEMENT, TYPE=S8R, ELSET=SHELL")
    eid = 0
    pressured: list[int] = []
    for ej in range(nv):
        for ei in range(nu):
            i, j = 2 * ei, 2 * ej
            ring = [
                (i, j), (i + 2, j), (i + 2, j + 2), (i, j + 2),
                (i + 1, j), (i + 2, j + 1), (i + 1, j + 2), (i, j + 1),
            ]  # fmt: skip
            eid += 1
            if pressure is not None and pressure[0](ei, ej):
                pressured.append(eid)
            lines.append(
                f"{eid}, " + ", ".join(str(node_id(a, b, nu, periodic=periodic_u)) for a, b in ring)
            )
    if pressured:
        lines.append("*ELSET, ELSET=PRESSURED")
        lines += [f"{e}," for e in pressured]

    lines += [
        "*MATERIAL, NAME=MAT",
        "*ELASTIC",
        f"{young:.12e}, {poisson:.12e}",
    ]
    if density:
        lines += ["*DENSITY", f"{density:.12e}"]
    lines.append(f"*SHELL SECTION, ELSET=SHELL, MATERIAL=MAT\n{thickness:.12e}")

    lines.append("*BOUNDARY")
    for nid, xyz in coords.items():
        for dof in sorted(set(boundary(*map(float, xyz)))):
            lines.append(f"{nid}, {dof}, {dof}")

    if nlgeom is None:
        lines.append("*STEP\n*STATIC")
    else:
        lines.append(f"*STEP, NLGEOM, INC={20 * nlgeom}\n*STATIC\n{1.0 / nlgeom:.6e}, 1.0, 1e-8, {1.0 / nlgeom:.6e}")
    if gravity is not None:
        g = float(np.linalg.norm(gravity))
        d = np.asarray(gravity, dtype=float) / g
        lines.append(f"*DLOAD\nSHELL, GRAV, {g:.12e}, {d[0]:.12e}, {d[1]:.12e}, {d[2]:.12e}")
    if pressure is not None and pressured:
        lines.append(f"*DLOAD\nPRESSURED, P, {pressure[1]:.12e}")
    if nodal_loads:
        lines.append("*CLOAD")
        for nid, (dof, value) in nodal_loads.items():
            lines.append(f"{nid}, {dof}, {value:.12e}")
    printed = sorted(set(print_nodes))
    if printed:
        lines.append("*NSET, NSET=PRINTED")
        lines += [f"{nid}," for nid in printed]
        lines += ["*NODE PRINT, NSET=PRINTED", "U"]
    lines.append("*END STEP")
    inp_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return coords


def read_dat_displacements(dat_path: Path) -> dict[int, np.ndarray]:
    """Displacements ``{node: (ux, uy, uz)}`` from the ``*NODE PRINT`` block of a ``.dat``."""
    out: dict[int, np.ndarray] = {}
    in_block = False
    for line in dat_path.read_text(encoding="utf-8", errors="replace").splitlines():
        if "displacements" in line:
            in_block = True
            continue
        if not in_block:
            continue
        parts = line.split()
        if len(parts) == 4:
            try:
                out[int(parts[0])] = np.asarray([float(x) for x in parts[1:]])
            except ValueError:
                in_block = False
        elif parts:
            in_block = False
    return out
