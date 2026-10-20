"""Optional reference-surface transforms for blade shell meshes."""

from collections import defaultdict

import numpy as np


def apply_midsurface_inset(mesh, props, *, verbose: bool = True) -> float:
    """Inset an OML blade mesh in place; return mean moved-node distance in metres.

    Use solver properties from ``build_rust_properties``. Each outer-skin node
    moves by the mean half-thickness of its adjacent skin laminates along the
    averaged skin normal projected into its section's XY plane. Nodes belonging
    to ``*_TE_FLAT`` elements stay on the OML to avoid trailing-edge inversions;
    interior/web-only nodes are untouched (shared skin/web boundary nodes move).

    This is a single OML-to-mid transform, not an incremental mesh operation:
    call it once on a fresh OML mesh, before assembly or other mesh consumers.
    """
    x0 = mesh.coords_array.copy()
    idx = mesh.node_id_to_index
    outer = {idx[n] for n in mesh.get_node_set("allOuterShellNods").nodes}
    thickness = {}
    for name in mesh.element_sets_names:
        laminate = props.get(name)
        if laminate is None or not hasattr(laminate, "total_thickness"):
            continue
        for element in mesh.get_element_set(name).elements:
            thickness[getattr(element, "id", element)] = laminate.total_thickness

    normal_sum = defaultdict(lambda: np.zeros(3))
    thickness_sum = defaultdict(float)
    count = defaultdict(int)
    for element in mesh.elements:
        ids = [idx[n] for n in element.node_ids]
        if not all(i in outer for i in ids) or element.id not in thickness:
            continue
        p = x0[ids]
        normal = (
            np.cross(p[2] - p[0], p[3] - p[1])
            if len(ids) == 4 else np.cross(p[1] - p[0], p[2] - p[0])
        )
        normal /= np.linalg.norm(normal)
        for i in ids:
            normal_sum[i] += normal
            thickness_sum[i] += thickness[element.id]
            count[i] += 1

    z = np.round(x0[:, 2], 6)
    sections = defaultdict(list)
    for i in outer:
        sections[z[i]].append(i)
    centers = {height: x0[nodes, :2].mean(0) for height, nodes in sections.items()}
    votes = sum(
        np.sign(np.dot(normal_sum[i][:2], centers[z[i]] - x0[i, :2]))
        for i in normal_sum
    )
    sign = 1.0 if votes > 0 else -1.0
    te_elements = set()
    for name in mesh.element_sets_names:
        if name.endswith("_TE_FLAT"):
            te_elements.update(
                getattr(e, "id", e) for e in mesh.get_element_set(name).elements
            )
    te_nodes = {
        idx[n] for element in mesh.elements if element.id in te_elements
        for n in element.node_ids
    }
    x1 = x0.copy()
    shifts = []
    for i in normal_sum:
        if i in te_nodes:
            continue
        normal = sign * normal_sum[i]
        normal[2] = 0.0
        normal /= np.linalg.norm(normal)
        half_thickness = 0.5 * thickness_sum[i] / count[i]
        x1[i] += half_thickness * normal
        shifts.append(half_thickness)
    mesh.coords_array = x1
    if verbose:
        print(
            f"inset: {len(shifts)} nodes, sign {sign:+.0f} (votes {votes:+.0f}), "
            f"mean half-thickness {np.mean(shifts) * 1e3:.1f} mm, "
            f"max {np.max(shifts) * 1e3:.1f} mm", flush=True,
        )
    return float(np.mean(shifts))
