"""The CCX writer's node/element labels must be correct for any entity-id scheme.

The entity id counters (``Node._id_counter``, ``MeshElement._id_counter``) are
process-global and never reset by ``MeshModel`` itself, so a mesh built after
another one has node and element ids that are neither 0-based nor contiguous.
The writer's ``.msh`` labels are 1-based *indices* into ``mesh.nodes`` /
``mesh.elements``; before this test the ``.nam`` element/node sets and the
element connectivity were labelled from the raw ids, so with drifted ids the
deck silently referenced nodes and elements that did not exist:

    *NODE, NSET=Nall            node labels 1..4
    *ELEMENT, TYPE=S4           connectivity 1001, 1002, 1004, 1003   <- absent
    *ELSET, ELSET=EPLATE        501                                   <- absent

The two tests here pin that writer behaviour: the exported deck must be
byte-identical whether the ids start at 0 or at an arbitrary offset, and every
label the deck writes must point at an entity that exists. They validate the
writer itself, not a physical quantity against an independent reference.
"""

from __future__ import annotations

from pathlib import Path

import pytest

pytest.importorskip("_aeroelast", reason="Rust backend not available")

import numpy as np  # noqa: E402
from tests.conftest import ccx_bin_or_skip  # noqa: E402
from tests.support.ccx_io import parse_frd_disp, run_ccx  # noqa: E402
from tests.support.ccx_plate import OFFSET, build_mesh, write_deck  # noqa: E402


def _parse_msh(msh_path: Path) -> tuple[set[int], set[int], list[list[int]]]:
    """Return (node labels, element labels, element connectivity labels)."""
    node_labels: set[int] = set()
    element_labels: set[int] = set()
    connectivity: list[list[int]] = []
    mode = None
    for line in msh_path.read_text().splitlines():
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.startswith("*"):
            upper = stripped.upper()
            if upper.startswith("*NODE"):
                mode = "node"
            elif upper.startswith("*ELEMENT"):
                mode = "element"
            else:
                mode = None
            continue
        parts = [p.strip() for p in stripped.split(",")]
        if mode == "node":
            node_labels.add(int(parts[0]))
        elif mode == "element":
            element_labels.add(int(parts[0]))
            connectivity.append([int(p) for p in parts[1:] if p])
    return node_labels, element_labels, connectivity


def _parse_nam(nam_path: Path) -> dict[str, set[int]]:
    """Return set name -> labels for every ELSET and NSET."""
    sets: dict[str, set[int]] = {}
    current: str | None = None
    for line in nam_path.read_text().splitlines():
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.startswith("*"):
            head = stripped.split(",")[1].strip() if "," in stripped else ""
            current = head.split("=", 1)[1].strip() if "=" in head else None
            if current:
                sets[current] = set()
            continue
        if current is None:
            continue
        sets[current].update(int(p) for p in stripped.split(",") if p.strip())
    return sets


def test_deck_is_id_offset_invariant(tmp_path: Path) -> None:
    """The exported .msh/.nam are byte-identical for 0-based and offset ids."""
    base_dir = tmp_path / "base"
    base_dir.mkdir()
    base = write_deck(base_dir, build_mesh(0), "m")
    shifted_dir = tmp_path / "shifted"
    shifted_dir.mkdir()
    shifted = write_deck(shifted_dir, build_mesh(OFFSET), "m")

    for suffix in (".msh", ".nam", ".inp"):
        base_text = base.with_suffix(suffix).read_text()
        shifted_text = shifted.with_suffix(suffix).read_text()
        assert base_text == shifted_text, (
            f"{suffix} differs between 0-based and offset ids; the writer is mixing "
            f"entity ids with 1-based indices"
        )


def test_deck_labels_are_internally_consistent(tmp_path: Path) -> None:
    """Every connectivity and set label must point at an entity that exists."""
    inp = write_deck(tmp_path, build_mesh(OFFSET), "m")
    node_labels, element_labels, connectivity = _parse_msh(inp.with_suffix(".msh"))
    named_sets = _parse_nam(inp.with_suffix(".nam"))

    referenced_nodes = {label for row in connectivity for label in row}
    assert referenced_nodes <= node_labels, (
        f"connectivity references absent nodes: {sorted(referenced_nodes - node_labels)}"
    )

    for set_name, labels in named_sets.items():
        if set_name.upper().startswith("N"):
            assert labels <= node_labels, f"{set_name} references absent nodes"
        else:
            assert labels <= element_labels, f"{set_name} references absent elements"

    # The section points at the element set, which must contain every element.
    plate = next(labels for name, labels in named_sets.items() if name.upper() == "EPLATE")
    assert plate == element_labels, (
        f"EPLATE={sorted(plate)} != all elements {sorted(element_labels)}"
    )


def test_offset_ids_give_the_same_ccx_result(tmp_path: Path) -> None:
    """Same model, two id schemes: CCX must return the same displacement."""
    ccx_bin = ccx_bin_or_skip()
    results = {}
    for name, offset in (("base", 0), ("shifted", OFFSET)):
        workdir = tmp_path / name
        workdir.mkdir()
        mesh = build_mesh(offset)
        inp = write_deck(workdir, mesh, "m")
        completed = run_ccx(inp, ccx_bin)
        if completed.returncode != 0:
            pytest.fail(f"CCX failed for {name}:\n{completed.stdout[-2000:]}")
        free_ids = sorted(
            mesh.node_id_to_index[n.id] + 1 for n in mesh.get_node_set("free_face").nodes.values()
        )
        disp = parse_frd_disp(inp.with_suffix(".frd"), free_ids)
        results[name] = np.mean([disp[nid][2] for nid in free_ids])

    assert results["base"] != 0.0
    rel = abs(results["base"] - results["shifted"]) / abs(results["base"])
    assert rel < 1e-9, (
        f"id scheme changed the result: base={results['base']}, shifted={results['shifted']}"
    )
