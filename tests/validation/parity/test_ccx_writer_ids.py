"""The CCX result must not depend on the entity-id scheme.

The entity id counters (``Node._id_counter``, ``MeshElement._id_counter``) are
process-global and never reset by ``MeshModel`` itself, so a mesh built after
another one has node and element ids that are neither 0-based nor contiguous.
This is a physical parity check: the same 2x2 quad plate is exported with a
0-based and an arbitrarily offset id scheme, both decks are solved by
CalculiX, and the free-face displacement must match. The comparison is against
a different code (CCX), not against our own writer.

The writer's own deck/id-scheme contract is checked separately under
``tests/software/contracts/test_ccx_writer_ids.py``.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from tests.conftest import ccx_bin_or_skip

pytest.importorskip("_aeroelast", reason="Rust backend not available")

from tests.support.ccx_io import parse_frd_disp, run_ccx  # noqa: E402
from tests.support.ccx_plate import OFFSET, build_mesh, write_deck  # noqa: E402


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
