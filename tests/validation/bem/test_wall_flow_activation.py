"""Guard for issue #16: the production coupling mesh cannot realise a ring section.

This is the activation guard for the BEM wall-flow moment realisation (roadmap
item P2 of #18).  It drives the **production** coupling-mesh builder
``aeroelast.cli.run_bem_fsi._build_mesh`` with
``coupling_node_set: allOuterShellNods`` - the node set every campaign and
smoke YAML uses - and then builds the production ``ForceProjector`` with the
official deck's element property map, exactly as the activated path will.

Measured today the coupling mesh is **node-only**: ``_build_mesh`` filters the
generator mesh to ``MeshModel(nodes=filtered_nodes)`` and drops every element
(``src/aeroelast/cli/run_bem_fsi.py``).  With no elements, ``ring_section``
finds no chordwise wall edge and every entry of
``ForceProjector._strip_ring_sections`` is ``None``, so
``_ring_section_is_realisable`` (``src/aeroelast/solvers/bem/force_projection.py``,
around line 1036) refuses the strip and it falls back to the minimum-norm
``_distribute`` field - the field ``moment_realization_over_delivers`` says
over-delivers rigid torsion.

The test records the two counts that authorise the change:

* ``len(mesh.elements)`` - the filtered coupling mesh's element count
  (**0** today), and
* ``usable of total ring sections usable`` - how many of the projector's
  precomputed physical rings pass the realisability gate (**0 today**).
"""

from __future__ import annotations

import pytest

pytest.importorskip("ccblade", reason="ccblade not installed (pip install -e '.[bem]')")
pytest.importorskip("_aeroelast", reason="Rust backend not available")

from aeroelast.cli.run_bem_fsi import _build_mesh  # noqa: E402
from aeroelast.models.blade.model import Blade  # noqa: E402
from aeroelast.solvers.bem.force_projection import ForceProjector  # noqa: E402

from tests.support.openfast_bem import build_blade_aero_from_aerodyn  # noqa: E402
from tests.support.paths import DATA_DIR  # noqa: E402

YAML = DATA_DIR / "IEA-15-240-RWT.yaml"
AIRFOILS = DATA_DIR / "airfoils"
AD_PRIMARY = DATA_DIR / "reference" / "iea15mw_openfast" / "case" / "IEA-15-240-RWT_AeroDyn15.dat"

#: Deliberately coarse: the point is the element count, not the campaign mesh.
#: ``element_size: 0.25`` on this deck is ~27 609 coupling nodes; 1.0 keeps the
#: test well under its time budget while still producing a real blade mesh.
ELEMENT_SIZE = 1.0
N_SAMPLES = 150

SPAN_DIR = [0.0, 0.0, 1.0]
NODE_SET = "allOuterShellNods"


def _production_cfg() -> dict:
    """The coupling-mesh section of the fluid YAML every campaign and smoke case uses."""
    return {
        "mesh": {
            "source": "generator",
            "generator": {
                "type": "BladeMesh",
                "params": {
                    "yaml_file": str(YAML),
                    "airfoil_dir": str(AIRFOILS),
                    "element_size": ELEMENT_SIZE,
                    "n_samples": N_SAMPLES,
                    "airfoil_spacing": "constant",
                    "span_grading": "chord",
                },
                "coupling_node_set": NODE_SET,
            },
        }
    }


@pytest.fixture(scope="module")
def production():
    """The production coupling mesh, the deck's property map and the real BladeAero."""
    for path in (YAML, AIRFOILS, AD_PRIMARY):
        if not path.exists():
            pytest.skip(f"required input not present: {path}")

    built = _build_mesh(_production_cfg(), YAML)
    # arity-stable: a later task may add the property map as a third value
    mesh, viz_mesh = built[0], built[1]
    assert viz_mesh is not None, "a filtered coupling mesh must keep the full mesh as viz_mesh"

    # The element property map and the mesh, built the way the generator builds them.
    model = Blade(str(YAML), element_size=ELEMENT_SIZE)
    model.generate_mesh()

    return {
        "mesh": mesh,
        "viz_mesh": viz_mesh,
        "node_set": viz_mesh.get_node_set(NODE_SET),
        "props": model.get_element_properties(),
        "blade_aero": build_blade_aero_from_aerodyn(AD_PRIMARY),
    }


def test_coupling_node_set_filter_keeps_the_node_list_and_the_skin_elements(production):
    """The filter is exactly the node set's nodes, in its order, plus the skin elements.

    This pins the #16 wiring: the previous filter returned
    ``MeshModel(nodes=list(ns.nodes.values()))``, so the node list and its order
    are the coupling vertex order and must not move, while the fully contained
    elements must survive or no ring section is realisable.
    """
    mesh = production["mesh"]
    viz_mesh = production["viz_mesh"]
    requested = list(production["node_set"].nodes.keys())

    assert [node.id for node in mesh.nodes] == requested
    assert [tuple(node.coords.tolist()) for node in mesh.nodes] == [
        tuple(viz_mesh.node_map[node_id].coords.tolist()) for node_id in requested
    ]

    requested_set = set(requested)
    assert mesh.elements, "the coupling mesh kept no element at all"
    for element in mesh.elements:
        assert set(element.node_ids) <= requested_set, (
            f"element {element.id} drags a node outside the coupling node set"
        )

    # The shear-web elements are the ones the filter leaves behind, and the full
    # mesh survives as viz_mesh for the surface VTU output.
    assert len(mesh.elements) < len(viz_mesh.elements)
    assert viz_mesh.get_node_set(NODE_SET).node_ids == production["node_set"].node_ids

    for name, element_set in mesh.element_sets.items():
        assert name in viz_mesh.element_sets, f"element set '{name}' is not the source mesh's"
        assert {element.id for element in element_set.elements} <= {
            element.id for element in mesh.elements
        }, f"element set '{name}' claims an element the filter dropped"


def test_production_coupling_mesh_realises_a_ring_section(production):
    """The production coupling mesh must keep its skin elements so a ring is realisable.

    Before the fix ``_build_mesh`` returned a nodes-only coupling mesh, so the
    first assertion failed with the measured ``0`` element count (0 of 185 ring
    sections usable at this mesh size).  It quotes the two counts the feature
    note records so a green run also carries the evidence.
    """
    mesh = production["mesh"]
    proj = ForceProjector(
        mesh,
        production["blade_aero"],
        span_direction=SPAN_DIR,
        element_properties=production["props"],
    )

    # One entry per physical ring of every strip; ``None`` means the realisability
    # gate refused it.  Flatten the per-strip lists into the ring-section count the
    # feature note records ("0 of 671 ring sections usable").
    sections = [s for per_strip in proj._strip_ring_sections for s in per_strip]
    usable = sum(1 for s in sections if s is not None)
    total = len(sections)
    counts = f"{usable} of {total} ring sections usable"

    assert len(mesh.elements) > 0, (
        "the production coupling mesh dropped every element: "
        f"len(mesh.elements) = {len(mesh.elements)}, "
        f"len(mesh.nodes) = {len(mesh.nodes)}, {counts}; "
        "_build_mesh must carry the skin elements so ring_section can find a "
        "closed cell and _ring_section_is_realisable stops falling back to "
        "_distribute"
    )
    assert total > 0, (
        f"the projector precomputed no ring section at all: "
        f"len(mesh.elements) = {len(mesh.elements)}, {counts}"
    )
    assert usable > 0, (
        "no ring section passes _ring_section_is_realisable: "
        f"len(mesh.elements) = {len(mesh.elements)}, {counts}"
    )
