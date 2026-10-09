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


def test_production_coupling_mesh_realises_a_ring_section():
    """The production coupling mesh must keep its skin elements so a ring is realisable.

    Today ``_build_mesh`` returns a nodes-only coupling mesh, so the first
    assertion fails with the measured ``0`` element count.  It quotes the two
    counts the feature note records so a green run also carries the evidence.
    """
    for path in (YAML, AIRFOILS, AD_PRIMARY):
        if not path.exists():
            pytest.skip(f"required input not present: {path}")

    cfg = {
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
                "coupling_node_set": "allOuterShellNods",
            },
        }
    }
    built = _build_mesh(cfg, YAML)
    # arity-stable: a later task may add a third value
    mesh, viz_mesh = built[0], built[1]  # noqa: F841

    # The element property map, built the way the mesh generator builds it.
    model = Blade(str(YAML), element_size=ELEMENT_SIZE)
    model.generate_mesh()
    props = model.get_element_properties()
    blade_aero = build_blade_aero_from_aerodyn(AD_PRIMARY)

    proj = ForceProjector(
        mesh,
        blade_aero,
        span_direction=SPAN_DIR,
        element_properties=props,
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
