"""The gmsh blade mesh as the BEM-FSI coupling mesh (odd/tasks/gmsh-official-mesher.md, T1 / G-C).

Drives the production coupling-mesh builder ``aeroelast.cli.run_bem_fsi._build_mesh`` with
``type: GmshBladeMesh`` on the windIO 2.0 IEA-15 definition (``coupling_node_set: allOuterShellNods``) and builds the
production ``ForceProjector`` with the deck's element property map.  The projector groups the coupling nodes into physical
rings and needs a realisable wall-flow section on each; the NuMAD mesh has constant-span rings by construction, the gmsh
skin mesh does not, so this is the integration question the one-way sensitivity tool cannot answer.
"""

from __future__ import annotations

import pytest

pytest.importorskip("ccblade", reason="ccblade not installed (pip install -e '.[bem]')")
pytest.importorskip("_aeroelast", reason="Rust backend not available")

from aeroelast.cli.run_bem_fsi import _build_mesh  # noqa: E402
from aeroelast.solvers.bem.force_projection import ForceProjector  # noqa: E402

from tests.support.openfast_bem import build_blade_aero_from_aerodyn  # noqa: E402
from tests.support.paths import DATA_DIR  # noqa: E402

YAML = DATA_DIR / "IEA15MW" / "reference" / "upstream" / "IEA-15-240-RWT.yaml"
AD_PRIMARY = DATA_DIR / "reference" / "iea15mw_openfast" / "case" / "IEA-15-240-RWT_AeroDyn15.dat"
NODE_SET = "allOuterShellNods"


@pytest.fixture(scope="module")
def coupling():
    for path in (YAML, AD_PRIMARY):
        if not path.exists():
            pytest.skip(f"required input not present: {path}")
    cfg = {
        "mesh": {
            "source": "generator",
            "generator": {
                "type": "GmshBladeMesh",
                "params": {"yaml_file": str(YAML), "element_size": 0.25},
                "coupling_node_set": NODE_SET,
            },
        }
    }
    built = _build_mesh(cfg, YAML)
    return {"mesh": built[0], "viz_mesh": built[1], "props": built[2] if len(built) > 2 else None}


def test_gmsh_generator_builds_the_coupling_mesh_with_deck_properties(coupling):
    mesh, viz = coupling["mesh"], coupling["viz_mesh"]
    assert viz is not None and len(mesh.nodes) > 0 and len(mesh.elements) > 0
    assert len(mesh.elements) < len(viz.elements)  # the web elements are left behind
    assert coupling["props"], "the generator's deck data did not give an element property map"


def test_projector_rings_are_realisable_on_the_gmsh_mesh(coupling):
    proj = ForceProjector(
        coupling["mesh"],
        build_blade_aero_from_aerodyn(AD_PRIMARY),
        span_direction=[0.0, 0.0, 1.0],
        element_properties=coupling["props"],
    )
    # the sections project() realises the wall flow on: one- and two-node stubs of the
    # unstructured ring labelling are not rings and take only the minimum-norm residual
    sections = [s for per_strip in proj._strip_wall_flow_sections for s in per_strip]
    stubs = sum(len(g) < 3 for per_strip in proj._strip_ring_groups for g in per_strip)
    usable = sum(1 for s in sections if s is not None)
    print(
        f"[gmsh-coupling] {len(coupling['mesh'].nodes)} nodes, {usable} of {len(sections)} ring sections usable, "
        f"{stubs} stub groups excluded"
    )
    assert sections, "the projector precomputed no ring section"
    assert usable == len(sections), f"{len(sections) - usable} of {len(sections)} ring sections are unrealisable"
