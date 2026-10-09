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

import numpy as np
import pytest

pytest.importorskip("ccblade", reason="ccblade not installed (pip install -e '.[bem]')")
pytest.importorskip("_aeroelast", reason="Rust backend not available")

from aeroelast.cli.run_bem_fsi import _build_mesh  # noqa: E402
from aeroelast.core.mesh.model import MeshModel  # noqa: E402
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


def test_factory_participant_realises_every_precomputed_ring(production):
    """The production factory must hand the deck's properties to the projector.

    The CLI now builds the deck's section-property map and forwards it through
    ``build_from_config`` to ``ForceProjector``; with it every precomputed ring
    of the element-bearing coupling mesh must be realisable, which is what makes
    the wall-flow realisation run in production (issue #16).
    """
    from aeroelast.solvers.bem.fsi_participant import build_from_config  # noqa: PLC0415

    cfg = {
        "blade_file": str(YAML),
        "bem": {
            "wind_speed": 10.59,
            "omega": 7.56 * 2.0 * np.pi / 60.0,
            "pitch": 0.0,
            "air_density": 1.225,
            "dynamic_viscosity": 1.81206e-5,
            "hub_height": 150.0,
            "shear_exp": 0.0,
            "default_re": 1e7,
        },
    }
    participant = build_from_config(
        production["mesh"],
        cfg,
        viz_mesh=production["viz_mesh"],
        element_properties=production["props"],
    )

    sections = [s for per_strip in participant._projector._strip_ring_sections for s in per_strip]
    total = len(sections)
    assert total > 0, "the factory-built projector precomputed no ring section at all"

    unusable = [i for i, section in enumerate(sections) if section is None]
    assert not unusable, (
        f"{len(unusable)} of {total} ring sections are unrealisable through the production "
        f"factory: indices {unusable[:5]}"
    )
    without_properties = [
        i
        for i, section in enumerate(sections)
        if section is not None and not section.from_element_properties
    ]
    assert not without_properties, (
        f"{len(without_properties)} of {total} ring sections fell back to the geometric "
        f"S = 1.0 split instead of the deck's laminates: indices {without_properties[:5]}"
    )


def test_wall_flow_leaves_the_applied_resultant_unchanged(production):
    """The invariance guard: the wall-flow split must not move the resultant.

    The wall-flow realisation changes the *distribution* of the applied load, not
    its resultant - ``moment_realization_over_delivers`` measured that the realised
    force and moment are exact.  This pins it: ``sum(F)`` and the total moment
    ``sum(r_j x f_j)`` about the origin must agree, relative to the wall-flow
    values, to ``1e-9``.  A larger difference is a finding, never a bound to widen.
    """
    from aeroelast.solvers.bem.engine import BEMSolver  # noqa: PLC0415

    mesh = production["mesh"]
    blade_aero = production["blade_aero"]
    result = BEMSolver(
        blade_aero,
        rho=1.225,
        mu=1.81206e-5,
        hub_height=150.0,
        shear_exp=0.0,
    ).compute(10.59, 7.56, 0.0)

    wall = ForceProjector(
        mesh,
        blade_aero,
        span_direction=SPAN_DIR,
        element_properties=production["props"],
    )
    fallback = ForceProjector(mesh, blade_aero, span_direction=SPAN_DIR)
    forces_wall = wall.project(result)
    forces_fallback = fallback.project(result)

    coords = mesh.coords_array

    def resultant(forces):
        return forces.sum(axis=0), np.cross(coords, forces).sum(axis=0)

    force_wall, moment_wall = resultant(forces_wall)
    force_fallback, moment_fallback = resultant(forces_fallback)
    delta_force = float(np.linalg.norm(force_wall - force_fallback))
    delta_moment = float(np.linalg.norm(moment_wall - moment_fallback))
    scale_force = float(np.linalg.norm(force_wall))
    scale_moment = float(np.linalg.norm(moment_wall))
    relative_force = delta_force / scale_force
    relative_moment = delta_moment / scale_moment

    assert relative_force <= 1e-9, (
        "the wall-flow split moved the applied resultant force: "
        f"relative difference {relative_force:.3e} "
        f"(|dF| = {delta_force:.6e} N, |F| = {scale_force:.6e} N)"
    )
    assert relative_moment <= 1e-9, (
        "the wall-flow split moved the applied resultant moment: "
        f"relative difference {relative_moment:.3e} "
        f"(|dM| = {delta_moment:.6e} N.m, |M| = {scale_moment:.6e} N.m)"
    )


def test_without_the_property_map_every_strip_falls_back(production):
    """A property-less projector must keep the pre-#16 minimum-norm fallback.

    The element-bearing coupling mesh makes ``ring_section`` return a geometric
    (``S = 1.0``) split even without properties, so those entries are not ``None``;
    it is ``from_element_properties`` that marks them non-physical and sends every
    strip to ``_distribute``.  The check that matters is therefore the observable
    fallback: the produced forces must equal the pre-#16 node-only projector's,
    which by construction has no realisable ring at all.
    """
    mesh = production["mesh"]
    blade_aero = production["blade_aero"]

    without_properties = ForceProjector(mesh, blade_aero, span_direction=SPAN_DIR)
    sections = [s for per_strip in without_properties._strip_ring_sections for s in per_strip]
    assert sections, "the coupling mesh carries no ring section to refuse"
    assert all(section is not None for section in sections), (
        "the element-bearing coupling mesh should still extract a geometric section; "
        "a None here means ring_section failed for another reason"
    )
    assert all(
        section is not None and not section.from_element_properties for section in sections
    ), "a property-less projector must mark every ring geometric-only"

    # The pre-#16 geometry was node-only, so its entries were ``None``; reproducing
    # it here shows the element-bearing property-less path falls back identically.
    node_only = ForceProjector(
        MeshModel(nodes=list(mesh.nodes)), blade_aero, span_direction=SPAN_DIR
    )
    node_only_sections = [
        s for per_strip in node_only._strip_ring_sections for s in per_strip
    ]
    assert node_only_sections and all(section is None for section in node_only_sections)

    from aeroelast.solvers.bem.engine import BEMSolver  # noqa: PLC0415

    result = BEMSolver(
        blade_aero,
        rho=1.225,
        mu=1.81206e-5,
        hub_height=150.0,
        shear_exp=0.0,
    ).compute(10.59, 7.56, 0.0)
    forces = without_properties.project(result)
    forces_pre_16 = node_only.project(result)
    assert np.array_equal(forces, forces_pre_16), (
        "the property-less element-bearing projector left the minimum-norm fallback: "
        f"max |dF| = {np.max(np.abs(forces - forces_pre_16)):.6e} N"
    )
