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
    # unstructured ring labelling are not rings and take only the minimum-norm residual; a real
    # (>= 3 node) ring with no usable section is skipped and takes only its minimum-norm share
    per_strip = proj._strip_wall_flow_sections
    sections = [s for strip_sections in per_strip for s in strip_sections]
    stubs = sum(len(g) < 3 for strip_groups in proj._strip_ring_groups for g in strip_groups)
    usable = sum(1 for s in sections if s is not None)
    skipped = len(sections) - usable
    fallback = [k for k, strip_sections in enumerate(per_strip) if strip_sections and not any(strip_sections)]
    print(
        f"[gmsh-coupling] {len(coupling['mesh'].nodes)} nodes, {usable} of {len(sections)} ring sections usable, "
        f"{skipped} unrealisable rings skipped, {stubs} stub groups excluded, "
        f"{len(fallback)} strips with no realisable ring (minimum-norm fallback): {fallback}"
    )
    assert sections, "the projector precomputed no ring section"
    # strict about strips: every strip that has a ring must realise the wall flow on at least one of them
    assert not fallback, f"strips {fallback} have no realisable ring and fall back to the minimum-norm field"


# --- mesh-size keys reach GmshBladeMesh (L5 T6): forwarded only when present, defaults untouched ---------------------------

class _Stop(Exception):
    pass


class _StubGmsh:
    """Records the constructor kwargs; ``generate`` stops the caller before any meshing."""

    calls: list = []

    def __init__(self, **kwargs):
        type(self).calls.append(kwargs)

    def generate(self, *args, **kwargs):
        raise _Stop



_MESH_KEYS = {"chord_divisions": 20, "span_aspect": 1.0, "h_range": [0.05, 0.25], "min_band_elements": 2,
              "min_band_fraction": 1.0}


def _cli_kwargs(monkeypatch, tmp_path, params):
    import aeroelast.core.mesh.gmsh_blade as gb

    _StubGmsh.calls = []
    monkeypatch.setattr(gb, "GmshBladeMesh", _StubGmsh)
    cfg = {"mesh": {"source": "generator", "generator": {"type": "GmshBladeMesh", "params": params}}}
    with pytest.raises(_Stop):
        _build_mesh(cfg, tmp_path / "case.yaml")
    assert len(_StubGmsh.calls) == 1
    return _StubGmsh.calls[0]


def _runner_kwargs(monkeypatch, params):
    from types import SimpleNamespace

    import aeroelast.solvers.fsi.runner as runner

    _StubGmsh.calls = []
    monkeypatch.setattr(runner, "GmshBladeMesh", _StubGmsh)
    runner.FSIRunner._gmsh_blade_generator(SimpleNamespace(config_path=None), params)
    assert len(_StubGmsh.calls) == 1
    return _StubGmsh.calls[0]


def test_cli_forwards_chord_sizing_keys_to_gmsh_blade_mesh(monkeypatch, tmp_path):
    kw = _cli_kwargs(monkeypatch, tmp_path, {"yaml_file": "/abs/blade.yaml", "skin_reference": "mid", **_MESH_KEYS})
    assert kw["chord_divisions"] == 20
    assert kw["span_aspect"] == 1.0
    assert tuple(kw["h_range"]) == (0.05, 0.25)
    assert kw["min_band_elements"] == 2
    assert kw["min_band_fraction"] == 1.0


def test_runner_forwards_chord_sizing_keys_to_gmsh_blade_mesh(monkeypatch):
    kw = _runner_kwargs(monkeypatch, {"yaml_file": "/abs/blade.yaml", "skin_reference": "mid", **_MESH_KEYS})
    assert kw["chord_divisions"] == 20
    assert kw["span_aspect"] == 1.0
    assert tuple(kw["h_range"]) == (0.05, 0.25)
    assert kw["min_band_elements"] == 2
    assert kw["min_band_fraction"] == 1.0


def test_cli_without_chord_keys_keeps_the_previous_kwargs(monkeypatch, tmp_path):
    kw = _cli_kwargs(monkeypatch, tmp_path, {"yaml_file": "/abs/blade.yaml"})
    assert kw == {"yaml_file": "/abs/blade.yaml", "element_size": 0.25, "skin_reference": "oml"}


def test_runner_without_chord_keys_keeps_the_previous_kwargs(monkeypatch):
    kw = _runner_kwargs(monkeypatch, {"yaml_file": "/abs/blade.yaml"})
    assert kw == {"yaml_file": "/abs/blade.yaml", "element_size": 0.25, "skin_reference": "mid"}
