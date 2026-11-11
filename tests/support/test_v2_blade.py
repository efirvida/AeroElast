"""Unit tests of tests/support/v2_blade.py (cache key, shape contract on a small eta patch, v1-consumer contract on
the cached full mesh). The full-mesh test skips when the es 0.25 cache has not been produced (sbatch, ~10 min)."""

import numpy as np
import pytest

pytest.importorskip("petsc4py", reason="PETSc not available")

from tests.support import v2_blade as vb  # noqa: E402

PATCH = (0.30, 0.32)  # a few bays, ~9 s to generate


def test_cache_key_depends_on_arguments_and_source(tmp_path, monkeypatch):
    base = vb.cache_key(0.25)
    assert base == vb.cache_key(0.25)
    assert base != vb.cache_key(0.2)
    assert base != vb.cache_key(0.25, collapse_ratio=0.15)
    assert base != vb.cache_key(0.25, eta_bounds=PATCH)
    assert base != vb.cache_key(0.25, min_strip_m=0.0)
    # YAML bytes
    other = tmp_path / "other.yaml"
    other.write_bytes(vb.DEFAULT_YAML.read_bytes() + b"\n# x\n")
    assert base != vb.cache_key(0.25, yaml_file=other)
    # mesher source digest
    monkeypatch.setattr(vb, "source_digest", lambda: "changed")
    assert base != vb.cache_key(0.25)


def test_element_size_range_enforced():
    for es in (0.1, 0.3, 0.5):
        with pytest.raises(ValueError, match="outside"):
            vb.build_v2_blade(es)


def test_bad_skin_reference_rejected():
    with pytest.raises(ValueError, match="skin_reference"):
        vb.mesh_from_raw({}, "bogus")


@pytest.fixture(scope="module")
def patch(tmp_path_factory):
    with pytest.MonkeyPatch.context() as mp:  # scoped: the full-mesh test must see the real cache dir
        mp.setenv("V2_BLADE_CACHE_DIR", str(tmp_path_factory.mktemp("v2cache")))
        yield (vb.build_v2_blade(0.25, eta_bounds=PATCH),
               vb.build_v2_blade(0.25, skin_reference="oml", eta_bounds=PATCH))


def test_patch_shape_contract(patch):
    (mesh, props), (mesh_oml, _) = patch
    names = set(mesh.element_sets_names)
    sec = sorted(n for n in names if n.startswith("SEC"))
    assert sec and set(sec) == set(props)  # properties keyed by SEC{i}, one per element set
    assert {"RootNodes", "TipNodes", "allOuterShellNods", "allOuterShellEls"} <= names | set(mesh.node_sets)
    root = mesh.get_node_set("RootNodes").nodes
    tip = mesh.get_node_set("TipNodes").nodes
    assert root and tip and not (set(root) & set(tip))
    z = mesh.coords_array[:, 2]
    idx = mesh.node_id_to_index
    assert z[[idx[n] for n in root]].max() < z[[idx[n] for n in tip]].min()
    assert mesh.coords_array.shape == (mesh.node_count, 3)
    assert all(hasattr(mesh, a) for a in ("numad_mesh_data", "v2_info"))
    assert {s["elementSet"] for s in mesh.numad_mesh_data["sections"]} == set(sec)
    # every element belongs to exactly one SEC set
    assert sum(len(mesh.get_element_set(n).elements) for n in sec) == mesh.elements_count


def test_mid_reference_moves_skin_inward_only(patch):
    (mid, _), (oml, _) = patch
    assert mid.coords_array.shape == oml.coords_array.shape
    # same rcm numbering on both (same connectivity), so rows are comparable
    d = np.linalg.norm(mid.coords_array - oml.coords_array, axis=1)
    skin_ids = set(mid.get_node_set("allOuterShellNods").nodes)
    skin = np.zeros(mid.node_count, dtype=bool)
    skin[[mid.node_id_to_index[i] for i in skin_ids]] = True
    assert d[skin].max() > 1e-3  # laminates are mm thick, not zero
    assert d[skin].max() < 0.1  # and nowhere near a web height
    web_only = ~skin
    if web_only.any():
        assert d[web_only].max() == 0.0


def test_memo_returns_same_objects(patch):
    (mesh, props), _ = patch
    again = vb.build_v2_blade(0.25, eta_bounds=PATCH)
    assert again[0] is mesh and again[1] is props


def test_full_cached_mesh_has_what_v1_consumers_read(monkeypatch):
    monkeypatch.delenv("V2_BLADE_CACHE_DIR", raising=False)  # the module-scoped patch fixture may still be active
    try:
        raw = vb.load_raw(0.25, generate=False)
    except FileNotFoundError:
        pytest.skip("es 0.25 full-blade cache absent (produce with the sbatch job)")
    mesh, props = vb.mesh_from_raw(raw, "mid")
    # V-02/S-2/S-4/S-1: RootNodes clamp, coords_array, properties keyed by element-set name
    root = sorted(mesh.get_node_set("RootNodes").nodes.keys())
    idx = mesh.node_id_to_index
    assert len(root) > 50 and all(isinstance(idx[n], (int, np.integer)) for n in root)
    z = mesh.coords_array[:, 2]
    assert z[[idx[n] for n in root]].max() < z.min() + 0.2
    tip = np.nonzero(z >= z.max() - 0.2)[0]
    assert len(tip) > 20  # S-2/S-7 tip ring found from coords_array
    assert set(props) == {n for n in mesh.element_sets_names if n.startswith("SEC")}
    for name, lam in props.items():
        assert isinstance(lam, dict) or hasattr(lam, "total_thickness")
    # blade length and mass of the right order (official blade: 117 m span, ~67 t BeamDyn)
    assert 116.0 < z.max() - z.min() < 118.0
    from aeroelast.models.blade.model import build_rust_properties  # noqa: F401  (properties came from it)

    assert set(mesh.numad_mesh_data) >= {"materials", "sections"}


def test_mesher_is_part_of_the_cache_key_and_gmsh_digest_covers_its_source(monkeypatch):
    plain = vb.cache_key(0.25)
    assert plain == vb.cache_key(0.25, mesher="windio2")
    assert plain != vb.cache_key(0.25, mesher="gmsh")
    gmsh_key = vb.cache_key(0.25, mesher="gmsh")
    assert gmsh_key == vb.cache_key(0.25, mesher="gmsh")
    assert vb.source_digest("gmsh") != vb.source_digest("windio2")
    # editing gmsh_skin.py invalidates the gmsh key only
    real = vb.source_digest
    monkeypatch.setattr(vb, "source_digest", lambda mesher="windio2": real(mesher) + ("x" if mesher == "gmsh" else ""))
    assert vb.cache_key(0.25) == plain
    assert vb.cache_key(0.25, mesher="gmsh") != gmsh_key


def test_gmsh_mesher_rejects_collapse_and_partial_blade_and_unknown_mesher():
    with pytest.raises(ValueError, match="no collapse"):
        vb.load_raw(0.25, 0.15, mesher="gmsh", generate=False)
    with pytest.raises(ValueError, match="no collapse"):
        vb.load_raw(0.25, None, PATCH, mesher="gmsh", generate=False)
    with pytest.raises(ValueError, match="mesher"):
        vb.load_raw(0.25, mesher="bogus", generate=False)


def test_mesher_environment_variable_selects_gmsh(monkeypatch):
    monkeypatch.setenv("V2_BLADE_MESHER", "gmsh")
    assert vb._resolve_mesher(None) == "gmsh"
    assert vb._resolve_mesher("windio2") == "windio2"
    monkeypatch.delenv("V2_BLADE_MESHER")
    assert vb._resolve_mesher(None) == "windio2"


def test_element_size_env_override_wins_and_is_range_checked(monkeypatch):
    """``V2_BLADE_ELEMENT_SIZE`` replaces the size a test hard-codes (campaign convergence runs); unset = unchanged."""
    from tests.support import v2_blade

    monkeypatch.delenv("V2_BLADE_ELEMENT_SIZE", raising=False)
    assert v2_blade._resolve_element_size(0.25) == 0.25
    monkeypatch.setenv("V2_BLADE_ELEMENT_SIZE", "0.177")
    assert v2_blade._resolve_element_size(0.25) == 0.177
    monkeypatch.setenv("V2_BLADE_ELEMENT_SIZE", "0.5")
    with pytest.raises(ValueError):
        v2_blade.load_raw(0.25, generate=False)
