"""Cheap unit tests of tests/support/v2_sectional.py on a small eta patch of the windIO v2 blade (~10 s to generate)
and on a synthetic ring stack.  The full-blade S-1 comparison is tests/validation/blade/v2/test_v2_iea15mw_s1_sectional.py."""

import numpy as np
import pytest

pytest.importorskip("petsc4py", reason="PETSc not available")

from tests.support import v2_blade as vb  # noqa: E402
from tests.support import v2_sectional as vs  # noqa: E402

PATCH = (0.30, 0.32)


@pytest.fixture(scope="module")
def patch(tmp_path_factory):
    with pytest.MonkeyPatch.context() as mp:
        mp.setenv("V2_BLADE_CACHE_DIR", str(tmp_path_factory.mktemp("v2cache")))
        mesh, props = vb.build_v2_blade(0.25, eta_bounds=PATCH)
        ring = vs.ring_ids_for_v2_mesh(mesh)
        ext = vs.V2SectionalExtractor(mesh, props, ring_ids=ring)
        yield mesh, props, ring, ext


def test_rings_are_ordered_in_span_and_stations_are_one_bay(patch):
    mesh, _, ring, ext = patch
    z = mesh.coords_array[:, 2]
    zr = np.array([z[ring == r].mean() for r in range(ring.max() + 1)])
    assert ring.max() >= 5 and np.all(np.diff(zr) > 0)  # ring index grows with span
    mem = ext.station_membership()
    assert sorted(mem) == list(range(ring.max()))  # n_rings - 1 bays, keys consecutive
    for r, m in mem.items():
        assert set(ring[m["nodes"]]) == {r, r + 1}  # a bay holds exactly its two rings
    # the root ring is the RootNodes set
    root = sorted(ext.mesh.node_id_to_index[n] for n in mesh.get_node_set("RootNodes").nodes)
    assert set(root) == set(np.nonzero(ring == 0)[0])


def test_mass_partition_is_conservative(patch):
    """sum of bay masses = total lumped mass minus half the root ring (v1 partition: the last bay has no upper face, so
    the tip ring counts fully there; the first bay's lower face is the root ring, half)."""
    _, _, ring, ext = patch
    mem = ext.station_membership()
    total = 0.0
    for r in mem:
        sl = ext._station_slice(r, mem)
        total += ext.section_mass(r, mem)["mass_dens"] * max(sl["z1"] - sl["z0"], 1e-9)
    end = 0.5 * ext.node_mass[ring == 0].sum()
    assert total == pytest.approx(ext.node_mass.sum() - end, rel=1e-9)


def test_total_lumped_mass_matches_layup_integral(patch):
    """Independent: sum over elements of area * sum(rho_k h_k) from the layup tables (not the assembled M)."""
    mesh, _, _, ext = patch
    nd = mesh.numad_mesh_data
    rho = {m["name"]: m["density"] for m in nd["materials"]}
    areal = {s["elementSet"]: sum(rho[mat] * h for mat, h, _ in s["layup"]) for s in nd["sections"]}
    x = ext.coords
    idx = mesh.node_id_to_index
    expect = 0.0
    for name, sets in areal.items():
        for e in mesh.get_element_set(name).elements:
            p = x[[idx[n] for n in e.node_ids]]
            cross = np.cross(p[2] - p[0], p[-1] - p[1]) / 2 if len(p) == 4 else np.cross(p[1] - p[0], p[2] - p[0]) / 2
            expect += sets * np.linalg.norm(cross)
    assert ext.node_mass.sum() == pytest.approx(expect, rel=2e-3)


def test_extract_selects_nearest_bays_and_mass_path_only(patch):
    _, _, _, ext = patch
    mem = ext.station_membership()
    zs = np.array([mem[r]["z"] for r in sorted(mem)])
    ref = np.array([zs[1] + 0.01, zs[4] - 0.02])
    sec = ext.extract(stations=ref)
    assert np.allclose(sec.z, zs[[1, 4]])
    assert np.all(sec.mass_dens > 0) and np.all(np.isnan(sec.ea)) and np.all(np.isnan(sec.ei_edge))
    assert len(ext.extract().z) == len(zs)  # stations=None keeps all


def test_ring_ids_from_z_splits_a_synthetic_stack():
    rng = np.random.default_rng(0)
    z = np.repeat([0.0, 0.25, 0.5, 0.9], 6) + rng.uniform(-0.02, 0.02, 24)
    xyz = np.column_stack([rng.normal(size=24), rng.normal(size=24), z])
    ring = vs.ring_ids_from_z(xyz, tol=0.1)
    assert ring.tolist() == np.repeat([0, 1, 2, 3], 6).tolist()


def test_bay_centres_are_the_49_grid_bays():
    c = vs.bay_centres_v2()
    assert c.shape == (49,) and c[0] == pytest.approx(1.1938775510204083) and c[-1] == pytest.approx(115.80612244897959)
    assert np.allclose(np.diff(c), 117.0 / 49)


def test_bad_ring_ids_rejected(patch):
    mesh, props, ring, _ = patch
    with pytest.raises(ValueError, match="one entry per mesh node"):
        vs.V2SectionalExtractor(mesh, props, ring_ids=ring[:-1])


class _Elem:
    def __init__(self, id, node_ids):
        self.id, self.node_ids = id, node_ids


class _Mesh:
    def __init__(self, elements, n):
        self.elements, self.node_count = elements, n
        self.node_id_to_index = {i: i for i in range(n)}


def _synthetic(elem_nodes, z, ring):
    """A V2SectionalExtractor shell around a synthetic mesh (bypasses the heavy base-class setup)."""
    cls = vs.make_v2_extractor_class()
    ext = cls.__new__(cls)
    ext.mesh = _Mesh([_Elem(k, n) for k, n in enumerate(elem_nodes)], len(z))
    ext.coords = np.column_stack([np.zeros(len(z)), np.zeros(len(z)), z])
    ext.ring_ids = np.asarray(ring)
    return ext


def test_element_spanning_three_rings_goes_to_the_bay_of_its_centroid():
    # 4 rings of 2 nodes (z 0, 1, 2.5, 3); element 2 touches rings 0,1,2 with centroid z = 0.875 -> bay 0
    z = [0, 0, 1, 1, 2.5, 2.5, 3, 3]
    ring = [0, 0, 1, 1, 2, 2, 3, 3]
    elems = [[0, 1, 3, 2], [2, 3, 5, 4], [0, 1, 3, 4], [4, 5, 7, 6]]
    mem = _synthetic(elems, z, ring).station_membership()
    assert sorted(mem) == [0, 1, 2]
    assert set(mem[0]["nodes"]) == {0, 1, 2, 3, 4}  # the 3-ring element joined bay 0 (centroid 0.875 < 1)
    assert set(mem[1]["nodes"]) == {2, 3, 4, 5}
    # a 3-ring element with its centroid in the upper bay goes up instead
    mem2 = _synthetic([[0, 1, 3, 2], [2, 3, 5, 4], [1, 4, 5, 3], [4, 5, 7, 6]], z, ring).station_membership()
    assert set(mem2[1]["nodes"]) == {1, 2, 3, 4, 5}  # centroid z = 1.5 -> bay 1


def test_membership_is_a_partition_of_the_elements_and_deterministic():
    rng = np.random.default_rng(1)
    z = np.repeat([0.0, 1.0, 1.4, 3.0, 4.0], 4) + rng.uniform(-0.1, 0.1, 20)
    ring = np.repeat(np.arange(5), 4)
    elems = [list(rng.choice(20, 4, replace=False)) for _ in range(60)]  # arbitrary ring spans, incl. 1 / 3+ rings
    a = _synthetic(elems, z, ring).station_membership()
    b = _synthetic(elems, z, ring).station_membership()
    assert all(np.array_equal(a[k]["nodes"], b[k]["nodes"]) and a[k]["z"] == b[k]["z"] for k in a)
    zr = np.array([z[ring == r].mean() for r in range(5)])
    counted = 0
    for el in elems:
        bay = [k for k in a if set(el) <= set(a[k]["nodes"])]
        assert bay  # every element lies inside at least one bay's node set (it was assigned to exactly one)
        counted += 1
    assert counted == len(elems) and set(a) <= {0, 1, 2, 3} and zr[0] < zr[-1]


def test_plain_two_ring_elements_unchanged_by_the_centroid_rule(patch):
    """Plain-mesh behaviour: the bay of every element is its lower ring, bit for bit the pre-fix definition."""
    mesh, _, ring, ext = patch
    idx = mesh.node_id_to_index
    old: dict = {}
    for e in mesh.elements:
        nodes = [idx[n] for n in e.node_ids]
        rs = ring[nodes]
        assert rs.max() - rs.min() == 1
        old.setdefault(int(rs.min()), set()).update(nodes)
    new = ext.station_membership()
    assert sorted(new) == sorted(old)
    for r in old:
        arr = np.array(sorted(old[r]), dtype=np.int64)
        assert np.array_equal(new[r]["nodes"], arr)
        assert new[r]["z"] == float(ext.coords[arr, 2].mean())


def test_substitution_is_restored():
    import aeroelast.postprocess.sectional as sec

    before = sec.SectionalExtractor
    with vs.v2_extractor_substituted():
        assert sec.SectionalExtractor is not before
    assert sec.SectionalExtractor is before
