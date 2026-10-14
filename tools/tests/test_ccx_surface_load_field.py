"""The S8R deck accepts a fixed vector load at each surface node."""

import pytest

from tests.support.ccx_surface import node_id, write_s8r_surface_deck


def _write(path, **kwargs):
    return write_s8r_surface_deck(
        path,
        lambda u, v: (u, v, 0.0),
        1,
        1,
        young=2e11,
        poisson=0.3,
        thickness=0.08,
        boundary=lambda x, y, z: [],
        **kwargs,
    )


def test_vector_dead_load_emits_all_components_without_pressure(tmp_path):
    deck = tmp_path / "vector.inp"
    nid = node_id(0, 0, 1)
    _write(deck, load_field={nid: (1.25, -2.5, 0.0)}, nlgeom=4)
    text = deck.read_text()
    assert "*STEP, NLGEOM" in text
    assert "*DLOAD" not in text
    assert f"*CLOAD\n{nid}, 1, 1.250000000000e+00\n{nid}, 2, -2.500000000000e+00\n" in text
    assert f"{nid}, 3," not in text.split("*CLOAD\n", 1)[1]


def test_unknown_vector_load_node_is_rejected_before_writing(tmp_path):
    deck = tmp_path / "invalid.inp"
    with pytest.raises(ValueError, match="unknown load node: 9999"):
        _write(deck, load_field={9999: (1.0, 0.0, 0.0)})
    assert not deck.exists()


def test_existing_scalar_nodal_load_keeps_its_deck(tmp_path):
    deck = tmp_path / "scalar.inp"
    nid = node_id(0, 0, 1)
    _write(deck, nodal_loads={nid: (2, 3.0)})
    text = deck.read_text()
    assert "*STEP\n*STATIC" in text
    assert "*STEP, NLGEOM" not in text
    assert f"*CLOAD\n{nid}, 2, 3.000000000000e+00\n" in text
