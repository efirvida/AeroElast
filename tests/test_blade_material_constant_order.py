"""Out-of-plane material constants reach the element in the element's own order.

windIO (and NuMAD's material dict) lists ``G = [G12, G13, G23]`` and ``nu = [nu12, nu13, nu23]``
with 1 the fibre axis. The Rust ``OrthotropicMaterial`` and the Python ``OrthotropicMaterial``
take ``(g12, g23, g13)`` and ``(nu12, nu23, nu31)``. Candidate D of the blade input audit
(``odd/tasks/blade-input-audit.md``, T5): the blade model passed the windIO order through, which
swapped g13/g23 (the transverse-shear moduli c44/c55) and put nu13 in the nu23 slot.
"""

import pytest

pytest.importorskip("_aeroelast")

from aeroelast.models.blade.model import build_rust_properties, material_factory  # noqa: E402

E = [40.0e9, 10.0e9, 9.0e9]
G = [4.0e9, 3.0e9, 2.0e9]  # windIO: G12, G13, G23
NU = [0.30, 0.25, 0.40]  # windIO: nu12, nu13, nu23
MAT = {"name": "uni", "density": 1900.0, "elastic": {"E": E, "nu": NU, "G": G}}


def _rust_material():
    data = {
        "materials": [MAT],
        "sections": [{"elementSet": "s", "layup": [["uni", 0.01, 0.0], ["uni", 0.02, 0.0]]}],
    }
    lam = build_rust_properties(data)["s"]
    return lam.plies_data()[0]["material"]


def test_rust_material_gets_element_order():
    m = _rust_material()
    assert m["g12"] == G[0]
    assert m["g13"] == G[1]
    assert m["g23"] == G[2]
    assert m["nu12"] == NU[0]
    assert m["nu23"] == NU[2]
    assert m["nu31"] == pytest.approx(NU[1] * E[2] / E[0], rel=1e-12)


def test_in_plane_constants_unchanged():
    m = _rust_material()
    assert (m["e1"], m["e2"], m["e3"], m["rho"]) == (E[0], E[1], E[2], 1900.0)


def test_python_material_gets_element_order():
    mat = material_factory(MAT)
    assert tuple(mat.G) == (G[0], G[2], G[1])
    assert mat.nu[0] == NU[0]
    assert mat.nu[1] == NU[2]
    assert mat.nu[2] == pytest.approx(NU[1] * E[2] / E[0], rel=1e-12)
