"""Guard for the composite laminate's public surface (SDD task 12.2).

The composite transverse shear is the one place where this repository deliberately
departs from the "no numerical factor" rule of the element formulation: the shell
element consumes the UNCORRECTED plain section integral of the shear stiffness, and
the classical ``5/6`` factor stays in the material layer where a user can see and
change it (ADRs in ``openspec/changes/archive/2026-09-29-mitc4plusd-faithful/design.md``). That split
only holds while the surface carrying it keeps its shape, which is what this test
pins:

* the Rust ``ShellConstitutive`` field set -- what the element reads;
* the Rust ``Laminate`` method set -- what produces it, including the uncorrected
  section integral ``shear_stiffness_uncorrected`` and ``applied_shear_correction_factor``;
* the Python ``Laminate``/``Ply`` public surface -- what a user calls, which is a
  DIFFERENT set from the Rust type of the same name (``get_ABD_matrix``, ``n_plies``,
  ``shear_correction_factor`` and the symmetric/balanced predicates are Python-only).

A change here is not automatically a bug: it is a contract change in the composite
path, and it should be made on purpose, with this list updated in the same commit.
"""

from __future__ import annotations

import re
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
MATERIALS_MOD = REPO / "crates" / "aeroelast-core" / "src" / "materials" / "mod.rs"
LAMINATE_RS = REPO / "crates" / "aeroelast-core" / "src" / "materials" / "laminate.rs"

#: Fields the element reads off `ShellConstitutive`.
SHELL_CONSTITUTIVE_FIELDS = frozenset({"cm", "cb_coupling", "cb", "cs", "cm_raw"})

#: Methods `impl Laminate` must keep exposing (Rust side).
LAMINATE_RS_METHODS = frozenset({
    "new",
    "abd_matrix_flat",
    "applied_shear_correction_factor",
    "shear_stiffness_uncorrected",
    "to_shell_constitutive",
    "to_shell_constitutive_with_offset",
    "is_symmetric",
    "is_balanced",
})

#: The user-facing surface of the Python wrapper (`aeroelast.core.laminate`).
LAMINATE_PY_PUBLIC = frozenset({
    "get_ABD_matrix",
    "get_equivalent_properties",
    "is_balanced",
    "is_symmetric",
    "n_plies",
    "shear_correction_factor",
    "shear_stiffness_uncorrected",
    "total_thickness",
})
PLY_PY_PUBLIC = frozenset({"angle_rad", "strength", "z_bottom", "z_top"})


def test_laminate_public_surface_unchanged():
    """The composite contract's three surfaces keep their exact shape."""
    mod_src = MATERIALS_MOD.read_text()
    struct = mod_src.split("pub struct ShellConstitutive", 1)
    assert len(struct) == 2, "ShellConstitutive is no longer defined in materials/mod.rs"
    body = struct[1].split("}", 1)[0]
    fields = set(re.findall(r"pub ([a-z_][a-z0-9_]*)\s*:", body))
    assert fields == SHELL_CONSTITUTIVE_FIELDS, (
        "the ShellConstitutive field set changed; the element reads these directly: "
        f"missing={sorted(SHELL_CONSTITUTIVE_FIELDS - fields)} unexpected={sorted(fields - SHELL_CONSTITUTIVE_FIELDS)}"
    )
    for field in sorted(SHELL_CONSTITUTIVE_FIELDS):
        assert f"pub {field}:" in body, f"ShellConstitutive lost the field {field!r}"

    lam_src = LAMINATE_RS.read_text()
    impl = lam_src.split("impl Laminate", 1)
    assert len(impl) == 2, "impl Laminate is no longer in laminate.rs"
    methods = set(re.findall(r"pub fn ([a-z_][a-z0-9_]*)", impl[1]))
    assert methods == LAMINATE_RS_METHODS, (
        "the Rust Laminate method set changed, and the composite shear contract is carried by "
        f"`shear_stiffness_uncorrected`: missing={sorted(LAMINATE_RS_METHODS - methods)} "
        f"unexpected={sorted(methods - LAMINATE_RS_METHODS)}"
    )

    from aeroelast.core.laminate import Laminate, Ply  # noqa: PLC0415 - import cost is the point

    py_laminate = {n for n in dir(Laminate) if not n.startswith("_")}
    assert py_laminate == LAMINATE_PY_PUBLIC, (
        "the Python Laminate surface changed: the Python wrapper is a DIFFERENT set from the "
        f"Rust type: missing={sorted(LAMINATE_PY_PUBLIC - py_laminate)} "
        f"unexpected={sorted(py_laminate - LAMINATE_PY_PUBLIC)}"
    )
    py_ply = {n for n in dir(Ply) if not n.startswith("_")}
    assert py_ply == PLY_PY_PUBLIC, (
        f"the Python Ply surface changed: missing={sorted(PLY_PY_PUBLIC - py_ply)} "
        f"unexpected={sorted(py_ply - PLY_PY_PUBLIC)}"
    )
