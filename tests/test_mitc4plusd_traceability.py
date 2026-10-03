"""Traceability gate for the MITC4+/D element documentation (SDD task 12.1).

Three scenarios, all of them reading files only -- no Rust build, no PETSc:

1. Every ingredient row of ``docs/formulations/shell-elements.md`` section 2.2 names
   a source, and every code symbol the row names EXISTS in
   ``crates/aeroelast-core/src/elements/mitc4.rs``. That is what makes this
   traceability instead of a second copy of the table: a renamed or deleted function
   fails here instead of leaving a stale claim in the documentation.
2. No forbidden, non-paper ingredient is present in the element's production code
   (the part before ``#[cfg(test)]``) or claimed by the table. The forbidden set is
   the one the spec enumerates: a drilling penalty, the ERC/``beta_w`` treatment,
   selective reduced integration of the in-plane shear, the rotation bubble, the
   classical ``5/6`` shear correction and hourglass scaffolding.
3. Every citation in the table resolves: author-year against the canonical
   bibliography ``docs/validation/references.yaml``, and the quoted equations against
   the paper extracts.

Each scenario carries a non-vacuity control, because a gate that scans nothing
passes for the wrong reason. The strongest one is in the forbidden-ingredient test:
the same detectors are run against ``mitc3.rs``, where ``k_drill`` legitimately
exists, and the test asserts the detector FIRES there.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest
import yaml

REPO = Path(__file__).resolve().parents[1]
DOC = REPO / "docs" / "formulations" / "shell-elements.md"
MODULE = REPO / "crates" / "aeroelast-core" / "src" / "elements" / "mitc4.rs"
MITC3 = REPO / "crates" / "aeroelast-core" / "src" / "elements" / "mitc3.rs"
# The canonical bibliography is the validation store. `docs/validation/references.yaml`
# is the rector; `docs/references.md` was a generated view of it and is gone.
BIBLIOGRAPHY = REPO / "docs" / "validation" / "references.yaml"
EXTRACTS = (
    REPO / "docs" / "formulations" / "mitc4plus-2017-extract.md",
    REPO / "docs" / "formulations" / "mitc4plusd-2025-extract.md",
)

TABLE_START = "### 2.2 Code → equation table"
TABLE_END = "### 2.3 "

#: The spec's forbidden, non-inherited ingredient set, each with the symbols that
#: would betray it in the code. Keys are the human names used in failure messages.
FORBIDDEN = {
    "the uncited drilling penalty": ("k_drill", "drilling_scale"),
    "the ERC drilling treatment and beta_w warping penalty": (
        "compute_ke_local_erc", "beta_w", "BETA_W", "erc_",
    ),
    "selective reduced integration of the in-plane shear": ("cm_normal",),
    "the 2-DOF rotation bubble": (
        "GpBubble", "b_kappa_bubble", "b_gamma_mitc4_plus", "bubble_function",
    ),
    "hourglass scaffolding": ("hg_stiffness_factor", "compute_hourglass", "h_orth"),
}
#: The classical shear correction factor is a number, not a symbol. The production
#: path must not apply it; the TESTS may compute it as the reference the element is
#: required to miss, which is why this is checked only outside `#[cfg(test)]`.
SHEAR_CORRECTION_LITERALS = ("5.0 / 6.0", "5.0/6.0", "5.0 / 6", "0.83333333")

_SYMBOL = re.compile(r"`([A-Za-z_][A-Za-z0-9_]*)`")
_CITATION = re.compile(r"([A-Z][A-Za-z]+)[^()]{0,80}?\((\d{4})\)")
_EQUATION = re.compile(r"Eqs?\.\s*\(?(\d+[a-z]?)")


def _table_rows() -> list[tuple[str, str, str]]:
    """The data rows of the section 2.2 ingredient table."""
    text = DOC.read_text()
    assert TABLE_START in text, f"{DOC.name} no longer has the {TABLE_START!r} heading"
    block = text.split(TABLE_START, 1)[1].split(TABLE_END, 1)[0]
    rows: list[tuple[str, str, str]] = []
    for line in block.splitlines():
        line = line.strip()
        if not line.startswith("|"):
            continue
        cells = [c.strip() for c in re.split(r"(?<!\\)\|", line)[1:-1]]
        if len(cells) != 3:
            continue
        if set(cells[0]) <= set("-: ") or cells[0].lower().startswith("code"):
            continue  # separator or header
        rows.append((cells[0], cells[1], cells[2]))
    return rows


def _production_source() -> str:
    """`mitc4.rs` up to the test module: the code that actually ships."""
    source = MODULE.read_text()
    marker = "#[cfg(test)]"
    assert marker in source, f"{MODULE.name} has no {marker} boundary to split at"
    return source.split(marker, 1)[0]


def test_traceability_identity_table_names_source_for_every_ingredient():
    rows = _table_rows()
    assert len(rows) >= 6, f"the ingredient table shrank to {len(rows)} rows: fix the gate, not the count"

    module = MODULE.read_text()
    checked = 0
    for code, paper, _verified in rows:
        symbols = _SYMBOL.findall(code)
        assert symbols, f"row {code!r} names no code symbol, so nothing ties it to the implementation"
        assert paper, f"row {code!r} names no source"
        for symbol in symbols:
            assert symbol in module, (
                f"{code!r} cites `{symbol}`, which `{MODULE.name}` does not define: "
                "the documentation has drifted away from the code"
            )
            checked += 1
    assert checked >= 6, f"only {checked} symbols were checked; the gate would be close to vacuous"


def test_traceability_no_forbidden_ingredient_claimed_or_present():
    production = _production_source()

    # (a) the production path must not contain any forbidden ingredient
    for name, symbols in FORBIDDEN.items():
        hits = [s for s in symbols if s in production]
        assert not hits, f"{name} is present in the production path of {MODULE.name}: {hits}"

    # (b) the classical shear correction factor must not be applied either
    hits = [lit for lit in SHEAR_CORRECTION_LITERALS if lit in production]
    assert not hits, f"the 5/6 shear correction factor is applied in the production path: {hits}"

    # (c) the table must not CLAIM any of them either, so the docs cannot assert an
    #     ingredient the element does not have
    forbidden_flat = [s for symbols in FORBIDDEN.values() for s in symbols]
    for code, _paper, _verified in _table_rows():
        claimed = [s for s in forbidden_flat if s in code]
        assert not claimed, f"the ingredient table claims the forbidden ingredient {claimed}"

    # (d) NON-VACUITY CONTROL: the detectors must be able to fire at all. `mitc3.rs`
    #     legitimately still has its own drilling penalty, so the same check there
    #     must find it. Without this, a typo in a symbol name would make (a) pass
    #     forever while checking nothing.
    mitc3 = MITC3.read_text()
    fired = [s for s in FORBIDDEN["the uncited drilling penalty"] if s in mitc3]
    assert fired, (
        "the forbidden-ingredient detectors did not fire on mitc3.rs either, which means they "
        "match nothing anywhere and scenario (a) is vacuous"
    )


def test_traceability_equation_citations_resolve_to_extract_or_paper():
    entries = yaml.safe_load(BIBLIOGRAPHY.read_text(encoding="utf-8"))["references"]
    extracts = "".join(path.read_text() for path in EXTRACTS)

    rows = _table_rows()
    resolved = 0
    for code, paper, _verified in rows:
        citations = _CITATION.findall(paper)
        assert citations, f"row {code!r} carries no author-year citation"
        for surname, year in citations:
            resolves = [
                entry
                for entry in entries
                if str(entry.get("year")) == str(year)
                and any(surname in author for author in entry.get("authors") or [])
            ]
            # Resolved against the store's entries rather than by grepping text: a year
            # appearing anywhere in the file is not a citation that resolves.
            assert resolves, (
                f"{surname} ({year}) does not resolve to an entry of "
                f"{BIBLIOGRAPHY.relative_to(REPO)}"
            )
            resolved += 1
        for equation in _EQUATION.findall(paper):
            assert f"({equation})" in extracts or f"({equation}" in extracts, (
                f"row {code!r} cites Eq. ({equation}), which appears in neither extract"
            )

    assert resolved >= 6, f"only {resolved} citations were resolved; the gate would be close to vacuous"


if __name__ == "__main__":  # pragma: no cover - manual diagnosis
    raise SystemExit(pytest.main([__file__, "-v"]))
