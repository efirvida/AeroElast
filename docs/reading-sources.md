# Reading sources

The single protocol for turning a held source into a quotable equation, table cell,
figure value or number. Every document that quotes a source points here. Do not restate
these rules elsewhere: restating them is how they drift apart, which is what this file
replaced.

Where sources live: `.sources/papers/` (gitignored, recoverable from git history — see
`docs/validation/references.yaml`). The bibliography that keys them is
`docs/validation/references.yaml`.

## Why this exists

The recovered PDFs include scans with **no text layer over their mathematics**.
`pdftotext` silently drops, reorders or interleaves equations, matrices and table cells.
A number obtained that way is not a reference; it is a transcription of an extraction
failure. Two of those failures are recorded in
`docs/formulations/mitc4plus-2017-extract.md`: `pdftotext` put Eq. (16) on a PDF page
where it is not, and the first revision of the validation matrix was assembled with
`pdftotext -layout` and had to be redone.

## Render the page, then read it with vision

```bash
pdfinfo .sources/papers/<file>.pdf | grep Pages      # how long is it
pdftoppm -png -r 200 -f <first> -l <last> .sources/papers/<file>.pdf /tmp/<slug>
```

Then read the rendered page image with vision.

- **200 dpi** reads body text and simple tables. **300 dpi** for subscripts, small type,
  dense tables and anything you must transcribe symbol by symbol.
- **Render only the pages that can hold the item.** `pdftotext -f p -l p` may be used to
  *locate* a page, and never to read a value: confirm visually before quoting.
- **PDF page number is not printed page number.** Record the printed one, and say so
  when the offset is not obvious.

## Never

- `pdftotext` for an equation, a matrix, a table cell, or any number you will quote.
- A value reconstructed from memory, from a similar paper, or from an earlier draft.
- A figure value read off a low-resolution render without stating how it was read.
- A citation whose page or cell you have not seen.

## Equations

Transcribe in the paper's own notation, symbol by symbol, with the equation numbers.
When the paper's notation differs from the repository's, record **both** and the mapping
between them; a silent change of notation is how an implementation diverges from its
source. Then validate the transcription against something independent: a limiting case,
a dimensional check, or a quantity the paper also tabulates. An equation that cannot be
validated is reported as unvalidated, never quoted as settled.

## Tables

A cell is only a reference together with its coordinates. Record, every time:

- the table number and its caption (including a normalization note such as "normalized
  against the Kirchhoff solution" — the same cell means different things normalized
  differently);
- the **column** header and the **row** header, verbatim;
- the mesh and parameters the cell belongs to, e.g. `Table 12, MITC4+ N=16,
  t/L = 0.0002667`;
- the element the column belongs to: these papers tabulate MITC4, MITC4+ and other
  codes side by side, and reading the wrong column is indistinguishable from a wrong
  result.

Quote the cell verbatim. Never round it, and never complete a truncated one.

## Figures and graphs

- Record axis labels, ranges and units before reading anything off the curve.
- Digitize against the tick grid, and state the reading uncertainty (for example
  "about ±0.5% of full scale"). A digitized value is weaker evidence than a tabulated
  one; when the paper also tabulates it, use the table.
- A curve or a mode shape carries no single number. Record the quantity *and* the point:
  `first flapwise mode, 0.526 Hz, Fig. 7(a)`.

## Record the provenance

Every quoted value carries three things, and a citation missing any of them is not
evidence: the **source key** (as in `docs/validation/references.yaml`), the **printed
page**, and the **locator** (`Eq. (9)-(11)`, `Table 12`, `Fig. 5`). Written out:

```text
Ko, Lee & Bathe 2017, p. 406, Eqs. (9)-(11)
Ko, Lee, Lee & Bathe 2017, Table 12, MITC4+ N=16
```

In the validation store this is the comparison's `reference.label`; in a formulation note
it sits next to the equation it justifies.

## Token economy

Reading a page with vision costs tokens, so:

- render the page you need, not the document;
- ask one targeted question per page ("what is the N=16, MITC4+ cell of Table 12?")
  instead of asking the reader to describe the page;
- **read the extract documents first.** `docs/formulations/*-extract.md` already record
  the equations and cells that were read with vision, with their page cites. Re-extract
  only what is missing, and add what you read so the next reader does not pay for it
  again;
- when a value is already in a store (`docs/validation/references.yaml`,
  `docs/validation/rows/*.yaml`), quote the store and cite its provenance rather than
  re-reading the paper.
