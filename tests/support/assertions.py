"""The suite's canonical comparison forms, and the shape the validation store reads.

A validation test is an ordinary pytest test: it runs, it asserts, and it fails when the physics is
wrong. What this module fixes is the *form* of one comparison, because the store has to read what a
test compares and against what.

Both forms name their reference. That name is not decoration: it is the only thing that distinguishes
our result against an independent reference from two of our own numbers, and no structural rule can
tell those apart -- which is why the store used to be given a per-test list of what to ignore instead
of reading the tests. Everything else a test asserts -- a symmetry, an invariant, a dominance
relation, a load ruler, a setup constant -- stays in the test exactly as it was and is simply not a
comparison.

Two forms, because the suite makes two kinds of claim, and the difference is real -- but both report
in one shape, so the store carries one residual pattern and a failure report needs one parse:

- :func:`assert_relative_error` -- our value against a reference *value*, compared by relative error.
- :func:`assert_residual_below` -- a residual the test already computed and printed against the bound
  the physics justifies. A residual is a relative difference already, so an absolute bound against a
  distance from a reference is what is being asserted, and the reference is what the residual was
  taken against.

Both print the residual they assert, in the form the group's residual pattern reads, because the
measured margin is transcribed from the test's own output and ``regression`` re-reads it there.
"""

from __future__ import annotations

# The four things a comparison may be against, as docs/validation-policy.md rule 6 defines them.
REFERENCE_KINDS = ("paper", "code", "analytical", "self")


def _check_kind(kind: str) -> str:
    if kind not in REFERENCE_KINDS:
        raise ValueError(
            f"reference kind {kind!r} is not one of {REFERENCE_KINDS}; rule 6 names what a "
            "reference may be"
        )
    return kind


def assert_relative_error(
    value: float,
    reference: float,
    *,
    tol: float,
    reference_name: str,
    kind: str,
    what: str,
) -> None:
    """Assert ``value`` is within ``tol`` of an independent ``reference``, by relative error.

    ``reference_name`` says what the reference *is* -- a paper's table, a code and its version, a
    closed form -- and reaches the store through the extractor, which is why the call is written in
    the test body and not inside a local helper.

    ``tol`` is authoritative: it is the bound the physics justifies, and this function adds no
    ceiling of its own.
    """
    if reference == 0:
        raise ValueError(
            "assert_relative_error needs a non-zero reference; a bound against zero is an absolute "
            "tolerance and belongs in assert_residual_below"
        )
    rel = abs(value - reference) / abs(reference)
    print(f"{what} vs {reference_name}: {value:.6g} vs {reference:.6g} (error: {rel:.4%}, bound: {tol:.4%})")
    assert rel < tol, (
        f"{what}: {rel:.4%} against {reference_name} {reference:.6e}, "
        f"above the {tol:.4%} bound (kind {_check_kind(kind)})"
    )


def assert_residual_below(
    residual: float,
    *,
    tol: float | None = None,
    atol: float | None = None,
    unit: str = "%",
    reference_name: str,
    kind: str,
    what: str,
) -> None:
    """Assert a residual the test already holds is below its bound, against a named reference.

    This is the suite's most common claim: a relative error against a published cell, a modal gap
    against another code, a distance from a closed form.

    ``tol`` is a fraction -- the residual is the difference divided by the reference -- and ``atol``
    is a quantity in ``unit``, for a residual that is not a fraction at all: an angle in degrees, a
    length, a force. Exactly one of the two is given, and the store takes the tolerance kind from
    which one it was, so an absolute bound is recorded as an absolute bound instead of being called
    a relative error.
    """
    if (tol is None) == (atol is None):
        raise ValueError("give exactly one of tol (a fraction) or atol (a quantity in unit)")
    limit = tol if tol is not None else atol
    if limit is None:  # the pair check above already refused this
        raise ValueError("give exactly one of tol (a fraction) or atol (a quantity in unit)")
    if tol is not None:
        shown = f"{residual:.4%}"
        bound = f"{tol:.4%}"
    else:
        shown = f"{residual:.6g} {unit}".strip()
        bound = f"{limit:.6g} {unit}".strip()
    print(f"{what} vs {reference_name}: {shown} (bound: {bound})")
    assert residual < limit, (
        f"{what}: {shown} against {reference_name}, above the {bound} bound "
        f"(kind {_check_kind(kind)})"
    )
