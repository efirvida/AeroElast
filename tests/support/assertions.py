"""The suite's one comparison form, and the shape the validation store reads.

A validation test is an ordinary pytest test: it runs, it asserts, and it fails when the physics is
wrong. What this module fixes is the *form* of one assertion, because the store has to read what a
test compares and against what.

A comparison is a call to :func:`assert_relative_error`, and it names its reference. That name is not
decoration: it is the only thing that distinguishes our result against an independent reference from
two of our own numbers, and no structural rule can tell those apart -- which is why the store used to
be given a per-test list of what to ignore instead of reading the tests. Everything else a test
asserts -- a symmetry, an invariant, a dominance relation, a load ruler, a setup constant -- stays in
the test exactly as it was and is simply not a comparison.

The function prints the residual it asserts, in the form the group's residual pattern reads, because
the measured margin is transcribed from the test's own output and ``regression`` re-reads it there.
"""

from __future__ import annotations


def assert_relative_error(
    value: float,
    reference: float,
    *,
    tol: float,
    reference_name: str,
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
            "tolerance and belongs in an explicit assertion of its own"
        )
    rel = abs(value - reference) / abs(reference)
    print(f"{what} vs {reference_name}: {value:.6g} vs {reference:.6g} (error: {rel:.4%})")
    assert rel < tol, (
        f"{what}: {value:.6e} against {reference_name} {reference:.6e} is {rel:.4%}, "
        f"above the {tol:.4%} bound"
    )
