# Certified zero and equality boundary

semialg distinguishes two different questions that must not be conflated.

- **Identity zero** asks whether an expression is the zero expression. `certified_zero` uses the package-wide `exprtest.zerotest` certification boundary for this purpose.
- **Pointwise zero** asks whether an expression is zero at every admissible parameter value relevant to the current decision. `certified_pointwise_zero` preserves `None` when a nonzero symbolic expression can still vanish on a parameter stratum.

This distinction matters for expressions such as `a`: it is not the zero polynomial, but it vanishes when `a = 0`. A negative identity result must therefore never by itself justify division by `a`, exclusion of the `a = 0` stratum, a regularity claim, or a scalar `a != 0` truth decision.

The certified-boundary audit removed all production decisions of the form `simplify(...) != 0`. The remaining `simplify(...) == 0` occurrences were reviewed as positive-only conservative recognizers: failure to recognize equality causes a fallback or missed simplification rather than certification of a false mathematical conclusion. They are tracked as a minor cleanup item in the release notes.
