# Certified symbolic decisions

semialg distinguishes identity questions from pointwise parameter questions.

- **Identity zero** asks whether an expression is identically zero. `certified_zero` uses the package-wide `exprtest.zerotest` certification boundary.
- **Pointwise zero** asks whether an expression is zero for the admissible parameter values relevant to a decision. `certified_pointwise_zero` preserves `None` when a nonzero symbolic expression can still vanish on a parameter stratum.
- **Sign** uses `certified_sign`; parameter-dependent expressions remain unresolved unless assumptions establish a pointwise sign.

For example, `a` is not the zero polynomial but vanishes on the stratum `a = 0`. A negative identity result therefore cannot justify division by `a`, exclusion of that stratum, a regularity claim, or a pointwise `a != 0` conclusion.

Exact representation properties such as `Poly.is_zero` are safe to query directly. General SymPy expressions use the shared certification helpers whenever the result controls a mathematical branch. Conservative syntactic recognizers may decline a simplification, but they must not turn failure to recognize an identity into a mathematical negation.
