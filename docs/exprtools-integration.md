# Exprtools integration (stage 6)

The shared zero boundary imports `exprtools.zerotest` directly. Identity decisions still distinguish pointwise nonzero proofs, and UNKNOWN remains unresolved.

The new package is developed in the `exprtools/` subproject of
[BhuvaneshBhatt/exprtest](https://github.com/BhuvaneshBhatt/exprtest/pull/1).
The legacy exprtest distribution is a compatibility wrapper sharing the same
implementations and proof caches; new runtime imports use exprtools.

## Validation

Against baseline `4256346da9f9bf40863341edf8060ebc812f9784`, the selected downstream suite
passed 41 tests before import changes and 43
after changes and integration regressions. Tests were run on Python 3.12.14 /
SymPy 1.14.0 / python-flint 0.9.0. These are focused migration suites, not a claim
that this project's entire test suite or every optional integration was run.

Selected files:

- `test_certified_zero_boundaries.py`
- `test_certified_decision_portfolio.py`
- `test_symbolic_assumption_simplification.py`
- `test_exact_arithmetic_and_symbol_resolution.py`
- `test_exprtools_integration.py`

A clean wheel environment containing exprtools, asymptotic, semialg, funcprops,
odeanalysis and probstats was tested outside all checkouts with exprtest absent
and its import actively blocked. It checks package provenance, dependency
metadata, certified nested radicals, independent refinement, domain holes,
parameter uncertainty, recurrence evidence, and exact posterior normalization.
A second smoke run installs the compatibility wrapper and checks shared classes,
functions and proof caches. Build artifacts and raw test results accompany the
migration report in the originating task.

## Release and development

No distributions were published. Before releasing this change, publish exprtools
first, then semialg and funcprops, then their consumers; choose new downstream
release versions at release time. The dependency name is `exprtools>=0.1.0`.
CI installs the immutable exprtools source commit
`b6c38bd2b5fc27d8f88e30026d52868050308e18` pending its first release.
The old exprtest wheel is optional for consumers that still use the legacy API.
