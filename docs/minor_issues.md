# Lower-priority findings

| ID | Finding | Current contract |
| --- | --- | --- |
| M1 | Constant-matrix definiteness on an empty domain | The constant-inertia shortcut requires a true domain. Universal definiteness on an empty domain returns true. |
| M2 | Empty loci and `is_zero_dimensional` | The predicate means Krull dimension exactly zero. Its documentation excludes the empty locus, which has dimension -1. |
| M3 | SymPy matrix columns as maps | Formula-level image/preimage helpers normalize columns as coordinate tuples and reject non-column matrices. |

The previously observed `exprtools` signal-timeout warning during profiling remains a dependency diagnostic; it was not turned into a broad exception handler in semialg. See [Limitations](limitations.md) for geometric scope boundaries.
