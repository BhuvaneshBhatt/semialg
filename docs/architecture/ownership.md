# Architecture and API ownership

`semialg` is organized by mathematical ownership instead of by the history of individual algorithms.  New code should extend the narrowest owning subsystem and should not create a second implementation of an existing mathematical decision.

| Subsystem | Ownership |
|---|---|
| `decision`, `qe`, `cad_algorithms` | first-order decisions, quantifier elimination, cylindrical decomposition |
| `algebraic` | exact algebraic numbers, signs, comparisons, ideals and certificates |
| `regions`, `topology` | region operations, incidence, connectedness, triangulation and homology |
| `optimization` and `_optimization_*` | extrema, ranges, KKT geometry and certificates |
| `decomposition.parametric`, `parameters` | parametric CAD, exceptional-boundary provenance, certified strata, representative fibers, and parameter feasibility |
| `region_integrate` and `_region_integrate_*` | ambient and intrinsic integration |
| `function_graph` | exact graph encodings used by image/range operations |

The root namespace is reserved for primary user operations and their result/certificate types.  Internal helpers belong to their owning subsystem.  In particular, topology implementation modules live under `semialg.topology`; new topology code should not add another top-level `*_topology.py` module.

## Algorithm-selection rule

Cheap exact structure is attempted before generic CAD/QE.  A specialized path may return a mathematical result only when its side conditions are certified.  Failure to certify is not evidence that the side condition is false: the operation falls through to another exact backend or reports unsupported computation.

This rule is especially important for bounds/signs, parameter representatives, graph algebraization, and integration charts.  It prevents easy polynomial or polyhedral problems from becoming unnecessarily expensive CAD problems without weakening exactness.

### Algebra implementation ownership

Large exact-algebra modules keep their public algorithm entry points stable while
private support code is divided by semantic responsibility. Immutable result and
replay-certificate records live in `_modular_types`, `_gtz_primary_types`, and
`_root_certificate_types`. Equality-only zero-dimensional quotient algebra, RUR
construction/solving, exact border bases, and their low-level linear-algebra kernels
are canonically owned by `algroots>=0.8.0`; Semialg keeps compatibility facades and
adds the real-sign, Thom, Boolean-formula, inequality-filtering, witness, and
decision/QE layers that are genuinely semialgebraic. New equality-only quotient
algorithms must be added to algroots rather than reimplemented under `semialg.algebraic`.

The remaining public algebraic algorithm modules retain hot loops and backend calls
directly for modular reconstruction, GTZ recursion, and root isolation/refinement.

### Exact algebra implementation boundaries

The large exact-algebra modules use semantic ownership instead of size-based splitting. Function-field representation/arithmetic remains in `algebraic_function_fields`, univariate factorization lives in `_function_field_factorization`, and primitive-element realization plus replay lives in `_function_field_compression`. Modular Gröbner realization remains in `algebraic.modular`, resultant/subresultant realization lives in `_modular_resultants`, and certificate replay lives in `_modular_certification`. GTZ recursion remains in `gtz_primary` while replay/schema validation lives in `_gtz_primary_certification`. Root isolation/refinement remains in `roots`; root-interval certificate construction/replay lives in `_root_certification`.

These boundaries use direct function re-exports instead of forwarding wrappers, so the split does not add calls inside computational loops.

## Algroots namespace contract

Semialg 1.3.2 requires algroots 0.8.0 or newer. Delegation uses `algroots.quotient`, `algroots.rational_univariate`, `algroots.border_basis`, and `algroots.errors`; advanced algroots names are no longer imported from its root. The Semialg-facing adapters retain their branch filtering, error translation and point enrichment.
