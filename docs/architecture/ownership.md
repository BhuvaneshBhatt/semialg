# Architecture and API ownership

`semialg` is organized by mathematical ownership rather than by the history of individual algorithms.  New code should extend the narrowest owning subsystem and should not create a second implementation of an existing mathematical decision.

| Subsystem | Ownership |
|---|---|
| `decision`, `qe`, `cad_algorithms` | first-order decisions, quantifier elimination, cylindrical decomposition |
| `algebraic` | exact algebraic numbers, signs, comparisons, ideals and certificates |
| `regions`, `topology` | region operations, incidence, connectedness, triangulation and homology |
| `optimization` and `_optimization_*` | extrema, ranges, KKT geometry and certificates |
| `parameter_stratification`, `parameters` | parameter exceptional sets, certified strata and representative fibers |
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
`_root_certificate_types`; border-basis monomial/order-ideal normalization lives
in `_border_basis_normalization`, while `_border_basis_types` owns diagnostics and
errors. The public algorithm modules retain hot loops and backend calls directly,
so the separation does not add wrapper calls to modular reconstruction, GTZ
recursion, root isolation/refinement, or Macaulay realization.

### Exact algebra implementation boundaries

The large exact-algebra modules use semantic ownership rather than size-based splitting. Function-field representation/arithmetic remains in `algebraic_function_fields`, univariate factorization lives in `_function_field_factorization`, and primitive-element realization plus replay lives in `_function_field_compression`. Modular Gröbner realization remains in `algebraic.modular`, resultant/subresultant realization lives in `_modular_resultants`, and certificate replay lives in `_modular_certification`. GTZ recursion remains in `gtz_primary` while replay/schema validation lives in `_gtz_primary_certification`. Root isolation/refinement remains in `roots`; root-interval certificate construction/replay lives in `_root_certification`.

These boundaries use direct function re-exports rather than forwarding wrappers, so the split does not add calls inside computational loops.
