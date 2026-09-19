# Code Quality Notes

## Public API

The root package exports concise names for the main workflows:

- `cad` and `cad_text` for cylindrical algebraic decomposition
- `generic_cad` and `generic_cad_text` for parameter-generic decomposition
- `component_instances` for connected-component samples
- `find_instance` and `find_instance_text` for satisfying assignments
- `reduce_text` and `resolve_text` for text-based symbolic queries.

Long-form decomposition names are not exported.

## Algorithmic notes

Reduced CAD paths use explicit side-condition reports and certification. If a reduced projection path cannot be certified, the complete Collins decomposition is used. Algebraic samples are stored as explicit rational or algebraic objects so exact and approximate display are separate concerns.

## Performance constraints

Maintainbility changes should not add dispatch or allocation overhead to CAD,
root isolation, exact sign determination, algebraic comparison, projection, or
RUR inner loops. Prefer direct module-level helpers for repeated normalization,
relation parsing, interval decomposition, canonical keys, and orchestration.

When a structural change affects a performance-sensitive public path, compare warmed
before/after timings and fresh-process import/API timings where relevant. Small
run-to-run differences in CAD-heavy tests should be treated as noise unless they
are reproducible across repeated measurements.


## Focused implementation modules

Large public modules should not accumulate unrelated implementation layers.
The region modules separate exact region algebra from numerical geometry:

```text
semialg.cad_region
    exact CAD-region algebra/signatures/combine/extend
    re-exports numerical geometry from:
        semialg.cad_algorithms.meshing
        semialg.cad_algorithms.numerical_boundaries

semialg.symbolic_regions
    SemialgebraicRegion core value and topology convenience methods
    re-exports coercion and predicates from:
        semialg.region_coercion
        semialg.region_predicates

semialg.function_analysis
    public function-property facade and shared result/context types
    lazily delegates focused implementation to:
        semialg._function_analysis_convexity
        semialg._function_analysis_monotonicity
        semialg._function_analysis_partitions
        semialg._function_analysis_properties

semialg.optimization
    public optimization/range facade and orchestration
    delegates certification and structural specializations to:
        semialg._optimization_certification
        semialg._optimization_specializations
        semialg._optimization_parametric
        semialg._optimization_range
        semialg._range_special_cases
        semialg._function_graph_image
        semialg._critical_loci

semialg.decision.api
    public decision facade
    delegates predicate and solve implementations to:
        semialg.decision._predicates
        semialg.decision._solve_api

semialg.region_integrate
    public integration facade
    delegates geometric and CAD reduction to:
        semialg._region_integrate_geometry
        semialg._region_integrate_cad
        semialg._region_integrate_intrinsic
        semialg._region_integrate_parametric
```

Public modules re-export the *same function/class objects* from focused modules,
so triangulation, root evaluation, membership, and relation predicates do not
acquire an extra call layer.

Likewise, the top-level package namespace is declared in
`semialg._public_api.PUBLIC_EXPORTS`. `semialg.__all__` and lazy `__getattr__`
are derived from that one registry so documentation and star-import surfaces do
not maintain separate lists that can drift.

## Exception-boundary policy

Internal fallbacks catch explicit expected failure classes such as
`EXACT_OPERATION_ERRORS`, narrower polynomial error tuples, or package-specific
strategy failures. Production code does not catch `Exception` or
`BaseException`: unexpected programming errors must surface rather than being
converted into ordinary unsupported/unknown results. Optional dependencies catch
their import errors specifically.

The source-quality verifier and Ruff enforce this rule across the full package.

## Source hygiene

Comments should explain invariants, mathematical assumptions, certification
conditions, performance-sensitive choices, or non-obvious failure behavior.
Comments that only restate the following line of code should be omitted.
Substantial algorithms should have docstrings describing their mathematical
contract and the conditions under which they decline to produce a result.

## Type annotations

Annotations define stable interfaces, dataclass fields, callback contracts, and
non-obvious container shapes. Local variables are normally left to inference
unless an empty container or union-valued state would otherwise obscure the
intended invariant. The codebase does not use annotations as decoration: a type
that cannot be stated more precisely than `object` is kept broad only at public
sympification boundaries or genuinely heterogeneous result seams.
