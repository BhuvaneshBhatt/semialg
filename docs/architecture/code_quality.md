# Code Quality Notes

## Public API

The root package exports concise names for the main workflows:

- `cad` and `cad_text` for cylindrical algebraic decomposition
- `generic_cad` and `generic_cad_text` for parameter-generic decomposition
- `component_instances` for connected-component samples
- `find_instance` and `find_instance_text` for satisfying assignments
- `reduce_text` and `resolve_text` for text-based symbolic queries.

Long-form decomposition names are intentionally not exported.

## Algorithmic notes

Reduced CAD paths use explicit side-condition reports and certification. If a reduced projection path cannot be certified, the complete Collins decomposition is used. Algebraic samples are stored as explicit rational or algebraic objects so exact and approximate display are separate concerns.

## Refactoring constraints

Maintainability changes should not add dispatch or allocation overhead to CAD,
root isolation, exact sign determination, algebraic comparison, projection, or
RUR inner loops. Prefer direct module-level helpers for repeated normalization,
relation parsing, interval decomposition, canonical keys, and orchestration.

When a refactor affects a performance-sensitive public path, compare warmed
before/after timings and fresh-process import/API timings where relevant. Small
run-to-run differences in CAD-heavy tests should be treated as noise unless they
are reproducible across repeated measurements.

