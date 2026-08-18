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
