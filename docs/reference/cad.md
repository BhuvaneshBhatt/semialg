# CAD and structured geometry reference

## `cad(...)`

Constructs or queries a cylindrical algebraic decomposition for supported semialgebraic input. Depending on options/output, CAD can support decision procedures, topology, structured cells, bounds, and integration adapters.

## `generic_cad(...)`

Builds parameter-aware/generic CAD information together with exceptional-set handling where supported. Specialization preserves the defining polynomial identity of algebraic root sections.

## Structured cell APIs

- `extract_structured_cad_cells`
- `structured_cad_cells_to_vertical_bounds_2d`
- `extract_vertical_bounds_from_cad_2d`
- `cylindrical_solution_from_structured`

These expose typed cylindrical structure rather than requiring callers to parse arbitrary Boolean formulas.

## Connectivity

- `build_cad_adjacency_graph`
- `extract_cad_connectivity`

These derive adjacency/components from CAD cell information.

## Bounds and certificates

Public CAD types include `CADBound`, `CADCellBoundsCertificate`, `StructuredCADLevel`, `StructuredCADCell`, `StructuredCADCellDecomposition`, `CADResult`, and related result/certificate objects.

Algebraic section bounds can be represented by certified root functions rather than approximate decimal endpoints.

## Exactness

Certified lifting, root ordering, sector sampling, and topology decisions use exact algebraic comparisons. Numerical approximations may be used for presentation/exploration only where explicitly allowed.

See [CAD and QE](../cad_qe.md) and [Performance](../guides/performance.md).
