# Implementation Notes

This document summarizes internal design choices for `semialg`.

## Core design

`semialg` exposes high-level APIs for real semialgebraic reasoning:

- `cad` for cylindrical algebraic decomposition
- `generic_cad` for generic parameter decompositions
- `component_instances` for connected-component samples
- `find_instance` for satisfying assignments
- `reduce_text` and `resolve_text` for text-based quantifier elimination and decision queries.

The implementation uses a conservative Collins-style CAD as the correctness baseline, with reduced McCallum, Lazard, and TTICAD paths accepted only when side conditions and certification checks succeed.

## Algebraic samples

CAD cells carry explicit rational or algebraic sample objects rather than untyped symbolic expressions. This keeps sample comparison, root isolation, sign evaluation, and exact/approximate display separated from formula simplification.

## Reduced CAD certification

Reduced projection backends may repair side-condition failures by adding local delineating factors. A reduced decomposition is accepted only after certification. Otherwise, the system uses the complete Collins decomposition.

## Public objects

`CADFunction`, `GenericCADFunction`, `CellSet`, and `InstanceResult` are intended to be stable user-facing objects. They preserve exact sample information and expose approximate forms for display.
