# Limitations

`semialg` is conservative by design. The public API is broader than the fully general algorithms currently implemented underneath it.

## General limitations

- CAD/QE can be expensive, especially in higher dimensions (doubly exponential)
- Variable ordering matters
- Formula simplification is semantic but not canonical
- Some algorithms are exact for recognized cases and conservative elsewhere.

## Region Integration limitations

- Full arbitrary CAD-cell-to-bounds conversion is not complete
- Intrinsic-dimensional integration is implemented for selected cases, not arbitrary semialgebraic strata
- Higher-dimensional exact region integration relies on standard-shape recognition or explicit reducible forms
- Symbolic integration may fail even when the region decomposition is valid.

## Function Range limitations

`function_range` supports polynomial, rational, and common semialgebraic expression graphs. Arbitrary transcendental expressions are not generally supported because they fall outside first-order real closed field reasoning.

## Optimization limitations

The current optimizer supports useful exact cases but is not yet a complete global polynomial optimizer. Future improvements should add objective-value projection, KKT/critical-point methods, and certificate backends.
