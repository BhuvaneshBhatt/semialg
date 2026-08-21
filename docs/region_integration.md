# Region integration and measure

`semialg` provides a layered, exact-first region integration framework.

## Public API

- `reduce_region_integral`
- `integrate_over_region`
- `semialgebraic_measure`

## Architecture

The integration engine is layered:

1. recognize standard shapes where exact moment formulas are available;
2. reduce supported semialgebraic regions to iterated-integral pieces;
3. evaluate symbolically, numerically, or automatically depending on `method`;
4. support intrinsic-dimensional integration for selected lower-dimensional regions.

## Reducing to iterated integrals

```python
import sympy as sp
from semialg import reduce_region_integral

x, y = sp.symbols("x y", real=True)

red = reduce_region_integral(1, x**2 + y**2 <= 1, [x, y])
```

A reduced integral consists of one or more pieces, each with an integrand, limits, and a sign. For the unit disk, the result is equivalent to:

```python
sp.Integral(1, (y, -sp.sqrt(1 - x**2), sp.sqrt(1 - x**2)), (x, -1, 1))
```

## Symbolic, numeric, and auto modes

```python
from semialg import integrate_over_region

integrate_over_region(1, x**2 + y**2 <= 1, [x, y], method="symbolic")
# pi

integrate_over_region(1, x**2 + y**2 <= 1, [x, y], method="numeric")
# approximate numeric value

integrate_over_region(1, x**2 + y**2 <= 1, [x, y], method="auto")
# symbolic if possible, numeric otherwise
```

`method="symbolic"` is the default. It requires exact symbolic evaluation and raises `NotImplementedError` if any reduced piece remains unevaluated.
 It also never substitutes numerical root finding for failed exact real-root isolation: if an exact boundary cannot be isolated, the exact reduction is declined rather than approximated.

## Measure

```python
from semialg import semialgebraic_measure

semialgebraic_measure(x**2 + y**2 <= 1, [x, y])
# pi

```


### Bound validation

`bounds=` is validated centrally. Every bound key must resolve to one of the declared
integration variables; duplicate variables, malformed pairs, and exactly reversed
endpoints raise `ValueError`. This prevents misspelled bounds from being silently
ignored and prevents reversed intervals from producing signed/negative measures.

`integrate_over_standard_region` follows the same symbol-identity rule as formula-based
integration: a string variable name is resolved against the actual SymPy symbol in the
integrand before a new real symbol is created. Boolean intersections of explicit
intervals and boxes compare exact algebraic endpoints without floating-point ordering.

## Intrinsic-dimensional measure

By default, measure and integration use ambient Lebesgue measure in the supplied variables.

```python
semialgebraic_measure(sp.Eq(x**2 + y**2, 1), [x, y])
# 0 under ambient 2D measure

semialgebraic_measure(sp.Eq(x**2 + y**2, 1), [x, y], measure_dimension="intrinsic")
# 2*pi in supported cases

integrate_over_region(x**2, sp.Eq(x**2 + y**2, 1), [x, y], measure_dimension=1)
# pi in supported cases
```

## Supported standard shapes and forms

The current exact layer supports many common cases, including:

- intervals;
- axis-aligned boxes;
- the standard 2D unit simplex;
- origin-centered disks and annuli;
- axis-aligned ellipses;
- selected vertical-slice regions;
- selected graph curves and circles for intrinsic one-dimensional measure.


## Symbol identity and string variable names

Public region APIs accept either SymPy symbols or string variable names. A string such as `"x"` is resolved against the symbols already present in the integrand, region formula, and bounds before any new symbol is created. This matters because `Symbol("x")` and `Symbol("x", real=True)` are distinct SymPy objects even though they print the same way. If two incompatible same-name symbols are genuinely present, the API raises `ValueError` rather than guessing which one the string denotes.

## Limitations

Arbitrary-dimensional CAD cells can be converted to typed nested cylindrical bounds. Variable-dependent algebraic sections are represented by certified `AlgebraicRootFunction` objects, and full-dimensional cells have direct iterated-integral adapters. Lower-dimensional cells use a separate intrinsic adapter based on the induced metric of verified triangular graph cells; singular or non-graph strata that cannot be certified still fail conservatively.

## Typed CAD bounds

`extract_cylindrical_solution(...)` preserves both expression bounds and typed bounds.
The typed form distinguishes explicit/infinite endpoints from delineable
`AlgebraicRootFunction` boundaries and retains open/closed sector information.
`verify_cad_cell_bounds(cell)` checks triangular variable dependence, section
root certificates, adjacent-root ordering, and sample containment.

For direct cell integration, use `full_dimensional_cell_integral(...)`.
Lower-dimensional cells are intentionally handled by the separate
`intrinsic_cell_integral(...)` adapter, which uses the induced metric
`sqrt(det(J.T*J))` for verified triangular graph cells.

## Regular/singular intrinsic stratification

`stratify_intrinsic_solution(...)` classifies cylindrical solution cells by certified regularity. `IntrinsicStratification.regular_strata` and `.singular_strata` keep the distinction explicit. Algebraic sections require a cell-wide `DelineabilityCertificate` with verified regularity; an algebraic root function without such a certificate is never silently treated as a regular manifold graph.

For a regular graph cell with mapping Jacobian `J`, intrinsic integration uses

```text
sqrt(det(J.T * J))
```

as the Hausdorff metric factor. Singular strata of the requested dimension cause verified intrinsic integration to decline; lower-dimensional singular strata remain inspectable and do not contribute to a higher-dimensional Hausdorff measure.

## Explicit and parametric region validation

Explicit standard-region objects validate geometric invariants at construction time.
`IntervalRegion` and `BoxRegion` reject bounds whose order is exactly known to be
reversed. Radius-based regions reject provably negative radii, and
`SphericalShellRegion` additionally requires the inner radius not to exceed the
outer radius when that ordering is exactly decidable. Constructors for compound
regions check ambient-coordinate dimensions before integration. Symbolic values
whose sign or order cannot be established exactly are not rejected merely because
they are undecidable.

`ParametricRegion` requires each declared parameter to have exactly one integration
limit, rejects undeclared or duplicate limit variables, and requires multiplicity
to be provably positive. Parametric integration resolves string ambient-variable
names against the symbols already present in the integrand and mapping, preserving
SymPy symbol identity and assumptions.

Boolean-region integration uses exact intersection semantics. In particular,
`RegionDifference(A, B)` integrates over `A \ B`, equivalently subtracting the
integral over `A ∩ B`; it does not assume that `B` is contained in `A`.
