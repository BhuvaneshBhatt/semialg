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

## General Boolean regions

For multidimensional formulas containing `Or` or `Not`, ambient-measure integration first evaluates the entire Boolean formula on an adapted CAD and selects the resulting disjoint full-dimensional cells. Each selected cell is converted to certified nested cylindrical bounds and integrated exactly once. This prevents overlap double-counting and supports general Boolean unions and bounded complements whenever the CAD cell bounds are integrable by the existing exact adapter.

Explicit `bounds=` are conjoined to the Boolean formula before CAD decomposition. This makes expressions such as the complement of a bounded hole integrable over a finite ambient box without incorrectly treating the complement as globally unbounded. `semialgebraic_measure` uses the same path because it delegates to `integrate_over_region` with integrand `1`.

## Parameter-dependent semialgebraic regions

This is distinct from `integrate_over_parametric_region`, which integrates over an explicitly parametrized curve, surface, or volume. `integrate_over_region(..., parameters=[...], return_stratified=True)` instead treats selected free symbols in the *region formula* as parameters and returns a certified `ParameterStratifiedResult`.

```python
a = sp.Symbol("a", real=True)

result = integrate_over_region(
    x,
    sp.And(x >= 0, x <= a),
    [x],
    parameters=[a],
    return_stratified=True,
)

result.select({a: 2})
# 2

result.select({a: -1})
# 0
```

The parameter-space CAD separates feasible fibers from empty fibers. The integration reducer then keeps the parameter symbols in certified symbolic bounds; empty fibers receive integral zero, so the result covers the full parameter space. The supported symbolic reductions remain valid over an entire parameter stratum, including parameter-dependent axis-aligned intervals/boxes. General nonlinear parameter-dependent CAD root-function bounds and intrinsic parameter-dependent measure are outside the supported integration fragment.

`semialgebraic_measure` accepts the same `parameters=` and `return_stratified=True` options because it delegates to region integration with integrand `1`.

## Supported standard shapes and forms

The current exact layer supports many common cases, including:

- intervals;
- axis-aligned boxes;
- the standard 2D unit simplex;
- origin-centered disks and annuli;
- axis-aligned ellipses;
- selected vertical-slice regions;
- arbitrary-dimensional full-dimensional regions whose selected CAD cells expose certified triangular bounds;
- automatic coordinate permutation when another CAD lifting order gives simpler certified iterated bounds;
- selected graph curves and circles for intrinsic one-dimensional measure.


## Symbol identity and string variable names

Public region APIs accept either SymPy symbols or string variable names. A string such as `"x"` is resolved against the symbols already present in the integrand, region formula, and bounds before any new symbol is created. This matters because `Symbol("x")` and `Symbol("x", real=True)` are distinct SymPy objects even though they print the same way. If two incompatible same-name symbols are genuinely present, the API raises `ValueError` rather than guessing which one the string denotes.

## Limitations

Arbitrary-dimensional CAD cells can be converted to typed nested cylindrical bounds. Variable-dependent algebraic sections are represented by certified `AlgebraicRootFunction` objects, and full-dimensional cells have direct iterated-integral adapters. Ambient integration performs a bounded variable-order search instead of assuming the caller's coordinate order: inexpensive explicit cylindrical permutations are considered first, followed by a small CAD-order set including the original order and Brown-style suggestions. Candidate decompositions are ranked by cell count and symbolic boundary complexity. Lower-dimensional cells use a separate intrinsic adapter based on the induced metric of verified triangular graph cells; singular or non-graph strata that cannot be certified still fail conservatively.

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


## Parameter-dependent algebraic endpoints

For one integration variable, parameter-dependent polynomial boundaries no longer have to reduce to explicit affine bounds. If the symbolic box/interval reducer declines, `integrate_over_region(..., parameters=..., return_stratified=True)` can build a complete CAD in `(parameters..., x)`. Parameter-level cells are induced by the full projection tower, so the real roots defining each fiber sector are delineable and retain a fixed order on the stratum. The exact antiderivative is then evaluated at the reconstructed algebraic root functions.

For example, `x**2 <= a` yields zero measure on `a <= 0` and an algebraic-root endpoint difference on `a > 0`, specializing exactly to `2*sqrt(a)` at concrete nonnegative algebraic/rational parameter values. Multidimensional nonlinear parametric cell bounds are outside the supported integration fragment.
