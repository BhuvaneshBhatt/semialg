# semialg

`semialg` is an experimental Python package for exact symbolic reasoning over real polynomial and semialgebraic conditions. It provides CAD/QE-backed tools for satisfiability, implication, solving, region operations, symbolic optimization, range computation, semialgebraic integration, and assumption-based simplification.

The package is intended for research, education, and experimentation. It is not yet a production-grade replacement for mature computer algebra systems.

## Installation

```bash
python -m pip install semialg
```

## Quick start

```python
import sympy as sp
from semialg import is_satisfiable, implies, equivalent

x, y = sp.symbols("x y", real=True)

is_satisfiable(sp.And(x**2 + y**2 <= 1, x > 0, y > 0), [x, y])
# True

implies(x > 1, x**2 > 1, [x])
# True

equivalent(x**2 <= 1, sp.And(x >= -1, x <= 1), [x])
# True
```

## Features

### CAD and quantifier elimination

- Cylindrical algebraic decomposition for real polynomial constraints,
- Complete real quantifier-elimination workflows in supported cases,
- Generic decompositions and exceptional-set workflows.
- Real witness search and sample-point generation.

### Decision and solving

- `is_satisfiable`
- `is_tautology`
- `implies`
- `equivalent`
- `solve_semialgebraic`
- `sample_point`, `sample_points`
- `sign_at`, `sign_vector`

The decision helpers preserve their simple boolean API by default, but can return structured decision result objects with witnesses or counterexamples:

```python
res = implies(x >= 0, x > 0, [x], return_result=True)
res.valid          # False
res.counterexample # for example {x: 0}
```

Available result classes include `SatisfiabilityResult`, `TautologyResult`,
`ImplicationResult`, and `EquivalenceResult`. Witnesses are validated against
the original formula before they are returned.

Exact finite-system dispatch is available in this layer. `solve_semialgebraic(...,
method="auto")`, `is_satisfiable(..., return_result=True)`, and the public sampling helpers opportunistically use the RUR backend for supported finite zero-dimensional equality branches before falling back to interval/CAD-style methods. This means finite systems with irrational witnesses can be handled exactly:

```python
finite = sp.And(sp.Eq(x**2 + y**2, 1), sp.Eq(x, y), x > 0)
sol = solve_semialgebraic(finite, [x, y], method="rur", count=10)
sol.method   # "rational_univariate"
sol.samples  # ({x: sqrt(2)/2, y: sqrt(2)/2},)
```
Sampling distinguishes certified representative samples from explicit grid, random, and CAD-cell sampling workflows. The default remains exact and validated, while numerical random sampling is opt-in via `exact=False`:

```python
# Exact representative sampling: RUR first for finite systems, then rational
# witnesses, then requested CAD representatives.
sample_point(sp.Eq(x**2, 2), [x], strategy="representative")

# Deterministic rational grid sampling over a plotting/inspection window.
sample_points(
    (x >= 0) & (y >= 0) & (x + y <= 1),
    [x, y],
    count=5,
    strategy="grid",
    bounds=[(0, 1), (0, 1)],
    grid_resolution=5,
)

# Seeded numerical sampling, intended for exploratory/visual workflows.
sample_points(
    x**2 + y**2 < 1,
    [x, y],
    count=10,
    strategy="random",
    bounds=[(-1, 1), (-1, 1)],
    seed=1,
    exact=False,
)
```

Supported public sampling strategies are `"representative"`/`"auto"`,
`"rational"`, `"grid"`, `"random"`, and `"cad_cells"`. Every point returned
by the public helpers is checked against the input formula before it is exposed.


Sign evaluation strengthens `sign_at` and `sign_vector`. They use exact algebraic sign decisions for rational values, SymPy algebraic numbers such as `sqrt(2)`/`RootOf`, semialg `AlgebraicRoot` samples, and RUR-backed points.
RUR point signs are evaluated by substituting the RUR coordinate polynomials, reducing the expression modulo the defining polynomial, and deciding the resulting univariate algebraic sign. Numeric fallback is available only when `exact=False`:

```python
from semialg.algebraic import compute_rational_univariate_representation, solve_rur_points

rep = compute_rational_univariate_representation([x**2 + y**2 - 1, x - y], [x, y])
pt = next(p for p in solve_rur_points(rep) if sign_at(x, p, variables=[x, y]) > 0)

sign_vector([x - y, x + y, x**2 + y**2 - 1], pt, variables=[x, y])
# (0, 1, 0)
```

Decision contract coverage documents and tests the public decision-layer contracts.
The main user-facing guide is `docs/user_guide/decision_sampling_signs.md`;
the quality checklist is `docs/quality/decision_contracts.md`. The covered contracts are: boolean API compatibility, structured results, validated witnesses/counterexamples, `count=0` solving semantics, explicit sampling strategies, exact sign evaluation, and warning-safe Boolean normalization.


### Parameters and roots

- `classify_real_roots`
- `solvability_conditions`
- `root_count_conditions`
- `root_of` and related algebraic-bound infrastructure


### Exact zero-dimensional solving with RUR

`semialg` includes an exact rational-univariate-representation backend for finite polynomial systems over the rationals. The public entry point is
`solve_zero_dimensional_system(equations, inequalities=None, vars=[...],
backend="rur")`. It first checks that the equality ideal is zero-dimensional, then constructs a Rouillier-style quotient-algebra RUR using traces, a separating linear form, and a squarefree univariate parameter polynomial.
Optional inequalities and Boolean combinations of relational atoms are applied as exact filters on the algebraic candidate points.

```python
import sympy as sp
from semialg import solve_zero_dimensional_system

x, y = sp.symbols("x y", real=True)

result = solve_zero_dimensional_system(
    [sp.Eq(x**2 + y**2, 1), sp.Eq(x - y, 0)],
    inequalities=x > 0,
    vars=[x, y],
)

result.points
# ((sqrt(2)/2, sqrt(2)/2),)

result.representation.defining_polynomial
# Poly(_rur_t**2 - 2, _rur_t, domain='QQ')
```

The RUR result keeps useful exact metadata: the quotient-algebra dimension, the number of distinct geometric solutions, the separating linear form, and the coordinate rational functions in the univariate parameter. Higher-level solving paths use this backend opportunistically for finite equality branches before falling back to CAD/virtual-substitution style workflows for positive-dimensional semalgebraic sets.

### Region operations

- `region_union`, `region_intersection`, `region_difference`, `region_complement`
- `region_closure`, `region_interior`, `region_boundary`
- `region_dimension`, `region_components`
- `region_subset`, `region_equal`, `region_disjoint`
- `region_bounded`, `region_closed`, `region_compact`

### Optimization and range computation

- `semialgebraic_minimize`
- `semialgebraic_maximize`
- `function_range`

`function_range` treats range computation as a semialgebraic image problem. For an expression `f(x)` over a domain `C(x)`, it constructs a graph relation and eliminates the original variables to obtain a condition on a value symbol.

```python
import sympy as sp
from semialg import function_range

x, t = sp.symbols("x t", real=True)

function_range(x, sp.Or(x <= -1, x >= 1), [x], value_symbol=t)
# Abs(t) >= 1

function_range(sp.sqrt(1 - x**2), True, [x], value_symbol=t)
# (t >= 0) & (t <= 1)
```

### Region integration and moments

- `reduce_region_integral`
- `integrate_over_region`
- `semialgebraic_measure`
- `region_moment`
- `region_centroid`
- `region_covariance`

Examples:

```python
from semialg import semialgebraic_measure, integrate_over_region

semialgebraic_measure(x**2 + y**2 <= 1, [x, y])
# pi

integrate_over_region(x**2 + y**2, x**2 + y**2 <= 1, [x, y])
# pi/2

semialgebraic_measure(sp.Eq(x**2 + y**2, 1), [x, y], measure_dimension="intrinsic")
# 2*pi
```

### Symbolic simplification

- `simplify_boole`
- `simplify_piecewise`
- `simplify_system`
- `simplify_under_assumptions`
- `prove_positive`, `prove_nonnegative`, `prove_negative`, `prove_nonpositive`

```python
from semialg import simplify_boole, simplify_under_assumptions

simplify_boole((x > 1) & (x >= 0), [x])
# x > 1

simplify_boole((x**2 <= 1) & (x >= 0), [x])
# (x >= 0) & (x <= 1)

simplify_boole(2*x - 2 >= 0, [x])
# x >= 1

simplify_under_assumptions(sp.Abs(x), x >= 0, [x])
# x

simplify_under_assumptions(sp.sqrt((x - 1)**2), x >= 1, [x])
# x - 1

simplify_under_assumptions(sp.sqrt(x**2*y**2), (x >= 0) & (y <= 0), [x, y])
# -x*y

simplify_under_assumptions(sp.log(x**2), x > 0, [x])
# 2*log(x)

simplify_under_assumptions((x**2 - 1)/(x - 1), x > 1, [x])
# x + 1

conditional = simplify_under_assumptions((x**2 - 1)/(x - 1), True, [x], return_conditions=True)
conditional.expression, conditional.conditions
# (x + 1, (Ne(x - 1, 0),))
```

`BooleanSimplificationResult` is available with `simplify_boole(..., return_result=True)` for diagnostics, chosen variables, assumptions, and the simplified formula. `AssumptionSimplificationResult` is available from `simplify_under_assumptions(..., return_result=True)` or `return_conditions=True`; it records side conditions for domain-changing rewrites such as rational cancellation.

## Documentation

Start with:

- [`docs/index.md`](docs/index.md) for overview
- [`docs/api_overview.md`](docs/api_overview.md) for the public API map
- [`docs/function_range.md`](docs/function_range.md) for range computation
- [`docs/region_integration.md`](docs/region_integration.md) for integration and measure
- [`docs/symbolic_simplification.md`](docs/symbolic_simplification.md) for CAD-backed simplification
- [`docs/implementation_notes.md`](docs/implementation_notes.md) for design notes
- [`docs/future_directions.md`](docs/future_directions.md) for the roadmap
- [`notebooks/semialg_demo.ipynb`](notebooks/semialg_demo.ipynb) for an expanded guided notebook covering algebraic-geometry background, CAD/QE theory, RUR, region workflows, plotting, optimization, integration, topology, applications, exact-vs-numeric workflows, and performance limitations.

## Current limitations

`semialg` is conservative by design: it returns exact answers for supported cases and raises `NotImplementedError` rather than silently producing unreliable results.

Important current limitations:

- Full arbitrary CAD-cell-to-bounds conversion is not complete
- General intrinsic Hausdorff-measure integration over arbitrary semialgebraic strata is not complete
- Higher-dimensional CAD can become expensive (doubly exponential)
- Some region-integration paths currently rely on recognized shapes, vertical-slice decompositions, or specific parametrizations
- `function_range` supports polynomial, rational, and common semialgebraic expression graphs such as `Abs`, `sqrt`, `Min`, `Max`, and simple `Piecewise`, but not arbitrary transcendental expressions
- Formula simplifiers aim to remove contradictions, unreachable branches, and provably redundant conditions; they do not yet compute a unique minimal canonical semialgebraic formula in every case.

## Development

```bash
ruff format .
ruff check .
pytest
python scripts/clean_artifacts.py
```

### Subresultant PRS utilities

SymPy already exposes subresultant polynomial sequences through `sympy.subresultants` / `sympy.polys.polytools.subresultants`. `semialg` wraps that exact functionality in a small CAD/QE-oriented result object:

```python
import sympy as sp
from semialg import subresultant_prs

x = sp.symbols("x")
prs = subresultant_prs(x**3 - 2*x + 1, x**2 - 1, x, domain=sp.QQ)
prs.polynomials                 # exact Poly entries
prs.principal_coefficients      # leading/principal coefficients used by PRS decisions
prs.resultant                   # resultant of the two inputs
prs.source                      # normally "sympy.subresultants"
```

The wrapper gives internal algorithms a stable package-level interface for resultants, gcd-degree tests, projection diagnostics, and future subresultant-based Cauchy-index/root-counting code.

### Exact border-basis utilities

SymPy has mature Groebner-basis support, but it does not currently expose a
public border-basis algorithm comparable to its Gröbner APIs. `semialg`
therefore includes a small exact border-basis layer for zero-dimensional
rational ideals:

```python
import sympy as sp
from semialg import compute_border_basis, compute_border_basis_linear

x, y = sp.symbols("x y")
bb = compute_border_basis([x**2 - 1, y - x], [x, y])
bb_linear = compute_border_basis([x**2 - 1, y - x], [x, y], algorithm="linear")
# equivalent convenience helper:
bb_linear = compute_border_basis_linear([x**2 - 1, y - x], [x, y])

bb.order_monomials       # quotient basis/order ideal O
bb.border_monomials      # border dO
bb.as_exprs()            # border relations sigma - sum c_m*m
bb.normal_form(y)       # x
bb.coordinates(y)       # coordinate column in O
bb.multiplication_matrix(x)
bb.multiplication_matrix(x + y)
bb.commutation_certificate
bb.diagnostics
bb.has_commuting_multiplication_matrices()
```

The result object also exposes construction diagnostics useful for release
audits and backend debugging: quotient-basis rank, border-polynomial rank, exact commutator matrices, and failure messages. By default invalid order ideals or positive-dimensional inputs raise `BorderBasisError`; passing `strict=False` returns a diagnostic result instead when a supporting Groebner basis could be computed.

Two exact construction paths are available. The default `algorithm="groebner"`
path is Gröbner-derived: it uses a zero-dimensional Gröbner basis to obtain a quotient order ideal, reduces each border monomial into that quotient basis, and verifies the commuting multiplication-matrix criterion. The exact Macaulay
`algorithm="linear"` path forms exact Macaulay matrices, row-reduces polynomial multiples over the rational domain, extracts a divisor-closed quotient order ideal from nonpivot monomials, and reads border relations directly from the row space. A supporting Groebner basis is still used to certify zero-dimensionality and quotient dimension, but not to reduce the border monomials in the linear constructor.

This gives CAD/QE and zero-dimensional solving code a stable symbolic
border-basis data structure plus a native exact linear-algebra construction
route. It is not yet the numerical AVI/SVD variant used for approximate
vanishing ideals; that remains a natural future backend for floating-point or
empirical data.


### Piecewise and system simplification

`simplify_piecewise` simplifies branch expressions under their branch assumptions, removes unreachable branches after earlier conditions cover the domain, and can return a `PiecewiseSimplificationResult` with diagnostics. `simplify_system` records simple equality substitutions in `SimplifiedSystem.substitutions`, supports `eliminate_equalities=True`, and can return either a formula, a tuple of constraints, or the structured result via `output=`.

### Structured sign proofs

`prove_positive`, `prove_nonnegative`, `prove_negative`, and
`prove_nonpositive` still return booleans by default. Pass
`return_result=True` to get a `SignProofResult` with a proof method,
certificate, and counterexample when the claimed strict sign fails.

```python
prove_nonnegative(x**2 + y**2, [x, y])
# True

prove_positive(x**2, [x], return_result=True).counterexample
# {x: 0}
```
