# Solving and sampling reference
## Family contract

**Mathematical return.** Solving APIs return exact witnesses, finite algebraic solution representations, or samples from semialgebraic solution sets.

**Exactness and certification.** Exact sampling/witness paths use algebraic values and certified sign checks. Explicit numerical sampling modes are inexact by design and are not proof substitutes.

**Algorithm.** Depending on structure, solving may use virtual-substitution witnesses, RUR/zero-dimensional solving, or CAD-derived cells.

**Complexity and limitations.** A witness is evidence of satisfiability, not a complete description of a positive-dimensional solution set. Use region/CAD APIs when the full set matters.




## Primary API overview

This table is the substantive coverage target for the primary APIs assigned to this reference page. Each entry states the API's primary role; the family contract and detailed sections below explain shared algorithms, exactness guarantees, and limitations. It is maintained together with `docs/reference/primary_api_manifest.toml`, and documentation tests require every root-level primary API to map here rather than merely appearing in the generated public index.

| API | Kind | Role / return |
|---|---|---|
| `discretize_region_geometry` | function | Return lightweight geometry for explicit standard-region objects. |
| `discretize_solution` | function | Return a small plotting/discretization representation for a solution. |
| `find_instance` | function | Find satisfying assignments for a formula. |
| `RealAlgebraicFeasibilityResult` | class | Certified SAT/UNSAT/incomplete result for real polynomial equality feasibility. |
| `real_algebraic_feasibility` | function | Decide rational polynomial equality feasibility with certified ARS polar reduction. |
| `solve_real_algebraic_set` | function | Return one exact real point or certified emptiness for a polynomial equality set. |
| `WitnessSearchResult` | class | One-sided exact-verified witness-search result. |
| `find_negative_witness_fast` | function | Try odd-degree and rational-line negative witnesses with exact verification. |
| `is_zero_dimensional` | function | Return whether rational polynomial equations define a finite complex set. |
| `plot_region_geometry` | function | Plot an explicit standard-region object using Matplotlib. |
| `plot_solution` | function | Plot a 1D/2D solution using Matplotlib when available. |
| `sample_point` | function | Return one satisfying real sample point for ``formula``, or ``None``. |
| `sample_points` | function | Return satisfying real sample points for a quantifier-free formula. |
| `sign_at` | function | Return the sign of a polynomial/expression at a point. |
| `sign_vector` | function | Return the signs of ``polys`` at ``point`` in input order. |
| `solve_zero_dimensional_system` | function | Solve a finite rational polynomial system exactly. |

## `solve_semialgebraic`

```text
solve_semialgebraic(
    constraints, variables=None, *, parameters=None, domain="reals",
    count=1, samples=None, sample_mode=None, strategy=None, method="auto",
    variable_order=None, projection_order=None, normalize_domains=True,
    return_formula=False, output=None
)
```

Solves a real semialgebraic system and returns a `SemialgebraicSolution` by default. Depending on `output`/`return_formula`, it can expose formulas, cells, samples, or other structured representations.

`method="auto"` may use specialized exact methods such as rational-univariate representation for finite zero-dimensional equality systems before falling back to broader semialgebraic methods.

## `sample_point` and `sample_points`

Representative sampling is exact by default. Supported workflows include representative/automatic, rational, grid, random, and CAD-cell sampling. Numerical random sampling is explicitly opt-in with `exact=False`.

Returned public samples are checked against the original formula.

## `sign_at` and `sign_vector`

Evaluate polynomial/expression signs at exact points, including algebraic and RUR-backed points. Certified exact paths do not use a hidden floating-point sign guess.

## Instance helpers

`find_instance` returns one exact instance mapping when ``count=1`` and a tuple of instance mappings when multiple instances are requested. Set ``return_result=True`` for status, method, approximate candidates, and diagnostics.

`find_instance_formula`, `find_instance_text`, and `component_instances` follow the same direct-return convention. ``component_instances`` returns the tuple of component sample mappings by default; ``return_result=True`` exposes the component objects and CAD diagnostics.

## Edge cases

- `count=0` requests no samples rather than one implicit sample.
- An empty feasible set is distinct from an unsupported solving strategy.
- String variables follow the shared symbol-resolution rules.

See [Solving](../solving.md) for deeper examples.

## Solver configuration, domains, and resolution

`SemialgOptions` stores high-level solver controls and `SolveDomain` names the
supported solving domains. `resolve_formula` evaluates a structured formula to a
Boolean result when its truth can be decided exactly; `reduce_formula` instead
returns a simplified exact formula when symbolic structure should be retained.

### `function_domain(...)`

`function_domain(expr, variables=None)` computes recognized exact real-domain
conditions using the same semialgebraic function-graph machinery used by range
analysis whenever that graph is available.  The shared graph engine
handles rational expressions, `Abs`, `sign`, `Min`, `Max`, finite `Piecewise`
expressions, and rational powers.  Nested graph variables are projected away by
QE when useful, so domain constraints can be simplified back to conditions on
the requested input variables.

Rational powers deliberately follow SymPy's own semantics.  Ordinary
`x**Rational(p, q)` is a principal-branch `Pow`; for a noninteger rational
exponent its real-valued locus on a real base is therefore nonnegative (strictly
positive for negative exponents).  Explicit real roots should be written with
`sympy.real_root`.  SymPy canonicalizes an odd real root to a combination such
as `sign(x)*Abs(x)**(1/3)`, and semialg recognizes that representation as a real
root rather than reinterpreting ordinary `Pow`.

```python
import sympy as sp
from semialg import function_domain

x = sp.symbols("x", real=True)

function_domain(x ** sp.Rational(1, 3), [x])
# x >= 0

function_domain(sp.real_root(x, 3), [x])
# True

function_domain(sp.sqrt(1 - sp.sqrt(x)), [x])
# an equivalent polynomial condition for 0 <= x <= 1

function_domain(sp.sqrt(x - 1) / (x - 2), [x])
# (x >= 1) & Ne(x, 2), modulo exact simplification
```

The structural domain layer also recognizes denominator nonvanishing and
`log(arg) > 0`.  This is useful even though the graph of `log` itself is not
semialgebraic.  Thus the function's scope is: exact semialgebraic graph/domain
reasoning where available, plus selected nonalgebraic heads whose **real-domain
predicate** is semialgebraic.

This is still not a complete real-domain analyzer for arbitrary SymPy functions.
Unsupported domain-sensitive heads may require conditions that are not inferred.
A returned `True` therefore means either that the supported exact graph projects
to all requested real inputs, or that the structural fallback found no additional
restriction; it is not a blanket claim that every possible function head has
been analyzed.

`is_real_valued` combines these recognized domain conditions with assumptions
and exact semialgebraic implication checks where possible.

`is_zero_dimensional` detects finite polynomial solution sets before
`solve_zero_dimensional_system` performs exact solving and returns the tuple of exact solution-coordinate tuples by default. Set ``return_result=True`` to obtain `ZeroDimensionalSolveResult`, including the backend, RUR representation, status, and notes.
