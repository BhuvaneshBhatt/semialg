# Exact rules before CAD

Use the ordinary optimization and region interfaces for structured inputs. The default dispatch recognizes the cases below before building stationary systems, image CADs or closure-incidence CADs. These rules prove the answer from exact polynomial identities, root coverage or quadratic geometry. Parameters or an unavailable exact comparison send the problem to the general pipeline.

| Problem | Exact rule | Supported cases and limits |
| --- | --- | --- |
| Affine objective on an ellipsoid | Support in a positive-definite metric | One quadratic equality or strict/weak sublevel, real algebraic constant coefficients and positive-definite matrix; extra constraints and singular quadrics use other methods |
| Univariate polynomial optimization | Complete derivative/constraint root isolation and endpoint limits | Real algebraic constant coefficients, holes, isolated points, open endpoints and unbounded intervals; unavailable root isolation declines |
| Polynomial zero sets and complements | Closedness, density and empty interior of a nonzero polynomial zero set | Constant real coefficients; inequalities require additional sign geometry |
| Quadratic inequality topology | Exact symmetric inertia and the stationary value | Includes indefinite and singular semidefinite forms, empty sets and full space |
| Cartesian topology | Products of exact one-dimensional closures and interiors | Multivariate conjunctions grouped by a single coordinate, degree at most two and real algebraic coefficients; holes and thin factors are retained |
| Boundedness of a finite union | Coercive sublevels, subsets and finite unions | Positive-definite quadratics or a leading form equal to a positive multiple of an even power of squared norm; no negative conclusion from failed recognition |
| Intersection of principal ideals | The polynomial least common multiple | Polynomial rings over supported exact characteristic-zero fields; other ideals retain elimination |
| QE formula presentation | Existing CAD sign vectors and Boolean truth masks | Complete ambient sign CAD, bounded family size, and uniform projected truth on every realized sign vector; no extra CAD construction |

## Optimization and attainment

For `q(v) = v.T*A*v + b.T*v + c`, with positive-definite `A`, write its sublevel as `(v-center).T*A*(v-center) <= level`. An affine objective `a.T*v+d` has bounds `a.T*center+d ± sqrt(level*(a.T*A.inv()*a))`. Equality in Cauchy–Schwarz gives an exact optimizer. A strict nonconstant sublevel has the same infimum and supremum but neither is attained; a constant objective remains attained throughout a nonempty domain. Negative levels and zero-level strict sublevels are empty.

Univariate polynomial calculus isolates every constraint root and every derivative root. Constraint signs are constant between the cuts, while the cuts themselves are tested separately. Polynomial endpoint values and leading-term limits at infinity, together with stationary values, exhaust the extrema. Missing root coverage cannot produce a certified answer.

<!-- semialg-exec -->
```python
import sympy as sp
from semialg import semialgebraic_minimize

x, y = sp.symbols("x y", real=True)
closed = (x - 2) ** 2 + 4 * (y + 1) ** 2 <= 9
result = semialgebraic_minimize(3 * x + 4 * y + 5, closed, (x, y), return_result=True)
assert sp.simplify(result.value - (7 - 3 * sp.sqrt(13))) == 0
assert result.attained and result.certified
assert sp.simplify(closed.subs(result.point)) is sp.true
punctured = semialgebraic_minimize(x * x, sp.Ne(x, 0), (x,), return_result=True)
assert punctured.value == 0 and not punctured.attained and punctured.certified
```

Sign assumptions on optimization variables are lowered to explicit constraints in neutral real coordinates before dispatch. Result points use the original public symbols, and replay checks the same sign restrictions.

The [executable application](https://github.com/BhuvaneshBhatt/semialg/blob/main/examples/structural_rules.py) also checks open ellipsoids, an unbounded cubic and a three-dimensional box boundary. The ellipsoid support identity is an instance of [support-function calculus](https://www.stanford.edu/~boyd/cvxbook/bv_cvxbook.pdf).

## Topology and boundedness

A polynomial equality defines a closed set, but the zero set is not always the boundary of its weak inequality. For example, `(x-y)**4 >= 0` defines all of real space, with empty boundary. Quadratic inertia and the stationary value distinguish these cases, including thin weak sublevels and empty strict sublevels. Higher-degree inequalities without a structural proof keep semantic CAD as their default path.

One-dimensional interval topology keeps its existing specialized path, which is cheaper than symbolic inequality reduction for repeated closure queries. For Cartesian factors `S_i`, closure and interior commute with the finite product. Boundary is the union of `boundary(S_i)` times the closures of the other factors. This covers deleted interior points and lower-dimensional factors. Explicit `strategy="cad"` requests the general topology path for comparison; `strategy="syntactic"` retains the separate syntactic heuristic contract.

Coercive-sublevel boundedness is a positive proof. An intersection is bounded if one factor is bounded; a finite union is bounded when all branches are bounded. An unrecognized formula still receives the exact general boundedness analysis.

## Reuse without changing the proof

QE presentation can use the projection-polynomial signs already certified on a complete ambient CAD. It first checks that cells with identical complete sign vectors have identical projected truth. Only then can Boolean minimization replace root-function formulas with polynomial-sign formulas. Different truth on the same vector, truth-only reduced CADs and finite-variety decompositions retain their existing presentation path.

Exact algebraic-root presentations and canonical primary-decomposition constructor results use bounded caches with the package's normal lifetime and teardown. Root isolation, sign decisions and certificate replay retain their exact obligations. Clearing caches must change runtime only. See [performance and profiling](guides/performance.md) and [native proof replay](native_cad_proofs.md).

These rules do not remove CAD's worst-case complexity. Forced full CAD tests still construct their requested geometry, and generic high-degree or parameter-dependent problems can remain expensive.
