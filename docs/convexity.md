# Convexity

`semialg.is_convex(region, variables)` decides convexity of real semialgebraic
sets exactly.  The implementation uses progressively more general certified
methods and keeps the quantified segment definition as the complete fallback.

## Decision hierarchy

The current hierarchy is:

1. normalize the semialgebraic formula;
2. handle trivial false/true formulas and equality-defined singleton candidates;
3. decide one-dimensional convexity completely by exact CAD connectivity;
4. recognize intersections of affine equalities and inequalities;
5. recognize convex basic quadratic sets using exact constant-Hessian signs;
6. certify basic polynomial sublevel/superlevel intersections from global
   Hessian signs;
7. certify polynomial constraints relative to an already-certified convex
   domain;
8. reject cheaply certifiable disconnected sets using topology;
9. search CAD cell samples for an exactly verified segment counterexample;
10. decide the remaining cases by complete quantifier elimination of the
    segment definition.

The final fallback checks whether there exist `x`, `y` in the set and
`0 <= t <= 1` for which `t*x + (1-t)*y` leaves the set.  A negative answer to
that existential query is an exact convexity proof.

## Structured certificates

`convexity_certificate(region, variables)` returns a `ConvexityCertificate`
with the exact Boolean outcome, the method that decided it, optional witness
data, and method-specific details.  This is useful when downstream packages
need to distinguish a cheap structural proof from full quantified QE.

`quadratic_convexity_certificate(region, variables)` recognizes basic
quadratic intersections.  Equality constraints must be affine.  Quadratic
sublevels require a positive-semidefinite Hessian and quadratic superlevels a
negative-semidefinite Hessian.  If the representation is outside that class,
the routine returns an inconclusive certificate rather than guessing.

`polynomial_convexity_certificate(expr, variables, domain=..., sense=...)`
certifies the Hessian sign of a polynomial globally or on a semialgebraic
domain. Constant symmetric Hessians (in particular, quadratic polynomials) are
decided by exact congruence/LDL-style inertia computation, avoiding exponential
principal-minor enumeration. Nonconstant polynomial Hessians use the exact
principal-minor characterization; symbolically unresolved signs are sent to
complete CAD under the requested domain.

A `False` result from `polynomial_convexity_certificate` means that the
requested Hessian sign fails somewhere on the supplied domain.  It is a
statement about that Hessian certificate, not by itself a representation-
independent proof that every possible sublevel representation is nonconvex.
The Hessian condition is interpreted relative to the supplied domain, so it
holds vacuously on an empty domain.  On lower-dimensional domains a failed
ambient Hessian test can still be inconclusive about convexity of the restricted
function; callers should use `True` certificates as sufficient proofs and avoid
treating a failed Hessian sign as a set-convexity theorem.


## Function convexity

`function_convexity(expr, variables, domain=...)` classifies a supported real semialgebraic function as `affine`, `convex`, `concave`, `neither`, `nonconvex_domain`, or `unknown`. The natural real domain recognized by `function_domain` is intersected automatically with the explicit domain.

The implementation consumes the reusable backend primitives rather than duplicating them. It certifies domain convexity, uses affine-relative strict feasibility and safe affine presolve for lower-dimensional domains, uses `function_sign` and `matrix_definiteness` for Hessian analysis, and falls back to exact epi/hypograph or quantified Jensen reasoning through the semialgebraic graph layer. Parameterized calls return `ParameterStratifiedResult` branches with the same canonical classifications.

For a polynomial on a full-dimensional convex domain, the Hessian test is necessary and sufficient, so both positive and negative matrix-definiteness results are decisive. On lower-dimensional domains, a positive Hessian certificate is still sufficient, but a failed ambient Hessian sign is not treated as a nonconvexity proof; exact function-definition reasoning is used instead.

Examples:

```python
from sympy import Abs, symbols
from semialg import function_convexity

x, y, a = symbols("x y a", real=True)

assert function_convexity(x**2, [x]) == "strongly_convex"
assert function_convexity(-(x**2), [x]) == "strongly_concave"
assert function_convexity(2 * x + 1, [x]) == "affine"
assert function_convexity(Abs(x), [x]) == "convex"

restricted = function_convexity(-(y**2), [x, y], domain=(y == 0) & (x >= -1) & (x <= 1))
assert restricted == "affine"

parametric = function_convexity(a * x**2, [x], parameters=[a])
assert parametric.select({a: 2}) == "convex"
assert parametric.select({a: -2}) == "concave"
assert parametric.select({a: 0}) == "affine"
```

With `return_result=True`, nonparametric calls return `FunctionConvexityResult`, including the separate convex/concave decisions, domain certificate, matrix certificates where applicable, relative-interior/presolve diagnostics, and exact counterexample data when a definition-based failure is found. Unsupported transcendental graphs are reported as `unknown`; semialg does not guess from floating-point samples.

## Domain-relative Hessian certificates

Global convexity of a defining polynomial is stronger than necessary in many
intersections.  The staged domain-relative path begins with affine constraints,
which define a convex domain, and adds a nonlinear constraint only after its
Hessian sign has been proved throughout the domain already certified convex.
For example, the polynomial

```python
f = x**4 - x**2 + y**2
```

is not globally convex, but its Hessian is positive semidefinite on `x >= 1`.
Consequently an intersection such as

```python
And(x >= 1, f <= 10)
```

can be certified convex without constructing the full two-point segment query.
The order is proof-sensitive: a later constraint is never used to justify the
convexity of an earlier domain.

## Negative certificates

Convex nonempty sets are connected, so a certified disconnected topology is an
exact nonconvexity certificate.  The hierarchy uses cheap topology checks where
they are advantageous, especially factorized algebraic varieties.

For other cases, the segment search uses exact CAD cell sample points.  If two
verified feasible points have a verified infeasible midpoint, the returned
certificate records the endpoints, midpoint, and `t = 1/2`.  Failure to find
such a witness is never interpreted as convexity; the complete quantified
fallback remains authoritative.
