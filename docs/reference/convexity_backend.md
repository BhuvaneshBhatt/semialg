# Convexity backend primitives

These APIs expose exact semialgebraic reasoning primitives intended for higher-level packages such as symbolic optimization systems. They do not implement KKT logic, active-set policy, duality, or convex-expression composition rules.

## Family contract

**Mathematical return.** `function_sign`, `prove_zero`, and `prove_nonzero` decide scalar sign properties on semialgebraic regions. `function_sign` automatically intersects the supplied assumptions with the function's exact real domain. `matrix_definiteness`, `matrix_psd_on`, and `matrix_pd_on` decide exact symmetric-matrix definiteness by constant exact inertia or determinantal criteria. `matrix_rank_on` and `matrix_rank_stratification` reason about exact rank. `strict_feasible` and `affine_relative_interior_formula` provide affine-relative strict-feasibility semantics. `parametric_affine_reduction` branches on symbolic affine pivots instead of assuming them generically nonzero.

**Exactness and certification.** All decisions are exact. Failed nonparametric sign and matrix-definiteness claims retain exact counterexample points when the satisfiability backend produces one; matrix definiteness also attempts to recover an exact violating quadratic-form vector after specialization. Parameter-dependent queries return `ParameterStratifiedResult` objects with semialgebraic guards. No floating-point sign sampling is accepted as proof.

**Algorithm selection.** Constant symmetric matrices use exact congruence/LDL inertia. Polynomial symmetric matrices use principal minors for semidefiniteness and Sylvester leading minors for definiteness, with exact semialgebraic feasibility/QE on violations. Rank uses determinantal minors. Relative strict feasibility preserves explicit affine equalities and inequalities that are identically tight on that affine hull. Parametric affine reduction produces separate zero/nonzero pivot branches.

**Complexity and limitations.** Matrix minor counts grow combinatorially with dimension, and parameterized conditions may invoke complete CAD. Relative-interior support is restricted to conjunctive affine polynomial systems; arbitrary nonlinear relative interiors remain a separate problem. Parametric affine reduction handles affine pivots whose coefficients depend only on declared parameters and leaves more complicated pivots unreduced.

### Scalar sign queries

```python
from semialg import prove_nonnegative, prove_nonzero, prove_zero, function_sign

classification = function_sign(x**2, [x])
assert classification == "nonnegative"
assert prove_zero(x - x, [x])
assert prove_nonzero(x**2 + 1, [x])
```

The existing `prove_positive`, `prove_nonnegative`, `prove_negative`, and `prove_nonpositive` functions also accept `parameters=`. A parameterized proof returns a Boolean `ParameterStratifiedResult` rather than collapsing free parameters to generic assumptions.

### Symmetric matrix definiteness

```python
from semialg import matrix_definiteness, matrix_pd_on, matrix_psd_on

assert matrix_psd_on([[x**2, 0], [0, 1]], [x])
assert matrix_pd_on([[x**2 + 1, 0], [0, 1]], [x])
```

`matrix_definiteness` accepts `requested="positive_semidefinite"`, `"positive_definite"`, `"negative_semidefinite"`, or `"negative_definite"`. It is the reusable backend used by polynomial Hessian convexity certification.

### Rank and parameter strata

`matrix_rank_on` returns an integer only when the rank is constant over the requested region. `matrix_rank_stratification` partitions parameter space by exact determinantal rank and keeps exceptional determinant-zero loci separate.

### Relative strict feasibility

`affine_relative_interior_formula` turns genuine affine inequality boundaries into strict inequalities while preserving equality constraints and inequalities that are identically tight on the affine hull. `strict_feasible` decides nonemptiness of that relative-interior formula and can return an exact witness.

### Parametric affine reduction

`parametric_affine_reduction` never divides by an unresolved symbolic pivot unconditionally. For an equality such as `a*x + 1 == 0`, it returns a generic `a != 0` branch with the reversible substitution `x -> -1/a` and an exceptional `a == 0` branch.

### Reusable contexts

`SemialgebraicContext` exposes `function_sign`, `matrix_definiteness`, `strict_feasible`, `implies`, `add_constraints`, and `with_formula`. Derived contexts share the same `ExactComputationContext`, so related queries can reuse process-local exact caches.


## Function convexity

`function_convexity` is the high-level consumer of these primitives. It first intersects the explicit domain with the exact real domain recognized by `function_domain`, certifies that the resulting domain is convex, and then uses the cheapest exact route that is sufficient:

1. affine equality presolve and relative-interior reasoning for lower-dimensional affine domains;
2. exact Hessian sign analysis through `function_sign` and `matrix_definiteness` for polynomial and graph-supported smooth algebraic functions;
3. exact epi/hypograph convexity after domain-sensitive algebraization when that produces a direct semialgebraic set; and
4. the quantified Jensen definition with semialgebraic function graphs as the representation-independent fallback.

```python
from sympy import Abs, symbols
from semialg import function_convexity

x, a = symbols("x a", real=True)

assert function_convexity(x**2, [x]) == "strongly_convex"
assert function_convexity(-(x**2), [x]) == "strongly_concave"
assert function_convexity(3 * x + 1, [x]) == "affine"
assert function_convexity(Abs(x), [x]) == "convex"

parametric = function_convexity(a * x**2, [x])
assert parametric.select({a: 2}) == "strongly_convex"
assert parametric.select({a: -2}) == "strongly_concave"
assert parametric.select({a: 0}) == "affine"
```

The primary result is the strongest certified canonical classification. In addition to `affine`, `convex`, `concave`, `neither`, `nonconvex_domain`, and `unknown`, it can return `strictly_convex`, `strictly_concave`, `strongly_convex`, and `strongly_concave`. With `return_result=True`, nonparametric calls return `FunctionConvexityResult`, which records the separate convex/concave decisions, strict and strong curvature certificates, any certified strong-convexity modulus, domain certificate, Hessian certificates when used, relative-interior/presolve information, and any exact counterexample obtained by the definition fallback. Passing `properties="all"` additionally requests quasi-/strict-quasi-, pseudo-, and algebraically provable log-convex/log-concave properties. Univariate quasiconvexity uses the exact monotonicity partition; pseudoconvexity uses differentiability plus the exact stationary-global-minimum criterion; log curvature uses the algebraic matrix `f*H(f) - grad(f)*grad(f).T` after proving `f > 0`. Calls with remaining free symbols return `ParameterStratifiedResult` automatically (or equivalently when `parameters=` is supplied), with the same classifications as branch values.

A positive Hessian certificate is used as a sufficient proof on any certified convex domain. A negative Hessian result is treated as necessary only when the function is polynomial and the domain has nonempty ambient interior. This avoids incorrectly rejecting functions that are affine or convex only after restriction to a lower-dimensional domain.


## Function monotonicity

`function_monotonicity(expression, variable, domain=...)` classifies exact univariate monotonicity on the intersection of the supplied domain with the function's exact real domain. It uses derivative sign reasoning on convex one-dimensional domains, using rational zero-structure directly to distinguish weak from strict monotonicity when possible. Disconnected domains are analyzed by an exact monotonicity partition and cross-component certificates before any selective pairwise graph/QE fallback.

The strongest canonical outputs are `constant`, `strictly_increasing`, `increasing`, `strictly_decreasing`, `decreasing`, `nonmonotonic`, and `unknown`. In particular, isolated critical points do not incorrectly destroy strictness: `x**3` is certified `strictly_increasing` even though its derivative vanishes at zero.

```python
from sympy import Abs, symbols
from semialg import function_monotonicity

x, a = symbols("x a", real=True)

assert function_monotonicity(x**3, x) == "strictly_increasing"
assert function_monotonicity(x**2, x) == "nonmonotonic"
assert function_monotonicity(Abs(x), x, domain=x >= 0) == "strictly_increasing"

conditional = function_monotonicity(a * x, x)
assert conditional.select({a: 2}) == "strictly_increasing"
assert conditional.select({a: 0}) == "constant"
assert conditional.select({a: -2}) == "strictly_decreasing"
```

When symbolic quantities other than the monotonicity variable remain free, they are treated as parameters automatically and the result is an exact `ParameterStratifiedResult`. This is the same condition-generation policy used by `function_convexity`. The pairwise definition is used for parameter conditions, so strictness is expressed exactly rather than by the merely sufficient condition `f'(x) > 0`. With `return_result=True`, unconditional calls return `FunctionMonotonicityResult` containing the derivative/sign evidence and any exact pairwise counterexample.


## Function monotonic partition

`function_monotonic_partition(f, x, domain=...)` computes an exact univariate
partition using a sign-invariant CAD for `diff(f, x)`.  Sector cells are
classified as strictly increasing/decreasing or constant and isolated critical
sections are absorbed when they do not destroy the stronger neighboring
property.  Disconnected function domains remain separate.  Free symbols other
than `x` are treated as parameters automatically; a parameter-first CAD returns
a `ParameterStratifiedResult` whose values are exact partitions.

## Function convex partition

`function_convex_partition(f, x, domain=...)` is the analogous exact univariate
partition based on a sign-invariant CAD for `diff(f, x, 2)`.  Cells are
classified as convex, concave, or affine.  Supported semialgebraic/nonsmooth
expressions fall back through function-graph algebraization when their derivative
sign cannot be represented directly.  Parameter-dependent inflection geometry
is represented by parameter strata rather than generic assumptions.


## Function sign partition

`function_sign_partition(f, x, domain=...)` decomposes the exact real function
domain into connected regions on which the sign is respectively `positive`,
`zero`, or `negative`. It uses the same function-graph algebraization as
`function_sign`, so supported nonsmooth and algebraic expressions are handled
without numerical sampling. Disconnected natural-domain components remain
separate. Free symbols other than `x` are parameters automatically, and a
parameter-first CAD returns a `ParameterStratifiedResult` whose branch values
are exact sign partitions.

```python
from sympy import Abs, symbols
from semialg import function_sign_partition

x, a = symbols("x a", real=True)

assert tuple(kind for kind, _ in function_sign_partition(x * (x - 1), x)) == (
    "positive",
    "zero",
    "negative",
    "zero",
    "positive",
)
assert tuple(kind for kind, _ in function_sign_partition(Abs(x), x)) == (
    "positive",
    "zero",
    "positive",
)

conditional = function_sign_partition(a * x, x)
assert tuple(kind for kind, _ in conditional.select({a: 1})) == ("negative", "zero", "positive")
```


## Function smoothness

`function_smoothness(f, variables, domain=..., max_order=...)` reports exact continuity, the greatest certified `C^k` order, whether the function is smooth, and exact exceptional loci. Polynomial and rational functions are certified `C^∞` on their exact real function domains. Univariate `Abs`, `sign`, and finite `Piecewise` joins are checked by exact breakpoint analysis with one-sided branch limits; for example `Abs(x)` is continuous but not `C^1` exactly at `x = 0`. The structured result type is `FunctionSmoothnessResult`. For direct predicates, `is_function_continuous(...)` returns the certified continuity projection and `is_function_smooth(...)` returns the certified smoothness projection.

## Function mapping properties

`function_mapping_properties(mapping, variables, domain=..., codomain=...)` certifies injectivity, surjectivity onto the stated semialgebraic codomain, and bijectivity. It also records the exact image formula, collision witnesses for failed injectivity, and missing-value witnesses for failed surjectivity. Polynomial/rational maps reuse the expert formula-level `semialgebraic_image` primitive; supported algebraic maps use function-graph elimination as fallback. Full-space affine maps use exact Jacobian rank as a fast path. The structured result type is `FunctionMappingPropertiesResult`. The convenience predicates `is_injective(...)`, `is_surjective(...)`, and `is_bijective(...)` return the corresponding certified projections without discarding the aggregate analysis API.


## Function-property performance and cache reuse

The function-property layer shares lazy analysis products within each high-level query. Natural/effective domains, derivatives, gradients, Hessians, sign classifications, smoothness checks, and the univariate monotonicity partition are reused rather than recomputed by each derived property.

One-dimensional set convexity uses ordered CAD-cell contiguity rather than closure-connectivity. This matters for punctured domains: `Ne(x, 0)` is not convex, and consequently `function_monotonicity(1/x, x)` is correctly classified as `"nonmonotonic"` globally even though the function is strictly decreasing on each connected component.

Univariate convexity first classifies the exact sign of the second derivative when that derivative is semialgebraic. Polynomial/rational sign partitions use one CAD over numerator/denominator boundaries. Strong one-dimensional curvature uses the exact range of the second derivative where needed. `properties="all"` propagates implication chains and computes at most one monotonicity partition for all quasi-curvature questions. Scalar univariate mapping properties reuse strict monotonicity as an injectivity certificate before pairwise collision QE. These are proof-preserving short-circuits; unsupported transcendental cases such as `exp(x)` remain `"unknown"` in the semialgebraic core.
