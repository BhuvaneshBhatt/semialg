# Applications

`semialg.applications` contains thin, domain-oriented workflows built on the certified core algorithms. The application layer translates a practical question into semialgebraic decision, quantifier-elimination, range, or optimization problems and packages the result in a domain-friendly form.

Core mathematical operations such as `function_range`, `semialgebraic_measure`, `integrate_over_region`, `semialgebraic_minimize`, and CAD remain in the core package. They are reusable primitives rather than applications and are intentionally not duplicated under `semialg.applications`.

## Robust parameter and tolerance analysis

Use `robust_parameter_analysis()` when a model contains operating variables and design or uncertain parameters.

```python
import sympy as sp
from semialg.applications import robust_parameter_analysis

x, a, b = sp.symbols("x a b", real=True)

result = robust_parameter_analysis(
    x**2 + a*x + b >= 0,
    operating_variables=[x],
    parameters=[a, b],
)

result.robust_condition
# a**2 - 4*b <= 0
```

The result distinguishes three exact parameter regions:

- `feasible_condition`: at least one operating point satisfies the constraints;
- `robust_condition`: every operating point satisfies the constraints;
- `violation_condition`: at least one operating point violates the constraints.

For a single condition, `robust_parameter_region(..., quantifier="forall")` returns the robust parameter region directly. Use `quantifier="exists"` for existential feasibility.

For tolerance or bounded-operating analyses, pass `operating_domain`. Universal robustness then means that the requirement holds for **every operating point in that domain**, not for every real value. For example:

```python
x, a = sp.symbols("x a", real=True)
result = robust_parameter_analysis(
    x**2 <= a,
    [x],
    [a],
    operating_domain=sp.And(x >= -1, x <= 1),
)
result.robust_condition
# a >= 1
```

## Certified symbolic-math validation

The validation helpers use exact semialgebraic reasoning rather than numerical sampling.

```python
from semialg.applications import validate_identity

x = sp.symbols("x", real=True)
result = validate_identity((x + 1)**2, x**2 + 2*x + 1, [x])
assert result.valid
```

A false claim can include an exact counterexample when the decision layer can construct one:

```python
result = validate_identity(x**2, x, [x])
assert not result.valid
print(result.counterexample)
```

Available helpers include:

- `validate_identity()` for expression identities, optionally under assumptions;
- `validate_formula_equivalence()` for proposed solution sets or logical transformations;
- `validate_range()` for a proposed exact function-range formula.

These APIs are useful for validating CAS transformations, symbolic equation or inequality results, range computations, and test oracles for other symbolic systems.

## Exact optimization benchmark oracle

`exact_optimization_benchmark()` builds an exact reference result for a low-dimensional polynomial optimization problem.

```python
from semialg.applications import (
    exact_optimization_benchmark,
    validate_numeric_optimization,
)

x = sp.symbols("x", real=True)
benchmark = exact_optimization_benchmark(x**2, [x >= 1], [x])

benchmark.exact_value
# 1

check = validate_numeric_optimization(
    benchmark,
    numeric_value=1.000001,
    atol=1e-5,
)
assert check.within_tolerance
```

The benchmark retains the full `OptimizationResult`, including exact optimizer points, attainment, method diagnostics, and the global-certification flag. The numeric check is intentionally lightweight: it compares the reported numerical objective value against the exact optimum. Solver-specific feasibility and KKT diagnostics can be layered on top when a numerical solver exposes them.

## Why core operations do not live here

The application layer should remain thin. For example:

- `function_range()` is a general exact image/range primitive used by many applications;
- `semialgebraic_measure()` is a general geometric primitive;
- `semialgebraic_minimize()` and `semialgebraic_maximize()` are general optimization primitives;
- CAD and QE are foundational reasoning engines.

Duplicating these APIs in `semialg.applications` would create two public homes for the same operation and make documentation and maintenance less clear. Applications should compose the core APIs instead.

## Polynomial control stability regions

Use `polynomial_stability_analysis()` for strict continuous-time stability of a real polynomial characteristic equation. The application builds the Hurwitz matrix and returns the exact semialgebraic parameter condition under which every characteristic root lies in the open left half-plane.

```python
import sympy as sp
from semialg.applications import polynomial_stability_analysis

s, a, b = sp.symbols("s a b", real=True)
result = polynomial_stability_analysis(s**2 + a*s + b, s, [a, b])
result.condition
# (a > 0) & (a*b > 0), equivalent to (a > 0) & (b > 0)
```

`result.determinants` exposes the leading principal Hurwitz determinants. `polynomial_stability_region()` is the convenience form when only the parameter condition is needed.

The current scope is **strict Hurwitz stability** for real polynomial characteristic equations. Marginal stability, discrete-time Schur/Jury stability, delay systems, and nonlinear state-space stability are separate problems and are not implied by this API.

## Polynomial safety and invariant verification

`verify_polynomial_invariant()` checks an inductive invariant for a discrete polynomial transition system. The core consecution obligation is

\[
D(x) \land I(x) \Longrightarrow I(F(x)).
\]

Optional initial and unsafe conditions add the obligations

\[
D(x) \land I_0(x) \Longrightarrow I(x)
\]

and

\[
D(x) \land I(x) \Longrightarrow \neg U(x).
\]

```python
from semialg.applications import verify_polynomial_invariant

x = sp.symbols("x", real=True)
check = verify_polynomial_invariant(
    x >= 0,
    {x: x + 1},
    [x],
    initial_condition=sp.Eq(x, 0),
    unsafe_condition=x < 0,
)
assert check.valid
```

When an implication fails, the result retains an exact counterexample whenever the decision layer can construct one. This API verifies a **supplied** invariant; discovering invariants or barrier certificates is a different and substantially harder synthesis problem. The current workflow is state-only and rejects undeclared symbolic parameters rather than silently quantifying them.

## Polynomial response-surface analysis

`analyze_response_surface()` packages exact optimization and range analysis for polynomial surrogate or response-surface models.

```python
from semialg.applications import analyze_response_surface

x, y = sp.symbols("x y", real=True)
domain = sp.And(x >= -1, x <= 1, y >= -1, y <= 1)
analysis = analyze_response_surface(
    x**2 + y**2,
    [x, y],
    domain=domain,
    thresholds=[1],
)

analysis.minimum.value   # 0
analysis.maximum.value   # 2
analysis.gradient        # (2*x, 2*y)
```

The result contains exact minimum and maximum `OptimizationResult` objects, an exact `FunctionRangeResult`, the symbolic gradient and stationary condition, and exact superlevel-set formulas for requested thresholds. The application is intended for polynomial response surfaces with exact numeric/algebraic coefficients; undeclared symbolic coefficient parameters are rejected. Use the core parameter-stratified APIs when coefficient parameters must remain symbolic. Arbitrary statistical or machine-learning models should first be converted to a polynomial surrogate if an exact semialgebraic analysis is desired.


## Lyapunov-function verification

`verify_lyapunov_function()` checks a supplied polynomial Lyapunov function for a continuous-time polynomial vector field. It certifies that the candidate vanishes at the chosen equilibrium, that the equilibrium belongs to the requested domain, that the candidate is positive away from the equilibrium, and that its Lie derivative has the requested negative sign.

```python
import sympy as sp
from semialg.applications import verify_lyapunov_function

x = sp.Symbol("x", real=True)
result = verify_lyapunov_function(x**2, {x: -x}, [x])
assert result.valid
result.lie_derivative
# -2*x**2
```

By default the derivative condition is strict away from the equilibrium. Set `derivative_strict=False` to verify the weaker nonpositive derivative condition. The application verifies a supplied candidate; it does not synthesize Lyapunov functions. Undeclared symbolic parameters are rejected rather than silently quantified.

## Barrier-certificate verification

`verify_barrier_certificate()` verifies a polynomial barrier for a continuous-time polynomial vector field using the convention that `B <= 0` is the certified safe side. It proves three exact obligations:

\[
I(x) \Rightarrow B(x)\le 0,
\]

\[
U(x) \Rightarrow B(x)>0,
\]

and, on the barrier boundary,

\[
B(x)=0 \Rightarrow \dot B(x)\le0.
\]

```python
from semialg.applications import verify_barrier_certificate

x = sp.Symbol("x", real=True)
result = verify_barrier_certificate(
    x**2 - 1,
    {x: -x},
    [x],
    initial_condition=sp.Eq(x, 0),
    unsafe_condition=x**2 >= 4,
)
assert result.valid
```

Use `derivative_strict=True` when the desired sufficient condition requires a strictly negative boundary derivative. This API verifies a supplied barrier; barrier synthesis is outside its current scope.

## Certified polynomial sensitivity and monotonicity

`analyze_polynomial_sensitivity()` differentiates an exact polynomial model with respect to each predictor, computes the exact derivative range on the domain, and certifies derivative-sign conditions.

```python
from semialg.applications import analyze_polynomial_sensitivity

x, y = sp.symbols("x y", real=True)
domain = sp.And(x >= 0, x <= 2, y >= -1, y <= 1)
result = analyze_polynomial_sensitivity(x**2 + 3*y, [x, y], domain=domain)

result.directions[x].classification
# 'nondecreasing'
result.directions[y].classification
# 'strictly_increasing'
```

Classifications are `constant`, `strictly_increasing`, `strictly_decreasing`, `nondecreasing`, `nonincreasing`, or `mixed`. They are certified from derivative signs. On disconnected domains, interpret these as coordinate-wise derivative statements along line segments that remain inside the domain rather than as an ordering claim between arbitrary disconnected points.

## Constraint redundancy analysis

`analyze_constraint_redundancy()` determines whether each explicit constraint is implied by all the others.

```python
from semialg.applications import analyze_constraint_redundancy

x = sp.Symbol("x", real=True)
result = analyze_constraint_redundancy([x >= 2, x >= 1, x <= 5], [x])
result.redundant_indices
# (1,)
```

For a nonredundant constraint, the result retains a counterexample to the implication when one can be constructed. Indices refer to the input constraint sequence.

## Feasible-set diagnostics

`diagnose_feasible_set()` provides a higher-level explanation of a constraint system. For a feasible system it returns an exact witness plus redundancy information. For an infeasible system it can compute an inclusion-minimal conflicting subset of the supplied constraints.

```python
from semialg.applications import diagnose_feasible_set

x = sp.Symbol("x", real=True)
result = diagnose_feasible_set([x > 3, x < 2, x**2 <= 100], [x])
assert not result.feasible
result.conflict_indices
# (0, 1)
```

The reported conflict is **irreducible**: removing any reported member makes that selected conflict satisfiable. It is not guaranteed to have minimum cardinality. Conflict extraction uses repeated exact satisfiability checks and is intended primarily for modest diagnostic constraint sets.

## Application testing contract

The application namespace has a dedicated regression guard: every exported application function must be called directly by at least one behavioral test. Each application area has positive, negative/edge, and invalid-input or certification tests where those contracts apply. This complements the package-wide public-function coverage guard.

## Polynomial model comparison

`compare_polynomial_models()` compares two exact polynomial models over a semialgebraic domain. It reports the exact minimum and maximum of the signed difference, the maximum absolute discrepancy, dominance in either direction, and exact counterexamples to failed dominance or equivalence claims when available.

```python
import sympy as sp
from semialg.applications import compare_polynomial_models

x = sp.Symbol("x", real=True)
domain = sp.And(x >= 0, x <= 1)
result = compare_polynomial_models(x**2, x, [x], domain=domain)

result.first_le_second          # True
result.maximum_absolute_error  # 1/4
```

The maximum absolute error is obtained by maximizing the polynomial square

\[
(f-g)^2
\]

and then taking its exact nonnegative square root. This keeps the optimization problem polynomial rather than introducing `Abs` into the optimizer. The current application expects exact polynomial coefficients and rejects undeclared symbolic parameters.

## Parameter regime analysis

Parameter-regime analysis exposes exact qualitative changes as certified semialgebraic strata rather than sampled parameter values.

`analyze_parameter_regimes()` partitions parameter space by real solvability of a semialgebraic system:

```python
from semialg.applications import analyze_parameter_regimes

x, a = sp.symbols("x a", real=True)
regimes = analyze_parameter_regimes(sp.Eq(x**2 + a, 0), [x], [a])

regimes.select({a: -1})  # True
regimes.select({a: 1})   # False
```

`analyze_root_count_regimes()` partitions parameter space by the number of distinct real roots of a polynomial:

```python
from semialg.applications import analyze_root_count_regimes

root_regimes = analyze_root_count_regimes(x**2 + a, x, [a])
root_regimes.select({a: -1})  # 2
root_regimes.select({a: 0})   # 1
root_regimes.select({a: 1})   # 0
```

Both workflows wrap `ParameterStratifiedResult`; each branch has an exact guard, and `select()` evaluates the active regime under exact parameter substitution. The initial scope covers solvability and real-root-count regimes. More specialized regime quantities can be added without changing the result model.

## Polynomial probability

`polynomial_probability()` computes an exact normalized probability for a polynomial density on a semialgebraic support. The density need not already integrate to one.

```python
from semialg.applications import polynomial_probability

x = sp.Symbol("x", real=True)
result = polynomial_probability(
    x <= sp.Rational(1, 2),
    [x],
    density=2*x,
    bounds={x: (0, 1)},
)

result.normalizing_mass  # 1
result.event_mass        # 1/4
result.probability       # 1/4
```

Before integrating, semialg certifies that the polynomial density is nonnegative on the effective support. The support must have finite positive total mass. `geometric_probability()` is the uniform-density convenience API:

```python
from semialg.applications import geometric_probability

uniform = geometric_probability(
    x <= sp.Rational(1, 2),
    [x],
    support=sp.And(x >= 0, x <= 1),
)
uniform.probability  # 1/2
```

The probability application inherits the exact integration engine's scope. A polynomial density can still lead to a transcendental exact value such as `pi`; exact does not mean algebraic. Unsupported integrals are declined according to the core integration contract rather than silently sampled numerically.
