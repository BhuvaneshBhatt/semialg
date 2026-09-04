# Solving and sampling semialgebraic systems

`semialg` provides an exact solver for real systems of polynomial equations,
inequalities, and Boolean combinations of them.  `solve_semialgebraic` is the
main structured entry point.

## `solve_semialgebraic`

```python
import sympy as sp
from semialg import solve_semialgebraic

x, y = sp.symbols("x y", real=True)

sol = solve_semialgebraic(
    [x**2 + y**2 <= 1, x > 0, y > 0],
    [x, y],
)

sol.satisfiable
# True

sol.sample
# one satisfying sample point, when available
```

A `SemialgebraicSolution` contains:

- `formula`
- `variables`
- `satisfiable`
- `sample`
- `samples`
- `method`
- `diagnostics`

## Automatic exact method planner

`method="auto"` uses a structural planner.  The order is deliberately biased
toward cheap exact transformations and treats general CAD/QE as the final
fallback:

1. normalize algebraic/domain-sensitive constraints;
2. eliminate globally safe affine equalities and retain reconstruction rules;
3. split small explicit Boolean disjunctions branch-by-branch;
4. split conjunctive systems into independent variable-incidence blocks;
5. use the exact univariate interval/root reducer;
6. decide coupled affine systems by Fourier--Motzkin elimination;
7. recognize already-triangular/cylindrical bound descriptions;
8. compute exact equality-ideal dimension when a finite algebraic system is
   structurally possible; for zero-dimensional ideals, record the exact
   quotient-algebra dimension and use rational-univariate representation (RUR);
9. simplify zero-dimensional residual constraints by Groebner normal forms and
   exact radical/nonvanishing consequences before enumerating algebraic points;
10. use a bounded Groebner presolve to simplify positive-dimensional equality
    varieties and reduce inequalities modulo the equality ideal;
11. use the complete CAD/QE stack when none of the structural methods settles
    the problem.

The planner does not approximate and does not change the meaning of the
solution set.  Affine elimination is reconstructed in the returned coordinates.
For example:

```python
sol = solve_semialgebraic(
    [sp.Eq(x + y, 1), x >= 0, y >= 0],
    [x, y],
    count=0,
)

sol.formula
# Eq(x, 1 - y) & (y >= 0) & (y <= 1)
```

A coupled linear system is decided without CAD:

```python
z = sp.symbols("z", real=True)
linear = solve_semialgebraic(
    [x >= 0, y >= 0, z >= 0, x + y + z <= 1],
    [x, y, z],
    count=0,
)

linear.method
# 'linear_fourier_motzkin'
```

Independent blocks are solved separately.  Thus a system such as
`x**2 <= 1` together with `y**2 >= 4` is reduced by two one-dimensional
solves instead of one two-dimensional CAD.

### Planner diagnostics

`solution.diagnostics["planner_steps"]` records the methods considered by the
automatic planner, whether each was accepted, and why.  Nested
`child_plans`/`branch_plans` record incidence and Boolean decomposition.  The
system profile also records polynomial degree, equality/inequality counts,
linearity, candidate finite-dimensionality, variable blocks, and a suggested
sparse variable order.

These diagnostics are deterministic structural data; they are suitable for
performance-regression tests without depending on wall-clock thresholds.

### Equality-ideal analysis for finite systems

The finite-system backend shares an exact equality-ideal context.  A
grevlex Groebner basis supplies the leading monomial ideal, from which semialg
computes the exact Krull dimension.  For dimension zero, standard monomials
give the exact dimension of the quotient algebra.  The same context can lazily
compute per-coordinate univariate eliminants by FGLM, certify global polynomial
replacements of the form `x = p(other_variables)`, and test radical membership
with the Rabinowitsch construction.

These facts are used before real point filtering.  For example, `x**2 == 0`
together with `x != 0` is rejected because `x` belongs to the radical of the
equality ideal; no algebraic-root enumeration or CAD is required.  Remaining
inequalities are still evaluated exactly at RUR points, so the optimization
does not weaken solution verification.

The expert algebraic API exposes `EqualityIdealContext` and
`analyze_equality_ideal` from `semialg.algebraic`.

## Explicit methods

The public `method` selector can still force a solver family.  In particular:

- `method="linear"` requires an affine conjunctive real system and uses exact
  affine presolve plus Fourier--Motzkin feasibility;
- `method="rur"` requires a finite zero-dimensional algebraic system;
- `method="interval"` requires one solve variable;
- `method="cad"` / `"qe"` bypass the automatic structural solver and invoke
  the general exact reduction stack.

## Parameter conditions

When `parameters=[...]` is supplied, `solve_semialgebraic` returns exact
parameter-space solvability information through `parameter_conditions` and
`parameter_decomposition`.  The convenience output
`output="conditions"` returns just the parameter condition.

## Sampling and sign evaluation

```python
from semialg import sample_point, sample_points, sign_at, sign_vector

pt = sample_point(sp.And(x > 0, x < 1), [x])
sign_at(x - sp.Rational(1, 2), pt)

sign_vector([x, x - 1], {x: sp.Rational(1, 2)})
# (+, -) in the package's sign representation
```

## Current scope

The solver emphasizes exact real-semialgebraic solution sets, feasibility,
parameter conditions, structural decomposition, and representative samples.
It does not promise a globally human-minimal formula.  Positive-dimensional
systems can naturally remain as exact equations/inequalities or cylindrical
cells rather than being converted to point replacement rules.
