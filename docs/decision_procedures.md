# Decision procedures and inequality proving

The decision layer answers exact real-domain questions about formulas.

## Satisfiability and equivalence

```python
import sympy as sp
from semialg import is_satisfiable, is_tautology, implies, equivalent

x, y = sp.symbols("x y", real=True)

is_satisfiable(sp.And(x**2 + y**2 <= 1, x > 0, y > 0), [x, y])
# True

is_tautology(sp.Or(x < 0, x >= 0), [x])
# True

implies(x > 1, x**2 > 1, [x])
# True

equivalent(x**2 <= 1, sp.And(x >= -1, x <= 1), [x])
# True
```

## Inequality proving

The inequality provers reduce sign claims to satisfiability checks:

- `prove_nonnegative(f)` checks unsatisfiability of `f < 0`.
- `prove_positive(f)` checks unsatisfiability of `f <= 0`.
- `prove_nonpositive(f)` checks unsatisfiability of `f > 0`.
- `prove_negative(f)` checks unsatisfiability of `f >= 0`.

```python
from semialg import prove_positive, prove_nonnegative

prove_nonnegative((x - 1)**2, [x])
# True

prove_positive(x**2 + 1, [x])
# True

prove_nonnegative(x*y, [x, y], assumptions=sp.And(x >= 0, y >= 0))
# True
```

## Notes

These functions work over the real domain and are intended for polynomial and semialgebraic formulas. They are often used internally by simplification, region predicates, and optimization routines.

## Specialized algebraic feasibility and polynomial sign backends

For rational polynomial equality systems, `real_algebraic_feasibility` provides a
certified positive-dimensional alternative to CAD. It first obtains a replayable
certified equidimensional decomposition, constructs Aubry--Rouillier--Safey El Din
polar/critical systems with deterministic rational generic points, proves a strict
dimension drop at every recursive step, and delegates the resulting finite systems
to the exact rational-univariate solver. A genericity or decomposition failure is
reported as an incomplete result; it is never converted into an emptiness claim.
`solve_real_algebraic_set` is the convenience form that returns one exact witness,
`None` for certified emptiness, and raises when certification is incomplete.

The decomposition engine also uses triangular initial and separant splitters. Each
split is implemented by the exact identity
`V(I) = V(I + <h>) union V(I : h^infinity)` and is accepted only after radical
containment proves both branches strict. `initial_split` and `separant_split` are
recorded in the same replayable decomposition certificate as factor-incidence,
monomial, real-radical, and saturation strategies.

`zeng_negative_point` is the specialized exact polynomial-negativity backend. Fast
witness searches may prove that `f < 0` at an exact point. For the certified
coercive critical-value fragment of the Zeng semidefinite-polynomial strategy, the
backend proves attainment of the global minimum, solves the zero-dimensional
gradient locus exactly, and checks every critical value exactly. Its inexpensive
coercivity certificate is deliberately narrow: every top-degree monomial must
have an even exponent in every variable and a nonnegative rational coefficient,
and every variable must occur in a positive pure top-degree monomial. Mixed-odd
leading terms are never accepted by this certificate. Unsupported cases remain
incomplete. `polynomial_nonnegative` and `find_negative_point` provide
strict convenience APIs; the former is also consumed by global `prove_nonnegative`
and `prove_nonpositive` before general implication/QE.

`find_negative_witness_fast` contains the one-sided accelerators. Odd-degree/ray
searches and deterministic pseudo-random rational lines are allowed to establish
satisfiability only after exact substitution verifies the returned witness. Failure
to find a witness has no logical meaning and never proves infeasibility. The
univariate step chooses rational points from certified isolating intervals; it
does not manufacture witnesses from floating-point root approximations. Search
counts and coefficient bounds are validated at the public boundary.
