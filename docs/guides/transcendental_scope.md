# Transcendental and algebraic solving scope

semialg distinguishes three mechanisms that can look similar at the surface but have different proof obligations.

## Polynomial and semialgebraic solving

The core exact machinery works with polynomial equalities and inequalities over real closed fields: CAD, quantifier elimination, rational-univariate representations, exact algebraic roots, sign determination, and related certified operations. The [Certified algebraic roots](certified_algebraic_roots.md) guide describes the explicit polynomial interval tools.

## Exact algebraization

Some expressions containing functions such as absolute values, radicals, rational powers, or selected structured transcendental forms can be transformed into an equivalent semialgebraic problem by introducing auxiliary variables and exact side conditions. Algebraization is only valid when the transformation records enough constraints to preserve the original real semantics.

An algebraized problem is subsequently solved by semialgebraic machinery. The algebraization step is not itself a general-purpose transcendental root finder.

## Existing transcendental solver

`semialg.solve.transcendental` is a separate, existing solver subsystem for its explicitly supported transcendental fragment. It may use structure such as periodicity, monotonic pieces, or certified interval reasoning according to the contracts documented in [Transcendental methods](../transcendental_methods.md).

The package does not treat every expression containing `sin`, `cos`, `exp`, or `log` as automatically solvable. Unsupported forms should be declined rather than silently converted into an uncertified numerical answer.

<!-- semialg-exec -->
```python
import sympy as sp

from semialg.algebraic import certify_polynomial_root_interval

x = sp.Symbol("x")
try:
    certify_polynomial_root_interval(sp.sin(x) - x / 2, -1, 1, var=x)
except (TypeError, sp.PolynomialError):
    rejected_by_polynomial_api = True
else:
    rejected_by_polynomial_api = False

assert rejected_by_polynomial_api
```

This rejection does **not** mean the transcendental solver subpackage is absent. It means the polynomial certificate API keeps a narrower contract.

<!-- semialg-exec -->
```python
import importlib.util

assert importlib.util.find_spec("semialg.solve.transcendental") is not None
```

## Design rule

Choose the mechanism according to the mathematical representation:

1. use semialgebraic APIs when the problem is already polynomial/semialgebraic;
2. use exact algebraization only when an equivalence-preserving transformation is implemented;
3. use the existing transcendental solver only for its documented supported fragment;
4. otherwise report that the exact method is unsupported instead of using numerical evidence as certification.

This page is the canonical scope statement. Other documentation pages should link here rather than defining their own competing boundary.
