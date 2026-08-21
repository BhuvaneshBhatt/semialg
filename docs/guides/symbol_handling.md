# Symbol handling

SymPy symbols are identified by more than their printed names. Assumptions are part of their identity.

```python
import sympy as sp

x = sp.Symbol("x")
x_real = sp.Symbol("x", real=True)

x == x_real
# False
```

This distinction matters in exact symbolic algorithms.

## Prefer actual symbols

The clearest API style is:

```python
x, y = sp.symbols("x y", real=True)
formula = x**2 + y**2 <= 1

from semialg import is_satisfiable
is_satisfiable(formula, [x, y])
```

## String names are resolved, not recreated blindly

Many public APIs also accept names such as `"x"`. semialg resolves those names against the actual symbols already present in the formula, objective, integrand, region mapping, or other relevant context. It should preserve the original assumptions and identity.

This avoids creating a new `Symbol("x", real=True)` that is not the same object as an existing unassumed `Symbol("x")`.

## Ambiguous names are errors

If the context contains distinct symbols with the same printed name, a string cannot identify which one you mean. semialg rejects that ambiguity instead of guessing.

Pass the exact `Symbol` object in such cases.

## Parameters follow the same rule

Parameter APIs resolve string parameter names against the actual expression symbols. Returned parameter-stratified results should therefore use the same symbols you supplied in the mathematical problem.

## Bounds must name declared variables

Shared bound normalization rejects a bound whose key does not resolve to a declared variable. This catches misspellings and accidental extra symbols instead of silently ignoring them.

```python
# Conceptually invalid: y is not one of the declared integration variables.
# bounds={y: (0, 1)} with variables=[x]
```

## Do not compare formulas by printed text

Two expressions can print the same while containing assumption-distinct symbols. semialg's decision APIs use symbolic identity and exact decision procedures rather than `str()`/`sstr()` equality as proof of equivalence.

## Recommended practice

- Create symbols once and pass them through your workflow.
- Use strings for convenience only when names are unambiguous.
- Pass exact `Symbol` objects in reusable libraries or parameterized code.
- Keep assumptions intentional and consistent.
