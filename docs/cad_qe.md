# CAD and quantifier elimination

Many `semialg` operations reduce questions to first-order formulas over the real numbers. For example:

- satisfiability: `exists x. C(x)`;
- implication: `for all x. A(x) implies B(x)`;
- range: `exists x. C(x) and t = f(x)`;
- optimization level sets: `exists x. C(x) and f(x) <= t`;
- subset checks: there is no point satisfying `A(x) and not B(x)`.

Cylindrical algebraic decomposition and quantifier elimination provide exact answers for formulas involving polynomial equations and inequalities over real closed fields.

## Variable ordering and cost

CAD is sensitive to variable ordering and has doubly exponential worst-case complexity. `semialg` therefore uses CAD/QE as a correctness engine, but it also includes faster special cases for common low-dimensional problems, standard shapes, and univariate systems.

## Design implication

A useful pattern throughout the package is:

```text
syntactic preprocessing
  -> fast special case if recognized
  -> CAD/QE-backed semantic check or projection
  -> conservative failure if unsupported
```
