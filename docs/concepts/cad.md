# Cylindrical algebraic decomposition: the core idea

Cylindrical algebraic decomposition (CAD) partitions real space into finitely many connected semialgebraic cells on which the relevant polynomial signs, and therefore the truth of the target formula, are controlled.

## A two-dimensional example

Consider the unit disk

\[
x^2+y^2\le1.
\]

Projecting the boundary polynomial with respect to \(y\) identifies the critical \(x\)-values

\[
x=-1,\qquad x=1.
\]

The \(x\)-axis is decomposed into

```text
(-infinity,-1) | {-1} | (-1,1) | {1} | (1,infinity)
```

Over the middle interval, the boundary has two ordered root functions,

\[
y=-\sqrt{1-x^2},
\qquad
y=\sqrt{1-x^2}.
\]

Lifting over that base interval produces sectors and sections:

```text
y >  sqrt(1-x^2)        sector
y =  sqrt(1-x^2)        section
-sqrt(...) < y < sqrt(...)  sector  <-- disk interior
y = -sqrt(1-x^2)        section
y < -sqrt(1-x^2)        sector
```

At \(x=\pm1\), the two boundary branches coalesce. Exact root identity and multiplicity handling are therefore important during lifting and reconstruction.

## Sections and sectors

A **section** lies on a root of a projection/lifting polynomial.

A **sector** lies between adjacent roots, or above/below all roots.

Each cell has an exact representative sample. Truth/sign information at the sample is valid on the whole cell only because the projection/lifting construction establishes the required invariance.

## Projection and lifting

CAD has two broad stages:

1. **Projection** constructs lower-dimensional polynomials whose roots capture where higher-dimensional root behavior can change.
2. **Lifting** recursively decomposes fibers over lower-dimensional cells.

This is why variable ordering matters: changing which variable is projected first changes the projection polynomials and can change the cost dramatically. The effect is not merely cosmetic: different orders can change projection degrees, coefficient growth, the number of real roots that must be isolated, and the number of cells by very large factors. In Collins-style lifting order, the final variable is projected/eliminated first.

A useful order is therefore part of the algorithm, not just an output-format choice. `semialg` uses structural heuristics and bounded cost estimates when it is free to choose an order, and CAD-driven integration can search coordinate permutations to obtain simpler certified iterated bounds. Quantifier alternation constrains this freedom: variables may not be moved across `Exists`/`ForAll` blocks merely to improve performance because doing so changes the formula.

## Quantifier elimination

Once a truth-invariant CAD has been built, quantified formulas can be evaluated cellwise and truth can be propagated through quantifier blocks. The selected free-variable cells are then reconstructed as an equivalent quantifier-free semialgebraic formula.

## Reduced CAD and equational constraints

An equation known to be logically necessary can permit less projection/lifting work. `semialg` uses reduced/EC-aware paths only when the required conditions are certified. It does not treat an arbitrary resultant or projection polynomial as a logical equality constraint.

## CAD is exact but not cheap

General real quantifier elimination has severe worst-case complexity. `semialg` therefore tries presolve, quadratic virtual substitution, RUR/finite solving, and structural shortcuts before or around CAD when applicable.

See [How semialg chooses an algorithm](algorithm_selection.md) and the [Performance guide](../guides/performance.md).
