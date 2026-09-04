# Concepts and theory

This page summarizes the mathematical ideas behind semialgebraic computation. It is not intended as a substitute for a text in real algebraic geometry, but it gives the vocabulary needed to understand `semialg`'s algorithms and guarantees.

## Basic semialgebraic sets

A **basic closed semialgebraic set** has the form

\[
\{x\in\mathbb R^n:f_i(x)=0,\;g_j(x)\ge0\}.
\]

Allowing strict inequalities, disequalities, and finite Boolean combinations gives the general class of semialgebraic sets. Disjunctive normal form can express such a set as a finite union of basic pieces, although explicit DNF expansion can be exponentially expensive and should not be treated as a computational default.

## Real closed fields and quantifier elimination

The first-order theory of real closed fields admits **quantifier elimination**. Every formula constructed from polynomial relations, Boolean connectives, and real quantifiers is equivalent to a quantifier-free formula.

For example,

\[
\exists y\;(y^2=x)
\]

is equivalent over the reals to \(x\ge0\).

This gives decision procedures for statements such as:

- whether a region is empty;
- whether one region is contained in another;
- whether a polynomial is nonnegative on a region;
- whether a parameter value admits a solution;
- whether a proposed optimum is globally valid.

In `semialg`, `Exists`, `ForAll`, `reduce_formula`, and `qe_by_complete_cad` expose this logical layer.

## Tarski-Seidenberg theorem

The geometric form of quantifier elimination is the **Tarski-Seidenberg theorem**: the image of a semialgebraic set under a coordinate projection is semialgebraic.

If

\[
S\subseteq\mathbb R^{n+m},
\]

then

\[
\pi(S)=\{x\in\mathbb R^n:\exists y\in\mathbb R^m\;(x,y)\in S\}
\]

is semialgebraic.

This theorem underlies exact projection, polynomial images, parameter elimination, reachable-set calculations, and many optimization constructions. `SemialgebraicRegion.project()`, `.image()`, and `.preimage()` expose these operations directly.

## Cells and cylindrical decomposition

A **cell decomposition** partitions a semialgebraic set into finitely many simple pieces. In a cylindrical algebraic decomposition, cells are arranged recursively so that projection onto lower coordinates is either a section or sector over an existing cell.

Informally, over a base cell in variables \((x_1,\ldots,x_{k-1})\), the next coordinate \(x_k\) is either:

- a **section**, lying on a continuous algebraic root function; or
- a **sector**, lying strictly between two such root functions.

A sign-invariant CAD ensures that each relevant polynomial has constant sign on every cell. Consequently every Boolean formula built from those polynomial signs has constant truth value on each cell.

## Dimension

Semialgebraic dimension behaves much like geometric dimension. A point has dimension 0, a curve usually dimension 1, a surface dimension 2, and a solid subset of \(\mathbb R^3\) dimension 3.

A CAD gives a constructive route: each sector contributes one free coordinate, while each section does not. The maximum cell dimension in a selected region is its semialgebraic dimension.

Local dimension can vary. For a union of an isolated point and a curve, the isolated point has local dimension 0 while generic points on the curve have local dimension 1. `SemialgebraicRegion.local_dimension(point)` computes this through the exact CAD cell complex.

## Closure, interior, boundary, and regularity

Semialgebraic sets are stable under standard topological operations. Their closure, interior, and boundary are again semialgebraic. This allows exact computation of notions such as:

\[
\partial S=\overline S\setminus S^\circ.
\]

A set is **regular closed** when \(S=\overline{S^\circ}\), and **regular open** when \(S=(\overline S)^\circ\). These regularizations are useful in constructive solid geometry and in eliminating lower-dimensional artifacts.

## Connectedness and tame topology

Semialgebraic sets have finitely many connected components, and for semialgebraic sets connectedness and path connectedness coincide componentwise. CAD adjacency therefore provides a finite combinatorial object from which connectivity questions can be answered.

More generally, semialgebraic sets admit triangulations and finite stratifications. Their homotopy type is finite in a strong sense. This is the foundation for Euler characteristics, homology, Betti numbers, roadmaps, and eventually certified simplicial models.

## Euler characteristic

Given a finite cell decomposition with \(c_k\) cells of dimension \(k\), the compactly supported semialgebraic Euler characteristic is

\[
\chi_c(S)=\sum_k(-1)^k c_k.
\]

For compact semialgebraic sets this agrees with the ordinary Euler characteristic. Additivity is especially useful:

\[
\chi(A\cup B)=\chi(A)+\chi(B)-\chi(A\cap B).
\]

`CADCellComplex` and the region Euler API compute this invariant from exact CAD cells.

## Singular and regular points

For an algebraic hypersurface \(f=0\), a point is singular when

\[
f=0,\qquad \nabla f=0.
\]

For systems of equations, Jacobian-rank conditions generalize this test. Semialgebraic boundaries introduce additional subtleties: a corner formed by intersecting two smooth inequality boundaries is not necessarily an algebraic singularity of either boundary polynomial. `semialg` documents this distinction explicitly in its singular-locus API.

## Positivity and optimization

Questions of polynomial positivity over semialgebraic sets connect real algebraic geometry with optimization. A global minimum statement

\[
f(x)\ge m\quad\text{for all }x\in S
\]

is itself a quantified semialgebraic assertion. KKT systems, discriminants, algebraic critical points, CAD, and exact range computation provide complementary methods for proving such statements.

Algebraic positivity certificates such as Positivstellensatz representations form another major branch of the theory. `semialg` emphasizes exact decision/CAD and structural optimization certificates rather than a full sums-of-squares/Positivstellensatz engine.
