# Introduction to semialgebraic geometry

Semialgebraic geometry studies subsets of real Euclidean space that can be described by **finite Boolean combinations of polynomial equalities and inequalities**. A typical semialgebraic set is

\[
S=\{(x,y)\in\mathbb R^2: x^2+y^2\le 1,\; y\ge x\},
\]

but unions, complements, strict inequalities, disequalities, and quantified descriptions are allowed as well. For example,

\[
\{x\in\mathbb R:\exists y\;(y^2=x\land y\ge0)\}
\]

is semialgebraic. The fact that eliminating the quantified variable still gives a semialgebraic set is one of the central structural results of the subject.

`semialg` treats these objects computationally. Its basic language is a SymPy Boolean formula over real polynomial relations:

```python
import sympy as sp
from semialg import SemialgebraicRegion, is_satisfiable

x, y = sp.symbols("x y", real=True)
formula = (x**2 + y**2 <= 1) & (y >= x)
region = SemialgebraicRegion(formula, (x, y))

is_satisfiable(formula, (x, y))
# True
```

## Why this class of sets is special

Semialgebraic sets are broad enough to encode many exact geometric and optimization problems, but rigid enough to support strong finiteness theorems. Among their key properties:

- they are closed under finite union, intersection, difference, and complement;
- they are closed under Cartesian product;
- they are closed under polynomial maps and inverse images;
- most importantly, they are closed under coordinate projection;
- they have finitely many connected components;
- their dimension and topology are tame compared with arbitrary subsets of \(\mathbb R^n\);
- many logical questions about them are decidable over the real closed field \(\mathbb R\).

This combination of expressive power and tameness is what makes exact quantifier elimination, cylindrical algebraic decomposition, global polynomial optimization, and constructive topology possible.

## Algebraic versus semialgebraic sets

An **algebraic set** is cut out only by polynomial equations,

\[
V(f_1,\ldots,f_m)=\{x:f_1(x)=\cdots=f_m(x)=0\}.
\]

A semialgebraic set also permits inequalities and Boolean structure. Thus the unit circle

\[
x^2+y^2=1
\]

is algebraic, while the closed disk

\[
x^2+y^2\le1
\]

is semialgebraic but not algebraic.

Computationally this distinction matters. Algebraic techniques such as Gröbner bases, resultants, subresultants, rational univariate representations, and Jacobian rank tests are important building blocks, but inequalities require real-root ordering, sign determination, cells, and real quantifier elimination.

## Formulas and geometry

A semialgebraic formula can be viewed in two complementary ways:

1. **Logical:** a first-order formula over the ordered field of real numbers.
2. **Geometric:** a region consisting of all real points satisfying that formula.

`semialg` deliberately supports both viewpoints. Raw formulas are useful for one-off decisions and quantifier elimination; `SemialgebraicRegion` packages a formula with its ambient variables; `CADRegion` attaches a reusable cylindrical algebraic decomposition for repeated geometric and topological queries.

See [Region representations](../guides/region_representations.md) for the computational tradeoffs.

## A small example

Consider the annulus

\[
1\le x^2+y^2\le4.
\]

From this one description we can ask logically or geometrically different questions:

```python
from semialg import SemialgebraicRegion

annulus = SemialgebraicRegion(
    (x**2 + y**2 >= 1) & (x**2 + y**2 <= 4),
    (x, y),
)

annulus.contains((sp.Rational(3, 2), 0))
annulus.dimension()
annulus.euler_characteristic()
annulus.measure()
```

Projection, connected components, optimization, integration, and meshing are all operations on the same underlying semialgebraic object.

## Where to go next

- [Concepts and theory](theory.md) introduces quantifier elimination, Tarski-Seidenberg, dimension, stratification, and real algebraic topology.
- [Algorithms and techniques](algorithms.md) explains CAD, projection/lifting, virtual substitution, algebraic solving, topology, optimization, and integration.
- [Applications](applications.md) surveys geometry, robotics, verification, control, optimization, scientific modeling, and theorem-assisted computation.
- [Exactness and certification](../concepts/exactness_and_certification.md) explains what the package certifies and where numerical layers enter.
