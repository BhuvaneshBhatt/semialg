# Algorithms and techniques

Semialgebraic computation combines symbolic algebra, real-root algorithms, logic, geometry, and topology. `semialg` uses specialized methods first when their hypotheses are certified and retains complete CAD-based fallbacks where supported.

## Normalization and preprocessing

Before an expensive exact solve, formulas are normalized so mathematically equivalent input shapes expose the same structure. Useful preprocessing includes:

- canonical polynomial residuals for relations;
- primitive and square-free normalization;
- safe factor simplification;
- affine equality substitution;
- Boolean simplification;
- variable-incidence decomposition;
- detection of equational constraints;
- parameter exceptional-polynomial analysis.

The goal is not cosmetic simplification. Good preprocessing can reduce the dimension, number of active polynomials, projection burden, and number of CAD cells dramatically.

## Univariate real-root isolation

Many higher-dimensional algorithms ultimately rely on exact univariate root operations:

- counting distinct real roots;
- isolating algebraic roots in rational intervals;
- comparing algebraic numbers;
- determining polynomial signs at algebraic samples;
- preserving root order as parameters vary over a CAD base cell.

These operations support exact CAD sections, interval decomposition, root classification, and algebraic sample points.

## Resultants, subresultants, and discriminants

Eliminating a variable from polynomial equations often uses the resultant. Discriminants detect multiple-root events, while principal/subresultant coefficients capture degree and common-root changes. In CAD projection these quantities identify parameter values where fiber behavior may change.

`semialg` caches projection and polynomial-key computations because these exact algebraic operations can dominate repeated workloads.

## Cylindrical algebraic decomposition

CAD is a general-purpose exact engine for real polynomial logic.

### Projection

Starting with polynomials in \(n\) variables, a projection operator constructs lower-dimensional polynomials whose signs control root behavior in the eliminated coordinate. Different projection operators trade projection-set size against hypotheses such as well-orientedness.

### Base decomposition

The lowest-dimensional line is decomposed at the real roots of the projected polynomials.

### Lifting

Each lower-dimensional cell is lifted by isolating roots of the next-level polynomials over a sample point. Those roots define sections and sectors. Repeating this process yields a decomposition of \(\mathbb R^n\).

### Truth evaluation

When all relevant polynomial signs are invariant on each cell, the original formula is truth-invariant there. Quantifiers can then be evaluated by finite cell logic.

`CADRegion` retains the resulting decomposition so point location, topology, Boolean combination, integration, and meshing can reuse it.

## Equational constraints

If a conjunction contains an equation \(f=0\), solutions need only lie on that locus. Equational-constraint-aware projection can reduce projection polynomials and lifting work. Choosing a useful equation is itself an algorithmic problem; `semialg` scores candidates using projection burden rather than only syntactic degree.

## Virtual substitution

For supported low-degree quantified problems, especially quadratic ones, **virtual substitution** can eliminate variables without building a full CAD. It substitutes symbolic representatives of roots and boundary limits into the remaining conditions while preserving exact logical meaning.

This is often much cheaper than CAD and is therefore tried before complete decomposition when its applicability is certified.

## Zero-dimensional algebraic solving

When polynomial equations define finitely many complex/real points, algebraic solving can be substantially cheaper than CAD. Techniques include:

- Gröbner bases;
- rational univariate representations;
- shape-position methods;
- multiplication matrices;
- exact root isolation and sign filtering.

Inequalities are then checked exactly at the algebraic candidate points.

## Boolean reconstruction on a shared CAD

Once a sign-invariant CAD has been built, many Boolean simplifications can be phrased as truth masks over the same finite cells. `semialg` uses this for implication minimization and can generalize DNF branches to larger implicants whenever their CAD truth mask remains inside the target set. This avoids rebuilding pairwise CADs for every logical check.

## Exact optimization

The package uses a hierarchy rather than one universal optimizer:

1. box/polytope specializations;
2. separable Cartesian-product decomposition;
3. exact univariate critical-point calculus;
4. equality reduction;
5. KKT/active-set algebraic candidates;
6. positive-dimensional critical-locus handling;
7. exact CAD decision or complete function-range certification.

The important design principle is that a finite list of candidates is not automatically a global proof. Certification checks that no better feasible point exists.

## Region topology and cell complexes

An adapted CAD gives a finite exact combinatorial representation. From it one can construct:

- cell dimensions;
- closure incidence;
- adjacency graphs;
- connected components;
- local dimension;
- Euler characteristic.

`CADCellComplex` records unsigned incidence. Oriented boundary maps, homology, and Betti numbers are outside the supported topology API.

## Exact integration

A full-dimensional cylindrical CAD cell has nested bounds

\[
a_1<x_1<b_1,
\quad a_2(x_1)<x_2<b_2(x_1),\ldots
\]

which directly define an iterated integral. `semialg` combines such CAD extraction with cheaper recognizers for intervals, boxes, simplexes, radial regions, and ellipses. Intrinsic integration over lower-dimensional regular strata additionally requires metric/Jacobian factors.

## Numerical geometry over exact structure

Meshing and delineable curve/surface evaluation are numerical layers built from exact CAD branch data. Root indices and CAD adjacency come from exact computation; sampled coordinates are approximate.

This separation matters: a conforming sampled mesh is useful for visualization and downstream numerical work, but it is not automatically a proof of ambient isotopy to the exact region. See [Exact versus numerical regions](../guides/exact_vs_numerical_regions.md).

## Complexity

General quantifier elimination over real closed fields is expensive, and CAD has doubly exponential worst-case behavior in the number of variables. This is why variable order, decomposition, equational constraints, specialized solvers, caching, and conservative early exits matter so much in practice.

The package's strategy is therefore **exact-first, structure-aware, and fallback-driven**, not “always build a full CAD.”
