# Tame topology, triangulation, and roadmaps

A semialgebraic set can be defined by equations with very high degrees and complicated singularities. Yet, unlike arbitrary subsets of Euclidean space, it admits finite descriptions of its topological structure. This chapter explains the mathematical guarantees and distinguishes them from currently implemented `semialg` algorithms.

## Closure, interior, boundary, and local dimension

For \(S\subseteq\mathbb R^n\), the Euclidean closure \(\overline S\) adds limit points, the interior \(S^\circ\) consists of points with a neighborhood inside \(S\), and the boundary is

\[ \partial S=\overline S\setminus S^\circ. \]

All three are semialgebraic if \(S\) is. The closure is essential in this identity: for the open disk \(x^2+y^2<1\), the boundary is the circle \(x^2+y^2=1\), even though none of its points belongs to the original disk.

**Global dimension** is the largest cell dimension; **local dimension** may differ between points. The union of a line segment and a distant isolated point has global dimension 1 but local dimension 0 at the isolated point. Neither dimension should be confused with the dimension of the ambient coordinate space.

## Finitely many connected components

A semialgebraic set has finitely many connected components, and each component is semialgebraic. Moreover, connectedness implies existence of a semialgebraic path between any two points in a component. This is much stronger than arbitrary topological connectedness: wild infinite oscillations and pathological connected-but-not-path-connected examples are excluded from the semialgebraic category.

The component count need not equal the number of irreducible algebraic components. The algebraic set \(xy=0\) is reducible into two lines, but its real zero set is connected because the lines meet at the origin. Conversely, one irreducible real polynomial can have several connected components of its real zero set.

## Triangulation: what the theorem actually promises

The **semialgebraic triangulation theorem** states that semialgebraic sets can be modeled by finite collections of compatible simplices via semialgebraic homeomorphisms. For a **compact** set \(S\), one may take the realization \(|K|\) of a finite simplicial complex \(K\) and a semialgebraic homeomorphism \(|K|\to S\). For a general noncompact or nonclosed set, the correct formulation uses a suitable **union of open simplices**, not necessarily the entire realization of \(K\).

Why does the qualification matter? Every entire finite simplicial complex is compact. It cannot be homeomorphic to the open interval \((0,1)\), although that interval is semialgebraic and trivially triangulable as one open 1-simplex.

Compatibility means that finitely many specified semialgebraic subsets can be represented as unions of simplicial pieces. A circle, disk, punctured disk, and a union of arcs illustrate how topology can be converted to finite incidence data even when algebraic boundaries are curved.

**Important computational distinction:** The existence theorem is general, while the exact triangulation backend in `semialg` supports a narrower certified subset, including documented compact polyhedral and compact one-dimensional settings. General curved semialgebraic triangulation should not be inferred from a statement of the theorem; consult [theory](theory.md) and [feature matrix](../feature_matrix.md).

## Euler characteristic and inclusion–exclusion

For a finite semialgebraic cell decomposition, the **compactly supported Euler characteristic** is

\[ \chi_c(S)=\sum_{k\ge0}(-1)^k\,c_k, \]

where \(c_k\) is the number of dimension-\(k\) cells. It is additive on disjoint semialgebraic pieces. For compact semialgebraic sets, it agrees with the usual Euler characteristic.

Be careful with noncompact sets. One open interval is a single 1-cell, so \(\chi_c((0,1))=-1\), whereas the ordinary homotopy Euler characteristic of an interval is \(+1\). This is not a contradiction: they are different invariants outside the compact setting.

As a compact example, a closed disk has Euler characteristic 1. A closed annulus, homotopy-equivalent to a circle, has Euler characteristic 0. One can also use inclusion–exclusion when suitable sets have known characteristics:

\[ \chi(A\cup B)=\chi(A)+\chi(B)-\chi(A\cap B). \]

The package's CAD cell-complex and Euler interfaces use exact incidence data where their prerequisites are satisfied; see [region invariants](../guides/region_invariants.md).

## Roadmaps: a one-dimensional connectivity skeleton

A **roadmap** of a semialgebraic set is a curve-like subset designed to meet every connected component in a connected way and to support connecting other points to it. In geometric path planning, the roadmap turns a potentially high-dimensional connectivity query into the problem of navigating a lower-dimensional network.

Conceptually, a roadmap construction selects critical loci or slices of a projection, connects them by carefully chosen algebraic curves, and recurses until the resulting network certifies the components. Unlike a visualization mesh or a collection of sample points, a **certified roadmap** has a mathematical reachability property. The existence and complexity of roadmap algorithms are distinct from the general CAD theorem.

**Package scope:** `semialg` provides supported exact connectivity queries and cell-chain paths using its CAD/topology machinery, but does **not** claim a complete implementation of general algebraic roadmap algorithms. See [applications](applications.md) and [limitations](../limitations.md).

## Stratification and singularities

Many semialgebraic sets can be decomposed into finitely many smooth manifold-like **strata** with compatible frontier behavior. This separates regular pieces from singular points, corners, cusps, and lower-dimensional intersections.

For a hypersurface \(f=0\), an algebraic singular point satisfies \(f=0\) and \(\nabla f=0\). For a region given by inequalities, a sharp corner may have an irregular **boundary** even when each defining polynomial is nonsingular at that corner. For example, the quadrant \(x\ge0,y\ge0\) has a corner at the origin although the lines \(x=0\), \(y=0\) are smooth algebraic hypersurfaces. These distinctions matter in optimization, integration, and local geometry.

## Numerical geometry versus certified topology

A plotted sample of an algebraic curve can suggest connectivity, and a triangle mesh may approximate a surface well. Neither automatically proves topological equivalence. Reliable certification requires controlling roots, singularities, and incidence across cell boundaries. `semialg` deliberately distinguishes numerical geometry created from exact structure from a proof of ambient isotopy or a certified homeomorphism. See [exact versus numerical regions](../guides/exact_vs_numerical_regions.md).

## A semialg example

For a compact annulus, the existing example in the [introduction](introduction.md#a-small-example) builds a `SemialgebraicRegion` and asks for its dimension, Euler characteristic and measure. The [topology example](../examples/05_topology_and_components.md) demonstrates the applicable interfaces. These examples should be interpreted according to the documented backend's scope, not as claims of a general triangulation or roadmap engine.

**References:** Bochnak, Coste, and Roy (1998); Basu, Pollack, and Roy (2006). See [bibliography](../references.md#mathematical-background-and-further-reading).

**Next:** [Optimization and positivity](optimization_positivity.md) · [Applications](applications.md).
