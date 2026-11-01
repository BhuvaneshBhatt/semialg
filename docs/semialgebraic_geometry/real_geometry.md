# Real closed fields and semialgebraic geometry

The defining feature of semialgebraic geometry is the ability to talk about both **equations** and **order**. This chapter explains why quantifier elimination and strong finiteness theorems are possible in that setting, before showing how they connect to `semialg`.

## What counts as a semialgebraic set?

A **basic semialgebraic set** is described by finitely many simultaneous polynomial relations such as \(f_i=0\), \(g_j>0\), and \(h_k\ge0\). A general semialgebraic set is a finite union of such sets, equivalently a set defined by a quantifier-free Boolean formula in polynomial relations.

For example, the graph of the absolute-value function is semialgebraic despite involving a nonpolynomial-looking notation:

\[ \{(x,y):y=|x|\}= \{x\ge0,\;y=x\}\cup\{x\le0,\;y=-x\}. \]

The graph \(y=\sin x\) on the entire real line is **not** semialgebraic: a semialgebraic subset of the line must be a finite union of points and intervals, whereas intersecting the sine graph with the algebraic line \(y=0\) gives infinitely many isolated points.

Semialgebraic sets are closed under complements, finite unions and intersections, products, inverse images of polynomial maps, and images of semialgebraic maps. The image statement is not obvious: it follows from quantifier elimination.

## Real closed fields and why they matter

A **real closed field** is an ordered field in which every positive element has a square root and every odd-degree polynomial has a root. The real numbers are the standard example; the field of real algebraic numbers is another. Over any real closed field, polynomial equalities and inequalities satisfy the same first-order sentences with rational coefficients.

Being real closed is quite different from being algebraically closed. The polynomial \(x^2+1\) has no real root; a real closed field does not contain \(i\) compatible with its order. Its quadratic extension obtained by adjoining \(i\) is algebraically closed.

This provides the algebraic backdrop for an exact decision procedure: coefficients, roots, and comparisons can be represented symbolically without relying on approximations to settle strict inequalities.

## Tarski–Seidenberg and quantifier elimination

Let \(S\subseteq\mathbb R^{n+m}\) be semialgebraic. Its coordinate projection

\[ \pi(S)=\{x\in\mathbb R^n:\exists y\in\mathbb R^m\ (x,y)\in S\} \]

is again semialgebraic. This is the **Tarski–Seidenberg theorem**. In logical language, every first-order formula over a real closed field, built from polynomial relations, Boolean operators, and quantifiers, is equivalent to a quantifier-free formula in its free variables.

A simple example is

\[ \exists y\;(y^2=x)\quad\Longleftrightarrow\quad x\ge0. \]

The equivalence captures the whole range of a polynomial map by eliminating a variable. Universal quantification can be handled using negation: \(\forall y\,P\) is equivalent to \(\neg\exists y\,\neg P\). This applies to nested or alternating quantifiers, although their computational cost can be substantial.

### Worked example: projection of a disk above a line

Consider

\[ P(x):=\exists y\ \bigl[x^2+y^2\le1\ \land\ y\ge x\bigr]. \]

For any \(x\in[-1,1]\), the disk allows \(-\sqrt{1-x^2}\le y\le\sqrt{1-x^2}\). The additional constraint permits a solution precisely when

\[ x\le\sqrt{1-x^2}. \]

If \(x\le0\), that inequality is automatic (subject to \(x\ge-1\)). If \(x>0\), squaring is safe because both sides are nonnegative, and gives \(2x^2\le1\). Thus

\[ P(x)\quad\Longleftrightarrow\quad-1\le x\le\frac1{\sqrt2}. \]

The example is elementary enough to solve by hand but illustrates exactly what CAD automates: split into sign-consistent cases, preserve strict/weak boundaries, and eliminate the hidden coordinate without introducing extraneous solutions.

For the package's quantified-input forms and API options, see [CAD and QE](../cad_qe.md) and [quantified expressions](../quantified_expressions.md). Be careful not to confuse a **mathematically valid** first-order formula with a **cheap** CAD computation.

## Finiteness and dimension

Every semialgebraic set has finitely many connected components, each semialgebraic. A semialgebraic connected component is also semialgebraically path connected. A definable curve cannot oscillate in infinitely many new ways near a finite point as \(\sin(1/x)\) does.

The dimension \(\dim S\) can be read from a finite cell decomposition: it is the maximum dimension among cells making up \(S\). For a nonempty open subset of \(\mathbb R^n\), the dimension is \(n\); for a finite nonempty set it is 0. A coordinate projection cannot increase dimension.

Local dimension is sensitive to where one stands. In

\[ S=\{(x,0):0\le x\le1\}\cup\{(3,0)\}, \]

the interval has local dimension 1, while the isolated point \((3,0)\) has local dimension 0. See [topology, stratification, and roadmaps](topology.md).

## Topological consequences

Semialgebraic closure, interior, and boundary are semialgebraic, with the standard boundary identity

\[ \partial S=\overline S\setminus S^\circ. \]

The closure matters. For example \(S=(0,1)\subseteq\mathbb R\) has boundary \(\{0,1\}\), neither of whose points belongs to \(S\).

Every semialgebraic set admits a finite **semialgebraic triangulation** in an appropriate formulation. When the set is compact, it can be represented by the entire realization of a finite simplicial complex under a semialgebraic homeomorphism. For noncompact or nonclosed sets, one works with suitable unions of open simplices (or compatible pairs), rather than claiming a homeomorphism from an entire compact finite complex.

These results justify finite combinatorial approaches to connectivity and Euler characteristic, but they do not imply that an arbitrary triangulation is inexpensive to construct or is available through the package.

## O-minimality: semialgebraic geometry in a larger setting

An expansion of the ordered real field is **o-minimal** when every definable subset of the real line is a finite union of points and intervals. Semialgebraic sets form the smallest such structure containing the real field operations and order. The concept explains the common thread behind finite connected components, well-behaved curves, and cell decomposition.

Some larger o-minimal structures admit additional functions, but that **does not** make their full first-order theories accessible to semialg's polynomial CAD algorithms. General elementary transcendental equations and arbitrary o-minimal structures are outside its exact real-algebraic decision guarantees; see [transcendental scope](../guides/transcendental_scope.md).

## When theory becomes computation

A theorem guarantees that an equivalent quantifier-free description **exists**. To actually find one, a program needs an effective representation of coefficients and a terminating method for the input class. `semialg` offers virtual substitution and other specialized methods on supported instances, partial CAD when logical information is enough, and complete CAD paths with explicit proof conditions. Failure due to a missing exact backend or exhausted resource budget is not a proof that the formula is false.

**References:** Bochnak, Coste, and Roy (1998); Basu, Pollack, and Roy (2006). See the [bibliography](../references.md#mathematical-background-and-further-reading).

**Next:** [Algebraic computation](algebraic_computation.md) · [A worked CAD](cad_worked_example.md).
