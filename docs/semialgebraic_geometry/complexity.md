# Complexity, decidability, and practical scope

The decisive advantage of semialgebraic geometry is not that every problem is easy. It is that **first-order polynomial questions over the reals are decidable**, and many useful subproblems admit specialized exact algorithms. This chapter explains what that statement does—and does not—promise.

## Decidable does not mean computationally inexpensive

Tarski's quantifier-elimination theorem provides an algorithm for deciding statements built from polynomial equalities, inequalities, Boolean operations, and quantifiers over real closed fields. In particular, one can ask whether a region is empty, whether two exact semialgebraic formulas are equivalent, or whether a polynomial has a prescribed global lower bound.

The underlying input size includes much more than the number of variables: polynomial degrees, coefficient bit lengths, number of polynomials, Boolean structure, number and alternation of quantifier blocks, and details of the desired output. Two formulas in the same number of variables can differ dramatically in runtime.

For general CAD, worst-case complexity is **doubly exponential in the number of variables**: symbolic projection and repeated lifting can create extremely large families of polynomials and cells. A typical qualitative bound is of the form \((sd)^{2^{O(n)}}\) for \(s\) input polynomials of degree bounded by \(d\), suppressing important coefficient and output-size factors. This does not say that every instance is expensive or that a particular input will have that many cells.

Specialized algorithms for existential real feasibility and selected topological questions can have substantially better asymptotic bounds than full sign-invariant CAD. Those algorithms are not interchangeable with complete quantified formula elimination, particularly with unrestricted alternations of existential and universal quantifiers.

## The existential theory of the reals

The **existential theory of the reals** asks whether a finite system of polynomial equations and inequalities has a real solution. It can be formulated as the truth of a sentence of the form

\[ \exists x_1\cdots\exists x_n\;\Phi(x_1,\ldots,x_n), \]

with quantifier-free polynomial Boolean formula \(\Phi\). Its associated complexity class is written \(\exists\mathbb R\). This class describes computational problems polynomial-time reducible to real feasibility; it lies between NP and PSPACE under the usual reductions and has applications including geometric realization and constraint problems.

A finite satisfying witness need not have rational coordinates. Even \(x^2=2\) has no rational solution. Exact real-algebraic numbers, isolating intervals, and sign certificates are therefore natural components of a real-feasibility solver.

The existential theory of the **integers** is a very different subject: unrestricted polynomial Diophantine solvability over the integers is undecidable. Results about the real closed field must not be transferred to integers without extra hypotheses.

## Why algebraic structure matters

Several common shortcuts illustrate the principle of solving only the needed problem:

- **Linear or low-degree formulas:** linear elimination and supported virtual substitution can avoid a full CAD.
- **Equational constraints:** projection may need to preserve signs only on a necessary equation's zero set, provided the required hypotheses are certified.
- **Univariate problems:** Sturm sequences and algebraic root isolation replace high-dimensional decomposition.
- **Zero-dimensional algebraic systems:** Gröbner bases, rational univariate representations, and exact root filtering can be effective.
- **Recognized geometric regions:** boxes, polytopes, products, and ellipsoids may have direct exact optimization or integration algorithms.
- **Repeated queries:** caching and reuse of a sign-invariant CAD can amortize construction cost.

Specialization is safe only when the shortcut's applicability conditions are proved or otherwise soundly enforced. A failed heuristic candidate search is not a proof of nonexistence.

## Variable ordering and projection choices

CAD is cylindrical **relative to a variable order**. Changing that order can change the degrees and number of projection polynomials, and thus the number of cells dramatically. The tradeoff is constrained by logical semantics: quantified variables and free variables cannot be permuted arbitrarily without preserving the intended quantifier prefix and output variables.

Similarly, replacing full projection by reduced equational-constraint projection may eliminate many polynomials but introduces proof obligations about degree, delineability, and root order. Complete Lazard projection gives a more validated certified route for some exceptional fibers. Read the [CAD proof guide](../native_cad_proofs.md) for the concrete backends supported by `semialg`.

## Roadmaps versus CAD for topology

If the only question is whether two points lie in the same connected component, building a full sign-invariant CAD may do more work than is needed. Classical roadmap algorithms target a one-dimensional skeleton of the set; they are important precisely because connectivity need not require a full ambient decomposition. But **having a roadmap theorem** is not the same as **having an implementation of a complete roadmap algorithm**. `semialg` currently documents its exact connectivity and supported cell-chain paths separately from general roadmap construction.

For many real scientific problems, a good practical strategy is to first reduce degrees of freedom, eliminate simple equalities, exploit product or convex structure, and only then invoke a general exact solver. Numerical computation can suggest cases and candidate answers, but the exact decision must not depend solely on a floating-point tolerance.

## What a certified answer means

An exact calculation should distinguish the following outcomes:

1. **Established true:** a supported algorithm completes and certifies the requested property.
2. **Established false:** an exact counterexample or a complete decision proof refutes it.
3. **Unsupported or incomplete:** a backend is inapplicable, a verification condition is unresolved, or the resource budget is exhausted.

The third case is neither a positive nor a negative mathematical answer. A useful proof object must state the property it verifies: checking a stored polynomial identity is not automatically the same as independently verifying the geometry of a CAD, and truth-invariance for one formula is not sign-invariance for all polynomials.

`semialg` emphasizes explicit failure and provenance where the guaranteed computation cannot be completed. Consult [exactness and certification](../concepts/exactness_and_certification.md), [resource controls](../cad_qe.md), [certificate model](../concepts/certificate_model.md), and [limitations](../limitations.md).

## A practical problem-selection checklist

Before running a general CAD or QE task, determine:

1. Are all coefficients exact and real, and are all free parameters declared correctly?
2. Is the problem feasibility, projection, formula equivalence, optimization, connectivity, or an actual full decomposition?
3. Can affine substitution, low-degree elimination, exact univariate methods, or structural optimization answer it first?
4. Is a truth-invariant CAD sufficient, or will later queries require global sign information?
5. What happens at leading-coefficient degree drops, singular fibers, and excluded endpoints?
6. What exact certificate or counterexample is expected if the computation succeeds?
7. What are the memory, timeout, and cell budgets, and how will an inconclusive result be reported?

These questions often matter more than choosing the fastest projected polynomial routine for a poorly specified task.

**References:** Basu, Pollack, and Roy (2006); Collins (1975). See [bibliography](../references.md#mathematical-background-and-further-reading) and [performance guide](../guides/performance.md).

**Next:** [Exercises and worked answers](exercises.md) · [Applications](applications.md).
