# Cylindrical algebraic decomposition: a worked example

This tutorial constructs a complete two-variable CAD by hand, explains its exact logical guarantee, then relates the result to the implementation in `semialg`. Readers unfamiliar with CAD should begin with [real closed fields and quantifier elimination](real_geometry.md) and [Sturm sequences](algebraic_computation.md).

## What a CAD does

A **cylindrical algebraic decomposition** partitions \(\mathbb R^n\) into finitely many connected semialgebraic **cells** arranged recursively with respect to an ordered list of variables. Over each cell in the first \(k-1\) coordinates, the \(k\)-th coordinate is either

- a **section**: a point on a continuous algebraic root function; or
- a **sector**: a region strictly between consecutive root functions (including unbounded sectors).

A **sign-invariant** CAD gives every designated input polynomial constant sign \(+,0,-\) on each cell. A **truth-invariant** CAD need only make the chosen Boolean formula constant on each cell. The latter is often cheaper and is not interchangeable with global sign invariance when reusing a CAD to answer new questions.

The classical algorithm proceeds in two directions: **project** polynomial information to lower dimensions, then **lift** cells back up by isolating real roots over certified samples.

## Example: \(y^2=x^3-x\)

Study the real plane curve

\[ C=\{(x,y)\in\mathbb R^2:f(x,y)=0\},\qquad f(x,y)=y^2-x^3+x. \]

Order coordinates \(x<y\). The polynomial is quadratic in \(y\), with discriminant

\[ \operatorname{disc}_y(f)=4(x^3-x)=4x(x-1)(x+1). \]

Hence the potential changes in the number or coincidence of real \(y\)-roots occur at \(x=-1,0,1\). A projection to the \(x\)-axis must separate these points.

### Step 1: split the real line into base cells

The real roots of the projection polynomial yield seven ordered cells:

\[ (-\infty,-1),\ \{-1\},\ (-1,0),\ \{0\},\ (0,1),\ \{1\},\ (1,\infty). \]

We choose one sample point in each interval and use the exact defining point for each singleton. The important property is that the number and order of roots in \(y\) remain stable over an entire interval cell; this requires the projection theorem, not merely successful substitution at a single sample.

### Step 2: lift each base cell

Rewrite the equation as \(y^2=x^3-x\). Each interval can have zero or two distinct real roots, and each special \(x\)-point has a double root at \(y=0\).

| Base cell | Sign of \(x^3-x\) | Real \(y\)-roots | Cells in the vertical stack |
|---|---|---|---:|
| \(x<-1\) | Negative | None | 1 |
| \(x=-1\) | Zero | \(0\) (double) | 3 |
| \(-1<x<0\) | Positive | \(\pm\sqrt{x^3-x}\) | 5 |
| \(x=0\) | Zero | \(0\) (double) | 3 |
| \(0<x<1\) | Negative | None | 1 |
| \(x=1\) | Zero | \(0\) (double) | 3 |
| \(x>1\) | Positive | \(\pm\sqrt{x^3-x}\) | 5 |

For example, at \(x=-2\), \(x^3-x=-6\), so

\[ f(-2,y)=y^2+6>0\quad\text{for every real }y. \]

The entire vertical fiber belongs to one positive-sign sector. At \(x=-\tfrac12\), \(x^3-x=\tfrac38\), and the roots are \(\pm\sqrt{3/8}\). Lifting produces five cells: below the lower root, on the lower root, between roots, on the upper root, and above it. Their \(f\)-signs are respectively \(+\), \(0\), \(-\), \(0\), \(+\).

At \(x=-1,0,1\), \(f=y^2\): the three cells are \(y<0\), \(y=0\), \(y>0\), with signs \(+,0,+\).

The total number of two-dimensional CAD cells (of various dimensions) is

\[ 1+3+5+3+1+3+5=\boxed{21}. \]

### Step 3: interpret the cells geometrically

The zero set consists of two connected components:

1. A **bounded oval** over \(-1\le x\le0\), where the two \(y\)-branches meet at \((-1,0)\) and \((0,0)\).
2. An **unbounded component** over \(x\ge1\), where the branches meet at \((1,0)\).

Despite the tangency-looking double roots in the *vertical fibers*, the plane curve itself is nonsingular. Indeed

\[ \nabla f=(1-3x^2,\;2y). \]

A singular point would require \(y=0\), \(x^2=1/3\), and \(x^3-x=0\) simultaneously; there is no such real or complex point. Thus this is **not** a nodal cubic. A repeated root in one projection direction does not necessarily mean a singularity of the curve.

The sign-invariant CAD describes not only the curve but the positive and negative regions around it. For example, the entire strip \(0<x<1\) is a positive-sign cell for \(f\), whereas for \(-1<x<0\), the region between the two graph branches has negative sign.

## Step 4: quantifier elimination from the same idea

If the goal is only to project the real curve onto the \(x\)-axis, we ask

\[ \exists y\;\bigl[y^2=x^3-x\bigr]. \]

A real square exists exactly when its argument is nonnegative. Therefore the quantified formula is equivalent to

\[ (x\in[-1,0])\ \lor\ (x\in[1,\infty)), \]

or, algebraically, \(x^3-x\ge0\). The projected solution is a **union of intervals and boundary points**, all semialgebraic. This also shows why real projection is more precise than simply eliminating \(y\) with a complex algebraic resultant: the latter can forget the inequality determining whether real roots exist.

A second example is [projection of a disk above a line](real_geometry.md#worked-example-projection-of-a-disk-above-a-line), which has a different quantifier-free answer.

## Step 5: distinguish CAD refinements

A full sign-invariant CAD may perform more work than the question requires.

**Partial CAD.** For \(\exists y\,\Phi(x,y)\), a single true witness in a suitably certified stack can settle an entire base cell; for \(\forall y\,\Phi\), one false cell can settle it. A complete vertical stack is unnecessary when the formula's truth has already been established.

**Equational constraints (ECs).** If a formula requires \(f=0\), reduced projection can concentrate on sections where \(f\) vanishes, instead of maintaining every polynomial's sign throughout the ambient space. At successive projection levels, resultants can provide **implied** ECs. Correctness still requires the applicable degree, root-order, and nullification hypotheses.

**Truth-table invariant CAD (TTICAD).** Instead of making all polynomial signs invariant, one may preserve the truth values of a selected family of formulas. Different formulas can have different relevant ECs.

**Lazard projection and valuation-aware lifting.** Exceptional parameter fibers, including degree loss and positive-dimensional nullification, must not be ignored. The full Lazard route uses valuation information and has its own proof obligations. A complete CAD path may restart with a more validated operator when a reduced method cannot be certified.

These are different guarantees and strategies, not merely alternate names for the same decomposition. See [native CAD proofs](../native_cad_proofs.md), [CAD and QE](../cad_qe.md), and [algorithm selection](../concepts/algorithm_selection.md).

## A small reproducible algebraic check

The following uses SymPy only, so the mathematical projection computation can be checked without a semialg installation:

```python
import sympy as sp

x, y = sp.symbols("x y", real=True)
f = y**2 - x**3 + x
projection = sp.factor(sp.discriminant(f, y))
assert projection == 4 * x * (x - 1) * (x + 1)
for value, expected in [(-2, 0), (-sp.Rational(1, 2), 2), (sp.Rational(1, 2), 0), (2, 2)]:
    p = sp.Poly(f.subs(x, value), y)
    assert p.count_roots(-sp.oo, sp.oo) == expected
```

For **native semialg CAD construction and independent replay**, use the tested API in [native CAD proofs](../native_cad_proofs.md). The [CAD tutorial](../tutorials/cad_factorized_partial.md) illustrates reductions and partial lifting. These references, instead of the hand calculation, define which backends and coefficient fields are supported by the package.

## Exercises

1. Decompose the \(x\)-axis for \(g(x,y)=x^2+y^2-1\) projected in \(y\). How many stack cells occur over \(-1<x<1\) and \(x>1\)?
2. Eliminate \(y\) from \(\exists y\,[x^2+y^2=1]\) and include the endpoints.
3. Explain why \(\partial f/\partial y=0\) at a point does **not** suffice to show the plane curve \(f=0\) is singular.
4. For a formula \((f=0)\land(g>0)\), explain why a truth-invariant decomposition need not preserve \(g\)'s sign off \(f=0\).

**References:** Collins (1975); Arnon, Collins, and McCallum (1984); McCallum, Parusiński, and Păunescu (2019). See [algorithm references](../references.md#cylindrical-algebraic-decomposition).

**Next:** [Topology and roadmaps](topology.md) · [Complexity and scope](complexity.md).
