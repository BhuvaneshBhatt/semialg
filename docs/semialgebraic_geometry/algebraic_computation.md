# Algebraic algorithms: elimination, roots, and signs

This chapter connects polynomial algebra to the primitives used in real quantifier elimination. The mathematical examples can be followed without installing `semialg`; API links point to the implementations and their actual supported scope.

## Why univariate root counting is central

In one variable, a nonzero polynomial has finitely many roots. Exact algorithms can isolate each distinct real root in an interval with rational endpoints and decide the sign of another polynomial at that root. In multiple variables, CAD repeatedly reduces questions to precisely these operations over algebraic sample points.

Numerical approximations are useful for proposing isolating intervals, but close roots, multiplicities, and equality at algebraic points demand an exact certificate. See [certified algebraic roots](../guides/certified_algebraic_roots.md).

## Sturm sequences

Let \(f\in\mathbb Q[x]\) be squarefree. A **Sturm sequence** begins with \(f_0=f\), \(f_1=f'\) and continues by negative Euclidean remainders:

\[ f_{i+1}=-\operatorname{rem}(f_{i-1},f_i). \]

At a real number where \(f\ne0\), define \(V(t)\) as the number of sign changes in the evaluated sequence, ignoring zero entries. If neither endpoint is a root, **Sturm's theorem** gives the number of distinct roots in \((a,b)\) as \(V(a)-V(b)\). Endpoint roots require an explicit convention or one-sided evaluation; simply treating a closed interval as an open interval loses boundary roots.

### Worked example: three real roots of a cubic

Take \(f=x^3-x\). Its Sturm sequence can be chosen as

\[ f_0=x^3-x,\quad f_1=3x^2-1,\quad f_2=\tfrac23x,\quad f_3=1. \]

At \(x=-2\), the signs are \((- ,+,-,+)\), giving three changes. At \(x=2\), all four signs are positive, giving zero changes. Hence the polynomial has \(3-0=3\) distinct real roots in \((-2,2)\); they are \(-1,0,1\).

We can also count roots in an interval **without solving for them**. For example, at \(-\tfrac12\) the signs are \((+,-,-,+)\), giving two changes. At \(\tfrac12\) they are \((-, -,+,+)\), giving one. Therefore the interval \((-\tfrac12,\tfrac12)\) contains exactly one root.

SymPy can check the counts for this pure-mathematics example:

```python
import sympy as sp

x = sp.symbols("x")
f = x**3 - x
assert sp.Poly(f, x).count_roots(-2, 2) == 3
assert sp.Poly(f, x).count_roots(-sp.Rational(1, 2), sp.Rational(1, 2)) == 1
```

A squarefree part can be obtained by dividing a polynomial by its gcd with the derivative; counting **distinct** roots must not accidentally count repeated roots multiple times.

## Resultants, discriminants, and degree drops

The **resultant** \(\operatorname{Res}_y(f,g)\) is a polynomial in the remaining variables that vanishes when specialized polynomials in \(y\) have a common root, provided their degree and leading-coefficient behavior is interpreted correctly. A discriminant detects multiple-root events.

For example, for \(f(x,y)=y^2-x\), the discriminant in \(y\) is \(4x\). The real fibers have

- no real root if \(x<0\),
- one repeated root at \(x=0\), and
- two distinct roots if \(x>0\).

The simple parameter partition \(x<0\), \(x=0\), \(x>0\) is a small instance of a CAD projection-and-lifting pattern.

A caution about **degree drops**: the specialization of a resultant may vanish identically because a leading coefficient disappears, not solely because two ordinary finite roots collide. Safe CAD projection retains the coefficients, subresultant information, or valuation-aware data needed for such exceptional fibers.

## Deriving an equational constraint by elimination

Suppose a formula requires both \(z-y=0\) and \(z+y=0\). Their resultant with respect to \(z\) is \(2y\), so their simultaneous vanishing implies \(y=0\). That equation can become an **implied equational constraint** at the next projection level. However, merely knowing an equation holds on a solution locus does not justify dropping every other projection polynomial: reduced CAD must establish the appropriate root-order and degree conditions. See [native CAD proofs](../native_cad_proofs.md).

## Gröbner bases and rational univariate representations

A Gröbner basis is especially useful for solving **zero-dimensional** polynomial systems: after suitable elimination, one can often reduce the solutions to roots of a single polynomial. A **rational univariate representation** (RUR) encodes coordinates as rational functions of an algebraic parameter \(t\):

\[ q(t)=0,\qquad x_i=\frac{a_i(t)}{b_i(t)}. \]

To obtain real solutions, one isolates the real roots of \(q\), checks denominator nonvanishing, and evaluates inequality signs **exactly**. A complex polynomial solution is not automatically a real feasible point. Positive-dimensional ideals require different techniques or decomposition; there need not be a finite list of isolated solutions.

## Algebraic numbers and Thom information

An exact real algebraic number can be represented by a defining polynomial together with data identifying one real root, such as an isolating interval or appropriate derivative signs (a **Thom encoding**). The crucial operations are comparisons, ordering, equality, and sign determination. These support sample points and root correspondence between adjacent CAD cells.

For a polynomial fiber that changes as parameters vary, it is not sufficient to choose the numerically nearest root after a collision: one must preserve mathematically certified branch identity where available, or report the ambiguity. The [native proof guide](../native_cad_proofs.md) describes the supported branch and lifting certificates.

## Division of responsibilities in the project

`semialg` uses exact arithmetic and selected algebraic routines as part of its real-geometry workflows. Some general algebraic geometry and zero-dimensional machinery lives in related packages. The [package ownership guide](../package_ownership.md) and [algebraic API reference](../reference/algebraic.md) identify those boundaries. Statements about Sturm's theorem, RUR, or Gröbner bases here are mathematical explanations, not promises of a particular high-level `semialg` function name.

**References:** Cox, Little, and O'Shea (2015); Basu, Pollack, and Roy (2006); Collins (1975). See [bibliography](../references.md#mathematical-background-and-further-reading) and [algorithm references](../references.md#cylindrical-algebraic-decomposition).

**Next:** [Worked cylindrical algebraic decomposition](cad_worked_example.md).
