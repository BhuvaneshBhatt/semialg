# Polynomial optimization and positivity certificates

Optimizing a polynomial over polynomial constraints is a central problem in real algebraic geometry. A numerical optimizer may find an excellent point; an **exact** method must additionally prove that no better feasible point exists, or explain why the question falls outside its certified scope.

## Formulating optimization as a logical problem

For a semialgebraic set \(S\) and polynomial \(f\), the statement that \(m\) is a lower bound is

\[ \forall x\in S\quad f(x)\ge m, \]

which is a first-order formula over the real numbers when \(S\) is described by polynomial equations and inequalities. Conversely, a strict improvement exists precisely when

\[ \exists x\in S\quad f(x)<m. \]

Quantifier elimination can therefore certify the lower bound by proving that the improvement formula is false. This is mathematically general but often computationally expensive.

A different question is **attainment**: does there exist \(x\in S\) with \(f(x)=\inf_S f\)? For example, over \(S=(0,1)\), the function \(f(x)=x\) has infimum 0, but no minimum. Over \(S=[0,1]\), the same infimum is attained at \(x=0\).

The distinction persists in higher dimensions, especially for open or unbounded domains. If \(S\) is nonempty and compact, a continuous polynomial attains its minimum and maximum; without compactness, finite bounds need not be attained.

## Worked optimization example: an ellipsoid

Let \(A\) be a symmetric positive-definite matrix, \(c\) a fixed vector, and \(r>0\). Consider

\[ S=\{x:(x-c)^\mathsf TA(x-c)\le r^2\},\qquad \ell(x)=u^\mathsf Tx+\beta. \]

Using \(z=A^{1/2}(x-c)\), the constraint becomes \(\|z\|_2\le r\), and Cauchy–Schwarz gives

\[
|u^\mathsf T(x-c)|
\le r\sqrt{u^\mathsf T A^{-1}u}. \]

Thus the extrema are exactly

\[ \boxed{\min_S\ell=\beta+u^\mathsf Tc-r\sqrt{u^\mathsf TA^{-1}u}}, \]

\[ \boxed{\max_S\ell=\beta+u^\mathsf Tc+r\sqrt{u^\mathsf TA^{-1}u}}. \]

Both are attained on the closed ellipsoid (if \(u=0\), the objective is constant). This explicit formula provides a **certificate-producing structural shortcut** without enumerating a general CAD. It also explains why optimizing a linear objective on the corresponding **open** ellipsoid gives the same infimum and supremum but, when \(u\ne0\), neither endpoint is attained.

The [structural rules guide](../structural_rules.md) documents the implemented ellipsoid specialization, its assumptions, and its result objects. Do not infer support for arbitrary matrix inequality or nonlinear ellipsoidal expressions from this mathematical formula.

## Univariate extrema and critical points

For a polynomial \(p(x)\) over a closed bounded interval \([a,b]\), every global minimum occurs at an endpoint or at an interior critical point where \(p'(x)=0\). Exact real-root isolation and sign determination suffice to enumerate and compare these candidates.

For example, let \(p(x)=x^3-3x\) on \([-2,2]\). Its derivative is \(3x^2-3\), with critical points \(-1,1\). Evaluating gives

\[ p(-2)=-2,\quad p(-1)=2,\quad p(1)=-2,\quad p(2)=2. \]

The minimum is \(-2\), attained at \(x=-2\) and \(x=1\); the maximum is 2, attained at \(x=-1\) and \(x=2\). On punctured or open intervals, endpoint limits must be considered, without treating excluded points as optimizers.

For multivariate problems, Lagrange multipliers and Karush–Kuhn–Tucker systems produce possible candidates under suitable constraint qualifications. They do **not** by themselves prove global optimality. Degeneracies, singular feasible sets, and positive-dimensional critical loci need additional analysis; a final exact decision or range certificate can establish that no better point exists.

## Positivity, sums of squares, and real algebra

An **SOS polynomial** is a sum of squares of polynomials,

\[ p(x)=q_1(x)^2+\cdots+q_k(x)^2. \]

It is automatically nonnegative. The converse fails: there exist nonnegative polynomials that are not polynomial sums of squares. The classical **Motzkin polynomial** is

\[ M(x,y)=x^4y^2+x^2y^4-3x^2y^2+1. \]

It is nonnegative by the arithmetic–geometric mean inequality applied to \(x^4y^2\), \(x^2y^4\), and \(1\), yet it is not an SOS of polynomials. This prevents one from treating failure to find an SOS representation as proof that the polynomial becomes negative somewhere.

Artin's solution of Hilbert's 17th problem says a globally nonnegative real polynomial is a sum of squares of **rational functions**. More general **Positivstellensätze** characterize semialgebraic feasibility or positivity by algebraic identities involving the defining constraints and sums of squares. Some useful SOS certificates require assumptions such as an Archimedean quadratic module or strict positivity; a bounded-degree SOS search need not succeed even when the inequality is true.

## Gram matrices and semidefinite programming

Let \(z(x)\) be a vector of monomials of a chosen degree. A polynomial identity

\[ p(x)=z(x)^\mathsf T Q z(x),\qquad Q\succeq0 \]

is an SOS certificate when \(Q\) is real symmetric positive semidefinite. The unknown coefficients of \(Q\) are constrained linearly by matching polynomial coefficients, while \(Q\succeq0\) is a semidefinite constraint. This makes SOS search amenable to semidefinite programming (SDP).

An SDP solver generally returns approximate floating-point values. For an **exact mathematical certificate**, one must reconstruct suitable exact coefficients and independently verify the polynomial identity and positive semidefiniteness. A small numerical residual or a nearly positive matrix is not enough to establish an exact polynomial inequality.

The package includes **selected** SOS and Gram-certificate functionality; it does not claim a general arbitrary-degree Positivstellensatz decision engine or the whole Lasserre moment/SOS hierarchy. See the [SOS tutorial](../tutorials/sos_certification.md), [algorithm references](../references.md#sums-of-squares-and-gram-matrices), and [limitations](../limitations.md).

## Comparing three strategies

| Strategy | Strength | Limitation |
|---|---|---|
| Structural exact rule | Very fast for a recognized geometry or objective | Requires precise hypotheses |
| Algebraic critical points / KKT | Generates exact candidates in many constrained problems | Candidates need a completeness/globality proof |
| SOS / SDP certificate | Can provide concise proof identities for inequalities | May need nontrivial degree and assumption choices; numerical output needs exact replay |
| CAD / quantifier elimination | General logical route for supported exact real polynomial formulas | Can be prohibitively expensive |

A good solver uses the cheapest **certified** route and falls back validatedly when a shortcut is inapplicable. A resource-limit failure is not an inequality counterexample.

## Exercises

1. Prove that the infimum of \(x^2\) on \((-1,0)\) is 0 but is not attained. How does the answer change on \([-1,0]\)?
2. Apply the ellipsoid formula to \(x^2+4y^2\le1\) and objective \(x+y\). Determine both extrema.
3. Use AM–GM to prove \(M(x,y)\ge0\) for the Motzkin polynomial.
4. Explain why a KKT point need not be a global minimum, and why an approximate PSD matrix is not automatically an exact SOS certificate.

**References:** Basu, Pollack, and Roy (2006); Parrilo (2003). See [bibliography](../references.md#mathematical-background-and-further-reading).

**Next:** [Complexity and scope](complexity.md) · [Applications](applications.md).
