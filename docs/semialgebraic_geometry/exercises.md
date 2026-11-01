# Exercises, checks, and further reading

These exercises accompany the mathematical overview and do not assume undocumented `semialg` APIs. Several have short answers below. For executable package workflows, see the [examples catalog](../examples/index.md) and [demo notebook](https://github.com/BhuvaneshBhatt/semialg/blob/main/notebooks/semialg_demo.ipynb) in the repository (the notebook is not necessarily included in the documentation-only ZIP).

## I. Algebra and algebraic sets

**1. The real-versus-complex zero set.** Describe the real and complex zero sets of \(x^2+y^2=0\). Explain why a complex elimination method alone cannot report the dimension of the real solution set.

**2. Ideal membership.** Let \(I=(xy-1,x-y)\). Show directly that \(y^2-1\in I\), and use that result to list the real solutions. [Hint: multiply the second generator by \(y\).]

**3. Real radical.** Verify that the real zero set of \(x^2+y^2\) is only the origin. Explain why \(x\) and \(y\) belong to the real radical of \((x^2+y^2)\) but not to its ordinary radical.

## II. Projection and quantifiers

**4. Square roots and quantifier elimination.** Eliminate \(y\) from \(\exists y\,[y^2=x\land y\le2]\). Pay attention to negative and positive values of \(x\).

**5. Projected disk.** Prove that the projection onto the \(x\)-axis of \(\{(x,y):x^2+y^2\le1,\;y\ge x\}\) is \([-1,1/\sqrt2]\). Explain why squaring \(x\le\sqrt{1-x^2}\) without first checking the sign of \(x\) can lead to mistakes.

**6. Absolute value.** Express the graph of \(y=|x|\) using only polynomial equations, inequalities and Boolean operations. Is its image under the \(x\)-projection semialgebraic?

## III. Roots and CAD

**7. A Sturm count.** Construct the Sturm sequence of \(x^3-x\) and count roots in \((-2,2)\) and \((-1/2,1/2)\). What happens if an interval endpoint is a root?

**8. Circle CAD.** For \(f(x,y)=x^2+y^2-1\), project in \(y\), partition the \(x\)-axis, and count the vertical CAD cells over \(x<-1\), \(x=-1\), \(-1<x<1\), \(x=1\), and \(x>1\). Compare with the [worked cubic CAD](cad_worked_example.md).

**9. Implied EC.** Show that \(z-y=0\) and \(z+y=0\) imply \(y=0\). Calculate the resultant with respect to \(z\). Why is such an implication not, on its own, a certificate that every reduced projection step is valid?

## IV. Topology and optimization

**10. Boundaries.** Compute \(S^\circ\), \(\overline S\) and \(\partial S\) for \(S=(0,1)\). Explain the error in the formula \(\partial S=S\setminus S^\circ\).

**11. Euler characteristic.** Compare the compactly supported Euler characteristic of an open interval with that of a closed interval. Why does compactness matter when comparing it to the ordinary homotopy Euler characteristic?

**12. An ellipsoid.** Minimize and maximize \(x+y\) over \(x^2+4y^2\le1\), then explain what changes when the inequality is strict. [Hint: \(A=\operatorname{diag}(1,4)\).]

**13. SOS versus nonnegativity.** Prove \(x^4y^2+x^2y^4+1\ge3x^2y^2\) by AM–GM. What does that establish about the Motzkin polynomial, and what does it *not* establish about polynomial SOS representations?

**14. Model scope.** A program returns no answer after an exact CAD computation hits a cell budget. Explain why that is logically different from returning false. What metadata would help a caller decide what to try next?

## Selected answers and checks

- **4.** There is a nonpositive solution \(y=-\sqrt x\le2\) for every \(x\ge0\), so the result is simply \(x\ge0\). The upper bound \(y\le2\) does not restrict \(x\), because the negative branch remains available.
- **8.** Projection uses \(1-x^2\) (up to a nonzero factor). The stacks have \(1,3,5,3,1\) cells, respectively, for a total of 13 CAD cells.
- **9.** The resultant is \(2y\) up to an irrelevant nonzero constant/sign. A valid reduced CAD also needs hypotheses ensuring correct root and sign behavior during lifting.
- **10.** \(S^\circ=(0,1)\), \(\overline S=[0,1]\), \(\partial S=\{0,1\}\). The incorrect formula gives the empty set because it excludes boundary points not in \(S\).
- **11.** \(\chi_c((0,1))=-1\), while \(\chi_c([0,1])=+1\). The usual homotopy Euler characteristic of each interval is \(+1\), so the open case demonstrates the distinction.
- **12.** \(\max=\sqrt5/2\), \(\min=-\sqrt5/2\); both are attained on the boundary when the ellipsoid is closed. For the open ellipsoid, they become an unattained supremum and infimum.
- **14.** A budget failure leaves truth unknown. Report the budget, backend, relevant proof obligation, and any validated intermediate certificates; consider specialization, variable order, or increased resources rather than reporting a mathematical counterexample.

## Suggested reading path

If you are new to the subject, read [foundations](foundations.md), [real closed fields](real_geometry.md), [algebraic algorithms](algebraic_computation.md), and the [worked CAD](cad_worked_example.md) in that order. Readers mainly interested in geometric applications can proceed to [topology](topology.md) and [optimization](optimization_positivity.md) after the introductory real-geometry chapter.

The [bibliography](../references.md#mathematical-background-and-further-reading) provides verified book references; the existing [algorithm references](../references.md) provide the papers behind individual `semialg` backends.
