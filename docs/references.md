# Algorithm references

`semialg` combines several classical algorithms with conservative certification and fallback layers.  The citations below identify the mathematical provenance of the named backends.  A name such as *ARS-style* or *McCallum-style* means that the implementation uses the cited construction but may add `semialg`-specific exact verification, specialization, or fallback steps; it should not be read as a line-for-line reimplementation.

## Real solving and polar varieties

- **[ARS2002]** P. Aubry, F. Rouillier, and M. Safey El Din, “Real Solving for Positive Dimensional Systems,” *Journal of Symbolic Computation* 34(6), 543–560 (2002). DOI: https://doi.org/10.1006/jsco.2002.0563.  This is the primary reference for the positive-dimensional critical-point/polar reduction used by `real_algebraic_feasibility`.

## Polynomial semidefiniteness and global minima

- **[ZZ2004]** G. Zeng and X. Zeng, “An Effective Decision Method for Semidefinite Polynomials,” *Journal of Symbolic Computation* 37(1), 83–99 (2004). DOI: https://doi.org/10.1016/S0747-7171(03)00073-7.  This is the canonical source motivating the package's `zeng_*` naming.
- **[ZX2012]** G. Zeng and S. Xiao, “Global Minimization of Multivariate Polynomials Using Nonstandard Methods,” *Journal of Global Optimization* 53, 391–415 (2012). DOI: https://doi.org/10.1007/s10898-011-9718-x.  The current coercivity/attained-global-minimum/critical-value implementation is closer in spirit to this global-minimization line than to a literal implementation of [ZZ2004].

Accordingly, `zeng_negative_point` should be understood as a **Zeng-family inspired certified critical-value backend**, not as a verbatim implementation of the 2004 algorithm.  `semialg` additionally uses exact leading-form coercivity proofs and delegates positive-dimensional critical loci to its ARS-style solver.

## Regular chains and equidimensional decomposition

- **[Kalkbrener1993]** M. Kalkbrener, “A Generalized Euclidean Algorithm for Computing Triangular Representations of Algebraic Varieties,” *Journal of Symbolic Computation* 15(2), 143–167 (1993). DOI: https://doi.org/10.1006/jsco.1993.1011.  This is a primary source for regular-chain/triangular unmixed-dimensional decomposition.
- **[Lazard1991]** D. Lazard, “A New Method for Solving Algebraic Systems of Positive Dimension,” *Discrete Applied Mathematics* 33, 147–160 (1991).  Squarefree triangular-set and saturation ideas used by the regular-chain literature originate in this line of work.
- **[ChenEtAl2013]** C. Chen, J. H. Davenport, J. P. May, M. Moreno Maza, B. Xia, and R. Xiao, “Triangular Decomposition of Semi-algebraic Systems,” *Journal of Symbolic Computation* 49, 3–26 (2013). DOI: https://doi.org/10.1016/j.jsc.2011.12.014.  This is a useful modern reference connecting regular chains to real/semi-algebraic solving.

The current `equidimensional_decomposition` uses exact factor/ideal splitting plus squarefree-regular-chain certification.  The separate `primary_decomposition` API now certifies prime, principal, monomial, and zero-dimensional primary decompositions, but does **not** yet claim a full arbitrary positive-dimensional GTZ decomposition.

## Algebraic function fields and norm factorization

- **[Trager1976]** B. M. Trager, *Algorithms for Manipulating Algebraic Functions*, Master's thesis, MIT (1976).  The monogenic function-field factorizer in `semialg.algebraic_function_fields` uses the norm/resultant descent associated with Trager: shift by a multiple of the extension generator, compute the norm to the lower field, factor there, and recover factors by exact gcd upstairs.  The implementation characterizes exceptional shifts exactly by vanishing of the norm polynomial's discriminant in the main variable; squarefreeness is tested equivalently by an exact polynomial gcd, so no finite shift budget enters correctness or completeness. Arbitrary input polynomials are first decomposed into squarefree parts with exact multiplicities over the tower using characteristic-zero gcd/Yun decomposition, so the squarefree condition is internal to the Trager core rather than a public restriction.

`semialg` applies this construction recursively over towers whose base is a rational-function field.  If norm descent is inconclusive, triangular decomposition retains its quotient-ring/saturation refinement as a one-sided exact fallback: a discovered factor identity may split a branch, but failure to discover one is never interpreted as irreducibility.

## Sums of squares and Gram matrices

- **[Parrilo2000]** P. A. Parrilo, *Structured Semidefinite Programs and Semialgebraic Geometry Methods in Robustness and Optimization*, PhD thesis, California Institute of Technology (2000). https://thesis.caltech.edu/1647/
- **[Parrilo2003]** P. A. Parrilo, “Semidefinite Programming Relaxations for Semialgebraic Problems,” *computer algebra systeml Programming* 96, 293–320 (2003). DOI: https://doi.org/10.1007/s10107-003-0387-5.

These references underlie the Gram formulation used by `SOSCertificate`: for a suitable monomial vector `z`, an SOS polynomial admits `p = z.T*Q*z` with `Q` positive semidefinite.  `semialg` treats numerical SDP output only as a search hint and independently verifies the exact identity and exact PSD condition.

## Cylindrical algebraic decomposition

- **[Collins1975]** G. E. Collins, “Quantifier Elimination for Real Closed Fields by Cylindrical Algebraic Decomposition,” *LNCS* 33, 134–183 (1975). DOI: https://doi.org/10.1007/3-540-07407-4_17.
- **[ACM1984]** D. S. Arnon, G. E. Collins, and S. McCallum, “Cylindrical Algebraic Decomposition I: The Basic Algorithm,” *SIAM Journal on Computing* 13(4), 865–877 (1984). DOI: https://doi.org/10.1137/0213054.
- **[McCallum1988]** S. McCallum, “An Improved Projection Operation for Cylindrical Algebraic Decomposition of Three-dimensional Space,” *Journal of Symbolic Computation* 5, 141–161 (1988). DOI: https://doi.org/10.1016/S0747-7171(88)80010-5.
- **[Brown2001]** C. W. Brown, “Improved Projection for Cylindrical Algebraic Decomposition,” *Journal of Symbolic Computation* 32(5), 447–465 (2001). DOI: https://doi.org/10.1006/jsco.2001.0463.
- **[Lazard1994]** D. Lazard, “An Improved Projection for Cylindrical Algebraic Decomposition,” in *Algebraic Geometry and its Applications*, 467–476 (1994). DOI: https://doi.org/10.1007/978-1-4612-2628-4_29.
- **[BDEMW2016]** R. Bradford, J. H. Davenport, M. England, S. McCallum, and D. Wilson, “Truth Table Invariant Cylindrical Algebraic Decomposition,” *Journal of Symbolic Computation* 76, 1–35 (2016). DOI: https://doi.org/10.1016/j.jsc.2015.11.002.

### What `semialg` actually implements

The **complete fallback** is Collins-style projection/lifting.  At each projection level `semialg` squarefree-normalizes the active family and carries lower-variable constraints downward; for active polynomials it projects content, all coefficients, discriminants, and pairwise resultants.  The resulting tower is lifted with exact real-root isolation and sign/truth evaluation.  This is the backend whose metadata reports `projection="collins"` / `collins-complete`.

The McCallum, Lazard, and TTICAD paths are **reduced, equational-constraint-shaped accelerators**, not independent claims of a complete literal implementation of every cited projection theory.  With a designated equational constraint they retain its leading coefficient, discriminant/subresultant-discriminant information, and resultants against the other active polynomials.  McCallum-style lifting checks nullification/well-orientedness; Lazard-style lifting additionally uses valuation-aware specialization.  Formula-oriented McCallum/Lazard results with a formula-wide equational constraint are accepted directly when the reduced-projection theorem's nullification/well-orientedness conditions are certified; a complete Collins CAD is then not constructed. Other reduced results use refinement certification or fall back to complete Collins CAD when required.  TTICAD uses the family-aware safe driver inspired by [BDEMW2016].

## Primitive elements and certified modular factorization

The algebraic-function-field layer also supports exact primitive-element
compression for finite characteristic-zero monogenic towers over
`QQ(parameters)`.  Adjacent extensions are compressed by choosing
`theta = alpha + c*beta`.  Rather than relying on a bounded or theorem-only
parameter search, semialg constructs the exact power-basis determinant guard
`B(c)`: its columns are the coordinates of `1, theta, ..., theta^(d-1)` in the
old product basis.  The binomial theorem assembles this matrix without adjoining
`c` to the algebraic tower.  Thus `B(c) != 0` is an explicit exact certificate
that `theta` generates the full compositum.  A deterministic integer outside
the finite zero set is chosen, finite-dimensional linear algebra reconstructs
the minimal polynomial and old generators, and all old defining relations are
re-evaluated exactly.  Compression certificates store the guard, chosen
parameter, generator images, and primitive expression and can be independently
replayed.  Bidirectional transport lets selected deep-tower factorization paths
work in the compressed field and then map factors back for exact reconstruction
in the source tower; a conservative degree/depth cost heuristic prevents
automatic compression when it is unlikely to reduce arithmetic growth.

At the rational leaf, parameter-free univariate factorization uses SymPy's
modular integer factorization machinery after exact denominator/content
normalization.  The reconstructed rational factors are multiplied back over
`QQ`; only exact equality with the input polynomial marks the modular result
certified.  Modular output therefore accelerates search but never substitutes
for exact reconstruction.

## Modular Groebner bases

- **[Arnold2003]** E. A. Arnold, “Modular algorithms for computing Groebner
  bases,” *Journal of Symbolic Computation* 35 (2003), 403–419.  `semialg` uses
  the standard multimodular pattern of lucky-prime images, Chinese-remainder
  accumulation, and rational reconstruction, but adds an explicit exact proof
  boundary: reconstructed bases are accepted only after exact Groebner and
  two-way ideal-membership verification over `QQ`.


## Modular resultants and subresultants

- **[BrownTraub1971]** W. S. Brown and J. F. Traub, “On Euclid's Algorithm and
  the Theory of Subresultants,” *Journal of the ACM* 18(4), 505–514 (1971).

`semialg` uses modular images and CRT/rational reconstruction as an
acceleration layer for resultants and subresultant sequences.  The resultant
certificate is independent of the modular proposal: exact fixed-degree
Sylvester determinants on a degree-determined interpolation grid prove the
reconstructed parameter polynomial.  Subresultant candidates additionally
undergo exact characteristic-zero normalization replay, because a generic PRS
identity alone would leave scale ambiguous.

## Primary decomposition and associated primes

- **[GTZ1988]** P. Gianni, B. Trager, and G. Zacharias, “Gröbner Bases and
  Primary Decomposition of Polynomial Ideals,” *Journal of Symbolic
  Computation* 6 (1988), 149–167.

`semialg` implements the certified localization and contraction infrastructure used by the full
positive-dimensional algorithm.  `certified_independent_localization` chooses
a maximal independent variable set from the exact leading monomial ideal and
certifies that extension to `QQ(U)[X]` is zero-dimensional.
`contract_localized_ideal` contracts a localized ideal by first computing its
monic Groebner basis over `QQ(U)`, clearing coefficient denominators, and
saturating by the clearing multiplier.  `saturation_stabilization` computes
the ascending colon chain `I:h^m` with no finite search bound until exact
Groebner equality proves stabilization; it independently verifies both the
Rabinowitsch saturation and the GTZ identity
`I = (I:h^infinity) intersect (I + <h^m>)`.

The zero-dimensional coefficient-field layer provides `zero_dimensional_primary_decomposition` over `QQ` and finite
monogenic algebraic extensions of `QQ`.  It certifies each coefficient-field
modulus as irreducible, lifts the zero-dimensional ideal to an ordinary
polynomial ideal over `QQ`, isolates every maximal local component by exact
separator saturation, and transports the components back to the coefficient
field.  Replay checks the field tower, minimal primes, saturation steps, exact
intersection, and transport identities independently.

The recursive GTZ driver over `QQ` uses the localized finite
`QQ(U)`-algebra is handled by exact multiplication matrices and the trace form.
In characteristic zero the kernel of the trace pairing is the nilradical; the
trace rank is therefore the dimension of the reduced algebra.  A deterministic
linear-form search stops only when the squarefree characteristic-polynomial
degree equals that rank, certifying a primitive element of the reduced quotient.
Irreducible characteristic factors separate the localized maximal ideals, exact
saturation isolates their primary factors, and both primary ideals and radicals
are contracted back to the polynomial ring.  The lower-dimensional or
same-dimensional residual branch is obtained from a certified parameter
annihilator and the GTZ saturation identity, then decomposed recursively.

The recursive proof object stores this decomposition as a tree.  Independent replay
verifies localizations, Artinian trace data, contractions, saturation exponents,
recursive source/companion identities, dimensions and degrees, same-associated-
prime merging, irredundancy, and final exact ideal reconstruction without
rerunning the search heuristic.
