# Package ownership

Use the package that owns the mathematical question. These packages share exact arithmetic and keep their dependency direction acyclic.

| Question | Canonical package | Main interfaces |
| --- | --- | --- |
| Finite quotient arithmetic, exact coefficient fields, zero-dimensional solving | [algroots](https://github.com/BhuvaneshBhatt/algroots) | `quotient`, `domains`, `finite_fields`, `finite_decomposition`, root solvers |
| Reusable coefficient-field and plane-curve arithmetic | [algroots](https://github.com/BhuvaneshBhatt/algroots) | `subfields`, `function_fields`, `curves` |
| General coordinate rings, schemes, normalization and function decomposition | [alggeom](https://github.com/BhuvaneshBhatt/alggeom) | `ideals`, `rings`, `schemes`, `normalization`, `decomposition` |
| Real loci, inequalities, parameter conditions and exact real graph verification | [semialg](https://github.com/BhuvaneshBhatt/semialg) | `decision`, `qe`, `conditional`, `function_graph`, `function_verification` |

algroots has no dependency on alggeom or semialg. semialg depends on algroots; it interprets exact algebra through real signs, inequalities and witnesses. alggeom depends on algroots and uses semialg through its optional real extra. Cross-package applications belong in alggeom's real profile. The backends do not import the application package.

A positive-dimensional coordinate ring has a Krull dimension and generally infinite dimension over its coefficient field. `algroots.QuotientAlgebra.dimension` is a finite vector-space dimension. A real branch condition belongs to a semialgebraic graph, not to an algebraic ideal or its Zariski closure. Finite fields do not have a real ordering.

Function-decomposition search, coefficient reconstruction, chains and their branch-aware adapters are supplied by `alggeom.decomposition`. The backend field arithmetic is reusable independently of a decomposition search. Generic real graph existence, uniqueness and composition checks remain available through `semialg.function_verification`.

Rational-coordinate membership and coefficient rebasing have one implementation in `algroots.rational_coordinates`. semialg's GTZ certificates use algroots for finite quotient arithmetic while retaining their proof records and replay. Generic real-graph comparison lives in `semialg.function_verification`; alggeom supplies decomposition-specific graphs and guards.

## Function analysis and optimization consumers

| Question | Owner | Relationship to semialg |
| --- | --- | --- |
| Which outputs does a supported analytic function attain? | [funcprops](https://github.com/BhuvaneshBhatt/funcprops/blob/main/docs/range-analysis.md) | Reuses semialg for exact algebraic images and real side conditions |
| What is an exact semialgebraic infimum or optimizer? | semialg | Native exact real optimization, including algebraic parameter conditions |
| What is an infimum, optimizer or optimization certificate in the general solver portfolio? | [symbopt](https://github.com/BhuvaneshBhatt/symbopt/blob/main/docs/guide/ranges-and-extrema.md) | Uses optional exact real projection, optimization and inequality decisions |
| What are a limit, growth relation or remainder bound? | [asymptotic](https://github.com/BhuvaneshBhatt/asymptotic) | The local package uses semialg and funcprops; a future image-tail adapter must preserve their dependency direction |

semialg's `function_range` exact modes own the supported real algebraic image formula; its explicit bounds mode owns an interval enclosure. A consumer must retain that distinction. funcprops' critical-root and value facts can inform optimization without being optimization certificates. The semialg backend does not import funcprops or symbopt. See [range semantics](range-images.md).
