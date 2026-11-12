# Known issues and implementation boundaries

This page distinguishes reproducible defects from unsupported algorithms. The [limitations](limitations.md) describe the broader mathematical scope.

## Fixed defects with permanent regressions

| Defect | Correct behavior | Regression |
| --- | --- | --- |
| Native-fiber multiplicity was recorded as one after specialisation | `y**2-x**2+2` at `x=sqrt(2)` is `y**2`, with one distinct root of multiplicity two; nested towers and higher powers use exact derivatives | [Native multiplicity](https://github.com/BhuvaneshBhatt/semialg/blob/main/tests/test_native_multiplicity.py) |
| A Collins-to-Lazard restart reset the cell budget | Work performed by both attempts counts toward the cooperative budget | [Budget restarts](https://github.com/BhuvaneshBhatt/semialg/blob/main/tests/test_cad_budget_restarts.py) |
| Lazard could return or cache a result after its deadline elapsed during final verification | A final budget check rejects it before caching | [Budget deadlines](https://github.com/BhuvaneshBhatt/semialg/blob/main/tests/test_cad_budget_restarts.py) |
| Native CAD objects and public CAD functions containing read-only mappings could not be serialized | Constructor fields survive pickle, mappings remain immutable, and proof replay works after caches are cleared | [Proof lifecycle](https://github.com/BhuvaneshBhatt/semialg/blob/main/tests/test_native_proof_lifecycle.py) |
| Preimages lost map domains when target predicates simplified away | Rational poles and real-root domains are retained for every output coordinate, including unused coordinates and Piecewise branches | [Preimage domains](https://github.com/BhuvaneshBhatt/semialg/blob/main/tests/test_preimage_domains.py) |
| Witness adapters treated absent samples or incomplete backends as infeasibility | Unknown searches, incomplete CAD and non-finite un-enumerated sets stay unknown; empty finite solution sets remain certified unsatisfiable | [Instance status evidence](https://github.com/BhuvaneshBhatt/semialg/blob/main/tests/test_instance_status_evidence.py) |
| SymPy tuple solutions were unpacked as scalar values | Complex multivariate witnesses retain every coordinate and are not discarded as solving failures | [Complex tuple regression](https://github.com/BhuvaneshBhatt/semialg/blob/main/tests/test_instance_status_evidence.py) |
| Concave polygons were replaced by convex hulls under embeddings/conversions | Exact lazy affine images retain concavities; convex-only conversion rejects concave sources | [Polygon images](https://github.com/BhuvaneshBhatt/semialg/blob/main/tests/test_affine_polygon_images.py) |
| Geometric identity and shortcuts dropped construction conditions | Simplex/polygon identity includes validity conditions, unresolved shortcuts defer, and affine images retain the source conditions | [Condition identity](https://github.com/BhuvaneshBhatt/semialg/blob/main/tests/test_region_condition_identity.py) |
| Critical images lost components through output-symbol collisions or unavailable computations | Default output coordinates avoid collisions; unavailable locus/image steps raise instead of returning an incomplete union as an exact set | [Critical image completeness](https://github.com/BhuvaneshBhatt/semialg/blob/main/tests/test_critical_image_completeness.py) |

Exact existential affine-block substitution also precedes symbolic-region virtual substitution, avoiding the observed barycentric/embedding expression blow-up. Universal blocks and vanishing parameter coefficients retain their existing exact elimnation paths.

## Unsupported or partial functionality

The high-level `cad` and `quantifier_eliminate` interfaces share cooperative cell/deadline budgets across nested constructions. Structured non-strict CAD retains an unknown outcome on exhaustion; direct outputs and QE raise. The preprocessing auxiliary-variable limit is also enforced. These checks cannot interrupt a single polynomial arithmetic operation. Use process isolation for a hard wall timeout.

Automatic QE uses partial lifting after specialist rules, including free parameters and alternating blocks. Every free cell retains truth information. Projection is built up front, so its growth remains a cost even when lifting prunes most branches.

Branch correspondence handles the supported collision and complex-birth cases, including the retained torus regressions. One-parameter degree loss includes finite and signed infinite limits. Ambiguous higher-dimensional continuation and multibase degree loss can raise `ExactComputationUnavailable`. See the [failure-aware closure tutorial](tutorials/failure_aware_closure.md).

Numerical meshing does not provide an ambient-isotopy certificate. Its conformity and provenance checks concern a numerical representation, not an exact topological equivalence.

Lower-priority findings are retained in [Minor issues](minor_issues.md). The [capability priorities](capability_priorities.md) distinguish tractable extensions from major algorithmic projects.
