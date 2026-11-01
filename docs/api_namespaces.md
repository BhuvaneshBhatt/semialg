# API namespaces

The package root contains broad mathematical operations and common region/function APIs. Low-level algorithms, certificates, planners, diagnostics, result records, and specialized workflows live in their owning namespaces.

Specialist APIs include Thom encodings (`semialg.algebraic`), roadmap construction (`semialg.roadmaps`), Hardt trivialization, dimension/component stratification and triangulation/Betti workflows (`semialg.topology.semialgebraic`), parameter-space covers (`semialg.parametric_geometry`), root-count conditions (`semialg.parameters`), and parametric map degree (`semialg.map_degree`).

Broad convenience operations such as connected-component counting, topology summaries, solvability conditions, and common geometry, optimization, and integration functions remain at the package root.

Generic graph-composition checks and independent totality and single-valuedness checks live in `semialg.function_verification`. See [Real function verification](function_verification.md) for contracts, examples and limits.

Decomposition search, fixed-degree reconstruction, parameter guards and decomposition-specific verification live in [alggeom.decomposition](https://github.com/BhuvaneshBhatt/alggeom/blob/main/docs/guide.md). These workflows use semialg for real graph and branch conditions.

Exact real normalization-image coverage lives in `semialg.curve_parametrization`. See [Real curve parametrization](curve_parametrization.md) for exceptional and isolated real points and real constants embeddings.
