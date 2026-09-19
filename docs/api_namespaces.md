# API namespaces

The package root contains broad mathematical operations and common region/function APIs. Low-level algorithms, certificates, planners, diagnostics, result records, and specialized workflows live in their owning namespaces.

Specialist APIs include Thom encodings (`semialg.algebraic`), roadmap construction (`semialg.roadmaps`), Hardt trivialization, dimension/component stratification and triangulation/Betti workflows (`semialg.topology.semialgebraic`), parameter-space covers (`semialg.parametric_geometry`), root-count conditions (`semialg.parameters`), and parametric map degree (`semialg.map_degree`).

Broad convenience operations such as connected-component counting, topology summaries, solvability conditions, and common geometry, optimization, and integration functions remain at the package root.
