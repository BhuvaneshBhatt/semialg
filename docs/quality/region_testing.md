# Region testing strategy

Region tests are organized by mathematical contract under `tests/regions/` rather than by private implementation module.

- `test_equivalent_formulations.py`: representation independence.
- `test_projection_and_maps.py`: projection/image/preimage identities.
- `test_point_location.py`: exact tower signs and CAD reuse.
- `test_cell_complex.py`: dimension grouping and incidence duality.
- `test_known_region_corpus.py`: curated regions with known invariants.
- `test_singularities.py`: smooth, cusp, corner, and local-dimension distinctions.
- `test_meshing.py`: simplex validity, shared vertices, and sampled conformity.
- `test_measure.py`: integral/measure identities and reuse.
- `test_api_interoperability.py`: raw/standard/symbolic/CAD entry-point interoperability.
- `test_failure_modes.py`: conservative rejection of unsupported or unresolved input.

Metamorphic/Hypothesis tests live in `tests/properties/test_region_metamorphic.py`. Reuse contracts live with the performance tests because they protect algorithm selection and prevent accidental decomposition rebuilding without relying on fragile wall-clock thresholds.

The curated corpus includes both compact and noncompact/open examples because compact-support Euler characteristic differs from ordinary homotopy Euler characteristic on open cells.

## Long-session isolation

Long-session validation is run in a single controlled pytest process, with randomized-order CI used separately to expose process-state and cache-order bugs. The suite does not clear SymPy's global expression cache between tests: doing so would both distort normal library behavior and add substantial overhead. Performance investigations should therefore distinguish package-owned cache growth from external test-runner/process leakage, and test harnesses must terminate interrupted pytest subprocesses before starting another validation run.
