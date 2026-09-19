# High-level API default contract

These defaults are part of semialg's user-facing behavior and are guarded by tests so documentation cannot silently drift from Python signatures.

| API | Documented defaults |
|---|---|
| `cad` | `output="formula"`, `strategy="auto"`, `domain="reals"`, `diagnostics=True`, `strict=False`, `return_result=False`, `formula_form="nested"`, `max_formula_terms=512`, `max_preprocess_aux_vars=3` |
| `function_range` | `constraints=None`, `variables=None`, `domain="reals"`, `method="qe"`, `return_result=False`, `parameters=None`, `return_stratified=False`, `eliminate_quantifiers=False` |
| `semialg.parameters.root_count_conditions` | `parameters=None`, `return_result=False`, `return_stratified=False` |
| `classify_real_roots` | `parameters=None` |
| `semialgebraic_minimize` / `semialgebraic_maximize` | `constraints=None`, `variables=None`, `domain="reals"`, `return_result=False`, `certification="auto"`, `range_cost_limit=2500`, `recursion_limit=4`, `max_boolean_branches=32`, `parameters=None`, `return_stratified=False`, `eliminate_quantifiers=False` |
| `integrate_over_region` | `bounds=None`, `method="symbolic"`, `precision=50`, `measure_dimension="ambient"`, `return_result=False`, `parameters=None`, `return_stratified=False` |

When a default is changed intentionally, update this page and the corresponding contract test in the same change.
