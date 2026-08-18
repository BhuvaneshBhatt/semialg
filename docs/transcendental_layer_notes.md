
# Transcendental package layer (starter port)

This package introduces a separate `semialg.solve.transcendental` layer with:

- a problem-state object carrying formula, variable partition, domains, order, notes, and metadata
- a preprocessing pipeline with Piecewise simplification and family-based auxiliary replacement
- univariate transcendental root isolation scaffolding
- periodic bounding helpers based on SymPy periodicity detection
- special-function family classification and rewrite hooks
- bounded system-roots fallback using `nonlinsolve` and seeded `nsolve`
- answer cleanup and finite-point solved-form reconstruction

This is an architectural starter layer, not full full transcendental.


## Stage 19 deepening

Added, in order:

1. stronger quantifier-aware preprocessing and dispatch
   - `QuantifierBlock`
   - `normalize_quantifier_blocks`
   - `QuantifierAwareQuantifierDispatchPlan`
   - `build_quantifier_plan`

2. richer function-family variable replacement
   - auxiliary replacement adds explicit equality constraints
   - family ordering is driven by detected family density
   - quantifier-aware mode restricts replacement scope to equations

3. harder univariate inequality decomposition
   - `decompose_univariate_inequality`
   - sampled sign-change decomposition over a search window

4. better periodic solution reconstruction
   - `reconstruct_periodic_solution_from_representatives`
   - `reconstruct_periodic_intervals_from_fundamental_domain`

5. stronger system-roots fallback and certification
   - `CertifiedPoint`
   - residual-based certification metadata for fallback points


## Stage 20 deepening

Added:

- real quantifier elimination over limited genuinely-univariate transcendental families
  - `quantifier_elimination.py`
  - `QuantifierEliminationResult`
  - `eliminate_leading_real_quantifier_block(...)`

- stronger certified interval/root isolation
  - `CertifiedIntervalRoot`
  - sign-change bracketing plus bisection-based midpoint certification
  - `certified_sign_change_isolation`
  - `certified_interval_decomposition`

- broader special-family handlers
  - explicit `statistical` family detection (`erf`, `erfc`)
  - explicit `arg` family detection (`arg`)
  - direct family rewrites for `erf`, `arg`, and retained `LambertW`

- stronger multivariate orchestration and proof-style completeness tracking
  - `CompletenessCertificate`
  - `orchestrate_transcendental_system_search(...)`
  - engine-level completeness certificates propagated in `TranscendentalReductionResult`
