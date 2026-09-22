# Root API usage and parameter semantics

This page deepens the per-function reference for public APIs whose signatures carry non-obvious mathematical or algorithmic choices. It complements the family reference pages rather than replacing their derivations and background. Examples below are deliberately small and are drawn from calls exercised by the regression suite; they are intended to expose semantics, not to benchmark performance.

**Exactness and certification.** Unless an API explicitly documents approximate presentation or randomized candidate search, mathematical decisions are exact. Candidate witnesses/certificates are checked before they are allowed to establish a result; resource exhaustion or an unsupported path is not silently converted to `False`.

**Algorithm/backend.** `strategy`, `method`, ordering, and resource-control parameters select *how* semialg attempts a computation. They must not change the mathematical meaning of a completed certified result. Structured-result modes expose backend/status information when it matters.

**Complexity and limitations.** CAD, quantifier elimination, exact optimization, topology, and algebraic decomposition can be intrinsically expensive. Resource guards may therefore produce an incomplete/unknown status or a documented exception. Geometry helpers also require compatible ambient variables and may decline unsupported symbolic cases rather than guess.

## real_algebraic_feasibility

`real_algebraic_feasibility(equations: 'Iterable[sp.Expr | sp.Equality]', variables: 'Sequence[sp.Symbol]', *, max_generic_attempts: 'int' = 64, max_pieces: 'int' = 128, return_result: 'bool' = False) -> 'bool | None | RealAlgebraicFeasibilityResult'`

Decide real algebraic feasibility.

**Parameter semantics.**

- `equations` — Polynomial equations defining the algebraic set.
- `variables` — Ordered real variables of the problem; their order is semantically significant for elimination, coordinates, or returned tuples.
- `max_generic_attempts` — Bound on generic choices attempted by the algebraic solver before reporting incompleteness.
- `max_pieces` — Bound on decomposition pieces retained during the computation.
- `return_result` — When true, return the structured result/certificate/diagnostic object instead of only the convenience value.

**Result semantics.** Return ``True``, ``False``, or ``None`` (inconclusive) by default. Set ``return_result=True`` for witnesses, terminal systems, and diagnostics.

**Representative regression-backed call.**

```python
real_algebraic_feasibility((x**2 - 1,), (x,))
```

This call is exercised by the regression contract “single answer decisions default to mathematical values”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## solve_real_algebraic_set

`solve_real_algebraic_set(equations: 'Iterable[sp.Expr | sp.Equality]', variables: 'Sequence[sp.Symbol]', **kwargs) -> 'dict[sp.Symbol, sp.Expr] | None'`

Return one exact real point, ``None`` for certified emptiness, or raise if incomplete.

**Parameter semantics.**

- `equations` — Polynomial equations defining the algebraic set.
- `variables` — Ordered real variables of the problem; their order is semantically significant for elimination, coordinates, or returned tuples.
- `kwargs` — Kwargs controlling this operation; see the signature type/default for the accepted representation.

**Result semantics.** The return contract is `dict[sp.Symbol, sp.Expr] | None`; see the family reference for structured-result details and certification status.

**Representative regression-backed call.**

```python
solve_real_algebraic_set((x * y,), (x, y))
```

This call is exercised by the regression contract “ars reducible positive dimensional set uses certified decomposition”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## find_negative_witness_fast

`find_negative_witness_fast(polynomial: 'sp.Expr', variables: 'Sequence[sp.Symbol]', *, random_lines: 'int' = 8, seed: 'int' = 1234) -> 'WitnessSearchResult'`

Run cheap one-sided negative-witness strategies in deterministic order.

**Parameter semantics.**

- `polynomial` — Polynomial whose sign, geometry, or algebraic structure is requested.
- `variables` — Ordered real variables of the problem; their order is semantically significant for elimination, coordinates, or returned tuples.
- `random_lines` — Number of deterministic-seeded line probes used by the fast witness-search stage.
- `seed` — Deterministic random seed used only by randomized search/generation stages.

**Result semantics.** The return contract is `WitnessSearchResult`; see the family reference for structured-result details and certification status.

**Representative regression-backed call.**

```python
find_negative_witness_fast(x**2 + y**2 - 2, (x, y))
```

This call is exercised by the regression contract “fast witness search finds simple negative ray”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## find_negative_point

`find_negative_point(polynomial: 'sp.Expr', variables: 'Sequence[sp.Symbol]', **kwargs) -> 'dict[sp.Symbol, sp.Expr] | None'`

Return a certified negative point or ``None`` for certified nonnegativity.

**Parameter semantics.**

- `polynomial` — Polynomial whose sign, geometry, or algebraic structure is requested.
- `variables` — Ordered real variables of the problem; their order is semantically significant for elimination, coordinates, or returned tuples.
- `kwargs` — Kwargs controlling this operation; see the signature type/default for the accepted representation.

**Result semantics.** The return contract is `dict[sp.Symbol, sp.Expr] | None`; see the family reference for structured-result details and certification status.

**Representative regression-backed call.**

```python
find_negative_point(x**2, (x,))
```

This call is exercised by the regression contract “find negative point handles boundary minimum and certified absence”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## polynomial_nonnegative

`polynomial_nonnegative(polynomial: 'sp.Expr', variables: 'Sequence[sp.Symbol]', *, strategy: 'str' = 'auto', sos_backend: 'str | Callable[[sp.Expr, list[sp.Symbol]], object]' = 'auto', return_result: 'bool' = False, random_lines: 'int' = 8, seed: 'int' = 1234) -> 'bool | CertifiedDecisionResult'`

Decide certified global polynomial nonnegativity.

**Parameter semantics.**

- `polynomial` — Polynomial whose sign, geometry, or algebraic structure is requested.
- `variables` — Ordered real variables of the problem; their order is semantically significant for elimination, coordinates, or returned tuples.
- `strategy` — Algorithm-selection policy. Use the default automatic policy unless a specific exact backend is being tested or diagnosed.
- `sos_backend` — SOS search backend or callable; returned certificates are still independently checked before certification.
- `return_result` — When true, return the structured result/certificate/diagnostic object instead of only the convenience value.
- `random_lines` — Number of deterministic-seeded line probes used by the fast witness-search stage.
- `seed` — Deterministic random seed used only by randomized search/generation stages.

**Result semantics.** The default return is the mathematical Boolean. Set ``return_result=True`` to obtain the structured certified portfolio trace, including the selected

**Representative regression-backed call.**

```python
polynomial_nonnegative(x**2 - 1, (x,))
```

This call is exercised by the regression contract “polynomial nonnegative distinguishes global signs”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## solvability_conditions

`solvability_conditions(constraints: 'FormulaLike | Iterable[FormulaLike]', variables: 'Sequence[sp.Symbol | str]', parameters: 'Sequence[sp.Symbol | str] | None' = None, *, domain: 'str' = 'reals', return_result: 'bool' = False, return_stratified: 'bool' = False) -> 'sp.Expr | SolvabilityConditionsResult | ParameterStratifiedResult'`

Return parameter conditions for real solvability of a constraint system.

**Parameter semantics.**

- `constraints` — Polynomial equalities/inequalities or a Boolean formula defining the feasible set.
- `variables` — Ordered real variables of the problem; their order is semantically significant for elimination, coordinates, or returned tuples.
- `parameters` — Symbols treated as free parameters rather than eliminated decision variables.
- `domain` — Mathematical domain selector; semialg public decision procedures are principally exact over the reals.
- `return_result` — When true, return the structured result/certificate/diagnostic object instead of only the convenience value.
- `return_stratified` — Request parameter-space strata rather than only their combined condition.

**Result semantics.** The return contract is `sp.Expr | SolvabilityConditionsResult | ParameterStratifiedResult`; see the family reference for structured-result details and certification status.

**Representative regression-backed call.**

```python
solvability_conditions(x**2 + a < 0, [x], [a])
```

This call is exercised by the regression contract “solvability conditions for strict negative quadratic parameter”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## prove_zero

`prove_zero(expr: 'sp.Expr', variables: 'Sequence[sp.Symbol | str] | None' = None, *, assumptions: 'FormulaLike | Iterable[FormulaLike]' = True, strategy: 'str | None' = None, parameters: 'Sequence[sp.Symbol | str] | None' = None, return_result: 'bool' = False) -> 'bool | SignProofResult | ParameterStratifiedResult'`

Return whether the expression is identically zero on the stated domain.

**Parameter semantics.**

- `expr` — Symbolic expression being analyzed.
- `variables` — Ordered real variables of the problem; their order is semantically significant for elimination, coordinates, or returned tuples.
- `assumptions` — Additional real-domain conditions under which the requested statement or computation is interpreted.
- `strategy` — Algorithm-selection policy. Use the default automatic policy unless a specific exact backend is being tested or diagnosed.
- `parameters` — Symbols treated as free parameters rather than eliminated decision variables.
- `return_result` — When true, return the structured result/certificate/diagnostic object instead of only the convenience value.

**Result semantics.** The return contract is `bool | SignProofResult | ParameterStratifiedResult`; see the family reference for structured-result details and certification status.

**Representative regression-backed call.**

```python
semialg.prove_zero(x - x, [x])
```

This call is exercised by the regression contract “sign primitives and parameter conditions”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## prove_nonzero

`prove_nonzero(expr: 'sp.Expr', variables: 'Sequence[sp.Symbol | str] | None' = None, *, assumptions: 'FormulaLike | Iterable[FormulaLike]' = True, strategy: 'str | None' = None, parameters: 'Sequence[sp.Symbol | str] | None' = None, return_result: 'bool' = False) -> 'bool | SignProofResult | ParameterStratifiedResult'`

Return whether the expression is everywhere nonzero on the stated domain.

**Parameter semantics.**

- `expr` — Symbolic expression being analyzed.
- `variables` — Ordered real variables of the problem; their order is semantically significant for elimination, coordinates, or returned tuples.
- `assumptions` — Additional real-domain conditions under which the requested statement or computation is interpreted.
- `strategy` — Algorithm-selection policy. Use the default automatic policy unless a specific exact backend is being tested or diagnosed.
- `parameters` — Symbols treated as free parameters rather than eliminated decision variables.
- `return_result` — When true, return the structured result/certificate/diagnostic object instead of only the convenience value.

**Result semantics.** The return contract is `bool | SignProofResult | ParameterStratifiedResult`; see the family reference for structured-result details and certification status.

**Representative regression-backed call.**

```python
semialg.prove_nonzero(x**2 + 1, [x])
```

This call is exercised by the regression contract “sign primitives and parameter conditions”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## equivalent

`equivalent(lhs: 'FormulaLike', rhs: 'FormulaLike', variables: 'Sequence[sp.Symbol | str] | None' = None, *, domain: 'str' = 'reals', strategy: 'str | None' = None, return_result: 'bool' = False) -> 'bool | EquivalenceResult'`

Return whether two semialgebraic formulas define the same real set.

**Parameter semantics.**

- `lhs` — Left formula in the comparison.
- `rhs` — Right formula in the comparison.
- `variables` — Ordered real variables of the problem; their order is semantically significant for elimination, coordinates, or returned tuples.
- `domain` — Mathematical domain selector; semialg public decision procedures are principally exact over the reals.
- `strategy` — Algorithm-selection policy. Use the default automatic policy unless a specific exact backend is being tested or diagnosed.
- `return_result` — When true, return the structured result/certificate/diagnostic object instead of only the convenience value.

**Result semantics.** With ``return_result=True``, a false result includes a counterexample from the symmetric difference and, when it can be determined cheaply, the failed

**Representative regression-backed call.**

```python
equivalent(fm, cad, [y])
```

This call is exercised by the regression contract “fourier motzkin matches complete qe on linear projection”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## solve_semialgebraic

`solve_semialgebraic(constraints: 'FormulaLike | Iterable[FormulaLike]', variables: 'Sequence[sp.Symbol | str] | None' = None, *, parameters: 'Sequence[sp.Symbol | str] | None' = None, domain: 'str' = 'reals', count: 'int' = 1, samples: 'int | str | None' = None, sample_mode: 'str | None' = None, strategy: 'str | None' = None, method: 'str' = 'auto', variable_order: 'Sequence[sp.Symbol | str] | None' = None, projection_order: 'Sequence[sp.Symbol | str] | None' = None, normalize_domains: 'bool' = True, return_formula: 'bool' = False, output: 'str | None' = None) -> 'SemialgebraicSolution | sp.Expr | tuple[object, ...] | bool | None'`

Reduce, sample, and summarize a semialgebraic system over the reals.

**Parameter semantics.**

- `constraints` — Polynomial equalities/inequalities or a Boolean formula defining the feasible set.
- `variables` — Ordered real variables of the problem; their order is semantically significant for elimination, coordinates, or returned tuples.
- `parameters` — Symbols treated as free parameters rather than eliminated decision variables.
- `domain` — Mathematical domain selector; semialg public decision procedures are principally exact over the reals.
- `count` — Requested number of samples/solutions when the API supports bounded enumeration.
- `samples` — Sampling request or sampling policy for the returned solution representation.
- `sample_mode` — Policy controlling whether samples are global, per cell, or per connected component.
- `strategy` — Algorithm-selection policy. Use the default automatic policy unless a specific exact backend is being tested or diagnosed.
- `method` — Computation-method selector; automatic mode chooses an applicable certified route.
- `variable_order` — Explicit variable order for decomposition/solving; changes work ordering, not mathematical meaning.
- `projection_order` — Explicit elimination/projection order used by CAD-style algorithms.
- `normalize_domains` — Whether domain predicates are normalized before the main exact computation.
- `return_formula` — Request the formula view instead of the richer structured solution object.
- `output` — Requested public representation/view of the computed result.

**Result semantics.** ``output`` may be used as a convenience selector for common views of the solution. The default ``None`` preserves structured result-object behavior.

**Representative regression-backed call.**

```python
solve_semialgebraic(x**2 <= 1, [x])
```

This call is exercised by the regression contract “solve semialgebraic reduces univariate interval”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## apply_quantifiers

`apply_quantifiers(formula: 'sp.Basic | bool', quantifiers: 'Iterable[tuple[str, sp.Symbol]]') -> 'Boolean'`

Wrap ``formula`` in a prenex quantifier prefix.

**Parameter semantics.**

- `formula` — Boolean/relational semialgebraic formula defining the set or proposition.
- `quantifiers` — Ordered quantifier prefix, from outermost to innermost.

**Result semantics.** The iterable is ordered from outermost to innermost, matching semialg's existing internal ``(name, symbol)`` quantifier representation.

**Representative regression-backed call.**

```python
apply_quantifiers(matrix, prefix)
```

This call is exercised by the regression contract “apply and split quantifiers round trip”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## prove_negative

`prove_negative(expr: 'sp.Expr', variables: 'Sequence[sp.Symbol | str] | None' = None, *, assumptions: 'FormulaLike | Iterable[FormulaLike]' = True, strategy: 'str | None' = None, parameters: 'Sequence[sp.Symbol | str] | None' = None, return_result: 'bool' = False) -> 'bool | SignProofResult | ParameterStratifiedResult'`

Return whether the expression is certified negative on the stated domain.

**Parameter semantics.**

- `expr` — Symbolic expression being analyzed.
- `variables` — Ordered real variables of the problem; their order is semantically significant for elimination, coordinates, or returned tuples.
- `assumptions` — Additional real-domain conditions under which the requested statement or computation is interpreted.
- `strategy` — Algorithm-selection policy. Use the default automatic policy unless a specific exact backend is being tested or diagnosed.
- `parameters` — Symbols treated as free parameters rather than eliminated decision variables.
- `return_result` — When true, return the structured result/certificate/diagnostic object instead of only the convenience value.

**Result semantics.** The return contract is `bool | SignProofResult | ParameterStratifiedResult`; see the family reference for structured-result details and certification status.

**Representative regression-backed call.**

```python
prove_negative(-(x**2), [x])
```

This call is exercised by the regression contract “boolean api preserved for sign provers”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## prove_nonnegative

`prove_nonnegative(expr: 'sp.Expr', variables: 'Sequence[sp.Symbol | str] | None' = None, *, assumptions: 'FormulaLike | Iterable[FormulaLike]' = True, strategy: 'str | None' = None, parameters: 'Sequence[sp.Symbol | str] | None' = None, return_result: 'bool' = False) -> 'bool | SignProofResult | ParameterStratifiedResult'`

Return whether the expression is certified nonnegative on the stated domain.

**Parameter semantics.**

- `expr` — Symbolic expression being analyzed.
- `variables` — Ordered real variables of the problem; their order is semantically significant for elimination, coordinates, or returned tuples.
- `assumptions` — Additional real-domain conditions under which the requested statement or computation is interpreted.
- `strategy` — Algorithm-selection policy. Use the default automatic policy unless a specific exact backend is being tested or diagnosed.
- `parameters` — Symbols treated as free parameters rather than eliminated decision variables.
- `return_result` — When true, return the structured result/certificate/diagnostic object instead of only the convenience value.

**Result semantics.** The return contract is `bool | SignProofResult | ParameterStratifiedResult`; see the family reference for structured-result details and certification status.

**Representative regression-backed call.**

```python
prove_nonnegative(x**2, [x])
```

This call is exercised by the regression contract “boolean api preserved for sign provers”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## prove_nonpositive

`prove_nonpositive(expr: 'sp.Expr', variables: 'Sequence[sp.Symbol | str] | None' = None, *, assumptions: 'FormulaLike | Iterable[FormulaLike]' = True, strategy: 'str | None' = None, parameters: 'Sequence[sp.Symbol | str] | None' = None, return_result: 'bool' = False) -> 'bool | SignProofResult | ParameterStratifiedResult'`

Return whether the expression is certified nonpositive on the stated domain.

**Parameter semantics.**

- `expr` — Symbolic expression being analyzed.
- `variables` — Ordered real variables of the problem; their order is semantically significant for elimination, coordinates, or returned tuples.
- `assumptions` — Additional real-domain conditions under which the requested statement or computation is interpreted.
- `strategy` — Algorithm-selection policy. Use the default automatic policy unless a specific exact backend is being tested or diagnosed.
- `parameters` — Symbols treated as free parameters rather than eliminated decision variables.
- `return_result` — When true, return the structured result/certificate/diagnostic object instead of only the convenience value.

**Result semantics.** The return contract is `bool | SignProofResult | ParameterStratifiedResult`; see the family reference for structured-result details and certification status.

**Representative regression-backed call.**

```python
prove_nonpositive(-(x**2), [x])
```

This call is exercised by the regression contract “inequality provers”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## prove_positive

`prove_positive(expr: 'sp.Expr', variables: 'Sequence[sp.Symbol | str] | None' = None, *, assumptions: 'FormulaLike | Iterable[FormulaLike]' = True, strategy: 'str | None' = None, parameters: 'Sequence[sp.Symbol | str] | None' = None, return_result: 'bool' = False) -> 'bool | SignProofResult | ParameterStratifiedResult'`

Return whether the expression is certified positive on the stated domain.

**Parameter semantics.**

- `expr` — Symbolic expression being analyzed.
- `variables` — Ordered real variables of the problem; their order is semantically significant for elimination, coordinates, or returned tuples.
- `assumptions` — Additional real-domain conditions under which the requested statement or computation is interpreted.
- `strategy` — Algorithm-selection policy. Use the default automatic policy unless a specific exact backend is being tested or diagnosed.
- `parameters` — Symbols treated as free parameters rather than eliminated decision variables.
- `return_result` — When true, return the structured result/certificate/diagnostic object instead of only the convenience value.

**Result semantics.** The return contract is `bool | SignProofResult | ParameterStratifiedResult`; see the family reference for structured-result details and certification status.

**Representative regression-backed call.**

```python
prove_positive(x**2, [x])
```

This call is exercised by the regression contract “boolean api preserved for sign provers”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## reduce_formula

`reduce_formula(parsed: 'ParsedPrenexFormula', config=None, *, domain: 'str | SolveDomain | None' = None, assumptions=None, return_result: 'bool' = False, strategy: 'str | None' = None)`

Reduce a parsed real formula using the selected exact decision strategy.

**Parameter semantics.**

- `parsed` — Parsed prenex formula produced by semialg parsing/normalization machinery.
- `config` — Optional solver/CAD configuration; omitted values use the package defaults.
- `domain` — Mathematical domain selector; semialg public decision procedures are principally exact over the reals.
- `assumptions` — Additional real-domain conditions under which the requested statement or computation is interpreted.
- `return_result` — When true, return the structured result/certificate/diagnostic object instead of only the convenience value.
- `strategy` — Algorithm-selection policy. Use the default automatic policy unless a specific exact backend is being tested or diagnosed.

**Result semantics.** The return contract is `the documented result type for this operation`; see the family reference for structured-result details and certification status.

**Representative regression-backed call.**

```python
reduce_formula(true_existential)
```

This call is exercised by the regression contract “reduce and resolve formula cover false universal and detailed result”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## resolve_formula

`resolve_formula(parsed, config=None, *, domain: 'str | SolveDomain | None' = None, assumptions=None, return_result: 'bool' = False, strategy: 'str | None' = 'lazy')`

Resolve a parsed formula and return its exact solution representation.

**Parameter semantics.**

- `parsed` — Parsed prenex formula produced by semialg parsing/normalization machinery.
- `config` — Optional solver/CAD configuration; omitted values use the package defaults.
- `domain` — Mathematical domain selector; semialg public decision procedures are principally exact over the reals.
- `assumptions` — Additional real-domain conditions under which the requested statement or computation is interpreted.
- `return_result` — When true, return the structured result/certificate/diagnostic object instead of only the convenience value.
- `strategy` — Algorithm-selection policy. Use the default automatic policy unless a specific exact backend is being tested or diagnosed.

**Result semantics.** The return contract is `the documented result type for this operation`; see the family reference for structured-result details and certification status.

**Representative regression-backed call.**

```python
resolve_formula(false_universal)
```

This call is exercised by the regression contract “reduce and resolve formula cover false universal and detailed result”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## cad

`cad(formula: 'sp.Expr | Formula', variables: 'Sequence[sp.Symbol | str]', *, output: 'CADOutput' = 'formula', operation: 'TopoOp | None' = None, strategy: 'str' = 'auto', domain: 'str' = 'reals', assumptions: 'Iterable[sp.Expr] | sp.Expr | None' = None, max_cells: 'int | None' = None, timeout: 'float | None' = None, diagnostics: 'bool' = True, strict: 'bool' = False, return_result: 'bool' = False, formula_form: 'FormulaForm' = 'nested', max_formula_terms: 'int' = 512, max_preprocess_aux_vars: 'int | None' = 3)`

Compute a cylindrical algebraic decomposition for a real formula.

**Parameter semantics.**

- `formula` — Boolean/relational semialgebraic formula defining the set or proposition.
- `variables` — Ordered real variables of the problem; their order is semantically significant for elimination, coordinates, or returned tuples.
- `output` — Requested public representation/view of the computed result.
- `operation` — Topological/set operation requested from CAD.
- `strategy` — Algorithm-selection policy. Use the default automatic policy unless a specific exact backend is being tested or diagnosed.
- `domain` — Mathematical domain selector; semialg public decision procedures are principally exact over the reals.
- `assumptions` — Additional real-domain conditions under which the requested statement or computation is interpreted.
- `max_cells` — Resource guard on CAD cell construction; hitting it must not be interpreted as a mathematical false result.
- `timeout` — Optional resource time limit; exhaustion is reported as incomplete rather than silently certified.
- `diagnostics` — Whether to retain diagnostic information in structured results.
- `strict` — Select strict failure behavior for incomplete/unsupported computation where the API provides it.
- `return_result` — When true, return the structured result/certificate/diagnostic object instead of only the convenience value.
- `formula_form` — Requested shape of reconstructed Boolean formulas.
- `max_formula_terms` — Guard against explosive reconstructed-formula size.
- `max_preprocess_aux_vars` — Limit on auxiliary variables introduced during preprocessing.

**Result semantics.** By default the representation selected by ``output`` is returned directly: a formula, :class:`CellSet`, :class:`CADFunction`, or CAD tree.  Set

**Representative regression-backed call.**

```python
cad(formula, (X,))
```

This call is exercised by the regression contract “unimplemented cad budgets are never silently ignored”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## implicitize_polynomial_map

`implicitize_polynomial_map(mapping, parameters, *, image_variables=None, domain_equations=(), return_result: 'bool' = False)`

Implicitize a polynomial map by exact Groebner elimination.

**Parameter semantics.**

- `mapping` — Polynomial or affine map whose image/local geometry is analyzed.
- `parameters` — Symbols treated as free parameters rather than eliminated decision variables.
- `image_variables` — Coordinate symbols used for the image space.
- `domain_equations` — Equations restricting the parameter/domain variety.
- `return_result` — When true, return the structured result/certificate/diagnostic object instead of only the convenience value.

**Result semantics.** ``domain_equations`` may restrict the parameter space to an algebraic subvariety.  Inequalities are intentionally not accepted here: an arbitrary

**Representative regression-backed call.**

```python
implicitize_polynomial_map((sp.sin(t),), (t,))
```

This call is exercised by the regression contract “nonpolynomial mapping is rejected”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## minimal_prime_intersections

`minimal_prime_intersections(equations, variables=None, *, max_pieces: 'int | None' = None, max_order: 'int | None' = None) -> 'tuple[MinimalPrimeIntersection, ...]'`

Return exact intersections among certified minimal-prime components.

**Parameter semantics.**

- `equations` — Polynomial equations defining the algebraic set.
- `variables` — Ordered real variables of the problem; their order is semantically significant for elimination, coordinates, or returned tuples.
- `max_pieces` — Bound on decomposition pieces retained during the computation.
- `max_order` — Maximum algebraic/moment order considered.

**Result semantics.** The intersection of varieties is represented by the sum of their prime ideals.  Reduced Groebner generators make the result deterministic.

**Representative regression-backed call.**

```python
minimal_prime_intersections((x * y * z,), (x, y, z))
```

This call is exercised by the regression contract “minimal prime intersections respect max order”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## local_dimension_strata

`local_dimension_strata(equations, variables=None, *, max_pieces: 'int | None' = None) -> 'tuple[LocalDimensionStratum, ...]'`

Partition a reduced variety into constructible constant-local-dimension loci.

**Parameter semantics.**

- `equations` — Polynomial equations defining the algebraic set.
- `variables` — Ordered real variables of the problem; their order is semantically significant for elimination, coordinates, or returned tuples.
- `max_pieces` — Bound on decomposition pieces retained during the computation.

**Result semantics.** At a point of a reduced affine variety, local dimension is the maximum of the dimensions of the irreducible components through that point.  The

**Representative regression-backed call.**

```python
semialg.local_dimension_strata((x * y,), (x, y))
```

This call is exercised by the regression contract “algebraic statistics secondary calls”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## local_branch_geometry

`local_branch_geometry(equations, point, variables=None, *, max_pieces: 'int | None' = None) -> 'LocalBranchGeometry'`

Return exact local branch geometry at a point of a reduced variety.

**Parameter semantics.**

- `equations` — Polynomial equations defining the algebraic set.
- `point` — Point in the ambient coordinate system.
- `variables` — Ordered real variables of the problem; their order is semantically significant for elimination, coordinates, or returned tuples.
- `max_pieces` — Bound on decomposition pieces retained during the computation.

**Result semantics.** Each minimal-prime branch through the point is analyzed independently. ``tangent_dimension_excess`` is the excess of Zariski tangent dimension over

**Representative regression-backed call.**

```python
local_branch_geometry((y**2,), {x: 0, y: 0}, (x, y))
```

This call is exercised by the regression contract “nonradical presentation does not inflate local geometry”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## stratified_singular_geometry

`stratified_singular_geometry(equations, variables=None, *, max_pieces: 'int | None' = None, max_intersection_order: 'int | None' = None) -> 'StratifiedSingularGeometry'`

Return a certified stratified description of singular algebraic geometry.

**Parameter semantics.**

- `equations` — Polynomial equations defining the algebraic set.
- `variables` — Ordered real variables of the problem; their order is semantically significant for elimination, coordinates, or returned tuples.
- `max_pieces` — Bound on decomposition pieces retained during the computation.
- `max_intersection_order` — Maximum intersection order considered.

**Result semantics.** The result separates four notions that a global Jacobian singular locus conflates: certified minimal-prime components, their mutual intersections,

**Representative regression-backed call.**

```python
stratified_singular_geometry((x**2,), (x,))
```

This call is exercised by the regression contract “nonradical input is reduced before stratification”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## certified_radicalization

`certified_radicalization(equations, variables=None, *, max_pieces: 'int | None' = None) -> 'ReducedAlgebraicVariety'`

Return a certified reduced presentation ``sqrt(<equations>)``.

**Parameter semantics.**

- `equations` — Polynomial equations defining the algebraic set.
- `variables` — Ordered real variables of the problem; their order is semantically significant for elimination, coordinates, or returned tuples.
- `max_pieces` — Bound on decomposition pieces retained during the computation.

**Result semantics.** The function refuses to return uncertified generators.  This makes the reduced ideal, rather than a possibly non-radical input presentation, the

**Representative regression-backed call.**

```python
certified_radicalization((x**2,), (x, y))
```

This call is exercised by the regression contract “certified radicalization removes nilpotent multiplicity”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## irreducible_components

`irreducible_components(equations, variables=None, *, max_pieces: 'int | None' = None) -> 'tuple[IrreducibleAlgebraicComponent, ...]'`

Return all certified irreducible components of an affine variety.

**Parameter semantics.**

- `equations` — Polynomial equations defining the algebraic set.
- `variables` — Ordered real variables of the problem; their order is semantically significant for elimination, coordinates, or returned tuples.
- `max_pieces` — Bound on decomposition pieces retained during the computation.

**Result semantics.** Components are returned only when semialg certifies a complete, irredundant minimal-prime decomposition of the radical ideal.

**Representative regression-backed call.**

```python
irreducible_components((x * y,), (x, y))
```

This call is exercised by the regression contract “irreducible components ignore nonzero scalar and multiplicity”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## reduced_component_singular_loci

`reduced_component_singular_loci(equations, variables=None, *, max_pieces: 'int | None' = None) -> 'tuple[ReducedComponentSingularLocus, ...]'`

Return singular loci computed separately on certified reduced components.

**Parameter semantics.**

- `equations` — Polynomial equations defining the algebraic set.
- `variables` — Ordered real variables of the problem; their order is semantically significant for elimination, coordinates, or returned tuples.
- `max_pieces` — Bound on decomposition pieces retained during the computation.

**Result semantics.** The return contract is `tuple[ReducedComponentSingularLocus, ...]`; see the family reference for structured-result details and certification status.

**Representative regression-backed call.**

```python
reduced_component_singular_loci((x**2 * y**2,), (x, y))
```

This call is exercised by the regression contract “component crossing is not component singularity”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## zariski_closure

`zariski_closure(mapping, parameters, *, image_variables=None, domain_equations=(), return_result: 'bool' = False)`

Return the Zariski closure of a polynomial-map image.

**Parameter semantics.**

- `mapping` — Polynomial or affine map whose image/local geometry is analyzed.
- `parameters` — Symbols treated as free parameters rather than eliminated decision variables.
- `image_variables` — Coordinate symbols used for the image space.
- `domain_equations` — Equations restricting the parameter/domain variety.
- `return_result` — When true, return the structured result/certificate/diagnostic object instead of only the convenience value.

**Result semantics.** The closure is represented by the elimination ideal of the graph of the map (and optional algebraic parameter-domain equations).  The default

**Representative regression-backed call.**

```python
zariski_closure((s, t), (s, t), image_variables=(x, y))
```

This call is exercised by the regression contract “dominant polynomial map has full ambient zariski closure”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## classify_real_roots

`classify_real_roots(polynomial: 'sp.Poly | sp.Expr', variable: 'sp.Symbol | str', *, parameters: 'Sequence[sp.Symbol | str] | None' = None) -> 'RootClassificationResult'`

Classify real roots of a univariate polynomial or polynomial family.

**Parameter semantics.**

- `polynomial` — Polynomial whose sign, geometry, or algebraic structure is requested.
- `variable` — Single variable singled out for range/fiber/elimination analysis.
- `parameters` — Symbols treated as free parameters rather than eliminated decision variables.

**Result semantics.** The current public implementation is exact for unparameterized polynomials, complete for linear, quadratic, and cubic parameter families (including

**Representative regression-backed call.**

```python
classify_real_roots(X**2 + 1, X)
```

This call is exercised by the regression contract “root classifier rejects nonpolynomials”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## zeng_negative_point

`zeng_negative_point(polynomial: 'sp.Expr', variables: 'Sequence[sp.Symbol]', *, random_lines: 'int' = 8, seed: 'int' = 1234) -> 'PolynomialNegativityResult'`

Decide whether a rational polynomial is negative somewhere when certified.

**Parameter semantics.**

- `polynomial` — Polynomial whose sign, geometry, or algebraic structure is requested.
- `variables` — Ordered real variables of the problem; their order is semantically significant for elimination, coordinates, or returned tuples.
- `random_lines` — Number of deterministic-seeded line probes used by the fast witness-search stage.
- `seed` — Deterministic random seed used only by randomized search/generation stages.

**Result semantics.** The implementation is a Zeng-family inspired exact critical-value reduction; see [ZZ2004] and [ZX2012] in ``docs/references.md``.  It is not a literal

**Representative regression-backed call.**

```python
zeng_negative_point(x**3 + y**2, (x, y))
```

This call is exercised by the regression contract “zeng odd degree fast path is one sided and exact”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## function_range

`function_range(expression: 'sp.Expr', constraints: 'FormulaLike | Iterable[FormulaLike] | None' = None, variables: 'Sequence[sp.Symbol | str] | None' = None, *, value_symbol: 'sp.Symbol | str | None' = None, domain: 'str' = 'reals', method: 'str' = 'qe', return_result: 'bool' = False, parameters: 'Sequence[sp.Symbol | str] | None' = None, return_stratified: 'bool' = False, eliminate_quantifiers: 'bool' = False) -> 'sp.Expr | FunctionRangeResult | object'`

Return a quantifier-free formula describing a real function range.

**Parameter semantics.**

- `expression` — Symbolic expression being analyzed.
- `constraints` — Polynomial equalities/inequalities or a Boolean formula defining the feasible set.
- `variables` — Ordered real variables of the problem; their order is semantically significant for elimination, coordinates, or returned tuples.
- `value_symbol` — Symbol used to represent an exact objective/range value.
- `domain` — Mathematical domain selector; semialg public decision procedures are principally exact over the reals.
- `method` — Computation-method selector; automatic mode chooses an applicable certified route.
- `return_result` — When true, return the structured result/certificate/diagnostic object instead of only the convenience value.
- `parameters` — Symbols treated as free parameters rather than eliminated decision variables.
- `return_stratified` — Request parameter-space strata rather than only their combined condition.
- `eliminate_quantifiers` — Whether quantifiers are eliminated rather than preserved in simplified form.

**Result semantics.** The preferred direct backend uses the semialgebraic image formulation ``exists variables. constraints and value_symbol == expression``. It first

**Representative regression-backed call.**

```python
function_range(1 / x, x > 0, [x], value_symbol=t)
```

This call is exercised by the regression contract “function range univariate rational positive ray”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## semialgebraic_maximize

`semialgebraic_maximize(objective: 'sp.Expr', constraints: 'FormulaLike | Iterable[FormulaLike] | None' = None, variables: 'Sequence[sp.Symbol | str] | None' = None, *, domain: 'str' = 'reals', return_result: 'bool' = False, certification: "Literal['auto', 'complete', 'candidate']" = 'auto', range_cost_limit: 'int' = 2500, recursion_limit: 'int' = 4, max_boolean_branches: 'int' = 32, parameters: 'Sequence[sp.Symbol | str] | None' = None, return_stratified: 'bool' = False, eliminate_quantifiers: 'bool' = False) -> 'list[object] | OptimizationResult | object'`

Return ``[maximum_or_supremum, optimizer_points]`` by default.

**Parameter semantics.**

- `objective` — Exact symbolic objective optimized over the feasible region.
- `constraints` — Polynomial equalities/inequalities or a Boolean formula defining the feasible set.
- `variables` — Ordered real variables of the problem; their order is semantically significant for elimination, coordinates, or returned tuples.
- `domain` — Mathematical domain selector; semialg public decision procedures are principally exact over the reals.
- `return_result` — When true, return the structured result/certificate/diagnostic object instead of only the convenience value.
- `certification` — Controls whether the operation must return/verify exact certification metadata.
- `range_cost_limit` — Resource guard for exact range computation.
- `recursion_limit` — Explicit recursion/decomposition guard.
- `max_boolean_branches` — Guard on Boolean branch expansion.
- `parameters` — Symbols treated as free parameters rather than eliminated decision variables.
- `return_stratified` — Request parameter-space strata rather than only their combined condition.
- `eliminate_quantifiers` — Whether quantifiers are eliminated rather than preserved in simplified form.

**Result semantics.** Set ``return_result=True`` for an :class:`OptimizationResult` with attainment, certification, method, and diagnostic metadata.

**Representative regression-backed call.**

```python
semialgebraic_maximize(x, [x >= 0, x <= 2], [x])
```

This call is exercised by the regression contract “optimization defaults to value and optimizer points”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## semialgebraic_minimize

`semialgebraic_minimize(objective: 'sp.Expr', constraints: 'FormulaLike | Iterable[FormulaLike] | None' = None, variables: 'Sequence[sp.Symbol | str] | None' = None, *, domain: 'str' = 'reals', return_result: 'bool' = False, certification: "Literal['auto', 'complete', 'candidate']" = 'auto', range_cost_limit: 'int' = 2500, recursion_limit: 'int' = 4, max_boolean_branches: 'int' = 32, parameters: 'Sequence[sp.Symbol | str] | None' = None, return_stratified: 'bool' = False, eliminate_quantifiers: 'bool' = False) -> 'list[object] | OptimizationResult | object'`

Return ``[minimum_or_infimum, optimizer_points]`` by default.

**Parameter semantics.**

- `objective` — Exact symbolic objective optimized over the feasible region.
- `constraints` — Polynomial equalities/inequalities or a Boolean formula defining the feasible set.
- `variables` — Ordered real variables of the problem; their order is semantically significant for elimination, coordinates, or returned tuples.
- `domain` — Mathematical domain selector; semialg public decision procedures are principally exact over the reals.
- `return_result` — When true, return the structured result/certificate/diagnostic object instead of only the convenience value.
- `certification` — Controls whether the operation must return/verify exact certification metadata.
- `range_cost_limit` — Resource guard for exact range computation.
- `recursion_limit` — Explicit recursion/decomposition guard.
- `max_boolean_branches` — Guard on Boolean branch expansion.
- `parameters` — Symbols treated as free parameters rather than eliminated decision variables.
- `return_stratified` — Request parameter-space strata rather than only their combined condition.
- `eliminate_quantifiers` — Whether quantifiers are eliminated rather than preserved in simplified form.

**Result semantics.** Set ``return_result=True`` for an :class:`OptimizationResult` with attainment, certification, method, and diagnostic metadata.

**Representative regression-backed call.**

```python
semialgebraic_minimize((x - 2) ** 2, variables=[x])
```

This call is exercised by the regression contract “optimization defaults to value and optimizer points”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## function_smoothness

`function_smoothness(expression, variables: 'Sequence[sp.Symbol | str] | sp.Symbol | str | None' = None, *, domain=True, max_order: 'int' = 3) -> 'FunctionSmoothnessResult'`

Report continuity, C^k order, smoothness, and exact exceptional loci.

**Parameter semantics.**

- `expression` — Symbolic expression being analyzed.
- `variables` — Ordered real variables of the problem; their order is semantically significant for elimination, coordinates, or returned tuples.
- `domain` — Mathematical domain selector; semialg public decision procedures are principally exact over the reals.
- `max_order` — Maximum algebraic/moment order considered.

**Result semantics.** The current exact implementation is strongest for univariate semialgebraic expressions and for polynomial/rational multivariate functions.  ``max_order``

**Representative regression-backed call.**

```python
function_smoothness(x**3, x)
```

This call is exercised by the regression contract “function smoothness polynomial abs sign piecewise”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## function_mapping_properties

`function_mapping_properties(mapping, variables: 'Sequence[sp.Symbol | str] | sp.Symbol | str | None' = None, *, domain=True, codomain=True, image_variables: 'Sequence[sp.Symbol | str] | None' = None) -> 'FunctionMappingPropertiesResult'`

Certify injectivity, surjectivity, and bijectivity of a semialgebraic map.

**Parameter semantics.**

- `mapping` — Polynomial or affine map whose image/local geometry is analyzed.
- `variables` — Ordered real variables of the problem; their order is semantically significant for elimination, coordinates, or returned tuples.
- `domain` — Mathematical domain selector; semialg public decision procedures are principally exact over the reals.
- `codomain` — Target/codomain description of a mapping.
- `image_variables` — Coordinate symbols used for the image space.

**Result semantics.** The return contract is `FunctionMappingPropertiesResult`; see the family reference for structured-result details and certification status.

**Representative regression-backed call.**

```python
function_mapping_properties(x, x)
```

This call is exercised by the regression contract “function mapping properties scalar maps”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## centroid

`centroid(condition: 'object', variables: 'Sequence[sp.Symbol | str]', *, bounds: 'Sequence[tuple[sp.Symbol | str, object, object]] | Mapping[sp.Symbol | str, tuple[object, object]] | None' = None, method: 'str' = 'symbolic', precision: 'int' = 50, measure_dimension: 'object' = 'ambient', return_result: 'bool' = False) -> 'Mapping[sp.Symbol, sp.Expr] | RegionCentroidResult'`

Return the centroid of a finite-measure semialgebraic region.

**Parameter semantics.**

- `condition` — Semialgebraic condition restricting the computation.
- `variables` — Ordered real variables of the problem; their order is semantically significant for elimination, coordinates, or returned tuples.
- `bounds` — Exact coordinate bounds restricting search or construction.
- `method` — Computation-method selector; automatic mode chooses an applicable certified route.
- `precision` — Presentation/numerical precision used only where approximation is explicitly part of the API.
- `measure_dimension` — Dimension of the measure used for moments/centroids; useful for intrinsic lower-dimensional sets.
- `return_result` — When true, return the structured result/certificate/diagnostic object instead of only the convenience value.

**Result semantics.** The return contract is `Mapping[sp.Symbol, sp.Expr] | RegionCentroidResult`; see the family reference for structured-result details and certification status.

**Representative regression-backed call.**

```python
box.centroid()
```

This call is exercised by the regression contract “box measure and centroid”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## covariance_matrix

`covariance_matrix(condition: 'object', variables: 'Sequence[sp.Symbol | str]', *, bounds: 'Sequence[tuple[sp.Symbol | str, object, object]] | Mapping[sp.Symbol | str, tuple[object, object]] | None' = None, method: 'str' = 'symbolic', precision: 'int' = 50, measure_dimension: 'object' = 'ambient', return_result: 'bool' = False) -> 'sp.Matrix | RegionCovarianceResult'`

Return the covariance matrix of the uniform measure on a region.

**Parameter semantics.**

- `condition` — Semialgebraic condition restricting the computation.
- `variables` — Ordered real variables of the problem; their order is semantically significant for elimination, coordinates, or returned tuples.
- `bounds` — Exact coordinate bounds restricting search or construction.
- `method` — Computation-method selector; automatic mode chooses an applicable certified route.
- `precision` — Presentation/numerical precision used only where approximation is explicitly part of the API.
- `measure_dimension` — Dimension of the measure used for moments/centroids; useful for intrinsic lower-dimensional sets.
- `return_result` — When true, return the structured result/certificate/diagnostic object instead of only the convenience value.

**Result semantics.** The return contract is `sp.Matrix | RegionCovarianceResult`; see the family reference for structured-result details and certification status.

**Representative regression-backed call.**

```python
covariance_matrix(disk, [x, y])
```

This call is exercised by the regression contract “disk moments centroid covariance”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## semialgebraic_measure

`semialgebraic_measure(condition: 'object', variables: 'Sequence[sp.Symbol | str]', *, bounds: 'Sequence[tuple[sp.Symbol | str, object, object]] | Mapping[sp.Symbol | str, tuple[object, object]] | None' = None, measure_dimension: 'object' = 'ambient', return_result: 'bool' = False, parameters: 'Sequence[sp.Symbol | str] | None' = None, return_stratified: 'bool' = False) -> 'sp.Expr | MeasureResult | object'`

Return the exact measure of a supported semialgebraic set.

**Parameter semantics.**

- `condition` — Semialgebraic condition restricting the computation.
- `variables` — Ordered real variables of the problem; their order is semantically significant for elimination, coordinates, or returned tuples.
- `bounds` — Exact coordinate bounds restricting search or construction.
- `measure_dimension` — Dimension of the measure used for moments/centroids; useful for intrinsic lower-dimensional sets.
- `return_result` — When true, return the structured result/certificate/diagnostic object instead of only the convenience value.
- `parameters` — Symbols treated as free parameters rather than eliminated decision variables.
- `return_stratified` — Request parameter-space strata rather than only their combined condition.

**Result semantics.** The measure implementation delegates to the structural region-integral reducer with integrand ``1``. This keeps ``semialgebraic_measure`` aligned

**Representative regression-backed call.**

```python
semialgebraic_measure(moved, [x])
```

This call is exercised by the regression contract “generated affine translation invariants”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## region_measure

`region_measure(region: 'object', variables: 'Sequence[sp.Symbol | str] | None' = None, *, measure_dimension: 'object' = None, return_result: 'bool' = False) -> 'sp.Expr | MeasureResult'`

Return exact Euclidean/Hausdorff measure of a region when supported.

**Parameter semantics.**

- `region` — Semialgebraic or standard region on which the operation is performed.
- `variables` — Ordered real variables of the problem; their order is semantically significant for elimination, coordinates, or returned tuples.
- `measure_dimension` — Dimension of the measure used for moments/centroids; useful for intrinsic lower-dimensional sets.
- `return_result` — When true, return the structured result/certificate/diagnostic object instead of only the convenience value.

**Result semantics.** With ``measure_dimension=None``, canonical geometry uses intrinsic measure while formula regions use ambient Lebesgue measure. Pass ``"intrinsic"``,

**Representative regression-backed call.**

```python
region_measure(box)
```

This call is exercised by the regression contract “box measure and centroid”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## region_moment

`region_moment(condition: 'object', variables: 'Sequence[sp.Symbol | str]', powers: 'Sequence[int] | sp.Expr | None' = None, *, integrand: 'object | None' = None, bounds: 'Sequence[tuple[sp.Symbol | str, object, object]] | Mapping[sp.Symbol | str, tuple[object, object]] | None' = None, method: 'str' = 'symbolic', precision: 'int' = 50, measure_dimension: 'object' = 'ambient', return_result: 'bool' = False) -> 'sp.Expr | RegionMomentResult'`

Return a raw moment integral over a semialgebraic region.

**Parameter semantics.**

- `condition` — Semialgebraic condition restricting the computation.
- `variables` — Ordered real variables of the problem; their order is semantically significant for elimination, coordinates, or returned tuples.
- `powers` — Moment exponents or polynomial powers requested.
- `integrand` — Exact symbolic integrand.
- `bounds` — Exact coordinate bounds restricting search or construction.
- `method` — Computation-method selector; automatic mode chooses an applicable certified route.
- `precision` — Presentation/numerical precision used only where approximation is explicitly part of the API.
- `measure_dimension` — Dimension of the measure used for moments/centroids; useful for intrinsic lower-dimensional sets.
- `return_result` — When true, return the structured result/certificate/diagnostic object instead of only the convenience value.

**Result semantics.** ``powers`` gives a monomial moment. For variables ``(x, y)`` and powers ``(2, 1)``, the integrated moment is ``x**2*y``. Alternatively, callers

**Representative regression-backed call.**

```python
region_moment(x**2 <= 1, [x])
```

This call is exercised by the regression contract “raw moments on interval”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## integrate_over_region

`integrate_over_region(integrand: 'object', condition: 'object', variables: 'Sequence[sp.Symbol | str]', *, bounds: 'Sequence[tuple[sp.Symbol | str, object, object]] | Mapping[sp.Symbol | str, tuple[object, object]] | None' = None, method: 'str' = 'symbolic', precision: 'int' = 50, measure_dimension: 'object' = 'ambient', return_result: 'bool' = False, parameters: 'Sequence[sp.Symbol | str] | None' = None, return_stratified: 'bool' = False) -> 'sp.Expr | RegionIntegralResult | object'`

Integrate ``integrand`` over a supported semialgebraic region.

**Parameter semantics.**

- `integrand` — Exact symbolic integrand.
- `condition` — Semialgebraic condition restricting the computation.
- `variables` — Ordered real variables of the problem; their order is semantically significant for elimination, coordinates, or returned tuples.
- `bounds` — Exact coordinate bounds restricting search or construction.
- `method` — Computation-method selector; automatic mode chooses an applicable certified route.
- `precision` — Presentation/numerical precision used only where approximation is explicitly part of the API.
- `measure_dimension` — Dimension of the measure used for moments/centroids; useful for intrinsic lower-dimensional sets.
- `return_result` — When true, return the structured result/certificate/diagnostic object instead of only the convenience value.
- `parameters` — Symbols treated as free parameters rather than eliminated decision variables.
- `return_stratified` — Request parameter-space strata rather than only their combined condition.

**Result semantics.** The region is first reduced to explicit iterated-integral pieces using ``reduce_region_integral``. The ``method`` option controls evaluation:

**Representative regression-backed call.**

```python
integrate_over_region(1, para, [x, y])
```

This call is exercised by the regression contract “polygon tetrahedron and parallelogram integrals”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## reduce_region_integral

`reduce_region_integral(integrand: 'object', condition: 'object', variables: 'Sequence[sp.Symbol | str]', *, bounds: 'Sequence[tuple[sp.Symbol | str, object, object]] | Mapping[sp.Symbol | str, tuple[object, object]] | None' = None, return_integrals: 'bool' = False, parameters: 'Sequence[sp.Symbol | str] | None' = None) -> 'ReducedRegionIntegral | tuple[sp.Integral, ...]'`

Reduce a supported region integral to explicit iterated integrals.

**Parameter semantics.**

- `integrand` — Exact symbolic integrand.
- `condition` — Semialgebraic condition restricting the computation.
- `variables` — Ordered real variables of the problem; their order is semantically significant for elimination, coordinates, or returned tuples.
- `bounds` — Exact coordinate bounds restricting search or construction.
- `return_integrals` — Return the underlying exact integrals in addition to the reduced summary.
- `parameters` — Symbols treated as free parameters rather than eliminated decision variables.

**Result semantics.** This is the structural layer used by ``integrate_over_region``. It does not call ``sympy.integrate`` unless callers later ask to evaluate the

**Representative regression-backed call.**

```python
reduce_region_integral(x + y, cond, [x, y])
```

This call is exercised by the regression contract “cylindrical bounds handle disjoint union”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## quantifier_eliminate

`quantifier_eliminate(formula: 'object', quantifiers: 'Sequence[tuple[str, sp.Symbol | str]] | None' = None, *, variables: 'Sequence[sp.Symbol | str] | None' = None, strategy: 'str' = 'auto', return_result: 'bool' = False) -> 'sp.Expr | QuantifierEliminationResult'`

Eliminate real quantifiers with a certified specialist-first dispatcher.

**Parameter semantics.**

- `formula` — Boolean/relational semialgebraic formula defining the set or proposition.
- `quantifiers` — Ordered quantifier prefix, from outermost to innermost.
- `variables` — Ordered real variables of the problem; their order is semantically significant for elimination, coordinates, or returned tuples.
- `strategy` — Algorithm-selection policy. Use the default automatic policy unless a specific exact backend is being tested or diagnosed.
- `return_result` — When true, return the structured result/certificate/diagnostic object instead of only the convenience value.

**Result semantics.** ``formula`` may contain semialg :class:`Exists`/:class:`ForAll` nodes or a quantifier-free matrix accompanied by an explicit prenex ``quantifiers``

**Representative regression-backed call.**

```python
quantifier_eliminate(ForAll(x, x**2 + a >= 0))
```

This call is exercised by the regression contract “universal quantifier”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## is_real_valued

`is_real_valued(expr: 'sp.Expr', variables: 'Sequence[sp.Symbol | str] | None' = None, *, assumptions: 'FormulaLike' = True) -> 'bool'`

Return whether supported domain conditions follow from assumptions.

**Parameter semantics.**

- `expr` — Symbolic expression being analyzed.
- `variables` — Ordered real variables of the problem; their order is semantically significant for elimination, coordinates, or returned tuples.
- `assumptions` — Additional real-domain conditions under which the requested statement or computation is interpreted.

**Result semantics.** The return contract is `bool`; see the family reference for structured-result details and certification status.

**Representative regression-backed call.**

```python
is_real_valued(sp.sqrt(x - 1), [x], assumptions=x < 1)
```

This call is exercised by the regression contract “function domain and real valuedness helpers”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## simplify_system

`simplify_system(constraints: 'FormulaLike | Iterable[FormulaLike]', variables: 'Sequence[sp.Symbol | str] | None' = None, *, assumptions: 'FormulaLike | Iterable[FormulaLike]' = True, return_result: 'bool' = False, strategy: 'str | None' = None, eliminate_equalities: 'bool' = False, output: 'str' = 'formula') -> 'sp.Expr | tuple[sp.Expr, ...] | SimplifiedSystem'`

Simplify a real semialgebraic system with CAD/QE-backed checks.

**Parameter semantics.**

- `constraints` — Polynomial equalities/inequalities or a Boolean formula defining the feasible set.
- `variables` — Ordered real variables of the problem; their order is semantically significant for elimination, coordinates, or returned tuples.
- `assumptions` — Additional real-domain conditions under which the requested statement or computation is interpreted.
- `return_result` — When true, return the structured result/certificate/diagnostic object instead of only the convenience value.
- `strategy` — Algorithm-selection policy. Use the default automatic policy unless a specific exact backend is being tested or diagnosed.
- `eliminate_equalities` — Whether equality elimination is attempted during simplification.
- `output` — Requested public representation/view of the computed result.

**Result semantics.** The routine focuses on dependable semantic simplifications: contradiction detection, duplicate removal, redundancy removal by implication,

**Representative regression-backed call.**

```python
simplify_system([x > 0, x >= 0], [x])
```

This call is exercised by the regression contract “simplify system removes redundant constraints”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## simplify_boole

`simplify_boole(expr: 'FormulaLike | Iterable[FormulaLike]', variables: 'Sequence[sp.Symbol | str] | None' = None, *, assumptions: 'FormulaLike | Iterable[FormulaLike] | None' = None, strategy: 'str | None' = None, form: 'str' = 'auto', semantic: 'bool' = True, return_result: 'bool' = False) -> 'sp.Expr | BooleanSimplificationResult'`

Simplify a semialgebraic Boolean formula over the real numbers.

**Parameter semantics.**

- `expr` — Symbolic expression being analyzed.
- `variables` — Ordered real variables of the problem; their order is semantically significant for elimination, coordinates, or returned tuples.
- `assumptions` — Additional real-domain conditions under which the requested statement or computation is interpreted.
- `strategy` — Algorithm-selection policy. Use the default automatic policy unless a specific exact backend is being tested or diagnosed.
- `form` — Requested output/formula form.
- `semantic` — Whether simplification should use semantic equivalence rather than syntax alone.
- `return_result` — When true, return the structured result/certificate/diagnostic object instead of only the convenience value.

**Result semantics.** This is a conservative public simplifier. It first performs ordinary SymPy Boolean simplification, then uses the CAD/QE-backed decision wrappers to

**Representative regression-backed call.**

```python
simplify_boole(expr, [x, y])
```

This call is exercised by the regression contract “simplify boole removes redundant disjunct”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## simplify_piecewise

`simplify_piecewise(expr: 'sp.Expr', variables: 'Sequence[sp.Symbol | str] | None' = None, *, assumptions: 'FormulaLike | Iterable[FormulaLike] | None' = None, strategy: 'str | None' = None, return_result: 'bool' = False) -> 'sp.Expr | PiecewiseSimplificationResult'`

Simplify a Piecewise expression using semialgebraic branch conditions.

**Parameter semantics.**

- `expr` — Symbolic expression being analyzed.
- `variables` — Ordered real variables of the problem; their order is semantically significant for elimination, coordinates, or returned tuples.
- `assumptions` — Additional real-domain conditions under which the requested statement or computation is interpreted.
- `strategy` — Algorithm-selection policy. Use the default automatic policy unless a specific exact backend is being tested or diagnosed.
- `return_result` — When true, return the structured result/certificate/diagnostic object instead of only the convenience value.

**Result semantics.** Branches whose effective condition is unsatisfiable are removed. A branch becomes unconditional when its effective condition is a tautology. Adjacent

**Representative regression-backed call.**

```python
simplify_piecewise(expr, [x])
```

This call is exercised by the regression contract “simplify piecewise removes impossible branch”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## convert_region

`convert_region(region, representation='canonical', *, variables=None, return_result=False, bounds=None, dimension=None, require_verified=True, slices=4, tolerance=1e-09, precision=30, require_conforming=True)`

Convert a region to a certified exact or explicitly lossy representation.

**Parameter semantics.**

- `region` — Semialgebraic or standard region on which the operation is performed.
- `representation` — Requested geometric/algebraic representation.
- `variables` — Ordered real variables of the problem; their order is semantically significant for elimination, coordinates, or returned tuples.
- `return_result` — When true, return the structured result/certificate/diagnostic object instead of only the convenience value.
- `bounds` — Exact coordinate bounds restricting search or construction.
- `dimension` — Requested intrinsic/ambient dimension where the constructor or analysis needs it.
- `require_verified` — Require verification rather than accepting an unverified candidate.
- `slices` — Slice specification for decomposition/construction.
- `tolerance` — Numerical presentation tolerance where an API explicitly allows approximate geometry.
- `precision` — Presentation/numerical precision used only where approximation is explicitly part of the API.
- `require_conforming` — Require a conforming decomposition/mesh.

**Result semantics.** The return contract is `the documented result type for this operation`; see the family reference for structured-result details and certification status.

**Representative regression-backed call.**

```python
convert_region(p, "formula")
```

This call is exercised by the regression contract “p4 conversion contract”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## polygonal_region_from_paths

`polygonal_region_from_paths(paths, *, fill_rule='nonzero')`

Construct the exact planar filled set selected by a winding fill rule.

**Parameter semantics.**

- `paths` — Closed polygonal paths interpreted according to the selected fill rule.
- `fill_rule` — Rule used to interpret nested/self-overlapping polygonal paths.

**Result semantics.** The return contract is `the documented result type for this operation`; see the family reference for structured-result details and certification status.

**Representative regression-backed call.**

```python
polygonal_region_from_paths([bow], fill_rule="nonzero")
```

This call is exercised by the regression contract “p3 crossing and fill rules”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## triangulate_polytope

`triangulate_polytope(polytope: 'Polytope | Sequence[Sequence[object]]', *, strategy: 'str' = 'pulling') -> 'PolytopeDecomposition'`

Triangulate a full-dimensional convex polytope exactly and deterministically.

**Parameter semantics.**

- `polytope` — Polytope to decompose/triangulate.
- `strategy` — Algorithm-selection policy. Use the default automatic policy unless a specific exact backend is being tested or diagnosed.

**Result semantics.** The return contract is `PolytopeDecomposition`; see the family reference for structured-result details and certification status.

**Representative regression-backed call.**

```python
triangulate_polytope(p)
```

This call is exercised by the regression contract “simplex pulling is one simplex”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## decompose_polytope

`decompose_polytope(polytope, *, strategy: 'str' = 'pulling', target: 'str' = 'simplices')`

Decompose a convex polytope into an exact requested cell representation.

**Parameter semantics.**

- `polytope` — Polytope to decompose/triangulate.
- `strategy` — Algorithm-selection policy. Use the default automatic policy unless a specific exact backend is being tested or diagnosed.
- `target` — Requested target object or result; for certificate APIs this is the relation to prove, while decomposition APIs use it to select the desired cell representation.

**Result semantics.** The return contract is `the documented result type for this operation`; see the family reference for structured-result details and certification status.

**Representative regression-backed call.**

```python
sa.decompose_polytope(poly)
```

This call is exercised by the regression contract “polytope decomposition preserves triangle area”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## tetrahedralize_cell

`tetrahedralize_cell(vertices: 'Sequence[Sequence[object]]', *, vertex_ids: 'Sequence[int] | None' = None) -> 'MixedCellTetrahedralization'`

Tetrahedralize one convex 3-cell using a global vertex ordering.

**Parameter semantics.**

- `vertices` — Vertices defining the polygon/polytope.
- `vertex_ids` — Indices identifying vertices in a combinatorial representation.

**Result semantics.** The return contract is `MixedCellTetrahedralization`; see the family reference for structured-result details and certification status.

**Representative regression-backed call.**

```python
tetrahedralize_cell(cube())
```

This call is exercised by the regression contract “p6 cell counts and shared face conformity”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## tetrahedralize_cells

`tetrahedralize_cells(vertices: 'Sequence[Sequence[object]]', cells: 'Sequence[Sequence[int]]') -> 'MixedMeshTetrahedralization'`

Tetrahedralize cells conformingly using shared global vertex identifiers.

**Parameter semantics.**

- `vertices` — Vertices defining the polygon/polytope.
- `cells` — Cells participating in the decomposition.

**Result semantics.** The return contract is `MixedMeshTetrahedralization`; see the family reference for structured-result details and certification status.

**Representative regression-backed call.**

```python
tetrahedralize_cells(pts, cells)
```

This call is exercised by the regression contract “adjacent cubes share same face triangulation”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## canonicalize_polygon

`canonicalize_polygon(region)`

Return the strongest certified canonical representation of polygonal geometry.

**Parameter semantics.**

- `region` — Semialgebraic or standard region on which the operation is performed.

**Result semantics.** The return contract is `the documented result type for this operation`; see the family reference for structured-result details and certification status.

**Representative regression-backed call.**

```python
sa.canonicalize_polygon(p)
```

This call is exercised by the regression contract “canonicalize polygon is invariant under cyclic vertex rotation”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## canonicalize_polyhedron

`canonicalize_polyhedron(region)`

Return a certified canonical polyhedral representation without unsafe convexification.

**Parameter semantics.**

- `region` — Semialgebraic or standard region on which the operation is performed.

**Result semantics.** The return contract is `the documented result type for this operation`; see the family reference for structured-result details and certification status.

**Representative regression-backed call.**

```python
canonicalize_polyhedron(tc)
```

This call is exercised by the regression contract “p2 nonconvex shell preserved”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## canonicalize_region

`canonicalize_region(region)`

Canonicalize a supported region using exact structural recognition.

**Parameter semantics.**

- `region` — Semialgebraic or standard region on which the operation is performed.

**Result semantics.** The return contract is `the documented result type for this operation`; see the family reference for structured-result details and certification status.

**Representative regression-backed call.**

```python
sa.canonicalize_region(region)
```

This call is exercised by the regression contract “canonicalize region is idempotent on polygon”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## verify_nonnegative_combination_certificate

`verify_nonnegative_combination_certificate(cert, variables) -> 'bool'`

Replay an exact nonnegative-combination and optional SOS certificate.

**Parameter semantics.**

- `cert` — Certificate object to verify.
- `variables` — Ordered real variables of the problem; their order is semantically significant for elimination, coordinates, or returned tuples.

**Result semantics.** The verifier is a trust boundary: malformed typed payloads are rejected cleanly before polynomial replay begins.

**Representative regression-backed call.**

```python
verify_nonnegative_combination_certificate(c, (x, y))
```

This call is exercised by the regression contract “nonnegative combination”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## nonnegative_combination_certificate

`nonnegative_combination_certificate(target, premises, variables, *, allow_sos=True, sos_certificate=None)`

Find an exact certificate target = sum(lambda_i premise_i) + SOS, lambda_i >= 0.

**Parameter semantics.**

- `target` — Requested target object or result; for certificate APIs this is the relation to prove, while decomposition APIs use it to select the desired cell representation.
- `premises` — Premise inequalities/equalities from which the target relation is to be certified.
- `variables` — Ordered real variables of the problem; their order is semantically significant for elimination, coordinates, or returned tuples.
- `allow_sos` — Permit SOS-based certification as one of the exact-verified proof routes.
- `sos_certificate` — Optional supplied SOS certificate to verify/use.

**Result semantics.** The return contract is `the documented result type for this operation`; see the family reference for structured-result details and certification status.

**Representative regression-backed call.**

```python
semialg.nonnegative_combination_certificate(2 * x, (x,), (x,), allow_sos=False)
```

This call is exercised by the regression contract “semialgebraic constraint secondary calls”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## implied_polynomial_inequality

`implied_polynomial_inequality(region, inequality, variables=None, *, use_certificates=True)`

Certify that ``region`` implies a polynomial inequality.

**Parameter semantics.**

- `region` — Semialgebraic or standard region on which the operation is performed.
- `inequality` — Polynomial inequality under analysis.
- `variables` — Ordered real variables of the problem; their order is semantically significant for elimination, coordinates, or returned tuples.
- `use_certificates` — Enable certificate-producing/verification paths when supported.

**Result semantics.** The return contract is `the documented result type for this operation`; see the family reference for structured-result details and certification status.

**Representative regression-backed call.**

```python
sa.implied_polynomial_inequality(region, 2 * X >= 0, (X,))
```

This call is exercised by the regression contract “implied and redundant polynomial inequality agree”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## redundant_polynomial_inequalities

`redundant_polynomial_inequalities(region, variables=None)`

Return inequalities implied by the other constraints in their DNF clause.

**Parameter semantics.**

- `region` — Semialgebraic or standard region on which the operation is performed.
- `variables` — Ordered real variables of the problem; their order is semantically significant for elimination, coordinates, or returned tuples.

**Result semantics.** The return contract is `the documented result type for this operation`; see the family reference for structured-result details and certification status.

**Representative regression-backed call.**

```python
sa.redundant_polynomial_inequalities(region, (X,))
```

This call is exercised by the regression contract “implied and redundant polynomial inequality agree”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## component_constraint_descriptions

`component_constraint_descriptions(region, variables=None, *, max_pieces=None)`

Describe the real model separately on each certified irreducible closure component.

**Parameter semantics.**

- `region` — Semialgebraic or standard region on which the operation is performed.
- `variables` — Ordered real variables of the problem; their order is semantically significant for elimination, coordinates, or returned tuples.
- `max_pieces` — Bound on decomposition pieces retained during the computation.

**Result semantics.** The return contract is `the documented result type for this operation`; see the family reference for structured-result details and certification status.

**Representative regression-backed call.**

```python
sa.component_constraint_descriptions(sp.Eq(X * Y, 0), (X, Y))
```

This call is exercised by the regression contract “component constraint descriptions cover each axis component”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## polynomial_constraints

`polynomial_constraints(region, variables: 'Sequence[sp.Symbol | str] | None' = None) -> 'PolynomialConstraintSystem'`

Return a DNF-preserving structured polynomial constraint description.

**Parameter semantics.**

- `region` — Semialgebraic or standard region on which the operation is performed.
- `variables` — Ordered real variables of the problem; their order is semantically significant for elimination, coordinates, or returned tuples.

**Result semantics.** Boolean disjunction is represented by separate clauses; conjunction is represented inside a clause. Negations are pushed to relational atoms.

**Representative regression-backed call.**

```python
semialg.polynomial_constraints(region, (x, y))
```

This call is exercised by the regression contract “constraints and active boundary”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## active_constraints

`active_constraints(region, point, variables=None) -> 'ActiveConstraintResult'`

Return polynomial constraints active at a feasible exact point.

**Parameter semantics.**

- `region` — Semialgebraic or standard region on which the operation is performed.
- `point` — Point in the ambient coordinate system.
- `variables` — Ordered real variables of the problem; their order is semantically significant for elimination, coordinates, or returned tuples.

**Result semantics.** The return contract is `ActiveConstraintResult`; see the family reference for structured-result details and certification status.

**Representative regression-backed call.**

```python
sa.active_constraints(region, (0, 0), (X, Y))
```

This call is exercised by the regression contract “active constraints change exactly at triangle boundary”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## relative_interior

`relative_interior(region, variables=None) -> 'sp.Expr'`

Return the interior relative to the certified affine hull of ``region``.

**Parameter semantics.**

- `region` — Semialgebraic or standard region on which the operation is performed.
- `variables` — Ordered real variables of the problem; their order is semantically significant for elimination, coordinates, or returned tuples.

**Result semantics.** The return contract is `sp.Expr`; see the family reference for structured-result details and certification status.

**Representative regression-backed call.**

```python
sa.relative_interior(segment, (X, Y))
```

This call is exercised by the regression contract “relative boundary and interior partition closed segment in affine hull”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## relative_boundary

`relative_boundary(region, variables=None) -> 'sp.Expr'`

Return the boundary relative to the certified affine hull of ``region``.

**Parameter semantics.**

- `region` — Semialgebraic or standard region on which the operation is performed.
- `variables` — Ordered real variables of the problem; their order is semantically significant for elimination, coordinates, or returned tuples.

**Result semantics.** The return contract is `sp.Expr`; see the family reference for structured-result details and certification status.

**Representative regression-backed call.**

```python
sa.relative_boundary(segment, (X, Y))
```

This call is exercised by the regression contract “relative boundary and interior partition closed segment in affine hull”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## semialgebraic_tangent_cone

`semialgebraic_tangent_cone(region, point, variables=None) -> 'sp.Expr'`

Return the exact Bouligand tangent cone at a point of a semialgebraic set.

**Parameter semantics.**

- `region` — Semialgebraic or standard region on which the operation is performed.
- `point` — Point in the ambient coordinate system.
- `variables` — Ordered real variables of the problem; their order is semantically significant for elimination, coordinates, or returned tuples.

**Result semantics.** The cone is the closure of secant velocities ``(x-p)/t`` for ``x`` in the set and ``t > 0``. Quantifier elimination makes the construction exact.

**Representative regression-backed call.**

```python
sa.semialgebraic_tangent_cone(X >= 0, (0,), (X,))
```

This call is exercised by the regression contract “semialgebraic tangent cone of halfline is itself”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## parameterization_geometry

`parameterization_geometry(mapping, domain=True, parameters: 'Sequence[sp.Symbol | str] | None' = None, *, image_variables=None) -> 'ParameterizationGeometry'`

Analyze generic rank, rank-drop locus, image, and generic fiber dimension.

**Parameter semantics.**

- `mapping` — Polynomial or affine map whose image/local geometry is analyzed.
- `domain` — Mathematical domain selector; semialg public decision procedures are principally exact over the reals.
- `parameters` — Symbols treated as free parameters rather than eliminated decision variables.
- `image_variables` — Coordinate symbols used for the image space.

**Result semantics.** The return contract is `ParameterizationGeometry`; see the family reference for structured-result details and certification status.

**Representative regression-backed call.**

```python
semialg.parameterization_geometry((u,), sp.Eq(v, 0), (u, v))
```

This call is exercised by the regression contract “lower dimensional domain uses intrinsic tangent rank”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## parameterization_critical_locus

`parameterization_critical_locus(mapping, domain=True, parameters=None) -> 'sp.Expr'`

Return the exact generic-rank-drop locus of a parameterization.

**Parameter semantics.**

- `mapping` — Polynomial or affine map whose image/local geometry is analyzed.
- `domain` — Mathematical domain selector; semialg public decision procedures are principally exact over the reals.
- `parameters` — Symbols treated as free parameters rather than eliminated decision variables.

**Result semantics.** The return contract is `sp.Expr`; see the family reference for structured-result details and certification status.

**Representative regression-backed call.**

```python
sa.parameterization_critical_locus((t**2, t**3), sp.true, (t,))
```

This call is exercised by the regression contract “parameterization critical locus and values for parabola”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## parameterization_critical_values

`parameterization_critical_values(mapping, domain=True, parameters=None, *, image_variables=None) -> 'sp.Expr'`

Return the exact image of the parameterization's critical locus.

**Parameter semantics.**

- `mapping` — Polynomial or affine map whose image/local geometry is analyzed.
- `domain` — Mathematical domain selector; semialg public decision procedures are principally exact over the reals.
- `parameters` — Symbols treated as free parameters rather than eliminated decision variables.
- `image_variables` — Coordinate symbols used for the image space.

**Result semantics.** The return contract is `sp.Expr`; see the family reference for structured-result details and certification status.

**Representative regression-backed call.**

```python
semialg.parameterization_critical_values((t**2, t**3), sp.true, (t,))
```

This call is exercised by the regression contract “parameterization and qe secondary calls”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## is_singular

`is_singular(equations, point, variables=None, *, codimension: 'int | None' = None) -> 'bool'`

Return whether ``point`` is singular on the polynomial variety.

**Parameter semantics.**

- `equations` — Polynomial equations defining the algebraic set.
- `point` — Point in the ambient coordinate system.
- `variables` — Ordered real variables of the problem; their order is semantically significant for elimination, coordinates, or returned tuples.
- `codimension` — Requested codimension used to select strata/components.

**Result semantics.** The return contract is `bool`; see the family reference for structured-result details and certification status.

**Representative regression-backed call.**

```python
is_singular((cusp,), (0, 0), (x, y))
```

This call is exercised by the regression contract “singularity tangent and thom contracts”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## is_smooth

`is_smooth(equations, variables=None, *, codimension: 'int | None' = None) -> 'bool'`

Return whether the real polynomial variety has empty singular locus.

**Parameter semantics.**

- `equations` — Polynomial equations defining the algebraic set.
- `variables` — Ordered real variables of the problem; their order is semantically significant for elimination, coordinates, or returned tuples.
- `codimension` — Requested codimension used to select strata/components.

**Result semantics.** The return contract is `bool`; see the family reference for structured-result details and certification status.

**Representative regression-backed call.**

```python
is_smooth([cusp], [x, y])
```

This call is exercised by the regression contract “smoothness and point singularity conveniences”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## singular_locus

`singular_locus(equations, variables=None, *, codimension: 'int | None' = None, reduce: 'bool' = True, max_pieces: 'int | None' = None) -> 'sp.Expr'`

Return equations defining the singular locus of an algebraic variety.

**Parameter semantics.**

- `equations` — Polynomial equations defining the algebraic set.
- `variables` — Ordered real variables of the problem; their order is semantically significant for elimination, coordinates, or returned tuples.
- `codimension` — Requested codimension used to select strata/components.
- `reduce` — Whether to reduce/simplify the constructed result.
- `max_pieces` — Bound on decomposition pieces retained during the computation.

**Result semantics.** The Jacobian criterion is applied to a certified radical presentation by default, so nilpotent multiplicities such as ``x**2 = 0`` cannot create

**Representative regression-backed call.**

```python
cusp.singular_locus()
```

This call is exercised by the regression contract “cusp detects origin and local dimension”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## tangent_dimension

`tangent_dimension(equations, point, variables=None) -> 'int'`

Return the exact Zariski tangent-space dimension at ``point``.

**Parameter semantics.**

- `equations` — Polynomial equations defining the algebraic set.
- `point` — Point in the ambient coordinate system.
- `variables` — Ordered real variables of the problem; their order is semantically significant for elimination, coordinates, or returned tuples.

**Result semantics.** The return contract is `int`; see the family reference for structured-result details and certification status.

**Representative regression-backed call.**

```python
tangent_dimension((cusp,), (0, 0), (x, y))
```

This call is exercised by the regression contract “singularity tangent and thom contracts”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## as_cad_region

`as_cad_region(region: 'object', variables: 'Sequence[sp.Symbol | str] | None' = None, *, strategy: 'str' = 'auto') -> 'CADRegion'`

Coerce a region/formula/CAD result to a reusable :class:`CADRegion`.

**Parameter semantics.**

- `region` — Semialgebraic or standard region on which the operation is performed.
- `variables` — Ordered real variables of the problem; their order is semantically significant for elimination, coordinates, or returned tuples.
- `strategy` — Algorithm-selection policy. Use the default automatic policy unless a specific exact backend is being tested or diagnosed.

**Result semantics.** The return contract is `CADRegion`; see the family reference for structured-result details and certification status.

**Representative regression-backed call.**

```python
as_cad_region(box)
```

This call is exercised by the regression contract “exact region measure and integration hooks”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## analyze_affine_map

`analyze_affine_map(mapping: 'Sequence[object] | sp.MatrixBase', variables: 'Sequence[sp.Symbol | str] | None' = None, *, offset: 'Sequence[object] | None' = None) -> 'AffineMapAnalysis'`

Analyze an affine expression map or an explicit matrix/offset pair.

**Parameter semantics.**

- `mapping` — Polynomial or affine map whose image/local geometry is analyzed.
- `variables` — Ordered real variables of the problem; their order is semantically significant for elimination, coordinates, or returned tuples.
- `offset` — Offset applied to a generated regular object.

**Result semantics.** With ``variables`` supplied, ``mapping`` is interpreted as a vector of affine expressions.  Without variables it is interpreted as a matrix and

**Representative regression-backed call.**

```python
analyze_affine_map([[sp.I]])
```

This call is exercised by the regression contract “affine maps reject nonreal data at boundary”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## convexity_certificate

`convexity_certificate(region, variables: 'Sequence[sp.Symbol | str] | None' = None, *, strategy: 'str | None' = None) -> 'ConvexityCertificate'`

Decide semialgebraic set convexity through an exact staged hierarchy.

**Parameter semantics.**

- `region` — Semialgebraic or standard region on which the operation is performed.
- `variables` — Ordered real variables of the problem; their order is semantically significant for elimination, coordinates, or returned tuples.
- `strategy` — Algorithm-selection policy. Use the default automatic policy unless a specific exact backend is being tested or diagnosed.

**Result semantics.** The stages are formula normalization, trivial/empty/singleton handling, complete one-dimensional classification, affine/polyhedral recognition,

**Representative regression-backed call.**

```python
convexity_certificate(region, (x, y))
```

This call is exercised by the regression contract “segment counterexample certificate replays from returned witness”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## is_convex

`is_convex(region, variables: 'Sequence[sp.Symbol | str] | None' = None, *, strategy: 'str | None' = None) -> 'bool'`

Return whether a semialgebraic region is convex, exactly.

**Parameter semantics.**

- `region` — Semialgebraic or standard region on which the operation is performed.
- `variables` — Ordered real variables of the problem; their order is semantically significant for elimination, coordinates, or returned tuples.
- `strategy` — Algorithm-selection policy. Use the default automatic policy unless a specific exact backend is being tested or diagnosed.

**Result semantics.** The return contract is `bool`; see the family reference for structured-result details and certification status.

**Representative regression-backed call.**

```python
is_convex(region, (x,))
```

This call is exercised by the regression contract “generated one dimensional interval rewrites remain convex”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## argmax_set

`argmax_set(expression, region=True, variables=None) -> 'sp.Expr'`

Return the exact global maximizer set as a semialgebraic formula.

**Parameter semantics.**

- `expression` — Symbolic expression being analyzed.
- `region` — Semialgebraic or standard region on which the operation is performed.
- `variables` — Ordered real variables of the problem; their order is semantically significant for elimination, coordinates, or returned tuples.

**Result semantics.** The return contract is `sp.Expr`; see the family reference for structured-result details and certification status.

**Representative regression-backed call.**

```python
argmax_set(x, interval, (x,))
```

This call is exercised by the regression contract “interval extrema and level sets have exact endpoint semantics”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## closest_points

`closest_points(left, right, variables=None)`

Return all exact closest point pairs when the distance is attained.

**Parameter semantics.**

- `left` — Left operand/region in the binary operation.
- `right` — Right operand/region in the binary operation.
- `variables` — Ordered real variables of the problem; their order is semantically significant for elimination, coordinates, or returned tuples.

**Result semantics.** The return contract is `the documented result type for this operation`; see the family reference for structured-result details and certification status.

**Representative regression-backed call.**

```python
closest_points(x <= 0, x >= 2, [x])
```

This call is exercised by the regression contract “coordinate range diameter nearest and support operations”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## coordinate_range

`coordinate_range(region, variable, variables=None, *, return_result: 'bool' = False)`

Return the exact range of one coordinate over a region.

**Parameter semantics.**

- `region` — Semialgebraic or standard region on which the operation is performed.
- `variable` — Single variable singled out for range/fiber/elimination analysis.
- `variables` — Ordered real variables of the problem; their order is semantically significant for elimination, coordinates, or returned tuples.
- `return_result` — When true, return the structured result/certificate/diagnostic object instead of only the convenience value.

**Result semantics.** The return contract is `the documented result type for this operation`; see the family reference for structured-result details and certification status.

**Representative regression-backed call.**

```python
coordinate_range(region, x, [x])
```

This call is exercised by the regression contract “coordinate range agrees with membership”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## diameter

`diameter(region, variables=None, *, return_result: 'bool' = False)`

Return the exact Euclidean diameter (supremal pairwise distance).

**Parameter semantics.**

- `region` — Semialgebraic or standard region on which the operation is performed.
- `variables` — Ordered real variables of the problem; their order is semantically significant for elimination, coordinates, or returned tuples.
- `return_result` — When true, return the structured result/certificate/diagnostic object instead of only the convenience value.

**Result semantics.** The return contract is `the documented result type for this operation`; see the family reference for structured-result details and certification status.

**Representative regression-backed call.**

```python
diameter(moved, [x])
```

This call is exercised by the regression contract “translation preserves diameter and shifts centroid”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## distance_set

`distance_set(left, right=None, variables=None, *, distance_symbol=None) -> 'sp.Expr'`

Return the exact set of Euclidean pairwise distances as a formula.

**Parameter semantics.**

- `left` — Left operand/region in the binary operation.
- `right` — Right operand/region in the binary operation.
- `variables` — Ordered real variables of the problem; their order is semantically significant for elimination, coordinates, or returned tuples.
- `distance_symbol` — Symbol used to represent an exact distance value in elimination formulas.

**Result semantics.** The return contract is `sp.Expr`; see the family reference for structured-result details and certification status.

**Representative regression-backed call.**

```python
distance_set(left_point, right_point, (x,))
```

This call is exercised by the regression contract “interval geometry queries cross check one another”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## extrema_set

`extrema_set(expression, region=True, variables=None) -> 'sp.Expr'`

Return the union of the exact global minimum and maximum sets.

**Parameter semantics.**

- `expression` — Symbolic expression being analyzed.
- `region` — Semialgebraic or standard region on which the operation is performed.
- `variables` — Ordered real variables of the problem; their order is semantically significant for elimination, coordinates, or returned tuples.

**Result semantics.** The return contract is `sp.Expr`; see the family reference for structured-result details and certification status.

**Representative regression-backed call.**

```python
extrema_set(x, interval, (x,))
```

This call is exercised by the regression contract “interval extrema and level sets have exact endpoint semantics”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## has_empty_interior

`has_empty_interior(region, variables=None) -> 'bool'`

Return whether the region has empty ambient interior.

**Parameter semantics.**

- `region` — Semialgebraic or standard region on which the operation is performed.
- `variables` — Ordered real variables of the problem; their order is semantically significant for elimination, coordinates, or returned tuples.

**Result semantics.** The return contract is `bool`; see the family reference for structured-result details and certification status.

**Representative regression-backed call.**

```python
has_empty_interior(point, (x,))
```

This call is exercised by the regression contract “dimension boundedness density and path contracts”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## inertia_tensor

`inertia_tensor(region, variables, **kwargs) -> 'sp.Matrix'`

Return the unit-density second moment-of-inertia tensor about the origin.

**Parameter semantics.**

- `region` — Semialgebraic or standard region on which the operation is performed.
- `variables` — Ordered real variables of the problem; their order is semantically significant for elimination, coordinates, or returned tuples.
- `kwargs` — Kwargs controlling this operation; see the signature type/default for the accepted representation.

**Result semantics.** The return contract is `sp.Matrix`; see the family reference for structured-result details and certification status.

**Representative regression-backed call.**

```python
inertia_tensor(square, [x, y])
```

This call is exercised by the regression contract “moment covariance identity on centered square”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## intersects

`intersects(left, right, variables=None, *, strategy: 'str | None' = None) -> 'bool'`

Return whether two semialgebraic regions have nonempty intersection.

**Parameter semantics.**

- `left` — Left operand/region in the binary operation.
- `right` — Right operand/region in the binary operation.
- `variables` — Ordered real variables of the problem; their order is semantically significant for elimination, coordinates, or returned tuples.
- `strategy` — Algorithm-selection policy. Use the default automatic policy unless a specific exact backend is being tested or diagnosed.

**Result semantics.** The return contract is `bool`; see the family reference for structured-result details and certification status.

**Representative regression-backed call.**

```python
intersects(left, right, [x])
```

This call is exercised by the regression contract “disjointness and intersection are complements”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## is_interior_disjoint

`is_interior_disjoint(left, right, variables=None, *, strategy: 'str | None' = None) -> 'bool'`

Return whether two regions have disjoint ambient interiors.

**Parameter semantics.**

- `left` — Left operand/region in the binary operation.
- `right` — Right operand/region in the binary operation.
- `variables` — Ordered real variables of the problem; their order is semantically significant for elimination, coordinates, or returned tuples.
- `strategy` — Algorithm-selection policy. Use the default automatic policy unless a specific exact backend is being tested or diagnosed.

**Result semantics.** Boundary contact is permitted: two closed regions that meet only along their boundaries are interior-disjoint.

**Representative regression-backed call.**

```python
semialg.is_interior_disjoint(left, right, (x,))
```

This call is exercised by the regression contract “new region boolean and topological relations are directly callable”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## is_bounded

`is_bounded(region, variables=None, *, strategy: 'str | None' = None) -> 'bool'`

Return whether the semialgebraic region is bounded.

**Parameter semantics.**

- `region` — Semialgebraic or standard region on which the operation is performed.
- `variables` — Ordered real variables of the problem; their order is semantically significant for elimination, coordinates, or returned tuples.
- `strategy` — Algorithm-selection policy. Use the default automatic policy unless a specific exact backend is being tested or diagnosed.

**Result semantics.** The return contract is `bool`; see the family reference for structured-result details and certification status.

**Representative regression-backed call.**

```python
region.is_bounded()
```

This call is exercised by the regression contract “known region basic invariants”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## is_closed

`is_closed(region, variables=None, *, strategy: 'str | None' = None) -> 'bool'`

Return whether the semialgebraic region is closed.

**Parameter semantics.**

- `region` — Semialgebraic or standard region on which the operation is performed.
- `variables` — Ordered real variables of the problem; their order is semantically significant for elimination, coordinates, or returned tuples.
- `strategy` — Algorithm-selection policy. Use the default automatic policy unless a specific exact backend is being tested or diagnosed.

**Result semantics.** The return contract is `bool`; see the family reference for structured-result details and certification status.

**Representative regression-backed call.**

```python
region.is_closed()
```

This call is exercised by the regression contract “known region basic invariants”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## is_compact

`is_compact(region, variables=None, *, strategy: 'str | None' = None) -> 'bool'`

Return whether the semialgebraic region is compact.

**Parameter semantics.**

- `region` — Semialgebraic or standard region on which the operation is performed.
- `variables` — Ordered real variables of the problem; their order is semantically significant for elimination, coordinates, or returned tuples.
- `strategy` — Algorithm-selection policy. Use the default automatic policy unless a specific exact backend is being tested or diagnosed.

**Result semantics.** The return contract is `bool`; see the family reference for structured-result details and certification status.

**Representative regression-backed call.**

```python
is_compact(region, [x])
```

This call is exercised by the regression contract “compact interval remains compact under affine bijection”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## is_connected

`is_connected(region, variables=None) -> 'bool'`

Decide connectedness; for semialgebraic sets this equals path connectedness.

**Parameter semantics.**

- `region` — Semialgebraic or standard region on which the operation is performed.
- `variables` — Ordered real variables of the problem; their order is semantically significant for elimination, coordinates, or returned tuples.

**Result semantics.** The return contract is `bool`; see the family reference for structured-result details and certification status.

**Representative regression-backed call.**

```python
is_connected(closed, [x])
```

This call is exercised by the regression contract “geometric property and set relation api”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## is_dense_in

`is_dense_in(subset, ambient, variables=None, *, strategy: 'str | None' = None) -> 'bool'`

Return whether ``subset`` is dense in ``ambient`` in the ambient Euclidean topology.

**Parameter semantics.**

- `subset` — Candidate subset whose relation to the ambient region is tested.
- `ambient` — Ambient coordinate description for the constructed/analyzed object.
- `variables` — Ordered real variables of the problem; their order is semantically significant for elimination, coordinates, or returned tuples.
- `strategy` — Algorithm-selection policy. Use the default automatic policy unless a specific exact backend is being tested or diagnosed.

**Result semantics.** The return contract is `bool`; see the family reference for structured-result details and certification status.

**Representative regression-backed call.**

```python
is_dense_in(sp.Ne(x, 0), sp.true, [x])
```

This call is exercised by the regression contract “dense linear image and distance set operations”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## is_disjoint

`is_disjoint(left, right, variables=None, *, strategy: 'str | None' = None) -> 'bool'`

Return whether two semialgebraic regions are disjoint.

**Parameter semantics.**

- `left` — Left operand/region in the binary operation.
- `right` — Right operand/region in the binary operation.
- `variables` — Ordered real variables of the problem; their order is semantically significant for elimination, coordinates, or returned tuples.
- `strategy` — Algorithm-selection policy. Use the default automatic policy unless a specific exact backend is being tested or diagnosed.

**Result semantics.** The return contract is `bool`; see the family reference for structured-result details and certification status.

**Representative regression-backed call.**

```python
is_disjoint(a, x > hi, [x])
```

This call is exercised by the regression contract “generated set algebra identities”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## is_empty

`is_empty(region, variables=None, *, strategy: 'str | None' = None) -> 'bool'`

Return whether the semialgebraic region is empty.

**Parameter semantics.**

- `region` — Semialgebraic or standard region on which the operation is performed.
- `variables` — Ordered real variables of the problem; their order is semantically significant for elimination, coordinates, or returned tuples.
- `strategy` — Algorithm-selection policy. Use the default automatic policy unless a specific exact backend is being tested or diagnosed.

**Result semantics.** The return contract is `bool`; see the family reference for structured-result details and certification status.

**Representative regression-backed call.**

```python
is_empty(closed, [x])
```

This call is exercised by the regression contract “geometric property and set relation api”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## is_equal

`is_equal(left, right, variables=None, *, strategy: 'str | None' = None) -> 'bool'`

Return whether two semialgebraic regions define the same set.

**Parameter semantics.**

- `left` — Left operand/region in the binary operation.
- `right` — Right operand/region in the binary operation.
- `variables` — Ordered real variables of the problem; their order is semantically significant for elimination, coordinates, or returned tuples.
- `strategy` — Algorithm-selection policy. Use the default automatic policy unless a specific exact backend is being tested or diagnosed.

**Result semantics.** The return contract is `bool`; see the family reference for structured-result details and certification status.

**Representative regression-backed call.**

```python
is_equal(left, right, [x])
```

This call is exercised by the regression contract “set relation algebra identity cases”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## is_full_dimensional

`is_full_dimensional(region, variables=None) -> 'bool'`

Return whether the region has full dimension in its ambient variables.

**Parameter semantics.**

- `region` — Semialgebraic or standard region on which the operation is performed.
- `variables` — Ordered real variables of the problem; their order is semantically significant for elimination, coordinates, or returned tuples.

**Result semantics.** The return contract is `bool`; see the family reference for structured-result details and certification status.

**Representative regression-backed call.**

```python
is_full_dimensional(closed, [x])
```

This call is exercised by the regression contract “geometric property and set relation api”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## is_open

`is_open(region, variables=None, *, strategy: 'str | None' = None) -> 'bool'`

Return whether the semialgebraic region is open.

**Parameter semantics.**

- `region` — Semialgebraic or standard region on which the operation is performed.
- `variables` — Ordered real variables of the problem; their order is semantically significant for elimination, coordinates, or returned tuples.
- `strategy` — Algorithm-selection policy. Use the default automatic policy unless a specific exact backend is being tested or diagnosed.

**Result semantics.** The return contract is `bool`; see the family reference for structured-result details and certification status.

**Representative regression-backed call.**

```python
is_open(opened, [x])
```

This call is exercised by the regression contract “geometric property and set relation api”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## is_subset

`is_subset(left, right, variables=None, *, strategy: 'str | None' = None) -> 'bool'`

Return whether one semialgebraic region is contained in another.

**Parameter semantics.**

- `left` — Left operand/region in the binary operation.
- `right` — Right operand/region in the binary operation.
- `variables` — Ordered real variables of the problem; their order is semantically significant for elimination, coordinates, or returned tuples.
- `strategy` — Algorithm-selection policy. Use the default automatic policy unless a specific exact backend is being tested or diagnosed.

**Result semantics.** The return contract is `bool`; see the family reference for structured-result details and certification status.

**Representative regression-backed call.**

```python
is_subset(a, b, [x])
```

This call is exercised by the regression contract “generated set algebra identities”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## level_set

`level_set(expression, value=0, region=True) -> 'sp.Expr'`

Return ``region ∩ {expression = value}``.

**Parameter semantics.**

- `expression` — Symbolic expression being analyzed.
- `value` — Level/threshold value defining a level, sublevel, or superlevel set.
- `region` — Semialgebraic or standard region on which the operation is performed.

**Result semantics.** The return contract is `sp.Expr`; see the family reference for structured-result details and certification status.

**Representative regression-backed call.**

```python
level_set(x**2, 1, region)
```

This call is exercised by the regression contract “level set constructors preserve region constraint”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## linear_image

`linear_image(region, matrix, variables=None)`

Return the exact linear image ``A*x``.

**Parameter semantics.**

- `region` — Semialgebraic or standard region on which the operation is performed.
- `matrix` — Matrix defining the linear/affine operation or quadratic data.
- `variables` — Ordered real variables of the problem; their order is semantically significant for elimination, coordinates, or returned tuples.

**Result semantics.** This is the zero-offset specialization of :func:`semialg.affine_image` and uses the same representation-dispatching contract.

**Representative regression-backed call.**

```python
sa.linear_image(interval, [[3]], (X,))
```

This call is exercised by the regression contract “linear image commutes with positive scalar interval bounds”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## minkowski_sum

`minkowski_sum(left, right, variables=None)`

Return the exact Minkowski sum of two semialgebraic regions.

**Parameter semantics.**

- `left` — Left operand/region in the binary operation.
- `right` — Right operand/region in the binary operation.
- `variables` — Ordered real variables of the problem; their order is semantically significant for elimination, coordinates, or returned tuples.

**Result semantics.** Canonical points, intervals, boxes, and filled balls retain structural representations.  General cases use exact semialgebraic image elimination.

**Representative regression-backed call.**

```python
minkowski_sum(left, right)
```

This call is exercised by the regression contract “polytope minkowski sum from vertex sums”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## moment_matrix

`moment_matrix(region, variables, **kwargs) -> 'sp.Matrix'`

Return the normalized raw second-moment matrix ``E[x x.T]``.

**Parameter semantics.**

- `region` — Semialgebraic or standard region on which the operation is performed.
- `variables` — Ordered real variables of the problem; their order is semantically significant for elimination, coordinates, or returned tuples.
- `kwargs` — Kwargs controlling this operation; see the signature type/default for the accepted representation.

**Result semantics.** The return contract is `sp.Matrix`; see the family reference for structured-result details and certification status.

**Representative regression-backed call.**

```python
moment_matrix(square, [x, y])
```

This call is exercised by the regression contract “moment covariance identity on centered square”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## nearest_point

`nearest_point(point, region, variables=None)`

Return all exact nearest points when the distance is attained.

**Parameter semantics.**

- `point` — Point in the ambient coordinate system.
- `region` — Semialgebraic or standard region on which the operation is performed.
- `variables` — Ordered real variables of the problem; their order is semantically significant for elimination, coordinates, or returned tuples.

**Result semantics.** The return contract is `the documented result type for this operation`; see the family reference for structured-result details and certification status.

**Representative regression-backed call.**

```python
nearest_point((3,), interval, [x])
```

This call is exercised by the regression contract “coordinate range diameter nearest and support operations”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## scale

`scale(region, factor, variables=None) -> 'sp.Expr'`

Scale a region about the origin by a scalar factor.

**Parameter semantics.**

- `region` — Semialgebraic or standard region on which the operation is performed.
- `factor` — Scale factor applied to the region or expression.
- `variables` — Ordered real variables of the problem; their order is semantically significant for elimination, coordinates, or returned tuples.

**Result semantics.** The return contract is `sp.Expr`; see the family reference for structured-result details and certification status.

**Representative regression-backed call.**

```python
scale(region, lam, [x])
```

This call is exercised by the regression contract “positive scaling scales diameter and width”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## squared_distance_range

`squared_distance_range(left, right=None, variables=None, *, return_result: 'bool' = False)`

Return the exact range of squared pairwise distances.

**Parameter semantics.**

- `left` — Left operand/region in the binary operation.
- `right` — Right operand/region in the binary operation.
- `variables` — Ordered real variables of the problem; their order is semantically significant for elimination, coordinates, or returned tuples.
- `return_result` — When true, return the structured result/certificate/diagnostic object instead of only the convenience value.

**Result semantics.** When ``right`` is omitted, both points range over ``left``.

**Representative regression-backed call.**

```python
squared_distance_range(left_point, right_point, (x,))
```

This call is exercised by the regression contract “interval geometry queries cross check one another”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## sublevel_set

`sublevel_set(expression, value=0, region=True, *, strict: 'bool' = False) -> 'sp.Expr'`

Return ``region ∩ {expression <= value}`` (or strict variant).

**Parameter semantics.**

- `expression` — Symbolic expression being analyzed.
- `value` — Level/threshold value defining a level, sublevel, or superlevel set.
- `region` — Semialgebraic or standard region on which the operation is performed.
- `strict` — Select strict failure behavior for incomplete/unsupported computation where the API provides it.

**Result semantics.** The return contract is `sp.Expr`; see the family reference for structured-result details and certification status.

**Representative regression-backed call.**

```python
sublevel_set(x**2, 1, region)
```

This call is exercised by the regression contract “level set constructors preserve region constraint”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## superlevel_set

`superlevel_set(expression, value=0, region=True, *, strict: 'bool' = False) -> 'sp.Expr'`

Return ``region ∩ {expression >= value}`` (or strict variant).

**Parameter semantics.**

- `expression` — Symbolic expression being analyzed.
- `value` — Level/threshold value defining a level, sublevel, or superlevel set.
- `region` — Semialgebraic or standard region on which the operation is performed.
- `strict` — Select strict failure behavior for incomplete/unsupported computation where the API provides it.

**Result semantics.** The return contract is `sp.Expr`; see the family reference for structured-result details and certification status.

**Representative regression-backed call.**

```python
superlevel_set(x**2, 1, region)
```

This call is exercised by the regression contract “level set constructors preserve region constraint”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## support_function

`support_function(region, direction, variables=None, *, return_result: 'bool' = False)`

Return ``sup(x·direction)`` over the region.

**Parameter semantics.**

- `region` — Semialgebraic or standard region on which the operation is performed.
- `direction` — Direction vector used by support/width or directional geometry.
- `variables` — Ordered real variables of the problem; their order is semantically significant for elimination, coordinates, or returned tuples.
- `return_result` — When true, return the structured result/certificate/diagnostic object instead of only the convenience value.

**Result semantics.** The return contract is `the documented result type for this operation`; see the family reference for structured-result details and certification status.

**Representative regression-backed call.**

```python
support_function(right, (1,), [x])
```

This call is exercised by the regression contract “minkowski support function identity for intervals”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## translate

`translate(region, vector, variables=None) -> 'sp.Expr'`

Translate a region by ``vector`` while preserving coordinate symbols.

**Parameter semantics.**

- `region` — Semialgebraic or standard region on which the operation is performed.
- `vector` — Translation/linear vector in the ambient coordinate system.
- `variables` — Ordered real variables of the problem; their order is semantically significant for elimination, coordinates, or returned tuples.

**Result semantics.** The return contract is `sp.Expr`; see the family reference for structured-result details and certification status.

**Representative regression-backed call.**

```python
translate(interval, (2,), [x])
```

This call is exercised by the regression contract “affine image translate scale and minkowski sum”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## width

`width(region, direction, variables=None) -> 'sp.Expr'`

Return exact directional width ``max u·x - min u·x``.

**Parameter semantics.**

- `region` — Semialgebraic or standard region on which the operation is performed.
- `direction` — Direction vector used by support/width or directional geometry.
- `variables` — Ordered real variables of the problem; their order is semantically significant for elimination, coordinates, or returned tuples.

**Result semantics.** The return contract is `sp.Expr`; see the family reference for structured-result details and certification status.

**Representative regression-backed call.**

```python
width(region, (1,), [x])
```

This call is exercised by the regression contract “width equals two sided support sum”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## bounding_box

`bounding_box(region, variables=None, *, return_result: 'bool' = False)`

Compute the exact axis-aligned bounding box by coordinate optimization.

**Parameter semantics.**

- `region` — Semialgebraic or standard region on which the operation is performed.
- `variables` — Ordered real variables of the problem; their order is semantically significant for elimination, coordinates, or returned tuples.
- `return_result` — When true, return the structured result/certificate/diagnostic object instead of only the convenience value.

**Result semantics.** The return contract is `the documented result type for this operation`; see the family reference for structured-result details and certification status.

**Representative regression-backed call.**

```python
bounding_box(interval, (x,))
```

This call is exercised by the regression contract “interval geometry queries cross check one another”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## critical_values

`critical_values(expression, region=True, variables=None) -> 'tuple[sp.Expr, ...]'`

Return exact objective values from isolated and constant KKT components.

**Parameter semantics.**

- `expression` — Symbolic expression being analyzed.
- `region` — Semialgebraic or standard region on which the operation is performed.
- `variables` — Ordered real variables of the problem; their order is semantically significant for elimination, coordinates, or returned tuples.

**Result semantics.** The return contract is `tuple[sp.Expr, ...]`; see the family reference for structured-result details and certification status.

**Representative regression-backed call.**

```python
critical_values(x**2, interval, (x,))
```

This call is exercised by the regression contract “interval extrema and level sets have exact endpoint semantics”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## critical_value_image

`critical_value_image(expression, region=True, variables=None, *, value_symbol: 'sp.Symbol | str | None' = None) -> 'CriticalValueImage'`

Return isolated values and exact images of positive-dimensional critical loci.

**Parameter semantics.**

- `expression` — Symbolic expression being analyzed.
- `region` — Semialgebraic or standard region on which the operation is performed.
- `variables` — Ordered real variables of the problem; their order is semantically significant for elimination, coordinates, or returned tuples.
- `value_symbol` — Symbol used to represent an exact objective/range value.

**Result semantics.** Unlike :func:`critical_values`, this operation does not discard a critical component merely because the objective varies on it.  Each supported KKT or

**Representative regression-backed call.**

```python
critical_value_image(x**2, sp.true, (x,))
```

This call is exercised by the regression contract “critical value image public pipeline handles discrete case”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## distance_between_regions

`distance_between_regions(left, right, variables=None, *, return_result: 'bool' = False)`

Compute exact Euclidean distance between two semialgebraic regions.

**Parameter semantics.**

- `left` — Left operand/region in the binary operation.
- `right` — Right operand/region in the binary operation.
- `variables` — Ordered real variables of the problem; their order is semantically significant for elimination, coordinates, or returned tuples.
- `return_result` — When true, return the structured result/certificate/diagnostic object instead of only the convenience value.

**Result semantics.** The return contract is `the documented result type for this operation`; see the family reference for structured-result details and certification status.

**Representative regression-backed call.**

```python
distance_between_regions(x <= 0, x >= 2, [x])
```

This call is exercised by the regression contract “distances are exact”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## distance_to_region

`distance_to_region(point, region, variables=None, *, return_result: 'bool' = False)`

Compute exact Euclidean distance from a point to a semialgebraic region.

**Parameter semantics.**

- `point` — Point in the ambient coordinate system.
- `region` — Semialgebraic or standard region on which the operation is performed.
- `variables` — Ordered real variables of the problem; their order is semantically significant for elimination, coordinates, or returned tuples.
- `return_result` — When true, return the structured result/certificate/diagnostic object instead of only the convenience value.

**Result semantics.** The return contract is `the documented result type for this operation`; see the family reference for structured-result details and certification status.

**Representative regression-backed call.**

```python
sa.distance_to_region((5,), region, (X,))
```

This call is exercised by the regression contract “nearest point realizes distance to region”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## fiber

`fiber(region, substitutions: 'Mapping[sp.Symbol | str, sp.Expr]') -> 'sp.Expr'`

Specialize a semialgebraic family at fixed parameter/coordinate values.

**Parameter semantics.**

- `region` — Semialgebraic or standard region on which the operation is performed.
- `substitutions` — Symbolic substitutions defining a parameterization or transformed representation.

**Result semantics.** String keys are resolved against the actual symbols in ``region`` and an ambiguous same-name symbol is rejected rather than guessed.

**Representative regression-backed call.**

```python
fiber(formula, {xr: 1})
```

This call is exercised by the regression contract “fiber rejects ambiguous same name symbol”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## is_path_connected

`is_path_connected(region, variables=None, *, max_pair_checks: 'int | None' = None) -> 'bool'`

Decide path connectedness via exact CAD connectivity.

**Parameter semantics.**

- `region` — Semialgebraic or standard region on which the operation is performed.
- `variables` — Ordered real variables of the problem; their order is semantically significant for elimination, coordinates, or returned tuples.
- `max_pair_checks` — Guard on pairwise comparisons/intersections.

**Result semantics.** Semialgebraic connected sets are semialgebraically path connected, so the existing certified CAD connected-component graph suffices for the decision.

**Representative regression-backed call.**

```python
is_path_connected(region, [x])
```

This call is exercised by the regression contract “path connectivity and path chain 1d”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## path_between

`path_between(region, start, end, variables=None, *, max_pair_checks: 'int | None' = None) -> 'CADPathResult'`

Return a certified CAD cell-chain connecting two points in a region.

**Parameter semantics.**

- `region` — Semialgebraic or standard region on which the operation is performed.
- `start` — Start point/value.
- `end` — End point/value.
- `variables` — Ordered real variables of the problem; their order is semantically significant for elimination, coordinates, or returned tuples.
- `max_pair_checks` — Guard on pairwise comparisons/intersections.

**Result semantics.** For a one-dimensional connected component the waypoint chain is an explicit piecewise-linear path. In higher dimensions the current implementation

**Representative regression-backed call.**

```python
path_between(region, (0,), (1,), [x])
```

This call is exercised by the regression contract “path connectivity and path chain 1d”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## semialgebraic_projection

`semialgebraic_projection(region, eliminate: 'Sequence[sp.Symbol | str]', variables: 'Sequence[sp.Symbol | str] | None' = None) -> 'sp.Expr'`

Project ``region`` by existentially eliminating the requested variables.

**Parameter semantics.**

- `region` — Semialgebraic or standard region on which the operation is performed.
- `eliminate` — Variables to project/eliminate from the input set.
- `variables` — Ordered real variables of the problem; their order is semantically significant for elimination, coordinates, or returned tuples.

**Result semantics.** The return contract is `sp.Expr`; see the family reference for structured-result details and certification status.

**Representative regression-backed call.**

```python
semialgebraic_projection(X >= 0, (Y,), (X,))
```

This call is exercised by the regression contract “geometry rejects mismatched coordinates”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## connected_component_count

`connected_component_count(region, variables=None) -> 'int'`

Return the exact number of semialgebraically connected components.

**Parameter semantics.**

- `region` — Semialgebraic or standard region on which the operation is performed.
- `variables` — Ordered real variables of the problem; their order is semantically significant for elimination, coordinates, or returned tuples.

**Result semantics.** The return contract is `int`; see the family reference for structured-result details and certification status.

**Representative regression-backed call.**

```python
connected_component_count(region, (x,))
```

This call is exercised by the regression contract “zero dimensional components match real roots”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## connected_component_samples

`connected_component_samples(region, variables=None)`

Return one exact CAD sample point from every connected component.

**Parameter semantics.**

- `region` — Semialgebraic or standard region on which the operation is performed.
- `variables` — Ordered real variables of the problem; their order is semantically significant for elimination, coordinates, or returned tuples.

**Result semantics.** The return contract is `the documented result type for this operation`; see the family reference for structured-result details and certification status.

**Representative regression-backed call.**

```python
connected_component_samples(region, (x,))
```

This call is exercised by the regression contract “zero dimensional components match real roots”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## topology_summary

`topology_summary(region, variables=None, *, compact_support: 'bool' = True) -> 'TopologySummary'`

Return exact component, Euler, and supported Betti-number information.

**Parameter semantics.**

- `region` — Semialgebraic or standard region on which the operation is performed.
- `variables` — Ordered real variables of the problem; their order is semantically significant for elimination, coordinates, or returned tuples.
- `compact_support` — Optional compact support restriction needed by operations that otherwise may be unbounded.

**Result semantics.** ``b_0`` is available in every dimension. For compact sets of dimension at most one, ``b_1`` follows exactly from ``chi = b_0 - b_1``. Higher Betti

**Representative regression-backed call.**

```python
topology_summary(disk, (x, y))
```

This call is exercised by the regression contract “compact convex disk has trivial positive betti numbers”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## betti_number

`betti_number(region, degree: 'int', variables=None) -> 'int'`

Return a certified Betti number when the current exact backend supports it.

**Parameter semantics.**

- `region` — Semialgebraic or standard region on which the operation is performed.
- `degree` — Polynomial/generation degree.
- `variables` — Ordered real variables of the problem; their order is semantically significant for elimination, coordinates, or returned tuples.

**Result semantics.** The return contract is `int`; see the family reference for structured-result details and certification status.

**Representative regression-backed call.**

```python
betti_number(box, 2)
```

This call is exercised by the regression contract “polyhedral triangulation drives exact higher betti numbers”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## clip_affine_subspace_to_box

`clip_affine_subspace_to_box(point: 'Sequence[object]', directions: 'Sequence[Sequence[object]]', bounds: 'Sequence[Sequence[object]]') -> 'AffineBoxClip'`

Clip a low-dimensional affine subspace to a box without CAD.

**Parameter semantics.**

- `point` — Point in the ambient coordinate system.
- `directions` — Directions controlling this operation; see the signature type/default for the accepted representation.
- `bounds` — Exact coordinate bounds restricting search or construction.

**Result semantics.** The specialization is limited to parameter dimension one or two, where enumerating active box facets is cheaper and simpler than a CAD.

**Representative regression-backed call.**

```python
semialg.clip_affine_subspace_to_box(
    (0, 0),
    ((1, 1),),
    ((0, 1), (0, 1)),
)
```

This call is exercised by the regression contract “affine box clip is obtained through clipping operation”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## local_dimension

`local_dimension(region: 'object', point, variables: 'Sequence[sp.Symbol] | None' = None) -> 'int'`

Exact local semialgebraic dimension at a point.

**Parameter semantics.**

- `region` — Semialgebraic or standard region on which the operation is performed.
- `point` — Point in the ambient coordinate system.
- `variables` — Ordered real variables of the problem; their order is semantically significant for elimination, coordinates, or returned tuples.

**Result semantics.** The local dimension is the maximum dimension of a selected CAD cell whose closure contains the point. Returns ``-1`` when the point is not in the

**Representative regression-backed call.**

```python
cusp.local_dimension((0, 0))
```

This call is exercised by the regression contract “cusp detects origin and local dimension”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## region_boundary_result

`region_boundary_result(region: 'object', variables: 'Sequence[sp.Symbol] | None' = None) -> 'RegionBoundaryResult'`

Return exact boundary cells, membership status, active residuals, and CAD.

**Parameter semantics.**

- `region` — Semialgebraic or standard region on which the operation is performed.
- `variables` — Ordered real variables of the problem; their order is semantically significant for elimination, coordinates, or returned tuples.

**Result semantics.** The decomposition is computed once and retained on the result so plotting, meshing, singularity analysis, and topology code can reuse the same exact

**Representative regression-backed call.**

```python
region_boundary_result(closed)
```

This call is exercised by the regression contract “region regularization and variable contracts”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## region_active_boundary_strata

`region_active_boundary_strata(region: 'object', variables: 'Sequence[sp.Symbol] | None' = None) -> 'tuple[ActiveBoundaryStratum, ...]'`

Stratify the exact boundary by realized active inequality constraints.

**Parameter semantics.**

- `region` — Semialgebraic or standard region on which the operation is performed.
- `variables` — Ordered real variables of the problem; their order is semantically significant for elimination, coordinates, or returned tuples.

**Result semantics.** This consumes :func:`region_boundary_result`, reusing its exact boundary CAD and active-residual metadata instead of launching an independent

**Representative regression-backed call.**

```python
region_active_boundary_strata(region, (x, y))
```

This call is exercised by the regression contract “disjunctive region reports only realized active boundary cells”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## region_nonsmooth_locus

`region_nonsmooth_locus(region: 'object', variables: 'Sequence[sp.Symbol] | None' = None) -> 'sp.Expr'`

Return the exact recognized nonsmooth/corner locus of a region boundary.

**Parameter semantics.**

- `region` — Semialgebraic or standard region on which the operation is performed.
- `variables` — Ordered real variables of the problem; their order is semantically significant for elimination, coordinates, or returned tuples.

**Result semantics.** Algebraic singularities come from :func:`region_singular_locus`. Corners and ridges are added only on realized active-boundary strata containing at

**Representative regression-backed call.**

```python
region_nonsmooth_locus(region, (x, y))
```

This call is exercised by the regression contract “nonsmooth locus uses realized active sets for disjunction corner”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## region_regular_locus

`region_regular_locus(region: 'object', variables: 'Sequence[sp.Symbol] | None' = None) -> 'sp.Expr'`

Return the part of ``region`` outside its algebraic boundary singular locus.

**Parameter semantics.**

- `region` — Semialgebraic or standard region on which the operation is performed.
- `variables` — Ordered real variables of the problem; their order is semantically significant for elimination, coordinates, or returned tuples.

**Result semantics.** The return contract is `sp.Expr`; see the family reference for structured-result details and certification status.

**Representative regression-backed call.**

```python
region_regular_locus(closed)
```

This call is exercised by the regression contract “region regularization and variable contracts”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## region_singular_locus

`region_singular_locus(region: 'object', variables: 'Sequence[sp.Symbol] | None' = None) -> 'sp.Expr'`

Return the exact reduced-real singular locus on the actual boundary.

**Parameter semantics.**

- `region` — Semialgebraic or standard region on which the operation is performed.
- `variables` — Ordered real variables of the problem; their order is semantically significant for elimination, coordinates, or returned tuples.

**Result semantics.** This formula-only convenience API is strict: if component- sensitive regularity cannot be certified, it raises ``NotImplementedError``

**Representative regression-backed call.**

```python
semialg.region_singular_locus(region)
```

This call is exercised by the regression contract “region analysis primary function contracts”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## region_singular_locus_result

`region_singular_locus_result(region: 'object', variables: 'Sequence[sp.Symbol] | None' = None) -> 'SingularLocusResult'`

Return singular-locus geometry together with completeness certification.

**Parameter semantics.**

- `region` — Semialgebraic or standard region on which the operation is performed.
- `variables` — Ordered real variables of the problem; their order is semantically significant for elimination, coordinates, or returned tuples.

**Result semantics.** Exact formulas are returned only when every algebraic decomposition needed by the component-relative Jacobian analysis is certified complete.  When a

**Representative regression-backed call.**

```python
semialg.region_singular_locus_result(region)
```

This call is exercised by the regression contract “region analysis primary function contracts”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## region_boundary

`region_boundary(region: 'FormulaLike | Iterable[FormulaLike]', variables: 'Sequence[sp.Symbol | str] | None' = None, *, strategy: 'str | None' = None) -> 'sp.Expr'`

Return the Euclidean boundary of a semialgebraic region.

**Parameter semantics.**

- `region` — Semialgebraic or standard region on which the operation is performed.
- `variables` — Ordered real variables of the problem; their order is semantically significant for elimination, coordinates, or returned tuples.
- `strategy` — Algorithm-selection policy. Use the default automatic policy unless a specific exact backend is being tested or diagnosed.

**Result semantics.** The default implementation is the CAD-semantic intersection ``closure(S) ∩ closure(complement(S))``.  This removes false boundaries at

**Representative regression-backed call.**

```python
region_boundary(region, [x])
```

This call is exercised by the regression contract “touching closed intervals have no internal boundary seam”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## region_complement

`region_complement(region: 'FormulaLike | Iterable[FormulaLike]')`

Return the complement of an implicit or unified semialgebraic region.

**Parameter semantics.**

- `region` — Semialgebraic or standard region on which the operation is performed.

**Result semantics.** The return contract is `the documented result type for this operation`; see the family reference for structured-result details and certification status.

**Representative regression-backed call.**

```python
region_complement(x > 0)
```

This call is exercised by the regression contract “boolean region operations”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## region_difference

`region_difference(lhs: 'FormulaLike | Iterable[FormulaLike]', rhs: 'FormulaLike | Iterable[FormulaLike]') -> 'sp.Expr'`

Return ``lhs`` minus ``rhs`` for implicit or unified regions.

**Parameter semantics.**

- `lhs` — Left formula in the comparison.
- `rhs` — Right formula in the comparison.

**Result semantics.** The return contract is `sp.Expr`; see the family reference for structured-result details and certification status.

**Representative regression-backed call.**

```python
region_difference(left, right)
```

This call is exercised by the regression contract “region operations accept unified and explicit regions”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## region_interior

`region_interior(region: 'FormulaLike | Iterable[FormulaLike]', variables: 'Sequence[sp.Symbol | str] | None' = None, *, strategy: 'str | None' = None) -> 'sp.Expr'`

Return the Euclidean interior of a semialgebraic region.

**Parameter semantics.**

- `region` — Semialgebraic or standard region on which the operation is performed.
- `variables` — Ordered real variables of the problem; their order is semantically significant for elimination, coordinates, or returned tuples.
- `strategy` — Algorithm-selection policy. Use the default automatic policy unless a specific exact backend is being tested or diagnosed.

**Result semantics.** By default this is computed as the CAD-semantic complement of the closure of the complement, so internal CAD sections are retained when their entire

**Representative regression-backed call.**

```python
region_interior(region, [x])
```

This call is exercised by the regression contract “touching closed intervals have no internal boundary seam”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## region_symmetric_difference

`region_symmetric_difference(lhs: 'FormulaLike | Iterable[FormulaLike]', rhs: 'FormulaLike | Iterable[FormulaLike]')`

Return the exact symmetric difference of two semialgebraic regions.

**Parameter semantics.**

- `lhs` — Left formula in the comparison.
- `rhs` — Right formula in the comparison.

**Result semantics.** The result contains points belonging to exactly one operand.  Unified and explicit region inputs are lowered through the same exact representation

**Representative regression-backed call.**

```python
sa.region_symmetric_difference(a, b)
```

This call is exercised by the regression contract “symmetric difference is commutative”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## affine_image

`affine_image(region, matrix, offset=None, variables=None)`

Return the exact affine image ``A*x+b`` of a region.

**Parameter semantics.**

- `region` — Semialgebraic or standard region on which the operation is performed.
- `matrix` — Matrix defining the linear/affine operation or quadratic data.
- `offset` — Offset applied to a generated regular object.
- `variables` — Ordered real variables of the problem; their order is semantically significant for elimination, coordinates, or returned tuples.

**Result semantics.** The function dispatches on representation.  Canonical geometry and :class:`HRepresentation` inputs preserve structure whenever possible.

**Representative regression-backed call.**

```python
affine_image(tri, A, (4, 5))
```

This call is exercised by the regression contract “affine image preserves point and simplex structure”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## affine_preimage

`affine_preimage(region, matrix, offset=None, variables=None, *, target_variables=None)`

Return the exact preimage under ``x -> A*x+b``.

**Parameter semantics.**

- `region` — Semialgebraic or standard region on which the operation is performed.
- `matrix` — Matrix defining the linear/affine operation or quadratic data.
- `offset` — Offset applied to a generated regular object.
- `variables` — Ordered real variables of the problem; their order is semantically significant for elimination, coordinates, or returned tuples.
- `target_variables` — Coordinate symbols used for the target representation.

**Result semantics.** Canonical inputs preserve structure for invertible square maps and otherwise lower to an exact symbolic region.  :class:`SemialgebraicRegion` inputs

**Representative regression-backed call.**

```python
affine_preimage(hrep, ((1, 0),))
```

This call is exercised by the regression contract “hrepresentation rectangular preimage is supported”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## region_image

`region_image(region, mapping, variables=None, *, image_variables=None, parameters=None)`

Return the exact image under a symbolic map.

**Parameter semantics.**

- `region` — Semialgebraic or standard region on which the operation is performed.
- `mapping` — Polynomial or affine map whose image/local geometry is analyzed.
- `variables` — Ordered real variables of the problem; their order is semantically significant for elimination, coordinates, or returned tuples.
- `image_variables` — Coordinate symbols used for the image space.
- `parameters` — Symbols treated as free parameters rather than eliminated decision variables.

**Result semantics.** Canonical inputs preserve canonical affine structure and otherwise return a lazy :class:`TransformedRegion`.  :class:`SemialgebraicRegion` inputs use

**Representative regression-backed call.**

```python
region_image(source, (x**2,), (x,))
```

This call is exercised by the regression contract “polynomial image lowers to semialgebraic formula”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## region_preimage

`region_preimage(target, mapping, variables=None, *, target_variables=None)`

Return an exact symbolic preimage, preserving invertible affine structure when possible.

**Parameter semantics.**

- `target` — Requested target object or result; for certificate APIs this is the relation to prove, while decomposition APIs use it to select the desired cell representation.
- `mapping` — Polynomial or affine map whose image/local geometry is analyzed.
- `variables` — Ordered real variables of the problem; their order is semantically significant for elimination, coordinates, or returned tuples.
- `target_variables` — Coordinate symbols used for the target representation.

**Result semantics.** The return contract is `the documented result type for this operation`; see the family reference for structured-result details and certification status.

**Representative regression-backed call.**

```python
region_preimage(target, (x**2,), (x,))
```

This call is exercised by the regression contract “polynomial preimage is exact”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## deduplicate_indexed_vertices

`deduplicate_indexed_vertices(vertices: 'Sequence[Sequence[object]]', cells: 'Sequence[Sequence[int]]') -> 'IndexedVertexData'`

Deduplicate exact vertices and remap zero-based cell indices.

**Parameter semantics.**

- `vertices` — Vertices defining the polygon/polytope.
- `cells` — Cells participating in the decomposition.

**Result semantics.** Equality is SymPy structural equality after sympification; no numerical tolerance is introduced.  The first occurrence of each coordinate is kept.

**Representative regression-backed call.**

```python
sa.deduplicate_indexed_vertices(first.vertices, first.cells)
```

This call is exercised by the regression contract “deduplicate indexed vertices is idempotent and reindexes cells”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## polygon_vertices

`polygon_vertices(region: 'Polygon | PolygonalSet') -> 'tuple[_PointData, ...]'`

Return unique polygon vertices in deterministic first-occurrence order.

**Parameter semantics.**

- `region` — Semialgebraic or standard region on which the operation is performed.

**Result semantics.** The return contract is `tuple[_PointData, ...]`; see the family reference for structured-result details and certification status.

**Representative regression-backed call.**

```python
polygon_vertices(region)
```

This call is exercised by the regression contract “polygonal set orients outer and hole boundaries”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## outer_polygons

`outer_polygons(region: 'Polygon | PolygonalSet') -> 'tuple[Polygon, ...]'`

Return all outer polygon boundaries.

**Parameter semantics.**

- `region` — Semialgebraic or standard region on which the operation is performed.

**Result semantics.** The return contract is `tuple[Polygon, ...]`; see the family reference for structured-result details and certification status.

**Representative regression-backed call.**

```python
outer_polygons(region)
```

This call is exercised by the regression contract “polygonal set orients outer and hole boundaries”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## inner_polygons

`inner_polygons(region: 'Polygon | PolygonalSet') -> 'tuple[Polygon, ...]'`

Return all polygon hole boundaries.

**Parameter semantics.**

- `region` — Semialgebraic or standard region on which the operation is performed.

**Result semantics.** The return contract is `tuple[Polygon, ...]`; see the family reference for structured-result details and certification status.

**Representative regression-backed call.**

```python
inner_polygons(region)
```

This call is exercised by the regression contract “polygonal set orients outer and hole boundaries”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## polyhedron_vertices

`polyhedron_vertices(region: 'PolyhedralShell | Polyhedron') -> 'tuple[_PointData, ...]'`

Return unique boundary vertices in deterministic first-occurrence order.

**Parameter semantics.**

- `region` — Semialgebraic or standard region on which the operation is performed.

**Result semantics.** The return contract is `tuple[_PointData, ...]`; see the family reference for structured-result details and certification status.

**Representative regression-backed call.**

```python
polyhedron_vertices(solid)
```

This call is exercised by the regression contract “polyhedron keeps components and cavities explicit”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## polyhedron_face_indices

`polyhedron_face_indices(region: 'PolyhedralShell | Polyhedron') -> 'tuple[tuple[tuple[int, ...], ...], ...]'`

Return face indices grouped by boundary shell.

**Parameter semantics.**

- `region` — Semialgebraic or standard region on which the operation is performed.

**Result semantics.** Each shell keeps its own local vertex indexing.  This avoids silently changing topology merely to manufacture one global coordinate array.

**Representative regression-backed call.**

```python
polyhedron_face_indices(shell)
```

This call is exercised by the regression contract “polyhedral shell requires closed consistently oriented manifold”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## outer_polyhedra

`outer_polyhedra(region: 'PolyhedralShell | Polyhedron') -> 'tuple[PolyhedralShell, ...]'`

Return all outer boundary shells.

**Parameter semantics.**

- `region` — Semialgebraic or standard region on which the operation is performed.

**Result semantics.** The return contract is `tuple[PolyhedralShell, ...]`; see the family reference for structured-result details and certification status.

**Representative regression-backed call.**

```python
outer_polyhedra(solid)
```

This call is exercised by the regression contract “polyhedron keeps components and cavities explicit”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## inner_polyhedra

`inner_polyhedra(region: 'PolyhedralShell | Polyhedron') -> 'tuple[PolyhedralShell, ...]'`

Return all cavity boundary shells.

**Parameter semantics.**

- `region` — Semialgebraic or standard region on which the operation is performed.

**Result semantics.** The return contract is `tuple[PolyhedralShell, ...]`; see the family reference for structured-result details and certification status.

**Representative regression-backed call.**

```python
inner_polyhedra(solid)
```

This call is exercised by the regression contract “polyhedron keeps components and cavities explicit”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## convex_hull

`convex_hull(points, *, canonical=True)`

Return the exact convex hull of a finite point set.

**Parameter semantics.**

- `points` — Input points, interpreted in the supplied coordinate order.
- `canonical` — Whether to canonicalize representation/order.

**Result semantics.** Redundant and interior input points are removed exactly.  Lower-dimensional point sets are solved in an exact intrinsic affine chart and mapped back.

**Representative regression-backed call.**

```python
semialg.convex_hull([(0, 0), (2, 0), (0, 2), (1, 1)])
```

This call is exercised by the regression contract “exact hull and seeded generation public calls”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## random_polygon

`random_polygon(*, vertex_count=6, coordinate_bound=10, seed=None)`

Generate a reproducible exact convex lattice polygon.

**Parameter semantics.**

- `vertex_count` — Number of generated vertices.
- `coordinate_bound` — Symmetric coordinate bound used by bounded search.
- `seed` — Deterministic random seed used only by randomized search/generation stages.

**Result semantics.** The return contract is `the documented result type for this operation`; see the family reference for structured-result details and certification status.

**Representative regression-backed call.**

```python
semialg.random_polygon(seed=4)
```

This call is exercised by the regression contract “random generators are seeded exact”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## random_polytope

`random_polytope(dimension=3, *, point_count=None, coordinate_bound=10, seed=None)`

Generate a reproducible exact full-dimensional random lattice polytope.

**Parameter semantics.**

- `dimension` — Requested intrinsic/ambient dimension where the constructor or analysis needs it.
- `point_count` — Number of generated/sample points.
- `coordinate_bound` — Symmetric coordinate bound used by bounded search.
- `seed` — Deterministic random seed used only by randomized search/generation stages.

**Result semantics.** The return contract is `the documented result type for this operation`; see the family reference for structured-result details and certification status.

**Representative regression-backed call.**

```python
semialg.random_polytope(3, seed=7)
```

This call is exercised by the regression contract “random generators are seeded exact”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## subdivide_triangular_faces

`subdivide_triangular_faces(source, *, levels=1)`

Subdivide every triangular face into four triangles using shared exact midpoints.

**Parameter semantics.**

- `source` — Source region/object being mapped.
- `levels` — Number or sequence of refinement levels.

**Result semantics.** The return contract is `the documented result type for this operation`; see the family reference for structured-result details and certification status.

**Representative regression-backed call.**

```python
semialg.subdivide_triangular_faces(ico, levels=1)
```

This call is exercised by the regression contract “triangular subdivision and geodesic projection”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## geodesic_refinement

`geodesic_refinement(source, *, levels=1, center=(0, 0, 0), radius=None)`

Refine a triangular shell and project its vertices exactly to a sphere.

**Parameter semantics.**

- `source` — Source region/object being mapped.
- `levels` — Number or sequence of refinement levels.
- `center` — Center of the generated/transformed object.
- `radius` — Radius parameter of the generated object.

**Result semantics.** The return contract is `the documented result type for this operation`; see the family reference for structured-result details and certification status.

**Representative regression-backed call.**

```python
semialg.geodesic_refinement(ico, levels=1)
```

This call is exercised by the regression contract “triangular subdivision and geodesic projection”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## polyhedral_intersection

`polyhedral_intersection(left, right)`

Return the exact canonical intersection of two full-dimensional convex polytopes.

**Parameter semantics.**

- `left` — Left operand/region in the binary operation.
- `right` — Right operand/region in the binary operation.

**Result semantics.** Returns ``None`` for an empty intersection and raises ``NotImplementedError`` when the exact structural backend does not cover the input dimensions.

**Representative regression-backed call.**

```python
sa.polyhedral_intersection(a, b)
```

This call is exercised by the regression contract “polyhedral intersection and boolean intersection agree”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## polyhedral_boolean

`polyhedral_boolean(operation, left, right)`

Exact structure-preserving convex-polyhedral Boolean fast path.

**Parameter semantics.**

- `operation` — Topological/set operation requested from CAD.
- `left` — Left operand/region in the binary operation.
- `right` — Right operand/region in the binary operation.

**Result semantics.** Intersection is constructed geometrically. Union/difference preserve a canonical operand when containment makes that exact without subdivision;

**Representative regression-backed call.**

```python
sa.polyhedral_boolean("intersection", a, b)
```

This call is exercised by the regression contract “polyhedral intersection and boolean intersection agree”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## RegularPolygon

`RegularPolygon(sides: 'int', center: 'Sequence[object]' = (0, 0), radius: 'object' = 1, *, rotation: 'object' = 0) -> 'Polygon'`

Return a canonical regular polygon as a :class:`Polygon`.

**Parameter semantics.**

- `sides` — Number of polygon/prism/pyramid sides.
- `center` — Center of the generated/transformed object.
- `radius` — Radius parameter of the generated object.
- `rotation` — Rotation applied to a generated regular object.

**Result semantics.** Vertices lie on the circumcircle with the first vertex at ``rotation`` radians from the positive x-axis.

**Representative regression-backed call.**

```python
RegularPolygon(4)
```

This call is exercised by the regression contract “regular polygon returns canonical polygon”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## Cube

`Cube(center: 'Sequence[object]' = (0, 0, 0), side: 'object' = 1) -> 'Parallelepiped'`

Return an axis-aligned cube as a canonical :class:`Parallelepiped`.

**Parameter semantics.**

- `center` — Center of the generated/transformed object.
- `side` — Side-length parameter.

**Result semantics.** The return contract is `Parallelepiped`; see the family reference for structured-result details and certification status.

**Representative regression-backed call.**

```python
Cube(side=0)
```

This call is exercised by the regression contract “named solid positive lengths are validated”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## Tetrahedron

`Tetrahedron(vertices: 'Sequence[Sequence[object]] | None' = None, *, center: 'Sequence[object]' = (0, 0, 0), edge: 'object' = 1) -> 'Simplex'`

Return a tetrahedron as a canonical :class:`Simplex`.

**Parameter semantics.**

- `vertices` — Vertices defining the polygon/polytope.
- `center` — Center of the generated/transformed object.
- `edge` — Edge identifier used by mesh/refinement operations.

**Result semantics.** With explicit ``vertices`` those vertices are used directly. Otherwise a regular tetrahedron centered at ``center`` with the requested edge length

**Representative regression-backed call.**

```python
Tetrahedron(vertices)
```

This call is exercised by the regression contract “explicit tetrahedron preserves vertices”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## Octahedron

`Octahedron(center: 'Sequence[object]' = (0, 0, 0), edge: 'object' = 1) -> 'Polytope'`

Return a regular octahedron as a canonical :class:`Polytope`.

**Parameter semantics.**

- `center` — Center of the generated/transformed object.
- `edge` — Edge identifier used by mesh/refinement operations.

**Result semantics.** The return contract is `Polytope`; see the family reference for structured-result details and certification status.

**Representative regression-backed call.**

```python
Octahedron()
```

This call is exercised by the regression contract “triangular subdivision and geodesic projection”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## Icosahedron

`Icosahedron(center: 'Sequence[object]' = (0, 0, 0), edge: 'object' = 1) -> 'Polytope'`

Return a regular icosahedron as a canonical :class:`Polytope`.

**Parameter semantics.**

- `center` — Center of the generated/transformed object.
- `edge` — Edge identifier used by mesh/refinement operations.

**Result semantics.** The return contract is `Polytope`; see the family reference for structured-result details and certification status.

**Representative regression-backed call.**

```python
Icosahedron(edge=3)
```

This call is exercised by the regression contract “regular platonic polytope edges”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## Dodecahedron

`Dodecahedron(center: 'Sequence[object]' = (0, 0, 0), edge: 'object' = 1) -> 'Polytope'`

Return a regular dodecahedron as a canonical :class:`Polytope`.

**Parameter semantics.**

- `center` — Center of the generated/transformed object.
- `edge` — Edge identifier used by mesh/refinement operations.

**Result semantics.** The return contract is `Polytope`; see the family reference for structured-result details and certification status.

**Representative regression-backed call.**

```python
Dodecahedron(edge=3)
```

This call is exercised by the regression contract “regular platonic polytope edges”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## Prism

`Prism(base: 'Sequence[Sequence[object]] | StandardRegion', vector: 'Sequence[object]') -> 'Polytope'`

Extrude a vertex-defined base by ``vector`` and return a polytope.

**Parameter semantics.**

- `base` — Base polygon/region for a prism or pyramid.
- `vector` — Translation/linear vector in the ambient coordinate system.

**Result semantics.** The return contract is `Polytope`; see the family reference for structured-result details and certification status.

**Representative regression-backed call.**

```python
Prism(base, (0, 0, 1))
```

This call is exercised by the regression contract “standard regions validate ambient dimensions”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## Pyramid

`Pyramid(base: 'Sequence[Sequence[object]] | StandardRegion', apex: 'Sequence[object]') -> 'Polytope'`

Join a vertex-defined base to an apex and return a polytope.

**Parameter semantics.**

- `base` — Base polygon/region for a prism or pyramid.
- `apex` — Apex point of a cone/pyramid-like construction.

**Result semantics.** The return contract is `Polytope`; see the family reference for structured-result details and certification status.

**Representative regression-backed call.**

```python
Pyramid(base, (0, 0, 1))
```

This call is exercised by the regression contract “standard regions validate ambient dimensions”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## as_semialgebraic_region

`as_semialgebraic_region(region: 'object', variables: 'Sequence[sp.Symbol | str] | None' = None) -> 'SemialgebraicRegion'`

Coerce a formula or explicit :class:`Geometry` to ``SemialgebraicRegion``.

**Parameter semantics.**

- `region` — Semialgebraic or standard region on which the operation is performed.
- `variables` — Ordered real variables of the problem; their order is semantically significant for elimination, coordinates, or returned tuples.

**Result semantics.** The return contract is `SemialgebraicRegion`; see the family reference for structured-result details and certification status.

**Representative regression-backed call.**

```python
region.as_semialgebraic_region()
```

This call is exercised by the regression contract “membership forms agree”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## interior_of_closure

`interior_of_closure(region: 'object', variables: 'Sequence[sp.Symbol | str] | None' = None, *, strategy: 'str | None' = None) -> 'SemialgebraicRegion'`

Return interior(closure(region)).

**Parameter semantics.**

- `region` — Semialgebraic or standard region on which the operation is performed.
- `variables` — Ordered real variables of the problem; their order is semantically significant for elimination, coordinates, or returned tuples.
- `strategy` — Algorithm-selection policy. Use the default automatic policy unless a specific exact backend is being tested or diagnosed.

**Result semantics.** The return contract is `SemialgebraicRegion`; see the family reference for structured-result details and certification status.

**Representative regression-backed call.**

```python
interior_of_closure(closed)
```

This call is exercised by the regression contract “region regularization and variable contracts”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## closure_of_interior

`closure_of_interior(region: 'object', variables: 'Sequence[sp.Symbol | str] | None' = None, *, strategy: 'str | None' = None) -> 'SemialgebraicRegion'`

Return closure(interior(region)).

**Parameter semantics.**

- `region` — Semialgebraic or standard region on which the operation is performed.
- `variables` — Ordered real variables of the problem; their order is semantically significant for elimination, coordinates, or returned tuples.
- `strategy` — Algorithm-selection policy. Use the default automatic policy unless a specific exact backend is being tested or diagnosed.

**Result semantics.** The return contract is `SemialgebraicRegion`; see the family reference for structured-result details and certification status.

**Representative regression-backed call.**

```python
closure_of_interior(opened)
```

This call is exercised by the regression contract “region regularization and variable contracts”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## region_variables

`region_variables(region: 'object', variables: 'Sequence[sp.Symbol | str] | None' = None, *, kind: 'str' = 'coordinates') -> 'tuple[sp.Symbol, ...]'`

Return coordinate variables, parameters, or all symbols of a region.

**Parameter semantics.**

- `region` — Semialgebraic or standard region on which the operation is performed.
- `variables` — Ordered real variables of the problem; their order is semantically significant for elimination, coordinates, or returned tuples.
- `kind` — Variant of the operation/result requested.

**Result semantics.** The return contract is `tuple[sp.Symbol, ...]`; see the family reference for structured-result details and certification status.

**Representative regression-backed call.**

```python
region_variables(closed)
```

This call is exercised by the regression contract “region regularization and variable contracts”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## is_regular_closed_region

`is_regular_closed_region(region: 'object', variables: 'Sequence[sp.Symbol | str] | None' = None, *, strategy: 'str | None' = None) -> 'bool'`

Return whether a region equals the closure of its interior.

**Parameter semantics.**

- `region` — Semialgebraic or standard region on which the operation is performed.
- `variables` — Ordered real variables of the problem; their order is semantically significant for elimination, coordinates, or returned tuples.
- `strategy` — Algorithm-selection policy. Use the default automatic policy unless a specific exact backend is being tested or diagnosed.

**Result semantics.** The return contract is `bool`; see the family reference for structured-result details and certification status.

**Representative regression-backed call.**

```python
is_regular_closed_region(closed)
```

This call is exercised by the regression contract “regular open and closed region operations”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## is_regular_open_region

`is_regular_open_region(region: 'object', variables: 'Sequence[sp.Symbol | str] | None' = None, *, strategy: 'str | None' = None) -> 'bool'`

Return whether a region equals the interior of its closure.

**Parameter semantics.**

- `region` — Semialgebraic or standard region on which the operation is performed.
- `variables` — Ordered real variables of the problem; their order is semantically significant for elimination, coordinates, or returned tuples.
- `strategy` — Algorithm-selection policy. Use the default automatic policy unless a specific exact backend is being tested or diagnosed.

**Result semantics.** The return contract is `bool`; see the family reference for structured-result details and certification status.

**Representative regression-backed call.**

```python
is_regular_open_region(opened)
```

This call is exercised by the regression contract “regular open and closed region operations”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.

## simplify_region

`simplify_region(region: 'object', variables: 'Sequence[sp.Symbol | str] | None' = None, *, exact: 'bool' = False) -> 'SemialgebraicRegion'`

Canonicalize a symbolic semialgebraic region formula.

**Parameter semantics.**

- `region` — Semialgebraic or standard region on which the operation is performed.
- `variables` — Ordered real variables of the problem; their order is semantically significant for elimination, coordinates, or returned tuples.
- `exact` — Request exact symbolic construction rather than an approximate presentation.

**Result semantics.** The return contract is `SemialgebraicRegion`; see the family reference for structured-result details and certification status.

**Representative regression-backed call.**

```python
simplify_region(expr, [x, y])
```

This call is exercised by the regression contract “canonical simplification is shared and idempotent”, so the example corresponds to behavior checked by the test suite rather than illustrative pseudocode.


## `parametric_cad`

`parametric_cad(formula, variables, *, parameters, ...)` is the canonical parameter-aware CAD operation. Parameters are ordered before ordinary variables in the underlying decomposition, so the result records a cylindrical partition of parameter space together with exact symbolic fibers. The default return is a `ParametricCADResult`, not merely the generic formula.

Use `parameters` for symbols whose regime changes should be exposed; `variables` are the fiber coordinates. `strategy` selects the CAD planner/backend, `assumptions` restrict the problem before decomposition, and `specialize_fibers=False` suppresses representative cylindrical-fiber extraction while retaining exact symbolic strata. `solvability_conditions` remains the lighter operation when only the feasible parameter condition is required.

```python
import sympy as sp
from semialg import parametric_cad

x, a = sp.symbols("x a", real=True)
result = parametric_cad(sp.Eq(a * x, 1), (x,), parameters=(a,))
assert result.parameter_condition != sp.false
assert result.generic_cases
assert result.exceptional_cases
```

`result.strata` preserves both full-dimensional generic and lower-dimensional exceptional parameter cells, including exceptional cells with empty fibers. `result.exceptional_analysis` records degree-drop, coefficient-sign, discriminant, and resultant polynomials discovered as provenance for possible stratum boundaries. The union of nonempty stratum conditions is the same parameter-feasibility question answered by `solvability_conditions`; the richer result additionally retains the decomposition, exceptional regimes, symbolic fibers, representative samples, and optional specialized cylindrical solutions.

For `parametric_cad`, `formula` is the exact semialgebraic relation to decompose and `variables` are the fiber coordinates. `parameters` names the symbolic parameter coordinates. `output="result"` returns the complete structured result; the explicit `formula`, `cases`, `cells`, and `function` selectors project convenience views. `strategy` controls CAD planning, `domain` must currently be real, `assumptions` are conjoined before decomposition, `strict=True` turns an unsupported domain into an immediate error, and `specialize_fibers` controls whether representative cylindrical fiber solutions are extracted.
