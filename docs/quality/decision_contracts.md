# Decision layer quality notes

These notes describe the public contracts for the decision, sampling, and sign-evaluation layer. The checks document expected behavior without adding a
separate solving backend.

## Public contracts under test

1. Boolean-returning calls remain boolean-returning by default
2. `return_result=True` returns structured objects with stable fields
3. Witnesses and counterexamples are validated before exposure
4. `solve_semialgebraic(..., count=0)` performs analysis without collecting samples
5. RUR dispatch is used opportunistically for supported finite systems
6. Sampling strategies have explicit meanings and reject unknown names
7. Exact sign evaluation handles rational, algebraic, and RUR-backed points
8. Boolean formulas are not routed through scalar simplification paths that can trigger SymPy `Mul(And(...))` deprecation warnings.
9. String variable names preserve the exact SymPy symbol objects already present in a problem; ambiguous same-name symbols with different assumptions are rejected.
10. Exact certificate paths do not promote fixed-precision numerical comparisons to exact truth/sign/order claims.

## Focused tests

The main regression files are:

- `tests/test_decision_boolean_contracts.py`
- `tests/test_decision_exact_dispatch.py`
- `tests/test_decision_sampling_contracts.py`
- `tests/test_decision_sign_evaluation.py`
- `tests/test_decision_public_contracts.py`

These tests are intended to be small, public-API-level contracts. They should
remain readable and should not depend on private implementation details except
where a backend method string is part of the documented behavior.

## Quality checklist

Before packaging or merging changes:

```bash
python -m compileall src tests
pytest -q
python scripts/verify_source_quality.py
```

Generated archives should not include `__pycache__`, `.pytest_cache`, scratch
`run_pytest_*` scripts, ad-hoc `check_*.py` files, or test logs.
