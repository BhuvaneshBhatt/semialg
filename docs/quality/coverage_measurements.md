# Coverage measurements

Coverage identifies missing behavior. Tests must check mathematical results, witnesses, endpoint semantics, degree drops or failure contracts; executing a line is insufficient.

The full non-slow baseline measured 3094 passing tests across 425 source modules: **80.37% statements, 68.30% branches and 77.02% combined**. Combined coverage weights statements and branches by their raw counts. The former 85% decision and 80% geometry subsystem figures were unmeasured targets and were not wired into CI; they did not describe achieved coverage.

CI now runs branch coverage and the subsystem checker. The evidence-backed floors leave a small margin below that baseline. Higher targets guide additional behavioral tests.

| Scope | Measured combined baseline | Enforced floor | Next target |
| --- | ---: | ---: | ---: |
| Overall | 77.02% | 76.0% | 78% |
| Decision, parameters and roots | 77.86% | 77.5% | 85% |
| Algebraic certification | 79.95% | 79.5% | 80% |
| Geometry and topology | 75.97% | 75.5% | 80% |

`tests/coverage_policy.toml` partitions all package modules into thirteen subsystems. Missing branch data, omitted modules, overlapping owners, unassigned modules and regressions below a floor fail the checker. Remaining subsystems are reported without an invented floor.

```sh
python -m pytest --cov=semialg --cov-branch --cov-report=json:build/coverage/coverage.json --cov-report=xml:build/coverage/coverage.xml --cov-report=html:build/coverage/html
python scripts/check_coverage.py build/coverage/coverage.json --markdown build/coverage/summary.md
```

Reports are uploaded even when a gate fails. The default coverage run includes every test marker. Scheduled seeded oracles and focused profiling supplement ordinary CI; their reports are kept separately from the complete-suite coverage measurement. A floor is raised after a fresh full measurement supports it.
