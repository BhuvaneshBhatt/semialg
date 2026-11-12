# Development tools

Repository-maintenance utilities belong here when they are needed. Runtime code does not depend on this directory.

## Isolated slow-test runner

`python tools/run_slow_tests.py` collects tests marked `slow` and executes each node in its own Python process, one at a time. These tests also run with plain `pytest`; the isolated runner is a diagnostic tool. It keeps a hang, leaked state, or abnormal process exit in one test from obscuring the rest of the suite. Results are written incrementally to `.test-reports/slow/` as JSON, Markdown, and one log per test.

Use `--soft-timeout` to flag performance outliers and `--hard-timeout` to stop a pathological test and continue with the next node. Collection paths can be supplied positionally, for example:

```bash
python tools/run_slow_tests.py tests/test_local_geometry_property_contracts.py \
  --soft-timeout 20 --hard-timeout 90
```

## Reproducible measurements

`python tools/profile_tests.py --output build/profile -- -m slow tests` runs the unchanged test bodies and records per-phase outcomes, failures and one cProfile file per call. `run.json` contains the exact command, seed, interpreter/dependency versions, module origin and source digest. A process timeout returns 124 and retains partial reports. Inspect a call with `python -m pstats build/profile/<digest>.prof`.

`python tools/benchmark_certification.py --repeats 3 --output benchmarks.json` checks the fixed inputs in `tests/data/performance_cases.json`, clears package and SymPy caches for each cold call, and immediately repeats it warm while retaining the cold result. Reports include raw times, medians and before/after exact cache/work counters. Timings are informational; ordinary tests enforce mathematical results and structural work contracts without machine-specific time thresholds.

`python tools/qe_oracles.py --generated 112 --seed 1729 --output oracle.json` compares four execution paths with exact sign-cell expectations derived from supplied factor roots. The default corpus is also a permanent pytest suite. Scheduled jobs run seeded expansions with minimum/latest dependencies and FLINT/Python ground arithmetic.

`python tools/check_certificate_audit.py certificate-audit.jsonl` reports inactive proof branches and distinguishes mathematical replay from constructor/schema rejection. `python tools/check_distribution.py dist/*.whl dist/*.tar.gz` verifies package bytes across source, installed wheel, sdist and a rebuilt wheel. Install the wheel first; source-tree imports are rejected by the distribution check. The sdist rebuild uses the host build environment (`python -m pip install build "setuptools>=77" wheel`), so CI installs these tools explicitly. The checker also compares packaged documentation, tests, examples and maintenance contracts, and rejects omissions before rebuilding.
