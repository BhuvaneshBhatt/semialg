# Development tools

Repository-maintenance utilities belong here when they are needed. Runtime code does not depend on this directory.

## Isolated slow-test runner

`python tools/run_slow_tests.py` collects tests marked `slow` and executes each
node in its own Python process, one at a time. This prevents a hang, leaked
state, or abnormal process exit in one expensive test from obscuring the rest
of the suite. Results are written incrementally to `.test-reports/slow/` as
JSON, Markdown, and one log per test.

Use `--soft-timeout` to flag performance outliers and `--hard-timeout` to stop
a pathological test and continue with the next node. Collection paths can be
supplied positionally, for example:

```bash
python tools/run_slow_tests.py tests/test_local_geometry_property_contracts.py \
  --soft-timeout 20 --hard-timeout 90
```
