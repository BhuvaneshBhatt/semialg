# Coverage measurements

Coverage is measured to find untested behavior, not to reward tests that merely execute
lines. Exact results, failure modes, witnesses, endpoint semantics, parameter degree
drops, and conservative refusal branches are higher-value coverage than artificial calls
whose only assertion is that no exception occurred.

## What is measured

CI runs the complete non-slow suite with branch measurement and writes JSON, XML, HTML,
and a Markdown subsystem summary under `build/coverage/`. The policy checker computes
three percentages from the same raw counts:

- **statement coverage** = covered statements / statements;
- **branch coverage** = covered branches / branches;
- **combined coverage** = (covered statements + covered branches) / (statements + branches).

The existing `coverage.py` `fail_under = 74.5` gate is the combined percentage when
branch measurement is enabled. A reported aggregate percentage should not be treated as a branch-only
figure unless it came from an explicitly branch-only calculation.

## Threshold progression

The current release tree retains the **74.5% combined gate**. Raising a number
without a fresh full report would create false precision. The intended progression is:

| Checkpoint | Combined target | Evidence expected before adopting it |
| --- | ---: | --- |
| Current | 74.5% | existing behavioral/property/failure-path suite |
| Decision/parameter/root expansion | 76% | critical decision branches, parameter strata, root degree drops and failures |
| Topology/optimization/integration matrices | 78% | exact boundary/attainment/measure behavior and important unsupported paths |
| Remaining thin API families | 80% | independent behavioral evidence rather than execution-only calls |

A target becomes an enforced floor only after the full non-slow report on the exact tree
is above it with reasonable margin. Removing obsolete unreachable code is preferable to
writing tests solely to cover dead branches.

## Subsystem ownership

`tests/coverage_policy.toml` partitions all package source files into thirteen reviewable
subsystems. The checker rejects overlap and rejects an unassigned file, so a new nested
subsystem cannot silently fall into a large catch-all group. Shared top-level plumbing is
the only deliberate fallback group.

The initial enforced combined floors are **85%** for public decision/parameter/root
layers and **80%** for algebraic certification/comparison and high-level
geometry/topology. These are the high-value layers for which the current test inventory
already provides dense behavioral evidence. Other subsystems are reported from the first
fresh baseline before a local floor is committed.

This prevents the three most important small subsystems from being hidden by a large
well-tested module while avoiding invented floors for areas that have not yet received a
fresh subsystem measurement.

## Running the checks

```sh
mkdir -p build/coverage
pytest -m "not slow" \
  --cov=semialg --cov-branch \
  --cov-report=xml:build/coverage/coverage.xml \
  --cov-report=json:build/coverage/coverage.json \
  --cov-report=html:build/coverage/html
python scripts/check_coverage.py \
  build/coverage/coverage.json \
  --markdown build/coverage/summary.md
```

The checker fails if the repository-wide enforced floor regresses, if a configured
subsystem floor regresses, if branch data is absent, or if the report is incomplete.
The CI artifact upload uses `if: always()` so available coverage diagnostics survive a
failed test or gate.
