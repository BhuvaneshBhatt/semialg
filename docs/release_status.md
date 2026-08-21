# Release status

The current package version is **0.2.0b1**, a beta pre-release of semialg 0.2.0.

This release is suitable for evaluation, testing, experimentation, and development against the documented APIs. It should not yet be treated as a final compatibility commitment. Solver coverage, heuristics, result metadata, and performance may change before 0.2.0 final.

## Correctness and certification

semialg distinguishes exact/certified results from incomplete, subset, superset, witness-only, and bounded-window results. Callers should inspect the documented result semantics rather than infer completeness from the presence of a formula or witness.

A beta designation does not weaken a result that is explicitly reported as exact or certified: those labels retain their documented meaning. The beta designation instead reflects that the package API and implementation are still being hardened and may receive corrections before the final release.

For critical applications, independently validate results and report reproducible cases that appear inconsistent with the documented contracts.

## Versioning

The beta uses the PEP 440 version **0.2.0b1**. A future final release will use **0.2.0** after beta feedback and release checks are complete.

When installing from PyPI, pre-releases are normally excluded unless explicitly requested. This helps prevent users who request ordinary stable releases from receiving the beta unintentionally.
