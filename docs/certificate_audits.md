# Certificate-field audits

`tests/certificate_fields.json` classifies the stored fields of each registered certificate family. New or renamed fields fail the audit until their role is recorded. No field is skipped merely because it is named method, source or modulus.

Proof fields have explicit mutations and replay must reject them. Non-proof fields require a reason. semialg additionally records conditional branches as inactive when absent, constructor/schema rejections separately from exact replay, and per-root overrides when a shared component type plays a different role in another certificate. Unexpected verifier crashes fail the tests.

Set `CERTIFICATE_AUDIT_REPORT=certificate-audit.jsonl` during pytest to retain field paths and outcomes. The report distinguishes mutated evidence, provenance, diagnostics and inactive branches. This is semantic certificate mutation coverage, not a source-code mutation score. The registry covers the exercised families; it does not establish that every possible malformed Python object or every certificate family in the package has been enumerated.

Independent arithmetic identities, coordinate transports and prescribed fibers complement replay tests so a constructor and verifier cannot serve as each other's only oracle. See the package testing guide for the registered workflows.

Native CAD tests mutate projection factors/provenance, EC implications and guards, Taylor valuations/residues, root ranks, cached expressions, selector contexts, isolating intervals, signs, bounds and stack coverage. Branch tests mutate both multiplicity spectra and Sturm continuation ranks. Statement coordinates and source-cell assumptions are classified as context; the caller must bind them to the query being proved. Replay outcome counters are derived data, not certificates.

## Aggregate branch coverage

The CI coverage job writes `certificate-audit.jsonl` and runs `tools/check_certificate_audit.py`. It merges fixtures by the declaring class and field, so an optional branch inactive in one test must become active in another. Each exercised proof field must have a rejection witness, and each exercised root family must include a well-formed mutation rejected by replay. Constructor and schema rejections are counted separately and cannot satisfy the family-level mathematical-replay requirement. The summary retains per-field activation counts and omissions.

`tests/certificate_audit_inventory.json` requires the maintained families and proof fields, so deleting all fixtures for a family also fails the aggregate gate. This inventory covers the exercised mutation families, not every dataclass in the package. `tests/certificate_fields.json` requires explicit roles for every traversed dataclass field and supports context-specific nested policies. Adding a new field without a policy fails immediately. The audit does not claim exhaustive resistance to arbitrary malicious certificates or mathematical independence when the verifier recomputes with construction code.
