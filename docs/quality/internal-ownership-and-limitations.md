# Internal ownership and practical limitations

This page records permanent internal ownership rules discovered in the maintained architecture.

## Exact cylindrical samples

`topology.sample_context.CylindricalSampleContext` owns dependency-aware closure of CAD sample towers. A coordinate may depend only on earlier cylindrical coordinates. Consumers should request a prefix assignment rather than independently substituting nested `root_of` expressions.

`NativeRegionAnalysis` owns reusable region CAD state. Expensive closed path formulas and exact sample assignments used by connectivity are cached with that analysis rather than reconstructed by topology callers.

## Algebraic modules

The large algebraic modules are intentionally split by mathematical responsibility rather than file size: `samples` owns sample representation/conversion, `comparison` owns exact ordering/equality, `roots` owns isolation/refinement, and `signs` owns certified sign determination. `algebraic_decomposition` remains the orchestration/certification layer for ideal and regular-chain decompositions. Cross-module helpers should move only when they have a single lower-level owner; splitting orchestration mechanically would increase call indirection without removing mathematical coupling.

## Exact strategy outcomes

Specialized exact algorithms distinguish success, inapplicability, and inability to certify. Internal dispatch code should not encode those three states as broad exceptions or collapse `UNKNOWN` to false. Programming errors continue to propagate.

## Function graphs

The exact graph fragment includes rational expressions and exact rational powers (including explicit real-root forms), plus the documented finite piecewise and wrapper operations. General exponential, logarithmic, trigonometric, irrational-power, and variable-exponent graphs are not semialgebraic and are declined rather than approximated.

## Topology

Known standard/polyhedral regions use their certified finite triangulation before general CAD/roadmap topology. General exact topology remains substantially more expensive. Native CAD connectivity reuses cached CAD paths and sample assignments.

## Integration

Automatic integration reports whether numerical evaluation was requested or was used because symbolic iterated integration remained unevaluated. Exact integration is intentionally incomplete: an integral of a semialgebraic function need not itself be algebraic or elementary.

## Resource diagnostics

CAD result diagnostics report the requested strategy, effective backend, cell/projection counts, and resource-limit state. `max_cells` and wall-clock `timeout` are currently not backend-enforced; requesting them therefore produces an explicit unsupported result rather than a false guarantee.
