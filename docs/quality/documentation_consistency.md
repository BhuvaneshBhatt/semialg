# Documentation consistency

The documentation is checked against the live public API and executable examples. The checks keep names, signatures, examples, navigation, and mathematical contracts synchronized with the implementation.

## What is checked

- README examples and links use the public API and remain valid when rendered outside the repository.
- Local-geometry documentation covers closure, germs, sign strata, components, curve selection, local bounds, parameter strata, local maps, path independence, and weighted blow-up charts.
- Gallery Markdown and Python examples stay synchronized and executable.
- Public functions appear in the API reference with usable signatures and behavioral examples.
- Exact computations are distinguished from numerical approximations and heuristic planning.
- Source distributions contain the documentation and examples required by the published navigation.

## Mathematical test coverage

The local-geometry property corpus compares incremental sign strata with an exhaustive small-system reference, checks permutation behavior, exercises symbol sign assumptions and multiple witnesses, verifies parameter-stratum disjointness, and covers exact preimage-fiber outcomes. Topology cases include nodes, cusps, tacnodes, and coordinate transformations because these stress CAD and algebraic-sign behavior.

The algebraic-core suite also checks that a cached CAD sector sample remains strictly between refined exact bounds or is replaced.

## Build checks

Normal validation checks relative links, navigation, executable examples, public-API documentation coverage, source-distribution contents, and formatting. A separate installed-wheel test executes representative public examples against the installed package. Strict site construction is part of the documentation build when MkDocs is available.
