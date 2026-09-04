# Exception hierarchy

All package-specific failures derive from `SemialgError`. The hierarchy distinguishes invalid input or unsupported exact fragments from failures of a particular exact strategy.

```text
SemialgError
├── CertificationFailure
├── DimensionMismatchError
├── FormulaNormalizationError
└── SemialgStrategyFailure
    ├── AlgebraicSolvingError
    ├── BackendFailure
    ├── ExactEvaluationFailure
    ├── QuantifierEliminationError
    ├── ReconstructionFailure
    ├── ResourceLimitError
    └── UnsupportedFragmentError
```

Some classes also inherit `ValueError` so ordinary Python/SymPy input-validation handlers can catch them naturally.

## `SemialgError`

Base class for package-specific errors. Catch this only when a caller genuinely wants to treat every `semialg` failure uniformly.

## `FormulaNormalizationError`

The input formula, variable list, bounds, or related structural data could not be normalized into the required semialgebraic representation.

## `DimensionMismatchError`

A point, coordinate list, mapping, matrix, or region dimension is inconsistent with the declared ambient variables.

## `UnsupportedFragmentError`

The requested exact strategy does not support the mathematical fragment supplied. This is different from a proof that the problem has no solution.

High-level planners may catch this internally and try another exact backend.

## `SemialgStrategyFailure`

Base class for failures of an exact computational strategy after the input has been accepted.

### `AlgebraicSolvingError`

Exact algebraic solving/RUR/root machinery could not complete the required operation.

### `BackendFailure`

A selected computational backend failed in a way that permits a higher-level caller to consider fallback.

### `ExactEvaluationFailure`

An exact sign, comparison, specialization, or evaluation could not be certified by the implemented method. The package should not replace this with a hidden finite-precision guess.

### `QuantifierEliminationError`

An exact QE operation failed to complete under the selected strategy.

### `ReconstructionFailure`

The internal exact decomposition was available but could not be reconstructed into the requested public representation.

### `ResourceLimitError`

A configured or implementation resource limit was reached. It does not imply that the mathematical problem is unsupported in principle.

## `CertificationFailure`

A candidate/result was produced, but the requested certification obligation could not be established.

## Catching exceptions

Prefer catching the narrowest class meaningful to the application:

```python
from semialg import ResourceLimitError, UnsupportedFragmentError

try:
    ...
except UnsupportedFragmentError:
    # choose a different exact formulation/backend
    ...
except ResourceLimitError:
    # simplify the problem or adjust resource policy
    ...
```

For the distinction between unsupported input, strategy failure, and conservative exact failure, see [Errors and failure modes](../guides/errors_and_failure_modes.md).
