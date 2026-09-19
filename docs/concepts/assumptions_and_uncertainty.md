# Exactness, assumptions, and uncertainty

Symbolic expressions can look similar while supporting different mathematical conclusions. semialg distinguishes these cases explicitly.

## Syntactic simplification

A simplifier rewrites an expression. Failure to rewrite an expression to `0` is not proof that it is nonzero. Recognizers therefore use certified equality and zero predicates when a branch depends on a mathematical decision.

## Certified identity

A certified identity establishes that two exact expressions denote the same algebraic object in the stated setting. This is stronger than matching expression trees and is why equivalent polynomial presentations should not change correctness.

## Pointwise truth and parameters

A symbolic expression can be nonzero as a polynomial yet vanish for particular parameter values. Parameter-dependent problems may therefore require conditions or strata rather than a single unconditional Boolean answer.

## Assumptions

Assumptions are inputs to a decision problem, not facts that semialg silently adds to its mathematical output. When an algorithm needs a sign, nonzeroness, or reality condition, that condition must be certified from exact data and the supplied assumptions.

## Failure to certify

`None`, an unsupported result, or an incomplete structured result means that the selected certified method did not establish the requested claim. It must not be interpreted as `False`. See the limitations and failure-mode guides for operation-specific behavior.
