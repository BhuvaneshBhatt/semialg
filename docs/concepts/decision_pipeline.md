# How semialg decides a problem

`semialg` is designed around exact mathematical conclusions rather than a single universal algorithm. A typical operation follows this progression:

1. **Normalize the input.** Polynomial relations, variables, assumptions, and Boolean structure are put into stable internal forms without changing their meaning.
2. **Recognize exact structure.** Cheap recognizers detect supported linear, low-degree, geometric, finite, or otherwise structured cases. Algebraically equivalent presentations are expected to select equivalent mathematical paths.
3. **Use a specialized exact algorithm.** When its hypotheses are certified, a specialized solver can avoid the cost of general CAD.
4. **Certify the result.** Candidate witnesses, signs, roots, decompositions, bounds, and identities are checked by exact predicates. Heuristics may choose work; they do not establish truth.
5. **Fall back to general exact machinery.** Supported problems that are not handled structurally can use algebraic elimination, quantifier elimination, or CAD.

A specialized method declining an input is not a negative mathematical result. It means that another certified path must be used or that the requested operation is outside the implemented fragment.

## Why use this architecture?

General CAD can be expensive even for problems with simple structure. Recognizing that structure first gives common problems much cheaper exact solutions while preserving a general certified route for the real-polynomial fragments that semialg supports. This also keeps numerical approximation separate from mathematical certification.
