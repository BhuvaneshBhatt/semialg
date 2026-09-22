# Certified decision boundaries

Symbolic mathematical decisions use shared certification helpers rather than Python truth tests. General expression equality and sign questions pass through the package's zero/equality/sign boundary; exact algebraic properties such as `Poly.is_zero` remain direct because they are properties of an exact polynomial representation.

Geometry validation follows the same rule. A parameter-dependent sign or equality is preserved as unresolved when it cannot be certified. Parameterized ellipsoids, for example, retain unresolved symmetry and positive-definiteness conditions instead of treating uncertainty as failure.

Boolean coercion is reserved for ordinary Python state. Mathematical APIs preserve three-valued outcomes where uncertainty is part of the contract. Backend fallbacks catch only failures that mean the backend cannot handle the input; malformed internal state and invariant violations are allowed to surface as errors.

`StandardRegion` is the shared protocol for explicit structured geometry. It currently includes composite representations because measurement, sampling, topology, transformations, integration, and structural operations dispatch on that protocol. Changing that hierarchy requires a replacement structured-geometry protocol rather than a class-only rename.
