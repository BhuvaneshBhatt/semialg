# Quantifier-elimination architecture

`quantifier_eliminate` is the common entry point for certified real quantifier elimination. It reuses the package's existing exact engines rather than duplicating their algebra.

The dispatcher normalizes the quantified problem, applies exact affine and equality-ideal preprocessing, detects independent variable-incidence blocks, and selects the narrowest certified backend that covers the resulting structure. Linear problems can use Fourier--Motzkin elimination; supported low-degree formulas can use virtual substitution; zero-dimensional equality structure can use Gröbner/RUR machinery; reduced or complete CAD handles the remaining supported polynomial formulas.

Independent conjunctive blocks are eliminated separately when their quantified variables are disjoint. Equality-ideal Gröbner consequences may be added as equivalent preprocessing information, but an elimination ideal is never treated as a real existential image: Zariski projection and real projection have different semantics.

Reusable parameterized CAD is represented by `ParametricCADFunction`. It stores certified cases, supports exact specialization through `specialize()`, and can recover the represented formula with `normal_form()`. Unresolved parameter conditions remain unresolved rather than being coerced through Python truth testing.

The dispatcher keeps backend roles separate: Gröbner bases provide exact equality-ideal consequences, RUR describes finite real algebraic solution geometry, and Fourier--Motzkin, virtual substitution, and CAD perform certified real elimination.
