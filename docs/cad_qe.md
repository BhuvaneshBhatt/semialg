# CAD and quantifier elimination

This page is the practical bridge between the conceptual and reference documentation.

- New to CAD? Read [CAD concepts](concepts/cad.md).
- Want to understand dispatch among presolve, virtual substitution, RUR, and CAD? Read [How semialg chooses an algorithm](concepts/algorithm_selection.md).
- Looking up functions and result types? Use the [Decision and QE reference](reference/decision_and_qe.md) and [CAD reference](reference/cad.md).
- Tuning a slow computation? See the [Performance guide](guides/performance.md).

## Partial CAD and full CAD

Partial CAD interleaves lifting with truth evaluation. An existential branch stops after a true child, and a universal branch stops after a false child. A completed truth decision does not require all ambient cells. This is the principle of [Collins and Hong's partial CAD](https://www.sciencedirect.com/science/article/pii/S0747717108801526).

semialg uses this distinction as follows:

| Request | CAD behavior when specialist methods do not settle it |
| --- | --- |
| Closed real sentence through `resolve`, or automatic `quantifier_eliminate` | Lazy quantifier-guided partial lifting, including alternating prefixes |
| Witness search | Lazy existential search with exact pruning and early success |
| Automatic free-parameter elimination | Complete free-coordinate prefix with lazy quantified lifting and a truth value on every free cell |
| Explicit `strategy="collins"` or `qe_by_complete_cad` | Full CAD construction and truth propagation |

The lazy driver still builds its projection set up front. It can restart with full Lazard projection while retaining partial lifting. Reduced projection and partial lifting describe different choices: a smaller projection family does not by itself imply lazy construction. `LazyResolveResult.status="complete"` means the logical decision is complete. Its traversal statistics do not describe a complete ambient CAD, and a sampled point is not a Skolem strategy for an alternating prefix.

Free-parameter QE uses the same quantifier-guided search separately over every free cell. `strategy="partial"` or `"lazy"` selects this route directly. `LazyQEResult` retains the complete free-cell truth table and visited suffix cells without claiming a complete ambient CAD. See [parameter-cell QE and resource limits](tutorials/parameter_cell_qe.md). Bounded linear-Boolean Fourier–Motzkin elimination runs before polynomial backends in automatic dispatch, including alternating prefixes; excessive expansion or parameter-dependent coefficients leave the problem to the general solver.

<!-- semialg-exec -->
```python
import sympy as sp
from semialg import quantifier_eliminate

x, y = sp.symbols("x y", real=True)
result = quantifier_eliminate(sp.Eq(y**2, x), (("forall", x), ("exists", y)), return_result=True)
assert result.formula is sp.false
assert result.method == "partial_cad_resolve"
assert result.backend_result.stats.stopped_early
condition = quantifier_eliminate(sp.Eq(y**2, x), (("exists", y),))
assert condition == (x >= 0)
```

Presolve retains polynomial Boolean syntax, including nested disjunctions, and eliminated quantified variables do not reappear as free coordinates. CAD constructors require every polynomial symbol in their coordinate tuple and reject explicitly nonreal coefficients. Zero-coordinate CADs accept exact real constants and retain those constant statements for replay.

## First-class quantified expressions

For programmatic formulas, prefer `Exists` and `ForAll` nodes instead of manually encoding quantifier-prefix tuples.

```python
import sympy as sp
from semialg import Exists, ForAll
from semialg.solve import reduce_complete_expr

x, y = sp.symbols("x y", real=True)

formula = ForAll(x, Exists(y, sp.Eq(x + y, 0)))
reduce_complete_expr(formula)
# True
```

Text interfaces remain supported, but programmatic quantifier nodes preserve Symbol identity and compose more safely.

## Virtual substitution before CAD

Quadratic virtual substitution is an exact backend for the supported low-degree fragment. The high-level planner can try it before CAD for QE and witness tasks. If it declines, that is not a mathematical failure; another exact backend can still solve the formula.

Constant polynomial relations remain valid inputs when SymPy has left them unevaluated. For example, conjoining `Eq(0, 0, evaluate=False)` with a quadratic constraint does not prevent elimination. Zero residuals are handled for equality, disequality, strict and weak inequalities, including infinity and infinitesimal substitutions.

## Variable ordering

CAD is highly order-sensitive. Automatic planning preserves quantifier-block semantics while reordering variables only where logically legal. Use `suggest_variable_order` / `suggest_cad_variable_order` for diagnostics instead of guessing an order solely from printed expression size.

## Finite equality varieties

When common polynomial equalities define a finite complex variety, automatic CAD may use the [Gröbner variety CAD](concepts/groebner_variety_cad.md) backend. It follows only compatible algebraic sections and reconstructs the exact finite solution set. Use `return_result=True` to inspect `result.cad.backend` and the `variety_only`, `equality_dimension`, and `quotient_dimension` diagnostics.

```python
import sympy as sp
from semialg import cad

x, y = sp.symbols("x y", real=True)
formula = sp.And(sp.Eq(x - y, 0), sp.Eq(x**2 + y**2, 2), x > 0)

auto = cad(formula, (x, y), return_result=True)
full = cad(formula, (x, y), strategy="collins", return_result=True)

assert auto.cad.backend == "groebner-variety"
assert full.cad.backend == "collins-complete"
assert auto.formula == sp.And(sp.Eq(x, 1), sp.Eq(y, 1))
```

The specialized backend has an exact fallback contract: inability to certify a fiber or algebraic sign causes a complete-CAD retry, not rejection of a candidate section. Pure existential QE may use this specialization; universal and mixed prefixes retain full-space CAD semantics.

Prefix projection through `CADFunction.project` or `project_cell_set` also preserves finite-set semantics. Projecting a finite set returns its exact coordinate points and preserves correlations between retained coordinates. It cannot fill the gaps between those points:

```python
finite = cad(sp.And(sp.Eq(y**2, 1), sp.Eq(x, y)), (x, y), output="function")
projected = finite.project((x,))
assert projected(-1) and projected(1)
assert not projected(0)
```

`qe_from_cad` rebuilds finite-variety inputs from the supplied matrix before elimination. A source variety may have discarded points using residual inequalities, and a later matrix may admit those points again. Rebuilding also supplies the ambient complement required by universal quantifiers. Other incomplete sources are rejected. Finite projected unions carry `CellUnion.finite_variety=True`; their reconstruction keeps the exact coordinate relations between points, even when every listed cell is selected.

```python
from semialg.formula import parse_formula
from semialg.qe.complete import qe_from_cad

equations = sp.And(sp.Eq(y**2, 1), sp.Eq(x, y))
positive_source = cad(sp.And(equations, y > 0), (x, y), return_result=True)
projection = qe_from_cad(positive_source.cad, (x, y), (("exists", y),), parse_formula(equations))
assert bool(projection.subs(x, -1))
assert bool(projection.subs(x, 1))
assert not bool(projection.subs(x, 0))
```

## Reduced projection and fallback

The `collins-complete` entry point starts with all-coefficient McCallum projection: contents, coefficients, discriminants and pairwise resultants of a squarefree basis. Isolated nullifications receive local delineating derivatives. If ordinary lifting cannot discharge its obligations, it restarts with the separate full Lazard backend. This includes positive-dimensional projection nullification; it is not an EC shortcut or a claim to implement Collins reducta/subresultant projection.

Use `cad(..., strategy="lazard-complete")` or `semialg.cad_algorithms.decomp_lazard_complete` to request full Lazard projection directly. That operator retains contents, leading and trailing coefficients, discriminants and pairwise resultants. Native Hasse derivatives compute the lexicographically minimal Taylor valuation over selected algebraic parent samples. The resulting mixed derivative remains a moving polynomial selector, so a section records its residue family instead of a polynomial that vanishes identically on the base.

Recursive formula EC lifting records implied resultant equations and the designated EC at each level. The highest reduced step preserves signs on its EC sections. Intermediate reduced steps additionally project leading coefficients and discriminants of every active polynomial, all non-EC pairwise resultants, and all coefficients of the designated EC. The driver checks the recorded nonzero leading-coefficient conditions on admissible lower EC sections. It can then certify truth by degree/order induction without constructing a reference CAD. A failed induction uses the ambient fallback. A product equation never implies an equation for an individual factor.

For example, `y*z + x - 1` has a degree drop on `y=0`; a nonzero value at one sample cannot justify omitting that stratum. Its exact existential projection is `Ne(y,0) | Eq(x,1)`. Positive Lazard valuations also do not justify accepting a reduced EC proof: the Taylor residue and the original equation have different zeros on a nullifying fiber.

`semialg.cad_algorithms.verify_cad` independently reconstructs projection closure and each native stack over verified parents. It checks local delineating families, root-selector identities and absolute ranks, isolating intervals, coverage, neighboring bounds and exact signs. Full Lazard valuation/residue claims and recursive EC implication/degree/order data are replayed. Completion flags and diagnostic counters are not proof inputs. Supply `formula=` for a recursive truth CAD and `polynomials=` to bind an ambient proof to an external statement. Finite-variety and reduced/TTICAD proof formats remain outside this geometric verifier. The checker reuses the exact arithmetic kernel; it is not a separate formal proof assistant.

See [Native CAD proofs](native_cad_proofs.md) for executable examples, the precise branch-continuation scope and the mapping to the projection theorems in [Algorithm references](references.md#cylindrical-algebraic-decomposition).

## Caching

Exact computation contexts and bounded process-local caches reuse projection, root, sign, specialization, comparison, RUR work, and immutable equality-ideal bases. Formula evaluation reuses the lifting engine's exact sign keys, including symbol identity; scaled products and repeated factors still receive their own exact sign evaluation when necessary. Refinement point location retains the previous cylindrical path to reuse common prefixes without accumulating queried points. Implication minimization shares truth bitsets for identical polynomial boundaries within one exact-computation context. Image certification canonicalizes exact diagonal-quadratic translations, transforming the objective and domain together before caching the complete image; attainment and missing image values are preserved.

## Exactness

Algebraic-coefficient root isolation can split a polynomial into squarefree factors of different multiplicities. Before assigning global root ranks, each factor interval must isolate a single distinct root of the full polynomial. Isolating one factor alone does not prove this: intervals for different factors may overlap. The isolator retains the single-factor fast path and separates mixed-multiplicity factors against their full family.

CAD/QE correctness depends on exact root ordering and sign/truth invariance. Fixed-precision numerical comparisons are not used as hidden substitutes on certified paths. See [Exactness and certification](concepts/exactness_and_certification.md).

Root-family rank certificates use absolute ordered-root indices. A partial list of selected roots does not determine the family's total real-root count; only complete isolation supplies that count. Sign-based root comparisons also check a genuine sign change across nonzero endpoints, since specialization can introduce even multiplicity and a ranked node's left endpoint can be a different root.

## Point location and cylindrical cells

A CAD cell stores one sample point, but its boundaries can move with preceding coordinates. Membership queries specialize those boundaries at the queried base point. This matters even for a half-plane:

```python
import sympy as sp
from semialg import cad

x, y = sp.symbols("x y", real=True)
region = cad(y > x, (x, y), output="function")
assert region(10, 11)
assert not region(10, 9)
assert not region(10, 10)
```

Typed algebraic boundaries use zero-based ranks among the distinct real roots of their defining polynomial. Multiplicity does not create extra sections, and roots from other polynomials in the same stack do not change that rank. A reconstructed `root_of` expression binds its fiber variable; its free symbols are the coefficient parameters, and substitutions preserve that binding.

Certificates used for a typed root's radical presentation, regularity or ordering must match its defining polynomial, fiber variable, root rank, base-variable order, base cell and stack position. Copying a bound with another root identity does not transfer the old certificate. Without a matching proof, the bound can still retain an exact ordered-root expression, but it cannot claim certified regularity or global ordering.

An explicit `base_point` requests a pointwise root comparison. A cell-wide order is not substituted for that request when a coordinate is missing or the requested fiber has no corresponding real root:

```python
from semialg.cad_algorithms.cells import extract_cylindrical_solution

strip = extract_cylindrical_solution(sp.And(x > 0, x * y**2 <= 1), (x, y))
level = strip.full_dimensional_cells[0].levels[1]
comparison = level.typed_lower.compare_certified(level.typed_upper, base_point={x: 1})
assert comparison.certified and comparison.relation == -1
assert comparison.scope == "point"
missing_fiber = level.typed_lower.compare_certified(level.typed_upper, base_point={x: 0})
assert not missing_fiber.certified
```

`extract_cylindrical_solution(..., selected_only=False)` and `extract_structured_cad_cells(..., selected_only=False)` include the ambient cells on which the input formula is false. If the input came from a finite-variety specialization, extraction builds the full ambient CAD before returning those cells. Each cell's `selected` flag records whether it belongs to the solution.

Solution extraction rejects an incomplete CAD result. Coverage certificates require a completed source CAD and all of its selected cells or all of its ambient cells, with no omissions or duplicates. Structured cells supplied without that provenance can be rendered and checked individually, but do not establish coverage of a formula.

Virtual-substitution prepasses handle uniform existential or universal prefixes. They leave alternating prefixes to CAD, including when partial elimination is requested, so an existential variable is never moved across a universal binder.

## Degenerate closure boundaries

The closure of a cylindrical cell may require more than replacing strict bounds by weak bounds. For example, the surface `x > 0, y = x*z` has the entire `z` axis in its closure at `x = y = 0`; substituting that base point into `z = y/x` loses the limiting fiber. CAD topology raises `ExactComputationUnavailable` when a vanishing boundary denominator, an unresolved multiple-root branch or an exact incidence evaluation requires further refinement. Such a failure does not establish non-incidence, emptiness or a completed closure computation. Regular closure cases and exterior computations continue to retain internal CAD sections whenever those sections belong to the resulting region.

For reconstruction calls with `closed=True`, merged interval blocks honor closed outer bounds just as single-cell reconstruction does. This is exact for regular intervals and boxes. The rendering option does not resolve the degenerate closure cases described above; use the topology operations when computing a Euclidean closure.

`CADFunction.project` and `project_cell_set` support nonempty coordinate prefixes. Use existential quantifier elimination for projection onto no coordinates, which is a feasibility question; the CAD methods raise `NotImplementedError` for this case.

## Reusing CAD regions safely

A Boolean combination can reuse a CAD only when its truth is constant on each existing cell. A callback that adds a boundary, such as adding `x < -1` to a region for `x > 0`, requires refinement:

```python
import sympy as sp
from semialg import as_cad_region
from semialg.cad_region import cad_combine

x = sp.Symbol("x", real=True)
positive = as_cad_region(x > 0, (x,))
combined = cad_combine(positive, op=lambda condition: sp.Or(condition, x < -1))
assert combined.result.as_function()(-2)
assert not combined.result.as_function()(-1)
```

A finite-variety CAD can be reused for combinations that remain false outside its equality variety, including unions and intersections of its selected cells. A complement needs an ambient decomposition. Symmetric difference constructs a common refinement when the operands do not already share one.

QE checks the scope of the source CAD's invariance proof before reusing its samples. A full Collins CAD certifies signs for its polynomial factors and projection factors, so nonzero scalar multiples, powers and Boolean combinations of those signs can reuse it. A reduced truth certificate applies to its recorded formula and its complement; it does not establish arbitrary sign conditions away from an equational constraint. Raw reduced projection towers remain incomplete until a safe driver certifies them.

When a new query lies outside the proved scope, `qe_from_cad` builds a decomposition for that matrix and reports the backend it actually used. For example, a CAD for `x = 0` does not include the boundary at `x = 1` needed here:

```python
from semialg.cad_algorithms.decomposition import decomp_collins_complete
from semialg.formula import parse_formula
from semialg.qe.complete import qe_from_cad

x = sp.Symbol("x", real=True)
source = decomp_collins_complete((x,), (x,))
answer = qe_from_cad(source, (x,), (("exists", x),), parse_formula(sp.And(x > 0, x < 1)))
assert answer is sp.true
```

## Zero-coordinate CADs

The real space with no coordinates has one point, the empty tuple. Its CAD contains one ambient cell; a constant true formula selects it and a false formula selects no cells. Relifting, sign verification, point location, complete QE and full coordinate restriction preserve this convention. Constant QE results record their selected cells and sample assignments just as nonconstant results do.

```python
from semialg import as_cad_region, cad

unit = cad(sp.true, (), return_result=True)
assert unit.as_function()()
location = as_cad_region(unit).locate_point(())
assert location.selected and location.dimension == 0
assert location.sign_vector == ()
assert unit.cad.verify_sign_invariance().ok
```

## Pruning and witness semantics

Lazy CAD pruning preserves an undecided algebraic equality as unknown. A failed simplification is insufficient to reject an equality or accept a disequality; lifting and exact sign tables can still decide the leaf. Constant sentences and zero-coordinate instance searches return their Boolean truth directly. If an existential search short-circuits on a prefix, its returned assignment supplies values for all requested coordinates.

`witness_from_qe_result` samples the free-coordinate solution set of a completed QE result. It preserves supplied coordinates, completes missing free coordinates, and identifies the source cell at the returned point. It does not reconstruct bound variables or a strategy for alternating quantifiers. `existential_wit_assign` checks supplied candidate samples after fixing the requested coordinates; it does not relift a new fiber and does not turn failure to find a candidate into an infeasibility proof.

A virtual-substitution witness callback returns an assignment or reports that it found none. The latter gives `status="unknown"`. A ground reduced formula that is exactly false establishes `status="unsat"`; a true ground formula needs no coordinate search. High-level instance finding preserves these statuses and can fall back to CAD when witness reconstruction is unresolved.

Structured and text instance queries conjoin explicit assumptions with the matrix before dispatch to RUR, virtual substitution or lazy CAD. Symbols introduced only by those assumptions are included in the returned coordinates. A request with `count <= 0` performs no search and returns an empty collection with `status="unknown"`; it supplies no infeasibility conclusion.

## Public resource controls

`cad` and `quantifier_eliminate` accept `max_cells` and `timeout`. A shared cooperative budget includes preprocessing, failed attempts, intermediate levels, cache returns and output construction. Structured non-strict CAD returns unknown when exhausted; direct outputs, strict CAD and QE raise `ResourceLimitError`. Read the [resource tutorial](tutorials/parameter_cell_qe.md) before relying on deadlines.
