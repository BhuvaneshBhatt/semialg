# Understanding result objects

Many `semialg` functions return more than a bare SymPy expression because exact computation often has several logically distinct outputs: a value, a witness, an attainment statement, a certificate, or parameter guards.

## Why not always return a SymPy expression?

Consider

\[
\inf_{-1<x<1}x=-1.
\]

The exact value is `-1`, but no feasible point attains it. A useful optimization result must therefore distinguish:

- the optimum/infimum value;
- whether it is attained;
- optimizer witnesses when they exist;
- whether global optimality has been certified.

Likewise, a parameter-dependent answer may have different exact values on different semialgebraic regions of parameter space.

## Common result concepts

### `value`

The mathematical value computed by the operation. For optimization this may be an infimum or supremum even when it is not attained.

### `attained`

Whether a feasible point realizes the optimum value.

### witness / optimizer points

Exact points demonstrating feasibility or attainment. A finite witness list is not the same thing as the complete optimizer set; use `argmin_set` or `argmax_set` when the entire locus matters.

### `certified`

Where exposed, this means the requested global mathematical claim has been established by an exact procedure. An exact candidate point is not automatically a certified global optimum.

### formula and quantifiers

Some exact parametric answers are naturally first-order relations. A result can therefore contain an exact formula plus explicit quantifiers instead of prematurely forcing another expensive QE pass.

### guarded branches / strata

A parameter-stratified result represents

\[
C_1(p)\Rightarrow v_1(p),\quad
C_2(p)\Rightarrow v_2(p),\ldots
\]

where each \(C_i\) is an exact semialgebraic guard.

## Selecting a parameter branch

Stratified results support selection by exact parameter assignment where applicable:

```python
result.select({a: 4})
```

String keys may also be accepted by APIs that resolve them against the result's actual parameter symbols. Exact Symbol keys are preferable when same-name symbols with different assumptions are in scope.

## `Piecewise` versus a stratified result

A SymPy `Piecewise` is convenient when every branch can be expressed as a simple symbolic value. A stratified result carries more information:

- explicit semialgebraic guards;
- branch-level values/results;
- coverage and disjointness information;
- certification metadata;
- potentially quantified exact relations.

Use `.as_piecewise()` only when the result type supports it and a plain symbolic expression is the representation you need.


## Metadata is certified when present

Cheap metadata fields are intentionally conservative. In particular, a numeric
`dimension` value means that the dimension has been established exactly by the
metadata path (for example by affine rank or a certified cylindrical
decomposition). `dimension=None` means that this inexpensive layer did not
certify a single global dimension; it does **not** mean that the set lacks a
dimension.

Parameter-dependent affine systems illustrate why this matters. For `a*x = 0`
the fiber in `x` has dimension 0 when `a != 0` and dimension 1 when `a = 0`.
Without first stratifying parameter space, reporting either integer globally
would be wrong, so cheap metadata correctly leaves the field unknown.

## CAD and algebraic result objects

CAD and exact algebraic solving expose structured objects because cell bounds, root identity, isolating data, and provenance matter for subsequent certified operations. Avoid converting these objects to floating point merely for convenience if they will feed another exact computation.


## Empty, partial, unsupported, and failed are different

These states should not be conflated:

- **empty exact result** — the requested set is mathematically empty and that has been certified;
- **partial result** — useful exact information was produced, but the API explicitly records that the requested representation/computation is incomplete;
- **unsupported fragment** — the implemented method does not cover the input;
- **certification failure** — a candidate or intermediate object exists, but the required exact proof could not be established;
- **resource limit** — a supported strategy was stopped by an explicit computational limit;
- **implementation error** — an unexpected programming failure, which should propagate rather than masquerade as mathematical emptiness.

When an API exposes fields such as `partial`, `complete`, `certified`, or `status`, consume those fields rather than inferring success merely from a nonempty symbolic object.

## Practical checklist

When consuming a structured result, ask:

1. What mathematical object does `value` or `formula` denote?
2. Is the optimum/value attained?
3. Is a returned point merely a witness, or the complete solution locus?
4. Is the global conclusion certified?
5. Does the answer depend on a parameter guard?
6. Is the relation still quantified?

For the precise distinction between exact representation and proof, see [Exactness and certification](exactness_and_certification.md).
