# Benchmark and conformance suite

`semialg.benchmarks` separates **mathematical conformance** from **performance measurement**. A benchmark result is useful only when the result is known to be correct; wall-clock time is therefore reported alongside deterministic structural counters rather than used as the correctness criterion.

## Established benchmark sources

There is no single universal benchmark suite covering all of real algebraic geometry. The strongest reusable collections differ by algorithm.

| Area | Established sources/families | What semialg records |
| --- | --- | --- |
| CAD | Wilson's *CAD Example Bank*; classical Collins/McCallum, Davenport--Heintz, Buchberger--Hong and later EC/TTICAD examples | cells, projection/lifting work where exposed, ordering, well-orientedness/certification, CPU and wall time |
| Quantifier elimination | CAD/QEPCAD literature examples and textbook decision/QE problems | quantifier structure, selected backend, cells visited/constructed, output size, certificate/equivalence checks, time |
| Roadmaps/topology | examples and constructions in Basu--Pollack--Roy and roadmap papers | component count, critical points, recursive subproblems, roadmap size, certificate checks, CAD work, time |
| Component instances | no widely adopted independent benchmark bank; use topology literature plus explicit conformance families | exact sample per component, component count/dimension, mixed-dimensional and singular cases |
| Gröbner bases | Cyclic, Katsura, ECO, Noon and related standard polynomial-system families | basis size/degrees, reductions/certificates where exposed, time |
| Real roots / RUR / regular chains / decomposition | standard polynomial-system families and published algorithm examples; no one cross-system corpus dominates | root/component counts, degrees, isolating data/certificates, algebraic work counters, time |

The catalog includes the complete Wilson CAD Example Bank v4, the associated TTICAD corpus, and parameterized Cyclic, Katsura, ECO, and Noon Gröbner families. Imported cases preserve source identifiers and mathematical formulations rather than copying another system's executable benchmark code.

## Case schema

Every `BenchmarkCase` records a domain, variables, formula or polynomial system, provenance, tags, difficulty metadata, and known expected results where available. `run_with()` measures both CPU and wall time and accepts algorithm-specific correctness checkers and structural-counter extractors. `run_groebner()` demonstrates the pattern by certifying that the computed basis reduces every input generator to zero.

```python
from semialg.benchmarks import all_cases, run_groebner

case = next(case for case in all_cases() if case.name == "groebner_cyclic_3")
result = run_groebner(case)
assert result.correct
print(result.metrics)
```

## Interpreting performance

Do not compare raw timings from different machines as if they were algorithmic invariants. Prefer cell counts, projection counts, basis size/degree, roots isolated, recursive subproblems, or other deterministic work measures. CPU and wall time are still retained for same-machine regression tracking.

## Provenance

The converted CAD examples cite D. Wilson, *Real Geometry and Connectedness via Triangular Description: CAD Example Bank*, DOI `10.15125/BATH-00069` (CC BY-SA 4.0). General real-algebraic-geometry conformance cases cite S. Basu, R. Pollack, and M.-F. Roy, *Algorithms in Real Algebraic Geometry*, 2nd edition.

### Published CAD corpora

`cad_cases()` contains all 68 entries of David Wilson's CAD Example Bank v4. Each
record preserves the source procedure/index identifier and the variable order in
the bank. The bank does not itself attach a canonical cell count to every entry,
so a missing reference count is not represented as zero or inferred from another
implementation.

`tticad_cases()` contains all 29 examples from the published Section 8.2 TTICAD
dataset (DOI 10.15125/BATH-00076). Formula grouping and equational constraints
are retained in `expected["formulae"]`, in addition to the flattened polynomial
list used by generic CAD runners. Reference cell counts belong to a particular
algorithm/configuration and should be stored with that algorithm label rather
than treated as an intrinsic property of the mathematical problem.

## Running the benchmarks

The benchmark package has a command-line runner. Run individual cases while developing; the complete published corpora can be expensive.

```bash
python -m semialg.benchmarks cad --case CMXYExamples:1
python -m semialg.benchmarks tticad --case Section82:1:IntA
python -m semialg.benchmarks groebner cyclic 5
python -m semialg.benchmarks groebner katsura 5
```

Omit `--case` to run the complete Wilson CAD bank or the complete TTICAD corpus:

```bash
python -m semialg.benchmarks cad
python -m semialg.benchmarks tticad
```

`run_cad()` constructs a Collins sign-invariant CAD in the prescribed source variable order and records full-dimensional cell count, counts by level, projection-polynomial counts, maximum stack size, exact-root isolation/comparison/cache counters, specialization-cache counters, CPU time, and wall time. Published cell counts are compared only when the stored reference names the same configuration; counts from reduced EC, TTICAD, QEPCAD, or other algorithms are not treated as Collins reference values.

`run_tticad()` is separate. It preserves the published formula-family and equational-constraint structure and runs semialg's certified family-aware TTICAD path. Its metrics report the effective backend and whether certification required a Collins fallback. The `published_tticad` reference count is the TTICAD column of Table 2 in Bradford et al.; it is a published reference measurement, not an expected equality for semialg's independently implemented algorithm.

The Gröbner runner accepts parameterized `cyclic`, `katsura`, `eco`, and `noon` families. For example, `groebner cyclic 7` constructs the size-7 family at run time rather than requiring a separate checked-in case for every size. Gröbner correctness is checked by reducing every input generator with the computed basis.
