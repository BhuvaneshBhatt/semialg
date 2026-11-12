"""Independent sign-cell expectations for exact, factored univariate QE inputs."""

from __future__ import annotations

import argparse
import json
import random
from dataclasses import dataclass
from functools import cmp_to_key
from pathlib import Path

import sympy as sp


@dataclass(frozen=True)
class OracleCase:
    """A formula and truth values obtained from its known real zero set."""

    index: int
    variable: sp.Symbol
    polynomials: tuple
    formula: sp.Basic
    exists: bool
    forall: bool


def make_cases(count=28, seed=28419):
    """Partition the line at supplied exact roots, independently of any QE backend."""
    variable = sp.Symbol("x", real=True)
    rng = random.Random(seed)
    root_sets = [
        (),
        (-2,),
        (-1, 0),
        (0, 1),
        (-2, 0, 2),
        (-sp.sqrt(2), sp.sqrt(2)),
        (sp.sqrt(2) - 1, sp.sqrt(2) + 1),
    ]

    def compare(a, b):
        return int(sp.sign(sp.simplify(a - b)))

    for index in range(count):
        roots = root_sets[index % len(root_sets)]
        other = root_sets[(index * 3 + 1) % len(root_sets)]
        first = (-1 if index % 2 else 1) * sp.prod(
            (variable - r) ** rng.choice((1, 2)) for r in roots
        )
        second = sp.prod(variable - r for r in other)
        left = (sp.Lt, sp.Le, sp.Eq, sp.Ne)[index % 4](first, 0)
        right = (sp.Ge, sp.Gt, sp.Ne, sp.Eq)[(index // 4) % 4](second, 0)
        formula = (sp.And if index % 3 else sp.Or)(left, right)
        cuts = sorted(set(roots + other), key=cmp_to_key(compare))
        points = list(cuts)
        if cuts:
            points.extend((cuts[0] - 1, cuts[-1] + 1))
            points.extend((a + b) / 2 for a, b in zip(cuts, cuts[1:], strict=False))
        else:
            points.append(sp.S.Zero)
        values = [bool(sp.simplify(formula.subs(variable, point))) for point in points]
        yield OracleCase(index, variable, (first, second), formula, any(values), all(values))


def check_case(case):
    """Compare automatic, lazy, full Collins and full Lazard decisions to the oracle."""
    from semialg import quantifier_eliminate
    from semialg.cad_algorithms import decomp_lazard_complete
    from semialg.formula import parse_formula
    from semialg.partial.qe import lazy_resolve_formula
    from semialg.qe.complete import qe_by_complete_cad, qe_from_cad

    matrix = parse_formula(case.formula)
    variables = (case.variable,)
    cad = decomp_lazard_complete(case.polynomials, variables)
    rows = []
    for quantifier in ("exists", "forall"):
        expected = getattr(case, quantifier)
        prefix = ((quantifier, case.variable),)
        actual = {
            "automatic": bool(quantifier_eliminate(case.formula, prefix)),
            "lazy": lazy_resolve_formula(variables, iter(prefix), matrix).truth_value,
            "full": qe_by_complete_cad(
                variables,
                prefix,
                matrix,
                use_presolve=False,
                allow_variety_cad=False,
                return_result=True,
            ).truth_value,
            "lazard": bool(qe_from_cad(cad, variables, prefix, matrix)),
        }
        if any(value != expected for value in actual.values()):
            raise AssertionError((case.index, quantifier, case.formula, expected, actual))
        rows.append(
            {
                "case": case.index,
                "quantifier": quantifier,
                "formula": str(case.formula),
                "expected": expected,
                "results": actual,
            }
        )
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--generated", type=int, default=28)
    parser.add_argument("--seed", type=int, default=28419)
    parser.add_argument("--output", type=Path, default=Path("qe-oracle-results.json"))
    args = parser.parse_args()
    if args.generated < 1:
        parser.error("--generated must be positive")
    records = [row for case in make_cases(args.generated, args.seed) for row in check_case(case)]
    args.output.write_text(
        json.dumps({"seed": args.seed, "cases": records}, indent=2) + "\n", encoding="utf-8"
    )
    print(f"{len(records)} oracle decisions; {4 * len(records)} backend comparisons passed")


if __name__ == "__main__":
    main()
