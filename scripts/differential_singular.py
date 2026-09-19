"""Optional exact-algebra differential corpus against Singular ``primdecGTZ``.

Singular is an independent oracle for development only.  Every semialg result
must first pass its own exact replay; oracle agreement never substitutes for a
semialg certificate.
"""

from __future__ import annotations

import argparse
import random
import re
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path

import sympy as sp

import semialg
from semialg.algebraic_decomposition import associated_primes, primary_decomposition


@dataclass(frozen=True)
class SingularComparison:
    dimension: int
    component_count: int
    reconstruction_ok: bool
    radicals_ok: bool
    primary_components_ok: bool
    multiplicities_ok: bool

    @property
    def agreed(self) -> bool:
        return all(
            (
                self.reconstruction_ok,
                self.radicals_ok,
                self.primary_components_ok,
                self.multiplicities_ok,
            )
        )


def _singular_expr(expr: sp.Expr) -> str:
    return sp.sstr(sp.expand(expr)).replace("**", "^")


def _ideal_decl(name: str, generators: tuple[sp.Expr, ...]) -> str:
    body = ",".join(_singular_expr(g) for g in generators) or "0"
    return f"ideal {name}={body};"


def _singular_compare(generators, variables, result) -> SingularComparison | None:
    executable = shutil.which("Singular")
    if executable is None:
        return None
    ring_vars = ",".join(str(v) for v in variables)
    lines = [
        f"ring r=0,({ring_vars}),dp;",
        'LIB "primdec.lib";',
        "proc idealEq(ideal A, ideal B) {",
        "  ideal RA=reduce(std(A),std(B));",
        "  ideal RB=reduce(std(B),std(A));",
        "  if (size(RA)==0 && size(RB)==0) { return(1); }",
        "  return(0);",
        "}",
        _ideal_decl("I", tuple(generators)),
        "list P=primdecGTZ(I);",
    ]
    for index, component in enumerate(result.components, start=1):
        lines.append(_ideal_decl(f"SQ{index}", tuple(component.equations)))
        lines.append(_ideal_decl(f"SP{index}", tuple(component.radical)))
    if result.components:
        lines.append("ideal SJ=SQ1;")
        for index in range(2, len(result.components) + 1):
            lines.append(f"SJ=intersect(SJ,SQ{index});")
    else:
        lines.append("ideal SJ=1;")
    lines.extend(
        [
            'print("SEMIALG_DIM="+string(dim(std(I))));',
            'print("SEMIALG_COMPONENTS="+string(size(P)));',
            'print("SEMIALG_RECON="+string(idealEq(I,SJ)));',
            "int radicals_ok=1; int primary_ok=1; int multiplicities_ok=1;",
        ]
    )
    for index, component in enumerate(result.components, start=1):
        lines.extend(
            [
                "int matched=0;",
                "for (int k=1; k<=size(P); k++) {",
                f"  if (idealEq(SP{index},P[k][2])==1) {{",
                "    matched=1;",
                f"    if (idealEq(SQ{index},P[k][1])==0) {{ primary_ok=0; }}",
                f"    if (deg(std(P[k][1]))!={int(component.degree)}) {{ multiplicities_ok=0; }}",
                "  }",
                "}",
                "if (matched==0) { radicals_ok=0; primary_ok=0; multiplicities_ok=0; }",
            ]
        )
    lines.extend(
        [
            'print("SEMIALG_RADICALS="+string(radicals_ok));',
            'print("SEMIALG_PRIMARY="+string(primary_ok));',
            'print("SEMIALG_MULTIPLICITIES="+string(multiplicities_ok));',
            "quit;",
        ]
    )
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / "case.sing"
        path.write_text("\n".join(lines) + "\n", encoding="utf-8")
        proc = subprocess.run(
            [executable, "-q", str(path)],
            check=True,
            capture_output=True,
            text=True,
            timeout=60,
        )

    def value(name: str) -> int:
        match = re.search(rf"{name}=(\d+)", proc.stdout)
        if match is None:
            raise RuntimeError(f"could not parse {name} from Singular output:\n{proc.stdout}")
        return int(match.group(1))

    return SingularComparison(
        dimension=value("SEMIALG_DIM"),
        component_count=value("SEMIALG_COMPONENTS"),
        reconstruction_ok=bool(value("SEMIALG_RECON")),
        radicals_ok=bool(value("SEMIALG_RADICALS")),
        primary_components_ok=bool(value("SEMIALG_PRIMARY")),
        multiplicities_ok=bool(value("SEMIALG_MULTIPLICITIES")),
    )


def _generated_cases(count: int, seed: int):
    rng = random.Random(seed)
    x, y = sp.symbols("x y")
    cases = [
        ((x * y,), (x, y)),
        (((x + y) ** 2, y * (x + y)), (x, y)),
        ((x**2, x * y), (x, y)),
        ((x**2 * (x - 1) ** 3,), (x,)),
    ]
    while len(cases) < count:
        a = rng.choice((-3, -2, -1, 1, 2, 3))
        b = rng.choice((-3, -2, -1, 1, 2, 3))
        e1 = rng.randint(1, 3)
        e2 = rng.randint(1, 3)
        if rng.random() < 0.55:
            cases.append((((x - a) ** e1 * (x - b) ** e2,), (x,)))
        else:
            linear = x + a * y
            cases.append(((linear**e1, y * linear), (x, y)))
    return tuple(cases[:count])


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--generated", type=int, default=12)
    parser.add_argument("--seed", type=int, default=1729)
    args = parser.parse_args(argv)
    if shutil.which("Singular") is None:
        print("Singular executable not found; optional differential corpus skipped.")
        return 0

    cases = _generated_cases(max(3, args.generated), args.seed)
    for case_index, (generators, variables) in enumerate(cases, start=1):
        result = primary_decomposition(generators, variables)
        if not result.complete or not semialg.replay_certificate(result).verified:
            raise AssertionError(f"semialg failed exact replay for case {case_index}: {generators}")
        primes = associated_primes(generators, variables)
        if not primes.complete:
            raise AssertionError(
                f"associated primes incomplete for case {case_index}: {generators}"
            )
        singular = _singular_compare(tuple(generators), tuple(variables), result)
        assert singular is not None
        semialg_dim = max(component.dimension for component in result.components)
        if singular.dimension != semialg_dim:
            raise AssertionError(
                f"dimension mismatch for {generators}: Singular={singular.dimension}, semialg={semialg_dim}"
            )
        if singular.component_count != len(primes.primes):
            raise AssertionError(
                f"associated-prime count mismatch for {generators}: "
                f"Singular={singular.component_count}, semialg={len(primes.primes)}"
            )
        if not singular.agreed:
            raise AssertionError(f"primary differential mismatch for {generators}: {singular}")
        print(f"ok {case_index}/{len(cases)}: {generators}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
