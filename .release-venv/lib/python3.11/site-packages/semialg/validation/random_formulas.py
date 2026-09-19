from __future__ import annotations

import random
from collections.abc import Sequence
from dataclasses import dataclass

import sympy as sp

from .corpus import ValidationCase


@dataclass(frozen=True)
class RandomFormulaConfig:
    seed: int = 0
    variables: tuple[str, ...] = ("x", "y")
    max_degree: int = 2
    max_terms: int = 4
    atom_count: int = 3
    quantifier_count: int = 1
    coefficient_min: int = -3
    coefficient_max: int = 3


def _monomial(rng: random.Random, symbols: Sequence[sp.Symbol], max_degree: int) -> sp.Expr:
    remaining = rng.randint(0, max_degree)
    powers = [0] * len(symbols)
    for _ in range(remaining):
        powers[rng.randrange(len(symbols))] += 1
    out = sp.Integer(1)
    for sym, power in zip(symbols, powers, strict=True):
        if power:
            out *= sym**power
    return out


def random_polynomial(
    rng: random.Random, symbols: Sequence[sp.Symbol], config: RandomFormulaConfig
) -> sp.Expr:
    terms: list[sp.Expr] = []
    for _ in range(rng.randint(1, config.max_terms)):
        coeff = 0
        while coeff == 0:
            coeff = rng.randint(config.coefficient_min, config.coefficient_max)
        terms.append(sp.Integer(coeff) * _monomial(rng, symbols, config.max_degree))
    poly = sp.expand(sum(terms, sp.Integer(0)))
    return poly if poly != 0 else sp.Integer(1)


def random_atom_text(
    rng: random.Random, symbols: Sequence[sp.Symbol], config: RandomFormulaConfig
) -> str:
    op = rng.choice(("=", "!=", ">", ">=", "<", "<="))
    poly = random_polynomial(rng, symbols, config)
    return f"({sp.sstr(poly)} {op} 0)"


def random_matrix_text(
    rng: random.Random, symbols: Sequence[sp.Symbol], config: RandomFormulaConfig
) -> str:
    atoms = [random_atom_text(rng, symbols, config) for _ in range(config.atom_count)]
    if len(atoms) == 1:
        return atoms[0]
    text = atoms[0]
    for atom in atoms[1:]:
        connective = rng.choice(("and", "or"))
        text = f"({text}) {connective} ({atom})"
    return text


def random_validation_cases(
    count: int, config: RandomFormulaConfig | None = None
) -> tuple[ValidationCase, ...]:
    cfg = config or RandomFormulaConfig()
    rng = random.Random(cfg.seed)
    symbols = tuple(sp.Symbol(name, real=True) for name in cfg.variables)
    cases: list[ValidationCase] = []
    for idx in range(count):
        qvars = cfg.variables[-cfg.quantifier_count :] if cfg.quantifier_count else ()
        quantifiers = tuple((rng.choice(("exists", "forall")), name) for name in qvars)
        cases.append(
            ValidationCase(
                name=f"random_{cfg.seed}_{idx}",
                formula_text=random_matrix_text(rng, symbols, cfg),
                variables=cfg.variables,
                quantifiers=quantifiers,
                tags=("random",),
                metadata={"seed": cfg.seed, "index": idx},
            )
        )
    return tuple(cases)


__all__ = ["RandomFormulaConfig", "random_polynomial", "random_validation_cases"]
