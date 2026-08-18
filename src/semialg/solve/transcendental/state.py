from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field, replace

import sympy as sp

from ..domains import SolveDomain, normalize_domain


@dataclass(frozen=True)
class QuantifierBlock:
    quantifier: str
    variables: tuple[sp.Symbol, ...]


@dataclass(frozen=True)
class TransProblemState:
    formula: sp.Expr
    free_variables: tuple[sp.Symbol, ...]
    quantified_variables: tuple[sp.Symbol, ...] = ()
    parameter_variables: tuple[sp.Symbol, ...] = ()
    quantifier_blocks: tuple[QuantifierBlock, ...] = ()
    variable_domains: dict[sp.Symbol, tuple[object, ...]] = field(default_factory=dict)
    default_domain: SolveDomain = SolveDomain.REALS
    variable_order: tuple[sp.Symbol, ...] = ()
    notes: tuple[str, ...] = ()
    metadata: dict = field(default_factory=dict)

    @property
    def all_variables(self) -> tuple[sp.Symbol, ...]:
        seen = []
        for seq in (self.parameter_variables, self.free_variables, self.quantified_variables):
            for v in seq:
                if v not in seen:
                    seen.append(v)
        return tuple(seen)

    @property
    def active_variable_order(self) -> tuple[sp.Symbol, ...]:
        if self.variable_order:
            return self.variable_order
        return self.parameter_variables + self.free_variables + self.quantified_variables

    @property
    def has_quantifiers(self) -> bool:
        return bool(self.quantified_variables or self.quantifier_blocks)

    def with_formula(
        self, formula: sp.Expr, *, note: str | None = None, metadata_updates: Mapping | None = None
    ) -> TransProblemState:
        notes = self.notes if note is None else self.notes + (note,)
        metadata = dict(self.metadata)
        if metadata_updates:
            metadata.update(metadata_updates)
        return replace(self, formula=sp.simplify(formula), notes=notes, metadata=metadata)

    def add_note(self, note: str) -> TransProblemState:
        return replace(self, notes=self.notes + (note,))

    def with_variable_domains(
        self, extra_domains: Mapping[sp.Symbol, Sequence[object]]
    ) -> TransProblemState:
        merged = dict(self.variable_domains)
        for var, doms in extra_domains.items():
            merged[var] = tuple(doms)
        return replace(self, variable_domains=merged)

    def with_quantifier_blocks(self, blocks: Sequence[QuantifierBlock]) -> TransProblemState:
        qvars = tuple(v for block in blocks for v in block.variables)
        return replace(self, quantifier_blocks=tuple(blocks), quantified_variables=qvars)


def norm_quant_blocks(
    quantifier_blocks: Sequence[QuantifierBlock | tuple[str, Sequence[sp.Symbol]]] | None,
    quantified_variables: Sequence[sp.Symbol],
) -> tuple[QuantifierBlock, ...]:
    if quantifier_blocks is None:
        if quantified_variables:
            return (QuantifierBlock("exists", tuple(quantified_variables)),)
        return ()
    out = []
    seen: set[sp.Symbol] = set()
    for block in quantifier_blocks:
        if isinstance(block, QuantifierBlock):
            q, vs = block.quantifier, tuple(block.variables)
        else:
            q, raw = block
            vs = tuple(raw)
        q = str(q).lower()
        if q not in {"exists", "forall"}:
            raise ValueError("quantifier must be 'exists' or 'forall'")
        for v in vs:
            if v in seen:
                raise ValueError(f"duplicate quantified variable: {v}")
            seen.add(v)
        if vs:
            out.append(QuantifierBlock(q, vs))
    return tuple(out)


def build_trans_state(
    formula: sp.Expr,
    free_variables: Sequence[sp.Symbol],
    *,
    quantified_variables: Sequence[sp.Symbol] = (),
    parameter_variables: Sequence[sp.Symbol] = (),
    quantifier_blocks: Sequence[QuantifierBlock | tuple[str, Sequence[sp.Symbol]]] | None = None,
    variable_domains: Mapping[sp.Symbol, Sequence[object]] | None = None,
    default_domain: str | SolveDomain = SolveDomain.REALS,
    variable_order: Sequence[sp.Symbol] | None = None,
    notes: Sequence[str] = (),
    metadata: Mapping | None = None,
) -> TransProblemState:
    free_variables = tuple(free_variables)
    quantified_variables = tuple(quantified_variables)
    parameter_variables = tuple(parameter_variables)
    blocks = norm_quant_blocks(quantifier_blocks, quantified_variables)
    qvars = (
        tuple(v for block in blocks for v in block.variables) if blocks else quantified_variables
    )
    order = (
        tuple(variable_order)
        if variable_order is not None
        else parameter_variables + free_variables + qvars
    )
    vdomains = {k: tuple(v) for k, v in (variable_domains or {}).items()}
    return TransProblemState(
        formula=sp.simplify(formula),
        free_variables=free_variables,
        quantified_variables=qvars,
        parameter_variables=parameter_variables,
        quantifier_blocks=blocks,
        variable_domains=vdomains,
        default_domain=normalize_domain(default_domain),
        variable_order=order,
        notes=tuple(notes),
        metadata=dict(metadata or {}),
    )


__all__ = [
    "QuantifierBlock",
    "TransProblemState",
    "norm_quant_blocks",
    "build_trans_state",
]
