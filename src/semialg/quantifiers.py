"""Small symbolic quantifier nodes used by semialg.

SymPy does not currently expose public ``Exists``/``ForAll`` expression
classes.  semialg needs quantified formulas in a few reconstruction paths, so
this module provides deliberately small, SymPy-compatible Boolean nodes.

The nodes are representation objects rather than quantifier-elimination
engines.  They understand lexical binding well enough to report free symbols
and perform capture-avoiding ``subs`` operations.  Their bound variables range
over the ambient domain of the surrounding semialg problem; narrower domains
should be stated explicitly in the body, e.g. ``Contains(k, S.Integers)``.
"""

from __future__ import annotations

from collections.abc import Iterable

import sympy as sp
from sympy.core.basic import Basic
from sympy.logic.boolalg import Boolean, as_Boolean


def _normalize_variables(variables: sp.Symbol | Iterable[sp.Symbol]) -> tuple[sp.Symbol, ...]:
    if isinstance(variables, sp.Symbol):
        result = (variables,)
    else:
        try:
            result = tuple(variables)
        except TypeError as exc:
            raise TypeError("quantifier variables must be a Symbol or iterable of Symbols") from exc

    if not result:
        return ()
    if any(not isinstance(variable, sp.Symbol) for variable in result):
        raise TypeError("quantifier variables must all be SymPy Symbols")
    if len(set(result)) != len(result):
        raise ValueError("quantifier variables must be distinct")
    return result


def _fresh_bound_symbol(variable: sp.Symbol) -> sp.Dummy:
    """Return a scope-local replacement preserving the variable assumptions."""

    assumptions = dict(variable.assumptions0)
    assumptions.pop("commutative", None)
    return sp.Dummy(variable.name, **assumptions)


class _Quantifier(Boolean):
    """Base class for semialg's lexical quantifier expressions."""

    __slots__ = ()
    is_Quantifier = True

    def __new__(
        cls,
        variables: sp.Symbol | Iterable[sp.Symbol],
        formula: sp.Basic | bool,
    ):
        bound = _normalize_variables(variables)
        body = as_Boolean(formula)

        if not bound:
            return body
        if body in (sp.true, sp.false):
            return body

        # Vacuous binders carry no information and complicate substitution.
        body_free = body.free_symbols
        bound = tuple(variable for variable in bound if variable in body_free)
        if not bound:
            return body

        return Basic.__new__(cls, sp.Tuple(*bound), body)

    @property
    def variables(self) -> tuple[sp.Symbol, ...]:
        return tuple(self.args[0])

    @property
    def bound_symbols(self) -> tuple[sp.Symbol, ...]:
        """Bound variables, following the naming convention of SymPy binders."""

        return self.variables

    @property
    def formula(self) -> Boolean:
        return self.args[1]

    @property
    def free_symbols(self) -> set[sp.Symbol]:
        return set(self.formula.free_symbols) - set(self.variables)

    def _eval_subs(self, old: Basic, new: Basic):
        """Perform capture-avoiding substitution beneath the binder."""

        bound = set(self.variables)
        old_free = set(getattr(old, "free_symbols", set()))

        # An occurrence involving a bound variable is local to this scope and
        # must not be targeted by an outer substitution.
        if old in bound or old_free & bound:
            return self

        variables = self.variables
        body = self.formula
        new_free = set(getattr(new, "free_symbols", set()))
        capture = bound & new_free

        if capture:
            renaming = {
                variable: _fresh_bound_symbol(variable)
                for variable in variables
                if variable in capture
            }
            renamed_body = body.xreplace(renaming)
            renamed_variables = tuple(renaming.get(variable, variable) for variable in variables)
        else:
            renamed_body = body
            renamed_variables = variables

        replaced_body = renamed_body._subs(old, new)
        if replaced_body == renamed_body:
            return self
        return self.func(renamed_variables, replaced_body)

    def _eval_simplify(self, **kwargs):
        simplified = sp.simplify(self.formula, **kwargs)
        return self.func(self.variables, simplified)

    def _sympystr(self, printer) -> str:
        variables = self.variables
        if len(variables) == 1:
            rendered_variables = printer.doprint(variables[0])
        else:
            rendered_variables = "(" + ", ".join(printer.doprint(v) for v in variables) + ")"
        return f"{self.func.__name__}({rendered_variables}, {printer.doprint(self.formula)})"

    def _latex(self, printer) -> str:
        symbol = r"\exists" if isinstance(self, Exists) else r"\forall"
        variables = ", ".join(printer.doprint(variable) for variable in self.variables)
        return rf"{symbol}_{{{variables}}}\,\left({printer.doprint(self.formula)}\right)"


class Exists(_Quantifier):
    """Existentially quantify one or more variables in a Boolean formula.

    ``Exists(x, phi)`` and ``Exists((x, y), phi)`` are supported.  The node is
    symbolic: it does not attempt quantifier elimination by itself.
    """


class ForAll(_Quantifier):
    """Universally quantify one or more variables in a Boolean formula."""


def apply_quantifiers(
    formula: sp.Basic | bool,
    quantifiers: Iterable[tuple[str, sp.Symbol]],
) -> Boolean:
    """Wrap ``formula`` in a prenex quantifier prefix.

    The iterable is ordered from outermost to innermost, matching semialg's
    existing internal ``(name, symbol)`` quantifier representation.
    """

    body = as_Boolean(formula)
    items = tuple(quantifiers)
    for name, variable in reversed(items):
        lowered = str(name).lower()
        if lowered == "exists":
            body = Exists(variable, body)
        elif lowered == "forall":
            body = ForAll(variable, body)
        else:
            raise ValueError(f"unsupported quantifier: {name!r}")
    return body


def split_quantifiers(
    formula: sp.Basic | bool,
) -> tuple[tuple[tuple[str, sp.Symbol], ...], Boolean]:
    """Return a prenex quantifier prefix and its quantifier-free matrix.

    Adjacent multi-variable quantifier nodes are flattened in lexical order.
    Only a leading prenex prefix is removed; quantified subformulas inside the
    matrix are intentionally preserved.
    """

    current = as_Boolean(formula)
    prefix: list[tuple[str, sp.Symbol]] = []
    while isinstance(current, _Quantifier):
        name = "exists" if isinstance(current, Exists) else "forall"
        prefix.extend((name, variable) for variable in current.variables)
        current = current.formula
    return tuple(prefix), current


__all__ = ["Exists", "ForAll", "apply_quantifiers", "split_quantifiers"]
