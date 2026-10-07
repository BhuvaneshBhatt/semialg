from __future__ import annotations

import json
import re
from collections.abc import Iterable
from dataclasses import dataclass
from functools import lru_cache
from importlib.resources import files
from typing import Any

import sympy as sp
from sympy import S, Symbol

from .conditions import parse_catalogue_expression
from .exactness import certify_identity_equal


class FunctionSlot(sp.Dummy):
    """A rule-local placeholder for one positional function argument.

    ``FunctionSlot`` deliberately behaves like a SymPy expression so registry
    conditions and results can be written naturally, while retaining an
    explicit positional binding.  Every construction is a fresh ``Dummy``; a
    human-readable name never creates a global-symbol dependency.
    """

    def __new__(cls, position: int, name: str | None = None):
        if isinstance(position, bool) or not isinstance(position, int) or position < 0:
            raise ValueError("FunctionSlot position must be a non-negative integer")
        display = name or f"slot{position}"
        obj = sp.Dummy.__new__(cls, display)
        obj._function_slot_position = position
        obj._function_slot_name = name
        return obj

    @property
    def position(self) -> int:
        return self._function_slot_position

    @property
    def display_name(self) -> str | None:
        return self._function_slot_name


def make_function_slots(*names: str) -> tuple[FunctionSlot, ...]:
    """Create fresh, positionally bound registry placeholders.

    Examples
    --------
    ``nu, z = make_function_slots("nu", "z")`` binds ``nu`` to argument slot 0 and
    ``z`` to argument slot 1 of whichever function head owns the rule.
    """

    return tuple(FunctionSlot(i, name) for i, name in enumerate(names))


def _replace_slots(value, substitutions: dict[FunctionSlot, Any]):
    if isinstance(value, sp.Basic):
        return value.xreplace(substitutions)
    if isinstance(value, tuple):
        return tuple(_replace_slots(v, substitutions) for v in value)
    if isinstance(value, list):
        return [_replace_slots(v, substitutions) for v in value]
    if isinstance(value, dict):
        return {
            _replace_slots(k, substitutions): _replace_slots(v, substitutions)
            for k, v in value.items()
        }
    return value


@dataclass(frozen=True)
class InstantiatedPropertyRule:
    """A native registry theorem instantiated for one concrete function call."""

    rule: "PropertyRule"
    expression: sp.Basic
    subject: tuple[Any, ...]
    conditions: Any
    result: Any
    substitutions: tuple[tuple[FunctionSlot, Any], ...]

    @property
    def substitution_map(self) -> dict[FunctionSlot, Any]:
        return dict(self.substitutions)


@dataclass(frozen=True)
class PropertyRule:
    """A localized, typed theorem about positional arguments of a function head.

    This type is responsible only for theorem matching and instantiation.
    Transport through compositions such as ``F(g(x))`` is handled separately
    by the registry transport layer.
    """

    head: Any
    property: str
    subject: tuple[FunctionSlot, ...]
    conditions: Any = S.true
    result: Any = None
    proof_kind: str | None = None
    verifier: Any = None
    source: str | None = None

    def __post_init__(self):
        subjects = (
            (self.subject,)
            if isinstance(self.subject, FunctionSlot)
            else tuple(self.subject)
        )
        if not subjects:
            raise ValueError(
                "PropertyRule subject must contain at least one FunctionSlot"
            )
        if not all(isinstance(slot, FunctionSlot) for slot in subjects):
            raise TypeError(
                "PropertyRule subject entries must be FunctionSlot instances"
            )
        if len(set(subjects)) != len(subjects):
            raise ValueError("PropertyRule subject cannot contain duplicate slots")
        object.__setattr__(self, "subject", subjects)
        object.__setattr__(self, "conditions", sp.sympify(self.conditions))

        used = set(subjects)
        for value in (self.conditions, self.result):
            if isinstance(value, sp.Basic):
                used.update(value.atoms(FunctionSlot))
        if used:
            positions = [slot.position for slot in used]
            if len(set(positions)) != len(positions):
                raise ValueError(
                    "A PropertyRule cannot contain distinct FunctionSlots for the same position"
                )

    @property
    def slots(self) -> tuple[FunctionSlot, ...]:
        used = set(self.subject)
        for value in (self.conditions, self.result):
            if isinstance(value, sp.Basic):
                used.update(value.atoms(FunctionSlot))
        return tuple(sorted(used, key=lambda slot: slot.position))

    def matches(self, expr) -> bool:
        expr = sp.sympify(expr)
        if getattr(expr, "func", None) != self.head:
            return False
        nargs = len(getattr(expr, "args", ()))
        return all(slot.position < nargs for slot in self.slots)

    def substitutions(self, expr) -> dict[FunctionSlot, Any]:
        expr = sp.sympify(expr)
        if not self.matches(expr):
            raise ValueError(f"Rule for {self.head} does not match {expr}")
        args = expr.args
        return {slot: args[slot.position] for slot in self.slots}

    def instantiate(self, expr) -> InstantiatedPropertyRule:
        expr = sp.sympify(expr)
        repl = self.substitutions(expr)
        return InstantiatedPropertyRule(
            rule=self,
            expression=expr,
            subject=tuple(repl[slot] for slot in self.subject),
            conditions=_replace_slots(self.conditions, repl),
            result=_replace_slots(self.result, repl),
            substitutions=tuple((slot, repl[slot]) for slot in self.slots),
        )


def _balanced_content(text: str, start: int, open_char="[", close_char="]"):
    depth = 0
    for i in range(start, len(text)):
        ch = text[i]
        if ch == open_char:
            depth += 1
        elif ch == close_char:
            depth -= 1
            if depth == 0:
                return text[start + 1 : i]
    return None


def _split_top_level(text: str):
    out = []
    start = 0
    depth = 0
    pairs = {"[": "]", "(": ")", "{": "}"}
    closes = set(pairs.values())
    stack = []
    for i, ch in enumerate(text):
        if ch in pairs:
            stack.append(pairs[ch])
            depth += 1
        elif ch in closes and stack and ch == stack[-1]:
            stack.pop()
            depth -= 1
        elif ch == "," and depth == 0:
            out.append(text[start:i].strip())
            start = i + 1
    out.append(text[start:].strip())
    return out


def _source_signature(lhs: str | None):
    """Extract the represented function pattern and property mode/variable."""
    if not lhs:
        return None, (), None, None
    lhs = lhs.split("/;", 1)[0].strip()
    first = lhs.find("[")
    if first < 0:
        return None, (), None, None
    outer = _balanced_content(lhs, first)
    if outer is None:
        return None, (), None, None
    outer_args = _split_top_level(outer)
    if not outer_args:
        return None, (), None, None
    fpat = outer_args[0].strip()
    if fpat.startswith("HoldPattern["):
        i = fpat.find("[")
        inner = _balanced_content(fpat, i)
        if inner is not None:
            fpat = inner.strip()
    hidx = fpat.find("[")
    if hidx < 0:
        return fpat, (), None, outer_args[1] if len(outer_args) > 1 else None
    head = fpat[:hidx].strip()
    inner = _balanced_content(fpat, hidx)
    args = tuple(_split_top_level(inner or "")) if inner is not None else ()
    prop_head = lhs[:first].strip()
    mode = None
    var_token = None
    if "SingularityProperties" in prop_head and len(outer_args) > 1:
        token = outer_args[1].strip()
        if token == "True":
            mode = "real"
        elif token == "False":
            mode = "complex"
        else:
            mode = "any"
    if "RealUnivariateProperties" in prop_head and len(outer_args) > 1:
        var_token = outer_args[1].strip()
    return head, args, mode, var_token


def _pattern_name(token: str):
    token = token.strip()
    m = re.match(r"([A-Za-z$][A-Za-z0-9$]*)_", token)
    if m:
        return m.group(1)
    m = re.match(r"([A-Za-z$][A-Za-z0-9$]*):", token)
    if m:
        return m.group(1)
    return None


def _literal_pattern_matches(token: str, actual) -> bool:
    token = token.strip()
    # Ordinary named patterns are handled by guards/assumptions.
    if "_" in token and re.match(r"^[A-Za-z$][A-Za-z0-9$]*_", token):
        return True
    # a:(2|3), a:3 style restricted alternatives/defaults.
    if ":" in token and re.match(r"^[A-Za-z$][A-Za-z0-9$]*:", token):
        rhs = token.split(":", 1)[1].strip()
        if rhs.startswith("(") and rhs.endswith(")") and "|" in rhs:
            vals = rhs[1:-1].split("|")
            return any(_literal_pattern_matches(v, actual) for v in vals)
        if rhs.startswith("_"):
            return True
        return _literal_pattern_matches(rhs, actual)
    # List/vector patterns and Blank/Repeated constructs are too rich to match
    # syntactically here; leave them to their explicit registry guard.
    if any(x in token for x in ("_", "Repeated[", "Vectors[", "Matrices[")):
        return True
    try:
        literal = parse_catalogue_expression(token)
        if isinstance(actual, sp.Basic) and isinstance(literal, sp.Basic):
            return certify_identity_equal(actual, literal) is True
        return actual == literal
    except (TypeError, ValueError):
        return False


@dataclass(frozen=True)
class RuleCase:
    value_source: str
    condition_source: str
    source_property_type: str | None = None
    source_lhs: str | None = None

    @property
    def condition(self):
        return parse_catalogue_expression(self.condition_source or "True")

    @property
    def value(self):
        return parse_catalogue_expression(self.value_source)

    @property
    def signature(self):
        return _source_signature(self.source_lhs)

    def substitutions(self, expr):
        _, patterns, _, _ = self.signature
        repl = {}
        for token, actual in zip(patterns, getattr(expr, "args", ()), strict=False):
            name = _pattern_name(token)
            if name:
                repl[Symbol(name)] = actual
        return repl

    def instantiate(self, expr):
        repl = self.substitutions(expr)
        value = self.value
        condition = self.condition
        if isinstance(value, sp.Basic):
            value = value.xreplace(repl)
        if isinstance(condition, sp.Basic):
            condition = condition.xreplace(repl)
        return value, condition

    def matches(self, expr, *, var=None, mode=None) -> bool:
        _, patterns, case_mode, var_token = self.signature
        args = getattr(expr, "args", ())
        if patterns and len(patterns) != len(args):
            return False
        for token, actual in zip(patterns, args, strict=False):
            if not _literal_pattern_matches(token, actual):
                return False
        if mode in ("real", "complex") and case_mode not in (None, "any", mode):
            return False
        if var is not None and var_token is not None:
            name = _pattern_name(var_token)
            if name:
                mapped = self.substitutions(expr).get(Symbol(name))
                if mapped is not None and mapped != var:
                    return False
        return True

    def specificity(self, *, mode=None) -> int:
        _, patterns, case_mode, _ = self.signature
        score = 0
        for token in patterns:
            t = token.strip()
            if ":" in t:
                score += 3
            elif re.search(r"_(Integer|Real|Complex|Positive|NonNegative)", t):
                score += 2
            elif "_" not in t:
                score += 4
        if mode in ("real", "complex") and case_mode == mode:
            score += 2
        return score


@dataclass(frozen=True)
class Family:
    raw: dict[str, Any]

    @property
    def id(self):
        return self.raw["family_id"]

    @property
    def sympy_names(self):
        return tuple(self.raw.get("sympy_names", ()))

    @property
    def mpmath_names(self):
        return tuple(self.raw.get("mpmath_names", ()))

    @property
    def catalogue_names(self):
        return tuple(self.raw.get("catalogue_names", ()))

    def cases(self, prop: str) -> tuple[RuleCase, ...]:
        p = self.raw.get("properties", {}).get(prop, {})
        return tuple(
            RuleCase(
                **{
                    k: c.get(k)
                    for k in (
                        "value_source",
                        "condition_source",
                        "source_property_type",
                        "source_lhs",
                    )
                }
            )
            for c in p.get("cases", ())
        )


class FunctionRegistry:
    def __init__(
        self,
        records: Iterable[dict[str, Any]],
        *,
        property_rules: Iterable[PropertyRule] = (),
        identity_rules: Iterable[Any] = (),
        derivative_rules: Iterable[Any] = (),
    ):
        self._families = tuple(Family(r) for r in records)
        self._by_name = {}
        for fam in self._families:
            for n in (
                fam.id,
                *fam.sympy_names,
                *fam.mpmath_names,
                *fam.catalogue_names,
            ):
                self._by_name.setdefault(n, fam)
        self._property_rules: list[PropertyRule] = []
        self._rules_by_head: dict[Any, list[PropertyRule]] = {}
        for rule in property_rules:
            self.register_property_rule(rule)
        self._identity_rules: list[Any] = []
        self._derivative_rules: list[Any] = []
        for rule in identity_rules:
            self.register_identity_rule(rule)
        for rule in derivative_rules:
            self.register_derivative_rule(rule)

    def __len__(self):
        return len(self._families)

    def __iter__(self):
        return iter(self._families)

    def get(self, name: str) -> Family | None:
        return self._by_name.get(name)

    def require(self, name: str) -> Family:
        f = self.get(name)
        if f is None:
            raise KeyError(name)
        return f

    def register_property_rule(self, rule: PropertyRule) -> PropertyRule:
        """Register one native localized property theorem.

        Catalogue-derived ``Family``/``RuleCase`` records remain independent in
        Catalogue entries can be adapted to this native rule API when their
        argument semantics are unambiguous.
        """

        if not isinstance(rule, PropertyRule):
            raise TypeError("rule must be a PropertyRule")
        self._property_rules.append(rule)
        self._rules_by_head.setdefault(rule.head, []).append(rule)
        return rule

    def register_identity_rule(self, rule):
        from .registry_derivatives import IdentityRule

        if not isinstance(rule, IdentityRule):
            raise TypeError("rule must be an IdentityRule")
        self._identity_rules.append(rule)
        return rule

    def register_derivative_rule(self, rule):
        from .registry_derivatives import DerivativeRule

        if not isinstance(rule, DerivativeRule):
            raise TypeError("rule must be a DerivativeRule")
        self._derivative_rules.append(rule)
        return rule

    def identity_rules(self, *, head=None):
        rules = tuple(self._identity_rules)
        return tuple(r for r in rules if head is None or r.head == head)

    def derivative_rules(self, *, head=None, position=None):
        rules = tuple(self._derivative_rules)
        return tuple(
            r
            for r in rules
            if (head is None or r.head == head)
            and (position is None or r.with_respect_to.position == position)
        )

    def match_identity_rules(self, expr):
        expr = sp.sympify(expr)
        return tuple(r for r in self.identity_rules(head=expr.func) if r.matches(expr))

    def match_derivative_rules(self, expr, *, position=None):
        expr = sp.sympify(expr)
        return tuple(
            r
            for r in self.derivative_rules(head=expr.func, position=position)
            if r.matches(expr)
        )

    def property_rules(
        self, *, head=None, property: str | None = None
    ) -> tuple[PropertyRule, ...]:
        """Return native rules, optionally filtered by head/property name."""

        rules = (
            tuple(self._rules_by_head.get(head, ()))
            if head is not None
            else tuple(self._property_rules)
        )
        if property is not None:
            rules = tuple(rule for rule in rules if rule.property == property)
        return rules

    def match_property_rules(
        self, expr, *, property: str | None = None
    ) -> tuple[PropertyRule, ...]:
        """Return native rules whose function head/slots match ``expr``."""

        expr = sp.sympify(expr)
        return tuple(
            rule
            for rule in self.property_rules(
                head=getattr(expr, "func", None), property=property
            )
            if rule.matches(expr)
        )

    def instantiate_property_rules(
        self, expr, *, property: str | None = None
    ) -> tuple[InstantiatedPropertyRule, ...]:
        """Instantiate all matching native rules for one concrete call."""

        return tuple(
            rule.instantiate(expr)
            for rule in self.match_property_rules(expr, property=property)
        )


@lru_cache(maxsize=1)
def registry() -> FunctionRegistry:
    p = files("funcprops").joinpath("data/catalogue.json")
    payload = json.loads(p.read_text())
    reg = FunctionRegistry(payload["families"])
    from .native_registry import install_native_rules

    install_native_rules(reg)
    return reg


def adapt_catalogue_rule_case(
    head,
    property_name: str,
    case: RuleCase,
    *,
    source: str | None = None,
) -> PropertyRule | None:
    """Conservatively adapt one catalogue case to a native ``PropertyRule``.

    Only ordinary positional named patterns are migrated. Literal/restricted,
    sequence, vector, and otherwise richer patterns deliberately remain in the
    catalogue representation because flattening them would lose theorem semantics.
    """

    _, patterns, _, _ = case.signature
    if len(patterns) != 1:
        return None
    names = []
    for token in patterns:
        token = token.strip()
        # Accept only simple named blanks, e.g. ``z_`` or ``nu_``.
        m = re.fullmatch(r"([A-Za-z$][A-Za-z0-9$]*)_", token)
        if m is None:
            return None
        names.append(m.group(1))
    slots = tuple(FunctionSlot(i, name) for i, name in enumerate(names))
    substitutions = {
        Symbol(name): slot for name, slot in zip(names, slots, strict=True)
    }
    value = case.value
    condition = case.condition
    if isinstance(value, sp.Basic):
        value = value.xreplace(substitutions)
    if isinstance(condition, sp.Basic):
        condition = condition.xreplace(substitutions)

    # Reject residual source-pattern symbols: the native rule would otherwise
    # appear localized while still depending on global catalogue names.
    residual = set()
    for obj in (value, condition):
        if isinstance(obj, sp.Basic):
            residual.update(obj.free_symbols)
    if residual & {Symbol(name) for name in names}:
        return None
    return PropertyRule(
        head,
        property_name,
        slots,
        conditions=condition,
        result=value,
        proof_kind="catalogue-adaptation",
        source=source or case.source_lhs or "catalogue",
    )
