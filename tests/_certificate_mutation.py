"""Table-driven mutation checks for replayable exact certificates.

Every populated dataclass field reached from a registered certificate is
mutated independently unless it is explicitly classified as search/provenance
metadata. Adding a field to a nested proof object therefore expands the test
surface automatically.
"""

from __future__ import annotations

from collections.abc import Callable, Iterator
from dataclasses import fields, is_dataclass, replace
from typing import Any

import sympy as sp

PROVENANCE_FIELDS = frozenset({"method", "primes", "modulus", "source", "notes"})


def _mutate_atom(value: Any) -> Any:
    if isinstance(value, bool):
        return not value
    if isinstance(value, int):
        return value + 1
    if isinstance(value, sp.Symbol):
        return sp.Dummy(f"tampered_{value}")
    if isinstance(value, sp.Poly):
        return sp.Poly(value.as_expr() + 1, *value.gens, domain=value.domain)
    if isinstance(value, sp.ImmutableMatrix):
        matrix = sp.MutableDenseMatrix(value)
        if matrix.rows and matrix.cols:
            matrix[0, 0] = sp.expand(matrix[0, 0] + 1)
        else:
            return sp.ImmutableMatrix([[1]])
        return sp.ImmutableMatrix(matrix)
    if isinstance(value, sp.Basic):
        return sp.expand(value + 1)
    if isinstance(value, str):
        if value == "lex":
            return "grevlex"
        if value == "grevlex":
            return "lex"
        return f"{value}_tampered"
    raise TypeError(f"no conservative certificate mutator for {type(value)!r}")


def _nested_mutations(value: Any, path: str) -> Iterator[tuple[str, Any]]:
    if is_dataclass(value) and not isinstance(value, type):
        for field in fields(value):
            if field.name in PROVENANCE_FIELDS:
                continue
            current = getattr(value, field.name)
            for child_path, changed in _nested_mutations(current, f"{path}.{field.name}"):
                yield child_path, replace(value, **{field.name: changed})
        return
    if isinstance(value, tuple):
        for index, item in enumerate(value):
            for child_path, changed in _nested_mutations(item, f"{path}[{index}]"):
                yield child_path, (*value[:index], changed, *value[index + 1 :])
        return
    if isinstance(value, frozenset):
        for item in value:
            for child_path, changed in _nested_mutations(item, f"{path}{{item}}"):
                yield child_path, frozenset((value - {item}) | {changed})
        return
    try:
        changed = _mutate_atom(value)
    except TypeError:
        return
    yield path, changed


def field_mutations(certificate: Any, *, ignored: frozenset[str] = frozenset()):
    """Yield one conservative mutation for every populated proof-field path."""
    if not is_dataclass(certificate):
        raise TypeError("certificate must be a dataclass instance")
    classified = ignored | PROVENANCE_FIELDS
    for field in fields(certificate):
        if field.name in classified:
            continue
        current = getattr(certificate, field.name)
        for path, changed in _nested_mutations(current, field.name):
            yield path, replace(certificate, **{field.name: changed})


def assert_certificate_fields_reject_mutation(
    certificate: Any,
    verifier: Callable[[Any], bool],
    *,
    ignored: frozenset[str] = frozenset(),
) -> None:
    """Require all generated proof mutations to be rejected by exact replay."""
    assert verifier(certificate), f"baseline {type(certificate).__name__} does not verify"
    tested = []
    for path, mutated in field_mutations(certificate, ignored=ignored):
        tested.append(path)
        try:
            accepted = bool(verifier(mutated))
        except (
            ArithmeticError,
            AttributeError,
            IndexError,
            TypeError,
            ValueError,
            sp.PolynomialError,
            sp.polys.polyerrors.CoercionFailed,
        ):
            accepted = False
        assert not accepted, (
            f"tampered proof field was accepted: {type(certificate).__name__}.{path}"
        )
    assert tested, f"no populated proof fields tested for {type(certificate).__name__}"
