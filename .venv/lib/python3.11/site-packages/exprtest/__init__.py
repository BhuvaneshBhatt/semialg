"""Fast symbolic expression predicates with exact, numeric, and heuristic methods."""

from importlib.metadata import PackageNotFoundError, version

from ._domains import (
    ElementOf,
    is_algebraic,
    is_integer,
    is_prime,
    is_rational,
    is_real,
    number_kind,
)
from ._result import Verdict, ZeroClassification
from .core import profile_zerotest, zerotest

try:
    __version__ = version("exprtest")
except PackageNotFoundError:
    __version__ = "0+unknown"

__all__ = [
    "ElementOf",
    "Verdict",
    "ZeroClassification",
    "is_algebraic",
    "is_integer",
    "is_prime",
    "is_rational",
    "is_real",
    "number_kind",
    "profile_zerotest",
    "zerotest",
]
