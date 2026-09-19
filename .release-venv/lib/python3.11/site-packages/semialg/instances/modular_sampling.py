from __future__ import annotations

import random
from collections.abc import Iterable
from itertools import product


def cartesian_power_tuples(values: Iterable[object], repeat: int) -> list[tuple[object, ...]]:
    pool = list(values)
    if repeat < 0:
        raise ValueError("repeat must be nonnegative")
    if repeat == 0:
        return [()]
    return list(product(pool, repeat=repeat))


def sample_modular_points(
    variable_count: int, modulus: int, sample_count: int, *, seed: int | None = None
) -> list[tuple[int, ...]]:
    if variable_count < 0:
        raise ValueError("variable_count must be nonnegative")
    if modulus <= 0:
        raise ValueError("modulus must be positive")
    if sample_count <= 0:
        return []
    universe_size = modulus**variable_count
    if universe_size <= max(1000, 10 * sample_count):
        pts = cartesian_power_tuples(range(modulus), variable_count)
        rng = random.Random(seed)
        rng.shuffle(pts)
        return pts[:sample_count]
    rng = random.Random(seed)
    seen: set[tuple[int, ...]] = set()
    out: list[tuple[int, ...]] = []
    while len(out) < sample_count:
        pt = tuple(rng.randrange(modulus) for _ in range(variable_count))
        if pt not in seen:
            seen.add(pt)
            out.append(pt)
    return out


__all__ = ["cartesian_power_tuples", "sample_modular_points"]
