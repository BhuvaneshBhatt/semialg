"""Small examples comparing semialg workflows with familiar CAD/QE tasks."""

from __future__ import annotations

import sympy as sp

from semialg import cad, find_instance, generic_cad


def main() -> None:
    x, y = sp.symbols("x y", real=True)
    disk = x**2 + y**2 <= 1
    print("CAD:", cad(disk, [x, y]))
    print("generic:", generic_cad(disk, [x, y]))
    print("witness:", find_instance(disk & (x > 0), [x, y]).first())


if __name__ == "__main__":
    main()
