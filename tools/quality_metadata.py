"""Environment and source identity shared by reproducible measurement tools."""

from __future__ import annotations

import hashlib
import importlib.metadata
import os
import platform
from pathlib import Path


def measurement_identity():
    """Record versions and a path-sensitive digest of the measured source tree."""
    import semialg

    measured_root = Path(semialg.__file__).resolve().parent
    modules = {
        "semialg/" + p.relative_to(measured_root).as_posix(): hashlib.sha256(
            p.read_bytes()
        ).hexdigest()
        for p in sorted(measured_root.rglob("*.py"))
    }
    digest = hashlib.sha256(
        "\n".join(f"{name}:{value}" for name, value in modules.items()).encode()
    ).hexdigest()
    return {
        "python": platform.python_version(),
        "platform": platform.platform(),
        "ground": os.environ.get("SYMPY_GROUND_TYPES", "auto"),
        "versions": {
            name: importlib.metadata.version(name)
            for name in ("sympy", "python-flint", "algroots", "semialg")
        },
        "module_origin": str(Path(semialg.__file__).resolve()),
        "source_sha256": digest,
        "module_sha256": modules,
    }
