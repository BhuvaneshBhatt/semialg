"""Run the permanent installed-wheel contract outside the source checkout."""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("wheel", type=Path)
    parser.add_argument("--source-root", type=Path, default=Path(__file__).resolve().parents[1])
    args = parser.parse_args(argv)

    wheel = args.wheel.resolve()
    source_root = args.source_root.resolve()
    contract = source_root / "tests" / "release" / "test_installed_wheel_contract.py"
    if not wheel.is_file():
        raise FileNotFoundError(wheel)

    with tempfile.TemporaryDirectory(prefix="semialg-wheel-contract-") as temp_name:
        temp = Path(temp_name)
        target = temp / "site-packages"
        target.mkdir()
        subprocess.run(
            [
                sys.executable,
                "-m",
                "pip",
                "install",
                "--no-deps",
                "--target",
                str(target),
                str(wheel),
            ],
            check=True,
        )
        copied = temp / contract.name
        shutil.copy2(contract, copied)
        env = os.environ.copy()
        env.update(
            {
                "PYTHONPATH": str(target),
                "SEMIALG_INSTALLED_WHEEL_TEST": "1",
                "SEMIALG_SOURCE_ROOT": str(source_root),
                "SEMIALG_TARGET_INSTALL": str(target),
            }
        )
        completed = subprocess.run(
            [sys.executable, "-m", "pytest", "-q", str(copied)],
            cwd=temp,
            env=env,
            check=False,
        )
        return completed.returncode


if __name__ == "__main__":
    raise SystemExit(main())
