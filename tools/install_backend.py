"""Install declared backend dependencies or their matching local source trees."""

import argparse
import subprocess
import sys
import tomllib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def main():
    """Install the declared core or real dependency profile."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profile", choices=("core", "real"), default="real")
    parser.add_argument("--source-root", type=Path, default=ROOT.parent)
    args = parser.parse_args()
    specifications = tomllib.loads((ROOT / "ci/backends.toml").read_text())
    names = [name for name in specifications if name != "semialg" or args.profile == "real"]
    targets = []
    for name in names:
        source = args.source_root / name
        targets.append(
            str(source.resolve())
            if (source / "pyproject.toml").is_file()
            else specifications[name]["requirement"]
        )
    subprocess.run([sys.executable, "-m", "pip", "install", *targets], check=True)


if __name__ == "__main__":
    main()
