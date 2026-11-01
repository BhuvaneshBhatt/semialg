"""Prepare and install the exact backend sources used by the test matrix."""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
import tarfile
import tomllib
from pathlib import Path, PurePosixPath

ROOT = Path(__file__).resolve().parents[1]


def source_digest(directory):
    """Hash the backend metadata and Python sources, independent of checkout location."""
    paths = [directory / "pyproject.toml", *sorted((directory / "src").rglob("*.py"))]
    digest = hashlib.sha256()
    for path in sorted(paths, key=lambda path: path.relative_to(directory).as_posix()):
        digest.update(path.relative_to(directory).as_posix().encode() + b"\0")
        digest.update(path.read_bytes() + b"\0")
    return digest.hexdigest()


def prepare_backend(name, specification, destination, source_root=None):
    """Restore a pinned baseline, apply its checked patch, and verify the resulting sources."""
    patch = ROOT / "ci" / specification["patch"]
    if hashlib.sha256(patch.read_bytes()).hexdigest() != specification["patch_sha256"]:
        raise ValueError(f"{name}: backend patch checksum does not match ci/backends.toml")
    directory = destination / name
    directory.mkdir(parents=True, exist_ok=False)
    git = ["git", "-c", "core.autocrlf=false", "-C", str(directory)]
    subprocess.run([*git, "init", "--quiet"], check=True)
    if source_root:
        source = str((source_root / name).resolve())
        subprocess.run(
            [*git, "fetch", "--quiet", "--depth=1", source, specification["commit"]], check=True
        )
        subprocess.run([*git, "checkout", "--quiet", "--detach", "FETCH_HEAD"], check=True)
        commit = subprocess.check_output([*git, "rev-parse", "HEAD"], text=True).strip()
    else:
        baseline = ROOT / "ci" / specification["source_archive"]
        if (
            hashlib.sha256(baseline.read_bytes()).hexdigest()
            != specification["source_archive_sha256"]
        ):
            raise ValueError(f"{name}: baseline archive checksum does not match ci/backends.toml")
        with tarfile.open(baseline) as archive:
            commit = archive.pax_headers.get("comment")
            for member in archive:
                if member.isdir():
                    continue
                path = PurePosixPath(member.name)
                if (
                    not member.isfile()
                    or path.is_absolute()
                    or ".." in path.parts
                    or "\\" in member.name
                ):
                    raise ValueError(f"{name}: invalid baseline archive member {member.name}")
                target = directory.joinpath(*path.parts)
                target.parent.mkdir(parents=True, exist_ok=True)
                target.write_bytes(archive.extractfile(member).read())
    if commit != specification["commit"]:
        raise ValueError(f"{name}: baseline commit does not match ci/backends.toml")
    subprocess.run([*git, "apply", "--check", str(patch)], check=True)
    subprocess.run([*git, "apply", str(patch)], check=True)
    digest = source_digest(directory)
    if digest != specification["source_sha256"]:
        raise ValueError(f"{name}: patched source checksum does not match ci/backends.toml")
    return directory, {"commit": commit, "source_sha256": digest}


def main():
    """Install core or real backends, with an offline source option for validation."""
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument("--destination", type=Path, default=ROOT / ".ci-backends")
    parser.add_argument(
        "--source-root",
        type=Path,
        help="use local Git repositories instead of bundled baseline archives",
    )
    parser.add_argument(
        "--prepare-only", action="store_true", help="verify sources without installing them"
    )
    args = parser.parse_args()
    specifications = tomllib.loads((ROOT / "ci/backends.toml").read_text(encoding="utf-8"))
    packages = ("algroots",)
    directories, manifest = [], {}
    for name in packages:
        directory, record = prepare_backend(
            name, specifications[name], args.destination.resolve(), args.source_root
        )
        directories.append(str(directory))
        manifest[name] = record
    (args.destination / "manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8", newline="\n"
    )
    if not args.prepare_only:
        subprocess.run([sys.executable, "-m", "pip", "install", *directories], check=True)
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
