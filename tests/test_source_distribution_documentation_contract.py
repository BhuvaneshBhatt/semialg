"""Source distributions must retain reproducible documentation-quality inputs."""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_documentation_quality_registries_and_assets_are_source_files():
    required = (
        "docs/reference/primary_api_manifest.toml",
        "docs/reference/root_function_documentation_adequacy.toml",
        "docs/assets/decision-flow.svg",
        "docs/assets/cad-flow.svg",
        "docs/assets/analysis-flow.svg",
        "scripts/generate_primary_api_manifest.py",
        "scripts/generate_root_api_documentation_adequacy.py",
    )
    assert [name for name in required if not (ROOT / name).is_file()] == []


def test_sdist_build_inputs():
    from setuptools import Distribution
    from setuptools.command.egg_info import manifest_maker

    manifest = manifest_maker(Distribution())
    manifest.initialize_options()
    manifest.finalize_options()
    manifest.template = str(ROOT / "MANIFEST.in")
    from setuptools.command.egg_info import FileList

    manifest.filelist = FileList()
    manifest.filelist.findall(str(ROOT))
    # Manifest commands match paths relative to the source root.
    manifest.filelist.allfiles = [
        str(Path(p).relative_to(ROOT)) for p in manifest.filelist.allfiles
    ]
    manifest.read_template()
    required = {
        "examples/catalog.json",
        "tools/check_examples.py",
        "tools/generate_specialist_api.py",
        "tools/install_backend.py",
        "ci/minimum.txt",
        "ci/backends.toml",
        "ci/sources/algroots.tar.gz",
        "ci/patches/algroots.patch",
        "docs/assets/linkage-workspace.svg",
        "tests/certificate_fields.json",
    }
    assert required <= set(manifest.filelist.files)
