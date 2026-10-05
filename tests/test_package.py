"""Package-level invariants."""

import importlib
import sys
from pathlib import Path

import pytest

import gdc2maf


def test_no_module_shadows_a_standard_library_module():
    # a module named e.g. select.py breaks `import select` for anything that
    # ends up with the package directory on sys.path
    package_dir = Path(gdc2maf.__file__).parent
    names = {p.stem for p in package_dir.glob("*.py")} - {"__init__"}
    assert not (names & sys.stdlib_module_names)


def test_everything_in_all_is_importable():
    for name in gdc2maf.__all__:
        assert hasattr(gdc2maf, name), name


def test_version_is_exposed():
    assert gdc2maf.__version__.count(".") == 2


def test_importing_the_package_configures_no_logging_handlers():
    import logging

    handlers = logging.getLogger("gdc2maf").handlers
    assert all(isinstance(h, logging.NullHandler) for h in handlers)


def test_every_module_imports_cleanly():
    for module in ["api", "attrition", "cases", "cli", "client", "clinical", "cohort",
                   "download", "files", "maf", "provenance", "reports", "selection",
                   "spec", "tcga"]:
        assert importlib.import_module(f"gdc2maf.{module}")


def read_pyproject_dependencies():
    """Return the distribution names in [project.dependencies], and the repo root."""
    import re

    # tomllib is 3.11+; the package supports 3.10, so skip rather than fail there
    tomllib = pytest.importorskip("tomllib")

    root = Path(gdc2maf.__file__).parent.parent.parent
    with open(root / "pyproject.toml", "rb") as f:
        deps = tomllib.load(f)["project"]["dependencies"]
    return {re.split(r"[<>=!~\s]", d)[0] for d in deps}, root


def requirement_lines(root):
    """Return the non-comment, non-blank lines of requirements.txt."""
    text = (root / "requirements.txt").read_text()
    return [
        line.strip()
        for line in text.splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ]


def test_requirements_txt_is_exactly_the_runtime_dependencies():
    # It must not grow a second copy of the dev or docs extras: those are
    # declared once in pyproject.toml, and a copy here would drift.
    deps, root = read_pyproject_dependencies()
    import re

    listed = {re.split(r"[<>=!~\s]", line)[0] for line in requirement_lines(root)}
    assert listed == deps, f"requirements.txt lists {listed}, pyproject lists {deps}"


def test_conda_environment_covers_every_runtime_dependency():
    deps, root = read_pyproject_dependencies()
    listed = (root / "environment.yml").read_text()
    missing = [d for d in deps if d not in listed]
    assert not missing, f"missing from environment.yml: {missing}"



README_MARKERS = ["intro-start", "intro-end", "quickstart-start", "quickstart-end"]


def test_readme_keeps_the_markers_the_docs_home_page_includes():
    # docs/index.md pulls these slices out of the README with MyST's include
    # directive, so the home page is never a second copy of the text. Renaming
    # or dropping a marker would silently empty the home page.
    root = Path(gdc2maf.__file__).parent.parent.parent
    readme = (root / "README.md").read_text()
    for marker in README_MARKERS:
        assert f"<!-- docs-include: {marker} -->" in readme, marker


def test_docs_home_page_includes_the_readme_rather_than_repeating_it():
    root = Path(gdc2maf.__file__).parent.parent.parent
    index = (root / "docs" / "index.md").read_text()
    assert "{include} ../README.md" in index
    for marker in README_MARKERS:
        assert marker in index, marker


def test_every_sphinx_static_path_exists():
    # Sphinx warns when an html_static_path entry is missing, and the docs CI
    # builds with -W. An empty directory cannot be tracked by git, so a path
    # that exists only on a developer's machine fails the build on a fresh
    # checkout; each one needs a tracked file to keep it in the repository.
    import ast

    root = Path(gdc2maf.__file__).parent.parent.parent
    docs = root / "docs"
    if not docs.is_dir():  # installed without the repository alongside
        return

    conf = ast.parse((docs / "conf.py").read_text())
    static_paths = [
        ast.literal_eval(node.value)
        for node in conf.body
        if isinstance(node, ast.Assign)
        and any(
            getattr(t, "id", None) == "html_static_path" for t in node.targets
        )
    ]
    for paths in static_paths:
        for entry in paths:
            path = docs / entry
            assert path.is_dir(), (
                f"docs/conf.py html_static_path entry missing: {entry}"
            )
            assert any(path.iterdir()), (
                f"docs/{entry} is empty, so git will not track it and the "
                "docs build will warn on a fresh checkout"
            )
