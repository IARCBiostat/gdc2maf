"""Package-level invariants."""

import importlib
import sys
from pathlib import Path

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
