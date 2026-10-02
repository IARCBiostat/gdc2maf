"""Sphinx configuration for the gdc2maf documentation.

Nothing here restates documentation that lives elsewhere. The home page
includes the README's own text, the command line page is generated from the
argparse parser, and the reference pages are generated from the docstrings, so
every piece of documentation has exactly one source.
"""

from importlib import metadata

project = "gdc2maf"
author = "Ali Farnudi"
copyright = "2026, Ali Farnudi"  # noqa: A001 - Sphinx expects this name
release = metadata.version("gdc2maf")
version = release

extensions = [
    "sphinx.ext.autodoc",       # pull in the docstrings
    "sphinx.ext.autosummary",   # tables of objects, with generated stub pages
    "sphinx.ext.napoleon",      # read numpydoc-style sections
    "sphinx.ext.intersphinx",   # link to pandas and the standard library
    "sphinx.ext.viewcode",      # a "source" link next to each object
    "myst_parser",              # Markdown, so the README can be included as-is
    "sphinxarg.ext",            # the CLI page, generated from the parser
    "sphinx_design",            # tab sets, so one page can show both interfaces
]

# Markdown everywhere; .rst stays available for anything that needs it.
source_suffix = {".md": "markdown", ".rst": "restructuredtext"}
myst_enable_extensions = ["colon_fence", "deflist", "substitution"]
myst_heading_anchors = 3

autosummary_generate = True
autodoc_member_order = "bysource"
autodoc_typehints = "description"
autodoc_default_options = {
    "members": True,
    "undoc-members": False,
    "show-inheritance": True,
}
napoleon_numpy_docstring = True
napoleon_google_docstring = False

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "pandas": ("https://pandas.pydata.org/docs", None),
}

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

html_theme = "furo"
html_static_path = ["_static"]
html_title = f"gdc2maf {release}"
html_theme_options = {
    "source_repository": "https://github.com/IARCBiostat/gdc2maf/",
    "source_branch": "main",
    "source_directory": "docs/",
}

# Warn about every cross-reference that does not resolve. The CI builds with
# -W, so a role pointing at a name that is not documented fails the build
# instead of quietly rendering as plain text.
nitpicky = True

# Types that appear in numpydoc "Parameters" sections and are not Python
# objects Sphinx can link to: prose descriptions, pandas spellings used in the
# docstrings, and the loose spellings numpydoc encourages.
nitpick_ignore_regex = [
    (r"py:class", r"pd\..*"),
    (r"py:class", r"np\..*"),
    (r"py:class", r"optional"),
    (r"py:class", r"array_like"),
    (r"py:class", r"file object"),
    (r"py:class", r"sequence"),
    (r"py:class", r"iterable"),
    (r"py:class", r"pairs"),
    (r"py:class", r".*\s.*"),  # anything with a space is prose, not a type
]
