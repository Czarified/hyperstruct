"""Sphinx configuration."""

project = "Hyperstruct"
author = "Benjamin Crews"
copyright = "2024, Benjamin Crews"
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx_click",
    "myst_parser",
    "sphinxcontrib.bibtex",
]
bibtex_bibfiles = ["references.bib"]
autodoc_typehints = "description"
html_theme = "furo"
