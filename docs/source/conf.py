# Configuration file for the Sphinx documentation builder.

# -- Project information

project = 'bowshockpy'
copyright = '2025, Guillermo Blazquez-Calero'
author = 'Guillermo Blazquez-Calero'

release =  "0.1.1"
version =  "0.1.1"

# -- General configuration

extensions = [
    'sphinx.ext.duration',
    'sphinx.ext.doctest',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.intersphinx',
    'nbsphinx'
]

intersphinx_mapping = {
    'python': ('https://docs.python.org/3/', None),
    'sphinx': ('https://www.sphinx-doc.org/en/master/', None),
}
intersphinx_disabled_domains = ['std']

autodoc_mock_imports = ["bowshockpy"]

templates_path = ['_templates']

# -- Options for HTML output

html_theme = 'sphinx_rtd_theme'

# -- Options for EPUB output
epub_show_urls = 'footnote'
