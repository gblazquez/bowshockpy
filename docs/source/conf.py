# Configuration file for the Sphinx documentation builder.

# -- Project information


# sys.path.insert(0, os.path.abspath(".."))
#
# import bowshockpy
# from bowshockpy import __version__



project = 'bowshockpy'
copyright = '2025, Guillermo Blazquez-Calero & contributors'
author = 'Guillermo Blazquez-Calero & contributors'

release =  "0.3.1"
version =  "0.3.1"

# release = __version__
# version = __version__

# -- General configuration

extensions = [
    'sphinx.ext.duration',
    'sphinx.ext.doctest',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.intersphinx',
    'sphinx_design',
    'nbsphinx'
]

intersphinx_mapping = {
    'python': ('https://docs.python.org/3/', None),
    'sphinx': ('https://www.sphinx-doc.org/en/master/', None),
}
intersphinx_disabled_domains = ['std']

# autodoc_mock_imports = ["bowshockpy"]

templates_path = ['_templates']

# -- Options for HTML output

# html_theme = 'pydata_sphinx_theme'
# html_theme = 'sphinx_rtd_theme'
# html_theme = "furo"
html_theme = 'sphinx_book_theme'
# html_theme = "press"

html_title = "BowshockPy"
# html_logo = "path/to/myimage.png"

html_theme_options = {
    # "github_url": "https://github.com/gblazquez/bowshockpy",
    "repository_url": "https://github.com/gblazquez/bowshockpy",
    "use_repository_button": True,
#    "use_source_button": True,
}

# -- Options for EPUB output
epub_show_urls = 'footnote'
