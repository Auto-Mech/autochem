# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'AutoChem'
copyright = '2024, Andreas V. Copan'
author = 'Andreas V. Copan'
release = '2024.6.5'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'myst_parser',
    'sphinx.ext.autodoc',
    'sphinx_autodoc_typehints',
    'sphinx.ext.mathjax',
    'sphinx.ext.intersphinx',
]

autodoc_default_options = {'members': True}
autoclass_content = 'class'

intersphinx_mapping = {
  'python': ('https://docs.python.org/3.12', None),
  'numpy': ('https://docs.scipy.org/doc/numpy/', None),
}

templates_path = ['_templates']
exclude_patterns = []



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_book_theme'
html_static_path = ['_static']
html_theme_options = {
  "show_toc_level":2
}
