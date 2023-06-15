# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

import os
import sys
import sphinx_rtd_theme

sys.path.insert(0, os.path.abspath('..'))

project = 'PolyCleaver'
copyright = '2023, Eric Mates-Torres'
author = 'Eric Mates-Torres'
release = '0.0.1'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ['sphinx.ext.viewcode', 'sphinx.ext.napoleon', 'sphinx.ext.autodoc', 'sphinx.ext.todo', 'sphinx_design']

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

html_logo = "_static/icon.png"
html_favicon = "_static/favicon.svg"
html_theme = 'sphinx_rtd_theme'
html_theme_options = {'display_version': True}
html_static_path = ['_static']
html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]
html_css_files = [
    'custom.css',
]




# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output


