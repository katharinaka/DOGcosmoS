# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'DOGcosmoS'
copyright = '2023, Katharina Kain'
author = 'Katharina Kain'
release = '1.0.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
	'sphinx.ext.autodoc',
	'sphinx.ext.napoleon'	
]

templates_path = ['_templates']
exclude_patterns = []

# -- Path ------------------------
import os
import sys
sys.path.insert(0, os.path.abspath('../../DOGcosmoS/'))
sys.path.insert(0, os.path.abspath('/home/fs71897/kkain/swiftgalaxy'))
sys.path.insert(0, os.path.abspath('/home/fs71897/kkain/.local/lib/python3.10/site-packages/'))

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']

# set master doc file --------------------
master_doc = 'index'
