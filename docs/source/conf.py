import os
import sys
sys.path.insert(0, os.path.abspath('../..'))

project = 'Optimization Projects'
copyright = '2024, Salar Nouri'
author = 'Salar Nouri'

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.viewcode',
    'sphinx.ext.githubpages',
]

templates_path = ['_templates']
exclude_patterns = []

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']

autodoc_member_order = 'bysource'
add_module_names = False 
