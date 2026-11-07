"""Sphinx settings for the ASPCAP user guide."""

from importlib.metadata import version as package_version

project = "ASPCAP"
release = package_version("aspcap")
extensions = ["sphinx.ext.autodoc", "sphinx.ext.napoleon"]
exclude_patterns = ["_build"]
html_theme = "alabaster"
