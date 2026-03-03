# This shim allows `pip install -e .` to work with older pip versions (<22).
# With pip >=22, pyproject.toml alone is sufficient.
from setuptools import setup
setup()
