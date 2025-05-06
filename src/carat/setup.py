"""Package setup for CarAT.

Defines metadata and dependencies for installing the CarAT library.
"""

from pathlib import Path

from setuptools import find_packages, setup

# Read dependencies
here = Path(__file__).parent
install_requires = here.joinpath("requirements.txt").read_text().splitlines()

setup(
    name="my-carat",
    version="0.1.0",
    package_dir={"": "src"},
    packages=["carat"] + [f"carat.{pkg}" for pkg in find_packages(where="src/carat")],
    install_requires=install_requires,
)
