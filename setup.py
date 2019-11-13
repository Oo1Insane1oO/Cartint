"""
Installer for cartint
"""
from setuptools import setup, find_packages

from setuptools.command.develop import develop
from setuptools.command.install import install

setup(
    name="cartint",
    version="0.1a",
    description="Solver for integrals over Gauss-Hermite functions in cartesian coordinates",
    author="Alocias Mariadason",
    package_dir={"": "interface"},
    packages=find_packages("interface"),
    entry_points={
        "console_scripts": [
            "cartrun=cartrun.main:run",
        ]
    },
    zip_safe=False,
)
