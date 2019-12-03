"""
Installer for cartint
"""
import os
import subprocess
from pathlib import Path
from setuptools import setup, find_packages, Extension

class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=''):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)

os.chdir(Path().cwd())

CMAKE_STATUS = subprocess.run(
    ["poetry", "run", "env", "CMAKE_BUILD_PARALLEL_LEVEL=$(nproc --all)", "cmake", ".", "-B", "build"],
    capture_output=True,
)

assert CMAKE_STATUS.returncode > -1, CMAKE_STATUS.stdout

MAKE_STATUS = subprocess.run(["poetry", "run", "make", "-C", "build"])

assert MAKE_STATUS.returncode > -1, MAKE_STATUS.stdout

setup(
    name="cartint",
    version="0.1a",
    description="Solver for integrals over Gauss-Hermite functions in cartesian coordinates",
    author="Alocias Mariadason",
    package_dir={"": "interface"},
    packages=find_packages("interface"),
    ext_modules=[CMakeExtension("integral", "binding")],
    entry_points={
        "console_scripts": [
            "cartrun=cartrun.main:run",
        ]
    },
    zip_safe=False,
)
