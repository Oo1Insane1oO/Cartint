"""
Read input yaml file with specs and run integration
"""
from typing import Union
from pathlib import Path

from .read_yaml import read_yaml

def run(filepath: Union[Path, str]):
    """
    Reads yaml file with specs, initializes pybind11 integral modules and calls specific integrate
    method.
    """
    parameters = read_yaml(filepath)
