"""
Read yaml-file and return dictionary with input parameters

example:

    >>> from tempfile import NamedTemporaryFile

    >>> with NamedTemporaryFile(mode="w+", suffix=".yaml") as tmpfile:
    ...     data = {
    ...         "foo" : "A",
    ...         "dim" : 2,
    ...         "basis" : [
    ...             [0, 0, 1, 1, -1, 0],
    ...             [0, 0, 1, 1,  1, 0]
    ...         ],
    ...         "quad" : {
    ...             "method" : "Chebychev",
    ...             "N" : 50,
    ...         }
    ...     }
    ...     yaml.dump(data, tmpfile)
    ...     _ = tmpfile.seek(0,0)
    ...     yaml_data = read_yaml(tmpfile.name)
    >>> for k in sorted(yaml_data.keys()):
    ...     print(f"{k} : {yaml_data[k]}")
    basis : [[0, 0, 1, 1, -1, 0], [0, 0, 1, 1, 1, 0]]
    dim : 2
    foo : A
    quad : {'N': 50, 'method': 'Chebychev'}
"""
from typing import Union, Dict, Any
from pathlib import Path

import yaml
try:
    from yaml import CSafeLoader as SafeLoader
except ImportError:
    from yaml import SafeLoader

def read_yaml(filepath: Union[Path, str]) -> Dict[str, Any]:
    """
    Read yaml file and return dictionary with input parameters

    Args:
        filepath:
            path to yaml file to read
    """
    with open(filepath, "r") as file_stream:
        return yaml.load(file_stream, Loader=SafeLoader)
