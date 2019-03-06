""" a z-matrix internal geometry module
"""
from ..constructors.zmatrix import from_data
from ._core import matrix
from ._core import values
from ._core import symbols
from ._core import key_matrix
from ._core import name_matrix
from ._core import value_matrix
from ._core import coordinate_key_matrix
from ._core import distance_names
from ._core import angle_names
from ._core import dihedral_names
from ._core import distance_keys
from ._core import angle_keys
from ._core import dihedral_keys
from ._core import set_names
from ._core import set_values
from ._io import from_zmat_string
from ._io import zmat_string
from ._io import matrix_block_string
from ._io import setval_block_string
from ._comp import almost_equal
from ._geom import geometry
from ._graph import connectivity_graph
# submodules
from . import tors

__all__ = [
    'from_data',
    'matrix',
    'values',
    'symbols',
    'key_matrix',
    'name_matrix',
    'value_matrix',
    'coordinate_key_matrix',
    'distance_names',
    'angle_names',
    'dihedral_names',
    'distance_keys',
    'angle_keys',
    'dihedral_keys',
    'set_names',
    'set_values',
    'from_zmat_string',
    'zmat_string',
    'matrix_block_string',
    'setval_block_string',
    'almost_equal',
    'geometry',
    'connectivity_graph',
    # submodules
    'tors',
]
