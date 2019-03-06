""" a z-matrix internal geometry module
"""
from ..constructors.zmatrix import from_data
from ._core import matrix
from ._core import values
from ._core import symbols
from ._core import key_matrix
from ._core import name_matrix
from ._core import distance_names
from ._core import angle_names
from ._core import dihedral_names
from ._core import set_names
from ._core import set_values
from ._math import almost_equal
from ._geom import geometry
from ._io import from_zmat_string
from ._io import zmat_string
from ._io import matrix_block_string
from ._io import setval_block_string

__all__ = [
    'from_data',
    'matrix',
    'values',
    'symbols',
    'key_matrix',
    'name_matrix',
    'distance_names',
    'angle_names',
    'dihedral_names',
    'set_names',
    'set_values',
    'almost_equal',
    'geometry',
    'from_zmat_string',
    'zmat_string',
    'matrix_block_string',
    'setval_block_string',
]
