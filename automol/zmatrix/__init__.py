""" a z-matrix internal geometry module
"""
from .._cnst.zmatrix import from_data
from ._core import from_matrices
from ._core import symbols
from ._core import distance_column
from ._core import angle_column
from ._core import torsion_column
from ._core import coordinate_matrix
from ._core import key_matrix
from ._core import value_matrix
from ._core import iter_
from ._core import iter_column_
from ._math import almost_equal
from ._io import from_zmat_string
from ._io import zmat_string
from ._io import zmat_string_matrix_block
from ._io import zmat_string_variable_block
from . import parse

__all__ = [
    'from_data',
    'from_matrices',
    'symbols',
    'distance_column',
    'angle_column',
    'torsion_column',
    'coordinate_matrix',
    'key_matrix',
    'value_matrix',
    'iter_',
    'iter_column_',
    'almost_equal',
    'from_zmat_string',
    'zmat_string',
    'zmat_string_matrix_block',
    'zmat_string_variable_block',
    'parse',
]
