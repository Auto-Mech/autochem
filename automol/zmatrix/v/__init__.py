""" variable z-matrix
"""
from automol.constructors.vmatrix import from_data
from automol.zmatrix.v._core import symbols
from automol.zmatrix.v._core import key_matrix
from automol.zmatrix.v._core import name_matrix
from automol.zmatrix.v._core import coordinate_key_matrix
from automol.zmatrix.v._core import coordinates
from automol.zmatrix.v._core import names
from automol.zmatrix.v._core import distance_names
from automol.zmatrix.v._core import central_angle_names
from automol.zmatrix.v._core import dihedral_angle_names
from automol.zmatrix.v._core import angle_names
from automol.zmatrix.v._core import set_names
from automol.zmatrix.v._core import standard_names
from automol.zmatrix.v._core import standard_form
from automol.zmatrix.v._core import is_valid
from automol.zmatrix.v._core import is_standard_form
from automol.zmatrix.v._io import from_string
from automol.zmatrix.v._io import string

__all__ = [
    'from_data',
    'symbols',
    'key_matrix',
    'name_matrix',
    'coordinate_key_matrix',
    'coordinates',
    'names',
    'distance_names',
    'central_angle_names',
    'dihedral_angle_names',
    'angle_names',
    'set_names',
    'standard_names',
    'standard_form',
    'is_valid',
    'is_standard_form',
    'from_string',
    'string',
]
