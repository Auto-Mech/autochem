""" a z-matrix internal geometry module
"""
from automol.constructors.zmatrix import from_data
from automol.zmatrix._core import var_
from automol.zmatrix._core import symbols
from automol.zmatrix._core import key_matrix
from automol.zmatrix._core import name_matrix
from automol.zmatrix._core import value_matrix
from automol.zmatrix._core import coordinate_key_matrix
from automol.zmatrix._core import coordinates
from automol.zmatrix._core import names
from automol.zmatrix._core import distance_names
from automol.zmatrix._core import central_angle_names
from automol.zmatrix._core import dihedral_angle_names
from automol.zmatrix._core import angle_names
from automol.zmatrix._core import values
from automol.zmatrix._core import set_names
from automol.zmatrix._core import set_values
from automol.zmatrix._core import standard_names
from automol.zmatrix._core import standard_form
from automol.zmatrix._core import is_valid
from automol.zmatrix._core import is_standard_form
from automol.zmatrix._io import from_string
from automol.zmatrix._io import string
from automol.zmatrix._geom import geometry
from automol.zmatrix._graph import connectivity_graph
from automol.zmatrix._comp import almost_equal
# submodules
from automol.zmatrix import v
from automol.zmatrix import tors

__all__ = [
    'from_data',
    'var_',
    'symbols',
    'key_matrix',
    'name_matrix',
    'value_matrix',
    'coordinate_key_matrix',
    'coordinates',
    'names',
    'distance_names',
    'central_angle_names',
    'dihedral_angle_names',
    'angle_names',
    'values',
    'set_names',
    'set_values',
    'standard_names',
    'standard_form',
    'is_valid',
    'is_standard_form',
    'from_string',
    'string',
    'geometry',
    'connectivity_graph',
    'almost_equal',
    # submodules
    'v',
    'tors',
]
