""" Level 2 Z-Matrix functions

Import hierarchy:
    _core       no dependencies
    _extra      dependencies: _core
"""

# common v-matrix/z-matrix functions:
from automol.vmat import symbols
from automol.vmat import key_matrix
from automol.vmat import name_matrix
from automol.vmat import count
from automol.vmat import atom_indices
from automol.vmat import coordinate_key_matrix
from automol.vmat import coordinates
from automol.vmat import coordinate
from automol.vmat import torsion_axis
from automol.vmat import names
from automol.vmat import distance_names
from automol.vmat import central_angle_names
from automol.vmat import dihedral_angle_names
from automol.vmat import angle_names
from automol.vmat import dummy_coordinate_names
from automol.vmat import standard_names
from automol.vmat import standard_name_matrix
from automol.vmat import distance_coordinate_name
from automol.vmat import central_angle_coordinate_name
from automol.vmat import dihedral_angle_coordinate_name
from automol.vmat import dummy_keys
from automol.vmat import dummy_source_dict
from automol.vmat import conversion_info
# core functions
# # constructors
from automol.zmat.base._core import from_data
from automol.zmat.base._core import from_geometry
# # getters
from automol.zmat.base._core import value_matrix
from automol.zmat.base._core import value_dictionary
from automol.zmat.base._core import value
# # setters
from automol.zmat.base._core import set_key_matrix
from automol.zmat.base._core import set_name_matrix
from automol.zmat.base._core import set_value_matrix
from automol.zmat.base._core import set_values_by_name
# # I/O
from automol.zmat.base._core import string
from automol.zmat.base._core import from_string
from automol.zmat.base._core import yaml_data
from automol.zmat.base._core import from_yaml_data
# # validation
from automol.zmat.base._core import is_valid
# # properties
from automol.zmat.base._core import torsion_coordinates
# # conversions
from automol.zmat.base._core import vmatrix
from automol.zmat.base._core import formula
# # relabelling
from automol.zmat.base._core import rename
from automol.zmat.base._core import standard_form
from automol.zmat.base._core import round_
# # add/remove atoms
from automol.zmat.base._core import add_atom
from automol.zmat.base._core import remove_atom
# # comparisons
from automol.zmat.base._core import almost_equal
# extra functions
from automol.zmat.base._extra import samples
from automol.zmat.base._extra import constraint_dict
from automol.zmat.base._extra import set_constraint_names


__all__ = [
    # common v-matrix/z-matrix functions:
    'symbols',
    'key_matrix',
    'name_matrix',
    'count',
    'atom_indices',
    'coordinate_key_matrix',
    'coordinates',
    'coordinate',
    'torsion_axis',
    'names',
    'distance_names',
    'central_angle_names',
    'dihedral_angle_names',
    'angle_names',
    'dummy_coordinate_names',
    'standard_names',
    'standard_name_matrix',
    # core functions
    # # constructors
    'from_data',
    'from_geometry',
    # # getters
    'value_matrix',
    'value_dictionary',
    'value',
    # # setters
    'set_key_matrix',
    'set_name_matrix',
    'set_value_matrix',
    'set_values_by_name',
    # # I/O
    'string',
    'from_string',
    'yaml_data',
    'from_yaml_data',
    # # validation
    'is_valid',
    # # properties
    'torsion_coordinates',
    # # conversions
    'vmatrix',
    'formula',
    # # relabelling
    'rename',
    'standard_form',
    'round_',
    # # add/remove atoms
    'add_atom',
    'remove_atom',
    # # comparisons
    'almost_equal',
    # # coordinate names
    'distance_coordinate_name',
    'central_angle_coordinate_name',
    'dihedral_angle_coordinate_name',
    # # dummy atom functions
    'dummy_keys',
    'dummy_source_dict',
    'conversion_info',
    # extra functions
    'samples',
    'constraint_dict',
    'set_constraint_names',
]
