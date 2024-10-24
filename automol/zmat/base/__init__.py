""" Level 2 Z-Matrix functions

Import hierarchy:
    _core       no dependencies
    _extra      dependencies: _core
"""

# common v-matrix/z-matrix functions:
from ...vmat import symbols
from ...vmat import key_matrix
from ...vmat import name_matrix
from ...vmat import count
from ...vmat import keys
from ...vmat import atom_indices
from ...vmat import coordinate_key_matrix
from ...vmat import coordinates
from ...vmat import distance_coordinates
from ...vmat import central_angle_coordinates
from ...vmat import dihedral_angle_coordinates
from ...vmat import coordinate
from ...vmat import torsion_axis
from ...vmat import names
from ...vmat import distance_names
from ...vmat import central_angle_names
from ...vmat import dihedral_angle_names
from ...vmat import angle_names
from ...vmat import standard_names
from ...vmat import standard_name_matrix
from ...vmat import distance_coordinate_name
from ...vmat import central_angle_coordinate_name
from ...vmat import dihedral_angle_coordinate_name
from ...vmat import dummy_keys
from ...vmat import dummy_source_dict
from ...vmat import conversion_info
from ...vmat import neighbor_keys
# core functions
# # constructors
from ._core import from_data
from ._core import from_geometry
# # getters
from ._core import value_matrix
from ._core import value_dictionary
from ._core import value
# # setters
from ._core import set_key_matrix
from ._core import set_name_matrix
from ._core import set_value_matrix
from ._core import set_values_by_name
# # I/O
from ._core import string
from ._core import from_string
from ._core import yaml_data
from ._core import from_yaml_data
# # validation
from ._core import is_valid
# # properties
from ._core import torsion_coordinates
# # conversions
from ._core import vmatrix
from ._core import formula
# # relabelling
from ._core import rename
from ._core import standard_form
from ._core import round_
# # add/remove atoms
from ._core import add_atom
from ._core import remove_atom
# # comparisons
from ._core import almost_equal
# extra functions
from ._extra import samples
from ._extra import constraint_dict
from ._extra import set_constraint_names


__all__ = [
    # common v-matrix/z-matrix functions:
    'symbols',
    'key_matrix',
    'name_matrix',
    'count',
    'keys',
    'atom_indices',
    'coordinate_key_matrix',
    'coordinates',
    'distance_coordinates',
    'central_angle_coordinates',
    'dihedral_angle_coordinates',
    'coordinate',
    'torsion_axis',
    'names',
    'distance_names',
    'central_angle_names',
    'dihedral_angle_names',
    'angle_names',
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
    'neighbor_keys',
    # extra functions
    'samples',
    'constraint_dict',
    'set_constraint_names',
]
