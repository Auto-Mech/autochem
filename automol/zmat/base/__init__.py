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
from automol.vmat import names
from automol.vmat import distance_names
from automol.vmat import central_angle_names
from automol.vmat import dihedral_angle_names
from automol.vmat import angle_names
from automol.vmat import dummy_coordinate_names
from automol.vmat import standard_names
from automol.vmat import standard_name_matrix
# core functions
# # constructors
from automol.zmat.base._core import from_data
from automol.zmat.base._core import from_geometry
# # getters
from automol.zmat.base._core import value_matrix
from automol.zmat.base._core import value_dictionary
# # setters
from automol.zmat.base._core import set_key_matrix
from automol.zmat.base._core import set_name_matrix
from automol.zmat.base._core import set_value_matrix
from automol.zmat.base._core import set_values_by_name
# # I/O
from automol.zmat.base._core import string
from automol.zmat.base._core import from_string
# # validation
from automol.zmat.base._core import is_valid
# # conversions
from automol.zmat.base._core import vmatrix
from automol.zmat.base._core import formula
# # relabelling
from automol.zmat.base._core import rename
from automol.zmat.base._core import standard_form
# # add/remove atoms
from automol.zmat.base._core import add_atom
from automol.zmat.base._core import remove_atom
# # comparisons
from automol.zmat.base._core import almost_equal
# # coordinate names
from automol.zmat.base._core import distance_coordinate_name
from automol.zmat.base._core import central_angle_coordinate_name
from automol.zmat.base._core import dihedral_angle_coordinate_name
from automol.zmat.base._core import dihedral_axis_name
# # dummy atom functions
from automol.zmat.base._core import dummy_keys
from automol.zmat.base._core import dummy_key_dictionary
from automol.zmat.base._core import dummy_neighbor_keys
from automol.zmat.base._core import linear_atom_keys
from automol.zmat.base._core import shift_down
from automol.zmat.base._core import shift_up
# extra functions
from automol.zmat.base._extra import samples
from automol.zmat.base._extra import torsional_sampling_ranges
from automol.zmat.base._extra import constraint_dct
from automol.zmat.base._extra import set_constraint_names
from automol.zmat.base._extra import coord_idxs
from automol.zmat.base._extra import bond_key_from_idxs
from automol.zmat.base._extra import get_babs1
from automol.zmat.base._extra import get_babs2


__all__ = [
    # common v-matrix/z-matrix functions:
    'symbols',
    'key_matrix',
    'name_matrix',
    'count',
    'atom_indices',
    'coordinate_key_matrix',
    'coordinates',
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
    # # setters
    'set_key_matrix',
    'set_name_matrix',
    'set_value_matrix',
    'set_values_by_name',
    # # I/O
    'string',
    'from_string',
    # # validation
    'is_valid',
    # # conversions
    'vmatrix',
    'formula',
    # # relabelling
    'rename',
    'standard_form',
    # # add/remove atoms
    'add_atom',
    'remove_atom',
    # # comparisons
    'almost_equal',
    # # coordinate names
    'distance_coordinate_name',
    'central_angle_coordinate_name',
    'dihedral_angle_coordinate_name',
    'dihedral_axis_name',
    # # dummy atom functions
    'dummy_keys',
    'dummy_key_dictionary',
    'dummy_neighbor_keys',
    'linear_atom_keys',
    'shift_down',
    'shift_up',
    # extra functions
    'samples',
    'torsional_sampling_ranges',
    'constraint_dct',
    'set_constraint_names',
    'coord_idxs',
    'bond_key_from_idxs',
    'get_babs1',
    'get_babs2',
]
