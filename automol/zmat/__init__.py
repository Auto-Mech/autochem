""" Z-Matrices

Level 2 Z-Matrix functions belong in geom/base/*.py. These do not require
converion to other basic types (geom, graph, zmat, inchi).

Level 4 Z-Matrix functions belong in geom/*.py. These **do** require
converion to other basic types (geom, graph, zmat, inchi).

Import hierarchy:
    _conv       no dependencies
    _extra      dependencies: _conv
"""

# L2
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
from automol.zmat.base._core import round_
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
from automol.zmat.base._core import dummy_parent_keys
from automol.zmat.base._core import dummy_conversion
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
# L4
# conversion functions
# # conversions
from automol.zmat._conv import graph
from automol.zmat._conv import graph_without_stereo
from automol.zmat._conv import geometry
from automol.zmat._conv import geometry_with_conversion_info
from automol.zmat._conv import rdkit_molecule
from automol.zmat._conv import py3dmol_view
from automol.zmat._conv import display
# # derived properties
from automol.zmat._conv import distance
from automol.zmat._conv import central_angle
from automol.zmat._conv import dihedral_angle
# # torsions
from automol.zmat._conv import torsion_coordinate_name
from automol.zmat._conv import torsion_leading_atom
# extra functions:
from automol.zmat._extra import is_atom_closest_to_bond_atom
# ring functions:
from automol.zmat._ring import all_rings_atoms
from automol.zmat._ring import all_rings_distances
from automol.zmat._ring import all_rings_dihedrals
from automol.zmat._ring import all_rings_dct
from automol.zmat._ring import all_rings_distances_reasonable
from automol.zmat._ring import ring_distances
from automol.zmat._ring import ring_dihedrals
from automol.zmat._ring import ring_distances_reasonable


__all__ = [
    # L2
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
    'dihedral_axis_name',
    # # dummy atom functions
    'dummy_keys',
    'dummy_key_dictionary',
    'dummy_parent_keys',
    'dummy_conversion',
    'linear_atom_keys',
    'shift_down',
    'shift_up',
    # extra base functions
    'samples',
    'torsional_sampling_ranges',
    'constraint_dct',
    'set_constraint_names',
    'coord_idxs',
    'bond_key_from_idxs',
    # L4
    # conversion functions
    # # conversions
    'graph',
    'graph_without_stereo',
    'geometry',
    'geometry_with_conversion_info',
    'rdkit_molecule',
    'py3dmol_view',
    'display',
    # # derived properties
    'distance',
    'central_angle',
    'dihedral_angle',
    # # torsions
    'torsion_coordinate_name',
    'torsion_leading_atom',
    # extra functions:
    'is_atom_closest_to_bond_atom',
    # ring functions
    'all_rings_atoms',
    'all_rings_distances',
    'all_rings_dihedrals',
    'all_rings_dct',
    'all_rings_distances_reasonable',
    'ring_distances',
    'ring_dihedrals',
    'ring_distances_reasonable'
]
