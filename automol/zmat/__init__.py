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
# # repulsion energy
from automol.zmat._conv import has_low_relative_repulsion_energy
from automol.zmat._conv import total_repulsion_energy
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
    # extra base functions
    'samples',
    'constraint_dict',
    'set_constraint_names',
    'coordinate',
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
    # # repulsion energy
    'has_low_relative_repulsion_energy',
    'total_repulsion_energy',
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
