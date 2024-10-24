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
from ..vmat import symbols
from ..vmat import key_matrix
from ..vmat import name_matrix
from ..vmat import count
from ..vmat import keys
from ..vmat import atom_indices
from ..vmat import coordinate_key_matrix
from ..vmat import coordinates
from ..vmat import distance_coordinates
from ..vmat import central_angle_coordinates
from ..vmat import dihedral_angle_coordinates
from ..vmat import coordinate
from ..vmat import torsion_axis
from ..vmat import names
from ..vmat import distance_names
from ..vmat import central_angle_names
from ..vmat import dihedral_angle_names
from ..vmat import angle_names
from ._conv import dummy_coordinate_names
from ..vmat import standard_names
from ..vmat import standard_name_matrix
from ..vmat import distance_coordinate_name
from ..vmat import central_angle_coordinate_name
from ..vmat import dihedral_angle_coordinate_name
from ..vmat import dummy_keys
from ..vmat import dummy_source_dict
from ..vmat import conversion_info
from ..vmat import neighbor_keys
# core functions
# # constructors
from .base._core import from_data
from .base._core import from_geometry
# # getters
from .base._core import value_matrix
from .base._core import value_dictionary
from .base._core import value
# # setters
from .base._core import set_key_matrix
from .base._core import set_name_matrix
from .base._core import set_value_matrix
from .base._core import set_values_by_name
# # I/O
from .base._core import string
from .base._core import from_string
from .base._core import yaml_data
from .base._core import from_yaml_data
# # validation
from .base._core import is_valid
# # properties
from .base._core import torsion_coordinates
# # conversions
from .base._core import vmatrix
from .base._core import formula
# # relabelling
from .base._core import rename
from .base._core import standard_form
from .base._core import round_
# # add/remove atoms
from .base._core import add_atom
from .base._core import remove_atom
# # comparisons
from .base._core import almost_equal
# extra functions
from .base._extra import samples
from .base._extra import constraint_dict
from .base._extra import set_constraint_names
# L4
# conversion functions
# # conversions
from ._conv import graph
from ._conv import graph_without_stereo
from ._conv import geometry
from ._conv import geometry_with_conversion_info
from ._conv import rdkit_molecule
from ._conv import py3dmol_view
from ._conv import display
# # derived properties
from ._conv import distance
from ._conv import central_angle
from ._conv import dihedral_angle
# # torsions
from ._conv import torsion_coordinate_name
from ._conv import torsion_leading_atom
# # repulsion energy
from ._conv import has_low_relative_repulsion_energy
from ._conv import total_repulsion_energy
# ring functions:
from ._ring import all_rings_atoms
from ._ring import all_rings_distances
from ._ring import all_rings_dihedrals
from ._ring import all_rings_dct
from ._ring import all_rings_distances_reasonable
from ._ring import ring_distances
from ._ring import ring_dihedrals
from ._ring import ring_distances_reasonable
from ._ring import complete_ring_dihedrals
from ._ring import samples_avg_dih



__all__ = [
    # L2
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
    'neighbor_keys',
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
    'complete_ring_dihedrals',
    'ring_distances_reasonable',
    'samples_avg_dih',
]
