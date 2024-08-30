""" Level 2 geometry functions (no dependencies on extern or other types)

Import hierarchy:
    _core       no dependencies
    _comp       dependencies: _core
"""

# core functions
# # constructors
from ._0core import from_data
from ._0core import subgeom
# # getters
from ._0core import symbols
from ._0core import coordinates
# # setters
from ._0core import set_coordinates
# # I/O
from ._0core import string
from ._0core import xyz_string
from ._0core import xyz_trajectory_string
from ._0core import from_string
from ._0core import from_xyz_string
from ._0core import xyz_string_comment
from ._0core import from_xyz_trajectory_string
from ._0core import yaml_data
from ._0core import from_yaml_data
# # validation
from ._0core import is_valid
# # conversions
from ._0core import formula
from ._0core import formula_string
# # properties
from ._0core import count
from ._0core import is_atom
from ._0core import is_diatomic
from ._0core import is_linear
from ._0core import atom_count
from ._0core import electron_count
from ._0core import atom_indices
from ._0core import dummy_atom_indices
from ._0core import masses
from ._0core import mass_weight_vector
from ._0core import total_mass
from ._0core import center_of_mass
from ._0core import mass_centered
from ._0core import aligned_to_principal_axes
from ._0core import reduced_mass
from ._0core import inertia_tensor
from ._0core import principal_axes
from ._0core import moments_of_inertia
from ._0core import rotational_constants
from ._0core import rotational_analysis
from ._0core import translational_normal_modes
from ._0core import rotational_normal_modes
# # geometric measurements
from ._0core import distance
from ._0core import central_angle
from ._0core import dihedral_angle
from ._0core import measure
from ._0core import zmatrix_row_values
# # binary functions
from ._0core import join
from ._0core import join_sequence
from ._0core import minimum_distance
from ._0core import permutation
# # adding/removing atoms
from ._0core import insert
from ._0core import remove
from ._0core import without_dummy_atoms
from ._0core import reorder
from ._0core import move_atom
from ._0core import swap_coordinates
# # transformations
from ._0core import round_
from ._0core import translate
from ._0core import translate_along_matrix
from ._0core import perturb
from ._0core import rotate
from ._0core import transform
from ._0core import transform_by_matrix
from ._0core import reflect_coordinates
from ._0core import shift_atom_position
# comparison functions
# # properties used for comparisons
from ._1comp import coulomb_spectrum
from ._1comp import distance_matrix
# # comparisons
from ._1comp import almost_equal
from ._1comp import almost_equal_coulomb_spectrum
from ._1comp import argunique_coulomb_spectrum
from ._1comp import almost_equal_dist_matrix
from ._1comp import minimum_volume_geometry
# intermolecular interactions
from ._2intmol import has_low_relative_repulsion_energy
from ._2intmol import total_repulsion_energy
from ._2intmol import repulsion_energy


__all__ = [
    # core functions
    # # constructors
    'from_data',
    'subgeom',
    # # getters
    'symbols',
    'coordinates',
    # # setters
    'set_coordinates',
    # # I/O
    'string',
    'xyz_string',
    'xyz_trajectory_string',
    'from_string',
    'from_xyz_string',
    'xyz_string_comment',
    'from_xyz_trajectory_string',
    'yaml_data',
    'from_yaml_data',
    # # validation
    'is_valid',
    # # conversions
    'formula',
    'formula_string',
    # # properties
    'count',
    'is_atom',
    'is_diatomic',
    'is_linear',
    'atom_count',
    'electron_count',
    'atom_indices',
    'dummy_atom_indices',
    'masses',
    'mass_weight_vector',
    'total_mass',
    'center_of_mass',
    'mass_centered',
    'aligned_to_principal_axes',
    'reduced_mass',
    'inertia_tensor',
    'principal_axes',
    'moments_of_inertia',
    'rotational_constants',
    'rotational_analysis',
    'translational_normal_modes',
    'rotational_normal_modes',
    # # geometric measurements
    'distance',
    'central_angle',
    'dihedral_angle',
    'measure',
    'zmatrix_row_values',
    # # binary functions
    'join',
    'join_sequence',
    'minimum_distance',
    'permutation',
    # # adding/removing atoms
    'insert',
    'remove',
    'without_dummy_atoms',
    'reorder',
    'move_atom',
    'swap_coordinates',
    # # transformations
    'round_',
    'translate',
    'translate_along_matrix',
    'perturb',
    'rotate',
    'transform',
    'transform_by_matrix',
    'reflect_coordinates',
    'shift_atom_position',
    # comparison functions
    # # properties used for comparisons
    'coulomb_spectrum',
    'distance_matrix',
    # # comparisons
    'almost_equal',
    'almost_equal_coulomb_spectrum',
    'argunique_coulomb_spectrum',
    'almost_equal_dist_matrix',
    'minimum_volume_geometry',
    # intermolecular interactions
    'has_low_relative_repulsion_energy',
    'total_repulsion_energy',
    'repulsion_energy',
]
