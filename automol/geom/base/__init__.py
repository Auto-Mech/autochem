""" Level 2 geometry functions (no dependencies on extern or other types)

Import hierarchy:
    _core       no dependencies
    _comp       dependencies: _core
"""

# core functions
# # constructors
from automol.geom.base._core import from_data
from automol.geom.base._core import from_subset
# # getters
from automol.geom.base._core import symbols
from automol.geom.base._core import coordinates
# # setters
from automol.geom.base._core import set_coordinates
# # I/O
from automol.geom.base._core import string
from automol.geom.base._core import xyz_string
from automol.geom.base._core import xyz_trajectory_string
from automol.geom.base._core import from_string
from automol.geom.base._core import from_xyz_string
from automol.geom.base._core import xyz_string_comment
from automol.geom.base._core import from_xyz_trajectory_string
# # validation
from automol.geom.base._core import is_valid
# # conversions
from automol.geom.base._core import formula
from automol.geom.base._core import formula_string
# # properties
from automol.geom.base._core import count
from automol.geom.base._core import is_atom
from automol.geom.base._core import is_diatomic
from automol.geom.base._core import is_linear
from automol.geom.base._core import atom_count
from automol.geom.base._core import atom_indices
from automol.geom.base._core import dummy_atom_indices
from automol.geom.base._core import masses
from automol.geom.base._core import total_mass
from automol.geom.base._core import center_of_mass
from automol.geom.base._core import mass_centered
from automol.geom.base._core import reduced_mass
from automol.geom.base._core import inertia_tensor
from automol.geom.base._core import principal_axes
from automol.geom.base._core import moments_of_inertia
from automol.geom.base._core import rotational_constants
# # geometric measurements
from automol.geom.base._core import distance
from automol.geom.base._core import central_angle
from automol.geom.base._core import dihedral_angle
from automol.geom.base._core import zmatrix_row_values
# # binary functions
from automol.geom.base._core import join
from automol.geom.base._core import minimum_distance
from automol.geom.base._core import permutation
# # adding/removing atoms
from automol.geom.base._core import insert
from automol.geom.base._core import remove
from automol.geom.base._core import without_dummy_atoms
from automol.geom.base._core import reorder
from automol.geom.base._core import move_atom
from automol.geom.base._core import swap_coordinates
# # transformations
from automol.geom.base._core import round_
from automol.geom.base._core import translate
from automol.geom.base._core import translate_along_matrix
from automol.geom.base._core import perturb
from automol.geom.base._core import rotate
from automol.geom.base._core import euler_rotate
from automol.geom.base._core import transform
from automol.geom.base._core import transform_by_matrix
from automol.geom.base._core import reflect_coordinates
from automol.geom.base._core import shift_atom_position
# comparison functions
# # properties used for comparisons
from automol.geom.base._comp import coulomb_spectrum
from automol.geom.base._comp import distance_matrix
# # comparisons
from automol.geom.base._comp import almost_equal
from automol.geom.base._comp import almost_equal_coulomb_spectrum
from automol.geom.base._comp import argunique_coulomb_spectrum
from automol.geom.base._comp import almost_equal_dist_matrix
from automol.geom.base._comp import minimum_volume_geometry


__all__ = [
    # core functions
    # # constructors
    'from_data',
    'from_subset',
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
    'atom_indices',
    'dummy_atom_indices',
    'masses',
    'total_mass',
    'center_of_mass',
    'mass_centered',
    'reduced_mass',
    'inertia_tensor',
    'principal_axes',
    'moments_of_inertia',
    'rotational_constants',
    # # geometric measurements
    'distance',
    'central_angle',
    'dihedral_angle',
    'zmatrix_row_values',
    # # binary functions
    'join',
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
    'euler_rotate',
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
]
