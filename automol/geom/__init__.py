""" cartesian geometries

Level 2 geometry functions belong in geom/base/*.py. These do not require
converion to other basic types (geom, graph, zmat, inchi).

Level 4 geometry functions belong in geom/*.py. These **do** require
converion to other basic types (geom, graph, zmat, inchi).

Import hierarchy:
    _conv       dependencies: automol.graph
    _extra      dependencies: automol.graph, _conv
    ts          dependencies: automol.graph
"""

# L2
# core functions
# # constructors
from automol.geom.base._0core import from_data
from automol.geom.base._0core import subgeom
# # getters
from automol.geom.base._0core import symbols
from automol.geom.base._0core import coordinates
# # setters
from automol.geom.base._0core import set_coordinates
# # I/O
from automol.geom.base._0core import string
from automol.geom.base._0core import xyz_string
from automol.geom.base._0core import xyz_trajectory_string
from automol.geom.base._0core import from_string
from automol.geom.base._0core import from_xyz_string
from automol.geom.base._0core import xyz_string_comment
from automol.geom.base._0core import from_xyz_trajectory_string
from automol.geom.base._0core import yaml_data
from automol.geom.base._0core import from_yaml_data
# # validation
from automol.geom.base._0core import is_valid
# # conversions
from automol.geom.base._0core import formula
from automol.geom.base._0core import formula_string
# # properties
from automol.geom.base._0core import count
from automol.geom.base._0core import is_atom
from automol.geom.base._0core import is_diatomic
from automol.geom.base._0core import is_linear
from automol.geom.base._0core import atom_count
from automol.geom.base._0core import atom_indices
from automol.geom.base._0core import dummy_atom_indices
from automol.geom.base._0core import masses
from automol.geom.base._0core import total_mass
from automol.geom.base._0core import center_of_mass
from automol.geom.base._0core import mass_centered
from automol.geom.base._0core import reduced_mass
from automol.geom.base._0core import inertia_tensor
from automol.geom.base._0core import principal_axes
from automol.geom.base._0core import moments_of_inertia
from automol.geom.base._0core import rotational_constants
# # geometric measurements
from automol.geom.base._0core import distance
from automol.geom.base._0core import central_angle
from automol.geom.base._0core import dihedral_angle
from automol.geom.base._0core import zmatrix_row_values
# # binary functions
from automol.geom.base._0core import join
from automol.geom.base._0core import minimum_distance
from automol.geom.base._0core import permutation
# # adding/removing atoms
from automol.geom.base._0core import insert
from automol.geom.base._0core import remove
from automol.geom.base._0core import without_dummy_atoms
from automol.geom.base._0core import reorder
from automol.geom.base._0core import move_atom
from automol.geom.base._0core import swap_coordinates
# # transformations
from automol.geom.base._0core import round_
from automol.geom.base._0core import translate
from automol.geom.base._0core import translate_along_matrix
from automol.geom.base._0core import perturb
from automol.geom.base._0core import rotate
from automol.geom.base._0core import transform
from automol.geom.base._0core import transform_by_matrix
from automol.geom.base._0core import reflect_coordinates
from automol.geom.base._0core import shift_atom_position
# comparison functions
# # properties used for comparisons
from automol.geom.base._1comp import coulomb_spectrum
from automol.geom.base._1comp import distance_matrix
# # comparisons
from automol.geom.base._1comp import almost_equal
from automol.geom.base._1comp import almost_equal_coulomb_spectrum
from automol.geom.base._1comp import argunique_coulomb_spectrum
from automol.geom.base._1comp import almost_equal_dist_matrix
from automol.geom.base._1comp import minimum_volume_geometry
# intermolecular interactions
from automol.geom.base._2intmol import has_low_relative_repulsion_energy
from automol.geom.base._2intmol import total_repulsion_energy
from automol.geom.base._2intmol import repulsion_energy
# L4
# MolSym interface
from automol.geom._molsym import point_group_from_geometry
# conversion functions:
# # conversions
from automol.geom._conv import graph
from automol.geom._conv import graph_without_stereo
from automol.geom._conv import connectivity_graph_deprecated
from automol.geom._conv import zmatrix
from automol.geom._conv import zmatrix_with_conversion_info
from automol.geom._conv import update_zmatrix
from automol.geom._conv import amchi
from automol.geom._conv import amchi_with_numbers
from automol.geom._conv import inchi
from automol.geom._conv import inchi_with_numbers
from automol.geom._conv import chi
from automol.geom._conv import chi_with_sort
from automol.geom._conv import smiles
from automol.geom._conv import rdkit_molecule
from automol.geom._conv import py3dmol_view
from automol.geom._conv import display
# # derived properties
from automol.geom._conv import is_connected
from automol.geom._conv import linear_atoms
from automol.geom._conv import closest_unbonded_atom_distances
from automol.geom._conv import could_be_forming_bond
from automol.geom._conv import ts_reacting_electron_direction
from automol.geom._conv import external_symmetry_factor
# # derived operations
from automol.geom._conv import apply_zmatrix_conversion
from automol.geom._conv import undo_zmatrix_conversion
from automol.geom._conv import set_distance
from automol.geom._conv import set_central_angle
from automol.geom._conv import set_dihedral_angle
# extra functions:
from automol.geom._extra import are_torsions_same
from automol.geom._extra import is_unique
from automol.geom._extra import hydrogen_bonded_structure
from automol.geom._extra import hydrogen_bonded_idxs
# align
from automol.geom._align import align
# ring functions
from automol.geom._ring import all_rings_angles_reasonable
from automol.geom._ring import ring_angles_reasonable
from automol.geom._ring import ring_fragments_geometry


__all__ = [
    # L2
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
    # L4
    # MolSym interface
    "point_group_from_geometry",
    # conversion functions:
    # # conversions
    'graph',
    'graph_without_stereo',
    'connectivity_graph_deprecated',
    'zmatrix',
    'zmatrix_with_conversion_info',
    'update_zmatrix',
    'amchi',
    'amchi_with_numbers',
    'inchi',
    'inchi_with_numbers',
    'chi',
    'chi_with_sort',
    'smiles',
    'rdkit_molecule',
    'py3dmol_view',
    'display',
    # # derived properties
    'is_connected',
    'linear_atoms',
    'closest_unbonded_atom_distances',
    'could_be_forming_bond',
    'ts_reacting_electron_direction',
    'external_symmetry_factor',
    # # derived operations
    'apply_zmatrix_conversion',
    'undo_zmatrix_conversion',
    'set_distance',
    'set_central_angle',
    'set_dihedral_angle',
    # extra functions:
    'are_torsions_same',
    'is_unique',
    'hydrogen_bonded_structure',
    'hydrogen_bonded_idxs',
    # align
    'align',
    # ring functions
    'all_rings_angles_reasonable',
    'ring_angles_reasonable',
    'ring_fragments_geometry',
]
