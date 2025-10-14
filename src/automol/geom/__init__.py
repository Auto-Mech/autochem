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
from .base._0core import from_data
from .base._0core import subgeom
# # getters
from .base._0core import symbols
from .base._0core import coordinates
# # setters
from .base._0core import set_coordinates
# # I/O
from .base._0core import string
from .base._0core import xyz_string
from .base._0core import xyz_trajectory_string
from .base._0core import from_string
from .base._0core import from_xyz_string
from .base._0core import from_xyz_string_with_mode
from .base._0core import xyz_string_comment
from .base._0core import from_xyz_trajectory_string
from .base._0core import parse_xyz_block
from .base._0core import yaml_data
from .base._0core import from_yaml_data
# # validation
from .base._0core import is_valid
# # conversions
from .base._0core import formula
from .base._0core import formula_string
# # properties
from .base._0core import count
from .base._0core import is_atom
from .base._0core import is_diatomic
from .base._0core import is_linear
from .base._0core import atom_count
from .base._0core import electron_count
from .base._0core import atom_indices
from .base._0core import dummy_atom_indices
from .base._0core import masses
from .base._0core import mass_weight_vector
from .base._0core import total_mass
from .base._0core import center_of_mass
from .base._0core import mass_centered
from .base._0core import aligned_to_principal_axes
from .base._0core import reduced_mass
from .base._0core import inertia_tensor
from .base._0core import principal_axes
from .base._0core import moments_of_inertia
from .base._0core import rotational_constants
from .base._0core import rotational_analysis
from .base._0core import translational_normal_modes
from .base._0core import rotational_normal_modes
# # geometric measurements
from .base._0core import distance
from .base._0core import central_angle
from .base._0core import dihedral_angle
from .base._0core import measure
from .base._0core import zmatrix_row_values
# # binary functions
from .base._0core import join
from .base._0core import join_sequence
from .base._0core import minimum_distance
from .base._0core import permutation
# # adding/removing atoms
from .base._0core import insert
from .base._0core import remove
from .base._0core import without_dummy_atoms
from .base._0core import reorder
from .base._0core import move_atom
from .base._0core import swap_coordinates
# # transformations
from .base._0core import round_
from .base._0core import displace
from .base._0core import translate
from .base._0core import translate_along_matrix
from .base._0core import perturb
from .base._0core import rotate
from .base._0core import transform
from .base._0core import transform_by_matrix
from .base._0core import reflect_coordinates
from .base._0core import shift_atom_position
from .base._0core import align_to_atoms
# comparison functions
# # properties used for comparisons
from .base._1comp import coulomb_spectrum
from .base._1comp import distance_matrix
# # comparisons
from .base._1comp import almost_equal
from .base._1comp import almost_equal_coulomb_spectrum
from .base._1comp import argunique_coulomb_spectrum
from .base._1comp import almost_equal_dist_matrix
from .base._1comp import minimum_volume_geometry
# intermolecular interactions
from .base._2intmol import has_low_relative_repulsion_energy
from .base._2intmol import total_repulsion_energy
from .base._2intmol import repulsion_energy
# L4
# MolSym interface
from ._0molsym import point_group_from_geometry
# conversion functions:
# # conversions
from ._1conv import graph
from ._1conv import graph_without_stereo
from ._1conv import connectivity_graph_deprecated
from ._1conv import zmatrix
from ._1conv import zmatrix_with_conversion_info
from ._1conv import update_zmatrix
from ._1conv import amchi
from ._1conv import amchi_with_numbers
from ._1conv import inchi
from ._1conv import inchi_with_numbers
from ._1conv import chi
from ._1conv import chi_with_sort
from ._1conv import smiles
from ._1conv import rdkit_molecule
from ._1conv import py3dmol_view
from ._1conv import display
from ._1conv import display_trajectory
# # derived properties
from ._1conv import is_connected
from ._1conv import linear_atoms
from ._1conv import closest_unbonded_atom_distances
from ._1conv import could_be_forming_bond
from ._1conv import ts_reacting_electron_direction
from ._1conv import external_symmetry_factor
# # derived operations
from ._1conv import apply_zmatrix_conversion
from ._1conv import undo_zmatrix_conversion
from ._1conv import set_distance
from ._1conv import set_central_angle
from ._1conv import set_dihedral_angle
# # interfaces
from ._1conv import ase_atoms
from ._1conv import from_ase_atoms
# vibrational analysis
from ._2vib import vibrational_analysis
# extra functions:
from ._extra import are_torsions_same
from ._extra import is_unique
from ._extra import hydrogen_bonded_structure
from ._extra import hydrogen_bonded_idxs
# align
from ._align import align
# ring functions
from ._ring import all_rings_angles_reasonable
from ._ring import ring_angles_reasonable
from ._ring import ring_fragments_geometry
from ._ring import ring_only_geometry
#adl cremer-pople paramters calculation
from ._ring import translate_to_ring_center
from ._ring import mean_ring_plane
from ._ring import normal_to_ring_plane
from ._ring import get_displacement
from ._ring import cremer_pople_params
from ._ring import checks_with_crest
from ._ring import dbscan

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
    'from_xyz_string_with_mode',
    'xyz_string_comment',
    'from_xyz_trajectory_string',
    'parse_xyz_block',
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
    'displace',
    'translate',
    'translate_along_matrix',
    'perturb',
    'rotate',
    'transform',
    'transform_by_matrix',
    'reflect_coordinates',
    'shift_atom_position',
    'align_to_atoms',
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
    # vibrational analysis
    'vibrational_analysis',
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
    'display_trajectory',
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
    # # interfaces
    'ase_atoms',
    'from_ase_atoms',
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
    'ring_only_geometry',
    #madl cremer-pople
    'translate_to_ring_center',
    'mean_ring_plane',
    'normal_to_ring_plane',
    'get_displacement',
    'cremer_pople_params',
    'checks_with_crest',
    'dbscan',
]
