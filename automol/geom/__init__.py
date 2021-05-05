""" cartesian geometry
"""

# Base geometry functions
# constructors
from automol.geom._base import from_subset
# getters
from automol.geom._base import symbols
from automol.geom._base import coordinates
from automol.geom._base import count
from automol.geom._base import atom_count
from automol.geom._base import atom_indices
from automol.geom._base import dummy_atom_indices
from automol.geom._base import components_graph
from automol.geom._base import connected
# validation
from automol.geom._base import is_valid
# setters
from automol.geom._base import set_coordinates
from automol.geom._base import without_dummy_atoms
# I/O
from automol.geom._base import from_string
from automol.geom._base import from_xyz_string
from automol.geom._base import string
from automol.geom._base import xyz_string
from automol.geom._base import xyz_trajectory_string
from automol.geom._base import from_xyz_trajectory_string
# geometric properties
from automol.geom._base import distance
from automol.geom._base import central_angle
from automol.geom._base import dihedral_angle
from automol.geom._base import zmatrix_row_values
from automol.geom._base import linear_atoms
from automol.geom._base import closest_unbonded_atoms
# chemical properties
from automol.geom._prop import is_atom
from automol.geom._prop import masses
from automol.geom._prop import total_mass
from automol.geom._prop import center_of_mass
from automol.geom._prop import reduced_mass
from automol.geom._prop import inertia_tensor
from automol.geom._prop import principal_axes
from automol.geom._prop import moments_of_inertia
from automol.geom._prop import rotational_constants
from automol.geom._prop import is_linear
from automol.geom._prop import permutation
# transformations and operations
from automol.geom._trans import remove_coordinates
from automol.geom._trans import join
from automol.geom._trans import reorder_coordinates
from automol.geom._trans import swap_coordinates
from automol.geom._trans import insert
from automol.geom._trans import insert_dummies_on_linear_atoms
from automol.geom._trans import displace
from automol.geom._trans import translate
from automol.geom._trans import transform
from automol.geom._trans import transform_by_matrix
from automol.geom._trans import rotate
from automol.geom._trans import euler_rotate
from automol.geom._trans import reflect_coordinates
from automol.geom._trans import mass_centered
# comparisons
from automol.geom._comp import distance_matrix
from automol.geom._comp import coulomb_spectrum
from automol.geom._comp import almost_equal_dist_matrix
from automol.geom._comp import almost_equal
from automol.geom._comp import minimum_distance
from automol.geom._comp import almost_equal_coulomb_spectrum
from automol.geom._comp import argunique_coulomb_spectrum
from automol.geom._comp import is_unique
from automol.geom._comp import are_torsions_same
from automol.geom._comp import are_torsions_same2
# extras
from automol.geom._extra import end_group_symmetry_factor
from automol.geom._extra import rot_permutated_geoms

# constructors
import automol.create.geom
# conversions
import automol.convert.geom

# submodules
from automol.geom import ts


def from_data(symbs, xyzs, angstrom=False):
    """ geometry data structure from symbols and coordinates
    """
    return automol.create.geom.from_data(
        symbols=symbs, coordinates=xyzs, angstrom=angstrom)


# conversions
def zmatrix(geo, ts_bnds=()):
    """ geometry => z-matrix
    """
    zma, _, _ = automol.convert.geom.zmatrix(geo, ts_bnds)
    return zma


def zmatrix_torsion_coordinate_names(geo, ts_bnds=()):
    """ z-matrix torsional coordinate names
    """
    return automol.convert.geom.zmatrix_torsion_coordinate_names(geo, ts_bnds)


def zmatrix_atom_ordering(geo, ts_bnds=()):
    """ z-matrix atom ordering
    """
    return automol.convert.geom.zmatrix_atom_ordering(geo, ts_bnds)


def external_symmetry_factor(geo):
    """ obtain external symmetry factor for a geometry using x2z
    """
    return automol.convert.geom.external_symmetry_factor(geo)


def graph(geo, stereo=True):
    """ geometry => graph
    """
    return automol.convert.geom.graph(
        geo, stereo=stereo)


def connectivity_graph(geo,
                       rqq_bond_max=3.5, rqh_bond_max=2.6, rhh_bond_max=1.9):
    """ geometry => connectivity graph
    """
    gra = automol.convert.geom.connectivity_graph(
        geo, rqq_bond_max=rqq_bond_max,
        rqh_bond_max=rqh_bond_max, rhh_bond_max=rhh_bond_max)
    return gra


def inchi(geo, stereo=True):
    """ geometry => inchi
    """
    return automol.convert.geom.inchi(geo, stereo=stereo)


def smiles(geo, stereo=True):
    """ geometry => inchi
    """
    ich = inchi(geo, stereo=stereo)
    return automol.convert.inchi.smiles(ich)


def formula(geo):
    """ geometry => formula
    """
    return automol.convert.geom.formula(geo)


def formula_string(geo):
    """ geometry => formula_string
    """
    return automol.convert.geom.formula_string(geo)


__all__ = [
    # Base geometry functions
    'from_subset',
    'symbols',
    'coordinates',
    'count',
    'atom_count',
    'atom_indices',
    'dummy_atom_indices',
    'components_graph',
    'connected',
    'is_valid',
    'set_coordinates',
    'without_dummy_atoms',
    'from_string',
    'from_xyz_string',
    'string',
    'xyz_string',
    'xyz_trajectory_string',
    'from_xyz_trajectory_string',
    'distance',
    'central_angle',
    'dihedral_angle',
    'zmatrix_row_values',
    'linear_atoms',
    'closest_unbonded_atoms',
    # chemical properties
    'is_atom',
    'masses',
    'total_mass',
    'center_of_mass',
    'reduced_mass',
    'inertia_tensor',
    'principal_axes',
    'moments_of_inertia',
    'rotational_constants',
    'is_linear',
    'permutation',
    # transformations and operations
    'remove_coordinates',
    'join',
    'reorder_coordinates',
    'swap_coordinates',
    'insert',
    'insert_dummies_on_linear_atoms',
    'displace',
    'translate',
    'transform',
    'transform_by_matrix',
    'rotate',
    'euler_rotate',
    'reflect_coordinates',
    'mass_centered',
    # comparisons
    'distance_matrix',
    'coulomb_spectrum',
    'almost_equal',
    'minimum_distance',
    'almost_equal_coulomb_spectrum',
    'argunique_coulomb_spectrum',
    'almost_equal_dist_matrix',
    'is_unique',
    'are_torsions_same',
    'are_torsions_same2',
    # extras
    'end_group_symmetry_factor',
    'rot_permutated_geoms',
    # submodules
    'ts'
]
