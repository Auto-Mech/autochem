""" cartesian geometry
"""
# constructors
from automol.geom._geom import from_subset
# getters
from automol.geom._geom import symbols
from automol.geom._geom import coordinates
from automol.geom._geom import count
from automol.geom._geom import atom_count
from automol.geom._geom import atom_indices
from automol.geom._geom import dummy_atom_indices
from automol.geom._geom import components_graph
from automol.geom._geom import connected
# validation
from automol.geom._geom import is_valid
# setters
from automol.geom._geom import set_coordinates
from automol.geom._geom import without_dummy_atoms
# I/O
from automol.geom._geom import from_string
from automol.geom._geom import from_xyz_string
from automol.geom._geom import string
from automol.geom._geom import xyz_string
from automol.geom._geom import xyz_trajectory_string
from automol.geom._geom import from_xyz_trajectory_string
# representations
from automol.geom._geom import coulomb_spectrum
# comparisons
from automol.geom._geom import almost_equal
from automol.geom._geom import minimum_distance
from automol.geom._geom import almost_equal_coulomb_spectrum
from automol.geom._geom import argunique_coulomb_spectrum
# operations
from automol.geom._geom import remove
from automol.geom._geom import join
from automol.geom._geom import reorder
from automol.geom._geom import swap_coordinates
from automol.geom._geom import insert
from automol.geom._geom import insert_dummies_on_linear_atoms
from automol.geom._geom import displace
from automol.geom._geom import translate
from automol.geom._geom import invert
from automol.geom._geom import transform
from automol.geom._geom import transform_by_matrix
from automol.geom._geom import rotate
from automol.geom._geom import euler_rotate
from automol.geom._geom import move_coordinates
from automol.geom._geom import reflect_coordinates
# geometric properties
from automol.geom._geom import distance
from automol.geom._geom import central_angle
from automol.geom._geom import dihedral_angle
from automol.geom._geom import zmatrix_row_values
from automol.geom._geom import linear_atoms
from automol.geom._geom import distance_matrix
from automol.geom._geom import closest_unbonded_atoms
# geometric properties
from automol.geom._geom import is_atom
from automol.geom._geom import masses
from automol.geom._geom import total_mass
from automol.geom._geom import center_of_mass
from automol.geom._geom import reduced_mass
from automol.geom._geom import mass_centered
from automol.geom._geom import inertia_tensor
from automol.geom._geom import principal_axes
from automol.geom._geom import moments_of_inertia
from automol.geom._geom import rotational_constants
from automol.geom._geom import is_linear
from automol.geom._geom import permutation
# extras
from automol.geom._extra import end_group_sym_factor
from automol.geom._extra import rot_permutated_geoms
from automol.geom._extra import almost_equal_dist_matrix
from automol.geom._extra import find_xyzp_using_internals

# constructors
import automol.create.geom
# conversions
import automol.convert.geom

# submodules
from automol.geom import ts


def from_data(syms, xyzs, angstrom=False):
    """ geometry data structure from symbols and coordinates
    """
    return automol.create.geom.from_data(
        symbols=syms, coordinates=xyzs, angstrom=angstrom)


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


def graph(geo, remove_stereo=False):
    """ geometry => graph
    """
    return automol.convert.geom.graph(
        geo, remove_stereo=remove_stereo)


def connectivity_graph(geo, dummy_bonds=True,
                       rqq_bond_max=3.5, rqh_bond_max=2.6, rhh_bond_max=1.9):
    """ geometry => connectivity graph
    """
    gra = automol.convert.geom.connectivity_graph(
        geo, dummy_bonds=dummy_bonds, rqq_bond_max=rqq_bond_max,
        rqh_bond_max=rqh_bond_max, rhh_bond_max=rhh_bond_max)
    return gra


def inchi(geo, remove_stereo=False):
    """ geometry => inchi
    """
    return automol.convert.geom.inchi(geo, remove_stereo=remove_stereo)


def smiles(geo, remove_stereo=False):
    """ geometry => inchi
    """
    ich = inchi(geo, remove_stereo=remove_stereo)
    return automol.convert.inchi.smiles(ich)


def formula(geo):
    """ geometry => formula
    """
    return automol.convert.geom.formula(geo)


__all__ = [
    # constructors
    'from_subset',
    # getters
    'symbols',
    'coordinates',
    'count',
    'atom_count',
    'atom_indices',
    'dummy_atom_indices',
    'components_graph',
    'connected',
    # validation
    'is_valid',
    # setters
    'set_coordinates',
    'without_dummy_atoms',
    # I/O
    'from_string',
    'from_xyz_string',
    'string',
    'xyz_string',
    'xyz_trajectory_string',
    'from_xyz_trajectory_string',
    # representations
    'coulomb_spectrum',
    # comparisons
    'almost_equal',
    'minimum_distance',
    'almost_equal_coulomb_spectrum',
    'argunique_coulomb_spectrum',
    # operations
    'remove',
    'join',
    'reorder',
    'swap_coordinates',
    'insert',
    'insert_dummies_on_linear_atoms',
    'displace',
    'translate',
    'invert',
    'transform',
    'transform_by_matrix',
    'rotate',
    'euler_rotate',
    'move_coordinates',
    'reflect_coordinates',
    # geometric properties
    'distance',
    'central_angle',
    'dihedral_angle',
    'zmatrix_row_values',
    'linear_atoms',
    'distance_matrix',
    'closest_unbonded_atoms',
    # geometric properties
    'is_atom',
    'masses',
    'total_mass',
    'center_of_mass',
    'reduced_mass',
    'mass_centered',
    'inertia_tensor',
    'principal_axes',
    'moments_of_inertia',
    'rotational_constants',
    'is_linear',
    'permutation',
    # extras
    'end_group_sym_factor',
    'rot_permutated_geoms',
    'almost_equal_dist_matrix',
    'find_xyzp_using_internals',
    # submodules
    'ts',
]
