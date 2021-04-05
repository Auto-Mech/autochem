""" z-matrix
"""
# constructors
from automol.zmat._zmat import from_geometry
# converters
from automol.zmat._zmat import geometry
# getters
from automol.zmat._zmat import count
from automol.zmat._zmat import symbols
from automol.zmat._zmat import atom_indices
from automol.zmat._zmat import key_matrix
from automol.zmat._zmat import name_matrix
from automol.zmat._zmat import value_matrix
from automol.zmat._zmat import value_dictionary
from automol.zmat._zmat import dummy_keys
from automol.zmat._zmat import dummy_neighbor_keys
from automol.zmat._zmat import linear_atom_keys
from automol.zmat._zmat import distance
from automol.zmat._zmat import central_angle
from automol.zmat._zmat import dihedral_angle
# setters
from automol.zmat._zmat import set_key_matrix
from automol.zmat._zmat import set_name_matrix
from automol.zmat._zmat import set_value_matrix
from automol.zmat._zmat import set_values_by_name
from automol.zmat._zmat import standard_name_matrix
from automol.zmat._zmat import standard_form
from automol.zmat._zmat import rename
from automol.vmat import standard_names
# operations
from automol.zmat._zmat import add_atom
from automol.zmat._zmat import remove_atom
from automol.zmat._zmat import join_replace_one
# comparisons
from automol.zmat._zmat import almost_equal
# I/O
from automol.zmat._zmat import from_string
from automol.zmat._zmat import string
# validators
from automol.zmat._zmat import is_valid
# extra names we need to change
from automol.zmat._zmat import set_constraint_names

# extra things to be added into the module
from automol.zmat._extra import torsion_coordinate_name
from automol.zmat._extra import torsion_leading_atom
from automol.zmat._extra import distance_coordinate_name
from automol.zmat._extra import central_angle_coordinate_name
from automol.zmat._extra import dihedral_angle_coordinate_name

# constructors
import automol.create.zmat
# conversions
import automol.convert.zmat

# to be moved to permanent file upon cleanup
from automol.zmat._new import is_atom_closest_to_bond_atom

# to be depreciated
from automol.zmat._zmat import samples
from automol.zmat._zmat import torsional_symmetry_numbers
from automol.zmat._zmat import torsional_sampling_ranges
from automol.zmat._zmat import torsional_scan_linspaces
from automol.zmat._zmat import dummy_coordinate_names
from automol.zmat._zmat import coordinates
from automol.zmat._zmat import dihedral_angle_names
from automol.zmat._zmat import dihedral_axis_name
from automol.zmat._zmat import coord_idxs
from automol.zmat._zmat import bond_key_from_idxs
from automol.zmat._zmat import get_babs1
from automol.zmat._zmat import get_babs2
from automol.zmat._zmat import constraint_dct
from automol.zmat._zmat import distance_names



def from_data(syms, key_mat, val_mat, name_mat=None,
              one_indexed=False, angstrom=False, degree=False):
    """ z-matrix constructor

    :param syms: atomic symbols
    :type syms: tuple[str]
    :param key_mat: key/index columns of the z-matrix, zero-indexed
    :type key_mat: tuple[tuple[int, int or None, int or None]]
    :param val_mat: z-matrix coordinate values
    :type val_mat: tuple[tuple[float, float or None, float or None]]
    :param name_mat: coordinate name columns of the z-matrix
    :type name_mat; tuple[tuple[str, str or None, str or None]]
    """
    return automol.create.zmat.from_data(
        symbols=syms, key_matrix=key_mat, value_matrix=val_mat,
        name_matrix=name_mat, one_indexed=one_indexed, angstrom=angstrom,
        degree=degree)


def torsion_coordinate_names(zma):
    """ z-matrix torsional coordinate names

    (currently assumes torsional coordinates generated through x2z)
    """
    name_dct = standard_names(zma)
    inv_name_dct = dict(map(reversed, name_dct.items()))
    geo = automol.geom.without_dummy_atoms(geometry(zma))
    tors_names = automol.convert.geom.zmatrix_torsion_coordinate_names(geo)
    tors_names = tuple(map(inv_name_dct.__getitem__, tors_names))
    return tors_names


def graph(zma, stereo=True, dummy=False):
    """ z-matrix => graph
    """
    return automol.convert.zmat.graph(zma, stereo=stereo, dummy=dummy)


def connectivity_graph(zma,
                       rqq_bond_max=3.5, rqh_bond_max=2.6, rhh_bond_max=1.9):
    """ z-matrix => connectivity graph
    """
    gra = automol.convert.zmat.connectivity_graph(
        zma, rqq_bond_max=rqq_bond_max, rqh_bond_max=rqh_bond_max,
        rhh_bond_max=rhh_bond_max)
    return gra


def formula(zma):
    """ zmatrix => formula
    """
    return automol.convert.zmat.formula(zma)


__all__ = [
    # constructors
    'from_geometry',
    # converters
    'geometry',
    # getters
    'count',
    'symbols',
    'key_matrix',
    'name_matrix',
    'value_matrix',
    'value_dictionary',
    'dummy_keys',
    'dummy_neighbor_keys',
    'linear_atom_keys',
    'distance',
    'central_angle',
    'dihedral_angle',
    # setters
    'set_key_matrix',
    'set_name_matrix',
    'set_value_matrix',
    'set_values_by_name',
    'standard_name_matrix',
    'rename',
    'standard_form',
    'standard_names',
    # operations
    'add_atom',
    'remove_atom',
    'join_replace_one',
    # comparisons
    'almost_equal',
    # I/O
    'from_string',
    'string',
    # validator
    'is_valid',
    # extra stuff 
    'set_constraint_names',

    # extras
    'torsion_coordinate_name',
    'torsion_leading_atom',
    'distance_coordinate_name',
    'central_angle_coordinate_name',
    'dihedral_angle_coordinate_name',

    # constructors
    'from_data',
    # conversions,
    'geometry',
    'graph',
    'formula',

    # to be moved to permanent file upon cleanup
    'is_atom_closest_to_bond_atom',

    # to be depreciated
    'samples',
    'torsional_symmetry_numbers',
    'torsional_sampling_ranges',
    'torsional_scan_linspaces',
    'dummy_coordinate_names',
    'coordinates',
    'dihedral_angle_names',
    'dihedral_axis_name',
    'coord_idxs',
    'bond_key_from_idxs',
    'get_babs1',
    'get_babs2',
    'constraint_dct',
    'distance_names',
]

