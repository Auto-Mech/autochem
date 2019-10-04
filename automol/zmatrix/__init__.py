""" z-matrix
"""
# getters
from automol.zmatrix._zmatrix import var_
from automol.zmatrix._zmatrix import count
from automol.zmatrix._zmatrix import symbols
from automol.zmatrix._zmatrix import key_matrix
from automol.zmatrix._zmatrix import name_matrix
from automol.zmatrix._zmatrix import value_matrix
from automol.zmatrix._zmatrix import coordinate_key_matrix
from automol.zmatrix._zmatrix import coordinates
from automol.zmatrix._zmatrix import bond_idxs
from automol.zmatrix._zmatrix import get_babs1
from automol.zmatrix._zmatrix import get_babs2
from automol.zmatrix._zmatrix import names
from automol.zmatrix._zmatrix import distance_names
from automol.zmatrix._zmatrix import central_angle_names
from automol.zmatrix._zmatrix import dihedral_angle_names
from automol.zmatrix._zmatrix import angle_names
from automol.zmatrix._zmatrix import dummy_coordinate_names
from automol.zmatrix._zmatrix import values
# validation
from automol.zmatrix._zmatrix import is_valid
# setters
from automol.zmatrix._zmatrix import set_keys
from automol.zmatrix._zmatrix import set_names
from automol.zmatrix._zmatrix import set_values
from automol.zmatrix._zmatrix import shift_row_to_end
from automol.zmatrix._zmatrix import standard_names
from automol.zmatrix._zmatrix import standard_form
# operations
from automol.zmatrix._zmatrix import join
from automol.zmatrix._zmatrix import insert_dummy_atom
# misc
from automol.zmatrix._zmatrix import is_standard_form
# I/O
from automol.zmatrix._zmatrix import from_string
from automol.zmatrix._zmatrix import string
# comparisons
from automol.zmatrix._zmatrix import almost_equal
# random sampling
from automol.zmatrix._zmatrix import samples
# z-matrix torsional degrees of freedom
from automol.zmatrix._zmatrix import torsional_symmetry_numbers
from automol.zmatrix._zmatrix import torsional_sampling_ranges
from automol.zmatrix._zmatrix import torsional_scan_linspaces

# submodules
from automol.zmatrix import ts

# constructors
import automol.create.zmatrix
# conversions
import automol.convert.zmatrix


def from_data(syms, key_mat, name_mat, val_dct,
              one_indexed=False, angstrom=False, degree=False):
    """ z-matrix constructor

    :param syms: atomic symbols
    :type syms: tuple[str]
    :param key_mat: key/index columns of the z-matrix, zero-indexed
    :type key_mat: tuple[tuple[float, float or None, float or None]]
    :param name_mat: coordinate name columns of the z-matrix
    :type name_mat; tuple[tuple[str, str or None, str or None]]
    :param val_dct: coordinate values, by coordinate name
    :type val_dct: dict
    """
    return automol.create.zmatrix.from_data(
        symbols=syms, key_matrix=key_mat, name_matrix=name_mat, values=val_dct,
        one_indexed=one_indexed, angstrom=angstrom, degree=degree)


def geometry(zma, remove_dummy_atoms=None):
    """ z-matrix => geometry
    """
    return automol.convert.zmatrix.geometry(
        zma, remove_dummy_atoms=remove_dummy_atoms)


def torsion_coordinate_names(zma):
    """ z-matrix torsional coordinate names

    (currently assumes torsional coordinates generated through x2z)
    """
    name_dct = standard_names(zma)
    inv_name_dct = dict(map(reversed, name_dct.items()))
    geo = automol.geom.without_dummy_atoms(geometry(zma))
    # line below should always break if zma variable has a dummy atom
    # assert var_(standard_form(zma)) == var_(automol.convert.geom.zmatrix(geo))
    tors_names = automol.convert.geom.zmatrix_torsion_coordinate_names(geo)
    tors_names = tuple(map(inv_name_dct.__getitem__, tors_names))
    return tors_names


def graph(zma, remove_stereo=False):
    """ z-matrix => graph
    """
    return automol.convert.zmatrix.graph(zma, remove_stereo=remove_stereo)


def formula(zma):
    """ zmatrix => formula
    """
    return automol.convert.zmatrix.formula(zma)


__all__ = [
    # constructors
    'from_data',

    # getters
    'var_',
    'count',
    'symbols',
    'key_matrix',
    'name_matrix',
    'value_matrix',
    'coordinate_key_matrix',
    'coordinates',
    'bond_idxs',
    'get_babs1',
    'get_babs2',
    'names',
    'distance_names',
    'central_angle_names',
    'dihedral_angle_names',
    'angle_names',
    'dummy_coordinate_names',
    'torsion_coordinate_names',
    'values',
    # validation
    'is_valid',
    # setters
    'set_keys',
    'set_names',
    'set_values',
    'shift_row_to_end',
    'standard_names',
    'standard_form',
    # operations
    'join',
    'insert_dummy_atom',
    # misc
    'is_standard_form',
    # I/O
    'from_string',
    'string',
    # comparisons
    'almost_equal',
    # random sampling
    'samples',
    # z-matrix torsional degrees of freedom
    'torsional_symmetry_numbers',
    'torsional_sampling_ranges',
    'torsional_scan_linspaces',

    # submodules
    'ts',

    # conversions,
    'geometry',
    'graph',
    'formula',
]
