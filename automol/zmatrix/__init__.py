""" z-matrix
"""
# constructors
from automol.zmatrix._zmatrix import from_geometry
# getters
from automol.zmatrix._zmatrix import var_
from automol.zmatrix._zmatrix import count
from automol.zmatrix._zmatrix import atom_count
from automol.zmatrix._zmatrix import symbols
from automol.zmatrix._zmatrix import key_matrix
from automol.zmatrix._zmatrix import name_matrix
from automol.zmatrix._zmatrix import value_matrix
from automol.zmatrix._zmatrix import coordinate_key_matrix
from automol.zmatrix._zmatrix import coordinates
from automol.zmatrix._zmatrix import coord_idxs
from automol.zmatrix._zmatrix import bond_key_from_idxs
from automol.zmatrix._zmatrix import get_babs1
from automol.zmatrix._zmatrix import get_babs2
from automol.zmatrix._zmatrix import names
from automol.zmatrix._zmatrix import distance_names
from automol.zmatrix._zmatrix import central_angle_names
from automol.zmatrix._zmatrix import dihedral_angle_names
from automol.zmatrix._zmatrix import new_distance_name
from automol.zmatrix._zmatrix import new_central_angle_name
from automol.zmatrix._zmatrix import new_dihedral_angle_name
from automol.zmatrix._zmatrix import angle_names
from automol.zmatrix._zmatrix import dummy_coordinate_names
from automol.zmatrix._zmatrix import atom_indices
from automol.zmatrix._zmatrix import atom_specifiers
from automol.zmatrix._zmatrix import atom_dependents
from automol.zmatrix._zmatrix import dummy_atom_anchors
from automol.zmatrix._zmatrix import values
from automol.zmatrix._zmatrix import _dihedral_edge_keys as dihedral_edge_keys

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
from automol.zmatrix._zmatrix import append
from automol.zmatrix._zmatrix import join
from automol.zmatrix._zmatrix import insert_dummy_atom
from automol.zmatrix._zmatrix import convert
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
from automol.zmatrix._util import shifted_standard_zmas_graphs

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
    tors_names = automol.convert.geom.zmatrix_torsion_coordinate_names(geo)
    tors_names = tuple(map(inv_name_dct.__getitem__, tors_names))
    return tors_names


def graph(zma, remove_stereo=False):
    """ z-matrix => graph
    """
    return automol.convert.zmatrix.graph(zma, remove_stereo=remove_stereo)


def connectivity_graph(zma,
                       rqq_bond_max=3.5, rqh_bond_max=2.6, rhh_bond_max=1.9):
    """ z-matrix => connectivity graph
    """
    gra = automol.convert.zmatrix.connectivity_graph(
        zma, rqq_bond_max=rqq_bond_max, rqh_bond_max=rqh_bond_max,
        rhh_bond_max=rhh_bond_max)
    return gra

def chain_between(zma, atm1, atm2):
    """ returns the chain of atoms connecting two given atoms
    """              
    chain1 = _chain_between(zma, atm1, atm2)
    chain2 = _chain_between(zma, atm2, atm1)
    if len(chain2) > 0 and len(chain2) < len(chain1):
        chain1 = chain2
    if len(chain1) < 1:
        chain1 = chain2
    return chain1

def _chain_between(zma, atm1, atm2):
    """ returns the chain of atoms connecting two given atoms
    """              
    def _loop(atm_nbh_keys_dct, nbh_atms, chain, atm2, iters, checked_atms):
        if iters < 800:
            new_atm_nbhs = []
            for atm in nbh_atms:
                if atm not in checked_atms:
                    new_atm_nbhs.append(atm)
            for nbh in new_atm_nbhs: 
                chain.append(nbh)
                #print('neigh atom', nbh)
                #print('chain', chain)
                #print('iter', iters)
                checked_atms.append(nbh)
                iters += 1
                if nbh == atm2:           
                    iters = 1000
                    return chain, nbh_atms, iters, checked_atms
                else:
                    nbh_atms = atm_nbh_keys_dct[nbh]
                    chain, atm_nbhs, iters, checked_atms = _loop(atm_nbh_keys_dct, nbh_atms, chain, atm2, iters, checked_atms)
                    if iters > 900:
                        break
                chain = chain[:-1]
                nbh_atms= atm_nbh_keys_dct[chain[-1]]
        return chain, nbh_atms, iters, checked_atms

    cgra = connectivity_graph(zma)
    atm_nbh_keys_dct = automol.graph.atom_neighbor_keys(cgra)
    nbh_atms = atm_nbh_keys_dct[atm1]
    chain = [atm1]
    checked_atms = [atm1]
    iters = 0
    return_chain, nbh_atms, iters, checked_atms =_loop(atm_nbh_keys_dct, nbh_atms, chain, atm2, iters, checked_atms)
    if iters < 900:
        return_chain = []
    return return_chain    


def formula(zma):
    """ zmatrix => formula
    """
    return automol.convert.zmatrix.formula(zma)


__all__ = [
    # constructors
    'from_data',
    'from_geometry',

    # getters
    'var_',
    'count',
    'atom_count',
    'symbols',
    'key_matrix',
    'name_matrix',
    'value_matrix',
    'coordinate_key_matrix',
    'coordinates',
    'coord_idxs',
    'bond_key_from_idxs',
    'get_babs1',
    'get_babs2',
    'names',
    'distance_names',
    'central_angle_names',
    'dihedral_angle_names',
    'new_distance_name',
    'new_central_angle_name',
    'new_dihedral_angle_name',
    'angle_names',
    'dummy_coordinate_names',
    'atom_indices',
    'atom_specifiers',
    'atom_dependents',
    'dummy_atom_anchors',
    'torsion_coordinate_names',
    'values',
    'dihedral_edge_keys',
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
    'append',
    'join',
    'insert_dummy_atom',
    'convert',
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
    'shifted_standard_zmas_graphs',

    # conversions,
    'geometry',
    'graph',
    'formula',
]
