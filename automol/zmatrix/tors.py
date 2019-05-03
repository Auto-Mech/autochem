""" functions for working with torsional degrees of freedom
"""
import numpy
from automol.graph import bond_symmetry_numbers as _bond_symmetry_numbers
from automol.zmatrix._graph import connectivity_graph as _connectivity_graph
from automol.zmatrix._core import set_values as _set_values
from automol.zmatrix._core import coordinates as _coordinates
from automol.zmatrix._core import (dihedral_angle_names as
                                   _dihedral_angle_names)


def symmetry_numbers(zma, tors_names):
    """ symmetry numbers for torsional dihedrals
    """
    dih_edg_key_dct = _dihedral_edge_keys(zma)
    assert set(tors_names) <= set(dih_edg_key_dct.keys())
    edg_keys = tuple(map(dih_edg_key_dct.__getitem__, tors_names))

    gra = _connectivity_graph(zma)
    bnd_sym_num_dct = _bond_symmetry_numbers(gra)
    assert set(edg_keys) <= set(bnd_sym_num_dct.keys())

    tors_sym_nums = tuple(map(bnd_sym_num_dct.__getitem__, edg_keys))
    return tors_sym_nums


def sampling_ranges(zma, tors_names):
    """ sampling ranges for torsional dihedrals
    """
    sym_nums = symmetry_numbers(zma, tors_names)
    return tuple((0, 2*numpy.pi/sym_num) for sym_num in sym_nums)


def samples(zma, nsamp, tors_range_dct):
    """ randomly sample over torsional dihedrals
    """
    tors_names = tuple(tors_range_dct.keys())
    tors_ranges = tuple(tors_range_dct.values())
    tors_vals_lst = _sample_over_ranges(tors_ranges, nsamp)

    zmas = tuple(_set_values(zma, dict(zip(tors_names, tors_vals)))
                 for tors_vals in tors_vals_lst)
    return zmas


def _dihedral_edge_keys(zma):
    """ dihedral bonds, by name
    """
    coo_dct = _coordinates(zma)
    dih_names = _dihedral_angle_names(zma)
    dih_keys_lst = tuple(map(coo_dct.__getitem__, dih_names))
    dih_edg_key_dct = {dih_name: frozenset(dih_key[1:3])
                       for dih_name, dih_keys in zip(dih_names, dih_keys_lst)
                       for dih_key in dih_keys}
    return dih_edg_key_dct


def _sample_over_ranges(rngs, nsamp):
    """ randomly sample over several ranges
    """
    nrng = len(rngs)
    samp_mat = numpy.random.rand(nsamp, nrng)
    for i, (start, stop) in enumerate(rngs):
        samp_mat[:, i] = samp_mat[:, i] * (stop - start) + start
    return tuple(map(tuple, samp_mat))
