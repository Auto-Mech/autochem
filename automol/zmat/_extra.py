""" extra Level 4 Z-Matrix functions
"""

import numpy
import automol.graph
import automol.geom
from automol.zmat._conv import graph
from automol.zmat._conv import geometry
from automol.zmat.base import coordinates
from automol.zmat.base import dihedral_angle_names


def torsional_symmetry_numbers(zma, tors_names,
                               frm_bnd_key=None, brk_bnd_key=None):
    """ symmetry numbers for torsional dihedrals
    """
    dih_edg_key_dct = _dihedral_edge_keys(zma)
    assert set(tors_names) <= set(dih_edg_key_dct.keys())
    edg_keys = tuple(map(dih_edg_key_dct.__getitem__, tors_names))

    gra = graph(zma, stereo=False)
    bnd_sym_num_dct = automol.graph.bond_symmetry_numbers(
        gra, frm_bnd_key, brk_bnd_key)
    tors_sym_nums = []
    for edg_key in edg_keys:
        if edg_key in bnd_sym_num_dct.keys():
            sym_num = bnd_sym_num_dct[edg_key]
        else:
            sym_num = 1
        tors_sym_nums.append(sym_num)

    tors_sym_nums = tuple(tors_sym_nums)
    return tors_sym_nums


def torsional_scan_linspaces(zma, tors_names, increment=0.5,
                             frm_bnd_key=None, brk_bnd_key=None):
    """ scan grids for torsional dihedrals
    """
    sym_nums = torsional_symmetry_numbers(
        zma, tors_names, frm_bnd_key=frm_bnd_key, brk_bnd_key=brk_bnd_key)
    intervals = tuple(2*numpy.pi/sym_num - increment for sym_num in sym_nums)
    npoints_lst = tuple(
        (int(interval / increment)+1) for interval in intervals)
    return tuple((0, interval, npoints)
                 for interval, npoints in zip(intervals, npoints_lst))


def is_atom_closest_to_bond_atom(zma, idx_rad, bond_dist):
    """ Check to see whether the radical atom is still closest to the bond
        formation site.
    """
    geo = geometry(zma)
    atom_closest = True
    for idx, _ in enumerate(geo):
        if idx < idx_rad:
            dist = automol.geom.distance(geo, idx, idx_rad)
            if dist < bond_dist-0.01:
                atom_closest = False
    return atom_closest


# helpers
def _dihedral_edge_keys(zma):
    """ dihedral bonds, by name
    """
    coo_dct = coordinates(zma)
    dih_names = dihedral_angle_names(zma)
    dih_keys_lst = tuple(map(coo_dct.__getitem__, dih_names))
    dih_edg_key_dct = {dih_name: frozenset(dih_key[1:3])
                       for dih_name, dih_keys in zip(dih_names, dih_keys_lst)
                       for dih_key in dih_keys}
    return dih_edg_key_dct
