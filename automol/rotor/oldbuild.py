"""
  Functions handling hindered rotor model calculations

    rotors = (
        rotor_dct1 = {
            rotor_tors_name1: {tors_inf: tors_bj},
            rotor_tors_name2: {tors_inf: tors_bj},
        },
        rotor_dct2 = {
            rotor_tors_name1: {tors_inf: tors_bj},
            rotor_tors_name2: {tors_inf: tors_bj},
        }
    )

"""

import itertools
from automol.rotor import _inf as inf
from automo.rotor.par import Torsion


# Build rotor and torsion objects
def from_zmatrix(zma, tors_names, scan_increment,
                 frm_bnd_keys=(), brk_bnd_keys=()),
                 model='1dhr':
    """ Build the rotor object consisting of all of the torsions
        and their info.

        :param zma: zmatrix
        :type zma:
        :param rotor_names: names of torsions in rotors
            [ ['D1'], ['D2', 'D3'], ['D4'] ]
        :type rotor_names: tuple(tuple(str))
    """

    rotor_dct = {}
    for tors_name in tors_names:
        rotor_dct[tors_name] = torsion(
            zma, tors_name, scan_increment,
            frm_bnd_keys=frm_bnd_keys, brk_bnd_keys=brk_bnd_keys
        )

    return rotor_dct


def from_geometry(geo, tors_names, scan_increment,
                  frm_bnd_keys=(), brk_bnd_keys=()):
    """ Build from a geometry
    """


def from_zmatrix(zma, tors_names, scan_increment,
                 frm_bnd_keys=(), brk_bnd_keys=()):


# Generate torsional names
def _geometry_torsional_names(geo,
                              frm_bnd_keys=(), brk_bnd_keys=()):
    """ build the torsional names
    """
    tors_names = automol.geom.zmatrix_torsion_coordinate_names(
        geo, ts_bnds=(frm_bnd_keys, brk_bnd_keys))


# from zmatrix __init__    
def _zmatrix_torsion_coordinate_names(zma):
    """ z-matrix torsional coordinate names
    (currently assumes torsional coordinates generated through x2z)
    """
    name_dct = standard_names(zma)
    inv_name_dct = dict(map(reversed, name_dct.items()))
    geo = automol.geom.without_dummy_atoms(geometry(zma))
    tors_names = automol.convert.geom.zmatrix_torsion_coordinate_names(geo)
    tors_names = tuple(map(inv_name_dct.__getitem__, tors_names))
    return tors_names


# Helper functions
def torsional_names(rotor_names):
    """ Build a flat list of all the torsions 
    """

