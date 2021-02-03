"""
 Handle rotor objects

 Rotors: (Rotor1, Rotor2, ..., RotorN)
 Rotor: (tors_obj_1, tors_obj_2, ..., tors_obj_N)
"""

from itertools import chain
import yaml
import numpy
import automol.zmat
import automol.pot
from automol.rotor import _tors as tors
from automol.rotor._name import group_torsions_into_rotors
from automol.rotor._util import graph_with_keys
from automol.rotor._util import sort_tors_names


# constructors
def from_zma(zma, tors_names=None, multi=False):
    """ Construct a list-of-lists of torsion objects
    """

    # Build a graph that is used to get torsion object info
    gra, lin_keys = graph_with_keys(zma)

    # Build the torsion objects
    tors_lst = tors.torsion_lst(zma, gra, lin_keys)

    # Place the torsions into order based on rotors
    rotors = group_torsions_into_rotors(
        tors_lst, name_grps=tors_names, multi=multi)

    return rotors


def from_data(zma, tors_inf_dct, tors_names=None):
    """ Build the rotors objects from existing data
    """

    tors_lst = ()
    for name, dct in tors_inf_dct.items():
        assert set(dct.keys()) >= {'symmetry', 'axis', 'groups'}, (
            'must have symmetry, axis, and groups in dct tp build torsions'
        )
        tors_lst += (tors.Torsion(zma, name, **dct),)

    rotors = group_torsions_into_rotors(tors_lst, name_grps=tors_names)

    return rotors


# Getters
def names(rotor_lst, flat=False):
    """ Get a flat list of names for list of rotors
    """
    _names = tuple(tuple(torsion.name for torsion in rotor)
                   for rotor in rotor_lst)
    if flat:
        _names = tuple(chain(*_names))
    return _names


def axes(rotor_lst):
    """ Build a list of list of axes(
    """
    return tuple(tuple(torsion.axis for torsion in rotor)
                 for rotor in rotor_lst)


def groups(rotor_lst):
    """ Build a list of list of axes(
    """
    return tuple(tuple(torsion.groups for torsion in rotor)
                 for rotor in rotor_lst)


def symmetries(rotor_lst):
    """ Build a list of list of axes(
    """
    return tuple(tuple(torsion.symmetry for torsion in rotor)
                 for rotor in rotor_lst)


def grids(rotor_lst, span=2.0*numpy.pi, increment=0.523599, flat=False):
    """ Build a list of list of grids
    """

    rotor_lst_grids = ()
    for rotor in rotor_lst:
        rotor_grids = ()
        for torsion in rotor:
            rotor_grids += (
                automol.pot.grid(
                    torsion.zma, torsion.name,
                    span, torsion.symmetry, increment, from_equilibrium=True),
            )
        rotor_lst_grids += (grids,)
    if flat:
        rotor_lst_grids = tuple(chain(*rotor_lst_grids))

    return rotor_lst_grids


# Manipulate the torsion objects
def relabel_for_geometry(rotor_lst):
    """ relabel the torsion objec tto correspond with a geometry converted
        from a z-matrix
    """
    geo = automol.zmat.geometry(rotor_lst[0][0].zma)
    geo_rotor_lst = tuple(
        tuple(tors.relabel_for_geometry(torsion) for torsion in rotor)
        for rotor in rotor_lst)

    return geo, geo_rotor_lst


# I/O
def string(rotor_lst):
    """ Write a list torsions to a string
    """

    def _encode_idxs(idxs):
        if len(idxs) == 1:
            idx_str = idxs[0]+1
        else:
            idx_str = '-'.join(str(val+1) for val in idxs)
        return idx_str

    tors_dct = {}
    for rotor in rotor_lst:
        for torsion in rotor:
            _axis = tuple(torsion.axis)
            _grps = torsion.groups
            tors_dct[torsion.name] = {
                    'axis1': _axis[0],
                    'group1': _encode_idxs(_grps[0]),
                    'axis2': _axis[1],
                    'group2': _encode_idxs(_grps[1]),
                    'symmetry': torsion.symmetry,
            }

    sort_tors_dct = {}
    tors_names = sort_tors_names(list(tors_dct.keys()))
    for name in tors_names:
        sort_tors_dct[name] = tors_dct[name]

    tors_str = yaml.dump(sort_tors_dct, sort_keys=False)

    return tors_str


def from_string(tors_str):
    """ read the transformation from a string
    """

    def _decode_idxs(idxs_str):
        if isinstance(idxs_str, int):
            idxs = int(idxs_str)-1
        else:
            idxs = tuple(map(int, idxs_str.split('-')))
            idxs = (val-1 for val in idxs)
        return idxs

    inf_dct = {}

    tors_dct = yaml.load(tors_str, Loader=yaml.FullLoader)
    for name, dct in tors_dct.items():
        _axis = frozenset({dct['axis1'], dct['axis2']})
        _grps = (_decode_idxs(dct['group1']), _decode_idxs(dct['group2']))
        symm = dct['symmetry']

        inf_dct[name] = {'axis': _axis, 'groups': _grps, 'symmetry': symm}
        # torsions += (Torsion('', name, _axis, _grps, symm),)

    return inf_dct
