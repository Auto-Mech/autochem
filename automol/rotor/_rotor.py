"""
 Handle rotor objects

 Rotors: (Rotor1, Rotor2, ..., RotorN)
 Rotor: (tors_obj_1, tors_obj_2, ..., tors_obj_N)
"""

from itertools import chain
from automol.rotor import _tors as tors
from automol.rotor._name import group_torsions_into_rotors
from automol.rotor._util import graph_with_keys


# constructors
def from_zma(zma, tors_names=None, tors_model=None):
    """ Construct a list-of-lists of torsion objects
    """

    # Build a graph that is used to get torsion object info
    gra, lin_keys = graph_with_keys(zma)

    # Build the torsion objects
    tors_lst = tors.torsion_lst(zma, gra, lin_keys)

    # Place the torsions into order based on rotors
    rotors = group_torsions_into_rotors(
        tors_lst, name_grps=tors_names, tors_model=tors_model)

    return rotors


def from_data(zma, tors_inf_dct, tors_names=None):
    """ Build the rotors objects from existing data
    """

    tors_lst = ()
    for name, dct in tors_inf_dct.items():
        assert set(dct.keys) >= {'symmetry', 'axis', 'groups'}
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
        _names = chain(*_names)
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
