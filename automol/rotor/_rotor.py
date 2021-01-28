"""
 Handle rotor objects

 Rotors: (Rotor1, Rotor2, ..., RotorN)
 Rotor: (tors_obj_1, tors_obj_2, ..., tors_obj_N)
"""

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
    tors_dct = tors.torsion_dct(zma, gra, lin_keys)

    # Place the torsions into order based on rotors
    rotors = group_torsions_into_rotors(tors_dct, tors_names=tors_names)

    return rotors


def from_data(tors_dct, tors_names=None):
    """ Build the rotors objects from existing data
    """

    # assert set(tors_dct.keys) >= {'sym', 'span', 'group'}

    # Place the torsions into order based on rotors
    rotors = group_torsions_into_rotors(tors_dct, tors_names=tors_names)

    return rotors


# Getters
def names(rotor):
    """ Build a list of list of names
    """
    return tuple(torsion.name for torsion in rotor)


def axes(rotor):
    """ Build a list of list of axes(
    """
    return tuple(torsion.axis for torsion in rotor)


def groups(rotor):
    """ Build a list of list of axes(
    """
    return tuple(torsion.groups for torsion in rotor)


def symmetries(rotor):
    """ Build a list of list of axes(
    """
    return tuple(torsion.symmetry for torsion in rotor)


def names_of_list(rotor_lst):
    """ Get a flat list of names for list of rotors
    """
    _names = ()
    for rotor in rotor_lst:
        _names += tuple(torsion.name for torsion in rotor)
    return _names
